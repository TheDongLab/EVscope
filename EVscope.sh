#!/usr/bin/env bash
# ==============================================================================
# EVscope.sh: Modular RNA-seq Analysis Pipeline
# Version: 1.0.0
# Description:
#   RNA-seq pipeline with comprehensive error handling, strict
#   variable scoping, robust input validation, and production-grade QC. Features
#   include circRNA detection (CIRCexplorer2/CIRI2), contamination screening,
#   tissue deconvolution, parallel processing, and MultiQC integration.
# Architecture:
#   - Strict POSIX-compliant with Bash 4.0+ extensions
#   - Modular step-based execution with dependency tracking
#   - Idempotent operations with checkpoint recovery
#   - Comprehensive logging with structured output
# Requirements:
#   - Bash >= 4.0 (for associative arrays, BASH_REMATCH)
#   - GNU coreutils, gawk, Python 3.8+, R 4.0+
#   - Conda/Mamba for environment management
#   - See check_dependencies() for complete tool list
# Author: Yiyong Zhao, Xianjun Dong
# License: Creative Commons Attribution 4.0 International (CC BY 4.0)
# Repository: https://github.com/YiyongZhao/EVscope
# Changelog:
#   v1.0.0 - Initial release with 27-step modular workflow
# ==============================================================================
# ------------------------------------------------------------------------------
# SHELL OPTIONS AND SAFETY SETTINGS
# ------------------------------------------------------------------------------
# -E: Inherit ERR trap in functions, command substitutions, and subshells
# -e: Exit immediately on non-zero exit status (errexit)
# -u: Treat unset variables as errors (nounset) - CRITICAL for safety
# -o pipefail: Return value of pipeline is status of last failed command
# Note: -u requires explicit handling of optional variables with ${VAR:-default}
set -Eeuo pipefail
# ------------------------------------------------------------------------------
# SCRIPT METADATA AND CONSTANTS
# ------------------------------------------------------------------------------
readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
if command -v realpath &>/dev/null; then
    readonly SCRIPT_DIR="$(realpath -e "$(dirname "${BASH_SOURCE[0]}")")"
else
    readonly SCRIPT_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
fi
readonly VERSION="1.0.0"
readonly PIPELINE_BASE_DIR="$SCRIPT_DIR"
readonly MIN_BASH_VERSION="4.0"
readonly MAX_PIPELINE_STEPS=27
readonly DEFAULT_TEMP_DIR="/tmp"
# ------------------------------------------------------------------------------
# GLOBAL STATE VARIABLES
# ------------------------------------------------------------------------------
declare -i verbosity=2
declare -i thread_count=1
declare -- sample_name=""
declare -- run_steps="all"
declare -- skip_steps=""
declare -- circ_tool="both"
declare -- read_count_mode="uniq"
declare -- gDNA_correction="no"
declare -- strand="unstrand"
declare -- config_file=""
declare -- output_dir=""
declare -- fastq_read1=""
declare -- fastq_read2=""
declare -- is_paired_end="false"
declare -i featurecounts_strand=0
declare -- featurecounts_paired=""
declare -- early_log_file=""
declare -a CLEANUP_FILES=()
declare -a CLEANUP_DIRS=()
# ------------------------------------------------------------------------------
# TERMINAL COLOR DEFINITIONS
# ------------------------------------------------------------------------------
if [[ -t 2 ]]; then
    readonly C_RESET='\033[0m'
    readonly C_RED='\033[0;31m'
    readonly C_GREEN='\033[0;32m'
    readonly C_YELLOW='\033[0;33m'
    readonly C_BLUE='\033[0;34m'
    readonly C_BOLD='\033[1m'
    readonly C_DIM='\033[2m'
else
    readonly C_RESET=''
    readonly C_RED=''
    readonly C_GREEN=''
    readonly C_YELLOW=''
    readonly C_BLUE=''
    readonly C_BOLD=''
    readonly C_DIM=''
fi
# ==============================================================================
# SECTION: UTILITY FUNCTIONS
# ==============================================================================
cleanup_on_exit() {
    local exit_code=$?
    set +e
    local file
    for file in "${CLEANUP_FILES[@]:-}"; do
        [[ -f "$file" ]] && rm -f "$file" 2>/dev/null
    done
    local dir
    for dir in "${CLEANUP_DIRS[@]:-}"; do
        [[ -d "$dir" ]] && rm -rf "$dir" 2>/dev/null
    done
    set -e
    exit "$exit_code"
}
trap cleanup_on_exit EXIT INT TERM HUP
error_handler() {
    local line_num="${1:-unknown}"
    local command="${2:-unknown}"
    local exit_status="${3:-1}"
    if declare -F log &>/dev/null; then
        log 5 "FATAL" "Unexpected error at line ${line_num}: '${command}' exited with status ${exit_status}"
    else
        echo -e "${C_RED}[FATAL] Unexpected error at line ${line_num}: '${command}' exited with status ${exit_status}${C_RESET}" >&2
    fi
}
trap 'error_handler "${LINENO}" "${BASH_COMMAND}" "$?"' ERR
log() {
    local level_num="${1:-2}"
    local level_name="${2:-INFO}"
    shift 2
    local message="${*:-}"
    if (( level_num < verbosity )); then
        return 0
    fi
    local color_prefix=""
    case "$level_name" in
        DEBUG) color_prefix="$C_BLUE" ;;
        INFO)  color_prefix="$C_GREEN" ;;
        WARN)  color_prefix="$C_YELLOW" ;;
        ERROR|FATAL) color_prefix="$C_RED" ;;
        *)     color_prefix="" ;;
    esac
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    echo -e "${color_prefix}[${timestamp}] ${level_name}:${C_RESET} ${message}" >&2
    if [[ -n "${output_dir:-}" && -d "$output_dir" ]]; then
        echo "[${timestamp}] ${level_name}: ${message}" >> "${output_dir}/EVscope_pipeline.log"
    elif [[ -n "${early_log_file:-}" && -f "$early_log_file" ]]; then
        echo "[${timestamp}] ${level_name}: ${message}" >> "$early_log_file"
    fi
    return 0
}
print_help() {
    cat << EOF
${C_BOLD}EVscope RNA-seq Analysis Pipeline v${VERSION}${C_RESET}
${C_DIM}Modular workflow for comprehensive EV RNA-seq analysis${C_RESET}
${C_BOLD}USAGE:${C_RESET}
    bash ${SCRIPT_NAME} [OPTIONS] --sample_name <n> --input_fastqs <file1> [file2]
${C_BOLD}REQUIRED ARGUMENTS:${C_RESET}
    --sample_name <n>         Unique sample identifier (alphanumeric, underscore, hyphen, dot)
    --input_fastqs <files>    Input FASTQ file(s): one for SE, two for PE (space-separated)
${C_BOLD}OPTIONAL ARGUMENTS:${C_RESET}
    --threads <int>           CPU threads for parallel operations (default: 1)
    --run_steps <spec>        Steps to execute: "all" or comma-separated list/ranges
    --skip_steps <spec>       Steps to skip (same format as --run_steps)
    --circ_tool <tool>        circRNA detection tool (default: both)
    --read_count_mode <mode>  Read counting strategy (default: uniq)
    --gDNA_correction <bool>  Apply genomic DNA correction (default: no)
    --strand <strand>         Library strand orientation (default: unstrand)
    --config <path>           Path to configuration file
    -V, --verbosity <level>   Logging verbosity: 1=DEBUG to 5=FATAL
    --dry-run                 Validate inputs and show execution plan
    -h, --help                Display this help message
    -v, --version             Display version information
EOF
    exit 0
}
print_version() {
    echo "${SCRIPT_NAME} Version: ${VERSION}"
    echo "Bash Version: ${BASH_VERSION}"
    echo "Platform: $(uname -s) $(uname -r) $(uname -m)"
    exit 0
}
check_bash_version() {
    local current_version="${BASH_VERSINFO[0]}.${BASH_VERSINFO[1]}"
    local min_major min_minor current_major current_minor
    min_major="${MIN_BASH_VERSION%%.*}"
    min_minor="${MIN_BASH_VERSION##*.}"
    current_major="${BASH_VERSINFO[0]}"
    current_minor="${BASH_VERSINFO[1]}"
    if (( current_major < min_major )) || \
       (( current_major == min_major && current_minor < min_minor )); then
        echo -e "${C_RED}[FATAL] Bash version ${MIN_BASH_VERSION}+ required (found: ${current_version})${C_RESET}" >&2
        exit 1
    fi
    return 0
}
sanitize_string() {
    local input="${1:-}"
    echo "${input//[^a-zA-Z0-9_.-]/_}"
}
get_absolute_path() {
    local path="${1:-}"
    if [[ -z "$path" ]]; then
        echo ""
        return 1
    fi
    if command -v realpath &>/dev/null; then
        realpath -e "$path" 2>/dev/null || echo ""
    elif [[ -d "$path" ]]; then
        (cd -P "$path" 2>/dev/null && pwd -P) || echo ""
    elif [[ -f "$path" ]]; then
        local dir file
        dir="$(dirname "$path")"
        file="$(basename "$path")"
        (cd -P "$dir" 2>/dev/null && echo "$(pwd -P)/${file}") || echo ""
    else
        echo ""
    fi
}
check_dependencies() {
    log 2 "INFO" "Checking software dependencies..."
    local missing_deps=0
    local -a core_dependencies=(
        "fastqc" "umi_tools" "cutadapt" "STAR" "samtools" "seqtk" "bwa" "perl" "python" "Rscript"
        "featureCounts" "ribodetector_cpu" "conda"
    )
    local cmd
    for cmd in "${core_dependencies[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            log 4 "ERROR" "Required tool not found in PATH: '${cmd}'"
            missing_deps=$((missing_deps + 1))
        else
            local version_info
            version_info="$("$cmd" --version 2>&1 | head -n1 || echo "version unknown")"
            log 1 "DEBUG" "Found ${cmd}: ${version_info}"
        fi
    done
    if [[ "$circ_tool" == "CIRCexplorer2" || "$circ_tool" == "both" ]]; then
        if ! command -v CIRCexplorer2 &>/dev/null; then
            log 4 "ERROR" "CIRCexplorer2 not found (required for circ_tool=${circ_tool})"
            missing_deps=$((missing_deps + 1))
        fi
    fi
    local -a config_scripts=(
        "${BBSPLIT_SCRIPT:-}"
        "${CIRI2_PERL_SCRIPT:-}"
        "${RSEM_CALC_EXPR:-}"
    )
    local script_path
    for script_path in "${config_scripts[@]}"; do
        if [[ -n "$script_path" && ! -f "$script_path" ]]; then
            log 4 "ERROR" "Config-defined script not found: '${script_path}'"
            missing_deps=$((missing_deps + 1))
        fi
    done
    if (( missing_deps > 0 )); then
        log 5 "FATAL" "Missing ${missing_deps} required dependencies. Install them and retry."
        exit 1
    fi
    log 2 "INFO" "All software dependencies verified."
    return 0
}
check_conda_envs() {
    log 2 "INFO" "Verifying Conda environments..."
    local existing_envs
    existing_envs=$(conda env list | awk '{print $1}')
    local missing=0
    local envs_to_check=()
    [[ -n "${PICARD_ENV:-}" ]] && envs_to_check+=("$PICARD_ENV")
    [[ -n "${KRAKEN2_ENV:-}" ]] && envs_to_check+=("$KRAKEN2_ENV")
    [[ -n "${MULTIQC_ENV:-}" ]] && envs_to_check+=("$MULTIQC_ENV")
    for env in "${envs_to_check[@]}"; do
        if ! echo "$existing_envs" | grep -qw "$env"; then
            log 4 "ERROR" "Conda environment not found: $env"
            missing=$((missing + 1))
        fi
    done
    if (( missing > 0 )); then
        log 5 "FATAL" "Missing required Conda environments."
        exit 1
    fi
    log 2 "INFO" "Conda environments verification passed."
}

validate_config_vars() {
    log 2 "INFO" "Validating configuration variables..."
    local missing_vars=0
    local -a required_vars=(
        "EVscope_PATH" "HUMAN_GENOME_FASTA" "STAR_INDEX" "TOTAL_GENE_GTF"
        "PICARD_ENV" "KRAKEN2_ENV" "JAVA_MEM"
        "GENCODE_V45_REFFLAT" "GENCODE_V45_GTF" "TOTAL_GENEID_META"
    )
    local var
    for var in "${required_vars[@]}"; do
        if [[ -z "${!var:-}" ]]; then
            log 4 "ERROR" "Required configuration variable not set: '${var}'"
            missing_vars=$((missing_vars + 1))
        fi
    done
    if (( missing_vars > 0 )); then
        log 5 "FATAL" "Missing ${missing_vars} required configuration variable(s)."
        log 5 "FATAL" "Check your configuration file: ${config_file}"
        exit 1
    fi
    local -a critical_paths=(
        "$EVscope_PATH" "$HUMAN_GENOME_FASTA" "$STAR_INDEX" "$TOTAL_GENE_GTF"
    )
    local path
    for path in "${critical_paths[@]}"; do
        if [[ ! -e "$path" ]]; then
            log 4 "ERROR" "Configuration path does not exist: '${path}'"
            missing_vars=$((missing_vars + 1))
        elif [[ ! -r "$path" ]]; then
            log 4 "ERROR" "Configuration path is not readable: '${path}'"
            missing_vars=$((missing_vars + 1))
        fi
    done
    if (( missing_vars > 0 )); then
        log 5 "FATAL" "Configuration validation failed. Fix paths and retry."
        exit 1
    fi
    log 2 "INFO" "All configuration variables validated."
    return 0
}
assert_file_exists() {
    local filepath="${1:-}"
    local description="${2:-file}"
    if [[ -z "$filepath" ]]; then
        log 5 "FATAL" "Empty path provided for ${description}"
        exit 1
    fi
    if [[ ! -f "$filepath" ]]; then
        log 5 "FATAL" "Required ${description} not found: ${filepath}"
        exit 1
    fi
    if [[ ! -r "$filepath" ]]; then
        log 5 "FATAL" "Required ${description} is not readable: ${filepath}"
        exit 1
    fi
    return 0
}
assert_dir_exists() {
    local dirpath="${1:-}"
    local description="${2:-directory}"
    if [[ -z "$dirpath" ]]; then
        log 5 "FATAL" "Empty path provided for ${description}"
        exit 1
    fi
    if [[ ! -d "$dirpath" ]]; then
        log 5 "FATAL" "Required ${description} not found: ${dirpath}"
        exit 1
    fi
    if [[ ! -x "$dirpath" ]]; then
        log 5 "FATAL" "Required ${description} is not accessible: ${dirpath}"
        exit 1
    fi
    return 0
}
validate_fastq_file() {
    local fastq_path="${1:-}"
    assert_file_exists "$fastq_path" "FASTQ file"    
    local first_line
    if [[ "$fastq_path" == *.gz ]]; then
        first_line="$(zcat "$fastq_path" 2>/dev/null | head -n1)" || true
        if [[ -z "$first_line" ]]; then
            log 5 "FATAL" "FASTQ file cannot be read (corrupted gzip?): ${fastq_path}"
            exit 1
        fi
    else
        first_line="$(head -n1 "$fastq_path")"
    fi

    if [[ ! "$first_line" =~ ^@ ]]; then
        log 5 "FATAL" "FASTQ file does not start with @ header: ${fastq_path}"
        exit 1
    fi
    
    log 1 "DEBUG" "FASTQ validation passed: ${fastq_path}"
    return 0
}
run_python() {
    if [[ -n "${CORE_PYTHON_ENV:-}" ]]; then
        conda run -n "$CORE_PYTHON_ENV" python "$@"
    else
        python "$@"
    fi
}
check_system_resources() {
    log 2 "INFO" "Checking system resources..."
    if [[ -f /proc/meminfo ]]; then
        local available_mem_kb
        available_mem_kb="$(awk '/MemAvailable/{print $2}' /proc/meminfo 2>/dev/null || echo "0")"
        local available_mem_gb=$(( available_mem_kb / 1024 / 1024 ))
        if (( available_mem_gb < 8 )); then
            log 3 "WARN" "Low available memory: ${available_mem_gb}GB (recommended: 8GB+)"
        else
            log 1 "DEBUG" "Available memory: ${available_mem_gb}GB"
        fi
        if [[ -z "${JAVA_MEM:-}" ]]; then
            JAVA_MEM="$(( available_mem_kb / 1024 / 2 ))m"
            log 1 "DEBUG" "Auto-configured JAVA_MEM: ${JAVA_MEM}"
        fi
    fi
    local output_parent
    output_parent="$(dirname "${output_dir:-/tmp}")"
    if [[ -d "$output_parent" ]]; then
        local available_space_kb
        available_space_kb="$(df -k "$output_parent" 2>/dev/null | awk 'NR==2{print $4}' || echo "0")"
        local available_space_gb=$(( available_space_kb / 1024 / 1024 ))
        if (( available_space_gb < 50 )); then
            log 3 "WARN" "Low disk space: ${available_space_gb}GB (recommended: 50GB+ for RNA-seq)"
        else
            log 1 "DEBUG" "Available disk space: ${available_space_gb}GB"
        fi
    fi
    local current_ulimit
    current_ulimit="$(ulimit -n 2>/dev/null || echo "256")"
    log 1 "DEBUG" "Current ulimit -n: ${current_ulimit}"
    if (( current_ulimit < 4096 )); then
        if ulimit -n 4096 2>/dev/null; then
            log 2 "INFO" "Increased file descriptor limit: ${current_ulimit} -> $(ulimit -n)"
        else
            log 3 "WARN" "Could not increase file descriptor limit (current: ${current_ulimit})"
        fi
    fi
    return 0
}
# ==============================================================================
# SECTION: PIPELINE STEP MANAGEMENT
# ==============================================================================
declare -A STEP_DESCRIPTIONS=(
    [1]="Raw FASTQ quality control using FastQC"
    [2]="UMI motif analysis and ratio calculation"
    [3]="UMI labeling and adapter/quality trimming"
    [4]="Quality control of trimmed FASTQs"
    [5]="Bacterial contamination detection (E. coli, Mycoplasma)"
    [6]="STAR Two-Pass Alignment (Initial + Refined for circRNA)"
    [7]="Library strand detection"
    [8]="CIRCexplorer2 circRNA detection"
    [9]="CIRI2 circRNA detection"
    [10]="Merge CIRCexplorer2 and CIRI2 circRNA results"
    [11]="RNA-seq metrics collection (Picard)"
    [12]="featureCounts quantification (unique-mapping mode)"
    [13]="gDNA-corrected featureCounts quantification"
    [14]="RSEM quantification (multi-mapping mode)"
    [15]="featureCounts-based expression matrix and RNA distribution plots"
    [16]="gDNA-corrected expression matrix and RNA distribution plots"
    [17]="RSEM-based expression matrix and RNA distribution plots"
    [18]="Genomic region read mapping analysis (3'UTR, 5'UTR, etc.)"
    [19]="Taxonomic classification using Kraken2"
    [20]="Tissue deconvolution for featureCounts results"
    [21]="Tissue deconvolution for gDNA-corrected results"
    [22]="Tissue deconvolution for RSEM results"
    [23]="rRNA detection using ribodetector"
    [24]="MultiQC comprehensive QC summary"
    [25]="Coverage analysis and BigWig generation (EMapper)"
    [26]="Coverage density plots (RNA types and meta-gene regions)"
    [27]="Final HTML report generation"
)
print_pipeline_steps() {
    log 2 "INFO" "EVscope Pipeline Steps (Version: ${VERSION})"
    log 2 "INFO" "========================================"
    log 2 "INFO" "Step  | Description"
    log 2 "INFO" "------|----------------------------------------------------------------"
    local step
    for step in $(seq 1 "$MAX_PIPELINE_STEPS"); do
        log 2 "INFO" "$(printf '%-5s | %s' "$step" "${STEP_DESCRIPTIONS[$step]:-Unknown step}")"
    done
    log 2 "INFO" "========================================"
}
parse_steps() {
    local input="${1:-}"
    [[ -z "$input" ]] && return 0
    input="${input//[[:space:]]/}"
    input="${input//[\[\]]/}"
    log 1 "DEBUG" "Parsing step specification: '${input}'"
    if [[ "$input" == "all" ]]; then
        seq 1 "$MAX_PIPELINE_STEPS"
        return 0
    fi
    local -a steps=()
    local IFS=','
    local part
    for part in $input; do
        [[ -z "$part" ]] && continue
        if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            local start="${BASH_REMATCH[1]}"
            local end="${BASH_REMATCH[2]}"
            if (( start > end )); then
                local temp="$start"
                start="$end"
                end="$temp"
            fi
            if (( start < 1 || end > MAX_PIPELINE_STEPS )); then
                log 5 "FATAL" "Step range out of bounds: ${part} (valid: 1-${MAX_PIPELINE_STEPS})"
                exit 1
            fi
            local i
            for ((i = start; i <= end; i++)); do
                steps+=("$i")
            done
        elif [[ "$part" =~ ^[0-9]+$ ]]; then
            if (( part < 1 || part > MAX_PIPELINE_STEPS )); then
                log 5 "FATAL" "Step number out of bounds: ${part} (valid: 1-${MAX_PIPELINE_STEPS})"
                exit 1
            fi
            steps+=("$part")
        else
            log 5 "FATAL" "Invalid step specification format: '${part}'"
            exit 1
        fi
    done
    printf '%s\n' "${steps[@]}" | sort -n | uniq
}
run_step() {
    local step_dir="${1:-}"
    local ignore_err="${2:-false}"
    shift 2
    if [[ -z "$step_dir" ]]; then
        log 5 "FATAL" "run_step called without step directory"
        exit 1
    fi
    if ! mkdir -p "$step_dir"; then
        log 5 "FATAL" "Failed to create step directory: ${step_dir}"
        exit 1
    fi
    if [[ -f "${step_dir}/step.done" ]]; then
        log 2 "INFO" "Step '$(basename "$step_dir")' already completed. Skipping."
        return 0
    fi
    local step_name
    step_name="$(basename "$step_dir")"
    log 2 "INFO" "==> Running Step: ${step_name} <=="
    local start_time
    start_time="$(date +%s)"
    local log_file="${step_dir}/step.stderr.log"
    if ( set -eo pipefail; "$@" ) 2> "$log_file"; then
        touch "${step_dir}/step.done"
        local end_time duration
        end_time="$(date +%s)"
        duration=$((end_time - start_time))
        log 2 "INFO" "Completed step: ${step_name} in ${duration} seconds"
        return 0
    else
        local exit_status=$?
        log 4 "ERROR" "Step FAILED: ${step_name} (exit code: ${exit_status})"
        log 4 "ERROR" "Check stderr log: ${log_file}"
        if [[ -s "$log_file" ]]; then
            log 4 "ERROR" "=== Last 50 lines of stderr ==="
            tail -50 "$log_file" >&2
            log 4 "ERROR" "=== End of stderr ==="
        fi
        if [[ "$ignore_err" == "true" ]]; then
            log 3 "WARN" "Ignoring error and continuing (non-critical step)"
            touch "${step_dir}/step.done"
            return 0
        else
            log 5 "FATAL" "Pipeline stopped due to critical error in ${step_name}"
            exit "$exit_status"
        fi
    fi
}
get_circ_expr_matrix() {
    local circ_matrix=""
    case "$circ_tool" in
        both)
            circ_matrix="${output_dir}/Step_10_circRNA_Merge/${sample_name}_combined_CIRCexplorer2_CIRI2.tsv"
            ;;
        CIRCexplorer2)
            circ_matrix="${output_dir}/Step_08_CIRCexplorer2_circRNA/${sample_name}_CIRCexplorer2_dedup_junction_readcounts_CPM.tsv"
            ;;
        CIRI2)
            circ_matrix="${output_dir}/Step_09_CIRI2_circRNA/${sample_name}_CIRI2_dedup_junction_readcounts_CPM.tsv"
            ;;
    esac
    echo "$circ_matrix"
}
# ==============================================================================
# SECTION: INDIVIDUAL STEP IMPLEMENTATIONS (1-27)
# ==============================================================================
_step_1_impl() {
    local step_dir="${output_dir}/Step_01_Raw_QC"
    if [[ "$is_paired_end" == "true" ]]; then
        fastqc -o "$step_dir" -t "$thread_count" "$fastq_read1" "$fastq_read2"
    else
        fastqc -o "$step_dir" -t "$thread_count" "$fastq_read1"
    fi
}
run_step_1() {
    local step_dir="${output_dir}/Step_01_Raw_QC"
    run_step "$step_dir" "false" _step_1_impl
}
_step_2_impl() {
    local step_dir="${output_dir}/Step_02_UMI_Analysis"
    local input_fq="${fastq_read2:-${fastq_read1}}"
    run_python "${EVscope_PATH}/bin/Step_02_plot_fastq2UMI_motif.py" \
        -head "${sample_name}" -fq "${input_fq}" -n 14 -r 1000000 -o "$step_dir"
    run_python "${EVscope_PATH}/bin/Step_02_calculate_ACC_motif_fraction.py" \
        --input_fastq "${input_fq}" --positions 9-11 \
        --output_tsv "${step_dir}/${sample_name}_ACC_motif_fraction.tsv"
}
run_step_2() {
    local step_dir="${output_dir}/Step_02_UMI_Analysis"
    run_step "$step_dir" "false" _step_2_impl
}
_step_3_impl() {
    local step_dir="${output_dir}/Step_03_UMI_Adaptor_Trim"
    local r1_umi="${step_dir}/${sample_name}_R1_umi_tools.fq.gz"
    local r2_umi="${step_dir}/${sample_name}_R2_umi_tools.fq.gz"
    local r1_trim="${step_dir}/${sample_name}_R1_adapter_trimmed.fq.gz"
    local r2_trim="${step_dir}/${sample_name}_R2_adapter_trimmed.fq.gz"
    local r1_umi_trim="${step_dir}/${sample_name}_R1_adapter_UMI_trimmed.fq.gz"
    local r2_umi_trim="${step_dir}/${sample_name}_R2_adapter_UMI_trimmed.fq.gz"
    local r1_clean="${step_dir}/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${step_dir}/${sample_name}_R2_clean.fq.gz"
    if [[ "$is_paired_end" == "true" ]]; then
        umi_tools extract --bc-pattern='NNNNNNNNNNNNNN' \
            --stdin="${fastq_read2}" --stdout="$r2_umi" \
            --read2-in="${fastq_read1}" --read2-out="$r1_umi" \
            --log="${step_dir}/UMI_extract.log" --umi-separator='_'
        cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC --overlap 3 --minimum-length 10 \
            -j "$thread_count" -o "$r1_trim" -p "$r2_trim" "$r1_umi" "$r2_umi"
        run_python "${EVscope_PATH}/bin/Step_03_UMIAdapterTrimR1.py" \
            --input_R1_fq "$r1_trim" --input_R2_fq "$r2_trim" \
            --output_R1_fq "$r1_umi_trim" --output_R2_fq "$r2_umi_trim" \
            --output_tsv "${step_dir}/${sample_name}_R1_readthrough_UMI_trimming.log" \
            --min-overlap 3 --min-length 10 --chunk-size 100000 --error-rate 0.1
        cutadapt -q 20 --minimum-length 10 -j "$thread_count" \
            -o "$r1_clean" -p "$r2_clean" "$r1_umi_trim" "$r2_umi_trim"
        run_python "${EVscope_PATH}/bin/Step_03_plot_fastq_read_length_dist.py" \
            --input_fastqs "$r1_trim" "$r1_clean" "$r2_clean" \
            --titles "R1 Adapter-Trimmed" "R1 Clean" "R2 Clean" \
            --output_pdf "${step_dir}/${sample_name}_read_length_distribution.pdf" \
            --output_png "${step_dir}/${sample_name}_read_length_distribution.png"
    else
        umi_tools extract --bc-pattern='NNNNNNNNNNNNNN' \
            --stdin="${fastq_read1}" --stdout="$r1_umi" \
            --log="${step_dir}/UMI_extract.log" --umi-separator='_'
        cutadapt -a AGATCGGAAGAGC --overlap 3 --minimum-length 10 \
            -j "$thread_count" -o "$r1_trim" "$r1_umi"
        run_python "${EVscope_PATH}/bin/Step_03_UMIAdapterTrimR1.py" \
            --input_R1_fq "$r1_trim" --output_R1_fq "$r1_umi_trim" \
            --output_tsv "${step_dir}/${sample_name}_R1_readthrough_UMI_trimming.log" \
            --min-overlap 3 --min-length 10 --chunk-size 100000 --error-rate 0.1
        cutadapt -q 20 --minimum-length 10 -j "$thread_count" -o "$r1_clean" "$r1_umi_trim"
        run_python "${EVscope_PATH}/bin/Step_03_plot_fastq_read_length_dist.py" \
            --input_fastqs "$r1_trim" "$r1_clean" "$r1_clean" \
            --titles "R1 Adapter-Trimmed" "R1 Clean" "R1 Clean (dup)" \
            --output_pdf "${step_dir}/${sample_name}_read_length_distribution.pdf" \
            --output_png "${step_dir}/${sample_name}_read_length_distribution.png"
    fi
}
run_step_3() {
    local step_dir="${output_dir}/Step_03_UMI_Adaptor_Trim"
    run_step "$step_dir" "false" _step_3_impl
}
_step_4_impl() {
    local step_dir="${output_dir}/Step_04_Trimmed_QC"
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    assert_file_exists "$r1_clean" "R1 clean FASTQ from Step 3"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_clean" "R2 clean FASTQ from Step 3"
        fastqc -o "$step_dir" -t "$thread_count" "$r1_clean" "$r2_clean"
    else
        fastqc -o "$step_dir" -t "$thread_count" "$r1_clean"
    fi
}
run_step_4() {
    local step_dir="${output_dir}/Step_04_Trimmed_QC"
    run_step "$step_dir" "false" _step_4_impl
}
_step_5_impl() {
    local step_dir="${output_dir}/Step_05_Bacterial_Filter"
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    assert_file_exists "$r1_clean" "R1 clean FASTQ"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_clean" "R2 clean FASTQ"
        bash "$BBSPLIT_SCRIPT" build=1 threads="$thread_count" \
            in1="$r1_clean" in2="$r2_clean" \
            ref="${ECOLI_GENOME_FASTA:-},${MYCOPLASMA_GENOME_FASTA:-}" \
            basename="${step_dir}/${sample_name}_%_R#.fq.gz" ambiguous=best path="$step_dir"
    else
        bash "$BBSPLIT_SCRIPT" build=1 threads="$thread_count" \
            in1="$r1_clean" ref="${ECOLI_GENOME_FASTA:-},${MYCOPLASMA_GENOME_FASTA:-}" \
            basename="${step_dir}/${sample_name}_%_R#.fq.gz" ambiguous=best path="$step_dir"
    fi
}
run_step_5() {
    local step_dir="${output_dir}/Step_05_Bacterial_Filter"
    run_step "$step_dir" "false" _step_5_impl
}
_step_6_initial_impl() {
    local step_dir="${output_dir}/Step_06_Alignment_Initial"
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local initial_bam="${step_dir}/${sample_name}_Aligned.sortedByCoord.out.bam"
    local dedup_bam="${step_dir}/${sample_name}_Aligned.sortedByCoord_umi_dedup.out.bam"
    local r1_dedup="${step_dir}/${sample_name}_R1_umi_dedup.clean.fq.gz"
    local r2_dedup="${step_dir}/${sample_name}_R2_umi_dedup.clean.fq.gz"
    assert_file_exists "$r1_clean" "R1 clean FASTQ"
    mkdir -p "$step_dir"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_clean" "R2 clean FASTQ"
        STAR --genomeDir "$STAR_INDEX" --readFilesIn "$r1_clean" "$r2_clean" \
             --outFileNamePrefix "${step_dir}/${sample_name}_" \
             --runThreadN "$thread_count" --twopassMode Basic --runMode alignReads \
             --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
             --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
             --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN "$thread_count"
        samtools index -@ "$thread_count" "$initial_bam"
        umi_tools dedup -I "$initial_bam" -S "$dedup_bam" \
            --log="${step_dir}/${sample_name}_umi_dedup.log" --extract-umi-method=read_id --paired
        samtools index -@ "$thread_count" "$dedup_bam"
        samtools fastq -@ "$thread_count" -1 "$r1_dedup" -2 "$r2_dedup" \
            -0 /dev/null -s /dev/null -n "$dedup_bam"
    else
        STAR --genomeDir "$STAR_INDEX" --readFilesIn "$r1_clean" \
             --outFileNamePrefix "${step_dir}/${sample_name}_" \
             --runThreadN "$thread_count" --twopassMode Basic --runMode alignReads \
             --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
             --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
             --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN "$thread_count"
        samtools index -@ "$thread_count" "$initial_bam"
        umi_tools dedup -I "$initial_bam" -S "$dedup_bam" \
            --log="${step_dir}/${sample_name}_umi_dedup.log" --extract-umi-method=read_id
        samtools index -@ "$thread_count" "$dedup_bam"
        samtools fastq -@ "$thread_count" "$dedup_bam" | gzip > "$r1_dedup"
    fi
}
_step_6_refined_impl() {
    local step_dir="${output_dir}/Step_06_Alignment_Refined"
    local r1_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_umi_dedup.clean.fq.gz"
    local r2_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_umi_dedup.clean.fq.gz"
    local final_bam="${step_dir}/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$r1_dedup" "R1 dedup FASTQ from Step 6 initial"
    mkdir -p "$step_dir"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_dedup" "R2 dedup FASTQ from Step 6 initial"
        STAR --genomeDir "$STAR_INDEX" --readFilesIn "$r1_dedup" "$r2_dedup" \
             --outFileNamePrefix "${step_dir}/${sample_name}_STAR_umi_dedup_" \
             --runThreadN "$thread_count" --twopassMode Basic --runMode alignReads \
             --quantMode GeneCounts --readFilesCommand zcat --outFilterMultimapNmax 100 \
             --winAnchorMultimapNmax 100 --outSAMtype BAM SortedByCoordinate \
             --chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimScoreMin 1 \
             --chimOutType Junctions WithinBAM --outBAMsortingThreadN "$thread_count"
    else
        STAR --genomeDir "$STAR_INDEX" --readFilesIn "$r1_dedup" \
             --outFileNamePrefix "${step_dir}/${sample_name}_STAR_umi_dedup_" \
             --runThreadN "$thread_count" --twopassMode Basic --runMode alignReads \
             --quantMode GeneCounts --readFilesCommand zcat --outFilterMultimapNmax 100 \
             --winAnchorMultimapNmax 100 --outSAMtype BAM SortedByCoordinate \
             --chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimScoreMin 1 \
             --chimOutType Junctions WithinBAM --outBAMsortingThreadN "$thread_count"
    fi
    samtools index -@ "$thread_count" "$final_bam"
    samtools flagstat -@ "$thread_count" "$final_bam" > \
        "${step_dir}/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.flagstat"
}
_step_6_impl() {
    mkdir -p "${output_dir}/Step_06_Alignment_Initial"
    mkdir -p "${output_dir}/Step_06_Alignment_Refined"
    _step_6_initial_impl
    _step_6_refined_impl
}
run_step_6() {
    local step_dir="${output_dir}/Step_06_Alignment_Refined"
    run_step "$step_dir" "false" _step_6_impl
}
_step_7_impl() {
    local step_dir="${output_dir}/Step_07_Strand_Detection"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    run_python "${EVscope_PATH}/bin/Step_07_bam2strand.py" \
        --input_bam "$final_bam" --bed "${GENCODE_V45_non_overlapping_exon_BED:-}" \
        --test_read_num 100000000 --output_dir "$step_dir"
}
run_step_7() {
    local step_dir="${output_dir}/Step_07_Strand_Detection"
    run_step "$step_dir" "false" _step_7_impl
}
_step_8_impl() {
    local step_dir="${output_dir}/Step_08_CIRCexplorer2_circRNA"
    local chimeric="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Chimeric.out.junction"
    local bsj_bed="${step_dir}/${sample_name}_back_spliced_junction.bed"
    local known_circs="${step_dir}/${sample_name}_circularRNA_known.txt"
    local final_out="${step_dir}/${sample_name}_CIRCexplorer2_dedup_junction_readcounts_CPM.tsv"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$chimeric" "Chimeric junctions from Step 6"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    CIRCexplorer2 parse -t STAR -b "$bsj_bed" "$chimeric"
    CIRCexplorer2 annotate -r "$GENCODE_V45_REFFLAT" -g "$HUMAN_GENOME_FASTA" -b "$bsj_bed" -o "$known_circs"
    run_python "${EVscope_PATH}/bin/Step_08_convert_CIRCexplorer2CPM.py" \
        --CIRCexplorer2_result "$known_circs" --input_bam "$final_bam" \
        --GeneID_meta_table "$TOTAL_GENEID_META" --output "$final_out"
}
run_step_8() {
    if [[ "$circ_tool" == "CIRCexplorer2" || "$circ_tool" == "both" ]]; then
        local step_dir="${output_dir}/Step_08_CIRCexplorer2_circRNA"
        run_step "$step_dir" "false" _step_8_impl
    else
        log 2 "INFO" "Skipping Step 8 (CIRCexplorer2) - circ_tool is '${circ_tool}'"
    fi
}
_step_9_impl() {
    local step_dir="${output_dir}/Step_09_CIRI2_circRNA"
    local r1_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_umi_dedup.clean.fq.gz"
    local r2_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_umi_dedup.clean.fq.gz"
    local bwa_sam="${step_dir}/${sample_name}_umi_dedup.bwa.sam"
    local ciri2_out="${step_dir}/${sample_name}_CIRI2_out.tsv"
    local final_out="${step_dir}/${sample_name}_CIRI2_dedup_junction_readcounts_CPM.tsv"
    assert_file_exists "$r1_dedup" "R1 dedup FASTQ from Step 6"
    trap 'rm -f "$bwa_sam"' RETURN
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_dedup" "R2 dedup FASTQ from Step 6"
        bwa mem -t "$thread_count" -T 19 "${BWA_INDEX:-}" "$r1_dedup" "$r2_dedup" > "$bwa_sam"
    else
        bwa mem -t "$thread_count" -T 19 "${BWA_INDEX:-}" "$r1_dedup" > "$bwa_sam"
    fi
    perl "$CIRI2_PERL_SCRIPT" -T "$thread_count" -I "$bwa_sam" \
        -O "$ciri2_out" -F "$HUMAN_GENOME_FASTA" -A "$GENCODE_V45_GTF"
    run_python "${EVscope_PATH}/bin/Step_09_convert_CIRI2CPM.py" \
        --CIRI2_result "$ciri2_out" --input_sam "$bwa_sam" \
        --output "$final_out" --GeneID_meta_table "$TOTAL_GENEID_META"
}
run_step_9() {
    if [[ "$circ_tool" == "CIRI2" || "$circ_tool" == "both" ]]; then
        local step_dir="${output_dir}/Step_09_CIRI2_circRNA"
        run_step "$step_dir" "false" _step_9_impl
    else
        log 2 "INFO" "Skipping Step 9 (CIRI2) - circ_tool is '${circ_tool}'"
    fi
}
_step_10_impl() {
    local step_dir="${output_dir}/Step_10_circRNA_Merge"
    local cirexp="${output_dir}/Step_08_CIRCexplorer2_circRNA/${sample_name}_CIRCexplorer2_dedup_junction_readcounts_CPM.tsv"
    local ciri2="${output_dir}/Step_09_CIRI2_circRNA/${sample_name}_CIRI2_dedup_junction_readcounts_CPM.tsv"
    assert_file_exists "$cirexp" "CIRCexplorer2 results from Step 8"
    assert_file_exists "$ciri2" "CIRI2 results from Step 9"
    run_python "${EVscope_PATH}/bin/Step_10_circRNA_merge.py" \
        --CIRCexplorer2 "$cirexp" --CIRI2 "$ciri2" \
        --output_matrix "${step_dir}/${sample_name}_combined_CIRCexplorer2_CIRI2.tsv" \
        --out_venn "${step_dir}/${sample_name}_Venn_diagram_of_circRNAs_identified_between_CIRCexplorer2_CIRI2.png"
}
run_step_10() {
    if [[ "$circ_tool" == "both" ]]; then
        local step_dir="${output_dir}/Step_10_circRNA_Merge"
        run_step "$step_dir" "false" _step_10_impl
    else
        log 2 "INFO" "Skipping Step 10 (circRNA Merge) - circ_tool is '${circ_tool}'"
    fi
}
_step_11_impl() {
    local step_dir="${output_dir}/Step_11_RNA_Metrics"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    local picard_out="${step_dir}/${sample_name}_picard_metrics.tsv"
    local insert_out="${step_dir}/${sample_name}_insert_size_metrics.tsv"
    local insert_pdf="${step_dir}/${sample_name}_insert_size_histogram.pdf"
    local insert_png="${step_dir}/${sample_name}_insert_size_histogram.png"
    local picard_strand="NONE"
    case "$strand" in
        reverse) picard_strand="SECOND_READ_TRANSCRIPTION_STRAND" ;;
        forward) picard_strand="FIRST_READ_TRANSCRIPTION_STRAND" ;;
    esac
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    conda run -n "$PICARD_ENV" picard -Xmx${JAVA_MEM} CollectRnaSeqMetrics \
        I="$final_bam" O="$picard_out" REF_FLAT="$GENCODE_V45_REFFLAT" \
        STRAND="$picard_strand" RIBOSOMAL_INTERVALS="${HUMAN_RRNA_INTERVAL:-}"
    conda run -n "$PICARD_ENV" picard -Xmx${JAVA_MEM} CollectInsertSizeMetrics \
        I="$final_bam" O="$insert_out" H="$insert_pdf"
    if command -v convert &>/dev/null; then
        conda run -n "$PICARD_ENV" convert -density 300 -background white -alpha remove \
            "$insert_pdf" "$insert_png" || true
    fi
}
run_step_11() {
    local step_dir="${output_dir}/Step_11_RNA_Metrics"
    run_step "$step_dir" "false" _step_11_impl
}
_step_12_impl() {
    local step_dir="${output_dir}/Step_12_featureCounts_Quant"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    if [[ "$is_paired_end" == "true" ]]; then
        featureCounts -a "$TOTAL_GENE_GTF" -o "${step_dir}/${sample_name}_featureCounts.tsv" \
            -p --countReadPairs -T "$thread_count" -s "$featurecounts_strand" \
            -g gene_id -t exon -B -C "$final_bam"
    else
        featureCounts -a "$TOTAL_GENE_GTF" -o "${step_dir}/${sample_name}_featureCounts.tsv" \
            -T "$thread_count" -s "$featurecounts_strand" -g gene_id -t exon -B -C "$final_bam"
    fi
}
run_step_12() {
    local step_dir="${output_dir}/Step_12_featureCounts_Quant"
    run_step "$step_dir" "false" _step_12_impl
}
_step_13_impl() {
    local step_dir="${output_dir}/Step_13_gDNA_Corrected_Quant"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    local forward_counts="${step_dir}/${sample_name}_featureCounts_forward.tsv"
    local reverse_counts="${output_dir}/Step_12_featureCounts_Quant/${sample_name}_featureCounts.tsv"
    local corrected="${step_dir}/${sample_name}_gDNA_corrected_counts.tsv"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    assert_file_exists "$reverse_counts" "featureCounts results from Step 12"
    if [[ "$is_paired_end" == "true" ]]; then
        featureCounts -a "$TOTAL_GENE_GTF" -o "$forward_counts" \
            -p --countReadPairs -T "$thread_count" -s 1 -g gene_id -t exon -B -C "$final_bam"
    else
        featureCounts -a "$TOTAL_GENE_GTF" -o "$forward_counts" \
            -T "$thread_count" -s 1 -g gene_id -t exon -B -C "$final_bam"
    fi
    run_python "${EVscope_PATH}/bin/Step_13_gDNA_corrected_featureCounts.py" \
        --strand "$strand" --forward_featureCounts_table "$forward_counts" \
        --reverse_featureCounts_table "$reverse_counts" --output "$corrected"
}
run_step_13() {
    if [[ "$gDNA_correction" == "yes" ]]; then
        local step_dir="${output_dir}/Step_13_gDNA_Corrected_Quant"
        run_step "$step_dir" "false" _step_13_impl
    else
        log 2 "INFO" "Skipping Step 13 - gDNA_correction is 'no'"
    fi
}
_step_14_impl() {
    local step_dir="${output_dir}/Step_14_RSEM_Quant"
    local r1_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_umi_dedup.clean.fq.gz"
    local r2_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_umi_dedup.clean.fq.gz"
    assert_file_exists "$r1_dedup" "R1 dedup FASTQ from Step 6"
    mkdir -p "${step_dir}/tmp"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_dedup" "R2 dedup FASTQ from Step 6"
        perl "$RSEM_CALC_EXPR" --paired-end --bowtie2 --strandedness "$strand" \
            --bowtie2-k 2 -p "$thread_count" --no-bam-output --seed 12345 \
            "$r1_dedup" "$r2_dedup" "${RSEM_BOWTIE2_INDEX:-}" \
            "${step_dir}/${sample_name}_RSEM" --temporary-folder "${step_dir}/tmp"
    else
        perl "$RSEM_CALC_EXPR" --bowtie2 --strandedness "$strand" \
            --bowtie2-k 2 -p "$thread_count" --no-bam-output --seed 12345 \
            "$r1_dedup" "${RSEM_BOWTIE2_INDEX:-}" \
            "${step_dir}/${sample_name}_RSEM" --temporary-folder "${step_dir}/tmp"
    fi
}
run_step_14() {
    local step_dir="${output_dir}/Step_14_RSEM_Quant"
    run_step "$step_dir" "false" _step_14_impl
}
_run_expression_analysis_impl() {
    local step_dir="$1"
    local input_counts="$2"
    local script="$3"
    local suffix="$4"
    local base_matrix="${step_dir}/${sample_name}_Gene_readcounts_normalized_expression_matrix_${suffix}.tsv"
    local circ_matrix
    circ_matrix="$(get_circ_expr_matrix)"
    local combined="${step_dir}/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_${suffix}.tsv"
    assert_file_exists "$input_counts" "Input counts file"
    run_python "${EVscope_PATH}/bin/${script}" \
        --featureCounts_out "$input_counts" --GeneID_meta_table "$TOTAL_GENEID_META" --output "$base_matrix"
    if [[ -n "$circ_matrix" && -f "$circ_matrix" ]]; then
        run_python "${EVscope_PATH}/bin/Step_15_combine_total_RNA_expr_matrix.py" \
            --gene_expr "$base_matrix" --circRNA_expr "$circ_matrix" --out_matrix "$combined"
    else
        log 3 "WARN" "circRNA matrix not found, using gene expression only"
        cp "$base_matrix" "$combined"
    fi
    run_python "${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_1subplot.py" \
        --sample_name "${sample_name}" --Expr_matrix "$combined" \
        --out_plot "${step_dir}/${sample_name}_RNA_type_composition_1subplot.pdf"
    run_python "${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_2subplots.py" \
        --sample_name "${sample_name}" --Expr_matrix "$combined" \
        --out_plot "${step_dir}/${sample_name}_RNA_type_composition_2subplots.pdf"
    run_python "${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_20subplots.py" \
        --sample_name "${sample_name}" --Expr_matrix "$combined" \
        --out_plot "${step_dir}/${sample_name}_RNA_type_composition_20subplots.pdf"
    run_python "${EVscope_PATH}/bin/Step_15_plot_top_expressed_genes.py" \
        --input_gene_expr_matrix "$combined" --gene_num_per_type 5 --total_gene_num 100 \
        --output_pdf "${step_dir}/${sample_name}_bar_plot_for_top_100_highly_expressed_genes.pdf" \
        --output_png "${step_dir}/${sample_name}_bar_plot_for_top_100_highly_expressed_genes.png"
}
_step_15_impl() {
    _run_expression_analysis_impl \
        "${output_dir}/Step_15_featureCounts_Expression" \
        "${output_dir}/Step_12_featureCounts_Quant/${sample_name}_featureCounts.tsv" \
        "Step_15_featureCounts2TPM.py" "featureCounts"
}
run_step_15() {
    local step_dir="${output_dir}/Step_15_featureCounts_Expression"
    run_step "$step_dir" "false" _step_15_impl
}
_step_16_impl() {
    _run_expression_analysis_impl \
        "${output_dir}/Step_16_gDNA_Corrected_Expression" \
        "${output_dir}/Step_13_gDNA_Corrected_Quant/${sample_name}_gDNA_corrected_counts.tsv" \
        "Step_15_featureCounts2TPM.py" "gDNA_correction"
}
run_step_16() {
    if [[ "$gDNA_correction" == "yes" ]]; then
        local step_dir="${output_dir}/Step_16_gDNA_Corrected_Expression"
        run_step "$step_dir" "false" _step_16_impl
    else
        log 2 "INFO" "Skipping Step 16 - gDNA_correction is 'no'"
    fi
}
_step_17_impl() {
    local step_dir="${output_dir}/Step_17_RSEM_Expression"
    local rsem="${output_dir}/Step_14_RSEM_Quant/${sample_name}_RSEM.genes.results"
    local base_matrix="${step_dir}/${sample_name}_Gene_readcounts_normalized_expression_matrix_RSEM.tsv"
    local circ_matrix
    circ_matrix="$(get_circ_expr_matrix)"
    local combined="${step_dir}/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_RSEM.tsv"
    assert_file_exists "$rsem" "RSEM results from Step 14"
    run_python "${EVscope_PATH}/bin/Step_17_RSEM2expr_matrix.py" \
        --RSEM_out "$rsem" --GeneID_meta_table "$TOTAL_GENEID_META" --output "$base_matrix"
    if [[ -n "$circ_matrix" && -f "$circ_matrix" ]]; then
        run_python "${EVscope_PATH}/bin/Step_15_combine_total_RNA_expr_matrix.py" \
            --gene_expr "$base_matrix" --circRNA_expr "$circ_matrix" --out_matrix "$combined"
    else
        log 3 "WARN" "circRNA matrix not found, using gene expression only"
        cp "$base_matrix" "$combined"
    fi
    run_python "${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_1subplot.py" \
        --sample_name "${sample_name}" --Expr_matrix "$combined" \
        --out_plot "${step_dir}/${sample_name}_RNA_type_composition_1subplot.pdf"
    run_python "${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_2subplots.py" \
        --sample_name "${sample_name}" --Expr_matrix "$combined" \
        --out_plot "${step_dir}/${sample_name}_RNA_type_composition_2subplots.pdf"
    run_python "${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_20subplots.py" \
        --sample_name "${sample_name}" --Expr_matrix "$combined" \
        --out_plot "${step_dir}/${sample_name}_RNA_type_composition_20subplots.pdf"
    run_python "${EVscope_PATH}/bin/Step_15_plot_top_expressed_genes.py" \
        --input_gene_expr_matrix "$combined" --gene_num_per_type 5 --total_gene_num 100 \
        --output_pdf "${step_dir}/${sample_name}_bar_plot_for_top_100_highly_expressed_genes.pdf" \
        --output_png "${step_dir}/${sample_name}_bar_plot_for_top_100_highly_expressed_genes.png"
}
run_step_17() {
    local step_dir="${output_dir}/Step_17_RSEM_Expression"
    run_step "$step_dir" "false" _step_17_impl
}
# Step 18: Parallel featureCounts for genomic regions - OPTIMIZED
_step_18_impl() {
    local step_dir="${output_dir}/Step_18_Genomic_Regions"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    declare -A region_files=()
    local -a pids=()
    local saf_path saf_basename current_s_opt
    for saf_path in "${STEP18_SAF_LIST[@]:-}"; do
        if [[ ! -f "$saf_path" ]]; then
            log 3 "WARN" "SAF file not found, skipping: ${saf_path}"
            continue
        fi
        saf_basename="$(basename "$saf_path" .saf)"
        local output_file="${step_dir}/${sample_name}_${saf_basename}_featureCounts.tsv"
        region_files["$saf_basename"]="$output_file"
        current_s_opt="$featurecounts_strand"
        if [[ "$saf_basename" == *"intergenic"* || "$saf_basename" == *"blacklist"* || \
              "$saf_basename" == *"INTERGENIC"* || "$saf_basename" == *"BLACKLIST"* ]]; then
            current_s_opt=0
        fi
        (
            if [[ "$is_paired_end" == "true" ]]; then
                featureCounts -F SAF -a "$saf_path" -o "$output_file" \
                    -p --countReadPairs -T 2 -s "$current_s_opt" -B -C "$final_bam"
            else
                featureCounts -F SAF -a "$saf_path" -o "$output_file" \
                    -T 2 -s "$current_s_opt" -B -C "$final_bam"
            fi
        ) &
        pids+=($!)
# Calc max concurrent jobs (featureCounts uses 2 threads per job)
        local max_jobs=$(( thread_count / 2 ))
        [[ "$max_jobs" -lt 1 ]] && max_jobs=1

        if (( ${#pids[@]} >= max_jobs )); then
            wait "${pids[0]}"
            pids=("${pids[@]:1}")
        fi
    done
    local pid
    for pid in "${pids[@]:-}"; do
        wait "$pid" || log 3 "WARN" "A parallel featureCounts job failed"
    done
    local -a plot_args=("--sampleName" "${sample_name}" "--output_dir" "$step_dir")
    local -A region_patterns=(
        ["5UTR"]="--input_5UTR_readcounts"
        ["exon"]="--input_exon_readcounts"
        ["3UTR"]="--input_3UTR_readcounts"
        ["intron"]="--input_intron_readcounts"
        ["promoter"]="--input_promoters_readcounts"
        ["downstream"]="--input_downstream_2Kb_readcounts"
        ["intergenic"]="--input_intergenic_readcounts"
        ["blacklist"]="--input_ENCODE_blacklist_readcounts"
    )
    local pattern flag matched_file saf_bn
    for pattern in "${!region_patterns[@]}"; do
        flag="${region_patterns[$pattern]}"
        matched_file=""
        for saf_bn in "${!region_files[@]}"; do
            if [[ "$saf_bn" == *"$pattern"* ]]; then
                matched_file="${region_files[$saf_bn]}"
                break
            fi
        done
        if [[ -n "$matched_file" && -f "$matched_file" ]]; then
            plot_args+=("$flag" "$matched_file")
        fi
    done
    if (( ${#plot_args[@]} > 4 )); then
        run_python "${EVscope_PATH}/bin/Step_18_plot_reads_mapping_stats.py" "${plot_args[@]}"
    else
        log 3 "WARN" "Insufficient region files for plotting in Step 18"
    fi
}
run_step_18() {
    local step_dir="${output_dir}/Step_18_Genomic_Regions"
    run_step "$step_dir" "false" _step_18_impl
}
_step_19_impl() {
    local step_dir="${output_dir}/Step_19_Taxonomy"
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local r1_ds="${step_dir}/${sample_name}_R1_downsampled.fq.gz"
    local r2_ds="${step_dir}/${sample_name}_R2_downsampled.fq.gz"
    local kraken_report="${step_dir}/${sample_name}_report.tsv"
    local krona_input="${step_dir}/${sample_name}_krona_input.tsv"
    local krona_html="${step_dir}/${sample_name}_krona.html"
    assert_file_exists "$r1_clean" "R1 clean FASTQ"
    local random_seed=100
    local ds_count=100000
    seqtk sample -s"$random_seed" "$r1_clean" "$ds_count" | gzip > "$r1_ds"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_clean" "R2 clean FASTQ"
        seqtk sample -s"$random_seed" "$r2_clean" "$ds_count" | gzip > "$r2_ds"
        conda run -n "$KRAKEN2_ENV" kraken2 --db "${KRAKEN_DB:-}" --threads "$thread_count" \
            --report "$kraken_report" --paired --gzip-compressed "$r1_ds" "$r2_ds"
    else
        conda run -n "$KRAKEN2_ENV" kraken2 --db "${KRAKEN_DB:-}" --threads "$thread_count" \
            --report "$kraken_report" --gzip-compressed "$r1_ds"
    fi
    run_python "${KRAKEN_TOOLS_DIR:-}/kreport2krona.py" -r "$kraken_report" -o "$krona_input"
    conda run -n "$KRAKEN2_ENV" ktImportText "$krona_input" -o "$krona_html"
}
run_step_19() {
    local step_dir="${output_dir}/Step_19_Taxonomy"
    run_step "$step_dir" "true" _step_19_impl
}
_run_deconvolution_impl() {
    local step_dir="$1"
    local source_matrix="$2"
    local step_num="$3"
    local tpm_matrix="${step_dir}/${sample_name}_TPM_matrix.csv"
    if [[ ! -f "$source_matrix" ]]; then
        log 3 "WARN" "Expression matrix not found for Step ${step_num} deconvolution: ${source_matrix}"
        return 0
    fi
    awk -F'\t' '{print $1","$5}' "$source_matrix" > "$tpm_matrix"
    local -a ref_files=("${GTEX_SMTS_REF:-}" "${Monaco2020_ImmuneCell_REF:-}" "${BRAIN_SC_REF:-}")
    local ref_file
    for ref_file in "${ref_files[@]}"; do
        if [[ -n "$ref_file" && -f "$ref_file" ]]; then
            log 2 "INFO" "Running deconvolution with reference: $(basename "$ref_file")"
            run_python "${EVscope_PATH}/bin/Step_22_run_RNA_deconvolution_ARIC.py" \
                --input_expr_file "$tpm_matrix" --ref_expr_file "$ref_file" \
                --output_dir "$step_dir" --sex None || log 3 "WARN" "Deconvolution failed for $(basename "$ref_file")"
        elif [[ -n "$ref_file" ]]; then
            log 3 "WARN" "Reference file not found: '${ref_file}'. Skipping."
        fi
    done
}
_step_20_impl() {
    _run_deconvolution_impl "${output_dir}/Step_20_featureCounts_Deconvolution" \
        "${output_dir}/Step_15_featureCounts_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_featureCounts.tsv" "20"
}
run_step_20() {
    local step_dir="${output_dir}/Step_20_featureCounts_Deconvolution"
    run_step "$step_dir" "true" _step_20_impl
}
_step_21_impl() {
    _run_deconvolution_impl "${output_dir}/Step_21_gDNA_Corrected_Deconvolution" \
        "${output_dir}/Step_16_gDNA_Corrected_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_gDNA_correction.tsv" "21"
}
run_step_21() {
    if [[ "$gDNA_correction" == "yes" ]]; then
        local step_dir="${output_dir}/Step_21_gDNA_Corrected_Deconvolution"
        run_step "$step_dir" "true" _step_21_impl
    else
        log 2 "INFO" "Skipping Step 21 - gDNA_correction is 'no'"
    fi
}
_step_22_impl() {
    _run_deconvolution_impl "${output_dir}/Step_22_RSEM_Deconvolution" \
        "${output_dir}/Step_17_RSEM_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_RSEM.tsv" "22"
}
run_step_22() {
    local step_dir="${output_dir}/Step_22_RSEM_Deconvolution"
    run_step "$step_dir" "true" _step_22_impl
}
_step_23_impl() {
    local step_dir="${output_dir}/Step_23_rRNA_Detection"
    local input_r1="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local input_r2="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local output_log="${step_dir}/${sample_name}_ribodetector_summary.log"
    local rrna_r1="${step_dir}/${sample_name}_rRNA_R1.fq.gz"
    local rrna_r2="${step_dir}/${sample_name}_rRNA_R2.fq.gz"
    assert_file_exists "$input_r1" "R1 clean FASTQ"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$input_r2" "R2 clean FASTQ"
        ribodetector_cpu -t "$thread_count" -l 100 -i "$input_r1" "$input_r2" \
            -e rrna --chunk_size 800 --log "$output_log" \
            -o /dev/null /dev/null -r "$rrna_r1" "$rrna_r2" || {
            log 3 "WARN" "ribodetector may have failed. Creating empty output files."
            touch "$rrna_r1" "$rrna_r2"
        }
    else
        ribodetector_cpu -t "$thread_count" -l 100 -i "$input_r1" \
            -e rrna --chunk_size 800 --log "$output_log" \
            -o /dev/null -r "$rrna_r1" || {
            log 3 "WARN" "ribodetector may have failed. Creating empty output file."
            touch "$rrna_r1"
        }
    fi
}
run_step_23() {
    local step_dir="${output_dir}/Step_23_rRNA_Detection"
    run_step "$step_dir" "true" _step_23_impl
}
_step_24_impl() {
    local step_dir="${output_dir}/Step_24_MultiQC_Summary"
    local -a multiqc_dirs=()
    local -a qc_args=("--output" "${step_dir}/${sample_name}_QC_summary.tsv")
    local -a potential_dirs=(
        "${output_dir}/Step_01_Raw_QC"
        "${output_dir}/Step_04_Trimmed_QC"
        "${output_dir}/Step_06_Alignment_Initial"
        "${output_dir}/Step_06_Alignment_Refined"
        "${output_dir}/Step_11_RNA_Metrics"
        "${output_dir}/Step_12_featureCounts_Quant"
        "${output_dir}/Step_14_RSEM_Quant"
        "${output_dir}/Step_19_Taxonomy"
    )
    local dir
    for dir in "${potential_dirs[@]}"; do
        [[ -d "$dir" ]] && multiqc_dirs+=("$dir")
    done
    
    if (( ${#multiqc_dirs[@]} > 0 )); then
        local cmd_prefix=""
        if [[ -n "${MULTIQC_ENV:-}" ]]; then
            cmd_prefix="conda run -n ${MULTIQC_ENV}"
        fi
        $cmd_prefix multiqc \
            --title "${sample_name} EVscope QC Report" \
            --filename "${sample_name}_multiqc_report" \
            --outdir "$step_dir" --force --export --flat \
            "${multiqc_dirs[@]}" || log 3 "WARN" "MultiQC encountered issues"
    else
        log 3 "WARN" "No QC directories found for MultiQC"
    fi

    local -a raw_fqc=()
    while IFS= read -r -d '' f; do
        raw_fqc+=("$f")
    done < <(find "${output_dir}/Step_01_Raw_QC" -name "*fastqc.zip" -print0 2>/dev/null || true)
    (( ${#raw_fqc[@]} > 0 )) && qc_args+=("--raw_fastqc_zips" "${raw_fqc[@]}")
    
    local -a trimmed_fqs=()
    while IFS= read -r -d '' f; do
        trimmed_fqs+=("$f")
    done < <(find "${output_dir}/Step_03_UMI_Adaptor_Trim" -name "*_clean.fq.gz" -print0 2>/dev/null || true)
    (( ${#trimmed_fqs[@]} > 0 )) && qc_args+=("--trimmed_fastqs" "${trimmed_fqs[@]}")

    local -a ecoli_fqs=()
    while IFS= read -r -d '' f; do
        ecoli_fqs+=("$f")
    done < <(find "${output_dir}/Step_05_Bacterial_Filter" -name "*Escherichia_coli*" -name "*.fq.gz" -print0 2>/dev/null || true)
    (( ${#ecoli_fqs[@]} > 0 )) && qc_args+=("--ecoli_fastqs" "${ecoli_fqs[@]}")

    local -a myco_fqs=()
    while IFS= read -r -d '' f; do
        myco_fqs+=("$f")
    done < <(find "${output_dir}/Step_05_Bacterial_Filter" -name "*Mycoplasma*" -name "*.fq.gz" -print0 2>/dev/null || true)
    (( ${#myco_fqs[@]} > 0 )) && qc_args+=("--myco_fastqs" "${myco_fqs[@]}")

    local -a ribo_fqs=()
    while IFS= read -r -d '' f; do
        ribo_fqs+=("$f")
    done < <(find "${output_dir}/Step_23_rRNA_Detection" -name "*_rRNA_*.fq.gz" -print0 2>/dev/null || true)
    (( ${#ribo_fqs[@]} > 0 )) && qc_args+=("--ribo_fastqs" "${ribo_fqs[@]}")
    
    local star_log_refined="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Log.final.out"
    [[ -f "$star_log_refined" ]] && qc_args+=("--STAR_log" "$star_log_refined")
    
    local star_log_initial="${output_dir}/Step_06_Alignment_Initial/${sample_name}_Log.final.out"
    [[ -f "$star_log_initial" ]] && qc_args+=("--STAR_log_initial" "$star_log_initial")

    local strand_file="${output_dir}/Step_07_Strand_Detection/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out_bam2strandness.tsv"
    [[ -f "$strand_file" ]] && qc_args+=("--bam2strand_file" "$strand_file")
    
    local picard_insert="${output_dir}/Step_11_RNA_Metrics/${sample_name}_insert_size_metrics.tsv"
    [[ -f "$picard_insert" ]] && qc_args+=("--picard_insert_file" "$picard_insert")
    
    local picard_rnaseq="${output_dir}/Step_11_RNA_Metrics/${sample_name}_picard_metrics.tsv"
    [[ -f "$picard_rnaseq" ]] && qc_args+=("--picard_rnaseq_file" "$picard_rnaseq")
    
    local acc_motif="${output_dir}/Step_02_UMI_Analysis/${sample_name}_ACC_motif_fraction.tsv"
    [[ -f "$acc_motif" ]] && qc_args+=("--ACC_motif_fraction" "$acc_motif")
    
    local expr_matrix="${output_dir}/Step_15_featureCounts_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_featureCounts.tsv"
    [[ -f "$expr_matrix" ]] && qc_args+=("--expression_matrix" "$expr_matrix")
    
    local kraken_report="${output_dir}/Step_19_Taxonomy/${sample_name}_report.tsv"
    [[ -f "$kraken_report" ]] && qc_args+=("--kraken_report" "$kraken_report")

    local step18_dir="${output_dir}/Step_18_Genomic_Regions"
    local fc_3utr fc_5utr fc_downstream fc_exon fc_blacklist fc_intergenic fc_intron fc_promoter
    fc_3utr=$(find "$step18_dir" -name "*3UTR*featureCounts.tsv" 2>/dev/null | head -1)
    fc_5utr=$(find "$step18_dir" -name "*5UTR*featureCounts.tsv" 2>/dev/null | head -1)
    fc_downstream=$(find "$step18_dir" -name "*downstream*featureCounts.tsv" 2>/dev/null | head -1)
    fc_exon=$(find "$step18_dir" -name "*exon*featureCounts.tsv" 2>/dev/null | head -1)
    fc_blacklist=$(find "$step18_dir" -name "*blacklist*featureCounts.tsv" 2>/dev/null | head -1)
    fc_intergenic=$(find "$step18_dir" -name "*intergenic*featureCounts.tsv" 2>/dev/null | head -1)
    fc_intron=$(find "$step18_dir" -name "*intron*featureCounts.tsv" 2>/dev/null | head -1)
    fc_promoter=$(find "$step18_dir" -name "*promoter*featureCounts.tsv" 2>/dev/null | head -1)
    [[ -f "$fc_3utr" ]] && qc_args+=("--featureCounts_3UTR" "$fc_3utr")
    [[ -f "$fc_5utr" ]] && qc_args+=("--featureCounts_5UTR" "$fc_5utr")
    [[ -f "$fc_downstream" ]] && qc_args+=("--featureCounts_downstream_2kb" "$fc_downstream")
    [[ -f "$fc_exon" ]] && qc_args+=("--featureCounts_exon" "$fc_exon")
    [[ -f "$fc_blacklist" ]] && qc_args+=("--featureCounts_ENCODE_blacklist" "$fc_blacklist")
    [[ -f "$fc_intergenic" ]] && qc_args+=("--featureCounts_intergenic" "$fc_intergenic")
    [[ -f "$fc_intron" ]] && qc_args+=("--featureCounts_intron" "$fc_intron")
    [[ -f "$fc_promoter" ]] && qc_args+=("--featureCounts_promoter_1500_500bp" "$fc_promoter")
    
    local qc_script="${EVscope_PATH}/bin/Step_24_generate_QC_matrix.py"
    local qc_output="${step_dir}/${sample_name}_QC_summary.tsv"
    
    if [[ -f "$qc_script" ]]; then
        run_python "$qc_script" "${qc_args[@]}" || log 3 "WARN" "QC summary generation encountered issues"
    fi
    
    if [[ ! -f "$qc_output" ]]; then
        log 5 "FATAL" "QC summary file was not generated: ${qc_output}"
        exit 1
    fi
}
run_step_24() {
    local step_dir="${output_dir}/Step_24_MultiQC_Summary"
    run_step "$step_dir" "true" _step_24_impl
}
_step_25_emapper_impl() {
    local step_dir="${output_dir}/Step_25_EMapper_BigWig_Quantification/EMapper_output"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    mkdir -p "$step_dir"
    run_python "${EVscope_PATH}/bin/Step_25_EMapper.py" \
        --num_threads "$thread_count" --input_bam "$final_bam" \
        --prefix "${sample_name}" --output_dir "$step_dir" \
        --strandness "$strand" --reference_fasta "$HUMAN_GENOME_FASTA" \
        --gtf "$TOTAL_GENE_GTF" --positional_assignment EM \
        --gene_disambiguation EM --no-cleanup --polyA_tail_detection
}
_step_25_expression_impl() {
    local step_dir="${output_dir}/Step_25_EMapper_BigWig_Quantification/bigwig2expression"
    local emapper_dir="${output_dir}/Step_25_EMapper_BigWig_Quantification/EMapper_output"
    mkdir -p "$step_dir"
    local f1r2_bw="${emapper_dir}/${sample_name}_F1R2.bw"
    local f2r1_bw="${emapper_dir}/${sample_name}_F2R1.bw"
    local combined_bw="${emapper_dir}/${sample_name}_unstranded.bw"
    if [[ "$strand" == "unstrand" ]]; then
        assert_file_exists "$combined_bw" "Unstranded BigWig from EMapper"
        run_python "${EVscope_PATH}/bin/Step_25_bigWig2Expression.py" \
            --input_combined_bw "$combined_bw" --gtf "$GENCODE_V45_GTF" \
            --output "${step_dir}/${sample_name}_gene_expression.tsv"
    else
        assert_file_exists "$f1r2_bw" "F1R2 BigWig from EMapper"
        assert_file_exists "$f2r1_bw" "F2R1 BigWig from EMapper"
        run_python "${EVscope_PATH}/bin/Step_25_bigWig2Expression.py" \
            --input_F1R2_bw "$f1r2_bw" --input_F2R1_bw "$f2r1_bw" \
            --gtf "$GENCODE_V45_GTF" --output "${step_dir}/${sample_name}_gene_expression.tsv"
    fi
}
_step_25_impl() {
    mkdir -p "${output_dir}/Step_25_EMapper_BigWig_Quantification"
    _step_25_emapper_impl
    _step_25_expression_impl
}
run_step_25() {
    local step_dir="${output_dir}/Step_25_EMapper_BigWig_Quantification"
    run_step "$step_dir" "false" _step_25_impl
}
_build_bed_list_string() {
    local -n arr=$1
    local result=""
    local first=true
    local item
    for item in "${arr[@]}"; do
        if [[ "$first" == "true" ]]; then
            result="$item"
            first=false
        else
            result="${result},${item}"
        fi
    done
    echo "[$result]"
}
_step_26_impl() {
    local step_dir="${output_dir}/Step_26_BigWig_Density_Plot"
    local bigwig="${output_dir}/Step_25_EMapper_BigWig_Quantification/EMapper_output/${sample_name}_unstranded.bw"
    assert_file_exists "$bigwig" "Unstranded BigWig from Step 25"
    local bed_rna bed_meta
    bed_rna="$(_build_bed_list_string STEP26_RNATYPE_BEDS)"
    local bed_labels_rna="[${STEP26_RNATYPE_LABELS:-}]"
    bed_meta="$(_build_bed_list_string STEP26_METAGENE_BEDS)"
    local bed_labels_meta="[${STEP26_METAGENE_LABELS:-}]"
    local num_beds num_labels
    if [[ -n "${STEP26_RNATYPE_BEDS[@]+cloud}" ]]; then
    num_beds="${#STEP26_RNATYPE_BEDS[@]}"
    else
    num_beds=0
    fi
    num_labels="$(echo "${STEP26_RNATYPE_LABELS:-}" | tr ',' '\n' | grep -c . || echo 0)"
    if (( num_beds > 0 && num_labels > 0 && num_beds != num_labels )); then
        log 3 "WARN" "STEP26_RNATYPE_BEDS count (${num_beds}) != labels count (${num_labels})"
    fi
    bash "${EVscope_PATH}/bin/Step_26_density_plot_over_RNA_types.sh" \
        --input_bw_file "$bigwig" --input_bed_files "$bed_rna" \
        --input_bed_labels "$bed_labels_rna" --output_dir "${step_dir}/RNA_types" \
        --random_tested_row_num_per_bed 100000
    bash "${EVscope_PATH}/bin/Step_26_density_plot_over_meta_gene.sh" \
        --input_bw_file "$bigwig" --input_bed_files "$bed_meta" \
        --input_bed_labels "$bed_labels_meta" --output_dir "${step_dir}/meta_gene" \
        --blackListFileName "${ENCODE_BLACKLIST_BED:-}" --random_tested_row_num_per_bed 100000
}
run_step_26() {
    local step_dir="${output_dir}/Step_26_BigWig_Density_Plot"
    run_step "$step_dir" "false" _step_26_impl
}
_step_27_impl() {
    local step_dir="${output_dir}/Step_27_HTML_Report"
    local abs_output_dir
    abs_output_dir="$(get_absolute_path "$output_dir")"
    local abs_bin_dir="${EVscope_PATH}/bin"
    local abs_figures_dir="${EVscope_PATH}/figures"
    local abs_step_dir="${abs_output_dir}/Step_27_HTML_Report"
    local rmd_input="${abs_bin_dir}/Step_27_html_report.Rmd"
    local report_file="${sample_name}_final_report.html"
    local local_rmd="${abs_step_dir}/${sample_name}_report.Rmd"
    assert_file_exists "$rmd_input" "R Markdown template"
    mkdir -p "$abs_step_dir"
    cp "$rmd_input" "$local_rmd"
    export EVscope_OUTPUT_DIR="$abs_output_dir"
    export EVscope_FIGURES_DIR="$abs_figures_dir"
    cd "$abs_step_dir" || exit 1
    Rscript -e "rmarkdown::render(input = '${local_rmd}', output_file = '${report_file}', output_dir = '${abs_step_dir}', intermediates_dir = '${abs_step_dir}')"
}
run_step_27() {
    local step_dir="${output_dir}/Step_27_HTML_Report"
    run_step "$step_dir" "false" _step_27_impl
}
# ==============================================================================
# SECTION: MAIN EXECUTION LOGIC
# ==============================================================================
main() {
    early_log_file="$(mktemp /tmp/evscope_early_XXXXXX.log)"
    CLEANUP_FILES+=("$early_log_file")
    check_bash_version
    config_file="${SCRIPT_DIR}/EVscope.conf"
    local -a input_fastqs=()
    local dry_run="false"
    while (( $# > 0 )); do
        case "$1" in
            --sample_name)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--sample_name requires a value"; exit 1; }
                sample_name="$2"; shift 2 ;;
            --threads)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--threads requires a value"; exit 1; }
                [[ ! "$2" =~ ^[0-9]+$ ]] && { log 5 "FATAL" "--threads must be a positive integer"; exit 1; }
                thread_count="$2"; shift 2 ;;
            --run_steps)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--run_steps requires a value"; exit 1; }
                run_steps="$2"; shift 2 ;;
            --skip_steps)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--skip_steps requires a value"; exit 1; }
                skip_steps="$2"; shift 2 ;;
            --circ_tool)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--circ_tool requires a value"; exit 1; }
                [[ ! "$2" =~ ^(CIRCexplorer2|CIRI2|both)$ ]] && { log 5 "FATAL" "--circ_tool must be 'CIRCexplorer2', 'CIRI2', or 'both'"; exit 1; }
                circ_tool="$2"; shift 2 ;;
            --read_count_mode)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--read_count_mode requires a value"; exit 1; }
                [[ ! "$2" =~ ^(uniq|multi)$ ]] && { log 5 "FATAL" "--read_count_mode must be 'uniq' or 'multi'"; exit 1; }
                read_count_mode="$2"; shift 2 ;;
            --gDNA_correction)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--gDNA_correction requires a value"; exit 1; }
                [[ ! "$2" =~ ^(yes|no)$ ]] && { log 5 "FATAL" "--gDNA_correction must be 'yes' or 'no'"; exit 1; }
                gDNA_correction="$2"; shift 2 ;;
            --strand)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--strand requires a value"; exit 1; }
                [[ ! "$2" =~ ^(reverse|forward|unstrand)$ ]] && { log 5 "FATAL" "--strand must be 'reverse', 'forward', or 'unstrand'"; exit 1; }
                strand="$2"; shift 2 ;;
            --config)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--config requires a value"; exit 1; }
                config_file="$2"; shift 2 ;;
            -V|--verbosity)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--verbosity requires a value"; exit 1; }
                [[ ! "$2" =~ ^[1-5]$ ]] && { log 5 "FATAL" "--verbosity must be 1-5"; exit 1; }
                verbosity="$2"; shift 2 ;;
            --input_fastqs)
                shift
                while (( $# > 0 )) && [[ "$1" != -* ]]; do
                    input_fastqs+=("$1"); shift
                done ;;
            --dry-run)
                dry_run="true"; shift ;;
            -h|--help)
                print_help ;;
            -v|--version)
                print_version ;;
            *)
                log 4 "ERROR" "Unknown option: $1"
                echo "Use --help for usage information" >&2; exit 1 ;;
        esac
    done
# Logic Check: gDNA correction requires stranded data
    if [[ "$gDNA_correction" == "yes" && "$strand" == "unstrand" ]]; then
        log 5 "FATAL" "Logic Error: --gDNA_correction yes requires --strand forward/reverse."
        exit 1
    fi
    if [[ -z "$sample_name" ]]; then
        log 5 "FATAL" "--sample_name is required"
        print_help
    fi
    if (( ${#input_fastqs[@]} < 1 || ${#input_fastqs[@]} > 2 )); then
        log 5 "FATAL" "--input_fastqs requires exactly 1 (SE) or 2 (PE) files"
        log 5 "FATAL" "Provided: ${#input_fastqs[@]} files"
        exit 1
    fi
    local sample_name_sanitized
    sample_name_sanitized="$(sanitize_string "$sample_name")"
    if [[ "$sample_name" != "$sample_name_sanitized" ]]; then
        log 3 "WARN" "Sample name sanitized: '${sample_name}' -> '${sample_name_sanitized}'"
    fi
    sample_name="$sample_name_sanitized"
    output_dir="${sample_name}_EVscope_output"
    if ! mkdir -p "$output_dir"; then
        log 5 "FATAL" "Cannot create output directory: ${output_dir}"
        exit 1
    fi
    touch "${output_dir}/EVscope_pipeline.log"
    if [[ -f "$early_log_file" && -s "$early_log_file" ]]; then
        cat "$early_log_file" >> "${output_dir}/EVscope_pipeline.log"
    fi
    config_file="$(get_absolute_path "$config_file")"
    if [[ -z "$config_file" || ! -f "$config_file" ]]; then
        log 5 "FATAL" "Configuration file not found: ${config_file:-<empty>}"
        exit 1
    fi
    # shellcheck source=/dev/null
    source "$config_file"
    log 2 "INFO" "Loaded configuration from: ${config_file}"
    validate_config_vars
    check_dependencies
    check_conda_envs
    check_system_resources
    local available_cores
    available_cores="$(nproc 2>/dev/null || echo 1)"
    if (( thread_count > available_cores )); then
        log 3 "WARN" "Requested threads (${thread_count}) > available cores (${available_cores})"
        thread_count="$available_cores"
    fi
    local fastq
    for fastq in "${input_fastqs[@]}"; do
        validate_fastq_file "$fastq"
    done
    fastq_read1="$(get_absolute_path "${input_fastqs[0]}")"
    fastq_read2=""
    is_paired_end="false"
    if (( ${#input_fastqs[@]} == 2 )); then
        is_paired_end="true"
        fastq_read2="$(get_absolute_path "${input_fastqs[1]}")"
    fi
    local run_steps_list
    run_steps_list="$(parse_steps "$run_steps")"
    local skip_steps_list=""
    if [[ -n "$skip_steps" ]]; then
        skip_steps_list="$(parse_steps "$skip_steps")"
    fi
    local final_steps_list
    if [[ -n "$skip_steps_list" ]]; then
        final_steps_list="$(echo "$run_steps_list" | grep -vFx -f <(echo "$skip_steps_list") || true)"
    else
        final_steps_list="$run_steps_list"
    fi
    local -a steps_to_run=()
    local step
    while IFS= read -r step; do
        [[ -z "$step" ]] && continue
        if [[ "$gDNA_correction" == "no" ]]; then
            case "$step" in
                13|16|21) log 1 "DEBUG" "Skipping step ${step} (requires gDNA_correction=yes)"; continue ;;
            esac
        fi
        if ! declare -F "run_step_${step}" &>/dev/null; then
            log 5 "FATAL" "Step ${step} implementation not found (run_step_${step})"
            exit 1
        fi
        steps_to_run+=("$step")
    done <<< "$final_steps_list"
    if (( ${#steps_to_run[@]} == 0 )); then
        log 5 "FATAL" "No steps to run after applying filters"
        exit 1
    fi
    case "$strand" in
        reverse) featurecounts_strand=2 ;;
        forward) featurecounts_strand=1 ;;
        *)       featurecounts_strand=0 ;;
    esac
    if [[ "$is_paired_end" == "true" ]]; then
        featurecounts_paired="-p --countReadPairs"
    else
        featurecounts_paired=""
    fi
    local start_time
    start_time="$(date +%s)"
    log 2 "INFO" "============================================================"
    log 2 "INFO" "EVscope Pipeline v${VERSION}"
    log 2 "INFO" "============================================================"
    log 2 "INFO" "Sample name:      ${sample_name}"
    log 2 "INFO" "Output directory: ${output_dir}"
    log 2 "INFO" "Paired-end mode:  ${is_paired_end}"
    log 2 "INFO" "Thread count:     ${thread_count}"
    log 2 "INFO" "circRNA tool:     ${circ_tool}"
    log 2 "INFO" "Strand:           ${strand}"
    log 2 "INFO" "gDNA correction:  ${gDNA_correction}"
    log 2 "INFO" "Steps to execute: ${steps_to_run[*]}"
    log 2 "INFO" "============================================================"
    if [[ "$dry_run" == "true" ]]; then
        log 2 "INFO" "DRY-RUN MODE: Execution plan validated successfully"
        print_pipeline_steps
        log 2 "INFO" "To run the pipeline, remove --dry-run flag"
        exit 0
    fi
    print_pipeline_steps
    for step in "${steps_to_run[@]}"; do
        "run_step_${step}"
    done
    local end_time total_time
    end_time="$(date +%s)"
    total_time=$((end_time - start_time))
    local hours minutes seconds
    hours=$((total_time / 3600))
    minutes=$(((total_time % 3600) / 60))
    seconds=$((total_time % 60))
    log 2 "INFO" "============================================================"
    log 2 "INFO" "Pipeline completed successfully!"
    log 2 "INFO" "Sample:     ${sample_name}"
    log 2 "INFO" "Duration:   ${hours}h ${minutes}m ${seconds}s (${total_time} seconds)"
    log 2 "INFO" "Output:     ${output_dir}"
    log 2 "INFO" "============================================================"
    return 0
}
# ==============================================================================
# SCRIPT ENTRY POINT
# ==============================================================================
main "$@"