#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# ==============================================================================
# EVscope.sh: Modular RNA-seq Analysis Pipeline
# Version: 1.1.0
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
# License: MIT
# Repository: https://github.com/TheDongLab/EVscope
# Changelog:
#   v1.0.0 - Initial release with 27-step modular workflow
#   v1.1.0 - Add kit-aware SMARTer smRNA-seq preprocessing and no-UMI mode
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
readonly VERSION="1.1.0"
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
declare -- kit="smarter_total_rna"
declare -- alignment_tag="umi_dedup"
declare -- result_dedup_tag="dedup"
declare -- config_file=""
declare -- output_dir=""
declare -- fastq_read1=""
declare -- fastq_read2=""
declare -- is_paired_end="false"
declare -i featurecounts_strand=0
declare -- featurecounts_paired=""
declare -- early_log_file=""
declare -- resume_mode="false"
declare -- force_mode="false"
declare -a CLEANUP_FILES=()
declare -a CLEANUP_DIRS=()
declare -a STEPS_TO_RUN=()
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
    --input_fastqs <files>    Input FASTQ(s): Total RNA and smRNA accept SE or paired R1/R2
${C_BOLD}OPTIONAL ARGUMENTS:${C_RESET}
    --threads <int>           CPU threads for parallel operations (default: 1)
    --kit <kit>               Library kit: smarter_total_rna (default) or smarter_smrna
    --run_steps <spec>        Steps to execute: "all" or comma-separated list/ranges
    --skip_steps <spec>       Steps to skip (same format as --run_steps)
    --circ_tool <tool>        circRNA detection tool (default: both)
    --read_count_mode <mode>  Read counting strategy (default: uniq)
    --gDNA_correction <bool>  Apply genomic DNA correction (default: no)
    --strand <strand>         Library strand orientation (default: unstrand; smRNA: forward)
    --config <path>           Path to configuration file
    --resume                  Resume an existing output directory after metadata checks
    --force                   Re-run steps and allow overwriting an existing output directory
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
normalize_kit() {
    local value="${1:-}"
    value="${value,,}"
    value="${value//-/_}"
    case "$value" in
        smarter_total_rna|smartseq_total_rna|smartseq_total|total_rna)
            printf '%s\n' "smarter_total_rna"
            ;;
        smarter_smrna|smarter_smrna_seq|smrna|smrna_seq)
            printf '%s\n' "smarter_smrna"
            ;;
        *)
            return 1
            ;;
    esac
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
step_is_selected() {
    local target="${1:-}"
    local step
    for step in "${STEPS_TO_RUN[@]:-}"; do
        [[ "$step" == "$target" ]] && return 0
    done
    return 1
}
any_step_selected() {
    local step
    for step in "$@"; do
        step_is_selected "$step" && return 0
    done
    return 1
}
add_unique_item() {
    local -n arr_ref=$1
    local item="${2:-}"
    [[ -z "$item" ]] && return 0
    local existing
    for existing in "${arr_ref[@]:-}"; do
        [[ "$existing" == "$item" ]] && return 0
    done
    arr_ref+=("$item")
}
python_command_available() {
    command -v python &>/dev/null || command -v python3 &>/dev/null
}
python_interpreter() {
    command -v python3 2>/dev/null || command -v python 2>/dev/null || return 1
}
conda_command() {
    local candidate
    if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE:-}" ]]; then
        printf '%s\n' "$CONDA_EXE"
        return 0
    fi
    if command -v conda >/dev/null 2>&1; then
        command -v conda
        return 0
    fi
    for candidate in \
        "${EVscope_PATH:-}/../soft/Miniforge_3/bin/conda" \
        "${HOME:-}/miniforge3/bin/conda" \
        "${HOME:-}/miniconda3/bin/conda" \
        "${HOME:-}/anaconda3/bin/conda"; do
        if [[ -n "$candidate" && -x "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}
run_conda() {
    local conda_bin
    conda_bin="$(conda_command)" || { log 5 "FATAL" "Conda executable not found. Set CONDA_EXE in EVscope.conf or initialize conda in PATH."; exit 1; }
    "$conda_bin" "$@"
}
build_rscript_command() {
    local -n cmd_ref=$1
    if [[ -n "${REPORT_R_ENV:-}" ]]; then
        local conda_bin
        conda_bin="$(conda_command)" || { log 4 "ERROR" "REPORT_R_ENV=${REPORT_R_ENV} is set, but conda executable was not found. Set CONDA_EXE in EVscope.conf."; return 1; }
        cmd_ref=("$conda_bin" "run" "-n" "$REPORT_R_ENV" "Rscript")
    else
        cmd_ref=("Rscript")
    fi
}
count_csv_values() {
    local csv="${1:-}"
    csv="${csv#[}"
    csv="${csv%]}"
    csv="${csv//\"/}"
    tr ',' '\n' <<< "$csv" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' | grep -c . || true
}
check_python_modules() {
    local py
    py="$(python_interpreter)" || { log 4 "ERROR" "Required Python interpreter not found in PATH: expected 'python3' or 'python'"; return 1; }
    "$py" - <<'PYMOD'
import importlib.util
import sys
missing = [name for name in ("matplotlib", "numpy") if importlib.util.find_spec(name) is None]
if missing:
    raise SystemExit("Missing Python modules: " + ", ".join(missing))
PYMOD
}
check_r_packages() {
    local -a packages=("knitr" "rmarkdown" "readr" "dplyr" "DT" "tools" "stringr" "ggplot2" "plotly" "tidyr" "xfun" "htmltools")
    local pkg_expr
    local -a rscript_cmd=()
    build_rscript_command rscript_cmd || return 1
    pkg_expr="c($(printf '"%s",' "${packages[@]}" | sed 's/,$//'))"
    "${rscript_cmd[@]}" -e "pkgs <- ${pkg_expr}; missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop('Missing R packages: ', paste(missing, collapse=', '))" \
        || { log 4 "ERROR" "Missing R packages required for Step 27 report rendering"; return 1; }
}
check_dependencies() {
    log 2 "INFO" "Checking software dependencies for selected steps..."
    local missing_deps=0
    local -a core_dependencies=()

    # Python is used by most custom EVscope scripts. run_python falls back to python3.
    if any_step_selected 2 3 7 8 9 10 15 16 17 18 20 21 22 24 25; then
        if [[ -z "${CORE_PYTHON_ENV:-}" ]] && ! python_command_available; then
            log 4 "ERROR" "Required Python interpreter not found in PATH: expected 'python' or 'python3'"
            missing_deps=$((missing_deps + 1))
        fi
    fi
    if step_is_selected 26; then
        if ! python_command_available; then
            log 4 "ERROR" "Required Python interpreter not found in PATH for Step 26: expected 'python3' or 'python'"
            missing_deps=$((missing_deps + 1))
        else
            check_python_modules || missing_deps=$((missing_deps + 1))
        fi
    fi

    any_step_selected 1 4 && add_unique_item core_dependencies "fastqc"
    step_is_selected 3 && add_unique_item core_dependencies "cutadapt"
    if [[ "$kit" == "smarter_total_rna" ]] && any_step_selected 3 6; then
        add_unique_item core_dependencies "umi_tools"
    fi
    step_is_selected 6 && { add_unique_item core_dependencies "STAR"; add_unique_item core_dependencies "samtools"; }
    any_step_selected 7 8 11 12 13 18 25 && add_unique_item core_dependencies "samtools"
    step_is_selected 9 && { add_unique_item core_dependencies "bwa"; add_unique_item core_dependencies "perl"; }
    step_is_selected 14 && add_unique_item core_dependencies "perl"
    any_step_selected 12 13 18 && add_unique_item core_dependencies "featureCounts"
    step_is_selected 19 && add_unique_item core_dependencies "seqtk"
    step_is_selected 23 && add_unique_item core_dependencies "ribodetector_cpu"
    step_is_selected 24 && [[ -z "${MULTIQC_ENV:-}" ]] && add_unique_item core_dependencies "multiqc"
    step_is_selected 26 && { add_unique_item core_dependencies "computeMatrix"; add_unique_item core_dependencies "plotProfile"; }
    step_is_selected 27 && [[ -z "${REPORT_R_ENV:-}" ]] && add_unique_item core_dependencies "Rscript"

    local needs_conda="false"
    if [[ -n "${CORE_PYTHON_ENV:-}" ]] && any_step_selected 2 3 7 8 9 10 15 16 17 18 20 21 22 24 25; then
        needs_conda="true"
    fi
    step_is_selected 11 && [[ -n "${PICARD_ENV:-}" ]] && needs_conda="true"
    step_is_selected 19 && [[ -n "${KRAKEN2_ENV:-}" ]] && needs_conda="true"
    step_is_selected 24 && [[ -n "${MULTIQC_ENV:-}" ]] && needs_conda="true"
    step_is_selected 27 && [[ -n "${REPORT_R_ENV:-}" ]] && needs_conda="true"
    if [[ "$needs_conda" == "true" ]] && ! conda_command >/dev/null 2>&1; then
        log 4 "ERROR" "Conda executable not found. Set CONDA_EXE in EVscope.conf or initialize conda in PATH."
        missing_deps=$((missing_deps + 1))
    fi

    local cmd
    for cmd in "${core_dependencies[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            log 4 "ERROR" "Required tool not found in PATH: '${cmd}'"
            missing_deps=$((missing_deps + 1))
        else
            local version_info
            version_info="$($cmd --version 2>&1 | head -n1 || echo "version unknown")"
            log 1 "DEBUG" "Found ${cmd}: ${version_info}"
        fi
    done

    if step_is_selected 8 && [[ "$circ_tool" == "CIRCexplorer2" || "$circ_tool" == "both" ]]; then
        if ! command -v CIRCexplorer2 &>/dev/null; then
            log 4 "ERROR" "CIRCexplorer2 not found (required for circ_tool=${circ_tool})"
            missing_deps=$((missing_deps + 1))
        fi
    fi

    local -a config_scripts=()
    step_is_selected 5 && config_scripts+=("${BBSPLIT_SCRIPT:-}")
    step_is_selected 9 && config_scripts+=("${CIRI2_PERL_SCRIPT:-}")
    step_is_selected 14 && config_scripts+=("${RSEM_CALC_EXPR:-}")
    local script_path
    for script_path in "${config_scripts[@]}"; do
        if [[ -n "$script_path" && ! -f "$script_path" ]]; then
            log 4 "ERROR" "Config-defined script not found: '${script_path}'"
            missing_deps=$((missing_deps + 1))
        fi
    done

    if step_is_selected 27; then
        check_r_packages || missing_deps=$((missing_deps + 1))
    fi

    if (( missing_deps > 0 )); then
        log 5 "FATAL" "Missing ${missing_deps} required dependencies for selected steps. Install them and retry."
        exit 1
    fi
    log 2 "INFO" "Selected-step software dependencies verified."
    return 0
}
check_conda_envs() {
    log 2 "INFO" "Verifying Conda environments for selected steps..."
    local -a envs_to_check=()
    if [[ -n "${CORE_PYTHON_ENV:-}" ]] && any_step_selected 2 3 7 8 9 10 15 16 17 18 20 21 22 24 25; then
        envs_to_check+=("$CORE_PYTHON_ENV")
    fi
    step_is_selected 11 && [[ -n "${PICARD_ENV:-}" ]] && envs_to_check+=("$PICARD_ENV")
    step_is_selected 19 && [[ -n "${KRAKEN2_ENV:-}" ]] && envs_to_check+=("$KRAKEN2_ENV")
    step_is_selected 24 && [[ -n "${MULTIQC_ENV:-}" ]] && envs_to_check+=("$MULTIQC_ENV")
    step_is_selected 27 && [[ -n "${REPORT_R_ENV:-}" ]] && envs_to_check+=("$REPORT_R_ENV")
    (( ${#envs_to_check[@]} == 0 )) && { log 1 "DEBUG" "No selected-step Conda environments to verify."; return 0; }
    local conda_bin existing_envs
    conda_bin="$(conda_command)" || { log 5 "FATAL" "Conda executable not found. Set CONDA_EXE in EVscope.conf or initialize conda in PATH."; exit 1; }
    existing_envs=$("$conda_bin" env list | awk '{print $1}')
    local missing=0 env
    for env in "${envs_to_check[@]}"; do
        if ! echo "$existing_envs" | grep -qw "$env"; then
            log 4 "ERROR" "Conda environment not found: $env"
            missing=$((missing + 1))
        fi
    done
    if (( missing > 0 )); then
        log 5 "FATAL" "Missing required Conda environments for selected steps."
        exit 1
    fi
    log 2 "INFO" "Selected-step Conda environment verification passed."
}
validate_config_vars() {
    log 2 "INFO" "Validating configuration variables for selected steps..."
    local missing_vars=0
    local -a required_vars=("EVscope_PATH")
    local -a path_vars=("EVscope_PATH")
    local var path

    add_required_var() { add_unique_item required_vars "$1"; }
    add_required_path_var() { add_unique_item required_vars "$1"; add_unique_item path_vars "$1"; }

    step_is_selected 5 && add_required_path_var "BBSPLIT_SCRIPT"
    step_is_selected 6 && add_required_path_var "STAR_INDEX"
    step_is_selected 7 && add_required_path_var "GENCODE_V45_non_overlapping_exon_BED"
    step_is_selected 8 && { add_required_path_var "GENCODE_V45_REFFLAT"; add_required_path_var "HUMAN_GENOME_FASTA"; add_required_path_var "TOTAL_GENEID_META"; }
    step_is_selected 9 && { add_required_var "BWA_INDEX"; add_required_path_var "CIRI2_PERL_SCRIPT"; add_required_path_var "HUMAN_GENOME_FASTA"; add_required_path_var "GENCODE_V45_GTF"; add_required_path_var "TOTAL_GENEID_META"; }
    step_is_selected 11 && { add_required_var "PICARD_ENV"; add_required_var "JAVA_MEM"; add_required_path_var "GENCODE_V45_REFFLAT"; }
    any_step_selected 12 13 && add_required_path_var "TOTAL_GENE_GTF"
    step_is_selected 14 && { add_required_path_var "RSEM_CALC_EXPR"; add_required_var "RSEM_BOWTIE2_INDEX"; }
    any_step_selected 15 16 17 && add_required_path_var "TOTAL_GENEID_META"
    step_is_selected 18 && { add_required_path_var "STEP18_MERGED_SAF"; add_required_var "STEP18_REGION_LABELS"; }
    step_is_selected 19 && { add_required_var "KRAKEN2_ENV"; add_required_path_var "KRAKEN_DB"; add_required_path_var "KRAKEN_TOOLS_DIR"; }
    step_is_selected 25 && { add_required_path_var "HUMAN_GENOME_FASTA"; add_required_path_var "TOTAL_GENE_GTF"; add_required_path_var "GENCODE_V45_GTF"; }
    step_is_selected 26 && { add_required_var "STEP26_RNATYPE_LABELS"; add_required_var "STEP26_METAGENE_LABELS"; }

    for var in "${required_vars[@]}"; do
        if [[ -z "${!var:-}" ]]; then
            log 4 "ERROR" "Required configuration variable not set for selected steps: '${var}'"
            missing_vars=$((missing_vars + 1))
        fi
    done
    if (( missing_vars > 0 )); then
        log 5 "FATAL" "Missing ${missing_vars} required configuration variable(s)."
        log 5 "FATAL" "Check your configuration file: ${config_file}"
        exit 1
    fi

    if step_is_selected 6; then
        local declaration
        for var in STAR_COMMON_ARGS STAR_TOTAL_RNA_ARGS STAR_SMRNA_SE_ARGS STAR_SMRNA_PE_ARGS \
                   STAR_EXTRA_ARGS STAR_TOTAL_RNA_EXTRA_ARGS STAR_SMRNA_EXTRA_ARGS; do
            declaration="$(declare -p "$var" 2>/dev/null || true)"
            if [[ "$declaration" != "declare -a "* ]]; then
                log 4 "ERROR" "STAR configuration variable must be a Bash indexed array: '${var}'"
                missing_vars=$((missing_vars + 1))
            fi
        done
    fi

    for var in "${path_vars[@]}"; do
        path="${!var:-}"
        if [[ -z "$path" ]]; then
            continue
        fi
        if [[ ! -e "$path" ]]; then
            log 4 "ERROR" "Configuration path does not exist (${var}): '${path}'"
            missing_vars=$((missing_vars + 1))
        elif [[ ! -r "$path" ]]; then
            log 4 "ERROR" "Configuration path is not readable (${var}): '${path}'"
            missing_vars=$((missing_vars + 1))
        fi
    done
    if step_is_selected 26; then
        local -a step26_arrays=("STEP26_RNATYPE_BEDS" "STEP26_METAGENE_BEDS")
        local -a step26_label_vars=("STEP26_RNATYPE_LABELS" "STEP26_METAGENE_LABELS")
        local idx arr_name labels_var labels_count bed_path
        for idx in 0 1; do
            arr_name="${step26_arrays[$idx]}"
            labels_var="${step26_label_vars[$idx]}"
            if ! declare -p "$arr_name" &>/dev/null; then
                log 4 "ERROR" "Required Step 26 BED array is not defined: ${arr_name}"
                missing_vars=$((missing_vars + 1))
                continue
            fi
            local -n beds_ref="$arr_name"
            if (( ${#beds_ref[@]} == 0 )); then
                log 4 "ERROR" "Required Step 26 BED array is empty: ${arr_name}"
                missing_vars=$((missing_vars + 1))
            fi
            labels_count="$(count_csv_values "${!labels_var:-}")"
            if (( labels_count == 0 )); then
                log 4 "ERROR" "Required Step 26 label list is empty: ${labels_var}"
                missing_vars=$((missing_vars + 1))
            elif (( ${#beds_ref[@]} != labels_count )); then
                log 4 "ERROR" "Step 26 BED/label count mismatch: ${arr_name} has ${#beds_ref[@]} BEDs but ${labels_var} has ${labels_count} labels"
                missing_vars=$((missing_vars + 1))
            fi
            for bed_path in "${beds_ref[@]}"; do
                if [[ ! -f "$bed_path" ]]; then
                    log 4 "ERROR" "Step 26 BED file not found (${arr_name}): '${bed_path}'"
                    missing_vars=$((missing_vars + 1))
                elif [[ ! -r "$bed_path" ]]; then
                    log 4 "ERROR" "Step 26 BED file is not readable (${arr_name}): '${bed_path}'"
                    missing_vars=$((missing_vars + 1))
                fi
            done
            unset -n beds_ref
        done
    fi

    if (( missing_vars > 0 )); then
        log 5 "FATAL" "Configuration validation failed for selected steps. Fix paths and retry."
        exit 1
    fi
    log 2 "INFO" "Selected-step configuration variables validated."
    return 0
}

set_star_config_defaults() {
    # Backward compatibility for custom configuration files created before v1.1.0.
    declare -p STAR_COMMON_ARGS &>/dev/null || STAR_COMMON_ARGS=(
        --readFilesCommand zcat --outFilterMultimapNmax 100
        --winAnchorMultimapNmax 100 --outSAMtype BAM SortedByCoordinate
    )
    declare -p STAR_TOTAL_RNA_ARGS &>/dev/null || STAR_TOTAL_RNA_ARGS=(
        --twopassMode Basic --chimSegmentMin 10 --chimJunctionOverhangMin 10
        --chimScoreMin 1 --chimOutType Junctions WithinBAM
    )
    declare -p STAR_SMRNA_SE_ARGS &>/dev/null || STAR_SMRNA_SE_ARGS=(
        --twopassMode None --alignEndsType EndToEnd --alignIntronMax 1
        --outFilterMismatchNmax 1 --outFilterMultimapScoreRange 0
        --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0
        --outFilterScoreMin 12 --outFilterMatchNmin 14
        --chimSegmentMin 10 --chimJunctionOverhangMin 10
        --chimScoreMin 1 --chimOutType Junctions WithinBAM
    )
    declare -p STAR_SMRNA_PE_ARGS &>/dev/null || STAR_SMRNA_PE_ARGS=(
        --twopassMode None --alignEndsType EndToEnd --alignIntronMax 1
        --outFilterMismatchNmax 1 --outFilterMultimapScoreRange 0
        --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0
        --outFilterScoreMin 24 --outFilterMatchNmin 28
        --chimSegmentMin 10 --chimJunctionOverhangMin 10
        --chimScoreMin 1 --chimOutType Junctions WithinBAM
    )
    declare -p STAR_EXTRA_ARGS &>/dev/null || STAR_EXTRA_ARGS=()
    declare -p STAR_TOTAL_RNA_EXTRA_ARGS &>/dev/null || STAR_TOTAL_RNA_EXTRA_ARGS=()
    declare -p STAR_SMRNA_EXTRA_ARGS &>/dev/null || STAR_SMRNA_EXTRA_ARGS=()
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
assert_nonempty_file() {
    local filepath="${1:-}"
    local description="${2:-file}"
    assert_file_exists "$filepath" "$description"
    if [[ ! -s "$filepath" ]]; then
        log 5 "FATAL" "Required ${description} is empty: ${filepath}"
        exit 1
    fi
    return 0
}
assert_valid_gzip() {
    local filepath="${1:-}"
    local description="${2:-gzip file}"
    assert_nonempty_file "$filepath" "$description"
    if ! gzip -t "$filepath"; then
        log 5 "FATAL" "Invalid gzip ${description}: ${filepath}"
        exit 1
    fi
    return 0
}
assert_valid_bam() {
    local filepath="${1:-}"
    local description="${2:-BAM file}"
    assert_nonempty_file "$filepath" "$description"
    if ! samtools quickcheck "$filepath"; then
        log 5 "FATAL" "Invalid ${description}: ${filepath}"
        exit 1
    fi
    return 0
}
link_or_copy_file() {
    local source_path="${1:-}"
    local destination_path="${2:-}"
    assert_file_exists "$source_path" "source file for no-dedup pass-through"
    rm -f -- "$destination_path"
    if ! ln -- "$source_path" "$destination_path" 2>/dev/null; then
        cp -- "$source_path" "$destination_path"
    fi
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
json_escape() {
    local value="${1:-}"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    value="${value//$'\n'/\\n}"
    printf '%s' "$value"
}
file_fingerprint() {
    local path="${1:-}"
    if [[ -z "$path" || ! -e "$path" ]]; then
        printf 'missing:%s' "$path"
        return 0
    fi
    local abs="" size="" mtime=""
    abs="$(get_absolute_path "$path" 2>/dev/null || printf '%s' "$path")"
    size="$(stat -c '%s' "$path" 2>/dev/null || printf 'NA')"
    mtime="$(stat -c '%Y' "$path" 2>/dev/null || printf 'NA')"
    printf '%s|%s|%s' "$abs" "$size" "$mtime"
}
current_run_signature() {
    local step_name="${1:-unknown}"
    local config_hash="NA"
    [[ -n "${config_file:-}" && -f "$config_file" ]] && config_hash="$(sha256sum "$config_file" | awk '{print $1}')"
    {
        printf 'version=%s\n' "$VERSION"
        printf 'step=%s\n' "$step_name"
        printf 'sample=%s\n' "$sample_name"
        printf 'kit=%s\n' "$kit"
        printf 'is_paired_end=%s\n' "$is_paired_end"
        printf 'strand=%s\n' "$strand"
        printf 'gDNA_correction=%s\n' "$gDNA_correction"
        printf 'circ_tool=%s\n' "$circ_tool"
        printf 'read_count_mode=%s\n' "$read_count_mode"
        printf 'threads=%s\n' "$thread_count"
        printf 'config=%s\n' "$config_file"
        printf 'config_sha256=%s\n' "$config_hash"
        printf 'fastq_read1=%s\n' "$(file_fingerprint "$fastq_read1")"
        printf 'fastq_read2=%s\n' "$(file_fingerprint "$fastq_read2")"
    } | sha256sum | awk '{print $1}'
}
read_step_signature() {
    local meta_file="${1:-}"
    [[ -f "$meta_file" ]] || return 1
    sed -n 's/^[[:space:]]*"signature_sha256"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$meta_file" | head -n1
}
write_step_meta() {
    local step_dir="${1:-}"
    local step_name="${2:-}"
    local signature="${3:-}"
    local meta_file="${step_dir}/step.meta.json"
    local timestamp
    timestamp="$(date '+%Y-%m-%dT%H:%M:%S%z')"
    {
        printf '{\n'
        printf '  "version": "%s",\n' "$(json_escape "$VERSION")"
        printf '  "step": "%s",\n' "$(json_escape "$step_name")"
        printf '  "sample_name": "%s",\n' "$(json_escape "$sample_name")"
        printf '  "signature_sha256": "%s",\n' "$(json_escape "$signature")"
        printf '  "timestamp": "%s",\n' "$(json_escape "$timestamp")"
        printf '  "config_file": "%s",\n' "$(json_escape "$config_file")"
        printf '  "fastq_read1": "%s",\n' "$(json_escape "$fastq_read1")"
        printf '  "fastq_read2": "%s",\n' "$(json_escape "$fastq_read2")"
        printf '  "parameters": {"kit": "%s", "strand": "%s", "gDNA_correction": "%s", "circ_tool": "%s", "read_count_mode": "%s", "threads": "%s"}\n' \
            "$(json_escape "$kit")" "$(json_escape "$strand")" "$(json_escape "$gDNA_correction")" "$(json_escape "$circ_tool")" "$(json_escape "$read_count_mode")" "$(json_escape "$thread_count")"
        printf '}\n'
    } > "$meta_file"
}
find_first_sample_file() {
    local dir="${1:-}"
    local pattern="${2:-}"
    [[ -d "$dir" ]] || return 0
    local -a matches=()
    while IFS= read -r -d '' f; do
        matches+=("$f")
    done < <(find "$dir" -maxdepth 1 -type f -name "${sample_name}${pattern}" -print0 2>/dev/null | sort -z)
    if (( ${#matches[@]} > 1 )); then
        log 3 "WARN" "Multiple files matched ${sample_name}${pattern} in ${dir}; using $(basename "${matches[0]}")"
    fi
    (( ${#matches[@]} > 0 )) && printf '%s' "${matches[0]}"
}

first_fastq_id() {
    local fastq_path="${1:-}"
    local header=""
    if [[ "$fastq_path" == *.gz ]]; then
        header="$(zcat "$fastq_path" 2>/dev/null | head -n1 || true)"
    else
        header="$(head -n1 "$fastq_path" 2>/dev/null || true)"
    fi
    header="${header#@}"
    header="${header%%[[:space:]]*}"
    header="${header%/1}"
    header="${header%/2}"
    printf '%s' "$header"
}
validate_fastq_pair_order() {
    local r1="${1:-}"
    local r2="${2:-}"
    local b1 b2 id1 id2
    b1="$(basename "$r1")"
    b2="$(basename "$r2")"
    if [[ "$b1" =~ (^|[^A-Za-z0-9])R?2([^A-Za-z0-9]|$) && "$b2" =~ (^|[^A-Za-z0-9])R?1([^A-Za-z0-9]|$) ]]; then
        log 5 "FATAL" "FASTQ mate order appears reversed: first file looks like R2 (${b1}), second file looks like R1 (${b2})"
        exit 1
    fi
    if [[ ! "$b1" =~ (^|[^A-Za-z0-9])R?1([^A-Za-z0-9]|$) || ! "$b2" =~ (^|[^A-Za-z0-9])R?2([^A-Za-z0-9]|$) ]]; then
        log 3 "WARN" "Could not confidently infer R1/R2 from FASTQ basenames (${b1}, ${b2}); validating first read IDs only."
    fi
    id1="$(first_fastq_id "$r1")"
    id2="$(first_fastq_id "$r2")"
    if [[ -n "$id1" && -n "$id2" && "$id1" != "$id2" ]]; then
        log 5 "FATAL" "FASTQ mate headers do not match between R1 and R2: '${id1}' vs '${id2}'"
        exit 1
    fi
}

run_python() {
    if [[ -n "${CORE_PYTHON_ENV:-}" ]]; then
        run_conda run -n "$CORE_PYTHON_ENV" python "$@"
    elif command -v python &>/dev/null; then
        python "$@"
    else
        python3 "$@"
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
    [27]="Detailed and concise standalone HTML report generation"
)
apply_kit_profile() {
    local strand_was_set="${1:-false}"
    case "$kit" in
        smarter_total_rna)
            alignment_tag="umi_dedup"
            result_dedup_tag="dedup"
            ;;
        smarter_smrna)
            alignment_tag="no_dedup"
            result_dedup_tag="no_dedup"
            if [[ "$strand_was_set" != "true" ]]; then
                strand="forward"
            fi
            STEP_DESCRIPTIONS[2]="SMARTer smRNA library-structure profile (no UMI)"
            STEP_DESCRIPTIONS[3]="SMARTer smRNA 5-prime/poly(A) and quality trimming"
            STEP_DESCRIPTIONS[6]="STAR alignment without UMI deduplication"
            ;;
        *)
            log 5 "FATAL" "Internal error: unsupported normalized kit '${kit}'"
            exit 1
            ;;
    esac
}
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
                log 5 "FATAL" "Step range must be ascending: ${part}"
                exit 1
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
    local step_name
    step_name="$(basename "$step_dir")"
    local signature existing_signature
    signature="$(current_run_signature "$step_name")"
    if [[ "$force_mode" != "true" && -f "${step_dir}/step.done" && ! -f "${step_dir}/step.failed" ]]; then
        existing_signature="$(read_step_signature "${step_dir}/step.meta.json" || true)"
        if [[ -z "$existing_signature" ]]; then
            log 5 "FATAL" "Step '${step_name}' has step.done but no metadata; use --force to regenerate or clean the step directory."
            exit 1
        fi
        if [[ "$existing_signature" != "$signature" ]]; then
            log 5 "FATAL" "Step '${step_name}' metadata does not match current inputs/config/parameters; use --force to regenerate."
            exit 1
        fi
        log 2 "INFO" "Step '${step_name}' already completed with matching metadata. Skipping."
        return 0
    fi
    if [[ -f "${step_dir}/step.done" && -f "${step_dir}/step.failed" ]]; then
        log 3 "WARN" "Step '${step_name}' has both step.done and step.failed; rerunning."
    elif [[ -f "${step_dir}/step.failed" ]]; then
        log 3 "WARN" "Step '${step_name}' previously failed; rerunning."
    fi
    rm -f "${step_dir}/step.done" "${step_dir}/step.failed"
    log 2 "INFO" "==> Running Step: ${step_name} <=="
    local start_time
    start_time="$(date +%s)"
    local log_file="${step_dir}/step.stderr.log"
    if ( set -eo pipefail; "$@" ) 2> "$log_file"; then
        rm -f "${step_dir}/step.failed"
        write_step_meta "$step_dir" "$step_name" "$signature"
        touch "${step_dir}/step.done"
        local end_time duration
        end_time="$(date +%s)"
        duration=$((end_time - start_time))
        log 2 "INFO" "Completed step: ${step_name} in ${duration} seconds"
        return 0
    else
        local exit_status=$?
        rm -f "${step_dir}/step.done"
        touch "${step_dir}/step.failed"
        log 4 "ERROR" "Step FAILED: ${step_name} (exit code: ${exit_status})"
        log 4 "ERROR" "Check stderr log: ${log_file}"
        if [[ -s "$log_file" ]]; then
            log 4 "ERROR" "=== Last 50 lines of stderr ==="
            tail -50 "$log_file" >&2
            log 4 "ERROR" "=== End of stderr ==="
        fi
        if [[ "$ignore_err" == "true" ]]; then
            log 3 "WARN" "Ignoring error and continuing (non-critical step; step.failed retained)"
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
            circ_matrix="${output_dir}/Step_08_CIRCexplorer2_circRNA/${sample_name}_CIRCexplorer2_${result_dedup_tag}_junction_readcounts_CPM.tsv"
            ;;
        CIRI2)
            circ_matrix="${output_dir}/Step_09_CIRI2_circRNA/${sample_name}_CIRI2_${result_dedup_tag}_junction_readcounts_CPM.tsv"
            ;;
    esac
    echo "$circ_matrix"
}
# ==============================================================================
# SECTION: INDIVIDUAL STEP IMPLEMENTATIONS (1-27)
# ==============================================================================
_step_1_impl() {
    local step_dir="${output_dir}/Step_01_Raw_QC"
    local expected_count=1
    if [[ "$is_paired_end" == "true" ]]; then
        expected_count=2
        fastqc -o "$step_dir" -t "$thread_count" "$fastq_read1" "$fastq_read2"
    else
        fastqc -o "$step_dir" -t "$thread_count" "$fastq_read1"
    fi
    local zip_count html_count
    zip_count="$(find "$step_dir" -maxdepth 1 -type f -name '*_fastqc.zip' -size +0c | wc -l)"
    html_count="$(find "$step_dir" -maxdepth 1 -type f -name '*_fastqc.html' -size +0c | wc -l)"
    if (( zip_count < expected_count || html_count < expected_count )); then
        log 5 "FATAL" "FastQC output incomplete in ${step_dir}: expected ${expected_count} zip/html files, found zip=${zip_count}, html=${html_count}"
        exit 1
    fi
}
run_step_1() {
    local step_dir="${output_dir}/Step_01_Raw_QC"
    run_step "$step_dir" "false" _step_1_impl
}
_step_2_impl() {
    local step_dir="${output_dir}/Step_02_UMI_Analysis"
    if [[ "$kit" == "smarter_smrna" ]]; then
        local input_layout="Single-end Read 1"
        [[ "$is_paired_end" == "true" ]] && input_layout="Paired-end Read 1 / Read 2"
        {
            printf 'Parameter\tValue\n'
            printf 'Kit\tSMARTer smRNA-Seq Kit for Illumina\n'
            printf 'Input layout\t%s\n' "$input_layout"
            printf 'Read orientation\tSense\n'
            printf 'UMI\tNo\n'
            printf 'R1 structure\t3-nt template-switch sequence + RNA insert + poly(A)\n'
            printf '5-prime fixed trim\t3 nt\n'
            printf '3-prime adapter\tAAAAAAAAAA\n'
            if [[ "$is_paired_end" == "true" ]]; then
                printf 'R2 structure\t14-nt poly(T) + reverse-complement RNA insert + 3-nt template-switch sequence + Illumina adapter\n'
                printf 'R2 5-prime fixed trim\t14 nt\n'
                printf 'R2 3-prime adapter\tNNNAGATCGGAAGAGC\n'
            fi
            printf 'Minimum retained length\t15 nt\n'
            printf 'Molecular deduplication\tDisabled\n'
        } > "${step_dir}/${sample_name}_library_profile.tsv"
        return 0
    fi

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

    if [[ "$kit" == "smarter_smrna" ]]; then
        # Takara's recommended structure trimming is: -m 15 -u 3 -a AAAAAAAAAA.
        # Quality trimming remains a separate pass so the length plot can distinguish
        # protocol-derived trimming from the pipeline's general Q20 filtering.
        if [[ "$is_paired_end" == "true" ]]; then
            cutadapt -u 3 -U 14 \
                -a 'AAAAAAAAAA;min_overlap=3' \
                -A 'NNNAGATCGGAAGAGC;min_overlap=10' \
                --minimum-length 15 -j "$thread_count" \
                --json "${step_dir}/${sample_name}_structure_trimming.cutadapt.json" \
                -o "$r1_trim" -p "$r2_trim" "$fastq_read1" "$fastq_read2"
            cutadapt -q 20 -Q 20 --minimum-length 15 -j "$thread_count" \
                --json "${step_dir}/${sample_name}_quality_trimming.cutadapt.json" \
                -o "$r1_clean" -p "$r2_clean" "$r1_trim" "$r2_trim"
            run_python "${EVscope_PATH}/bin/Step_03_plot_fastq_read_length_dist.py" \
                --input_fastqs "$r1_trim" "$r1_clean" "$r2_clean" \
                --titles "R1 Structure-Trimmed" "R1 Clean" "R2 Clean" \
                --output_pdf "${step_dir}/${sample_name}_read_length_distribution.pdf" \
                --output_png "${step_dir}/${sample_name}_read_length_distribution.png"
            assert_valid_gzip "$r2_clean" "SMARTer smRNA clean R2 FASTQ"
        else
            cutadapt -u 3 -a AAAAAAAAAA --overlap 3 --minimum-length 15 \
                -j "$thread_count" --json "${step_dir}/${sample_name}_structure_trimming.cutadapt.json" \
                -o "$r1_trim" "$fastq_read1"
            cutadapt -q 20 --minimum-length 15 -j "$thread_count" \
                --json "${step_dir}/${sample_name}_quality_trimming.cutadapt.json" \
                -o "$r1_clean" "$r1_trim"
            run_python "${EVscope_PATH}/bin/Step_03_plot_fastq_read_length_dist.py" \
                --input_fastqs "$fastq_read1" "$r1_trim" "$r1_clean" \
                --titles "Raw R1" "Structure-Trimmed R1" "Clean R1" \
                --output_pdf "${step_dir}/${sample_name}_read_length_distribution.pdf" \
                --output_png "${step_dir}/${sample_name}_read_length_distribution.png"
        fi
        assert_valid_gzip "$r1_clean" "SMARTer smRNA clean R1 FASTQ"
        return 0
    fi

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
configure_star_protocol_args() {
    if [[ "$kit" == "smarter_smrna" ]]; then
        if [[ "$is_paired_end" == "true" ]]; then
            STAR_PROTOCOL_ARGS=("${STAR_COMMON_ARGS[@]}" "${STAR_SMRNA_PE_ARGS[@]}" \
                "${STAR_EXTRA_ARGS[@]}" "${STAR_SMRNA_EXTRA_ARGS[@]}")
        else
            STAR_PROTOCOL_ARGS=("${STAR_COMMON_ARGS[@]}" "${STAR_SMRNA_SE_ARGS[@]}" \
                "${STAR_EXTRA_ARGS[@]}" "${STAR_SMRNA_EXTRA_ARGS[@]}")
        fi
    else
        STAR_PROTOCOL_ARGS=("${STAR_COMMON_ARGS[@]}" "${STAR_TOTAL_RNA_ARGS[@]}" \
            "${STAR_EXTRA_ARGS[@]}" "${STAR_TOTAL_RNA_EXTRA_ARGS[@]}")
    fi
}
_step_6_initial_impl() {
    local step_dir="${output_dir}/Step_06_Alignment_Initial"
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local initial_bam="${step_dir}/${sample_name}_Aligned.sortedByCoord.out.bam"
    local processed_bam="${step_dir}/${sample_name}_Aligned.sortedByCoord_${alignment_tag}.out.bam"
    local r1_processed="${step_dir}/${sample_name}_R1_${alignment_tag}.clean.fq.gz"
    local r2_processed="${step_dir}/${sample_name}_R2_${alignment_tag}.clean.fq.gz"
    assert_file_exists "$r1_clean" "R1 clean FASTQ"
    mkdir -p "$step_dir"
    configure_star_protocol_args
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_clean" "R2 clean FASTQ"
        STAR --genomeDir "$STAR_INDEX" --readFilesIn "$r1_clean" "$r2_clean" \
             --outFileNamePrefix "${step_dir}/${sample_name}_" \
             --runThreadN "$thread_count" --runMode alignReads \
             --outBAMsortingThreadN "$thread_count" \
             "${STAR_PROTOCOL_ARGS[@]}"
        samtools index -@ "$thread_count" "$initial_bam"
        assert_valid_bam "$initial_bam" "initial STAR BAM"
        assert_nonempty_file "${initial_bam}.bai" "initial STAR BAM index"
        if [[ "$kit" == "smarter_total_rna" ]]; then
            umi_tools dedup -I "$initial_bam" -S "$processed_bam" \
                --log="${step_dir}/${sample_name}_umi_dedup.log" --extract-umi-method=read_id --paired
            samtools index -@ "$thread_count" "$processed_bam"
            samtools fastq -@ "$thread_count" -1 "$r1_processed" -2 "$r2_processed" \
                -0 /dev/null -s /dev/null -n "$processed_bam"
        else
            link_or_copy_file "$initial_bam" "$processed_bam"
            link_or_copy_file "${initial_bam}.bai" "${processed_bam}.bai"
            link_or_copy_file "$r1_clean" "$r1_processed"
            link_or_copy_file "$r2_clean" "$r2_processed"
        fi
    else
        STAR --genomeDir "$STAR_INDEX" --readFilesIn "$r1_clean" \
             --outFileNamePrefix "${step_dir}/${sample_name}_" \
             --runThreadN "$thread_count" --runMode alignReads \
             --outBAMsortingThreadN "$thread_count" \
             "${STAR_PROTOCOL_ARGS[@]}"
        samtools index -@ "$thread_count" "$initial_bam"
        assert_valid_bam "$initial_bam" "initial STAR BAM"
        assert_nonempty_file "${initial_bam}.bai" "initial STAR BAM index"
        if [[ "$kit" == "smarter_total_rna" ]]; then
            umi_tools dedup -I "$initial_bam" -S "$processed_bam" \
                --log="${step_dir}/${sample_name}_umi_dedup.log" --extract-umi-method=read_id
            samtools index -@ "$thread_count" "$processed_bam"
            samtools fastq -@ "$thread_count" "$processed_bam" | gzip > "$r1_processed"
        else
            link_or_copy_file "$initial_bam" "$processed_bam"
            link_or_copy_file "${initial_bam}.bai" "${processed_bam}.bai"
            link_or_copy_file "$r1_clean" "$r1_processed"
        fi
    fi
    assert_valid_bam "$processed_bam" "kit-processed BAM"
    assert_nonempty_file "${processed_bam}.bai" "kit-processed BAM index"
    assert_valid_gzip "$r1_processed" "R1 kit-processed FASTQ"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_valid_gzip "$r2_processed" "R2 kit-processed FASTQ"
    fi
}
_step_6_refined_impl() {
    local step_dir="${output_dir}/Step_06_Alignment_Refined"
    local r1_processed="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_${alignment_tag}.clean.fq.gz"
    local r2_processed="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_${alignment_tag}.clean.fq.gz"
    local final_bam="${step_dir}/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$r1_processed" "R1 processed FASTQ from Step 6 initial"
    mkdir -p "$step_dir"
    configure_star_protocol_args
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_processed" "R2 processed FASTQ from Step 6 initial"
        STAR --genomeDir "$STAR_INDEX" --readFilesIn "$r1_processed" "$r2_processed" \
             --outFileNamePrefix "${step_dir}/${sample_name}_STAR_${alignment_tag}_" \
             --runThreadN "$thread_count" --runMode alignReads \
             --quantMode GeneCounts --outBAMsortingThreadN "$thread_count" \
             "${STAR_PROTOCOL_ARGS[@]}"
    else
        STAR --genomeDir "$STAR_INDEX" --readFilesIn "$r1_processed" \
             --outFileNamePrefix "${step_dir}/${sample_name}_STAR_${alignment_tag}_" \
             --runThreadN "$thread_count" --runMode alignReads \
             --quantMode GeneCounts --outBAMsortingThreadN "$thread_count" \
             "${STAR_PROTOCOL_ARGS[@]}"
    fi
    samtools index -@ "$thread_count" "$final_bam"
    samtools flagstat -@ "$thread_count" "$final_bam" > \
        "${step_dir}/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.flagstat"
    assert_valid_bam "$final_bam" "refined STAR BAM"
    assert_nonempty_file "${final_bam}.bai" "refined STAR BAM index"
    assert_nonempty_file "${step_dir}/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.flagstat" "refined STAR flagstat"
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
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    local star_log_refined="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Log.final.out"
    local star_log_args=()
    [[ -f "$star_log_refined" ]] && star_log_args=(--star_log "$star_log_refined")
    run_python "${EVscope_PATH}/bin/Step_07_bam2strand.py" \
        --input_bam "$final_bam" --bed "${GENCODE_V45_non_overlapping_exon_BED:-}" \
        --test_read_num 100000000 --output_dir "$step_dir" "${star_log_args[@]}"
}
run_step_7() {
    local step_dir="${output_dir}/Step_07_Strand_Detection"
    run_step "$step_dir" "false" _step_7_impl
}
_step_8_impl() {
    local step_dir="${output_dir}/Step_08_CIRCexplorer2_circRNA"
    local chimeric="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Chimeric.out.junction"
    local bsj_bed="${step_dir}/${sample_name}_back_spliced_junction.bed"
    local known_circs="${step_dir}/${sample_name}_circularRNA_known.txt"
    local final_out="${step_dir}/${sample_name}_CIRCexplorer2_${result_dedup_tag}_junction_readcounts_CPM.tsv"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.bam"
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
    local r1_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_${alignment_tag}.clean.fq.gz"
    local r2_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_${alignment_tag}.clean.fq.gz"
    local bwa_sam="${step_dir}/${sample_name}_${alignment_tag}.bwa.sam"
    local ciri2_out="${step_dir}/${sample_name}_CIRI2_out.tsv"
    local final_out="${step_dir}/${sample_name}_CIRI2_${result_dedup_tag}_junction_readcounts_CPM.tsv"
    assert_valid_gzip "$r1_dedup" "R1 processed FASTQ from Step 6"
    trap 'rm -f "$bwa_sam"' RETURN
    if [[ "$is_paired_end" == "true" ]]; then
        assert_valid_gzip "$r2_dedup" "R2 processed FASTQ from Step 6"
        bwa mem -t "$thread_count" -T 19 "${BWA_INDEX:-}" "$r1_dedup" "$r2_dedup" > "$bwa_sam"
    else
        bwa mem -t "$thread_count" -T 19 "${BWA_INDEX:-}" "$r1_dedup" > "$bwa_sam"
    fi
    assert_nonempty_file "$bwa_sam" "BWA SAM for CIRI2"
    perl "$CIRI2_PERL_SCRIPT" -T "$thread_count" -I "$bwa_sam" \
        -O "$ciri2_out" -F "$HUMAN_GENOME_FASTA" -A "$GENCODE_V45_GTF"
    assert_nonempty_file "$ciri2_out" "CIRI2 raw output"
    run_python "${EVscope_PATH}/bin/Step_09_convert_CIRI2CPM.py" \
        --CIRI2_result "$ciri2_out" --input_sam "$bwa_sam" \
        --output "$final_out" --GeneID_meta_table "$TOTAL_GENEID_META"
    assert_nonempty_file "$final_out" "CIRI2 CPM output"
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
    local cirexp="${output_dir}/Step_08_CIRCexplorer2_circRNA/${sample_name}_CIRCexplorer2_${result_dedup_tag}_junction_readcounts_CPM.tsv"
    local ciri2="${output_dir}/Step_09_CIRI2_circRNA/${sample_name}_CIRI2_${result_dedup_tag}_junction_readcounts_CPM.tsv"
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
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.bam"
    local picard_out="${step_dir}/${sample_name}_picard_metrics.tsv"
    local insert_out="${step_dir}/${sample_name}_insert_size_metrics.tsv"
    local insert_pdf="${step_dir}/${sample_name}_insert_size_histogram.pdf"
    local insert_png="${step_dir}/${sample_name}_insert_size_histogram.png"
    local picard_strand="NONE"
    case "$strand" in
        reverse) picard_strand="SECOND_READ_TRANSCRIPTION_STRAND" ;;
        forward) picard_strand="FIRST_READ_TRANSCRIPTION_STRAND" ;;
    esac
    assert_valid_bam "$final_bam" "Aligned BAM from Step 6"
    run_conda run -n "$PICARD_ENV" picard -Xmx${JAVA_MEM} CollectRnaSeqMetrics \
        I="$final_bam" O="$picard_out" REF_FLAT="$GENCODE_V45_REFFLAT" \
        STRAND="$picard_strand" RIBOSOMAL_INTERVALS="${HUMAN_RRNA_INTERVAL:-}"
    assert_nonempty_file "$picard_out" "Picard RNA-seq metrics"
    if [[ "$is_paired_end" == "true" ]]; then
        run_conda run -n "$PICARD_ENV" picard -Xmx${JAVA_MEM} CollectInsertSizeMetrics \
            I="$final_bam" O="$insert_out" H="$insert_pdf"
        assert_nonempty_file "$insert_out" "Picard insert-size metrics"
        assert_nonempty_file "$insert_pdf" "Picard insert-size histogram PDF"
        if command -v convert &>/dev/null; then
            run_conda run -n "$PICARD_ENV" convert -density 300 -background white -alpha remove \
                "$insert_pdf" "$insert_png" || true
        fi
    else
        log 2 "INFO" "Skipping Picard CollectInsertSizeMetrics for single-end data"
    fi
}
run_step_11() {
    local step_dir="${output_dir}/Step_11_RNA_Metrics"
    run_step "$step_dir" "false" _step_11_impl
}
_step_12_impl() {
    local step_dir="${output_dir}/Step_12_featureCounts_Quant"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.bam"
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
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.bam"
    local forward_counts="${step_dir}/${sample_name}_featureCounts_s1_forward.tsv"
    local reverse_counts="${step_dir}/${sample_name}_featureCounts_s2_reverse.tsv"
    local corrected="${step_dir}/${sample_name}_gDNA_corrected_counts.tsv"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    if [[ "$is_paired_end" == "true" ]]; then
        featureCounts -a "$TOTAL_GENE_GTF" -o "$forward_counts" \
            -p --countReadPairs -T "$thread_count" -s 1 -g gene_id -t exon -B -C "$final_bam"
        featureCounts -a "$TOTAL_GENE_GTF" -o "$reverse_counts" \
            -p --countReadPairs -T "$thread_count" -s 2 -g gene_id -t exon -B -C "$final_bam"
    else
        featureCounts -a "$TOTAL_GENE_GTF" -o "$forward_counts" \
            -T "$thread_count" -s 1 -g gene_id -t exon -B -C "$final_bam"
        featureCounts -a "$TOTAL_GENE_GTF" -o "$reverse_counts" \
            -T "$thread_count" -s 2 -g gene_id -t exon -B -C "$final_bam"
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
    local r1_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_${alignment_tag}.clean.fq.gz"
    local r2_dedup="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_${alignment_tag}.clean.fq.gz"
    local rsem_strand="$strand"
    [[ "$rsem_strand" == "unstrand" ]] && rsem_strand="none"
    assert_file_exists "$r1_dedup" "R1 processed FASTQ from Step 6"
    mkdir -p "${step_dir}/tmp"
    if [[ "$is_paired_end" == "true" ]]; then
        assert_file_exists "$r2_dedup" "R2 processed FASTQ from Step 6"
        perl "$RSEM_CALC_EXPR" --paired-end --bowtie2 --strandedness "$rsem_strand" \
            --bowtie2-k 2 -p "$thread_count" --no-bam-output --seed 12345 \
            "$r1_dedup" "$r2_dedup" "${RSEM_BOWTIE2_INDEX:-}" \
            "${step_dir}/${sample_name}_RSEM" --temporary-folder "${step_dir}/tmp"
    else
        perl "$RSEM_CALC_EXPR" --bowtie2 --strandedness "$rsem_strand" \
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
        run_python "${EVscope_PATH}/bin/Step_15_combine_total_RNA_expr_matrix.py" \
            --gene_expr "$base_matrix" --out_matrix "$combined"
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
        run_python "${EVscope_PATH}/bin/Step_15_combine_total_RNA_expr_matrix.py" \
            --gene_expr "$base_matrix" --out_matrix "$combined"
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
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"

    # Use pre-built merged SAF (7 regions, 801K features, mutually exclusive)
    local merged_saf="${STEP18_MERGED_SAF}"
    assert_file_exists "$merged_saf" "Merged meta-gene SAF from EVscope.conf"
    local merged_output="${step_dir}/${sample_name}_metagene_featureCounts.tsv"

    # featureCounts defaults: unique-mapping only (no -M), unstranded (-s 0)
    # PE: -p --countReadPairs (count fragments), -B (both ends mapped), -C (no chimera)
    log 1 "INFO" "Running featureCounts with merged meta-gene SAF (unique-mapping, unstranded)"
    if [[ "$is_paired_end" == "true" ]]; then
        featureCounts -F SAF -a "$merged_saf" -o "$merged_output" \
            -p --countReadPairs -B -C -T "$thread_count" -s 0 "$final_bam"
    else
        featureCounts -F SAF -a "$merged_saf" -o "$merged_output" \
            -T "$thread_count" -s 0 "$final_bam"
    fi

    # Split merged output into per-region files (for Step 24 and pie chart)
    declare -A region_files=()
    IFS=',' read -ra region_labels <<< "${STEP18_REGION_LABELS}"
    for region in "${region_labels[@]}"; do
        local region_output="${step_dir}/${sample_name}_HG38_${region}_noOverlap_featureCounts.tsv"
        head -2 "$merged_output" > "$region_output"
        grep "^${region}__" "$merged_output" >> "$region_output" || true
        region_files["$region"]="$region_output"
    done

    # Plot pie chart
    local -a plot_args=("--sampleName" "${sample_name}" "--output_dir" "$step_dir" "--saf_file" "${STEP18_MERGED_SAF}")
    local -A region_patterns=(
        ["5UTR"]="--input_5UTR_readcounts"
        ["exon"]="--input_exon_readcounts"
        ["3UTR"]="--input_3UTR_readcounts"
        ["intron"]="--input_intron_readcounts"
        ["promoter"]="--input_promoters_readcounts"
        ["downstream"]="--input_downstream_2Kb_readcounts"
        ["intergenic"]="--input_intergenic_readcounts"
    )
    local pattern flag
    for pattern in "${!region_patterns[@]}"; do
        flag="${region_patterns[$pattern]}"
        if [[ -n "${region_files[$pattern]+x}" && -f "${region_files[$pattern]}" ]]; then
            plot_args+=("$flag" "${region_files[$pattern]}")
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
        run_conda run -n "$KRAKEN2_ENV" kraken2 --db "${KRAKEN_DB:-}" --threads "$thread_count" \
            --report "$kraken_report" --paired --gzip-compressed "$r1_ds" "$r2_ds"
    else
        run_conda run -n "$KRAKEN2_ENV" kraken2 --db "${KRAKEN_DB:-}" --threads "$thread_count" \
            --report "$kraken_report" --gzip-compressed "$r1_ds"
    fi
    run_python "${KRAKEN_TOOLS_DIR:-}/kreport2krona.py" -r "$kraken_report" -o "$krona_input"
    run_conda run -n "$KRAKEN2_ENV" ktImportText "$krona_input" -o "$krona_html"
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
    local -a qc_args=("--output" "${step_dir}/${sample_name}_QC_summary.tsv" "--kit" "$kit")
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
            cmd_prefix="$(printf %q "$(conda_command)") run -n $(printf %q "${MULTIQC_ENV}")"
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
    done < <(find "${output_dir}/Step_01_Raw_QC" -maxdepth 1 -type f -name "*_fastqc.zip" -print0 2>/dev/null | sort -z || true)
    (( ${#raw_fqc[@]} > 0 )) && qc_args+=("--raw_fastqc_zips" "${raw_fqc[@]}")
    
    local -a trimmed_fqs=()
    while IFS= read -r -d '' f; do
        trimmed_fqs+=("$f")
    done < <(find "${output_dir}/Step_03_UMI_Adaptor_Trim" -maxdepth 1 -type f -name "${sample_name}_*_clean.fq.gz" -print0 2>/dev/null | sort -z || true)
    (( ${#trimmed_fqs[@]} > 0 )) && qc_args+=("--trimmed_fastqs" "${trimmed_fqs[@]}")

    local bbsplit_dir="${output_dir}/Step_05_Bacterial_Filter"
    local -a ecoli_fqs=()
    while IFS= read -r -d '' f; do
        ecoli_fqs+=("$f")
    done < <(find "$bbsplit_dir" -maxdepth 1 -type f -name "${sample_name}_*.fq.gz" \
        \( -iname "*escherichia*" -o -iname "*e.coli*" -o -iname "*ecoli*" -o -iname "*e_coli*" \) \
        -print0 2>/dev/null | sort -z || true)
    [[ -f "${bbsplit_dir}/step.done" ]] && qc_args+=("--ecoli_fastqs" "${ecoli_fqs[@]}")

    local -a myco_fqs=()
    while IFS= read -r -d '' f; do
        myco_fqs+=("$f")
    done < <(find "$bbsplit_dir" -maxdepth 1 -type f -name "${sample_name}_*.fq.gz" -iname "*mycoplasma*" -print0 2>/dev/null | sort -z || true)
    [[ -f "${bbsplit_dir}/step.done" ]] && qc_args+=("--myco_fastqs" "${myco_fqs[@]}")

    local ribo_dir="${output_dir}/Step_23_rRNA_Detection"
    local -a ribo_fqs=()
    while IFS= read -r -d '' f; do
        ribo_fqs+=("$f")
    done < <(find "$ribo_dir" -maxdepth 1 -type f -name "${sample_name}_rRNA_*.fq.gz" -print0 2>/dev/null | sort -z || true)
    [[ -f "${ribo_dir}/step.done" ]] && qc_args+=("--ribo_fastqs" "${ribo_fqs[@]}")
    
    local star_log_refined="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Log.final.out"
    [[ -f "$star_log_refined" ]] && qc_args+=("--STAR_log" "$star_log_refined")
    
    local star_log_initial="${output_dir}/Step_06_Alignment_Initial/${sample_name}_Log.final.out"
    [[ -f "$star_log_initial" ]] && qc_args+=("--STAR_log_initial" "$star_log_initial")

    local strand_file="${output_dir}/Step_07_Strand_Detection/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out_bam2strandness.tsv"
    [[ -f "$strand_file" ]] && qc_args+=("--bam2strand_file" "$strand_file")

    local length_strandness_file="${output_dir}/Step_07_Strand_Detection/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out_read_length_stratified_strandness.tsv"
    [[ -f "$length_strandness_file" ]] && qc_args+=("--length_strandness_file" "$length_strandness_file")
    
    local picard_insert="${output_dir}/Step_11_RNA_Metrics/${sample_name}_insert_size_metrics.tsv"
    [[ -f "$picard_insert" ]] && qc_args+=("--picard_insert_file" "$picard_insert")
    
    local picard_rnaseq="${output_dir}/Step_11_RNA_Metrics/${sample_name}_picard_metrics.tsv"
    [[ -f "$picard_rnaseq" ]] && qc_args+=("--picard_rnaseq_file" "$picard_rnaseq")
    
    local acc_motif="${output_dir}/Step_02_UMI_Analysis/${sample_name}_ACC_motif_fraction.tsv"
    [[ -f "$acc_motif" ]] && qc_args+=("--ACC_motif_fraction" "$acc_motif")
    
    local expr_matrix_combined="${output_dir}/Step_15_featureCounts_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_featureCounts.tsv"
    local expr_matrix_gene="${output_dir}/Step_15_featureCounts_Expression/${sample_name}_Gene_readcounts_normalized_expression_matrix_featureCounts.tsv"
    if [[ -f "$expr_matrix_combined" ]]; then
        qc_args+=("--expression_matrix" "$expr_matrix_combined")
    elif [[ -f "$expr_matrix_gene" ]]; then
        qc_args+=("--expression_matrix" "$expr_matrix_gene")
    fi
    
    local kraken_report="${output_dir}/Step_19_Taxonomy/${sample_name}_report.tsv"
    [[ -f "$kraken_report" ]] && qc_args+=("--kraken_report" "$kraken_report")

    local step18_dir="${output_dir}/Step_18_Genomic_Regions"
    local fc_3utr fc_5utr fc_downstream fc_exon fc_intergenic fc_intron fc_promoter
    fc_3utr="$(find_first_sample_file "$step18_dir" "*3UTR*featureCounts.tsv")"
    fc_5utr="$(find_first_sample_file "$step18_dir" "*5UTR*featureCounts.tsv")"
    fc_downstream="$(find_first_sample_file "$step18_dir" "*downstream*featureCounts.tsv")"
    fc_exon="$(find_first_sample_file "$step18_dir" "*exon*featureCounts.tsv")"
    fc_intergenic="$(find_first_sample_file "$step18_dir" "*intergenic*featureCounts.tsv")"
    fc_intron="$(find_first_sample_file "$step18_dir" "*intron*featureCounts.tsv")"
    fc_promoter="$(find_first_sample_file "$step18_dir" "*promoter*featureCounts.tsv")"
    [[ -f "$fc_3utr" ]] && qc_args+=("--featureCounts_3UTR" "$fc_3utr")
    [[ -f "$fc_5utr" ]] && qc_args+=("--featureCounts_5UTR" "$fc_5utr")
    [[ -f "$fc_downstream" ]] && qc_args+=("--featureCounts_downstream_2kb" "$fc_downstream")
    [[ -f "$fc_exon" ]] && qc_args+=("--featureCounts_exon" "$fc_exon")
    [[ -f "$fc_intergenic" ]] && qc_args+=("--featureCounts_intergenic" "$fc_intergenic")
    [[ -f "$fc_intron" ]] && qc_args+=("--featureCounts_intron" "$fc_intron")
    [[ -f "$fc_promoter" ]] && qc_args+=("--featureCounts_promoter_1500_500bp" "$fc_promoter")
    
    local qc_script="${EVscope_PATH}/bin/Step_24_generate_QC_matrix.py"
    local qc_output="${step_dir}/${sample_name}_QC_summary.tsv"
    rm -f "$qc_output"
    
    if [[ -f "$qc_script" ]]; then
        run_python "$qc_script" "${qc_args[@]}"
    else
        log 5 "FATAL" "QC summary script not found: ${qc_script}"
        exit 1
    fi
    
    if [[ ! -s "$qc_output" ]]; then
        log 5 "FATAL" "QC summary file was not generated or is empty: ${qc_output}"
        exit 1
    fi
}
run_step_24() {
    local step_dir="${output_dir}/Step_24_MultiQC_Summary"
    run_step "$step_dir" "true" _step_24_impl
}
_step_25_emapper_impl() {
    local step_dir="${output_dir}/Step_25_EMapper_BigWig_Quantification/EMapper_output"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_${alignment_tag}_Aligned.sortedByCoord.out.bam"
    assert_file_exists "$final_bam" "Aligned BAM from Step 6"
    mkdir -p "$step_dir"
    local -a polya_args=()
    if [[ "$kit" == "smarter_total_rna" ]]; then
        polya_args=(--polyA_tail_detection)
    fi
    run_python "${EVscope_PATH}/bin/Step_25_EMapper.py" \
        --num_threads "$thread_count" --input_bam "$final_bam" \
        --prefix "${sample_name}" --output_dir "$step_dir" \
        --strandness "$strand" --reference_fasta "$HUMAN_GENOME_FASTA" \
        --gtf "$TOTAL_GENE_GTF" --positional_assignment EM \
        --gene_disambiguation EM --no-cleanup "${polya_args[@]}"
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
_build_json_array_from_bash_array() {
    local -n arr=$1
    local result="["
    local first=true item
    for item in "${arr[@]}"; do
        if [[ "$first" == "true" ]]; then
            first=false
        else
            result+=","
        fi
        result+="\"$(json_escape "$item")\""
    done
    result+="]"
    echo "$result"
}
_build_json_array_from_csv() {
    local csv="${1:-}"
    csv="${csv#[}"
    csv="${csv%]}"
    csv="${csv//\"/}"
    local result="["
    local first=true item
    local IFS=','
    for item in $csv; do
        item="${item#${item%%[![:space:]]*}}"
        item="${item%${item##*[![:space:]]}}"
        [[ -z "$item" ]] && continue
        if [[ "$first" == "true" ]]; then
            first=false
        else
            result+=","
        fi
        result+="\"$(json_escape "$item")\""
    done
    result+="]"
    echo "$result"
}
_step_26_impl() {
    local step_dir="${output_dir}/Step_26_BigWig_Density_Plot"
    local bigwig="${output_dir}/Step_25_EMapper_BigWig_Quantification/EMapper_output/${sample_name}_unstranded.bw"
    assert_file_exists "$bigwig" "Unstranded BigWig from Step 25"

    local bed_rna bed_meta bed_labels_rna bed_labels_meta
    bed_rna="$(_build_json_array_from_bash_array STEP26_RNATYPE_BEDS)"
    bed_labels_rna="$(_build_json_array_from_csv "${STEP26_RNATYPE_LABELS:-}")"
    bed_meta="$(_build_json_array_from_bash_array STEP26_METAGENE_BEDS)"
    bed_labels_meta="$(_build_json_array_from_csv "${STEP26_METAGENE_LABELS:-}")"

    local rna_bed_count=0 meta_bed_count=0 rna_label_count meta_label_count
    declare -p STEP26_RNATYPE_BEDS &>/dev/null && rna_bed_count="${#STEP26_RNATYPE_BEDS[@]}"
    declare -p STEP26_METAGENE_BEDS &>/dev/null && meta_bed_count="${#STEP26_METAGENE_BEDS[@]}"
    rna_label_count="$(count_csv_values "${STEP26_RNATYPE_LABELS:-}")"
    meta_label_count="$(count_csv_values "${STEP26_METAGENE_LABELS:-}")"
    if (( rna_bed_count == 0 || meta_bed_count == 0 || rna_label_count == 0 || meta_label_count == 0 )); then
        log 5 "FATAL" "Step 26 requires non-empty BED arrays and labels for both RNA type and meta-gene plots."
        exit 1
    fi
    if (( rna_bed_count != rna_label_count )); then
        log 5 "FATAL" "STEP26_RNATYPE_BEDS count (${rna_bed_count}) != STEP26_RNATYPE_LABELS count (${rna_label_count})"
        exit 1
    fi
    if (( meta_bed_count != meta_label_count )); then
        log 5 "FATAL" "STEP26_METAGENE_BEDS count (${meta_bed_count}) != STEP26_METAGENE_LABELS count (${meta_label_count})"
        exit 1
    fi

    if [[ "$force_mode" == "true" ]]; then
        log 2 "INFO" "Force mode: clearing cached Step 26 matrix/profile outputs before recomputing."
        rm -rf -- "${step_dir}/RNA_types" "${step_dir}/meta_gene"
    fi

    bash "${EVscope_PATH}/bin/Step_26_density_plot_over_RNA_types.sh" \
        --input_bw_file "$bigwig" --input_bed_files "$bed_rna" \
        --input_bed_labels "$bed_labels_rna" --output_dir "${step_dir}/RNA_types" \
        --threads "$thread_count" --random_tested_row_num_per_bed 100000
    bash "${EVscope_PATH}/bin/Step_26_density_plot_over_meta_gene.sh" \
        --input_bw_file "$bigwig" --input_bed_files "$bed_meta" \
        --input_bed_labels "$bed_labels_meta" --output_dir "${step_dir}/meta_gene" \
        --blackListFileName "${ENCODE_BLACKLIST_BED:-}" --threads "$thread_count" --random_tested_row_num_per_bed 100000
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
    local summary_script="${abs_bin_dir}/Step_27_generate_summary_report.py"
    local summary_file="${sample_name}_summary_report.html"
    local local_rmd="${abs_step_dir}/${sample_name}_report.Rmd"
    assert_file_exists "$rmd_input" "R Markdown template"
    assert_file_exists "$summary_script" "concise summary report generator"
    mkdir -p "$abs_step_dir"
    cp "$rmd_input" "$local_rmd"

    local qc_summary="${abs_output_dir}/Step_24_MultiQC_Summary/${sample_name}_QC_summary.tsv"
    local expression_matrix="${abs_output_dir}/Step_15_featureCounts_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_featureCounts.tsv"
    if [[ ! -f "$expression_matrix" ]]; then
        expression_matrix="${abs_output_dir}/Step_15_featureCounts_Expression/${sample_name}_Gene_readcounts_normalized_expression_matrix_featureCounts.tsv"
    fi
    local -a summary_args=(
        --sample-name "$sample_name"
        --kit "$kit"
        --strand "$strand"
        --paired-end "$is_paired_end"
        --output "${abs_step_dir}/${summary_file}"
    )
    [[ -f "$qc_summary" ]] && summary_args+=(--qc-summary "$qc_summary")
    [[ -f "$expression_matrix" ]] && summary_args+=(--expression-matrix "$expression_matrix")
    run_python "$summary_script" "${summary_args[@]}"
    assert_nonempty_file "${abs_step_dir}/${summary_file}" "standalone concise summary report"

    export EVscope_OUTPUT_DIR="$abs_output_dir"
    export EVscope_FIGURES_DIR="$abs_figures_dir"
    export EVscope_KIT="$kit"
    export EVscope_RMD_INPUT="$local_rmd"
    export EVscope_REPORT_FILE="$report_file"
    export EVscope_REPORT_DIR="$abs_step_dir"
    local -a rscript_cmd=()
    build_rscript_command rscript_cmd || exit 1
    cd "$abs_step_dir" || exit 1
    "${rscript_cmd[@]}" -e 'rmarkdown::render(input = Sys.getenv("EVscope_RMD_INPUT"), output_file = Sys.getenv("EVscope_REPORT_FILE"), output_dir = Sys.getenv("EVscope_REPORT_DIR"), intermediates_dir = Sys.getenv("EVscope_REPORT_DIR"))'
    assert_nonempty_file "${abs_step_dir}/${report_file}" "final HTML report"
}
run_step_27() {
    local step_dir="${output_dir}/Step_27_HTML_Report"
    run_step "$step_dir" "false" _step_27_impl
}
# ==============================================================================
# SECTION: MAIN EXECUTION LOGIC
# ==============================================================================
main() {
    local early_tmp_dir="${TMPDIR:-/tmp}"
    mkdir -p "$early_tmp_dir"
    early_log_file="$(mktemp "${early_tmp_dir%/}/evscope_early_XXXXXX.log")"
    CLEANUP_FILES+=("$early_log_file")
    check_bash_version

    config_file="${SCRIPT_DIR}/EVscope.conf"
    local config_file_input="$config_file"
    local -a input_fastqs=()
    local dry_run="false"
    local strand_was_set="false"

    while (( $# > 0 )); do
        case "$1" in
            --sample_name)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--sample_name requires a value"; exit 1; }
                sample_name="$2"; shift 2 ;;
            --threads)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--threads requires a value"; exit 1; }
                [[ ! "$2" =~ ^[0-9]+$ ]] && { log 5 "FATAL" "--threads must be a positive integer"; exit 1; }
                (( $2 < 1 )) && { log 5 "FATAL" "--threads must be >= 1"; exit 1; }
                thread_count="$2"; shift 2 ;;
            --kit)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--kit requires a value"; exit 1; }
                local normalized_kit=""
                if ! normalized_kit="$(normalize_kit "$2")"; then
                    log 5 "FATAL" "--kit must be 'smarter_total_rna' or 'smarter_smrna' (received: '$2')"
                    exit 1
                fi
                kit="$normalized_kit"; shift 2 ;;
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
                strand="$2"; strand_was_set="true"; shift 2 ;;
            --config)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--config requires a value"; exit 1; }
                config_file="$2"; config_file_input="$2"; shift 2 ;;
            --resume)
                resume_mode="true"; shift ;;
            --force)
                force_mode="true"; shift ;;
            -V|--verbosity)
                [[ -z "${2:-}" ]] && { log 5 "FATAL" "--verbosity requires a value"; exit 1; }
                [[ ! "$2" =~ ^[1-5]$ ]] && { log 5 "FATAL" "--verbosity must be 1-5"; exit 1; }
                verbosity="$2"; shift 2 ;;
            --input_fastqs)
                shift
                while (( $# > 0 )) && [[ "$1" != -* ]]; do
                    if [[ "$1" == *,* ]]; then
                        local -a split_fastqs=()
                        local fq_piece
                        IFS=',' read -r -a split_fastqs <<< "$1"
                        for fq_piece in "${split_fastqs[@]}"; do
                            fq_piece="${fq_piece#${fq_piece%%[![:space:]]*}}"
                            fq_piece="${fq_piece%${fq_piece##*[![:space:]]}}"
                            [[ -n "$fq_piece" ]] && input_fastqs+=("$fq_piece")
                        done
                    else
                        input_fastqs+=("$1")
                    fi
                    shift
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

    apply_kit_profile "$strand_was_set"
    if [[ "$kit" == "smarter_smrna" && "$strand_was_set" == "true" && "$strand" != "forward" ]]; then
        log 3 "WARN" "SMARTer smRNA Read 1 is sense; explicit --strand ${strand} will override the recommended forward setting."
    fi
    if [[ "$gDNA_correction" == "yes" && "$strand" == "unstrand" ]]; then
        log 5 "FATAL" "Logic Error: --gDNA_correction yes requires --strand forward/reverse."
        exit 1
    fi
    if [[ -z "$sample_name" ]]; then
        log 5 "FATAL" "--sample_name is required"
        print_help
    fi
    if [[ ! "$sample_name" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
        log 5 "FATAL" "Invalid --sample_name '${sample_name}'. Use only letters, numbers, dot, underscore, and hyphen; the first character must be alphanumeric."
        exit 1
    fi
    if (( ${#input_fastqs[@]} < 1 || ${#input_fastqs[@]} > 2 )); then
        log 5 "FATAL" "--input_fastqs requires exactly 1 (SE) or 2 (PE) files"
        log 5 "FATAL" "Provided: ${#input_fastqs[@]} files"
        exit 1
    fi
    output_dir="${sample_name}_EVscope_output"
    if [[ -d "$output_dir" ]]; then
        if [[ -n "$(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" && "$resume_mode" != "true" && "$force_mode" != "true" ]]; then
            log 5 "FATAL" "Output directory already exists and is not empty: ${output_dir}. Use --resume to continue with metadata checks or --force to regenerate selected steps."
            exit 1
        fi
    elif ! mkdir -p "$output_dir"; then
        log 5 "FATAL" "Cannot create output directory: ${output_dir}"
        exit 1
    fi
    touch "${output_dir}/EVscope_pipeline.log"
    if [[ -f "$early_log_file" && -s "$early_log_file" ]]; then
        cat "$early_log_file" >> "${output_dir}/EVscope_pipeline.log"
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
        validate_fastq_pair_order "$fastq_read1" "$fastq_read2"
    fi

    local run_steps_list skip_steps_list final_steps_list step
    run_steps_list="$(parse_steps "$run_steps")"
    skip_steps_list=""
    if [[ -n "$skip_steps" ]]; then
        skip_steps_list="$(parse_steps "$skip_steps")"
    fi
    if [[ -n "$skip_steps_list" ]]; then
        final_steps_list="$(echo "$run_steps_list" | grep -vFx -f <(echo "$skip_steps_list") || true)"
    else
        final_steps_list="$run_steps_list"
    fi
    STEPS_TO_RUN=()
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
        STEPS_TO_RUN+=("$step")
    done <<< "$final_steps_list"
    if (( ${#STEPS_TO_RUN[@]} == 0 )); then
        log 5 "FATAL" "No steps to run after applying filters"
        exit 1
    fi

    config_file="$(get_absolute_path "$config_file")"
    if [[ -z "$config_file" || ! -f "$config_file" ]]; then
        log 5 "FATAL" "Configuration file not found: input='${config_file_input}' resolved='${config_file:-<empty>}'"
        exit 1
    fi
    # shellcheck source=/dev/null
    source "$config_file"
    log 2 "INFO" "Loaded configuration from: ${config_file}"
    set_star_config_defaults
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

    local abs_output_dir
    abs_output_dir="$(realpath -m "$output_dir" 2>/dev/null || printf '%s' "$output_dir")"
    log 2 "INFO" "============================================================"
    log 2 "INFO" "EVscope Pipeline v${VERSION}"
    log 2 "INFO" "============================================================"
    log 2 "INFO" "Sample name:      ${sample_name}"
    log 2 "INFO" "Library kit:      ${kit}"
    log 2 "INFO" "Output directory: ${output_dir} (${abs_output_dir})"
    log 2 "INFO" "Paired-end mode:  ${is_paired_end}"
    log 2 "INFO" "Thread count:     ${thread_count}"
    log 2 "INFO" "circRNA tool:     ${circ_tool}"
    log 2 "INFO" "Strand:           ${strand}"
    log 2 "INFO" "gDNA correction:  ${gDNA_correction}"
    log 2 "INFO" "Resume mode:      ${resume_mode}"
    log 2 "INFO" "Force mode:       ${force_mode}"
    log 2 "INFO" "Steps to execute: ${STEPS_TO_RUN[*]}"
    log 2 "INFO" "============================================================"

    if [[ "$dry_run" == "true" ]]; then
        log 2 "INFO" "DRY-RUN MODE: Execution plan validated successfully"
        print_pipeline_steps
        log 2 "INFO" "To run the pipeline, remove --dry-run flag"
        exit 0
    fi

    local start_time
    start_time="$(date +%s)"
    print_pipeline_steps
    for step in "${STEPS_TO_RUN[@]}"; do
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
    log 2 "INFO" "Output:     ${output_dir} (${abs_output_dir})"
    log 2 "INFO" "============================================================"
    return 0
}

# ==============================================================================
# SCRIPT ENTRY POINT
# ==============================================================================
main "$@"
