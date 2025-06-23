#!/usr/bin/env bash
#
# ==============================================================================
# EVscope.sh: A Comprehensive and Modular RNA-seq Analysis Pipeline
#
# Version: 2.4.0
#
# Description:
# This script provides a fully modular, user-driven pipeline for RNA-seq data
# analysis. It starts from raw FASTQ files and performs a series of steps
# including QC, trimming, alignment, quantification, and various downstream
# analyses like circular RNA detection, contamination screening, and tissue
# deconvolution. This version features a robust two-pass STAR alignment
# strategy critical for accurate circRNA detection.
#
# Author: [Your Name/Organization]
# Last Modified: [Date]
# ==============================================================================

# --- Strict Mode & Safety ---
# Exit immediately if a command exits with a non-zero status.
# Treat unset variables as an error when substituting.
# Pipelines (e.g., cmd1 | cmd2) return the exit status of the last command
# that failed, or zero if all commands in the pipeline succeed.
set -eo pipefail

# --- Global Declarations ---
# Using readonly for script constants to prevent accidental modification.
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
readonly VERSION="2.4.0"
# Set pipeline base directory early, making it available for the config file.
readonly PIPELINE_BASE_DIR=$SCRIPT_DIR

# ==============================================================================
# --- CORE UTILITIES: LOGGING, HELP, AND DEPENDENCY CHECKS ---
# ==============================================================================

# Default verbosity level (2: INFO). Can be changed by --verbosity option.
# 1=DEBUG, 2=INFO, 3=WARN, 4=ERROR, 5=FATAL
verbosity=2

# ANSI Color Codes for terminal output. They are only defined if stderr is a tty.
if [ -t 2 ]; then
    readonly C_RESET='\033[0m'
    readonly C_RED='\033[0;31m'
    readonly C_GREEN='\033[0;32m'
    readonly C_YELLOW='\033[0;33m'
    readonly C_BLUE='\033[0;34m'
    readonly C_BOLD='\033[1m'
fi

# Professional logging function.
# Usage: log <LEVEL_NUM> <LEVEL_NAME> "Message"
log() {
    local level_num=$1
    local level_name=$2
    shift 2
    local message="$*"

    # Only log messages at or above the current verbosity level.
    (( level_num < verbosity )) && return 0

    local color_prefix=""
    case $level_name in
        DEBUG) color_prefix=$C_BLUE ;;
        INFO)  color_prefix=$C_GREEN ;;
        WARN)  color_prefix=$C_YELLOW ;;
        ERROR|FATAL) color_prefix=$C_RED ;;
    esac

    local timestamp
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    # Print formatted message to stderr.
    echo -e "${color_prefix}[$timestamp] ${level_name}:${C_RESET} ${message}" >&2

    # Append plain message to the main log file if the output directory is set.
    if [[ -n "${output_dir:-}" && -d "$output_dir" ]]; then
        echo "[$timestamp] $level_name: $message" >> "${output_dir}/EVscope_pipeline.log"
    fi
}

# Displays the help message and exits.
print_help() {
    cat << EOF
Usage: bash $SCRIPT_NAME [options] --sample_name <name> --input_fastqs <files>

A fully modular RNA-seq analysis pipeline.

Required Options:
  --sample_name <name>    A unique name for the sample, used for output files and directories.
  --input_fastqs <files>  Input FASTQ file(s). For paired-end data, provide a single
                          string with paths separated by a comma (e.g., "r1.fq.gz,r2.fq.gz").

Optional Options:
  --threads <int>           Number of CPU threads to use (default: 1).
  --run_steps <list>        A comma-separated list of steps or ranges to run (e.g., "1,3,5-8", "all").
                            No dependency checks are performed; run prerequisite steps first.
  --skip_steps <list>       A comma-separated list of steps or ranges to skip from the --run_steps list.
  --circ_tool <tool>        circRNA detection tool to use: 'CIRCexplorer2', 'CIRI2', or 'both' (default: both).
  --read_count_mode <mode>  Read counting strategy: 'uniq' (featureCounts) or 'multi' (RSEM) (default: uniq).
  --gDNA_correction <yes|no> Apply genomic DNA contamination correction (default: no).
  --strandedness <strand>   Library strandedness: 'reverse', 'forward', or 'unstrand' (default: reverse).
  --config <path>           Path to a custom configuration file (default: EVscope.conf in script dir).
  -V, --verbosity <level>   Set logging verbosity: 1-DEBUG, 2-INFO (default), 3-WARN, 4-ERROR.
  -h, --help                Display this help message and exit.
  -v, --version             Display the pipeline version and exit.
EOF
    exit 0
}

# Displays the script version and exits.
print_version() {
    echo "$SCRIPT_NAME Version: $VERSION"
    exit 0
}

# Checks for the presence of essential command-line tools.
check_dependencies() {
    log 2 "INFO" "Checking for required software dependencies..."
    local missing_deps=0
    local dependencies=(
        "fastqc" "umi_tools" "cutadapt" "STAR" "samtools" "seqtk" "bwa" "perl" "python" "Rscript"
        "conda" "featureCounts" "ribodetector_cpu" "getopt" "awk" "tee" "date" "mkdir" "touch" "convert"
    )
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            log 4 "ERROR" "Dependency not found in PATH: '$cmd'."
            missing_deps=$((missing_deps + 1))
        fi
    done

    if [ "$missing_deps" -gt 0 ]; then
        log 5 "FATAL" "Please install all missing dependencies or ensure they are in your system's PATH."
        exit 1
    fi
    log 2 "INFO" "All basic software dependencies are satisfied."
}

# ==============================================================================
# --- PIPELINE STEP DEFINITIONS & WORKFLOW LOGIC ---
# ==============================================================================

# Prints a list of all available pipeline steps.
print_pipeline_steps() {
    log 2 "INFO" "EVscope Pipeline Steps (Version: $VERSION)"
    log 2 "INFO" "========================================"
    log 2 "INFO" "The pipeline consists of 27 steps. Use --run_steps and --skip_steps to customize execution."
    log 2 "INFO" "Step  | Description"
    log 2 "INFO" "------|----------------------------------------------------------------"
    log 2 "INFO" "1     | Raw FASTQ quality control using FastQC"
    log 2 "INFO" "2     | UMI motif analysis and ratio calculation"
    log 2 "INFO" "3     | UMI labeling and adapter/quality trimming"
    log 2 "INFO" "4     | Quality control of trimmed FASTQs"
    log 2 "INFO" "5     | Bacterial contamination detection (E. coli, Mycoplasma)"
    log 2 "INFO" "6     | STAR Two-Pass Alignment (Initial + Refined for circRNA)"
    log 2 "INFO" "7     | Library strandedness detection"
    log 2 "INFO" "8     | CIRCexplorer2 circRNA detection"
    log 2 "INFO" "9     | CIRI2 circRNA detection"
    log 2 "INFO" "10    | Merge CIRCexplorer2 and CIRI2 circRNA results"
    log 2 "INFO" "11    | RNA-seq metrics collection (Picard)"
    log 2 "INFO" "12    | featureCounts quantification (unique-mapping mode)"
    log 2 "INFO" "13    | gDNA-corrected featureCounts quantification"
    log 2 "INFO" "14    | RSEM quantification (multi-mapping mode)"
    log 2 "INFO" "15    | featureCounts-based expression matrix and RNA distribution plots"
    log 2 "INFO" "16    | gDNA-corrected expression matrix and RNA distribution plots"
    log 2 "INFO" "17    | RSEM-based expression matrix and RNA distribution plots"
    log 2 "INFO" "18    | Genomic region read mapping analysis (3'UTR, 5'UTR, etc.)"
    log 2 "INFO" "19    | Taxonomic classification using Kraken2"
    log 2 "INFO" "20    | Tissue deconvolution for featureCounts results"
    log 2 "INFO" "21    | Tissue deconvolution for gDNA-corrected results"
    log 2 "INFO" "22    | Tissue deconvolution for RSEM results"
    log 2 "INFO" "23    | rRNA detection using ribodetector"
    log 2 "INFO" "24    | Comprehensive QC summary generation"
    log 2 "INFO" "25    | Coverage analysis and BigWig generation (EMapper)"
    log 2 "INFO" "26    | Coverage density plots (RNA types and meta-gene regions)"
    log 2 "INFO" "27    | Final HTML report generation"
    log 2 "INFO" "========================================"
}

# Parses a user-provided string of steps (e.g., "1,3,5-8") into a sorted, unique list.
parse_steps() {
    local input; input=$(echo "$1" | tr -d '[] ')
    log 1 "DEBUG" "Parsing run_steps string: '$input'"
    if [[ "$input" == "all" ]]; then
        seq 1 27
        return
    fi

    local steps=()
    IFS=',' read -ra parts <<< "$input"
    for part in "${parts[@]}"; do
        if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            local start=${BASH_REMATCH[1]}
            local end=${BASH_REMATCH[2]}
            for ((i=start; i<=end; i++)); do steps+=("$i"); done
        elif [[ "$part" =~ ^[0-9]+$ ]]; then
            steps+=("$part")
        else
            log 5 "FATAL" "Invalid format in --run_steps or --skip_steps: '$part'"
            exit 1
        fi
    done
    # Sort numerically and remove duplicates
    printf "%s\n" "${steps[@]}" | sort -n | uniq
}

# Generic step execution wrapper. Creates a directory, checks for completion,
# runs the command, logs output, and handles errors.
# Usage: run_step <directory> <command_string> [ignore_errors_boolean]
run_step() {
    local dir="$1"
    local cmd="$2"
    local ignore_err="${3:-false}"

    if ! mkdir -p "$dir"; then
        log 5 "FATAL" "Failed to create step directory: $dir"
        exit 1
    fi

    # Check for a "done" file to allow skipping of completed steps.
    if [[ -f "$dir/step.done" ]]; then
        log 2 "INFO" "Step '$(basename "$dir")' already completed. Skipping."
        return 0
    fi

    log 2 "INFO" "==> Running Step: $(basename "$dir") <=="
    local start_time; start_time=$(date +%s)
    local log_file="$dir/step.stderr.log"

    # Execute the command, redirecting stderr to a log file.
    # The `eval` is necessary to correctly interpret the command string with quotes and variables.
    if (eval "$cmd") 2> "$log_file"; then
        touch "$dir/step.done"
        local end_time; end_time=$(date +%s)
        log 2 "INFO" "Finished step: $(basename "$dir") in $((end_time - start_time)) seconds."
    else
        local exit_status=$?
        log 4 "ERROR" "Step FAILED: $(basename "$dir") (Exit code: $exit_status)"
        log 4 "ERROR" "Check stderr log for details: $log_file"
        # Print error log to the console and main pipeline log.
        if [[ -s "$log_file" ]]; then
            cat "$log_file" >&2
            cat "$log_file" >> "${output_dir}/EVscope_pipeline.log"
        fi

        if [[ "$ignore_err" == "true" ]]; then
            log 3 "WARN" "Ignoring error as requested and continuing..."
            # Create a "done" file even on failure to prevent re-running.
            touch "$dir/step.done"
        else
            log 5 "FATAL" "Pipeline stopped due to a critical error."
            exit $exit_status
        fi
    fi
}

# ==============================================================================
# --- INDIVIDUAL STEP FUNCTIONS (1-27) ---
# ==============================================================================
# Each function defines the specific command for a pipeline step and calls the
# `run_step` wrapper for execution.

run_step_1() {
    local step_dir="${output_dir}/Step_01_Raw_QC"
    local cmd="fastqc -o \"$step_dir\" -t \"$thread_count\" ${input_fastqs[*]}"
    run_step "$step_dir" "$cmd"
}

run_step_2() {
    local step_dir="${output_dir}/Step_02_UMI_Analysis"
    local input_fq="${fastq_read2:-${fastq_read1}}" # Use R2 if paired, else R1
    local cmd="
        # 1. Plot UMI motif from the first 1M reads of R2.
        python \"${EVscope_PATH}/bin/Step_02_plot_fastq2UMI_motif.py\" \\
            -head \"${sample_name}\" \\
            -fq \"${input_fq}\" \\
            -n 14 -r 1000000 \\
            -o \"$step_dir\" && \\
        # 2. Calculate the fraction of 'ACC' motif at specific positions.
        python \"${EVscope_PATH}/bin/Step_02_calculate_ACC_motif_fraction.py\" \\
            --input_fastq \"${input_fq}\" \\
            --positions 9-11 \\
            --output_tsv \"${step_dir}/${sample_name}_ACC_motif_fraction.tsv\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_3() {
    local step_dir="${output_dir}/Step_03_UMI_Adaptor_Trim"
    # Define intermediate filenames for clarity
    local r1_umi_tools="${step_dir}/${sample_name}_R1_umi_tools.fq.gz"
    local r2_umi_tools="${step_dir}/${sample_name}_R2_umi_tools.fq.gz"
    local r1_adapter_trimmed="${step_dir}/${sample_name}_R1_adapter_trimmed.fq.gz"
    local r2_adapter_trimmed="${step_dir}/${sample_name}_R2_adapter_trimmed.fq.gz"
    local r1_adapter_umi_trimmed="${step_dir}/${sample_name}_R1_adapter_UMI_trimmed.fq.gz"
    local r2_adapter_umi_trimmed="${step_dir}/${sample_name}_R2_adapter_UMI_trimmed.fq.gz"
    local r1_clean="${step_dir}/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${step_dir}/${sample_name}_R2_clean.fq.gz"

    local cmd="
        # 1. Extract UMI from R2 and append it to read headers.
        umi_tools extract --bc-pattern='NNNNNNNNNNNNNN' \\
            --stdin=\"${fastq_read2:-${fastq_read1}}\" --stdout=\"$r2_umi_tools\" \\
            ${is_paired_end:+--read2-in=\"${fastq_read1}\" --read2-out=\"$r1_umi_tools\"} \\
            --log=\"${step_dir}/UMI_extract.log\" --umi-separator='_' && \\

        # 2. Trim standard adapters using cutadapt.
        cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC --overlap 3 --minimum-length 10 -j \"$thread_count\" \\
            -o \"$r1_adapter_trimmed\" ${is_paired_end:+-p \"$r2_adapter_trimmed\"} \\
            \"$r1_umi_tools\" ${is_paired_end:+\"$r2_umi_tools\"} && \\

        # 3. Custom trimming of read-through UMIs from R1.
        python \"${EVscope_PATH}/bin/Step_03_UMIAdapterTrimR1.py\" \\
            --input_R1_fq \"$r1_adapter_trimmed\" ${is_paired_end:+--input_R2_fq \"$r2_adapter_trimmed\"} \\
            --output_R1_fq \"$r1_adapter_umi_trimmed\" ${is_paired_end:+--output_R2_fq \"$r2_adapter_umi_trimmed\"} \\
            --output_tsv \"${step_dir}/${sample_name}_R1_readthrough_UMI_trimming.log\" \\
            --min-overlap 3 --min-length 10 --chunk-size 100000 --error-rate 0.1 && \\

        # 4. Final quality trimming (phred score >= 20).
        cutadapt -q 20 --minimum-length 10 -j \"$thread_count\" \\
            -o \"$r1_clean\" ${is_paired_end:+-p \"$r2_clean\"} \\
            \"$r1_adapter_umi_trimmed\" ${is_paired_end:+\"$r2_adapter_umi_trimmed\"} && \\
        
        # 5. Plot read length distribution of the final clean FASTQs.
        python \"${EVscope_PATH}/bin/Step_03_plot_fastq_read_length_dist.py\" \\
            --input_fastqs \"$r1_clean\" ${is_paired_end:+\"$r2_clean\"} \\
            --output_pdf \"${step_dir}/${sample_name}_read_length_distribution.pdf\" \\
            --output_png \"${step_dir}/${sample_name}_read_length_distribution.png\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_4() {
    local step_dir="${output_dir}/Step_04_Trimmed_QC"
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local cmd="fastqc -o \"$step_dir\" -t \"$thread_count\" \"$r1_clean\" ${is_paired_end:+\"$r2_clean\"}"
    run_step "$step_dir" "$cmd"
}

run_step_5() {
    local step_dir="${output_dir}/Step_05_Bacterial_Filter"
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local cmd="
        bash \"$BBSPLIT_SCRIPT\" build=1 threads=\"$thread_count\" \\
            in1=\"$r1_clean\" ${is_paired_end:+in2=\"$r2_clean\"} \\
            ref=\"$ECOLI_GENOME_FASTA,$MYCOPLASMA_GENOME_FASTA\" \\
            basename=\"${step_dir}/${sample_name}_%_R#.fq.gz\" \\
            ambiguous=best path=\"$step_dir\"
    "
    run_step "$step_dir" "$cmd"
}

# --- Step 6: Two-Pass Alignment is broken into two sub-functions for clarity ---
_run_step_6_initial() {
    local step_dir="${output_dir}/Step_06_Alignment_Initial"
    # Define file paths
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local initial_bam="${step_dir}/${sample_name}_Aligned.sortedByCoord.out.bam"
    local dedup_bam="${step_dir}/${sample_name}_Aligned.sortedByCoord_umi_dedup.out.bam"
    local readnames_txt="${step_dir}/${sample_name}_readnames.txt"
    local r1_dedup_fq="${step_dir}/${sample_name}_R1_umi_dedup.clean.fq.gz"
    local r2_dedup_fq="${step_dir}/${sample_name}_R2_umi_dedup.clean.fq.gz"

    local cmd="
        ulimit -n 65535;

        # 1. First-pass STAR alignment.
        STAR --genomeDir \"$STAR_INDEX\" \\
             --readFilesIn \"$r1_clean\" ${is_paired_end:+\"$r2_clean\"} \\
             --outFileNamePrefix \"${step_dir}/${sample_name}_\" \\
             --runThreadN \"$thread_count\" --twopassMode Basic --runMode alignReads \\
             --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \\
             --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \\
             --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN \"$thread_count\";

        # 2. Index the initial alignment BAM.
        samtools index -@ \"$thread_count\" \"$initial_bam\";

        # 3. Deduplicate reads based on UMI.
        umi_tools dedup -I \"$initial_bam\" -S \"$dedup_bam\" \\
            --log=\"${step_dir}/${sample_name}_umi_dedup.log\" \\
            --extract-umi-method=read_id ${is_paired_end:+--paired};

        # 4. Index the UMI-deduplicated BAM.
        samtools index -@ \"$thread_count\" \"$dedup_bam\";

        # 5. Extract read names from the deduplicated BAM for filtering FASTQs.
        samtools view -@ \"$thread_count\" \"$dedup_bam\" | cut -f1 | sort | uniq > \"$readnames_txt\";

        # 6. Create new FASTQ files containing only the deduplicated reads.
        if [ \"${is_paired_end}\" = true ]; then
            seqtk subseq \"$r1_clean\" \"$readnames_txt\" | gzip > \"$r1_dedup_fq\";
            seqtk subseq \"$r2_clean\" \"$readnames_txt\" | gzip > \"$r2_dedup_fq\";
        else
            seqtk subseq \"$r1_clean\" \"$readnames_txt\" | gzip > \"$r1_dedup_fq\";
        fi
    "
    # This sub-step is part of the larger step 6, so it doesn't create its own .done file.
    # The main `run_step_6` handles that.
    ( eval "$cmd" )
}

_run_step_6_refined() {
    local step_dir="${output_dir}/Step_06_Alignment_Refined"
    # Define file paths
    local r1_dedup_fq="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_umi_dedup.clean.fq.gz"
    local r2_dedup_fq="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_umi_dedup.clean.fq.gz"
    local final_bam="${step_dir}/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"

    local cmd="
        ulimit -n 65535;

        # 1. Second-pass STAR alignment using UMI-deduplicated FASTQs.
        STAR --genomeDir \"$STAR_INDEX\" \\
             --readFilesIn \"$r1_dedup_fq\" ${is_paired_end:+\"$r2_dedup_fq\"} \\
             --outFileNamePrefix \"${step_dir}/${sample_name}_STAR_umi_dedup_\" \\
             --runThreadN \"$thread_count\" --twopassMode Basic --runMode alignReads \\
             --quantMode GeneCounts --readFilesCommand zcat --outFilterMultimapNmax 100 \\
             --winAnchorMultimapNmax 100 --outSAMtype BAM SortedByCoordinate \\
             --chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimScoreMin 1 \\
             --chimOutType Junctions WithinBAM --outBAMsortingThreadN \"$thread_count\";

        # 2. Index the final analysis-ready BAM file.
        samtools index -@ \"$thread_count\" \"$final_bam\";

        # 3. Generate alignment statistics.
        samtools flagstat -@ \"$thread_count\" \"$final_bam\" > \"${step_dir}/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.flagstat\"
    "
    ( eval "$cmd" )
}

run_step_6() {
    # This function orchestrates the two-pass alignment and is wrapped by `run_step`.
    # It creates two directories, but the completion flag (`step.done`) is only
    # created in the final directory upon successful completion of both parts.
    local step_dir="${output_dir}/Step_06_Alignment_Refined"
    local cmd="
        # Create directories for both passes.
        mkdir -p \"${output_dir}/Step_06_Alignment_Initial\" && \\
        mkdir -p \"${output_dir}/Step_06_Alignment_Refined\" && \\
        # Execute the initial pass function.
        _run_step_6_initial && \\
        # Execute the refined pass function.
        _run_step_6_refined
    "
    run_step "$step_dir" "$cmd"
}
# --- End of Step 6 ---

run_step_7() {
    local step_dir="${output_dir}/Step_07_Strand_Detection"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    local cmd="
        python \"${EVscope_PATH}/bin/Step_07_bam2strand.py\" \\
            --input_bam \"$final_bam\" \\
            --bed \"$GENCODE_V45_BED\" \\
            --test_read_num 100000000 \\
            --output_dir \"$step_dir\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_8() {
    if [[ "$circ_tool" == "CIRCexplorer2" || "$circ_tool" == "both" ]]; then
        local step_dir="${output_dir}/Step_08_CIRCexplorer2_circRNA"
        local chimeric_junctions="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Chimeric.out.junction"
        local bsj_bed="${step_dir}/${sample_name}_back_spliced_junction.bed"
        local known_circs="${step_dir}/${sample_name}_circularRNA_known.txt"
        local final_output="${step_dir}/${sample_name}_CIRCexplorer2_dedup_junction_readcounts_CPM.tsv"
        local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
        local cmd="
            # 1. Parse STAR chimeric junctions to find back-spliced junctions (BSJ).
            CIRCexplorer2 parse -t STAR -b \"$bsj_bed\" \"$chimeric_junctions\";

            # 2. Annotate the BSJs using reference genome and gene models.
            CIRCexplorer2 annotate -r \"$GENCODE_V45_REFFLAT\" -g \"$HUMAN_GENOME_FASTA\" \\
                -b \"$bsj_bed\" -o \"$known_circs\";

            # 3. Convert raw counts to normalized CPM values.
            python \"${EVscope_PATH}/bin/Step_08_convert_CIRCexplorer2CPM.py\" \\
                --CIRCexplorer2_result \"$known_circs\" \\
                --input_bam \"$final_bam\" \\
                --GeneID_meta_table \"$TOTAL_GENEID_META\" \\
                --output \"$final_output\"
        "
        run_step "$step_dir" "$cmd"
    fi
}

run_step_9() {
    if [[ "$circ_tool" == "CIRI2" || "$circ_tool" == "both" ]]; then
        local step_dir="${output_dir}/Step_09_CIRI2_circRNA"
        local r1_dedup_fq="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_umi_dedup.clean.fq.gz"
        local r2_dedup_fq="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_umi_dedup.clean.fq.gz"
        local bwa_sam="${step_dir}/${sample_name}_umi_dedup.bwa.sam"
        local ciri2_out="${step_dir}/${sample_name}_CIRI2_out.tsv"
        local final_output="${step_dir}/${sample_name}_CIRI2_dedup_junction_readcounts_CPM.tsv"
        
        # Using a trap to ensure the large temporary SAM file is deleted even if the script fails.
        trap 'rm -f "$bwa_sam"' RETURN

        local cmd="
            # 1. Align deduplicated reads with BWA-MEM (required for CIRI2).
            bwa mem -t \"$thread_count\" -T 19 \"$BWA_INDEX\" \\
                \"$r1_dedup_fq\" ${is_paired_end:+\"$r2_dedup_fq\"} > \"$bwa_sam\" && \\

            # 2. Run CIRI2 to detect circRNAs from the BWA alignment.
            perl \"$CIRI2_PERL_SCRIPT\" -T \"$thread_count\" -I \"$bwa_sam\" \\
                -O \"$ciri2_out\" -F \"$HUMAN_GENOME_FASTA\" -A \"$GENCODE_V45_GTF\" && \\

            # 3. Convert raw CIRI2 counts to normalized CPM values.
            python \"${EVscope_PATH}/bin/Step_09_convert_CIRI2CPM.py\" \\
                --CIRI2_result \"$ciri2_out\" \\
                --input_sam \"$bwa_sam\" \\
                --output \"$final_output\" \\
                --GeneID_meta_table \"$TOTAL_GENEID_META\"
        "
        run_step "$step_dir" "$cmd"
    fi
}

run_step_10() {
    if [[ "$circ_tool" == "both" ]]; then
        local step_dir="${output_dir}/Step_10_circRNA_Merge"
        local cirexp_results="${output_dir}/Step_08_CIRCexplorer2_circRNA/${sample_name}_CIRCexplorer2_dedup_junction_readcounts_CPM.tsv"
        local ciri2_results="${output_dir}/Step_09_CIRI2_circRNA/${sample_name}_CIRI2_dedup_junction_readcounts_CPM.tsv"
        local cmd="
            python \"${EVscope_PATH}/bin/Step_10_circRNA_merge.py\" \\
                --CIRCexplorer2 \"$cirexp_results\" \\
                --CIRI2 \"$ciri2_results\" \\
                --output_matrix \"${step_dir}/${sample_name}_combined_CIRCexplorer2_CIRI2.tsv\" \\
                --out_venn \"${step_dir}/${sample_name}_Venn_diagram_of_circRNAs_identified_between_CIRCexplorer2_CIRI2.png\"
        "
        run_step "$step_dir" "$cmd"
    fi
}

run_step_11() {
    local step_dir="${output_dir}/Step_11_RNA_Metrics"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    local picard_metrics_tsv="${step_dir}/${sample_name}_picard_metrics.tsv"
    local insert_metrics_tsv="${step_dir}/${sample_name}_insert_size_metrics.tsv"
    local insert_histo_pdf="${step_dir}/${sample_name}_insert_size_histogram.pdf"
    local insert_histo_png="${step_dir}/${sample_name}_insert_size_histogram.png"

    local cmd="
        # 1. Run Picard CollectRnaSeqMetrics to get alignment and gene region statistics.
        conda run -n \"$PICARD_ENV\" picard -Xmx250g CollectRnaSeqMetrics \\
            I=\"$final_bam\" O=\"$picard_metrics_tsv\" REF_FLAT=\"$GENCODE_V45_REFFLAT\" \\
            STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=\"$HUMAN_RRNA_INTERVAL\" && \\

        # 2. Run Picard CollectInsertSizeMetrics to analyze fragment sizes.
        conda run -n \"$PICARD_ENV\" picard CollectInsertSizeMetrics \\
            I=\"$final_bam\" O=\"$insert_metrics_tsv\" H=\"$insert_histo_pdf\" && \\

        # 3. Convert the insert size PDF histogram to a PNG image for reports.
        conda run -n \"$PICARD_ENV\" convert -density 300 -background white -alpha remove \\
            \"$insert_histo_pdf\" \"$insert_histo_png\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_12() {
    local step_dir="${output_dir}/Step_12_featureCounts_Quant"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    local cmd="
        featureCounts -a \"$TOTAL_GENE_GTF\" \\
            -o \"${step_dir}/${sample_name}_featureCounts.tsv\" \\
            $featurecounts_paired -T \"$thread_count\" -s \"$featurecounts_strand\" \\
            -g gene_id -t exon -B -C \"$final_bam\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_13() {
    if [[ "$gDNA_correction" == "yes" ]]; then
        local step_dir="${output_dir}/Step_13_gDNA_Corrected_Quant"
        local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
        local forward_counts="${step_dir}/${sample_name}_featureCounts_forward.tsv"
        local reverse_counts="${output_dir}/Step_12_featureCounts_Quant/${sample_name}_featureCounts.tsv"
        local corrected_counts="${step_dir}/${sample_name}_gDNA_corrected_counts.tsv"
        local cmd="
            # 1. Count reads mapping to the sense strand (potential gDNA).
            featureCounts -a \"$TOTAL_GENE_GTF\" \\
                -o \"$forward_counts\" \\
                $featurecounts_paired -T \"$thread_count\" -s 1 -g gene_id -t exon -B -C \"$final_bam\" && \\

            # 2. Subtract sense counts from antisense counts to correct for gDNA.
            python \"${EVscope_PATH}/bin/Step_13_gDNA_corrected_featureCounts.py\" \\
                --strand \"$strandedness\" \\
                --forward_featureCounts_table \"$forward_counts\" \\
                --reverse_featureCounts_table \"$reverse_counts\" \\
                --output \"$corrected_counts\"
        "
        run_step "$step_dir" "$cmd"
    fi
}

run_step_14() {
    local step_dir="${output_dir}/Step_14_RSEM_Quant"
    local r1_dedup_fq="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R1_umi_dedup.clean.fq.gz"
    local r2_dedup_fq="${output_dir}/Step_06_Alignment_Initial/${sample_name}_R2_umi_dedup.clean.fq.gz"
    local cmd="
        # Run RSEM, which uses its own Bowtie2 alignment for quantification.
        # This handles multi-mapping reads differently than featureCounts.
        perl \"$RSEM_CALC_EXPR\" ${is_paired_end:+--paired-end} --bowtie2 \\
            --strandedness \"$strandedness\" --bowtie2-k 2 -p \"$thread_count\" \\
            --no-bam-output --seed 12345 \\
            \"$r1_dedup_fq\" ${is_paired_end:+\"$r2_dedup_fq\"} \\
            \"$RSEM_BOWTIE2_INDEX\" \\
            \"${step_dir}/${sample_name}_RSEM\" \\
            --temporary-folder \"${step_dir}/tmp\"
    "
    run_step "$step_dir" "$cmd"
}

# Generic function to run expression analysis and plotting for a given mode.
_run_expression_analysis() {
    local step_dir="$1"
    local input_counts_file="$2"
    local count_tool_script="$3" # e.g., Step_15_featureCounts2TPM.py
    local combined_matrix_suffix="$4"
    local base_expr_matrix="${step_dir}/${sample_name}_Gene_readcounts_normalized_expression_matrix_${combined_matrix_suffix}.tsv"
    local circ_expr_matrix="${output_dir}/Step_10_circRNA_Merge/${sample_name}_combined_CIRCexplorer2_CIRI2.tsv"
    local combined_matrix="${step_dir}/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_${combined_matrix_suffix}.tsv"

    local cmd="
        # 1. Convert raw counts to a standardized expression matrix (TPM, FPKM, etc.).
        python \"${EVscope_PATH}/bin/${count_tool_script}\" \\
            --featureCounts_out \"$input_counts_file\" --GeneID_meta_table \"$TOTAL_GENEID_META\" --output \"$base_expr_matrix\" && \\
        
        # 2. Combine linear RNA and circular RNA expression matrices.
        python \"${EVscope_PATH}/bin/Step_15_combine_total_RNA_expr_matrix.py\" \\
            --gene_expr \"$base_expr_matrix\" --circRNA_expr \"$circ_expr_matrix\" --out_matrix \"$combined_matrix\" && \\
        
        # 3. Generate RNA type distribution plots (1 subplot).
        python \"${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_1subplot.py\" \\
            --sample_name \"${sample_name}\" --Expr_matrix \"$combined_matrix\" --out_plot \"${step_dir}/${sample_name}_RNA_type_composition_1subplot.pdf\" && \\
        
        # 4. Generate RNA type distribution plots (2 subplots).
        python \"${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_2subplots.py\" \\
            --sample_name \"${sample_name}\" --Expr_matrix \"$combined_matrix\" --out_plot \"${step_dir}/${sample_name}_RNA_type_composition_2subplots.pdf\" && \\
        
        # 5. Generate RNA type distribution plots (20 subplots).
        python \"${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_20subplots.py\" \\
            --sample_name \"${sample_name}\" --Expr_matrix \"$combined_matrix\" --out_plot \"${step_dir}/${sample_name}_RNA_type_composition_20subplots.pdf\" && \\
        
        # 6. Generate a bar plot of the top 100 most highly expressed genes.
        python \"${EVscope_PATH}/bin/Step_15_plot_top_expressed_genes.py\" \\
            --input_gene_expr_matrix \"$combined_matrix\" --gene_num_per_type 5 --total_gene_num 100 \\
            --output_pdf \"${step_dir}/${sample_name}_bar_plot_for_top_100_highly_expressed_genes.pdf\" \\
            --output_png \"${step_dir}/${sample_name}_bar_plot_for_top_100_highly_expressed_genes.png\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_15() {
    _run_expression_analysis \
        "${output_dir}/Step_15_featureCounts_Expression" \
        "${output_dir}/Step_12_featureCounts_Quant/${sample_name}_featureCounts.tsv" \
        "Step_15_featureCounts2TPM.py" \
        "featureCounts"
}

run_step_16() {
    if [[ "$gDNA_correction" == "yes" ]]; then
        _run_expression_analysis \
            "${output_dir}/Step_16_gDNA_Corrected_Expression" \
            "${output_dir}/Step_13_gDNA_Corrected_Quant/${sample_name}_gDNA_corrected_counts.tsv" \
            "Step_15_featureCounts2TPM.py" \
            "gDNA_correction"
    fi
}

run_step_17() {
    # RSEM has a different input file format, so we use a specialized version of the helper function.
    local step_dir="${output_dir}/Step_17_RSEM_Expression"
    local rsem_results="${output_dir}/Step_14_RSEM_Quant/${sample_name}_RSEM.genes.results"
    local base_expr_matrix="${step_dir}/${sample_name}_Gene_readcounts_normalized_expression_matrix_RSEM.tsv"
    local circ_expr_matrix="${output_dir}/Step_10_circRNA_Merge/${sample_name}_combined_CIRCexplorer2_CIRI2.tsv"
    local combined_matrix="${step_dir}/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_RSEM.tsv"

    # Command is similar to the generic one, but uses a different script for step 1.
    local cmd="
        # 1. Convert RSEM results to a standardized expression matrix.
        python \"${EVscope_PATH}/bin/Step_17_RSEM2expr_matrix.py\" \\
            --RSEM_out \"$rsem_results\" --GeneID_meta_table \"$TOTAL_GENEID_META\" --output \"$base_expr_matrix\" && \\
        
        # 2. Combine linear and circRNA matrices.
        python \"${EVscope_PATH}/bin/Step_15_combine_total_RNA_expr_matrix.py\" \\
            --gene_expr \"$base_expr_matrix\" --circRNA_expr \"$circ_expr_matrix\" --out_matrix \"$combined_matrix\" && \\
        
        # 3-6. Generate all relevant RNA distribution and top expressed genes plots.
        python \"${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_1subplot.py\" --sample_name \"${sample_name}\" --Expr_matrix \"$combined_matrix\" --out_plot \"${step_dir}/${sample_name}_RNA_type_composition_1subplot.pdf\" && \\
        python \"${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_2subplots.py\" --sample_name \"${sample_name}\" --Expr_matrix \"$combined_matrix\" --out_plot \"${step_dir}/${sample_name}_RNA_type_composition_2subplots.pdf\" && \\
        python \"${EVscope_PATH}/bin/Step_15_plot_RNA_distribution_20subplots.py\" --sample_name \"${sample_name}\" --Expr_matrix \"$combined_matrix\" --out_plot \"${step_dir}/${sample_name}_RNA_type_composition_20subplots.pdf\" && \\
        python \"${EVscope_PATH}/bin/Step_15_plot_top_expressed_genes.py\" --input_gene_expr_matrix \"$combined_matrix\" --gene_num_per_type 5 --total_gene_num 100 --output_pdf \"${step_dir}/${sample_name}_bar_plot_for_top_100_highly_expressed_genes.pdf\" --output_png \"${step_dir}/${sample_name}_bar_plot_for_top_100_highly_expressed_genes.png\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_18() {
    local step_dir="${output_dir}/Step_18_Genomic_Regions"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    
    local cmd="
        # 1. Run featureCounts in a loop for each genomic region defined by a SAF file.
        for saf_path in \\
            \"\$HUMAN_3UTR_SAF\" \\
            \"\$HUMAN_5UTR_SAF\" \\
            \"\$HUMAN_DOWNSTREAM_2KB_SAF\" \\
            \"\$HUMAN_BLACKLIST_SAF\" \\
            \"\$HUMAN_EXON_SAF\" \\
            \"\$HUMAN_INTERGENIC_SAF\" \\
            \"\$HUMAN_INTRON_SAF\" \\
            \"\$HUMAN_PROMOTER_SAF\"
        do
            saf_basename=\$(basename \"\$saf_path\" .saf)
            featureCounts -F SAF -a \"\$saf_path\" \\
                -o \"${step_dir}/${sample_name}_\${saf_basename}_featureCounts.tsv\" \\
                \$featurecounts_paired -T \$thread_count -s \$featurecounts_strand -B -C \"$final_bam\";
        done && \\
        
        # 2. Aggregate the counts from all regions and generate a summary plot.
        python \"${EVscope_PATH}/bin/Step_18_plot_reads_mapping_stats.py\" \\
            --input_5UTR_readcounts \"${step_dir}/${sample_name}_HG38_5UTR_noOverlap_featureCounts.tsv\" \\
            --input_exon_readcounts \"${step_dir}/${sample_name}_HG38_exon_noOverlap_featureCounts.tsv\" \\
            --input_3UTR_readcounts \"${step_dir}/${sample_name}_HG38_3UTR_noOverlap_featureCounts.tsv\" \\
            --input_intron_readcounts \"${step_dir}/${sample_name}_HG38_intron_noOverlap_featureCounts.tsv\" \\
            --input_promoters_readcounts \"${step_dir}/${sample_name}_HG38_promoter_1500_500bp_noOverlap_featureCounts.tsv\" \\
            --input_downstream_2Kb_readcounts \"${step_dir}/${sample_name}_HG38_downstream_2kb_noOverlap_featureCounts.tsv\" \\
            --input_intergenic_readcounts \"${step_dir}/${sample_name}_HG38_intergenic_noOverlap.tsv\" \\
            --input_ENCODE_blacklist_readcounts \"${step_dir}/${sample_name}_HG38_ENCODE_blacklist_V2_featureCounts.tsv\" \\
            --sampleName \"${sample_name}\" \\
            --output_dir \"$step_dir\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_19() {
    local step_dir="${output_dir}/Step_19_Taxonomy"
    local r1_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local r2_clean="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local r1_downsampled="${step_dir}/${sample_name}_R1_downsampled.fq.gz"
    local r2_downsampled="${step_dir}/${sample_name}_R2_downsampled.fq.gz"
    local kraken_report="${step_dir}/${sample_name}_report.tsv"
    local krona_input="${step_dir}/${sample_name}_krona_input.tsv"
    local krona_html="${step_dir}/${sample_name}_krona.html"
    
    local cmd="
        # 1. Downsample clean FASTQs to 100k reads for faster classification.
        seqtk sample -s100 \"$r1_clean\" 100000 | gzip > \"$r1_downsampled\" && \\
        ${is_paired_end:+seqtk sample -s100 \"$r2_clean\" 100000 | gzip > \"$r2_downsampled\"}

        # 2. Run Kraken2 for taxonomic classification.
        conda run -n \"$KRAKEN2_ENV\" kraken2 --db \"$KRAKEN_DB\" --threads \"$thread_count\" \\
            --report \"$kraken_report\" ${is_paired_end:+--paired} --gzip-compressed \\
            \"$r1_downsampled\" ${is_paired_end:+\"$r2_downsampled\"} && \\
            
        # 3. Convert Kraken2 report to Krona input format.
        python \"${KRAKEN_TOOLS_DIR}/kreport2krona.py\" -r \"$kraken_report\" -o \"$krona_input\" && \\

        # 4. Generate an interactive Krona HTML plot.
        conda run -n \"$KRAKEN2_ENV\" ktImportText \"$krona_input\" -o \"$krona_html\"
    "
    # Ignore errors (true) because this step can fail if no non-human reads are found.
    run_step "$step_dir" "$cmd" "true"
}

# Generic function to run deconvolution analysis.
_run_deconvolution() {
    local step_dir="$1"
    local source_expr_matrix="$2"
    local step_num="$3"
    local tpm_matrix="${step_dir}/${sample_name}_TPM_matrix.csv"

    local cmd="
        # 1. Extract just the GeneID and TPM columns needed for deconvolution.
        awk -F'\t' '{print \$1\",\"\$5}' \"$source_expr_matrix\" > \"$tpm_matrix\";

        # 2. Loop through all reference expression profiles defined in the config file.
        for ref_file in \\
            \"\$GTEX_SMTS_REF\" \\
            \"\$GTEX_SMTSD_REF\" \\
            \"\$BRAIN_SC_REF\"
        do
            # Robustness check: ensure the reference file is defined and exists.
            if [[ -n \"\$ref_file\" && -f \"\$ref_file\" ]]; then
                log 2 'INFO' \"Running deconvolution for Step ${step_num} with reference: \$(basename \$ref_file)\"
                python \"${EVscope_PATH}/bin/Step_22_run_RNA_deconvolution_ARIC.py\" \\
                    --bulk_expr_file \"$tpm_matrix\" \\
                    --ref_expr_file \"\$ref_file\" \\
                    --output_dir \"$step_dir\";
            elif [[ -n \"\$ref_file\" ]]; then
                # Warn the user if a configured file is not found.
                log 3 'WARN' \"Reference file for Step ${step_num} not found: '\$ref_file'. Skipping this reference.\"
            fi
        done
    "
    # Ignore errors (true) as this is an auxiliary analysis and may fail if reference files are missing.
    run_step "$step_dir" "$cmd" "true"
}

run_step_20() {
    _run_deconvolution \
        "${output_dir}/Step_20_featureCounts_Deconvolution" \
        "${output_dir}/Step_15_featureCounts_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_featureCounts.tsv" \
        "20"
}

run_step_21() {
    if [[ "$gDNA_correction" == "yes" ]]; then
        _run_deconvolution \
            "${output_dir}/Step_21_gDNA_Corrected_Deconvolution" \
            "${output_dir}/Step_16_gDNA_Corrected_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_gDNA_correction.tsv" \
            "21"
    fi
}

run_step_22() {
    _run_deconvolution \
        "${output_dir}/Step_22_RSEM_Deconvolution" \
        "${output_dir}/Step_17_RSEM_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_RSEM.tsv" \
        "22"
}

run_step_23() {
    local step_dir="${output_dir}/Step_23_rRNA_Detection"
    local input_r1="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    local input_r2="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    local output_log="${step_dir}/${sample_name}_ribodetector_summary.log"
    local rrna_output_r1="${step_dir}/${sample_name}_rRNA_R1.fq.gz"
    local rrna_output_r2="${step_dir}/${sample_name}_rRNA_R2.fq.gz"

    local cmd="
        # Run ribodetector to identify and separate rRNA reads.
        # Non-rRNA reads are discarded to /dev/null; we only want the rRNA reads.
        ribodetector_cpu \\
            -t \"$thread_count\" \\
            -l 100 \\
            -i \"$input_r1\" ${is_paired_end:+\"$input_r2\"} \\
            -e rrna \\
            --chunk_size 800 \\
            --log \"$output_log\" \\
            -o /dev/null ${is_paired_end:+/dev/null} \\
            -r \"$rrna_output_r1\" ${is_paired_end:+-r \"$rrna_output_r2\"} || \\
        {
            # Fallback: If ribodetector exits with an error (e.g., no rRNA found),
            # create empty output files to prevent downstream steps from failing.
            log 3 'WARN' 'ribodetector may have failed because no rRNA was detected. Creating empty output files to ensure pipeline continuation.'
            touch \"$rrna_output_r1\"
            if [ \"$is_paired_end\" = true ]; then
                touch \"$rrna_output_r2\"
            fi
        }
    "
    # Ignore errors for robustness.
    run_step "$step_dir" "$cmd" "true"
}

run_step_24() {
    local step_dir="${output_dir}/Step_24_QC_Summary"
    local qc_script="${EVscope_PATH}/bin/Step_24_generate_QC_matrix.py"
    # Start building the command with the mandatory output argument.
    local qc_cmd="python \"$qc_script\" --output \"${step_dir}/${sample_name}_QC_matrix.tsv\""

    # Dynamically append optional arguments only if the corresponding input files exist.
    # This makes the step robust to which previous steps were run.
    local raw_fastqc_files=("${output_dir}/Step_01_Raw_QC/"*".zip")
    [ -e "${raw_fastqc_files[0]}" ] && qc_cmd+=" --raw_fastqc_zips ${raw_fastqc_files[*]}"
    
    local acc_motif_file="${output_dir}/Step_02_UMI_Analysis/${sample_name}_ACC_motif_fraction.tsv"
    [ -f "$acc_motif_file" ] && qc_cmd+=" --ACC_motif_fraction \"$acc_motif_file\""
    
    local trimmed_r1="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R1_clean.fq.gz"
    [ -f "$trimmed_r1" ] && qc_cmd+=" --trimmed_R1_fastq \"$trimmed_r1\""
    
    local trimmed_r2="${output_dir}/Step_03_UMI_Adaptor_Trim/${sample_name}_R2_clean.fq.gz"
    [ -f "$trimmed_r2" ] && qc_cmd+=" --trimmed_R2_fastq \"$trimmed_r2\""

    local ecoli_r1=("${output_dir}/Step_05_Bacterial_Filter/"*"coli"*"_R1.fq.gz")
    [ -e "${ecoli_r1[0]}" ] && qc_cmd+=" --ecoli_R1_fastq \"${ecoli_r1[0]}\""

    local myco_r1=("${output_dir}/Step_05_Bacterial_Filter/"*"myco"*"_R1.fq.gz")
    [ -e "${myco_r1[0]}" ] && qc_cmd+=" --myco_R1_fastq \"${myco_r1[0]}\""
    
    local star_log="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Log.final.out"
    [ -f "$star_log" ] && qc_cmd+=" --STAR_log \"$star_log\""
    
    local strand_file="${output_dir}/Step_07_Strand_Detection/bam2strand.xls"
    [ -f "$strand_file" ] && qc_cmd+=" --bam2strand_file \"$strand_file\""

    local picard_insert="${output_dir}/Step_11_RNA_Metrics/${sample_name}_insert_size_metrics.tsv"
    [ -f "$picard_insert" ] && qc_cmd+=" --picard_insert_file \"$picard_insert\""
    
    local picard_rnaseq="${output_dir}/Step_11_RNA_Metrics/${sample_name}_picard_metrics.tsv"
    [ -f "$picard_rnaseq" ] && qc_cmd+=" --picard_rnaseq_file \"$picard_rnaseq\""
    
    # Find the primary expression matrix, preferring the gDNA-corrected or standard featureCounts one.
    local expr_matrix=""
    if [ -f "${output_dir}/Step_16_gDNA_Corrected_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_gDNA_correction.tsv" ]; then
        expr_matrix="${output_dir}/Step_16_gDNA_Corrected_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_gDNA_correction.tsv"
    elif [ -f "${output_dir}/Step_15_featureCounts_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_featureCounts.tsv" ]; then
        expr_matrix="${output_dir}/Step_15_featureCounts_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_featureCounts.tsv"
    elif [ -f "${output_dir}/Step_17_RSEM_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_RSEM.tsv" ]; then
        expr_matrix="${output_dir}/Step_17_RSEM_Expression/${sample_name}_combined_expression_matrix_linearRNA_TPM_circRNA_CPM_RSEM.tsv"
    fi
    [ -n "$expr_matrix" ] && qc_cmd+=" --expression_matrix \"$expr_matrix\""

    # Add genomic region counts
    for region in 3UTR 5UTR downstream_2kb exon ENCODE_blacklist intergenic intron promoter_1500_500bp; do
        local fc_file_path
        fc_file_path=$(find "${output_dir}/Step_18_Genomic_Regions/" -name "${sample_name}*${region}*featureCounts.tsv" 2>/dev/null | head -n 1)
        [ -f "$fc_file_path" ] && qc_cmd+=" --featureCounts_${region} \"$fc_file_path\""
    done
    
    local kraken_report="${output_dir}/Step_19_Taxonomy/${sample_name}_report.tsv"
    [ -f "$kraken_report" ] && qc_cmd+=" --kraken_report \"$kraken_report\""

    local ribo_r1=("${output_dir}/Step_23_rRNA_Detection/"*"_rRNA_R1.fq.gz")
    [ -e "${ribo_r1[0]}" ] && qc_cmd+=" --ribo_R1_fastq \"${ribo_r1[0]}\""

    local downsampled_r1="${output_dir}/Step_19_Taxonomy/${sample_name}_R1_downsampled.fq.gz"
    [ -f "$downsampled_r1" ] && qc_cmd+=" --downsampled_trimmed_R1_fastq \"$downsampled_r1\""

    run_step "$step_dir" "$qc_cmd"
}

run_step_25() {
    local step_dir="${output_dir}/Step_25_EMapper_BigWig_Quantification"
    local final_bam="${output_dir}/Step_06_Alignment_Refined/${sample_name}_STAR_umi_dedup_Aligned.sortedByCoord.out.bam"
    local cmd="
        python \"${EVscope_PATH}/bin/Step_25_EMapper.py\" \\
            --num_threads 4 \\
            --input_bam \"$final_bam\" \\
            --sample_name \"${sample_name}\" \\
            --output_dir \"${step_dir}/EMapper_output\"
    "
    run_step "$step_dir" "$cmd"
}

run_step_26() {
    local step_dir="${output_dir}/Step_26_BigWig_Density_Plot"
    local bigwig_file="${output_dir}/Step_25_EMapper_BigWig_Quantification/EMapper_output/${sample_name}_final_unstranded.bw"
    # Format the long lists of BED files for readability.
    local bed_files_rna_types="[\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_protein_coding.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_tRNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_rRNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_Y_RNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_snoRNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_snRNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_scaRNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_vault_RNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_lncRNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_TEC_protein_coding.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_IG_genes.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_TR_genes.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_misc-sncRNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_miRNAs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_pseudogenes.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_ERVs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_LINEs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_SINEs.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_RNAtype_bed/HG38_piRNAs_gold_standard.bed]"
    local bed_labels_rna_types="[protein_coding,tRNAs,rRNAs,Y_RNAs,snoRNAs,snRNAs,scaRNAs,vault_RNAs,lncRNAs,TEC_protein_coding,IG_genes,TR_genes,misc-sncRNAs,miRNAs,pseudogenes,ERVs,LINEs,SINEs,piRNAs]"

    local bed_files_meta_gene="[\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_meta_bed/HG38_promoter_1500_500bp_noOverlap.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_meta_bed/HG38_5UTR_noOverlap.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_meta_bed/HG38_exon_noOverlap.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_meta_bed/HG38_intron_noOverlap.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_meta_bed/HG38_3UTR_noOverlap.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_meta_bed/HG38_downstream_2kb_noOverlap.bed,\\
        \${EVscope_PATH}/references/annotations_HG38/GeneCodeV45_meta_bed/HG38_intergenic_noOverlap.bed]"
    local bed_labels_meta_gene="[Promoters,5'UTR,Exon,Intron,3'UTR,Downstream,Intergenic]"

    local cmd="
        # 1. Generate density plots of signal coverage over different RNA biotypes.
        bash \"\${EVscope_PATH}/bin/Step_26_density_plot_over_RNA_types.sh\" \\
            --input_bw_file \"$bigwig_file\" \\
            --input_bed_files \"$bed_files_rna_types\" \\
            --input_bed_labels \"$bed_labels_rna_types\" \\
            --output_dir \"${step_dir}/RNA_types\" \\
            --random_tested_row_num_per_bed 50000 && \\

        # 2. Generate density plots of signal coverage over meta-gene regions (promoter, UTR, exon, etc.).
        bash \"\${EVscope_PATH}/bin/Step_26_density_plot_over_meta_gene.sh\" \\
            --input_bw_file \"$bigwig_file\" \\
            --input_bed_files \"$bed_files_meta_gene\" \\
            --input_bed_labels \"$bed_labels_meta_gene\" \\
            --output_dir \"${step_dir}/meta_gene\" \\
            --blackListFileName \"\${ENCODE_BLACKLIST_BED}\" \\
            --random_tested_row_num_per_bed 50000
    "
    run_step "$step_dir" "$cmd"
}

run_step_27() {
    local step_dir="${output_dir}/Step_27_HTML_Report"
    # Use `readlink -f` to get absolute paths, required for R/rmarkdown context.
    local abs_output_dir; abs_output_dir=$(readlink -f "${output_dir}")
    local abs_bin_dir; abs_bin_dir=$(readlink -f "${EVscope_PATH}/bin")
    local rmd_input="${abs_bin_dir}/Step_27_html_report.Rmd"
    local report_file="${sample_name}_final_report.html"
    
    # Construct the R command to be executed.
    local r_command="rmarkdown::render(input = '$rmd_input', output_file = '$report_file', output_dir = '$step_dir')"

    # The Rmd script needs an environment variable to find all the output files.
    # Exporting it makes it available to the Rscript process.
    local full_cmd="export EVscope_OUTPUT_DIR=\"$abs_output_dir\"; Rscript -e \"$r_command\""
    run_step "$step_dir" "$full_cmd"
}

# ==============================================================================
# --- MAIN EXECUTION LOGIC ---
# ==============================================================================

# Set a higher limit for open files, required by tools like STAR.
ulimit -n 655350

main() {
    # --- 1. Initialize Default Parameters and Parse Command-Line Arguments ---
    local config_file="${SCRIPT_DIR}/EVscope.conf"
    local thread_count=1
    local sample_name=""
    local run_steps="all" # Default is to run all steps.
    local skip_steps=""
    local circ_tool="both"
    local read_count_mode="uniq"
    local gDNA_correction="no"
    local strandedness="reverse"
    local -a input_fastqs=()

    # Define command-line options for `getopt`.
    local SHORT_OPTS="hvV:"
    local LONG_OPTS="sample_name:,input_fastqs:,threads:,run_steps:,skip_steps:,circ_tool:,read_count_mode:,gDNA_correction:,strandedness:,config:,verbosity:,help,version"
    
    local parsed_args
    parsed_args=$(getopt -o "$SHORT_OPTS" --long "$LONG_OPTS" -n "$SCRIPT_NAME" -- "$@")
    if [[ $? -ne 0 ]]; then print_help; fi
    eval set -- "$parsed_args"

    while true; do
        case "$1" in
            --sample_name) sample_name="$2"; shift 2 ;;
            --input_fastqs) IFS=',' read -r -a input_fastqs <<< "$2"; shift 2 ;;
            --threads) thread_count="$2"; shift 2 ;;
            --run_steps) run_steps="$2"; shift 2 ;;
            --skip_steps) skip_steps="$2"; shift 2 ;;
            --circ_tool) circ_tool="$2"; shift 2 ;;
            --read_count_mode) read_count_mode="$2"; shift 2 ;;
            --gDNA_correction) gDNA_correction="$2"; shift 2 ;;
            --strandedness) strandedness="$2"; shift 2 ;;
            --config) config_file="$2"; shift 2 ;;
            -V|--verbosity) verbosity="$2"; shift 2 ;;
            -h|--help) print_help ;;
            -v|--version) print_version ;;
            --) shift; break ;;
            *) echo "Internal error in argument parsing!"; exit 1 ;;
        esac
    done

    # --- 2. Set Up Logging, Output Directory, and Load Configuration ---
    if [[ -z "$sample_name" ]]; then
        log 5 "FATAL" "--sample_name is a required argument."
        print_help
    fi
    # Sanitize sample name to create a valid directory name.
    local sample_name_sanitized
    sample_name_sanitized=$(echo "$sample_name" | sed 's/[^a-zA-Z0-9_.-]/_/g')
    output_dir="${sample_name_sanitized}_EVscope_output"
    
    mkdir -p "$output_dir" || { echo -e "${C_RED}FATAL: Cannot create output directory: $output_dir${C_RESET}"; exit 1; }
    touch "${output_dir}/EVscope_pipeline.log" # Create log file early.

    if [[ ! -f "$config_file" ]]; then
        log 5 "FATAL" "Configuration file not found: '$config_file'"
        exit 1
    fi
    # Source the configuration file to load all paths and settings.
    # shellcheck source=/dev/null
    source "$config_file"
    log 2 "INFO" "Loaded configuration from '$config_file'"

    # --- 3. Pre-run Checks for Dependencies and Essential Files ---
    check_dependencies
    local essential_refs=("$HUMAN_GENOME_FASTA" "$STAR_INDEX" "$TOTAL_GENE_GTF")
    for ref in "${essential_refs[@]}"; do
        if [[ ! -e "$ref" ]]; then
            log 5 "FATAL" "Essential reference file or directory not found: $ref"
            exit 1
        fi
    done
    log 2 "INFO" "Essential reference files are present."

    # --- 4. Validate User Inputs ---
    if [[ ${#input_fastqs[@]} -eq 0 ]]; then
        log 5 "FATAL" "--input_fastqs is a required argument."
        print_help
    fi
    for fastq in "${input_fastqs[@]}"; do
        if [[ ! -f "$fastq" ]]; then
            log 5 "FATAL" "Input FASTQ file not found: $fastq"
            exit 1
        fi
    done

    # --- 5. Determine the Final List of Steps to Execute ---
    local run_steps_list
    run_steps_list=$(parse_steps "$run_steps")
    
    local skip_steps_list
    skip_steps_list=$(parse_steps "$skip_steps")
    
    # Filter the run list by removing any steps present in the skip list.
    local final_steps_array=()
    comm -23 <(echo "$run_steps_list") <(echo "$skip_steps_list") | while read -r step; do
        final_steps_array+=("$step")
    done

    # Further filter out steps that are conditional on other settings (e.g., gDNA correction).
    local steps_to_run=()
    for step in "${final_steps_array[@]}"; do
        if [[ "$gDNA_correction" == "no" && ( "$step" == "13" || "$step" == "16" || "$step" == "21" ) ]]; then
            log 1 "DEBUG" "Skipping step $step because --gDNA_correction is 'no'."
            continue
        fi
        steps_to_run+=("$step")
    done
    
    if [[ ${#steps_to_run[@]} -eq 0 && -n "$run_steps" ]]; then
        log 5 "FATAL" "No steps to run after applying filters. Check your --run_steps and --skip_steps options."
        exit 1
    fi

    # --- 6. Configure Run-time Variables Based on Inputs ---
    is_paired_end=false
    fastq_read1="${input_fastqs[0]}"; fastq_read2=""
    if [[ ${#input_fastqs[@]} -eq 2 ]]; then
        is_paired_end=true
        fastq_read2="${input_fastqs[1]}"
    fi

    # Convert user-friendly strandedness option to featureCounts' numeric code.
    featurecounts_strand=2 # Default: reverse
    if [[ "$strandedness" == "forward" ]]; then
        featurecounts_strand=1
    elif [[ "$strandedness" == "unstrand" ]]; then
        featurecounts_strand=0
    fi

    # Set paired-end flags for featureCounts.
    featurecounts_paired=""
    if [[ "$is_paired_end" == true ]]; then
        featurecounts_paired="-p --countReadPairs"
    fi

    # --- 7. Execute Pipeline ---
    local start_time; start_time=$(date +%s)
    log 2 "INFO" "Starting EVscope pipeline (Version: $VERSION) for sample: ${C_BOLD}${sample_name}${C_RESET}"
    log 1 "DEBUG" "Full command: $0 $*"
    log 2 "INFO" "Output will be saved to: ${C_BOLD}${output_dir}${C_RESET}"
    log 2 "INFO" "Final steps to be executed: ${steps_to_run[*]}"
    print_pipeline_steps

    for step in "${steps_to_run[@]}"; do
        "run_step_$step"
    done

    local total_time; total_time=$(( $(date +%s) - start_time ))
    log 2 "INFO" "Pipeline completed successfully for ${C_BOLD}${sample_name}${C_RESET} in ${total_time} seconds."
    log 2 "INFO" "Final output is in: ${C_BOLD}${output_dir}${C_RESET}"
}

# --- Script Entrypoint ---
# Pass all command-line arguments to the main function.
main "$@"