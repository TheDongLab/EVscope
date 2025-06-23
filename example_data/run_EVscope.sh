#!/usr/bin/env bash
#
# run_EVscope.sh
#
# ==============================================================================
# Description:
#   A wrapper script to execute the main EVscope pipeline (EVscope.sh) with
#   a pre-defined set of parameters for a specific sample. It automatically
#   handles the argument formatting for both single-end and paired-end
#   FASTQ data.
#
# Usage:
#   ./run_EVscope.sh
#
# Prerequisites:
#   - The main pipeline script (EVscope.sh) must be accessible via the
#     relative path defined in the `PIPELINE_SCRIPT` variable.
#   - Input FASTQ files must exist at the specified paths.
# ==============================================================================

# --- Script Behavior ---
# Exit immediately if a command exits with a non-zero status.
# Treat unset variables as an error when substituting.
# Pipelines return the exit status of the last command to fail.
set -euo pipefail

# --- Configuration: Script and File Paths ---
# Get the absolute directory of the current script to robustly locate the pipeline script.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
PIPELINE_SCRIPT="${SCRIPT_DIR}/../../EVscope.sh"

# --- Configuration: Sample Information ---
# The name of the sample to be processed.
SAMPLE_NAME="Example_Data"
# Full path to the forward (R1) or single-end FASTQ file.
FASTQ_R1="Example_Data/chr21_2000_reads_R1_001.fastq"
# Full path to the reverse (R2) FASTQ file. Leave empty ("") for single-end data.
FASTQ_R2="xample_Data/chr21_2000_reads_R2_001.fastq"

# --- Configuration: Pipeline Execution Parameters ---
# Number of threads to allocate for the pipeline run.
THREADS=20
# Specific steps of the pipeline to execute (e.g., "1-27", "1,3,5").
STEPS_TO_RUN="1-27"
# Verbosity level for logging: 1(DEBUG), 2(INFO), 3(WARN), 4(ERROR).
VERBOSITY_LEVEL=2

# --- Configuration: Pipeline Scientific Parameters ---
# Enable genomic DNA (gDNA) contamination correction. ("yes" or "no").
GDNA_CORRECTION="yes"
# Library strandedness ("forward", "reverse", or "unstranded").
STRANDEDNESS="reverse"


# ==============================================================================
# Main Execution Logic
#
# Do not edit below this line unless you know what you are doing.
# ==============================================================================
main() {
    # Validate that the main pipeline script exists and is executable.
    if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
        echo "Error: Pipeline script not found at the expected path: ${PIPELINE_SCRIPT}" >&2
        exit 1
    fi

    # This array will store all arguments to be passed to the pipeline script.
    # Using an array is a robust way to handle arguments that may contain spaces or special characters.
    local -a pipeline_args=()

    # Prepare the --input_fastqs argument based on whether R2 is provided.
    if [[ -n "$FASTQ_R2" ]]; then
        # Paired-end mode: Combine R1 and R2 paths into a single, comma-separated string.
        echo "INFO: Paired-end data detected. Preparing input FASTQs argument."
        pipeline_args+=(--input_fastqs "${FASTQ_R1},${FASTQ_R2}")
    else
        # Single-end mode.
        echo "INFO: Single-end data detected. Preparing input FASTQs argument."
        pipeline_args+=(--input_fastqs "$FASTQ_R1")
    fi

    # Add all other configured parameters to the arguments array.
    pipeline_args+=(
        --sample_name "$SAMPLE_NAME"
        --threads "$THREADS"
        --run_steps "$STEPS_TO_RUN"
        --gDNA_correction "$GDNA_CORRECTION"
        --strandedness "$STRANDEDNESS"
        --verbosity "$VERBOSITY_LEVEL"
    )

    echo "==================================================================="
    echo "Starting EVscope Pipeline for Sample: ${SAMPLE_NAME}"
    echo "==================================================================="

    # For logging and debugging, print the exact command that will be executed.
    # 'printf "%q "' is used to safely quote each argument for display.
    echo "Command to be executed:"
    printf "bash %q " "$PIPELINE_SCRIPT"
    printf "%q " "${pipeline_args[@]}"
    echo -e "\n-------------------------------------------------------------------\n"

    # Execute the main pipeline script with the prepared arguments.
    bash "$PIPELINE_SCRIPT" "${pipeline_args[@]}"

    echo -e "\n-------------------------------------------------------------------"
    echo "Pipeline execution finished for sample: ${SAMPLE_NAME}."
    echo "==================================================================="
}

# --- Script Entry Point ---
# The "$@" allows passing arguments from this script to the main function,
# though in this specific setup, no arguments are expected.
main "$@"