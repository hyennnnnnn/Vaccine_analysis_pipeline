#!/usr/bin/env bash
set -euo pipefail

# GC Pipeline (Parallel Mode) - Influenza Sequencing Data Analysis Pipeline
# Processes all samples step-by-step with parallel execution within each step
#
# Prerequisites: Run setup.sh first to create project directory structure

# Trap errors
trap 'echo "Error at line $LINENO" >&2' ERR

# === Conda Environment Activation ===
# Conda environment for the pipeline
# vac_env should contain:
#   - Python 3.9+
#   - pandas, numpy, matplotlib (for data analysis and plotting)
#   - biopython (for sequence analysis)
#   - pillow (for image processing)
# Use environment variable GC_CONDA_ENV if set, otherwise default to "vac_env"
CONDA_ENV="${GC_CONDA_ENV:-vac_env}"

# Initialize and activate conda environment
if command -v conda &> /dev/null; then
    # Initialize conda
    eval "$(conda shell.bash hook)" || {
        echo "❌ ERROR: Failed to initialize conda" >&2
        exit 1
    }

    # Activate conda environment
    echo "[SETUP] Activating conda environment: ${CONDA_ENV}" >&2
    conda activate "${CONDA_ENV}" || {
        echo "❌ ERROR: Failed to activate conda environment '${CONDA_ENV}'" >&2
        echo "Please ensure the environment exists: conda env list" >&2
        exit 1
    }

    # Verify activation
    if [[ "$CONDA_DEFAULT_ENV" == "$CONDA_ENV" ]]; then
        echo "[SETUP] ✓ Conda environment activated: ${CONDA_ENV}" >&2
        echo "[SETUP] ✓ Python: $(python3 --version 2>&1)" >&2
    else
        echo "❌ ERROR: Conda environment activation verification failed" >&2
        echo "Expected: ${CONDA_ENV}" >&2
        echo "Got: ${CONDA_DEFAULT_ENV}" >&2
        exit 1
    fi
else
    echo "❌ ERROR: conda command not found" >&2
    echo "Please ensure conda is installed and in PATH" >&2
    exit 1
fi

# === Configuration ===
VERBOSE=${VERBOSE:-1}
START_TIME=$(date +%s)
PARALLEL_JOBS=${PARALLEL_JOBS:-8}  # Number of samples to process in parallel

# === Logging Functions ===
log_step_start() {
    local step=$1
    local desc=$2
    echo "" >&2
    echo "========================================================================" >&2
    echo "[PIPELINE] $(date '+%Y-%m-%d %H:%M:%S') - Step $step: $desc" >&2
    echo "========================================================================" >&2
}

log_sample_step() {
    local step=$1
    local total=$2
    local desc=$3
    local sample=$4
    echo "[$(date '+%H:%M:%S')] ▶️  [$step/$total] $desc - Sample: $sample" >&2
}

log_error() {
    echo "[$(date '+%H:%M:%S')] ❌ Error: $1" >&2
}

log_skip() {
    echo "[$(date '+%H:%M:%S')] ⏭️  $1" >&2
}

log_success() {
    echo "[$(date '+%H:%M:%S')] ✅ $1" >&2
}

show_elapsed_time() {
    local current_time=$(date +%s)
    local elapsed=$((current_time - START_TIME))
    local minutes=$((elapsed / 60))
    local seconds=$((elapsed % 60))
    echo "⏱️  Total elapsed time: ${minutes}m ${seconds}s" >&2
}

# === Parse Arguments ===
SAMPLE_LIST_FILE=""
STRAIN=""
REGION=""
YEAR=""
PARENT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --samples)
            SAMPLE_LIST_FILE="$2"
            shift 2
            ;;
        --strain)
            STRAIN="$2"
            shift 2
            ;;
        --region)
            REGION="$2"
            shift 2
            ;;
        --year)
            YEAR="$2"
            shift 2
            ;;
        --jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        --parent-dir)
            PARENT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1" >&2
            exit 1
            ;;
    esac
done

# Use environment variable GC_PARENT_DIR if --parent-dir not specified
# Default to current directory if neither provided
PARENT_DIR="${PARENT_DIR:-${GC_PARENT_DIR:-$(pwd)}}"

# Usage check
if [[ -z "$REGION" || -z "$YEAR" || -z "$SAMPLE_LIST_FILE" || -z "$STRAIN" ]]; then
    echo "Usage: bash GC_pipeline_parallel.sh \\" >&2
    echo "  --region <region> \\" >&2
    echo "  --year <year> \\" >&2
    echo "  --strain <strain> \\" >&2
    echo "  --samples <sample_list_file> \\" >&2
    echo "  [--parent-dir <parent_directory>] \\" >&2
    echo "  [--jobs <parallel_jobs>]  # Default: 8" >&2
    echo "" >&2
    echo "Note: If --parent-dir is not specified, uses environment variable GC_PARENT_DIR" >&2
    echo "      or defaults to current directory" >&2
    echo "" >&2
    echo "Example:" >&2
    echo "  bash GC_pipeline_parallel.sh \\" >&2
    echo "    --region southern \\" >&2
    echo "    --year 2024 \\" >&2
    echo "    --strain A_H3N2 \\" >&2
    echo "    --parent-dir /data/projects \\" >&2
    echo "    --samples /data/projects/GreenCross_flu_southern_2024/A_H3N2/samples.txt \\" >&2
    echo "    --jobs 16" >&2
    echo "" >&2
    echo "Note: Run 'bash setup.sh --region <region> --year <year> --strain <strain> --parent-dir <parent_dir>' first" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR/scripts"

export BASE_PROJECT_DIR="${PARENT_DIR}/GreenCross_flu_${REGION}_${YEAR}"
export STRAIN_DIR="${BASE_PROJECT_DIR}/${STRAIN}"
PIPELINE_DIR="${STRAIN_DIR}/.pipeline"
mkdir -p "$PIPELINE_DIR"

# === Check if Setup was Completed ===
SETUP_MARKER="$PIPELINE_DIR/setup_${REGION}_${YEAR}_${STRAIN}_done"
if [[ ! -f "$SETUP_MARKER" ]]; then
    echo "========================================================================" >&2
    echo "❌ Error: Project not set up yet!" >&2
    echo "========================================================================" >&2
    echo "" >&2
    echo "Please run setup first:" >&2
    echo "  bash setup.sh --region ${REGION} --year ${YEAR} --strain ${STRAIN}" >&2
    echo "" >&2
    exit 1
fi

# === Source Configuration ===
if [[ ! -f "${SCRIPTS_DIR}/config.sh" ]]; then
    echo "Error: config.sh not found at ${SCRIPTS_DIR}/config.sh" >&2
    exit 1
fi
source "${SCRIPTS_DIR}/config.sh"

# === Pre-run Checks ===
if [ ! -d "${STRAIN_DIR}" ]; then
    echo "[PIPELINE] Error: Strain directory not found at ${STRAIN_DIR}" >&2
    echo "[PIPELINE] Please run 'bash setup.sh --region ${REGION} --year ${YEAR} --strain ${STRAIN}' first." >&2
    exit 1
fi

if [ ! -d "${FASTQ_INPUT_DIR}" ]; then
    echo "[PIPELINE] Error: FASTQ input directory not found at ${FASTQ_INPUT_DIR}" >&2
    exit 1
fi

if [ ! -f "${SAMPLE_LIST_FILE}" ]; then
    echo "[PIPELINE] Error: Sample list file not found at ${SAMPLE_LIST_FILE}" >&2
    exit 1
fi

# Check that all required scripts are executable
for script in 01_trimmomatic.sh 02_bwa_alignment.sh 03_bcftools_snp.sh 04_bcftools_indel.sh 05_consensus_translation_msa.sh 06_run_summary.sh; do
    if [ ! -x "${SCRIPTS_DIR}/pipeline/${script}" ]; then
        echo "[PIPELINE] Error: Script ${SCRIPTS_DIR}/pipeline/${script} is not executable." >&2
        echo "[PIPELINE] Please run: chmod +x ${SCRIPTS_DIR}/pipeline/${script}" >&2
        exit 1
    fi
done

# Create log directory
LOG_DIR="${STRAIN_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Define FastQC output directories
FASTQC_PRE_TRIM_DIR="${FASTQ_INPUT_DIR}/fastqc_pre_trim"
FASTQC_POST_TRIM_DIR="${TRIMMING_OUTPUT_DIR}/fastqc_post_trim"
mkdir -p "${FASTQC_PRE_TRIM_DIR}" "${FASTQC_POST_TRIM_DIR}"

# === Read Sample List ===
mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "${SAMPLE_LIST_FILE}" || true)
TOTAL_SAMPLES=${#SAMPLES[@]}

echo "========================================================================" >&2
echo "[PIPELINE] Starting parallel pipeline for ${TOTAL_SAMPLES} samples" >&2
echo "[PIPELINE] Parallel jobs per step: ${PARALLEL_JOBS}" >&2
echo "[PIPELINE] Sample list: ${SAMPLE_LIST_FILE}" >&2
echo "========================================================================" >&2

# === Step Processing Functions ===

# Function to run a step for a single sample
run_step1_sample() {
    local sampleID=$1
    local SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"
    mkdir -p "${SAMPLE_LOG_DIR}"

    local STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step1_fastqc_pre_done"
    if [[ -f "$STEP_MARKER" ]]; then
        log_skip "[1/8] Initial FastQC - Already completed for ${sampleID}"
        return 0
    fi

    log_sample_step 1 8 "Initial FastQC" "${sampleID}"
    if "${FASTQC_PATH}" -t "${THREADS}" -o "${FASTQC_PRE_TRIM_DIR}" \
        "${FASTQ_INPUT_DIR}/${sampleID}_1.fastq.gz" \
        "${FASTQ_INPUT_DIR}/${sampleID}_2.fastq.gz" \
        > "${SAMPLE_LOG_DIR}/01_fastqc_pre.log" 2>&1; then
        touch "$STEP_MARKER"
        log_success "Step 1 completed for ${sampleID}"
        return 0
    else
        log_error "Step 1 failed for ${sampleID}"
        return 1
    fi
}

run_step2_sample() {
    local sampleID=$1
    local SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"

    local STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step2_trimmomatic_done"
    if [[ -f "$STEP_MARKER" ]]; then
        log_skip "[2/8] Trimmomatic - Already completed for ${sampleID}"
        return 0
    fi

    # Check if previous step was completed
    if [[ ! -f "${PIPELINE_DIR}/sample_${sampleID}_step1_fastqc_pre_done" ]]; then
        log_error "Cannot run Step 2 for ${sampleID}: Step 1 not completed"
        return 1
    fi

    log_sample_step 2 8 "Trimmomatic" "${sampleID}"
    if "${SCRIPTS_DIR}/pipeline/01_trimmomatic.sh" "${sampleID}" \
        2>&1 | tee "${SAMPLE_LOG_DIR}/02_trimmomatic.log" > /dev/null; then
        touch "$STEP_MARKER"
        log_success "Step 2 completed for ${sampleID}"
        return 0
    else
        log_error "Step 2 failed for ${sampleID}"
        return 1
    fi
}

run_step3_sample() {
    local sampleID=$1
    local SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"

    local STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step3_fastqc_post_done"
    if [[ -f "$STEP_MARKER" ]]; then
        log_skip "[3/8] Post-Trimming FastQC - Already completed for ${sampleID}"
        return 0
    fi

    # Check if previous step was completed
    if [[ ! -f "${PIPELINE_DIR}/sample_${sampleID}_step2_trimmomatic_done" ]]; then
        log_error "Cannot run Step 3 for ${sampleID}: Step 2 not completed"
        return 1
    fi

    log_sample_step 3 8 "Post-Trimming FastQC" "${sampleID}"
    if "${FASTQC_PATH}" -t "${THREADS}" -o "${FASTQC_POST_TRIM_DIR}" \
        "${TRIMMING_OUTPUT_DIR}/${sampleID}/${sampleID}_1.trim.fastq.gz" \
        "${TRIMMING_OUTPUT_DIR}/${sampleID}/${sampleID}_2.trim.fastq.gz" \
        > "${SAMPLE_LOG_DIR}/03_fastqc_post.log" 2>&1; then
        touch "$STEP_MARKER"
        log_success "Step 3 completed for ${sampleID}"
        return 0
    else
        log_error "Step 3 failed for ${sampleID}"
        return 1
    fi
}

run_step4_sample() {
    local sampleID=$1
    local SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"

    local STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step4_bwa_alignment_done"
    if [[ -f "$STEP_MARKER" ]]; then
        log_skip "[4/8] BWA Alignment - Already completed for ${sampleID}"
        return 0
    fi

    # Check if previous step was completed
    if [[ ! -f "${PIPELINE_DIR}/sample_${sampleID}_step3_fastqc_post_done" ]]; then
        log_error "Cannot run Step 4 for ${sampleID}: Step 3 not completed"
        return 1
    fi

    log_sample_step 4 8 "BWA Alignment" "${sampleID}"
    if "${SCRIPTS_DIR}/pipeline/02_bwa_alignment.sh" "${sampleID}" "${STRAIN}" \
        > "${SAMPLE_LOG_DIR}/04_bwa_alignment.log" 2>&1; then
        touch "$STEP_MARKER"
        log_success "Step 4 completed for ${sampleID}"
        return 0
    else
        log_error "Step 4 failed for ${sampleID}"
        return 1
    fi
}

run_step5_sample() {
    local sampleID=$1
    local SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"

    local STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step5_bcftools_snp_done"
    if [[ -f "$STEP_MARKER" ]]; then
        log_skip "[5/8] BCFtools SNP Calling - Already completed for ${sampleID}"
        return 0
    fi

    # Check if previous step was completed
    if [[ ! -f "${PIPELINE_DIR}/sample_${sampleID}_step4_bwa_alignment_done" ]]; then
        log_error "Cannot run Step 5 for ${sampleID}: Step 4 not completed"
        return 1
    fi

    log_sample_step 5 8 "BCFtools SNP Calling" "${sampleID}"
    if "${SCRIPTS_DIR}/pipeline/03_bcftools_snp.sh" "${sampleID}" "${STRAIN}" \
        > "${SAMPLE_LOG_DIR}/05_bcftools_snp.log" 2>&1; then
        touch "$STEP_MARKER"
        log_success "Step 5 completed for ${sampleID}"
        return 0
    else
        log_error "Step 5 failed for ${sampleID}"
        return 1
    fi
}

run_step6_sample() {
    local sampleID=$1
    local SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"

    local STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step6_bcftools_indel_done"
    if [[ -f "$STEP_MARKER" ]]; then
        log_skip "[6/8] BCFtools INDEL Calling - Already completed for ${sampleID}"
        return 0
    fi

    # Check if previous step was completed
    if [[ ! -f "${PIPELINE_DIR}/sample_${sampleID}_step5_bcftools_snp_done" ]]; then
        log_error "Cannot run Step 6 for ${sampleID}: Step 5 not completed"
        return 1
    fi

    log_sample_step 6 8 "BCFtools INDEL Calling" "${sampleID}"
    if "${SCRIPTS_DIR}/pipeline/04_bcftools_indel.sh" "${sampleID}" "${STRAIN}" \
        > "${SAMPLE_LOG_DIR}/06_bcftools_indel.log" 2>&1; then
        touch "$STEP_MARKER"
        log_success "Step 6 completed for ${sampleID}"
        return 0
    else
        log_error "Step 6 failed for ${sampleID}"
        return 1
    fi
}

run_step7_sample() {
    local sampleID=$1
    local SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"

    local STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step7_msa_done"
    if [[ -f "$STEP_MARKER" ]]; then
        log_skip "[7/8] Consensus, Translation & MSA - Already completed for ${sampleID}"
        return 0
    fi

    # Check if previous step was completed
    if [[ ! -f "${PIPELINE_DIR}/sample_${sampleID}_step6_bcftools_indel_done" ]]; then
        log_error "Cannot run Step 7 for ${sampleID}: Step 6 not completed"
        return 1
    fi

    log_sample_step 7 8 "Consensus, Translation & MSA" "${sampleID}"
    if "${SCRIPTS_DIR}/pipeline/05_consensus_translation_msa.sh" "${sampleID}" "${STRAIN}" \
        > "${SAMPLE_LOG_DIR}/07_consensus_translation_msa.log" 2>&1; then
        touch "$STEP_MARKER"
        log_success "Step 7 completed for ${sampleID}"
        return 0
    else
        log_error "Step 7 failed for ${sampleID}"
        return 1
    fi
}

run_step8_sample() {
    local sampleID=$1
    local SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"

    local STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step8_summary_done"
    if [[ -f "$STEP_MARKER" ]]; then
        log_skip "[8/8] Generate Summary Report - Already completed for ${sampleID}"
        return 0
    fi

    # Check if previous step was completed
    if [[ ! -f "${PIPELINE_DIR}/sample_${sampleID}_step7_msa_done" ]]; then
        log_error "Cannot run Step 8 for ${sampleID}: Step 7 not completed"
        return 1
    fi

    log_sample_step 8 8 "Generate Summary Report" "${sampleID}"
    if "${SCRIPTS_DIR}/pipeline/06_run_summary.sh" "${sampleID}" \
        > "${SAMPLE_LOG_DIR}/08_summary.log" 2>&1; then
        touch "$STEP_MARKER"
        log_success "Step 8 completed for ${sampleID}"
        return 0
    else
        log_error "Step 8 failed for ${sampleID}"
        return 1
    fi
}

# Export functions and variables for parallel execution
export -f run_step1_sample run_step2_sample run_step3_sample run_step4_sample
export -f run_step5_sample run_step6_sample run_step7_sample run_step8_sample
export -f log_sample_step log_error log_skip log_success
export PIPELINE_DIR LOG_DIR FASTQ_INPUT_DIR FASTQC_PRE_TRIM_DIR FASTQC_POST_TRIM_DIR
export TRIMMING_OUTPUT_DIR SCRIPTS_DIR STRAIN FASTQC_PATH THREADS

# === Main Pipeline Execution ===

# Step 1: Initial FastQC
log_step_start 1 "Initial FastQC (all samples)"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "${PARALLEL_JOBS}" -I {} bash -c 'run_step1_sample "$@"' _ {}

# Step 2: Trimmomatic
log_step_start 2 "Trimmomatic (all samples)"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "${PARALLEL_JOBS}" -I {} bash -c 'run_step2_sample "$@"' _ {}

# Step 3: Post-Trimming FastQC
log_step_start 3 "Post-Trimming FastQC (all samples)"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "${PARALLEL_JOBS}" -I {} bash -c 'run_step3_sample "$@"' _ {}

# Step 4: BWA Alignment
log_step_start 4 "BWA Alignment (all samples)"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "${PARALLEL_JOBS}" -I {} bash -c 'run_step4_sample "$@"' _ {}

# Step 5: BCFtools SNP Calling
log_step_start 5 "BCFtools SNP Calling (all samples)"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "${PARALLEL_JOBS}" -I {} bash -c 'run_step5_sample "$@"' _ {}

# Step 6: BCFtools INDEL Calling
log_step_start 6 "BCFtools INDEL Calling (all samples)"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "${PARALLEL_JOBS}" -I {} bash -c 'run_step6_sample "$@"' _ {}

# Step 7: Consensus, Translation, and MSA
log_step_start 7 "Consensus, Translation & MSA (all samples)"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "${PARALLEL_JOBS}" -I {} bash -c 'run_step7_sample "$@"' _ {}

# Step 8: Generate Summary Report
log_step_start 8 "Generate Summary Report (all samples)"
printf '%s\n' "${SAMPLES[@]}" | xargs -P "${PARALLEL_JOBS}" -I {} bash -c 'run_step8_sample "$@"' _ {}

# === Mark Complete Samples ===
COMPLETED_COUNT=0
FAILED_COUNT=0

for sampleID in "${SAMPLES[@]}"; do
    SAMPLE_MARKER="${PIPELINE_DIR}/sample_${sampleID}_done"

    # Check if all 8 steps are completed
    ALL_STEPS_DONE=true
    for step in {1..8}; do
        STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step${step}_*_done"
        if ! ls $STEP_MARKER 1> /dev/null 2>&1; then
            ALL_STEPS_DONE=false
            break
        fi
    done

    if [[ "$ALL_STEPS_DONE" == "true" ]]; then
        touch "$SAMPLE_MARKER"
        ((COMPLETED_COUNT++))
    else
        ((FAILED_COUNT++))
    fi
done

# === Final Summary Report ===
echo "" >&2
echo "🎉 Pipeline execution completed!" >&2
show_elapsed_time
echo "" >&2
echo "================ Pipeline Execution Summary ================" >&2
echo "Total samples: ${TOTAL_SAMPLES}" >&2
echo "Successfully completed: ${COMPLETED_COUNT}" >&2
echo "Failed or incomplete: ${FAILED_COUNT}" >&2
echo "============================================================" >&2

if [[ $FAILED_COUNT -gt 0 ]]; then
    echo "" >&2
    echo "⚠️  Some samples failed or are incomplete. Check logs in: ${LOG_DIR}" >&2
    exit 1
fi

# === Generate Overall Summary (All Samples) ===
if [[ $COMPLETED_COUNT -gt 0 ]]; then
    echo "" >&2
    echo "========================================================================" >&2
    echo "[PIPELINE] $(date '+%Y-%m-%d %H:%M:%S') - Generating Overall Summary" >&2
    echo "========================================================================" >&2

    SUMMARY_DIR="${STRAIN_DIR}/5_summary"
    mkdir -p "${SUMMARY_DIR}"

    echo "[$(date '+%H:%M:%S')] 📊 Collecting data from all samples..." >&2

    if python3 "${SCRIPTS_DIR}/python/generate_overall_summary.py" \
        "${SAMPLE_LIST_FILE}" \
        "${STRAIN_DIR}" \
        "${SUMMARY_DIR}" \
        > "${LOG_DIR}/overall_summary.log" 2>&1; then
        echo "[$(date '+%H:%M:%S')] ✅ Overall summary generated successfully!" >&2
        echo "" >&2
        echo "📁 Summary files location:" >&2
        echo "   ${SUMMARY_DIR}/overall_summary_formatted.csv" >&2
        echo "   ${SUMMARY_DIR}/overall_summary_detailed.csv" >&2
        echo "   ${SUMMARY_DIR}/overall_summary_statistics.csv" >&2
        echo "   ${SUMMARY_DIR}/overall_read_counts.png" >&2
        echo "   ${SUMMARY_DIR}/overall_retention_percentages.png" >&2
        echo "   ${SUMMARY_DIR}/overall_variant_counts.png" >&2
    else
        echo "[$(date '+%H:%M:%S')] ⚠️  Warning: Overall summary generation failed. Check log: ${LOG_DIR}/overall_summary.log" >&2
    fi
    echo "" >&2
fi

echo "" >&2
echo "✅ All samples processed successfully!" >&2
