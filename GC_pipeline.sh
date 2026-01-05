#!/usr/bin/env bash
set -euo pipefail

# GC Pipeline - Influenza Sequencing Data Analysis Pipeline
# Analyzes influenza sequencing samples with progress tracking and resumability
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

# === Logging Functions ===
log_step() {
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
    echo "[$(date '+%H:%M:%S')] ⏭️  $1 - Sample already processed" >&2
}

log_success() {
    echo "[$(date '+%H:%M:%S')] ✅ $1" >&2
}

log_info() {
    echo "========================================================================" >&2
    echo "[PIPELINE] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
    echo "========================================================================" >&2
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
    echo "Usage: bash GC_pipeline.sh \\" >&2
    echo "  --region <region> \\" >&2
    echo "  --year <year> \\" >&2
    echo "  --strain <strain> \\" >&2
    echo "  --samples <sample_list_file> \\" >&2
    echo "  [--parent-dir <parent_directory>]" >&2
    echo "" >&2
    echo "Note: If --parent-dir is not specified, uses environment variable GC_PARENT_DIR" >&2
    echo "      or defaults to current directory" >&2
    echo "" >&2
    echo "Example:" >&2
    echo "  bash GC_pipeline.sh \\" >&2
    echo "    --region southern \\" >&2
    echo "    --year 2024 \\" >&2
    echo "    --strain A_H3N2 \\" >&2
    echo "    --parent-dir /data/projects \\" >&2
    echo "    --samples /data/projects/GreenCross_flu_southern_2024/A_H3N2/samples.txt" >&2
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

# === Sample Counting ===
TOTAL_SAMPLES=$(grep -c -v '^[[:space:]]*$' "${SAMPLE_LIST_FILE}" || true)
CURRENT_SAMPLE=0
PROCESSED_COUNT=0
SKIPPED_COUNT=0
FAILED_COUNT=0

log_info "Starting pipeline for ${TOTAL_SAMPLES} samples from ${SAMPLE_LIST_FILE}"

# Create log directory
LOG_DIR="${STRAIN_DIR}/logs"
mkdir -p "${LOG_DIR}"

# === Main Loop ===
while IFS= read -r sampleID || [[ -n "$sampleID" ]]; do
    # Skip empty lines
    if [ -z "$sampleID" ]; then
        continue
    fi

    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))

    # Check if sample was already successfully processed
    SAMPLE_MARKER="${PIPELINE_DIR}/sample_${sampleID}_done"
    if [[ -f "$SAMPLE_MARKER" ]]; then
        log_skip "Sample ${sampleID} (${CURRENT_SAMPLE}/${TOTAL_SAMPLES})"
        SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
        continue
    fi

    log_info "Processing Sample: ${sampleID} (${CURRENT_SAMPLE}/${TOTAL_SAMPLES})"

    # Create sample-specific log directory
    SAMPLE_LOG_DIR="${LOG_DIR}/${sampleID}"
    mkdir -p "${SAMPLE_LOG_DIR}"

    # Flag to track if all steps succeeded
    SAMPLE_SUCCESS=true

    # Define FastQC output directories
    FASTQC_PRE_TRIM_DIR="${FASTQ_INPUT_DIR}/fastqc_pre_trim"
    FASTQC_POST_TRIM_DIR="${TRIMMING_OUTPUT_DIR}/fastqc_post_trim"
    mkdir -p "${FASTQC_PRE_TRIM_DIR}" "${FASTQC_POST_TRIM_DIR}"

    # 1. Initial FastQC
    STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step1_fastqc_pre_done"
    if [[ -f "$STEP_MARKER" ]]; then
        echo "[$(date '+%H:%M:%S')] ⏭️  [1/8] Initial FastQC - Already completed for ${sampleID}" >&2
    else
        log_step 1 8 "Initial FastQC" "${sampleID}"
        if ! "${FASTQC_PATH}" -t "${THREADS}" -o "${FASTQC_PRE_TRIM_DIR}" \
            "${FASTQ_INPUT_DIR}/${sampleID}_1.fastq.gz" \
            "${FASTQ_INPUT_DIR}/${sampleID}_2.fastq.gz" \
            > "${SAMPLE_LOG_DIR}/01_fastqc_pre.log" 2>&1; then
            log_error "Initial FastQC failed for ${sampleID}"
            SAMPLE_SUCCESS=false
        else
            touch "$STEP_MARKER"
        fi
    fi

    if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
        # 2. Trimmomatic
        STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step2_trimmomatic_done"
        if [[ -f "$STEP_MARKER" ]]; then
            echo "[$(date '+%H:%M:%S')] ⏭️  [2/8] Trimmomatic - Already completed for ${sampleID}" >&2
        else
            log_step 2 8 "Trimmomatic" "${sampleID}"
            if ! "${SCRIPTS_DIR}/pipeline/01_trimmomatic.sh" "${sampleID}" \
                2>&1 | tee "${SAMPLE_LOG_DIR}/02_trimmomatic.log" > /dev/null; then
                log_error "Trimmomatic failed for ${sampleID}"
                SAMPLE_SUCCESS=false
            else
                touch "$STEP_MARKER"
            fi
        fi
    fi

    if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
        # 3. Post-Trimming FastQC
        STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step3_fastqc_post_done"
        if [[ -f "$STEP_MARKER" ]]; then
            echo "[$(date '+%H:%M:%S')] ⏭️  [3/8] Post-Trimming FastQC - Already completed for ${sampleID}" >&2
        else
            log_step 3 8 "Post-Trimming FastQC" "${sampleID}"
            if ! "${FASTQC_PATH}" -t "${THREADS}" -o "${FASTQC_POST_TRIM_DIR}" \
                "${TRIMMING_OUTPUT_DIR}/${sampleID}/${sampleID}_1.trim.fastq.gz" \
                "${TRIMMING_OUTPUT_DIR}/${sampleID}/${sampleID}_2.trim.fastq.gz" \
                > "${SAMPLE_LOG_DIR}/03_fastqc_post.log" 2>&1; then
                log_error "Post-trimming FastQC failed for ${sampleID}"
                SAMPLE_SUCCESS=false
            else
                touch "$STEP_MARKER"
            fi
        fi
    fi

    if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
        # 4. BWA Alignment
        STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step4_bwa_alignment_done"
        if [[ -f "$STEP_MARKER" ]]; then
            echo "[$(date '+%H:%M:%S')] ⏭️  [4/8] BWA Alignment - Already completed for ${sampleID}" >&2
        else
            log_step 4 8 "BWA Alignment" "${sampleID}"
            if ! "${SCRIPTS_DIR}/pipeline/02_bwa_alignment.sh" "${sampleID}" "${STRAIN}" \
                > "${SAMPLE_LOG_DIR}/04_bwa_alignment.log" 2>&1; then
                log_error "BWA Alignment failed for ${sampleID}"
                SAMPLE_SUCCESS=false
            else
                touch "$STEP_MARKER"
            fi
        fi
    fi

    if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
        # 5. BCFtools SNP Calling
        STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step5_bcftools_snp_done"
        if [[ -f "$STEP_MARKER" ]]; then
            echo "[$(date '+%H:%M:%S')] ⏭️  [5/8] BCFtools SNP Calling - Already completed for ${sampleID}" >&2
        else
            log_step 5 8 "BCFtools SNP Calling" "${sampleID}"
            if ! "${SCRIPTS_DIR}/pipeline/03_bcftools_snp.sh" "${sampleID}" "${STRAIN}" \
                > "${SAMPLE_LOG_DIR}/05_bcftools_snp.log" 2>&1; then
                log_error "BCFtools SNP calling failed for ${sampleID}"
                SAMPLE_SUCCESS=false
            else
                touch "$STEP_MARKER"
            fi
        fi
    fi

    if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
        # 6. BCFtools INDEL Calling
        STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step6_bcftools_indel_done"
        if [[ -f "$STEP_MARKER" ]]; then
            echo "[$(date '+%H:%M:%S')] ⏭️  [6/8] BCFtools INDEL Calling - Already completed for ${sampleID}" >&2
        else
            log_step 6 8 "BCFtools INDEL Calling" "${sampleID}"
            if ! "${SCRIPTS_DIR}/pipeline/04_bcftools_indel.sh" "${sampleID}" "${STRAIN}" \
                > "${SAMPLE_LOG_DIR}/06_bcftools_indel.log" 2>&1; then
                log_error "BCFtools INDEL calling failed for ${sampleID}"
                SAMPLE_SUCCESS=false
            else
                touch "$STEP_MARKER"
            fi
        fi
    fi

    if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
        # 7. Consensus Generation, Translation, and MSA
        STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step7_msa_done"
        if [[ -f "$STEP_MARKER" ]]; then
            echo "[$(date '+%H:%M:%S')] ⏭️  [7/8] Consensus, Translation & MSA - Already completed for ${sampleID}" >&2
        else
            log_step 7 8 "Consensus, Translation & MSA" "${sampleID}"
            if ! "${SCRIPTS_DIR}/pipeline/05_consensus_translation_msa.sh" "${sampleID}" "${STRAIN}" \
                > "${SAMPLE_LOG_DIR}/07_consensus_translation_msa.log" 2>&1; then
                log_error "Consensus, Translation & MSA failed for ${sampleID}"
                SAMPLE_SUCCESS=false
            else
                touch "$STEP_MARKER"
            fi
        fi
    fi

    if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
        # 8. Generate Summary Report
        STEP_MARKER="${PIPELINE_DIR}/sample_${sampleID}_step8_summary_done"
        if [[ -f "$STEP_MARKER" ]]; then
            echo "[$(date '+%H:%M:%S')] ⏭️  [8/8] Generate Summary Report - Already completed for ${sampleID}" >&2
        else
            log_step 8 8 "Generate Summary Report" "${sampleID}"
            if ! "${SCRIPTS_DIR}/pipeline/06_run_summary.sh" "${sampleID}" \
                > "${SAMPLE_LOG_DIR}/08_summary.log" 2>&1; then
                log_error "Summary report generation failed for ${sampleID}"
                SAMPLE_SUCCESS=false
            else
                touch "$STEP_MARKER"
            fi
        fi
    fi

    # Mark sample as complete if all steps succeeded
    if [[ "$SAMPLE_SUCCESS" == "true" ]]; then
        touch "$SAMPLE_MARKER"
        log_success "Successfully completed pipeline for sample: ${sampleID}"
        PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
    else
        log_error "Pipeline failed for sample: ${sampleID}. Check logs at: ${SAMPLE_LOG_DIR}"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi

    echo "" >&2

done < "${SAMPLE_LIST_FILE}"

# === Final Summary Report ===
echo "" >&2
echo "🎉 Pipeline execution completed!" >&2
show_elapsed_time
echo "" >&2
echo "================ Pipeline Execution Summary ================" >&2
echo "Total samples in list: ${TOTAL_SAMPLES}" >&2
echo "Successfully processed: ${PROCESSED_COUNT}" >&2
echo "Skipped (already done): ${SKIPPED_COUNT}" >&2
echo "Failed: ${FAILED_COUNT}" >&2
echo "============================================================" >&2

if [[ $FAILED_COUNT -gt 0 ]]; then
    echo "" >&2
    echo "⚠️  Some samples failed. Check logs in: ${LOG_DIR}" >&2
    exit 1
fi

if [[ $VERBOSE -eq 1 ]]; then
    echo "" >&2
    echo "📋 To disable verbose output, run with: VERBOSE=0 bash GC_pipeline.sh ..." >&2
fi

# === Generate Overall Summary (All Samples) ===
if [[ $PROCESSED_COUNT -gt 0 ]] || [[ $SKIPPED_COUNT -gt 0 ]]; then
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

log_info "All samples processed successfully!"
