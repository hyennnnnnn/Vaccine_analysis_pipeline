#!/usr/bin/env bash
set -euo pipefail

# GC Pipeline - Project Setup Script
# Creates project directory structure, copies data, and prepares for analysis

# Trap errors
trap 'echo "Error at line $LINENO" >&2' ERR

# === Parse Arguments ===
REGION=""
YEAR=""
RAW_DATA_DIR=""
REFERENCE_FILE=""
STRAIN=""
PARENT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --region)
            REGION="$2"
            shift 2
            ;;
        --year)
            YEAR="$2"
            shift 2
            ;;
        --raw-data)
            RAW_DATA_DIR="$2"
            shift 2
            ;;
        --reference)
            REFERENCE_FILE="$2"
            shift 2
            ;;
        --strain)
            STRAIN="$2"
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

# === Configuration ===
# Use environment variable GC_PARENT_DIR if --parent-dir not specified
# Default to current directory if neither provided
PARENT_DIR="${PARENT_DIR:-${GC_PARENT_DIR:-$(pwd)}}"

# === Usage Check ===
if [[ -z "$REGION" || -z "$YEAR" || -z "$STRAIN" ]]; then
    echo "Usage: bash setup.sh \\" >&2
    echo "  --region <region> \\" >&2
    echo "  --year <year> \\" >&2
    echo "  --strain <strain_name> \\" >&2
    echo "  [--parent-dir <parent_directory>] \\" >&2
    echo "  [--raw-data <raw_data_directory>] \\" >&2
    echo "  [--reference <reference_fasta_file>]" >&2
    echo "" >&2
    echo "Note: If --parent-dir is not specified, uses environment variable GC_PARENT_DIR" >&2
    echo "      or defaults to current directory" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "  # Minimal setup (manual data preparation required)" >&2
    echo "  bash setup.sh --region southern --year 2024 --strain A_H3N2 --parent-dir /data/projects" >&2
    echo "" >&2
    echo "  # Full automated setup" >&2
    echo "  bash setup.sh \\" >&2
    echo "    --region avianflu \\" >&2
    echo "    --year 2025 \\" >&2
    echo "    --strain A_H5N1 \\" >&2
    echo "    --parent-dir /data/projects \\" >&2
    echo "    --raw-data /path/to/raw_data \\" >&2
    echo "    --reference /path/to/2025_A_H5N1_8segments_spikein_ATGC.fasta" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR/scripts"

# === Setup Project Directory Structure ===
export BASE_PROJECT_DIR="${PARENT_DIR}/GreenCross_flu_${REGION}_${YEAR}"
export STRAIN_DIR="${BASE_PROJECT_DIR}/${STRAIN}"
PIPELINE_DIR="${STRAIN_DIR}/.pipeline"
mkdir -p "$PIPELINE_DIR"

SETUP_MARKER="$PIPELINE_DIR/setup_${REGION}_${YEAR}_${STRAIN}_done"

if [[ -f "$SETUP_MARKER" ]]; then
    echo "[SETUP] Project already set up at: ${STRAIN_DIR}" >&2
    echo "[SETUP] To re-run setup, delete marker: rm $SETUP_MARKER" >&2
    exit 0
fi

echo "========================================================================" >&2
echo "[SETUP] GC Pipeline - Project Setup" >&2
echo "========================================================================" >&2
echo "Region: ${REGION}" >&2
echo "Year: ${YEAR}" >&2
echo "Strain: ${STRAIN}" >&2
echo "Project directory: ${BASE_PROJECT_DIR}" >&2
echo "Strain directory: ${STRAIN_DIR}" >&2
echo "========================================================================" >&2
echo "" >&2

# === Source Configuration ===
if [[ ! -f "${SCRIPTS_DIR}/config.sh" ]]; then
    echo "Error: config.sh not found at ${SCRIPTS_DIR}/config.sh" >&2
    exit 1
fi
source "${SCRIPTS_DIR}/config.sh"

# === Validate Required Tools ===
echo "[SETUP] Step 1/6: Validating required tools..." >&2

validate_tool() {
    local tool_path=$1
    local tool_name=$2
    if [[ ! -x "$tool_path" ]]; then
        echo "  ❌ Error: $tool_name not found or not executable at: $tool_path" >&2
        return 1
    fi
    echo "  ✓ $tool_name" >&2
}

validate_tool "$BWA_PATH" "BWA" || exit 1
validate_tool "$SAMTOOLS_PATH" "SAMtools" || exit 1
validate_tool "$BCFTOOLS_PATH" "BCFtools" || exit 1

if [[ ! -f "$TRIMMOMATIC_JAR" ]]; then
    echo "  ❌ Error: Trimmomatic JAR not found at: $TRIMMOMATIC_JAR" >&2
    exit 1
fi
echo "  ✓ Trimmomatic" >&2

if [[ ! -f "$TRIMMOMATIC_ADAPTERS" ]]; then
    echo "  ❌ Error: Trimmomatic adapters not found at: $TRIMMOMATIC_ADAPTERS" >&2
    exit 1
fi
echo "  ✓ Trimmomatic adapters" >&2
echo "" >&2

# === Create Directory Structure ===
echo "[SETUP] Step 2/6: Creating project directory structure..." >&2
# Create base project directory
mkdir -p "${BASE_PROJECT_DIR}"
# Create strain-specific directory structure
"${SCRIPTS_DIR}/utils/setup_directories.sh" "${STRAIN_DIR}"
echo "  ✓ Directory structure created at ${STRAIN_DIR}" >&2
echo "" >&2

# === Copy Raw FASTQ Files ===
if [[ -n "$RAW_DATA_DIR" ]]; then
    echo "[SETUP] Step 3/6: Copying FASTQ files from ${RAW_DATA_DIR}..." >&2

    if [[ ! -d "$RAW_DATA_DIR" ]]; then
        echo "  ❌ Error: Raw data directory not found: $RAW_DATA_DIR" >&2
        exit 1
    fi

    # Count FASTQ files
    FASTQ_COUNT=$(find "$RAW_DATA_DIR" -maxdepth 1 -name "*.fastq.gz" | wc -l | tr -d ' ')

    if [[ $FASTQ_COUNT -eq 0 ]]; then
        echo "  ❌ Error: No .fastq.gz files found in $RAW_DATA_DIR" >&2
        exit 1
    fi

    echo "  Found ${FASTQ_COUNT} FASTQ files" >&2
    cp -v "$RAW_DATA_DIR"/*.fastq.gz "${STRAIN_DIR}/0_fastq/" >&2
    echo "  ✓ FASTQ files copied" >&2
    echo "" >&2

    # === Generate Sample List ===
    echo "[SETUP] Step 4/6: Generating sample list..." >&2

    SAMPLE_LIST="${STRAIN_DIR}/samples.txt"

    # Extract unique sample prefixes (before _1.fastq.gz or _2.fastq.gz)
    find "${STRAIN_DIR}/0_fastq" -name "*_1.fastq.gz" | \
        xargs -n1 basename | \
        sed 's/_1\.fastq\.gz$//' | \
        sort -u > "$SAMPLE_LIST"

    SAMPLE_COUNT=$(wc -l < "$SAMPLE_LIST" | tr -d ' ')

    if [[ $SAMPLE_COUNT -eq 0 ]]; then
        echo "  ⚠️  Warning: No paired-end samples found (expected *_1.fastq.gz and *_2.fastq.gz)" >&2
        rm "$SAMPLE_LIST"
    else
        echo "  ✓ Sample list created: ${SAMPLE_LIST}" >&2
        echo "  ✓ Detected ${SAMPLE_COUNT} samples:" >&2
        cat "$SAMPLE_LIST" | sed 's/^/    - /' >&2
    fi
    echo "" >&2
else
    echo "[SETUP] Step 3/6: Skipping FASTQ copy (--raw-data not provided)" >&2
    echo "" >&2
    echo "[SETUP] Step 4/6: Skipping sample list generation (--raw-data not provided)" >&2
    echo "" >&2
fi

# === Copy and Index Reference Genome ===
if [[ -n "$REFERENCE_FILE" ]]; then
    echo "[SETUP] Step 5/6: Copying reference genome..." >&2

    if [[ ! -f "$REFERENCE_FILE" ]]; then
        echo "  ❌ Error: Reference file not found: $REFERENCE_FILE" >&2
        exit 1
    fi

    REFERENCE_NAME=$(basename "$REFERENCE_FILE")
    REFERENCE_DEST="${STRAIN_DIR}/reference/${REFERENCE_NAME}"

    cp -v "$REFERENCE_FILE" "$REFERENCE_DEST" >&2
    echo "  ✓ Reference genome copied: ${REFERENCE_NAME}" >&2
    echo "" >&2

    # === Index Reference Genome ===
    echo "[SETUP] Step 6/6: Indexing reference genome with BWA..." >&2
    echo "  This may take several minutes..." >&2

    cd "${STRAIN_DIR}/reference"
    "${BWA_PATH}" index "$REFERENCE_NAME" > bwa_index.log 2>&1

    if [[ $? -eq 0 ]]; then
        echo "  ✓ BWA indexing completed" >&2
        echo "  ✓ Index files created:" >&2
        ls -1 "${REFERENCE_NAME}".* | sed 's/^/    - /' >&2
    else
        echo "  ❌ Error: BWA indexing failed. Check: ${STRAIN_DIR}/reference/bwa_index.log" >&2
        exit 1
    fi

    cd "$SCRIPT_DIR"
    echo "" >&2
else
    echo "[SETUP] Step 5/6: Skipping reference copy (--reference not provided)" >&2
    echo "" >&2
    echo "[SETUP] Step 6/6: Skipping BWA indexing (--reference not provided)" >&2
    echo "" >&2
fi

# === Mark Setup as Complete ===
touch "$SETUP_MARKER"

# === Display Summary and Next Steps ===
echo "========================================================================" >&2
echo "✅ Setup completed successfully!" >&2
echo "========================================================================" >&2
echo "" >&2

if [[ -n "$RAW_DATA_DIR" && -n "$REFERENCE_FILE" ]]; then
    # Full automated setup
    echo "🎉 Full automated setup complete! Ready to run pipeline." >&2
    echo "" >&2
    echo "Next step:" >&2
    echo "" >&2
    echo "  bash GC_pipeline.sh \\" >&2
    echo "    --region ${REGION} \\" >&2
    echo "    --year ${YEAR} \\" >&2
    echo "    --strain ${STRAIN} \\" >&2
    echo "    --samples ${STRAIN_DIR}/samples.txt" >&2
    echo "" >&2

elif [[ -n "$RAW_DATA_DIR" || -n "$REFERENCE_FILE" ]]; then
    # Partial automated setup
    echo "⚠️  Partial setup complete. Manual steps required:" >&2
    echo "" >&2

    if [[ -z "$RAW_DATA_DIR" ]]; then
        echo "1. Copy FASTQ files:" >&2
        echo "   cp /path/to/sample*.fastq.gz ${STRAIN_DIR}/0_fastq/" >&2
        echo "" >&2
        echo "2. Create sample list:" >&2
        echo "   ls ${STRAIN_DIR}/0_fastq/*_1.fastq.gz | \\" >&2
        echo "     xargs -n1 basename | sed 's/_1\.fastq\.gz$//' > ${STRAIN_DIR}/samples.txt" >&2
        echo "" >&2
    fi

    if [[ -z "$REFERENCE_FILE" ]]; then
        echo "3. Copy reference genome:" >&2
        echo "   cp /path/to/YYYY_STRAIN_8segments_spikein_ATGC.fasta ${STRAIN_DIR}/reference/" >&2
        echo "" >&2
        echo "4. Index reference genome:" >&2
        echo "   ${BWA_PATH} index ${STRAIN_DIR}/reference/YOUR_REFERENCE.fasta" >&2
        echo "" >&2
    fi

    echo "5. Run pipeline:" >&2
    echo "   bash GC_pipeline.sh \\" >&2
    echo "     --region ${REGION} \\" >&2
    echo "     --year ${YEAR} \\" >&2
    echo "     --strain ${STRAIN} \\" >&2
    echo "     --samples ${STRAIN_DIR}/samples.txt" >&2
    echo "" >&2

else
    # Minimal setup - all manual
    echo "📋 Manual data preparation required:" >&2
    echo "" >&2
    echo "1. Copy FASTQ files:" >&2
    echo "   cp /path/to/sample*.fastq.gz ${STRAIN_DIR}/0_fastq/" >&2
    echo "" >&2
    echo "2. Copy reference genome:" >&2
    echo "   cp /path/to/YYYY_${STRAIN}_8segments_spikein_ATGC.fasta ${STRAIN_DIR}/reference/" >&2
    echo "" >&2
    echo "3. Index reference genome:" >&2
    echo "   ${BWA_PATH} index ${STRAIN_DIR}/reference/YYYY_${STRAIN}_8segments_spikein_ATGC.fasta" >&2
    echo "" >&2
    echo "4. Create sample list:" >&2
    echo "   ls ${STRAIN_DIR}/0_fastq/*_1.fastq.gz | \\" >&2
    echo "     xargs -n1 basename | sed 's/_1\.fastq\.gz$//' > ${STRAIN_DIR}/samples.txt" >&2
    echo "" >&2
    echo "5. Run pipeline:" >&2
    echo "   bash GC_pipeline.sh \\" >&2
    echo "     --region ${REGION} \\" >&2
    echo "     --year ${YEAR} \\" >&2
    echo "     --strain ${STRAIN} \\" >&2
    echo "     --samples ${STRAIN_DIR}/samples.txt" >&2
    echo "" >&2
fi

echo "========================================================================" >&2
