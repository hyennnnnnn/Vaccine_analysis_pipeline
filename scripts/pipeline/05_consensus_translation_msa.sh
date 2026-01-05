#!/usr/bin/env bash
set -euo pipefail

# This script generates consensus sequences, translates to protein, and performs MSA.
# It sources a central configuration file for paths and settings.

# Source the configuration file
source "$(dirname "$0")/../config.sh" || { echo "Error: config.sh not found."; exit 1; }

# Check for arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <sampleID> <strain>"
    exit 1
fi
sampleID=$1
strain=$2

# Construct path to the reference genome by dynamically finding the file matching the strain
Reference=$(find "${REFERENCE_DIR}" -maxdepth 1 -type f -name "*_${strain}_*8segments_spikein_ATGC.fasta" | head -n 1)

if [ -z "${Reference}" ] || [ ! -f "${Reference}" ]; then
    echo "Error: Reference genome for strain '${strain}' not found in ${REFERENCE_DIR}."
    echo "Please ensure exactly one reference file matching '*_${strain}_*8segments_spikein_ATGC.fasta' exists."
    exit 1
fi

# Define input files
Sample_VC_Dir="${VARIANT_CALL_OUTPUT_DIR}/${sampleID}"
SNP_VCF="${Sample_VC_Dir}/${sampleID}.aln.snp.vcf"
INDEL_VCF="${Sample_VC_Dir}/${sampleID}.aln.indel.vcf"

# Define output directory
MSA_OUTPUT_DIR="${STRAIN_DIR}/4_msa/${sampleID}"
mkdir -p "${MSA_OUTPUT_DIR}" || { echo "Error: Failed to create output directory ${MSA_OUTPUT_DIR}"; exit 1; }

echo "Starting consensus generation, translation, and MSA for ${sampleID}..."

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Call the Python script
python3 "${SCRIPT_DIR}/../python/generate_consensus_translate_msa.py" \
    "${sampleID}" \
    "${Reference}" \
    "${SNP_VCF}" \
    "${INDEL_VCF}" \
    "${BCFTOOLS_PATH}" \
    "${SAMTOOLS_PATH}" \
    "${MAFFT_PATH}" \
    "${SEGMENTS_TO_ANALYZE}" \
    "${MSA_OUTPUT_DIR}" || {
    echo "Error: Python MSA script failed."
    exit 1
}

echo "${sampleID} - Consensus generation, translation, and MSA DONE!!!"
