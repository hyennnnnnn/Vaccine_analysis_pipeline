#!/usr/bin/env bash
set -euo pipefail

# This script calls SNPs from a sorted BAM file using BCFtools.
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
    echo "Error: Reference genome for strain '${strain}' not found or multiple found in ${REFERENCE_DIR}."
    echo "Please ensure exactly one reference file matching '*_${strain}_*8segments_spikein_ATGC.fasta' exists."
    exit 1
fi

# Define input and output directories
Sample_Align_Dir="${ALIGN_OUTPUT_DIR}/${sampleID}"
Sample_VC_Dir="${VARIANT_CALL_OUTPUT_DIR}/${sampleID}"
Sorted_BAM="${Sample_Align_Dir}/${sampleID}.aln.sort.bam"

# Create sample-specific output directory
mkdir -p "${Sample_VC_Dir}" || { echo "Error: Failed to create output directory ${Sample_VC_Dir}"; exit 1; }

# Change to the sample's variant calling directory
cd "${Sample_VC_Dir}" || exit 1

echo "Starting SNP calling for ${sampleID}..."

# BCFtools mpileup and call (EXACT same parameters as original script)
"${BCFTOOLS_PATH}" mpileup \
    -q30 -Ou \
    --max-depth 1000000 \
    -f "${Reference}" \
    "${Sorted_BAM}" \
| "${BCFTOOLS_PATH}" call -V indels -cv -Ob \
    -o "${sampleID}.aln.snp.bcf"

# Filter
"${BCFTOOLS_PATH}" view -i 'DP>=10000' "${sampleID}.aln.snp.bcf" > "${sampleID}.aln.snp.vcf" || {
    echo "Error: BCFtools filtering failed."
    exit 1
}

echo "${sampleID} - BCFtools SNP calling is DONE!!!"
