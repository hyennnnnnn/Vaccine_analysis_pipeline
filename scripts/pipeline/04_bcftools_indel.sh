#!/usr/bin/env bash
set -euo pipefail

# This script calls INDELs from a sorted BAM file using BCFtools.
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

echo "Starting INDEL calling for ${sampleID}..."

# BCFtools mpileup and call for INDELs (same style as SNP script)
"${BCFTOOLS_PATH}" mpileup \
    -q30 -Ou \
    -f "${Reference}" \
    "${Sorted_BAM}" \
| "${BCFTOOLS_PATH}" call -V snps -cv -Ob \
    -o "${sampleID}.aln.indel.bcf"

# Convert BCF to VCF (no filter for INDELs)
"${BCFTOOLS_PATH}" view "${sampleID}.aln.indel.bcf" \
> "${sampleID}.aln.indel.vcf"

# Convert BCF to VCF
"${BCFTOOLS_PATH}" view "${sampleID}.aln.indel.bcf" > "${sampleID}.aln.indel.vcf" || { echo "BCF to VCF conversion failed"; exit 1; }

echo "${sampleID} - BCFtools INDEL calling is DONE!!!"
