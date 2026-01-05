#!/usr/bin/env bash
set -euo pipefail

# This script aligns trimmed FASTQ files to a reference genome using BWA.
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
Sample_Trim_Dir="${TRIMMING_OUTPUT_DIR}/${sampleID}"
Sample_Align_Dir="${ALIGN_OUTPUT_DIR}/${sampleID}"

# Create sample-specific output directories
mkdir -p "${Sample_Align_Dir}/logs" || { echo "Error: Failed to create output directory ${Sample_Align_Dir}/logs"; exit 1; }

# Change to the sample's alignment directory to keep outputs organized
cd "${Sample_Align_Dir}" || exit 1

echo "Starting Alignment for ${sampleID}..."

# Define input files
Trimmed_R1="${Sample_Trim_Dir}/${sampleID}_1.trim.fastq.gz"
Trimmed_R2="${Sample_Trim_Dir}/${sampleID}_2.trim.fastq.gz"

# 1. mapping (global alignment)
echo "Step 1: Mapping with BWA..."
# 1st aln
"${BWA_PATH}" aln -t "${THREADS}" "${Reference}" "${Trimmed_R1}" -f "${sampleID}_1.sai" > "logs/${sampleID}_1.sai.log" 2>&1 || { echo "BWA aln R1 failed"; exit 1; }
# 2nd aln
"${BWA_PATH}" aln -t "${THREADS}" "${Reference}" "${Trimmed_R2}" -f "${sampleID}_2.sai" > "logs/${sampleID}_2.sai.log" 2>&1 || { echo "BWA aln R2 failed"; exit 1; }
# bwa sampe
"${BWA_PATH}" sampe "${Reference}" "${sampleID}_1.sai" "${sampleID}_2.sai" "${Trimmed_R1}" "${Trimmed_R2}" -f "${sampleID}.aln.sam" > "logs/${sampleID}.sam.log" 2>&1 || { echo "BWA sampe failed"; exit 1; }
echo "Done mapping, ${sampleID}"

# 2. transform SAM to BAM
echo "Step 2: Converting SAM to BAM..."
"${SAMTOOLS_PATH}" view -bS "${sampleID}.aln.sam" -o "${sampleID}.aln.bam" > "logs/${sampleID}.bam.log" 2>&1 || { echo "SAM to BAM failed"; exit 1; }
echo "Done transforming(make BAM), ${sampleID}"

# 3. sorting BAM file
echo "Step 3: Sorting BAM file..."
"${SAMTOOLS_PATH}" sort -@ "${THREADS}" "${sampleID}.aln.bam" -o "${sampleID}.aln.sort.bam" > "logs/${sampleID}.sort.log" 2>&1 || { echo "BAM sorting failed"; exit 1; }
echo "Done sorting, ${sampleID}"

# 4. indexing BAM file
echo "Step 4: Indexing BAM file..."
"${SAMTOOLS_PATH}" index "${sampleID}.aln.sort.bam" > "logs/${sampleID}.index.log" 2>&1 || { echo "BAM indexing failed"; exit 1; }
echo "Done indexing, ${sampleID}"

# 5. making flagstat file
echo "Step 5: Generating flagstat..."
mkdir -p ../flagstat
"${SAMTOOLS_PATH}" flagstat "${sampleID}.aln.sort.bam" > "${sampleID}.aln.flagstat" || { echo "Flagstat generation failed"; exit 1; }
cp "${sampleID}.aln.flagstat" ../flagstat/

# 6. making idxstat file
echo "Step 6: Generating idxstats..."
mkdir -p ../idxstat
"${SAMTOOLS_PATH}" idxstats "${sampleID}.aln.sort.bam" > "${sampleID}.aln.idxstat" || { echo "Idxstats generation failed"; exit 1; }
cp "${sampleID}.aln.idxstat" ../idxstat/

echo "${sampleID} - Global Alignment is Perfectly Done!!"
