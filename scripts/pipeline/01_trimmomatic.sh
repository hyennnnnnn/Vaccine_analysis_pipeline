#!/usr/bin/env bash
set -euo pipefail

# This script performs quality trimming on paired-end FASTQ files using Trimmomatic.
# It sources a central configuration file for paths and settings.

# Source the configuration file
source "$(dirname "$0")/../config.sh" || { echo "Error: config.sh not found."; exit 1; }

# Check for sampleID argument
if [ -z "$1" ]; then
    echo "Usage: $0 <sampleID>"
    exit 1
fi
sampleID=$1

# Define input and output files based on config variables
Input_R1="${FASTQ_INPUT_DIR}/${sampleID}_1.fastq.gz"
Input_R2="${FASTQ_INPUT_DIR}/${sampleID}_2.fastq.gz"
Sample_Trim_Dir="${TRIMMING_OUTPUT_DIR}/${sampleID}"
Output_R1_Paired="${Sample_Trim_Dir}/${sampleID}_1.trim.fastq.gz"
Output_R1_Unpaired="${Sample_Trim_Dir}/${sampleID}_1.unpaired.fastq.gz"
Output_R2_Paired="${Sample_Trim_Dir}/${sampleID}_2.trim.fastq.gz"
Output_R2_Unpaired="${Sample_Trim_Dir}/${sampleID}_2.unpaired.fastq.gz"
Log_File="${Sample_Trim_Dir}/${sampleID}_trim.log"

# Create sample-specific output directory
mkdir -p "${Sample_Trim_Dir}" || { echo "Error: Failed to create output directory ${Sample_Trim_Dir}"; exit 1; }

echo "Starting Trimming for ${sampleID}..."

java -jar "${TRIMMOMATIC_JAR}" PE \
    -threads "${THREADS}" \
    -phred33 \
    "${Input_R1}" \
    "${Input_R2}" \
    "${Output_R1_Paired}" "${Output_R1_Unpaired}" \
    "${Output_R2_Paired}" "${Output_R2_Unpaired}" \
    ILLUMINACLIP:"${TRIMMOMATIC_ADAPTERS}":2:30:10 \
    SLIDINGWINDOW:4:15 \
    LEADING:10 \
    TRAILING:10 \
    MINLEN:50 \
    > "${Log_File}" 2>&1

if [ $? -eq 0 ]; then
    echo "${sampleID} - Trimming DONE!!"
else
    echo "Error: Trimming failed for ${sampleID}. Check log file: ${Log_File}"
    exit 1
fi
