#!/usr/bin/env bash
# Note: Do not use 'set -euo pipefail' in config files that are sourced

# --- Configuration for Bioinformatics Pipeline ---

# This config file expects STRAIN_DIR to be set as an environment variable.
# STRAIN_DIR is the strain-specific directory (e.g., /path/to/project/GreenCross_flu_southern_2024/A_H5N1)
if [ -z "${STRAIN_DIR}" ]; then
  echo "Error: STRAIN_DIR environment variable is not set." >&2
  exit 1
fi

# --- Tool Paths ---
# These can be overridden by environment variables
# If not set, defaults to command name (assumes tools are in PATH)
export FASTQC_PATH="${FASTQC_PATH:-fastqc}"
export TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-/path/to/trimmomatic.jar}"
export BWA_PATH="${BWA_PATH:-bwa}"
export SAMTOOLS_PATH="${SAMTOOLS_PATH:-samtools}"
export BCFTOOLS_PATH="${BCFTOOLS_PATH:-bcftools}"
export TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS:-/path/to/trimmomatic/adapters/TruSeq3-PE-2.fa}"
export MAFFT_PATH="${MAFFT_PATH:-mafft}"

# --- Directory Structure Relative to STRAIN_DIR ---
# Input FASTQ files
export FASTQ_INPUT_DIR="${STRAIN_DIR}/0_fastq"
# Trimmomatic output
export TRIMMING_OUTPUT_DIR="${STRAIN_DIR}/1_trimming"
# BWA alignment output
export ALIGN_OUTPUT_DIR="${STRAIN_DIR}/2_align"
# Variant call output
export VARIANT_CALL_OUTPUT_DIR="${STRAIN_DIR}/3_variantCall"
# Reference genome directory (reference file itself will be constructed based on strain)
export REFERENCE_DIR="${STRAIN_DIR}/reference"

# --- Global Settings ---
export THREADS=8 # Number of threads to use for multi-threaded tools

# --- Segment Configuration for Translation and MSA ---
# Comma-separated list of segments to analyze (empty = analyze all)
# Common segments: PB2,PB1,PA,HA,NP,NA,M,NS
# Example: export SEGMENTS_TO_ANALYZE="HA,NA" (for selective analysis)
export SEGMENTS_TO_ANALYZE="HA,NA"  # Customize as needed
