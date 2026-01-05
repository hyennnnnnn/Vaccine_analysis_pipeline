#!/usr/bin/env bash
set -euo pipefail

# This script creates the standard directory structure for a new project.

# Check for argument
if [ -z "$1" ]; then
    echo "Usage: $0 <base_project_directory_path>"
    echo "  <base_project_directory_path>: The full path of the project directory to create."
    exit 1
fi

BASE_DIR=$1

echo "Creating directory structure in: ${BASE_DIR}"

# Create the main directories
mkdir -p "${BASE_DIR}/0_fastq/fastqc_pre_trim"
mkdir -p "${BASE_DIR}/1_trimming/fastqc_post_trim"
mkdir -p "${BASE_DIR}/2_align/flagstat"
mkdir -p "${BASE_DIR}/2_align/idxstat"
mkdir -p "${BASE_DIR}/3_variantCall"
mkdir -p "${BASE_DIR}/4_msa"
mkdir -p "${BASE_DIR}/5_summary"
mkdir -p "${BASE_DIR}/reference"

# Check if creation was successful
if [ -d "${BASE_DIR}/3_variantCall" ]; then
    echo "Directory structure created successfully."
else
    echo "Error: Failed to create directory structure."
    exit 1
fi
