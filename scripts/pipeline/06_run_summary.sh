#!/usr/bin/env bash
set -euo pipefail

# Individual sample summary generation has been disabled.
# Overall summary for all samples is generated after all samples are processed.

# 1. 인자값 파싱 및 기본값 설정
if [ -z "$1" ]; then
    echo "Usage: $0 <sampleID>"
    exit 1
fi
sampleID=$1

echo "${sampleID} - Individual summary skipped (Overall summary will be generated after all samples)"
