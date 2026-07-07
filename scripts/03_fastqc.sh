#!/usr/bin/env bash

set -euo pipefail
source 01_config.sh

export QC_MULTIQC="$PROJECT_DIR/qc/multiqc_reports"

THREADS=3

mkdir -p "$QC_MULTIQC"

echo "=== Running FastQC on raw reads ==="

fastqc \
    "$RAWDATA_DIR"/*.fastq.gz \
    --outdir "$QC_RAW_DIR" \
    --threads "$THREADS"

echo "=== Running MultiQC ==="

multiqc \
    "$QC_RAW_DIR" \
    --outdir "$QC_MULTIQC" \
    --filename raw_multiqc \
    --force
