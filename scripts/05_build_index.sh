#!/usr/bin/env bash
# scripts/05_build_index.sh
# Simple HISAT2 index — no splice-site extraction
# Perfectly valid for gene-level count matrix generation with featureCounts

set -euo pipefail
source 01_config.sh

GENOME_FA="$REF_GENOME_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa"
INDEX_PREFIX="$REF_INDEX_DIR/GRCm39"
THREADS=4   

echo "================================================="
echo "Building HISAT2 genome index "
echo "================================================="

mkdir -p "$REF_INDEX_DIR"

hisat2-build \
    -p "$THREADS" \
    "$GENOME_FA" \
    "$INDEX_PREFIX"

echo "================================================="
echo "Index complete."
ls -lh "$REF_INDEX_DIR"
