#!/usr/bin/env bash

set -euo pipefail
source 01_config.sh

# Download mouse reference genome (GRCm39) from Ensembl
# Mus musculus reference for RNA-seq alignment

ENSEMBL_RELEASE=116
SPECIES="mus_musculus"


echo "=========================================="
echo "Downloading Mouse Reference Genome"
echo "Species: Mus musculus (GRCm39)"
echo "Ensembl Release: $ENSEMBL_RELEASE"
echo "=========================================="

# --- Primary genome assembly (chromosome sequences) ---
wget -P "$REF_GENOME_DIR" \
  ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${SPECIES}/dna/\
Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# --- Decompress genome ---
gunzip -f "$REF_GENOME_DIR"/*.fa.gz

# --- GTF annotation file ---
wget -P "$REF_ANNOTATION_DIR" \
  ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${SPECIES}/\
Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz

# --- Decompress annotation ---
gunzip -f "$REF_ANNOTATION_DIR"/*.gtf.gz

echo "=========================================="
echo "Mouse reference download complete."
echo "Files saved in:"
echo "  $REF_GENOME_DIR"
echo "  $REF_ANNOTATION_DIR"
echo "=========================================="

ls -lh "$REF_GENOME_DIR" "$REF_ANNOTATION_DIR"
