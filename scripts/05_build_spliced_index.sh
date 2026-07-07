#!/usr/bin/env bash
# scripts/05_build_index.sh

set -euo pipefail
source 01_config.sh
 
GENOME_FA="$REF_GENOME_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa"
ANNOT_GTF="$REF_ANNOTATION_DIR/Mus_musculus.GRCm39.116.gtf"
INDEX_PREFIX="$REF_INDEX_DIR/GRCm39"
THREADS=3
 
# Extract splice sites and exons from GTF for splice-aware alignment
echo 'Extracting splice sites...'
hisat2_extract_splice_sites.py $ANNOT_GTF > $REF_INDEX_DIR/splicesites.txt
 
echo 'Extracting exon information...'
hisat2_extract_exons.py $ANNOT_GTF > $REF_INDEX_DIR/exons.txt
 
echo 'Building HISAT2 index (this may take 30-90 minutes)...'
hisat2-build \
    -p $THREADS \
    --ss $REF_INDEX_DIR/splicesites.txt \
    --exon $REF_INDEX_DIR/exons.txt \
    $GENOME_FA \
    $INDEX_PREFIX
 
echo 'HISAT2 index complete.'
ls -lh $REF_INDEX_DIR

