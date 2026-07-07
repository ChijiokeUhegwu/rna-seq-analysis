#!/usr/bin/env bash
# scripts/06_hisat2_align.sh
# Aligns raw paired-end reads directly — no trimming step
# Raw reads passed QC and require no trimming

set -euo pipefail
source 01_config.sh

THREADS=4
INDEX="$REF_INDEX_DIR/GRCm39"

echo "================================================="
echo "Starting HISAT2 alignment"
echo "Input: raw reads from $RAWDATA_DIR"
echo "Index: $INDEX"
echo "================================================="

for R1 in "$RAWDATA_DIR"/*_1.fastq.gz; do

    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="$RAWDATA_DIR/${SAMPLE}_2.fastq.gz"

    # Verify both files exist before attempting alignment
    if [[ ! -f "$R2" ]]; then
        echo "ERROR: R2 file not found for $SAMPLE — skipping"
        continue
    fi

    echo "-------------------------------"
    echo "Aligning: $SAMPLE"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "-------------------------------"

    # Pipe HISAT2 output directly into samtools sort
    # This avoids writing an intermediate SAM file (saves ~20-40 GB per sample)
    hisat2 \
        -p "$THREADS" \
        -x "$INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        --dta \
        --summary-file "$ALIGN_LOGS_DIR/${SAMPLE}_summary.txt" \
    | samtools sort \
        -@ "$THREADS" \
        -o "$ALIGN_BAM_DIR/${SAMPLE}.sorted.bam"

    samtools index "$ALIGN_BAM_DIR/${SAMPLE}.sorted.bam"

    samtools flagstat \
        "$ALIGN_BAM_DIR/${SAMPLE}.sorted.bam" \
        > "$ALIGN_LOGS_DIR/${SAMPLE}_flagstat.txt"

    echo "Done: $SAMPLE"
    grep "overall alignment rate" "$ALIGN_LOGS_DIR/${SAMPLE}_summary.txt"
    echo ""

done

echo "================================================="
echo "All alignments complete."
echo ""
echo "=== Alignment Rate Summary ==="
grep "overall alignment rate" "$ALIGN_LOGS_DIR"/*_summary.txt
