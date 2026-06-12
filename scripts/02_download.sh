#!/usr/bin/env bash
# scripts/02_download.sh
set -euo pipefail

# Pull in all folder paths
source 01_config.sh
 
ACCESSIONS="$METADATA_DIR/SRR_Acc_List.txt"
OUTDIR="$RAWDATA_DIR"
THREADS=3
 
mkdir -p $OUTDIR
 
while read SRR; do
    echo '=============================='
    echo 'Downloading: '$SRR
    echo '=============================='
 
    # prefetch downloads the .sra file to cache
    prefetch $SRR
 
    # fasterq-dump converts .sra to FASTQ
    # --split-files: separate R1 and R2 for paired-end data
    # --skip-technical: removes technical reads
    fasterq-dump $SRR \
        --outdir $OUTDIR \
        --split-files \
        --skip-technical \
        --threads $THREADS
 
    # Compress the output files to save disk space
    gzip $OUTDIR/${SRR}_1.fastq
    gzip $OUTDIR/${SRR}_2.fastq
 
    echo 'Done: '$SRR
done < $ACCESSIONS
 
echo 'All downloads complete.'
ls -lh $OUTDIR

