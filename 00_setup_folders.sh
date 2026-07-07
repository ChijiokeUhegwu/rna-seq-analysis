#!/bin/bash
# Set up project folders inside project directory
 
mkdir -p data/raw data/trimmed data/metadata
mkdir -p reference/genome reference/annotation reference/index
mkdir -p qc/raw_fastqc qc/trimmed_fastqc
mkdir -p alignment/sam alignment/bam alignment/logs
mkdir -p counts
mkdir -p scripts
mkdir -p R/scripts R/results/tables R/results/figures R/reports
mkdir -p environment docs
 
# Add .gitkeep to empty folders so git tracks them
find . -type d -empty -exec touch {}/.gitkeep \;
 
echo 'Folder structure created successfully.'
ls -R | head -60

