#!/usr/bin/env bash

# Define the main project directory as an environment variable
export PROJECT_DIR="/home/chijioke/rna-seq-analysis"

# 1. Data Directories
export RAWDATA_DIR="$PROJECT_DIR/data/raw"
export TRIMMED_DATA_DIR="$PROJECT_DIR/data/trimmed"
export METADATA_DIR="$PROJECT_DIR/data/metadata"

# 2. Reference Genome Directories
export REF_GENOME_DIR="$PROJECT_DIR/reference/genome"
export REF_ANNOTATION_DIR="$PROJECT_DIR/reference/annotation"
export REF_INDEX_DIR="$PROJECT_DIR/reference/index"

# 3. Quality Control (QC) Reports
export QC_RAW_DIR="$PROJECT_DIR/qc/raw_fastqc"
export QC_TRIMMED_DIR="$PROJECT_DIR/qc/trimmed_fastqc"
export QC_MULTIQC="$PROJECT_DIR/qc/multiqc_reports"

# 4. Alignment Outputs
export ALIGN_SAM_DIR="$PROJECT_DIR/alignment/sam"
export ALIGN_BAM_DIR="$PROJECT_DIR/alignment/bam"
export ALIGN_LOGS_DIR="$PROJECT_DIR/alignment/logs"

# 5. Quantification (Gene Counts)
export COUNTS_DIR="$PROJECT_DIR/counts"

# 6. Downstream R Analysis & Scripts
export SCRIPTS_DIR="$PROJECT_DIR/scripts"
export R_SCRIPTS_DIR="$PROJECT_DIR/R/scripts"
export R_TABLES_DIR="$PROJECT_DIR/R/results/tables"
export R_FIGURES_DIR="$PROJECT_DIR/R/results/figures"
export R_REPORTS_DIR="$PROJECT_DIR/R/reports"

# 7. Miscellaneous Documentation & Conda Environments
export ENV_DIR="$PROJECT_DIR/environment"
export DOCS_DIR="$PROJECT_DIR/docs"

# Automatically create all directories structurally if they don't exist
mkdir -p "$RAWDATA_DIR" "$TRIMMED_DATA_DIR" "$METADATA_DIR" \
         "$REF_GENOME_DIR" "$REF_ANNOTATION_DIR" "$REF_INDEX_DIR" \
         "$QC_RAW_DIR" "$QC_TRIMMED_DIR" \
         "$ALIGN_SAM_DIR" "$ALIGN_BAM_DIR" "$ALIGN_LOGS_DIR" \
         "$COUNTS_DIR" "$SCRIPTS_DIR" \
         "$R_SCRIPTS_DIR" "$R_TABLES_DIR" "$R_FIGURES_DIR" "$R_REPORTS_DIR" \
         "$ENV_DIR" "$DOCS_DIR"

echo "========================================================="
echo " Environment variables configured for RNA-Seq Pipeline"
echo " Project root: $PROJECT_DIR"
echo "========================================================="
