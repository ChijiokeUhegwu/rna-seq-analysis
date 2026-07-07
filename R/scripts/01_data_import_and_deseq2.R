# R/scripts/01_data_import_and_deseq2.R

# Step 1. Load required packages ----
library(DESeq2)
library(tidyverse)
library(ashr)


# Step 2. Import all .tabular files ----
# List all tabular files in the counts folder
files <- list.files(
  path    = "counts/",
  pattern = "\\.tabular$",
  full.names = TRUE
)

print(files) # confirm that all files were listed

# Read each count file
count_list <- lapply(files, function(f) {
  read.delim(f, header = TRUE)
})

# Merge all count tables by Geneid
count_matrix <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE),
                       count_list)

# Set gene IDs as row names
rownames(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL # remove the redundant Geneid column

# Convert to integer matrix (required by DESeq2)
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

# Quick check
view(count_matrix)
dim(count_matrix)       
head(count_matrix)
colSums(count_matrix)    # library sizes per sample

# Step 3: Load metadata from CSV ----
metadata <- read.csv("data/metadata/SraRunTable.csv",
                     stringsAsFactors = FALSE)

# Confirm it loaded correctly
print(metadata)

# Set Condition as a factor with Control as the reference level
metadata$Condition <- factor(metadata$Condition,
                             levels = c("Control", "FluD4", "FluD8"))

# Match metadata row order to count matrix column order
# This is critical — DESeq2 requires they align exactly
metadata <- metadata[match(colnames(count_matrix), metadata$SampleID), ]

# Set SampleID as rownames (DESeq2 uses this to match colData to countData)
rownames(metadata) <- metadata$SampleID
metadata$SampleID  <- NULL

# Confirm alignment before proceeding
cat("Count matrix columns:\n"); print(colnames(count_matrix))
cat("Metadata rows:\n");        print(rownames(metadata))
cat("Do they match?", all(colnames(count_matrix) == rownames(metadata)), "\n")


# Step 4: Create DESeq2 object ----
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = metadata,
  design    = ~ Condition       # three-level factor: Control, FluD4, FluD8
)

# Set the reference level (Control = baseline for comparison)
dds$Condition <- relevel(dds$Condition, ref = "Control")

# Step 5: Pre-filter low-count genes ----
# Remove genes with fewer than 10 counts across all samples combined
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")

# Step 6: Run DESeq2 ----
dds <- DESeq(dds)

# Step 7: Extract results — all pairwise comparisons ----
# FluD4 vs Control
res_D4 <- results(dds,
                  contrast = c("Condition", "FluD4", "Control"),
                  alpha    = 0.05)

# FluD8 vs Control
res_D8 <- results(dds,
                  contrast = c("Condition", "FluD8", "Control"),
                  alpha    = 0.05)

# FluD8 vs FluD4 (progression of infection)
res_D8vD4 <- results(dds,
                     contrast = c("Condition", "FluD8", "FluD4"),
                     alpha    = 0.05)

cat("=== FluD4 vs Control ===\n");   summary(res_D4)
cat("=== FluD8 vs Control ===\n");   summary(res_D8)
cat("=== FluD8 vs FluD4 ===\n");     summary(res_D8vD4)


# Step 8: LFC shrinkage and export for each comparison ----
resultsNames(dds)    # run this first to see exact coefficient names

shrink_and_export <- function(dds, contrast, filename_prefix) {
  
  # Shrinkage using ashr for contrasts 
  res_s <- lfcShrink(dds,
                     contrast = contrast,
                     type     = "ashr")
  
  # Full results table
  df_all <- as.data.frame(res_s) %>% 
    tibble::rownames_to_column("gene_id") %>% 
    arrange(padj)
  
  write_tsv(
    df_all,
    paste0("R/results/tables/", filename_prefix, "_all_genes.tsv")
  )
  
  # Significant DEGs
  df_sig <- df_all %>% 
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  
  write_tsv(
    df_sig,
    paste0("R/results/tables/", filename_prefix, "_sig_DEGs.tsv")
  )
  
  cat(filename_prefix, "— significant DEGs:", nrow(df_sig),
      "(up:", sum(df_sig$log2FoldChange > 0),
      "| down:", sum(df_sig$log2FoldChange < 0), ")\n")
  
  return(invisible(df_all))
}

# Run for each comparison
shrink_and_export(dds, c("Condition", "FluD4",  "Control"), "FluD4_vs_Control")
shrink_and_export(dds, c("Condition", "FluD8",  "Control"), "FluD8_vs_Control")
shrink_and_export(dds, c("Condition", "FluD8",  "FluD4"),   "FluD8_vs_FluD4")

saveRDS(dds, "R/results/tables/dds_object.rds")
cat("Analysis complete.\n")