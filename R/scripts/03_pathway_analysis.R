# R/scripts/03_pathway_analysis.R

# Step 1. Load required packages ----
library(clusterProfiler)
library(org.Mm.eg.db)      # Mouse annotation database
library(enrichplot)
library(tidyverse)


# Step 2. Load significant DEGs ----
sig_d4 <- read.delim("R/results/tables/FluD4_vs_Control_sig_DEGs.tsv")
sig_d8 <- read.delim("R/results/tables/FluD8_vs_Control_sig_DEGs.tsv")

# Pool all unique significant genes
all_sig <- unique(c(sig_d4$gene_id, sig_d8$gene_id))
cat("Unique significant genes for enrichment:", length(all_sig), "\n")

# Convert Ensembl IDs to Entrez IDs 
entrez_map <- bitr(all_sig,
                   fromType = "ENSEMBL",
                   toType   = c("ENTREZID", "SYMBOL"),
                   OrgDb    = org.Mm.eg.db)

cat("IDs mapped successfully:", nrow(entrez_map), "\n")
print(entrez_map)

# ── GO Enrichment (run but note: small gene list may return nothing) ───────
go_bp <- enrichGO(
  gene          = entrez_map$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

if (nrow(as.data.frame(go_bp)) > 0) {
  write.table(as.data.frame(go_bp),
              "R/results/tables/GO_BP_enrichment.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  pdf("R/results/figures/GO_BP_dotplot.pdf", width = 9, height = 7)
  print(dotplot(go_bp, showCategory = 15,
                title = "GO Biological Process — Significant DEGs"))
  dev.off()
  cat("GO enrichment results saved.\n")
} else {
  cat("No significant GO terms — gene list too small. Proceed with GSEA below.\n")
}

# ── GSEA — uses all tested genes ranked by log2FoldChange ─────────────────
# This is more appropriate than ORA when you have few significant DEGs
# Use FluD8 vs Control as it has the most DEGs

res_d8_all <- read.delim("R/results/tables/FluD8_vs_Control_all_genes.tsv")

# Remove NAs and map to Entrez
res_d8_clean <- res_d8_all |>
  filter(!is.na(log2FoldChange), !is.na(padj)) |>
  mutate(ensembl_id = gene_id)

entrez_all <- bitr(res_d8_clean$ensembl_id,
                   fromType = "ENSEMBL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Mm.eg.db)

res_d8_mapped <- res_d8_clean |>
  inner_join(entrez_all, by = c("ensembl_id" = "ENSEMBL"))

# Build named ranked vector: names = Entrez IDs, values = log2FoldChange
gene_ranks <- res_d8_mapped$log2FoldChange
names(gene_ranks) <- res_d8_mapped$ENTREZID

# Sort descending
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Remove duplicates (keep highest LFC per Entrez ID)
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

cat("Genes in ranked list:", length(gene_ranks), "\n")

# Run GSEA on GO Biological Process
gsea_bp <- gseGO(
  geneList      = gene_ranks,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = TRUE
)

if (nrow(as.data.frame(gsea_bp)) > 0) {
  write.table(as.data.frame(gsea_bp),
              "R/results/tables/GSEA_GO_BP_FluD8vsControl.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  pdf("R/results/figures/GSEA_dotplot.pdf", width = 10, height = 8)
  print(dotplot(gsea_bp, showCategory = 15, split = ".sign") +
          facet_grid(. ~ .sign) +
          ggtitle("GSEA — GO Biological Process",
                  subtitle = "FluD8 vs Control"))
  dev.off()
  cat("GSEA results saved.\n")
} else {
  cat("No significant GSEA terms at padj < 0.05.\n")
}

# Run GSEA on KEGG
gsea_kegg <- gseKEGG(
  geneList      = gene_ranks,
  organism      = "mmu",          # mmu = Mus musculus
  minGSSize     = 10,
  pvalueCutoff  = 0.05,
  verbose       = TRUE
)

if (nrow(as.data.frame(gsea_kegg)) > 0) {
  write.table(as.data.frame(gsea_kegg),
              "R/results/tables/GSEA_KEGG_FluD8vsControl.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  pdf("R/results/figures/GSEA_KEGG_dotplot.pdf", width = 10, height = 7)
  print(dotplot(gsea_kegg, showCategory = 15, split = ".sign") +
          facet_grid(. ~ .sign) +
          ggtitle("GSEA — KEGG Pathways", subtitle = "FluD8 vs Control"))
  dev.off()
  cat("GSEA KEGG results saved.\n")
} else {
  cat("No significant KEGG pathways at padj < 0.05.\n")
}