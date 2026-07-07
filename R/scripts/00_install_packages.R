# Install packages (run once) ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c(
  "DESeq2",
  "tidyverse",
  "ashr",
  "EnhancedVolcano",
  "pheatmap",
  "RColorBrewer",
  "clusterProfiler",
  "org.Mm.eg.db"
)

BiocManager::install(packages, ask = FALSE)

