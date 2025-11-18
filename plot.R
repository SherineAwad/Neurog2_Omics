#!/usr/bin/env Rscript

# ======================================
# Full functional heatmap with Neurog2-9SA
# ======================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ----------------------------
# User input
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Usage: Rscript plot_functional_heatmap.R <DE_csv_file>")
}

de_csv <- args[1]        # CSV file with DE results
output_png <- "Neurog2-S9A_heatmap.png"

# ----------------------------
# Define gene sets using exact CSV names
# ----------------------------
gene_dict <- list(
  Transgene = "Neurog2-9SA",
  MG_markers = c("Rlbp1", "Glul", "Sox9"),
  MGPC_progenitor = c("Ascl1", "Rgs13", "Dll1"),
  Neuronal_markers = c("Pou4f2", "Otx2", "Crx", "Vsx1")
)

all_genes <- unique(unlist(gene_dict))

# ----------------------------
# Load DE CSV
# ----------------------------
cat("Loading DE CSV:", de_csv, "\n")
de <- read.csv(de_csv, stringsAsFactors = FALSE)

# Check which genes exist in CSV
genes_found <- all_genes[all_genes %in% de$gene]
genes_missing <- setdiff(all_genes, genes_found)
if(length(genes_missing) > 0){
  warning("These genes were not found in DE CSV: ", paste(genes_missing, collapse = ", "))
}

# Subset DE to only genes found
heatmap_de <- de %>%
  filter(gene %in% genes_found) %>%
  select(gene, cluster, avg_log2FC) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)

# Convert to matrix
heatmap_matrix <- as.data.frame(heatmap_de)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix$gene <- NULL
heatmap_matrix <- as.matrix(heatmap_matrix)

# ----------------------------
# Row annotation (gene type)
# ----------------------------
row_annotation <- data.frame(
  GeneType = sapply(rownames(heatmap_matrix), function(g){
    if(g == gene_dict$Transgene) return("Transgene")
    if(g %in% gene_dict$MG_markers) return("MG_marker")
    if(g %in% gene_dict$MGPC_progenitor) return("MGPC_progenitor")
    if(g %in% gene_dict$Neuronal_markers) return("Neuronal_marker")
    return("Other")
  })
)
rownames(row_annotation) <- rownames(heatmap_matrix)

# ----------------------------
# Annotation colors
# ----------------------------
ann_colors <- list(
  GeneType = c(
    Transgene = "red",
    MG_marker = "blue",
    MGPC_progenitor = "purple",
    Neuronal_marker = "orange",
    Other = "grey"
  )
)

# ----------------------------
# Plot heatmap
# ----------------------------
cat("Plotting heatmap to:", output_png, "\n")
png(output_png, width = 1800, height = 1200, res = 150)
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  fontsize_col = 12,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  main = "Functional markers and Neurog2-9SA",
  border_color = "grey60"
)
dev.off()

cat("Heatmap saved successfully!\n")
cat("Check current working directory:", getwd(), "\n")

