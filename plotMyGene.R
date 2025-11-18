#!/usr/bin/env Rscript

# ============================
# Full DE heatmap script (No Seurat required)
# ============================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

# ----------------------------
# User inputs
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript plot_DE_heatmap.R <DE_csv_file> <gene_of_interest>")
}

de_csv <- args[1]        # CSV file with DE results
gene_query <- args[2]    # e.g., "Neurog2_S9A"

output_png <- paste0(gene_query, "TopMarkers_heatmap", ".png")

# ----------------------------
# Load DE CSV
# ----------------------------
cat("Loading DE CSV:", de_csv, "\n")
de <- read.csv(de_csv, stringsAsFactors = FALSE)

# Fix gene names if they contain special characters (optional)
de$gene <- gsub("-", "_", de$gene)
gene_query <- gsub("-", "_", gene_query)

if(!gene_query %in% de$gene){
  stop(paste("Gene", gene_query, "not found in DE CSV"))
}

# ----------------------------
# Select genes for heatmap
# ----------------------------
# Include transgene
heatmap_genes <- gene_query

# ==== CHANGE #2: Select top N using smallest p_val_adj ====
top_markers <- de %>%
  filter(gene != gene_query) %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%      # ← sorted by significance
  slice_head(n = 10) %>%        # ← top 5 most significant
  pull(gene)

heatmap_genes <- unique(c(heatmap_genes, top_markers))

cat("Genes included in heatmap:\n")
print(heatmap_genes)

# ----------------------------
# Prepare matrix: genes x clusters
# ----------------------------

# ==== CHANGE #1: enforce cluster order ====
cluster_order <- c("MG", "MGPC", "BC", "AC", "Rod", "Cones")

heatmap_df <- de %>%
  filter(gene %in% heatmap_genes) %>%
  select(gene, cluster, avg_log2FC) %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%  # enforce order
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)

# Convert to matrix
heatmap_matrix <- as.data.frame(heatmap_df)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix$gene <- NULL
heatmap_matrix <- as.matrix(heatmap_matrix)

# ----------------------------
# Plot heatmap and save PNG
# ----------------------------
cat("Saving heatmap to:", output_png, "\n")
png(output_png, width = 1500, height = 1200, res = 150)
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,   # do not reorder columns—your custom order stays
  fontsize_row = 10,
  fontsize_col = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = paste("Top markers +", gene_query),
  display_numbers = FALSE,
  border_color = "grey60"
)
dev.off()

cat("Heatmap saved successfully!\n")
cat("Check current working directory:", getwd(), "\n")

