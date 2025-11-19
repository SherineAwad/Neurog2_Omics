#!/usr/bin/env Rscript

# ============================
# DE Heatmap with Diagonal Gene Sorting
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
  stop("Usage: Rscript plot_DE_heatmap.R <DE_csv_file> <top_N_markers>")
}

de_csv <- args[1]        
top_n <- as.numeric(args[2])

if(is.na(top_n) || top_n <= 0){
  stop("Please provide a valid positive number for top_N_markers")
}

output_png <- paste0("Top",top_n, "_Markers_heatmap.png")

# ----------------------------
# Load DE CSV
# ----------------------------
cat("Loading DE CSV:", de_csv, "\n")
de <- read.csv(de_csv, stringsAsFactors = FALSE)

# Fix gene names (optional)
de$gene <- gsub("-", "_", de$gene)

# ----------------------------
# Fixed cluster order
# ----------------------------
cluster_order <- c('MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones')

# ----------------------------
# Select top N markers per cluster
# ----------------------------
top_genes <- de %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = top_n) %>%
  pull(gene) %>%
  unique()

cat("Number of genes selected:", length(top_genes), "\n")

# ----------------------------
# Prepare gene × cluster matrix
# ----------------------------
heatmap_df <- de %>%
  filter(gene %in% top_genes) %>%
  select(gene, cluster, avg_log2FC) %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)

# Convert to matrix
heatmap_matrix <- as.data.frame(heatmap_df)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix$gene <- NULL

# Add missing clusters (if any)
for(cl in cluster_order){
  if(!cl %in% colnames(heatmap_matrix)){
    heatmap_matrix[[cl]] <- 0
  }
}
heatmap_matrix <- heatmap_matrix[, cluster_order, drop = FALSE]

# ----------------------------
# DIAGONAL SORTING OF GENES
# ----------------------------

# Find cluster index (position of max logFC)
cluster_index_map <- setNames(seq_along(cluster_order), cluster_order)

gene_order_df <- heatmap_matrix %>%
  mutate(
    gene = rownames(.),
    best_cluster = apply(., 1, function(x) cluster_order[which.max(x)]),
    best_cluster_index = cluster_index_map[best_cluster],
    max_lfc = apply(., 1, max)
  ) %>%
  arrange(best_cluster_index, desc(max_lfc))

sorted_genes <- gene_order_df$gene

# Reorder matrix based on diagonal logic
heatmap_matrix <- heatmap_matrix[sorted_genes, , drop = FALSE]

cat("Diagonal sorting complete.\n")

# ----------------------------
# Plot heatmap
# ----------------------------
cat("Saving heatmap to:", output_png, "\n")

png(output_png, width = 1500, height = 1200, res = 150)
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = paste("Top", top_n, "Markers per Cluster (Diagonal Sorted)"),
  border_color = NA
)
dev.off()

cat("Heatmap saved successfully!\n")
cat("Working directory:", getwd(), "\n")

