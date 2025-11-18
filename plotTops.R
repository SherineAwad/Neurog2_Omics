#!/usr/bin/env Rscript

# ============================
# Full DE Heatmap Script (Top N markers per cluster)
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

de_csv <- args[1]        # CSV file with DE results
top_n <- as.numeric(args[2])

if(is.na(top_n) || top_n <= 0){
  stop("Please provide a valid number for top_N_markers")
}

output_png <- paste0("Top", top_n, "_Markers_heatmap.png")

# ----------------------------
# Load DE CSV
# ----------------------------
cat("Loading DE CSV:", de_csv, "\n")
de <- read.csv(de_csv, stringsAsFactors = FALSE)

# Fix gene names if they contain special characters (optional)
de$gene <- gsub("-", "_", de$gene)

# ----------------------------
# Select top N markers per cluster
# ----------------------------
top_markers <- de %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

cat("Number of unique genes included in heatmap:", length(top_markers), "\n")
cat("Genes included in heatmap (first 20 shown):\n")
print(head(top_markers, 20))

# ----------------------------
# Fixed cluster order
# ----------------------------
cluster_order <- c('MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones')

# ----------------------------
# Prepare matrix: genes x clusters
# ----------------------------
heatmap_df <- de %>%
  filter(gene %in% top_markers) %>%
  select(gene, cluster, avg_log2FC) %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)

# Convert to matrix
heatmap_matrix <- as.data.frame(heatmap_df)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix$gene <- NULL

# Ensure exact column order and fill missing clusters with 0
for(cl in cluster_order){
  if(!cl %in% colnames(heatmap_matrix)){
    heatmap_matrix[[cl]] <- 0
  }
}
heatmap_matrix <- heatmap_matrix[, cluster_order, drop = FALSE]

# ----------------------------
# Plot heatmap and save PNG
# ----------------------------
cat("Saving heatmap to:", output_png, "\n")
png(output_png, width = 1500, height = 1200, res = 150)
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # keep fixed column order
  fontsize_row = 10,
  fontsize_col = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = paste("Top", top_n, "markers per cluster"),
  display_numbers = FALSE,
  border_color = "grey60"
)
dev.off()

cat("Heatmap saved successfully!\n")
cat("Check current working directory:", getwd(), "\n")

