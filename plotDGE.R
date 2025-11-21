library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) stop("Usage: Rscript plot_DE_heatmap.R <DE_csv_file> <top_N_markers>")

de_csv <- args[1]
top_n <- as.numeric(args[2])
if(is.na(top_n) || top_n <= 0) stop("Please provide a valid positive number for top_N_markers")

output_png <- paste0("Top",top_n, "_Markers_heatmap.png")

de <- read.csv(de_csv, stringsAsFactors = FALSE)
de$gene <- gsub("-", "_", de$gene)

cluster_order <- c('MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones')

top_genes <- de %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = top_n) %>%
  pull(gene) %>%
  unique()

heatmap_df <- de %>%
  filter(gene %in% top_genes) %>%
  select(gene, cluster, avg_log2FC) %>%
  mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)

heatmap_matrix <- as.data.frame(heatmap_df)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix$gene <- NULL

for(cl in cluster_order){
  if(!cl %in% colnames(heatmap_matrix)) heatmap_matrix[[cl]] <- 0
}
heatmap_matrix <- heatmap_matrix[, cluster_order, drop = FALSE]

heatmap_matrix_z <- t(apply(heatmap_matrix, 1, function(x) {
  if(sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)
}))

cluster_index_map <- setNames(seq_along(cluster_order), cluster_order)
gene_order_df <- heatmap_matrix_z %>%
  as.data.frame() %>%
  mutate(
    gene = rownames(.),
    best_cluster = apply(., 1, function(x) cluster_order[which.max(x)]),
    best_cluster_index = cluster_index_map[best_cluster],
    max_lfc = apply(., 1, max)
  ) %>%
  arrange(best_cluster_index, desc(max_lfc))

sorted_genes <- gene_order_df$gene
heatmap_matrix_z <- heatmap_matrix_z[sorted_genes, , drop = FALSE]

# --- NEW COLOR PALETTE (only change) ---
custom_colors <- colorRampPalette(c(
  "#B5D1E1",  # light blue
  "#C0DAEA",  # very light blue
  "#FFFFFF",  # white
  "#FDFEFE",  # near-white
  "#E5A07E",  # light salmon
  "#C94832",  # medium red
  "#B5332A"   # deep red
))(200)

png(output_png, width = 1500, height = 1200, res = 150)
pheatmap(
  heatmap_matrix_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  color = custom_colors,
  main = paste("Top", top_n, "Markers per Cluster (Diagonal Sorted, Z-score)"),
  border_color = NA
)
dev.off()

