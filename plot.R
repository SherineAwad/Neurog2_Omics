#!/usr/bin/env Rscript

# ======================================
# Full functional heatmap with forced Neurog2-S9A
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
Transgene <- "Neurog2-S9A"   # force inclusion

MG_genes <- c("Rlbp1", "Glul", "Sox9", "Aqp4", "Slc1a3", "Vim")
MGPC_genes <- c("Ascl1", "Rgs13", "Dll1", "Nuak1", "Mapt", "Slc44a5")
Progenitor_genes <- c("Hes1", "Hes5", "Pax6", "Sox2", "Nestin", "Ccnd1", "Mki67")
Neuronal_genes <- c("Pou4f2", "Otx2", "Crx", "Vsx1", "Tubb3", "NeuN")

gene_dict <- list(
  Transgene = Transgene,
  MG_markers = MG_genes,
  MGPC_markers = MGPC_genes,
  Progenitor_markers = Progenitor_genes,
  Neuronal_markers = Neuronal_genes
)

all_genes <- unique(c(Transgene, MG_genes, MGPC_genes, Progenitor_genes, Neuronal_genes))

# ----------------------------
# Load DE CSV
# ----------------------------
cat("Loading DE CSV:", de_csv, "\n")
de <- read.csv(de_csv, stringsAsFactors = FALSE)

# ----------------------------
# Force Neurog2-S9A presence
# ----------------------------
# If Neurog2-S9A is missing, create a fake row with zeros for plotting
if(!Transgene %in% de$gene){
  warning("Neurog2-S9A not found in CSV. Adding it with zero expression.")
  clusters <- unique(de$cluster)
  new_row <- data.frame(
    gene = Transgene,
    cluster = rep(clusters, each = 1),
    avg_log2FC = 0
  )
  de <- bind_rows(de, new_row)
}

# ----------------------------
# Subset DE to only genes in our list
# ----------------------------
heatmap_de <- de %>%
  filter(gene %in% all_genes) %>%
  select(gene, cluster, avg_log2FC) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)

# Convert to matrix
heatmap_matrix <- as.data.frame(heatmap_de)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix$gene <- NULL
heatmap_matrix <- as.matrix(heatmap_matrix)

# ----------------------------
# Row annotation
# ----------------------------
row_annotation <- data.frame(
  GeneType = sapply(rownames(heatmap_matrix), function(g){
    if(g == Transgene) return("Neurog2-S9A")
    if(g %in% MG_genes) return("MG")
    if(g %in% MGPC_genes) return("MGPC")
    if(g %in% Progenitor_genes) return("Progenitor")
    if(g %in% Neuronal_genes) return("Neuronal")
    return("Other")
  })
)
rownames(row_annotation) <- rownames(heatmap_matrix)

# ----------------------------
# Force grouping by gene type
# ----------------------------
gene_type_levels <- c("Neurog2-S9A", "MG", "MGPC", "Progenitor", "Neuronal", "Other")
row_annotation$GeneType <- factor(row_annotation$GeneType, levels = gene_type_levels)
heatmap_matrix <- heatmap_matrix[order(row_annotation$GeneType), ]
row_annotation <- row_annotation[order(row_annotation$GeneType), , drop = FALSE]

# ----------------------------
# Annotation colors
# ----------------------------
ann_colors <- list(
  GeneType = setNames(
    c("red", "blue", "purple", "green", "orange", "grey"),
    gene_type_levels
  )
)

# ----------------------------
# Plot heatmap
# ----------------------------
cat("Plotting heatmap to:", output_png, "\n")
png(output_png, width = 1800, height = 1200, res = 150)
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,   # keep our ordered groups
  cluster_cols = TRUE,
  fontsize_row = 10,
  fontsize_col = 12,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  main = "Functional markers and Neurog2-S9A",
  border_color = "grey60"
)
dev.off()

cat("Heatmap saved successfully!\n")
cat("Check current working directory:", getwd(), "\n")

