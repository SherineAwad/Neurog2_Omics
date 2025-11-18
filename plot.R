#!/usr/bin/env Rscript

# ======================================
# Full functional heatmap with Neurog2-9SA and expanded gene lists
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
Transgene <- c("Neurog2-9SA")

MG_genes <- c(
  "Rlbp1", "Glul", "Sox9", "Aqp4", "Slc1a3", "Vim"
)

MGPC_genes <- c(
  "Ascl1", "Rgs13", "Dll1", "Nuak1", "Mapt", "Slc44a5"
)

Neuronal_genes <- c(
  "Pou4f2", "Otx2", "Crx", "Vsx1", "Tubb3", "NeuN"
)

Progenitor_genes <- c(
  "Hes1", "Hes5", "Pax6", "Sox2", "Nestin", "Ccnd1", "Mki67"
)

# Combine into dictionary
gene_dict <- list(
  Transgene = Transgene,
  MG_markers = MG_genes,
  MGPC_progenitor = MGPC_genes,
  Neuronal_markers = Neuronal_genes,
  Progenitor_markers = Progenitor_genes
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
    if(g %in% gene_dict$Transgene) return("Neurog2-S9A")
    if(g %in% gene_dict$MG_markers) return("MG")
    if(g %in% gene_dict$MGPC_progenitor) return("MGPC")
    if(g %in% gene_dict$Neuronal_markers) return("Neuronal")
    if(g %in% gene_dict$Progenitor_markers) return("Progenitor")
    return("Other")
  })
)
rownames(row_annotation) <- rownames(heatmap_matrix)

# ----------------------------
# Annotation colors
# ----------------------------
ann_colors <- list(
  GeneType = c(
    `Neurog2-S9A` = "red",
    MG = "blue",
    MGPC = "purple",
    Neuronal = "orange",
    Progenitor = "green", 
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

