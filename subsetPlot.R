library(Seurat)
library(pheatmap)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

# -------------------- PARAMETERS --------------------
top_n <- 30            # Number of genes to show in heatmap
logfc_threshold <- 0.5 # Minimum absolute log fold-change
# ----------------------------------------------------

# Load previously subsetted MG/MGPC Seurat object
subset_rds_file <- paste0(mysample, "_MG_MGPC_subset.rds")
subset_cells <- readRDS(subset_rds_file)
cat("Loaded subsetted Seurat object:", subset_rds_file, "\n")
cat("Number of cells:", ncol(subset_cells), "\n")

# Load DE results CSV
dge_file <- paste0(mysample, "_MG_MGPC_TH2_vs_TH1_DGE.csv")
dge_results <- read.csv(dge_file, row.names = 1)
cat("Loaded DE results:", dge_file, "\n")

# Filter DE genes based on thresholds
filtered_genes <- dge_results %>%
  filter(abs(avg_log2FC) >= logfc_threshold)

if (nrow(filtered_genes) == 0) {
  stop("No genes pass the thresholds!")
}

# Rank genes: first by adjusted p-value, then by abs(logFC)
ranked_genes <- filtered_genes %>%
  arrange(desc(abs(avg_log2FC)))

# Take top N genes
top_genes <- head(rownames(ranked_genes), top_n)

# Ensure only genes that exist in the Seurat object are used
assay_to_use <- if ("SCT" %in% names(subset_cells@assays)) "SCT" else "RNA"
existing_genes <- rownames(subset_cells[[assay_to_use]])
top_genes <- intersect(top_genes, existing_genes)
cat("Number of genes actually plotted (exist in assay):", length(top_genes), "\n")

if (length(top_genes) == 0) {
  stop("None of the top genes exist in the Seurat object assay!")
}

# Extract expression matrix
expr_matrix <- as.matrix(GetAssayData(subset_cells, assay = assay_to_use, slot = "data")[top_genes, ])
scaled_matrix <- t(scale(t(expr_matrix)))

# Column annotation for TH2/TH1
annotation_col <- data.frame(
  Sample = subset_cells$Sample,
  row.names = colnames(subset_cells)
)

# Plot heatmap
heatmap_file <- paste0(mysample, "_TH2_vs_TH1_heatmap_top", length(top_genes), ".png")
png(heatmap_file, width = 1200, height = 1000, res = 150)
pheatmap(
  scaled_matrix,
  annotation_col = annotation_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 10,
  main = paste0(mysample, ": Top ", length(top_genes), " DE Genes TH2 vs TH1"),
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
dev.off()
cat("Heatmap saved to:", heatmap_file, "\n")

