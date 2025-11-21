library(Seurat)
library(pheatmap)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

# -------------------- PARAMETERS --------------------
top_n <- 30            # Number of genes to show in heatmap
logfc_threshold <- 0.5 # Minimum absolute log fold-change
zscore_clip <- 2.5     # Max/min z-score for visualization
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

# Filter DE genes based on logFC threshold
filtered_genes <- dge_results %>%
  filter(abs(avg_log2FC) >= logfc_threshold)

if (nrow(filtered_genes) == 0) stop("No genes pass the thresholds!")

# Rank genes by descending absolute logFC
ranked_genes <- filtered_genes %>%
  arrange(desc(abs(avg_log2FC)))

# Take top N genes
top_genes <- head(rownames(ranked_genes), top_n)

# Ensure genes exist in the Seurat object
assay_to_use <- if ("SCT" %in% names(subset_cells@assays)) "SCT" else "RNA"
existing_genes <- rownames(subset_cells[[assay_to_use]])
top_genes <- intersect(top_genes, existing_genes)
cat("Number of genes actually plotted:", length(top_genes), "\n")
if (length(top_genes) == 0) stop("None of the top genes exist in assay!")

# Extract expression matrix (SCT or RNA)
expr_matrix <- as.matrix(GetAssayData(subset_cells, assay = assay_to_use, slot = "data")[top_genes, ])

# Row-wise z-score
scaled_matrix <- t(apply(expr_matrix, 1, function(x) {
  z <- (x - mean(x)) / sd(x)
  # Clip extreme values for better visualization
  pmin(pmax(z, -zscore_clip), zscore_clip)
}))

# Column annotation for TH2/TH1
annotation_col <- data.frame(
  Sample = subset_cells$Sample,
  row.names = colnames(subset_cells)
)

# ----------------------------------------------------
# CUSTOM COLOR PALETTE (ONLY CHANGE MADE)
# ----------------------------------------------------
custom_palette <- colorRampPalette(c(
  "#B5D1E1",  # light blue (low)
  "#C0DAEA",  # very light blue
  "#FFFFFF",  # white (zero)
  "#FDFEFE",  # near-white
  "#E5A07E",  # light salmon
  "#C94832",  # medium red
  "#B5332A"   # deep red (high)
))(100)

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
  color = custom_palette
)
dev.off()
cat("Heatmap saved to:", heatmap_file, "\n")

