library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(pheatmap)

set.seed(1234)
plan(sequential)
options(future.globals.maxSize = 80 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

input_rds <- paste0(mysample, "_annotated.rds")
myObject <- readRDS(input_rds)

cat("=== SUBSETTING MG AND MGPC CELLS ===\n")
subset_cells <- subset(myObject, idents = c("MG", "MGPC"))
cat("Number of cells after subset:", ncol(subset_cells), "\n")

# Save the subsetted Seurat object
subset_rds_file <- paste0(mysample, "_MG_MGPC_subset.rds")
saveRDS(subset_cells, file = subset_rds_file)
cat("Subsetted Seurat object saved to:", subset_rds_file, "\n")

# Set Sample as identity for TH1 vs TH2 comparison
Idents(subset_cells) <- "Sample"

# Prepare SCT if exists, else SCTransform
if ("SCT" %in% names(subset_cells@assays)) {
  DefaultAssay(subset_cells) <- "SCT"
  subset_cells <- PrepSCTFindMarkers(subset_cells)
} else {
  DefaultAssay(subset_cells) <- "RNA"
  subset_cells <- SCTransform(subset_cells, verbose = FALSE)
  subset_cells <- PrepSCTFindMarkers(subset_cells)
}

# Differential expression: TH1 vs TH2
cat("=== RUNNING DIFFERENTIAL EXPRESSION: TH1 vs TH2 ===\n")
dge_results <- FindMarkers(
  subset_cells,
  ident.1 = "TH1",
  ident.2 = "TH2",
  assay = "SCT",
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

output_csv <- paste0(mysample, "_TH1_vs_TH2_DGE.csv")
write.csv(dge_results, file = output_csv)
cat("DGE results saved to:", output_csv, "\n")

# Top 20 DE genes for heatmap
top_genes <- rownames(dge_results[order(dge_results$p_val_adj),])[1:20]

# Extract expression matrix for heatmap
expr_matrix <- as.matrix(GetAssayData(subset_cells, slot = "data")[top_genes, ])
scaled_matrix <- t(scale(t(expr_matrix)))

# Column annotation
annotation_col <- data.frame(
  Sample = subset_cells$Sample,
  row.names = colnames(subset_cells)
)

# Plot heatmap
heatmap_file <- paste0(mysample, "_TH1_vs_TH2_heatmap.png")
png(heatmap_file, width = 1500, height = 1200, res = 150)
pheatmap(
  scaled_matrix,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10,
  main = paste0(mysample, ": Top 20 DE Genes TH1 vs TH2")
)
dev.off()
cat("Heatmap saved to:", heatmap_file, "\n")

