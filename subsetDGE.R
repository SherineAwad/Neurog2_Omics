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

# Set Sample as identity for TH2 vs TH1 comparison
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

cat("=== RUNNING DIFFERENTIAL EXPRESSION: TH2 vs TH1 ===\n")
dge_results <- FindMarkers(
  subset_cells,
  ident.1 = "TH2",
  ident.2 = "TH1",
  assay = "SCT",
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  )

output_csv <- paste0(mysample, "_MG_MGPC_TH2_vs_TH1_DGE.csv")
write.csv(dge_results, file = output_csv)
cat("DGE results saved to:", output_csv, "\n")

