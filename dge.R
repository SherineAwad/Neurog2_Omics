#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(UpSetR)  # for UpSet plot
set.seed(1234)


library(future)

plan(sequential)   # stays the same (or 'multicore' if you want parallel)

# 80 GB limit
options(future.globals.maxSize = 80 * 1024^3)


args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

# --------------------------------------------------------------
# LOAD ANNOTATED OBJECT
# --------------------------------------------------------------
input_rds <- paste0(mysample, "_annotated.rds")
myObject <- readRDS(input_rds)

# --------------------------------------------------------------
# SETUP AND QUALITY CONTROL
# --------------------------------------------------------------
cat("=== DIFFERENTIAL GENE EXPRESSION ANALYSIS ===\n")
cat("Sample:", mysample, "\n")
cat("Number of cells:", ncol(myObject), "\n")
cat("Number of genes:", nrow(myObject), "\n")

celltypes <- levels(Idents(myObject))
cat("Available cell types:", paste(celltypes, collapse = ", "), "\n")

DefaultAssay(myObject) <- "RNA"

if ("SCT" %in% names(myObject@assays)) {
  cat("Using existing SCT assay...\n")
  DefaultAssay(myObject) <- "SCT"
  myObject <- PrepSCTFindMarkers(myObject)
} else {
  cat("Performing SCTransform normalization...\n")
  myObject <- SCTransform(myObject, verbose = FALSE)
  myObject <- PrepSCTFindMarkers(myObject)
}

# --------------------------------------------------------------
# DIFFERENTIAL EXPRESSION BETWEEN CELL TYPES
# --------------------------------------------------------------
cat("\n=== FINDING MARKERS BETWEEN CELL TYPES ===\n")

celltype_markers <- FindAllMarkers(
  myObject,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.1,      
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

write.csv(celltype_markers, paste0(mysample, "_all_celltype_markers.csv"), row.names = FALSE)

# --------------------------------------------------------------
# SAVE FINAL OBJECT
# --------------------------------------------------------------
cat("\n=== SAVING RESULTS ===\n")
myObject@misc$celltype_markers <- celltype_markers
saveRDS(myObject, file = paste0(mysample, "_with_DGE.rds"))

# --------------------------------------------------------------
# FINAL OUTPUT
# --------------------------------------------------------------
cat("\n=== ANALYSIS COMPLETE ===\n")

