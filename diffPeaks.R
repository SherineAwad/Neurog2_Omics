library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)

plan(sequential)   # stays the same (or 'multicore' if you want parallel)
options(future.globals.maxSize = 80 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]
myRDS <- paste0(mysample, "_with_DGE.rds")
myObject <- readRDS(myRDS)

# Set ATAC as default assay
DefaultAssay(myObject) <- "ATAC"

# Find differential peaks
myObject.atac.markers <- FindAllMarkers(
  myObject, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Save diff peaks table
file_name <- paste0(mysample, "_DiffPeaks.csv")
write.csv(myObject.atac.markers, file=file_name)

# Annotate peaks with nearest genes
Nearby_genes <- ClosestFeature(myObject, myObject.atac.markers$gene)
file_name <- paste0(mysample, "_annotatedDiffPeaks.csv")
write.csv(Nearby_genes, file=file_name)

# Store diff_peaks in the Seurat object for downstream scripts
myObject@misc$diff_peaks <- myObject.atac.markers

# Also store the annotations separately if needed
myObject@misc$peak_annotations <- Nearby_genes

# Save the updated Seurat object
myRDS_out <- paste0(mysample, "_DiffPeaks.rds")
saveRDS(myObject, file = myRDS_out)
message("Saved Seurat object with diff_peaks: ", myRDS_out)

