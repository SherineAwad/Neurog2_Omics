# Load required libraries
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)

# Set seed for reproducibility
set.seed(1234)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 9) {
  stop("Usage: Rscript filterS.R SAMPLE nATAC_max nRNA_max nATAC_min nRNA_min nFeature_RNA_min nucl_max percent_mt_max tss_min")
}

mysample <- args[1]
nATAC1 <- as.numeric(args[2])
nRNA1 <- as.numeric(args[3])
nATAC2 <- as.numeric(args[4])
nRNA2 <- as.numeric(args[5])
features <- as.numeric(args[6])
nucl <- as.numeric(args[7])
mt <- as.numeric(args[8])
tss <- as.numeric(args[9])  # ignored for now

cat("Filtering parameters:\n")
cat("nCount_ATAC <", nATAC1, "\n")
cat("nCount_RNA <", nRNA1, "\n")
cat("nCount_ATAC >", nATAC2, "\n")
cat("nCount_RNA >", nRNA2, "\n")
cat("nFeature_RNA >", features, "\n")
cat("nucleosome_signal <", nucl, "\n")
cat("percent.mt <", mt, "\n")
cat("TSS (ignored):", tss, "\n\n")

# Read Seurat object
myRDS <- paste0(mysample, "_preprocessed.rds")
if(!file.exists(myRDS)) stop("File not found: ", myRDS)
myObject <- readRDS(myRDS)
cat("Loaded Seurat object:", myRDS, "with", ncol(myObject), "cells\n\n")

# QC summary before filtering
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nucleosome_signal", "percent.mt")
print(summary(myObject@meta.data[, qc_metrics]))

# QC plot before filtering
png(paste0(mysample, "_QC_vlnplot_before_filtering.png"), width=2400, height=1200, res=150)
VlnPlot(myObject, features=qc_metrics, pt.size=0.1, ncol=5)
dev.off()
cat("QC plot before filtering saved.\n\n")

# Filter cells
myObject <- subset(
  x = myObject,
  subset = nCount_ATAC < nATAC1 &
           nCount_RNA < nRNA1 &
           nCount_ATAC > nATAC2 &
           nCount_RNA > nRNA2 &
           nFeature_RNA > features &
           nucleosome_signal < nucl &
           percent.mt < mt
)

if(ncol(myObject) == 0) {
  stop("Error: Filtering removed all cells. Relax thresholds.")
}
cat("Number of cells after filtering:", ncol(myObject), "\n\n")

# QC plot after filtering
png(paste0(mysample, "_QC_vlnplot_after_filtering.png"), width=2400, height=1200, res=150)
VlnPlot(myObject, features=qc_metrics, pt.size=0.1, ncol=5)
dev.off()
cat("QC plot after filtering saved.\n\n")

# Save filtered object
filtered_rds <- paste0(mysample, "_filtered.rds")
saveRDS(myObject, file=filtered_rds)
cat("Filtered Seurat object saved to:", filtered_rds, "\n")

