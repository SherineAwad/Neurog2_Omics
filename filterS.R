# Load required libraries
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)

# Set a random seed for reproducibility
set.seed(1234)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

mysample <- args[1]
nATAC1 <- as.double(args[2])
nRNA1 <- as.double(args[3])
nATAC2 <- as.double(args[4])
nRNA2 <- as.double(args[5])
features <- as.double(args[6])
nucl <- as.double(args[7])
tss <- as.double(args[8])
mt <- as.double(args[9])

# Print values (for debugging)
nATAC1; nRNA1; nATAC2; nRNA2; features; nucl; tss; mt

# Build file path to preprocessed Seurat object
myRDS <- paste0(mysample, "_preprocessed.rds")
myRDS

# Read Seurat object
myObject <- readRDS(myRDS)

# Example command:
# Rscript filter.R sample1 100000 30000 0 0 200 0.5 8 1

# Filter cells
# Original filtering line using TSS:
# myObject <- subset(x = myObject,
#                    subset = nCount_ATAC < nATAC1 &
#                             nCount_RNA < nRNA1 &
#                             nCount_ATAC > nATAC2 &
#                             nCount_RNA > nRNA2 &
#                             nFeature_RNA > features &
#                             nucleosome_signal < nucl &
#                             TSS.enrichment > tss &      # <- TSS-related line
#                             percent.mt < mt)

# We ignore TSS:
myObject <- subset(x = myObject,
                   subset = nCount_ATAC < nATAC1 &
                            nCount_RNA < nRNA1 &
                            nCount_ATAC > nATAC2 &
                            nCount_RNA > nRNA2 &
                            nFeature_RNA > features &
                            nucleosome_signal < nucl &
                            percent.mt < mt)

# -----------------------------
# QC Plots - PNG instead of PDF
# -----------------------------
figure_name <- paste0(mysample, "_QC_vlnplot_after_filtering.png")
png(file = figure_name, width = 2400, height = 1200, res = 150)
VlnPlot(
  object = myObject,
  features = c(
    "nCount_RNA", "nFeature_RNA",
    "nCount_ATAC", "nucleosome_signal", "percent.mt"
  ),
  pt.size = 0.1,
  ncol = 5
)
dev.off()
cat("QC plot saved to", figure_name, "\n")

# Save filtered object
myRDS <- paste0(mysample, "_filtered.rds")
saveRDS(myObject, file = myRDS)
cat("Filtered Seurat object saved to", myRDS, "\n")
