#!/usr/bin/env Rscript

library(Signac)
library(Seurat)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)

set.seed(1234)

# -------------------------------
# Arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: script.R <sample_prefix> <marker_prefix>")
}
mysample <- args[1]
marker_prefix <- args[2]

# -------------------------------
# Load Seurat object
# -------------------------------
myRDS <- paste0(mysample, "_analysed.rds")
if(!file.exists(myRDS)) stop(paste("File not found:", myRDS))
myObject <- readRDS(myRDS)
DefaultAssay(myObject) <- "SCT"

# -------------------------------
# Load marker genes
# -------------------------------
file_path <- paste0(marker_prefix, ".txt")
if(!file.exists(file_path)) stop(paste("Marker file not found:", file_path))
marker_genes <- readLines(file_path)

# -------------------------------
# Fix duplicate factor levels
# -------------------------------
myObject$cell_type <- as.character(Idents(myObject))

# Make unique factor levels automatically
unique_levels <- make.unique(myObject$cell_type)
myObject$cell_type <- factor(myObject$cell_type, levels = unique(unique_levels))
Idents(myObject) <- myObject$cell_type

# -------------------------------
# Validate genes
# -------------------------------
valid_genes <- marker_genes[marker_genes %in% rownames(myObject)]
if(length(valid_genes) == 0) stop("No valid genes found in dataset.")

# -------------------------------
# FeaturePlots (individual PNGs)
# -------------------------------
for(gene in valid_genes){
  out_png <- paste0(mysample, "_", gene, "_FeaturePlot.png")
  png(filename = out_png, width = 1600, height = 1400, res = 200)
  print(
    FeaturePlot(myObject,
                features = gene,
                order = TRUE,
                reduction = "umap.rna",
                cols = c("lightgrey", "red"),
                pt.size = 0.1) +
      ggtitle(gene)
  )
  dev.off()
}

message("✅ All plots saved successfully in the current directory!")

