#!/usr/bin/env Rscript

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(harmony)
set.seed(1234)

plan("multicore", workers = 8)
options(future.globals.maxSize = 20 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

# --------------------------------------------------------------
# LOAD PREVIOUS PROCESSED OBJECT
# --------------------------------------------------------------
input_rds <- paste(mysample, "_analysed.rds", sep="")
myObject <- readRDS(input_rds)

if (!"seurat_clusters" %in% colnames(myObject@meta.data)) {
  stop("seurat_clusters not found in metadata.")
}

# --------------------------------------------------------------
# REMOVE CLUSTERS 11 AND 19
# --------------------------------------------------------------
clusters_to_remove <- c("11", "19")
cells_remove <- WhichCells(myObject, idents = clusters_to_remove)
myObject <- subset(myObject, cells = setdiff(colnames(myObject), cells_remove))
Idents(myObject) <- "orig.ident"

# --------------------------------------------------------------
# REDUCED DIMENSIONS TO MAKE CLUSTERS CLOSER
# --------------------------------------------------------------

### RNA
DefaultAssay(myObject) <- "SCT"
myObject <- RunPCA(myObject, npcs = 10)

myObject <- RunHarmony(
  myObject,
  group.by.vars = "orig.ident",
  reduction = "pca",
  assay.use = "SCT",
  reduction.save = "harmony.rna.recluster",
  theta = 2,
  lambda = 2,
  sigma = 0.8,
  nclust = 50
)

### ATAC
DefaultAssay(myObject) <- "ATAC"
myObject <- RunSVD(myObject, n = 10)

myObject <- RunHarmony(
  myObject,
  group.by.vars = "orig.ident",
  reduction = "lsi",
  assay.use = "ATAC",
  project.dim = FALSE,
  reduction.save = "harmony.atac.recluster",
  theta = 2,
  lambda = 2,
  sigma = 0.8,
  nclust = 50
)

# --------------------------------------------------------------
# OPTIONAL: compress embeddings to bring clusters closer
# --------------------------------------------------------------
myObject[["harmony.rna.recluster"]]@cell.embeddings <- 
  myObject[["harmony.rna.recluster"]]@cell.embeddings * 0.5
myObject[["harmony.atac.recluster"]]@cell.embeddings <- 
  myObject[["harmony.atac.recluster"]]@cell.embeddings * 0.5

# --------------------------------------------------------------
# WNN
# --------------------------------------------------------------
myObject <- FindMultiModalNeighbors(
  myObject,
  reduction.list = list("harmony.rna.recluster", "harmony.atac.recluster"),
  dims.list = list(1:10, 2:10),
  knn.graph.name = "wknn_re",
  snn.graph.name = "wsnn_re",
  weighted.nn.name = "weighted.nn.re"
)

# --------------------------------------------------------------
# CLUSTER — LOWER RESOLUTION
# --------------------------------------------------------------
myObject <- FindClusters(
  myObject,
  graph.name = "wsnn_re",
  resolution = 0.3,  # merge small clusters → closer clusters
  algorithm = 3,
  verbose = FALSE
)

myObject$Sample <- myObject$orig.ident

# --------------------------------------------------------------
# RECOMPUTE UMAP — TIGHTER EMBEDDING
# --------------------------------------------------------------
myObject <- RunUMAP(
  myObject,
  nn.name = "weighted.nn.re",
  reduction.name = "umap.wnn.re",
  reduction.key = "wnnRe_",
  spread = 0.5,
  min.dist = 0.2
)

# --------------------------------------------------------------
# PLOTS
# --------------------------------------------------------------
figure_name <- paste(mysample, "_Recluster_Clusters.png", sep="")
png(file = figure_name, width = 1200, height = 800)
print(
  DimPlot(myObject, reduction = "umap.wnn.re",
          group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("Reclustered with Compressed & Tighter UMAP")
)
dev.off()

figure_name <- paste(mysample, "_Recluster_BySample.png", sep="")
png(file = figure_name, width = 1200, height = 800)
print(
  DimPlot(myObject, reduction = "umap.wnn.re",
          group.by = "Sample", repel = TRUE) +
  ggtitle("Reclustered Samples (Compressed & Tighter UMAP)")
)
dev.off()

# --------------------------------------------------------------
# SAVE RESULT
# --------------------------------------------------------------
output_rds <- paste(mysample, "_reclustered_tight.rds", sep="")
saveRDS(myObject, file = output_rds)

