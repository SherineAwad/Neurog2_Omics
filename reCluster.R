#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

# --------------------------------------------------------------
# LOAD ORIGINAL OBJECT
# --------------------------------------------------------------
input_rds <- paste0(mysample, "_analysed.rds")
myObject <- readRDS(input_rds)

# --------------------------------------------------------------
# CHECK REQUIRED REDUCTIONS AND GRAPHS
# --------------------------------------------------------------
if (!"umap.wnn.harmony" %in% names(myObject@reductions)) stop("umap.wnn.harmony missing!")
if (!"wsnn_harmony" %in% names(myObject@graphs)) stop("wsnn_harmony missing!")

# --------------------------------------------------------------
# REMOVE CLUSTERS 11 AND 19
# --------------------------------------------------------------
# Make sure Idents are set to seurat_clusters
Idents(myObject) <- "seurat_clusters"

clusters_to_remove <- c("11","19")
cells_remove <- WhichCells(myObject, idents = clusters_to_remove)
myObject <- subset(myObject, cells = setdiff(colnames(myObject), cells_remove))

# --------------------------------------------------------------
# RECLUSTER USING EXISTING SNN GRAPH
# --------------------------------------------------------------
myObject <- FindClusters(
  myObject,
  graph.name = "wsnn_harmony",
  resolution = 0.5,  #changed from 0.1 to 0.5 
  algorithm = 3,
  verbose = TRUE
)

# --------------------------------------------------------------
# UPDATE UMAP ON ORIGINAL HARMONY EMBEDDINGS
# --------------------------------------------------------------
# Keep the original Harmony embeddings: "harmony.rna" + "harmony.atac"
myObject <- RunUMAP(
  myObject,
  reduction = "harmony.rna",   # you can choose either rna, atac, or concatenate
  dims = 1:50,
  reduction.name = "umap.wnn.harmony",
  reduction.key = "wnnHarmony_",
  spread = 1,
  min.dist = 0.5
)

# --------------------------------------------------------------
# OPTIONAL: compress UMAP coordinates for closer clusters
# --------------------------------------------------------------
coords <- Embeddings(myObject, "umap.wnn.harmony")
coords <- coords * 0.6
myObject[["umap.wnn.harmony"]]@cell.embeddings <- coords

# --------------------------------------------------------------
# PLOTS
# --------------------------------------------------------------
figure_name <- paste0(mysample, "_Recluster_Clusters.png")
png(file=figure_name, width=1200, height=800)
print(
  DimPlot(myObject, reduction = "umap.wnn.harmony",
          group.by = "seurat_clusters", label=TRUE, repel=TRUE) +
    ggtitle("Reclustered (Removed 11 & 19) - Updated Harmony UMAP")
)
dev.off()

figure_name <- paste0(mysample, "_Recluster_BySample.png")
png(file=figure_name, width=1200, height=800)
print(
  DimPlot(myObject, reduction = "umap.wnn.harmony",
          group.by = "orig.ident", repel=TRUE) +
    ggtitle("Reclustered by Sample - Updated Harmony UMAP")
)
dev.off()

# --------------------------------------------------------------
# SAVE RESULT
# --------------------------------------------------------------
output_rds <- paste0(mysample, "_reClustered.rds")
saveRDS(myObject, file=output_rds)

