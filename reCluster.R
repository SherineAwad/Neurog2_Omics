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
Idents(myObject) <- "seurat_clusters"
clusters_to_remove <- c("11","19","26")
cells_remove <- WhichCells(myObject, idents = clusters_to_remove)
myObject <- subset(myObject, cells = setdiff(colnames(myObject), cells_remove))

# --------------------------------------------------------------
# RECLUSTER USING EXISTING SNN GRAPH
# --------------------------------------------------------------
myObject <- FindClusters(
  myObject,
  graph.name = "wsnn_harmony",
  resolution = 0.5,
  algorithm = 3,
  verbose = TRUE
)

# --------------------------------------------------------------
# UPDATE UMAP ON ORIGINAL HARMONY EMBEDDINGS
# --------------------------------------------------------------
myObject <- RunUMAP(
  myObject,
  reduction = "harmony.rna",
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
# 1. UMAP by cluster
figure_name <- paste0(mysample, "_Recluster_Clusters.png")
png(file=figure_name, width=1200, height=800)
print(
  DimPlot(myObject, reduction = "umap.wnn.harmony",
          group.by = "seurat_clusters", label=TRUE, repel=TRUE) +
    ggtitle("Reclustered (Removed 11,19, &26) - Updated Harmony UMAP")
)
dev.off()

# 2. UMAP by sample
figure_name <- paste0(mysample, "_Recluster_BySample.png")
png(file=figure_name, width=1200, height=800)
print(
  DimPlot(myObject, reduction = "umap.wnn.harmony",
          group.by = "orig.ident", repel=TRUE) +
    ggtitle("Reclustered by Sample - Updated Harmony UMAP")
)
dev.off()

# 3. Violin plot of UMI per cluster (detect poor-quality cells)
# Ensure cluster identity reflects Harmony WNN clusters
myObject$wnn_harmony_clusters <- Idents(myObject)  # store Harmony WNN clusters in a new column
Idents(myObject) <- "wnn_harmony_clusters"

# Violin plot of UMI per Harmony WNN cluster
DefaultAssay(myObject) <- "RNA"
figure_name <- paste0(mysample, "_UMI_Violin_HarmonyWNN.png")
png(file = figure_name, width = 1200, height = 800)
print(
  VlnPlot(myObject, features = "nCount_RNA", group.by = "wnn_harmony_clusters") +
    ggtitle("UMI Counts per Harmony WNN Cluster - Identify Poor-Quality Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

# --------------------------------------------------------------
# SAVE RESULT
# --------------------------------------------------------------
output_rds <- paste0(mysample, "_reClustered.rds")
saveRDS(myObject, file=output_rds)

