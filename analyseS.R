library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)

plan("multicore", workers = 8)
options(future.globals.maxSize = 20 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

myRDS <- paste(mysample, ".rds", sep="")
myObject <- readRDS(myRDS)

#### Analyze RNA part
DefaultAssay(myObject) <- "RNA"
myObject <- SCTransform(myObject, verbose = FALSE) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'UMAP_')

#### Analyze ATAC part
DefaultAssay(myObject) <- "ATAC"
myObject <- RunTFIDF(myObject)
myObject <- FindTopFeatures(myObject, min.cutoff = 20)
myObject <- RunSVD(myObject)
myObject <- RunUMAP(myObject, dims = 2:30, reduction = 'lsi', reduction.name = "umap.atac", reduction.key = "atacUMAP_")

figure_name <- paste(mysample, "_depth_cor.png", sep="")
png(file=figure_name, width=1200, height=800)
DepthCor(myObject)
dev.off()

### Combine RNA and ATAC with WNN analysis
myObject <- FindMultiModalNeighbors(myObject, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
myObject <- RunUMAP(myObject, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key ="wnnUMAP_")
myObject <- FindClusters(myObject, graph.name = "wsnn", resolution = 1.2, algorithm = 3, verbose = FALSE)

myObject$Sample <- myObject$orig.ident

# -----------------------------
# Original WNN plots by cluster
figure_name <- paste(mysample, "_Clusters.png", sep="")
png(file =figure_name, width =1200, height=800)
DimPlot(myObject, reduction = "umap.wnn", group.by = "orig.ident",  repel = TRUE) + ggtitle("WNN")
DimPlot(myObject, reduction = "umap.wnn", label=TRUE, repel = TRUE) + ggtitle("WNN")
dev.off()

p1 <- DimPlot(myObject, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(myObject, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE,  repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(myObject, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("WNN")

figure_name <- paste(mysample, "_Cluster_VlnPlot1.png", sep="")
png(file=figure_name, width=1200, height=800)
p1 + p3 & NoAxes()
dev.off()

figure_name <- paste(mysample, "_Cluster_VlnPlot2.png", sep="")
png(file=figure_name, width=1200, height=800)
p2 + p3 & NoAxes()
dev.off()

# -----------------------------
# NEW: WNN plot by sample for batch effect
figure_name <- paste(mysample, "_Sample_VlnPlot.png", sep="")
png(file=figure_name, width=1200, height=800)
DimPlot(myObject, reduction = "umap.wnn", group.by = "Sample", repel = TRUE) + ggtitle("WNN by Sample")
dev.off()

# -----------------------------
DefaultAssay(myObject) <- "ATAC"
myRDS <- paste(mysample, "_analysed.rds", sep="")
saveRDS(myObject, file = myRDS)

