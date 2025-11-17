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

myRDS <- paste(mysample, ".rds", sep="")
myObject <- readRDS(myRDS)

#### Analyze RNA part
DefaultAssay(myObject) <- "RNA"
myObject <- SCTransform(myObject, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'UMAP_')

#### Apply Harmony to RNA PCA
myObject <- RunHarmony(myObject, group.by.vars = "orig.ident", reduction = "pca", assay.use = "SCT", reduction.save = "harmony.rna")

#### Run UMAP on Harmony-corrected RNA
myObject <- RunUMAP(myObject, dims = 1:30, reduction = "harmony.rna", reduction.name = 'umap.rna.harmony', reduction.key = 'rnaHarmony_')

#### Analyze ATAC part
DefaultAssay(myObject) <- "ATAC"
myObject <- RunTFIDF(myObject)
myObject <- FindTopFeatures(myObject, min.cutoff = 20)
myObject <- RunSVD(myObject)
myObject <- RunUMAP(myObject, dims = 2:30, reduction = 'lsi', reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#### Apply Harmony to ATAC LSI - use project.dim = FALSE
myObject <- RunHarmony(myObject,
                       group.by.vars = "orig.ident",
                       reduction = "lsi",
                       assay.use = "ATAC",
                       project.dim = FALSE,
                       reduction.save = "harmony.atac")

#### Run UMAP on Harmony-corrected ATAC
myObject <- RunUMAP(myObject, dims = 2:30, reduction = "harmony.atac", reduction.name = "umap.atac.harmony", reduction.key = "atacHarmony_")

figure_name <- paste(mysample, "_depth_cor.png", sep="")
png(file=figure_name, width=1200, height=800)
DepthCor(myObject)
dev.off()

### Combine RNA and ATAC with WNN analysis using both original and Harmony-corrected reductions
# First do original WNN
myObject <- FindMultiModalNeighbors(
  myObject,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:20),
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn"
)

myObject <- RunUMAP(myObject, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_")
myObject <- FindClusters(myObject, graph.name = "wsnn", resolution = 1.2, algorithm = 3, verbose = FALSE)

# Now do Harmony-corrected WNN
myObject <- FindMultiModalNeighbors(
  myObject,
  reduction.list = list("harmony.rna", "harmony.atac"),
  dims.list = list(1:30, 2:20),
  knn.graph.name = "wknn_harmony",
  snn.graph.name = "wsnn_harmony",
  weighted.nn.name = "weighted.nn.harmony"
)

myObject <- RunUMAP(myObject, nn.name = "weighted.nn.harmony", reduction.name = "umap.wnn.harmony", reduction.key = "wnnHarmony_")
myObject <- FindClusters(myObject, graph.name = "wsnn_harmony", resolution = 1.2, algorithm = 3, verbose = FALSE)

myObject$Sample <- myObject$orig.ident

# -----------------------------
# Original WNN plots by cluster
figure_name <- paste(mysample, "_Clusters.png", sep="")
png(file = figure_name, width = 1200, height = 800)
p1 <- DimPlot(myObject, reduction = "umap.wnn", group.by = "orig.ident", repel = TRUE) + ggtitle("WNN")
p2 <- DimPlot(myObject, reduction = "umap.wnn", label = TRUE, repel = TRUE) + ggtitle("WNN")
print(p1 + p2)
dev.off()

# Add Harmony-corrected cluster plots
figure_name <- paste(mysample, "_Clusters_harmony.png", sep="")
png(file = figure_name, width = 1200, height = 800)
p1 <- DimPlot(myObject, reduction = "umap.wnn.harmony", group.by = "orig.ident", repel = TRUE) + ggtitle("WNN Harmony")
p2 <- DimPlot(myObject, reduction = "umap.wnn.harmony", label = TRUE, repel = TRUE) + ggtitle("WNN Harmony")
print(p1 + p2)
dev.off()

p1 <- DimPlot(myObject, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(myObject, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(myObject, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("WNN")

# Add Harmony-corrected versions
p1_harmony <- DimPlot(myObject, reduction = "umap.rna.harmony", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("RNA Harmony")
p2_harmony <- DimPlot(myObject, reduction = "umap.atac.harmony", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("ATAC Harmony")
p3_harmony <- DimPlot(myObject, reduction = "umap.wnn.harmony", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("WNN Harmony")

figure_name <- paste(mysample, "_Cluster_VlnPlot1.png", sep="")
png(file = figure_name, width = 1200, height = 800)
print(p1 + p3)
dev.off()

figure_name <- paste(mysample, "_Cluster_VlnPlot2.png", sep="")
png(file = figure_name, width = 1200, height = 800)
print(p2 + p3)
dev.off()

# Add Harmony comparison plots
figure_name <- paste(mysample, "_Cluster_VlnPlot1_harmony.png", sep="")
png(file = figure_name, width = 1200, height = 800)
print(p1_harmony + p3_harmony)
dev.off()

figure_name <- paste(mysample, "_Cluster_VlnPlot2_harmony.png", sep="")
png(file = figure_name, width = 1200, height = 800)
print(p2_harmony + p3_harmony)
dev.off()

p1 <- DimPlot(myObject, reduction = "umap.rna", group.by = "orig.ident", repel = TRUE) + ggtitle("RNA")
p3 <- DimPlot(myObject, reduction = "umap.wnn", group.by = "orig.ident", repel = TRUE) + ggtitle("WNN")

# Add Harmony versions
p1_harmony <- DimPlot(myObject, reduction = "umap.rna.harmony", group.by = "orig.ident", repel = TRUE) + ggtitle("RNA Harmony")
p3_harmony <- DimPlot(myObject, reduction = "umap.wnn.harmony", group.by = "orig.ident", repel = TRUE) + ggtitle("WNN Harmony")

figure_name <- paste(mysample, "_ident_VlnPlot.png", sep="")
png(file = figure_name, width = 1200, height = 800)
print(p1 + p3)
dev.off()

# Add Harmony ident plot
figure_name <- paste(mysample, "_ident_VlnPlot_harmony.png", sep="")
png(file = figure_name, width = 1200, height = 800)
print(p1_harmony + p3_harmony)
dev.off()

# -----------------------------
# WNN plot by sample for batch effect (before and after Harmony)
figure_name <- paste(mysample, "_WNN_by_sample.png", sep="")
png(file = figure_name, width = 1200, height = 800)
DimPlot(myObject, reduction = "umap.wnn", group.by = "Sample", repel = TRUE) + ggtitle("WNN by Sample")
dev.off()

figure_name <- paste(mysample, "_WNN_by_sample_harmony.png", sep="")
png(file = figure_name, width = 1200, height = 800)
DimPlot(myObject, reduction = "umap.wnn.harmony", group.by = "Sample", repel = TRUE) + ggtitle("WNN Harmony by Sample")
dev.off()

# -----------------------------
DefaultAssay(myObject) <- "ATAC"
myRDS <- paste(mysample, "_analysed.rds", sep="")
saveRDS(myObject, file = myRDS)
