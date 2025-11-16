library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1] 
myRDS <- paste(mysample, "_annotated.rds", sep="")
myRDS

myObject <- readRDS(myRDS)


# Finder all markers
DefaultAssay(myObject) <- "ATAC"
myObject.atac.markers <- FindAllMarkers(myObject, assay = "ATAC", test.use = "wilcox", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1) #was 0.25

head(myObject.atac.markers)

myObject.atac.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

file_name <- paste(mysample, "DiffPeaks.csv", sep="")
write.csv(myObject.atac.markers, file=file_name)


Nearby_genes <- ClosestFeature(myObject, myObject.atac.markers$gene)
file_name <- paste(mysample, "annotatedDiffPeaks.csv", sep="")


write.csv(Nearby_genes, file=file_name)

myRDS <- paste(mysample, "_diffPeaks.rds", sep="")
saveRDS(myObject, file = myRDS)



