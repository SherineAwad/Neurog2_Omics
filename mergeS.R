library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)
args <- commandArgs(trailingOnly = TRUE)

myRDS <- args[1]
mysample1 <- args[2]
mysample2 <- args[3]


myObject1 <- paste(mysample1,"_filtered.rds", sep="")
myObject2 <- paste(mysample2,"_filtered.rds", sep="") 

myObject1
myObject2
myRDS

myObject1 <- readRDS(myObject1)
myObject2 <- readRDS(myObject2)

myObject <- merge(myObject1, y = c(myObject2), add.cell.ids = c("Control", "OE"), project = myRDS)
myRDS <- paste(myRDS, ".rds", sep="")
saveRDS(myObject, file = myRDS)

