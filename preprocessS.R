# -----------------------------
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(rtracklayer)
library(Biostrings)
library(Rsamtools)
library(Matrix)
set.seed(1234)

# -----------------------------
# Arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

mysample <- args[1]
myh5 <- paste0(mysample, "_filtered_feature_bc_matrix.h5")
mytsv <- paste0(mysample, "_atac_fragments.tsv.gz")

cat("Sample:", mysample, "\n")
cat("RNA+ATAC matrix:", myh5, "\n")
cat("ATAC fragments:", mytsv, "\n\n")

# -----------------------------
# Load custom genome annotation
# -----------------------------
gtf_file <- "neurog2.gtf"
annotation <- rtracklayer::import(gtf_file)
annotation <- annotation[annotation$type == "gene"]
seqlevelsStyle(annotation) <- "UCSC"

# -----------------------------
# Ensure FASTA is indexed
# -----------------------------
fasta_file <- "neurog2.fa"
if(!file.exists(paste0(fasta_file, ".fai"))) {
  Rsamtools::indexFa(fasta_file)
}

# -----------------------------
# Load RNA + ATAC counts
# -----------------------------
counts <- Read10X_h5(myh5)

# -----------------------------
# Create Seurat object with RNA
# -----------------------------
myObject <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  project = mysample
)

# -----------------------------
# Calculate mitochondrial percentage RIGHT AFTER creating object
# -----------------------------
rna_genes <- rownames(myObject)
mt_genes <- rna_genes[grepl("^mt-", rna_genes)]
cat("Mitochondrial genes found:", length(mt_genes), "\n")

if(length(mt_genes) > 0){
  myObject[["percent.mt"]] <- PercentageFeatureSet(myObject, features = mt_genes, assay = "RNA")
  cat("percent.mt calculated - Range:", range(myObject$percent.mt), "\n")
} else {
  myObject[["percent.mt"]] <- 0
}

# -----------------------------
# Create ATAC assay using custom annotation
# -----------------------------
myObject[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = mytsv,
  annotation = annotation
)

# -----------------------------
# Quality control
# -----------------------------
DefaultAssay(myObject) <- "ATAC"
myObject <- NucleosomeSignal(myObject)

# -----------------------------
# QC Plots
# -----------------------------
figure_name <- paste0(mysample, "_QC_vlnplot.png")
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

# -----------------------------
# Save preprocessed Seurat object
# -----------------------------
myRDS <- paste0(mysample, "_preprocessed.rds")
saveRDS(myObject, file = myRDS)
cat("Seurat object saved to", myRDS, "\n")
