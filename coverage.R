#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
})

parser <- ArgumentParser()
parser$add_argument("-r", "--rds", required=TRUE)
parser$add_argument("-g", "--gene", required=TRUE)
parser$add_argument("-o", "--out", default=NULL) 
args <- parser$parse_args()

cat("Loading object...\n")
obj <- readRDS(args$rds)

# ---------------------------
# Extract gene coordinates
# ---------------------------
anno <- Annotation(obj)
gene <- args$gene

cat("Looking up gene:", gene, "\n")

gene_row <- anno[anno$gene_name == gene]

if (length(gene_row) == 0) {
  stop(paste("Gene not found in annotation:", gene))
}

seq <- as.character(seqnames(gene_row))[1]
start <- start(gene_row)[1]
end <- end(gene_row)[1]

region_string <- paste0(seq, "-", start, "-", end)

cat("Using region:", region_string, "\n")

# ---------------------------
# REMOVE ANNOTATION COMPLETELY
# ---------------------------
cat("Stripping annotation completely (fixing Signac bug)...\n")

empty_gr <- GRanges()   # empty genomics ranges
Annotation(obj) <- empty_gr

# ---------------------------
# Plot
# ---------------------------
cat("Plotting...\n")


p <- CoveragePlot(
  object = obj,
  region = region_string,
  peaks = TRUE,
  links = FALSE    # links need annotation, so OFF
)

# ---------------------------
# Save PNG
# ---------------------------
cat("Saving PNG:", args$out, "\n")

png(args$out, width=1600, height=800, res=150)
print(p)
dev.off()

cat("Done.\n")

