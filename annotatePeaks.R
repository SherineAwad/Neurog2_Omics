#!/usr/bin/env Rscript

# =======================================
# Annotate differential peaks with nearest genes (argparse version)
# =======================================

# -------------------------
# Install packages if needed
# -------------------------
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse", repos = "http://cran.us.r-project.org")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "http://cran.us.r-project.org")

# -------------------------
# Load libraries
# -------------------------
library(argparse)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

# -------------------------
# Command-line arguments
# -------------------------
parser <- ArgumentParser(description = "Annotate peaks with nearest gene")
parser$add_argument("-f", "--file", required = TRUE, help = "Input differential peaks CSV")
parser$add_argument("-g", "--gtf", required = TRUE, help = "Input genome GTF file")
parser$add_argument("-o", "--out", default = "annotated_peaks.csv", help = "Output CSV filename")
args <- parser$parse_args()

# -------------------------
# Load differential peaks
# -------------------------
dar <- read.csv(args$file, row.names = 1, stringsAsFactors = FALSE)
dar$peak <- rownames(dar)

# -------------------------
# Convert peak coordinates to GRanges
# -------------------------
coords <- strsplit(dar$peak, "-")
peak_gr <- GRanges(
  seqnames = sapply(coords, `[`, 1),
  ranges   = IRanges(
    start = as.numeric(sapply(coords, `[`, 2)),
    end   = as.numeric(sapply(coords, `[`, 3))
  ),
  peak_name = dar$peak
)

# -------------------------
# Load GTF
# -------------------------
gtf <- import(args$gtf)
genes <- gtf[gtf$type == "gene"]

# Extract gene metadata
gene_names <- mcols(genes)$gene_name
gene_types <- mcols(genes)$gene_type

# -------------------------
# Find nearest genes
# -------------------------
nearest_idx <- nearest(peak_gr, genes)
dar$nearest_gene <- gene_names[nearest_idx]
dar$gene_type    <- gene_types[nearest_idx]
dar$distance     <- distance(peak_gr, genes[nearest_idx])

# -------------------------
# Save annotated CSV
# -------------------------
write.csv(dar, args$out, row.names = FALSE)
cat("Annotated peaks saved to:", args$out, "\n")

