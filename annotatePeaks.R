#!/usr/bin/env Rscript

# =======================================
# Annotate differential peaks with nearest genes (argparse version)
# =======================================

if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse", repos = "http://cran.us.r-project.org")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "http://cran.us.r-project.org")

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
# Load differential peaks CSV
# -------------------------
dar <- read.csv(args$file, stringsAsFactors = FALSE)

# Your file has peaks in the first column (the one with empty header)
# The first column contains coordinates like "chr7-19628475-19629074"
peak_col <- 1  # use first column
peak_coords <- dar[[peak_col]]

# Remove quotes if present
peak_coords <- gsub('"', '', peak_coords)

# -------------------------
# Split coordinates safely
# -------------------------
coords <- strsplit(as.character(peak_coords), "-")

# Extract seqnames, start, end as plain vectors and remove malformed entries
seqnames_vec <- unlist(lapply(coords, `[`, 1))
start_vec    <- as.numeric(unlist(lapply(coords, `[`, 2)))
end_vec      <- as.numeric(unlist(lapply(coords, `[`, 3)))

valid_idx <- !is.na(seqnames_vec) & !is.na(start_vec) & !is.na(end_vec)
seqnames_vec <- seqnames_vec[valid_idx]
start_vec    <- start_vec[valid_idx]
end_vec      <- end_vec[valid_idx]
dar <- dar[valid_idx, ]
peak_coords <- peak_coords[valid_idx]

# -------------------------
# Create GRanges
# -------------------------
peak_gr <- GRanges(
  seqnames = seqnames_vec,
  ranges   = IRanges(start = start_vec, end = end_vec),
  peak_name = peak_coords
)

# -------------------------
# Load GTF
# -------------------------
gtf <- import(args$gtf)
genes <- gtf[gtf$type == "gene"]

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
