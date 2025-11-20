#!/usr/bin/env Rscript

# -------------------------
# Install packages if missing
# -------------------------
if(!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse", repos="http://cran.us.r-project.org")
if(!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap", repos="http://cran.us.r-project.org")
if(!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos="http://cran.us.r-project.org")

# -------------------------
# Load libraries
# -------------------------
library(argparse)
library(dplyr)
library(pheatmap)

# -------------------------
# Command-line parser
# -------------------------
parser <- ArgumentParser(description="Plot top N differential peaks heatmap from CSVs")
parser$add_argument("-d", "--diff", required=TRUE, help="Differential peaks CSV")
parser$add_argument("-a", "--annot", required=TRUE, help="Annotated peaks CSV")
parser$add_argument("-n", "--topN", type="integer", default=30, help="Top N peaks to plot")
parser$add_argument("-o", "--out", default="top_peaks_heatmap.png", help="Output heatmap PNG")
args <- parser$parse_args()

# -------------------------
# Load CSVs
# -------------------------
diff <- read.csv(args$diff, row.names=1, stringsAsFactors=FALSE)
annot <- read.csv(args$annot, stringsAsFactors=FALSE)

# Ensure peak columns
#if(!"query_region" %in% colnames(annot)) stop("Annotated CSV must have 'query_region' column")
#annot$peak <- annot$query_region

if(!"X" %in% colnames(annot)) stop("Annotated CSV must have 'X' column")
annot$peak <- annot$X



# -------------------------
# Intersect peaks
# -------------------------
common_peaks <- intersect(rownames(diff), annot$peak)
if(length(common_peaks) == 0) stop("No matching peaks between differential and annotation CSVs!")

diff_sub <- diff[common_peaks, , drop=FALSE]
annot_sub <- annot[match(common_peaks, annot$peak), ]

# -------------------------
# Add nearest gene column
# -------------------------
diff_sub$nearest_gene <- annot_sub$gene_name

# -------------------------
# Select top N peaks by abs(avg_log2FC)
# -------------------------
diff_top <- diff_sub %>% 
  mutate(abs_fc = abs(avg_log2FC)) %>% 
  arrange(desc(abs_fc)) %>% 
  head(args$topN)

# -------------------------
# Prepare matrix for heatmap
# -------------------------
heatmap_matrix <- as.matrix(diff_top[, c("avg_log2FC")])
rownames(heatmap_matrix) <- diff_top$nearest_gene

# -------------------------
# Plot heatmap
# -------------------------
png(args$out, width=1200, height=1000, res=150)
pheatmap(
  heatmap_matrix,
  cluster_rows=TRUE,
  cluster_cols=FALSE,
  show_rownames=TRUE,
  show_colnames=FALSE,
  fontsize_row=9,
  main=paste("Top", args$topN, "Differential Peaks by Nearest Gene")
)
dev.off()
cat("Heatmap saved as:", args$out, "\n")

