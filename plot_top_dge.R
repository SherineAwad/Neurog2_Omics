#!/usr/bin/env Rscript

# -------------------------
# Load required libraries (assume already installed)
# -------------------------
suppressMessages(library(argparse))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(pheatmap))
suppressMessages(library(readr))
suppressMessages(library(tibble))

# -------------------------
# Command-line parser
# -------------------------
parser <- ArgumentParser(description="Plot top N genes heatmap from CSV")
parser$add_argument("-f", "--file", required=TRUE, help="Input CSV file")
parser$add_argument("-n", "--topN", type="integer", default=20, help="Top N genes to plot")
parser$add_argument("-o", "--out", default="top_genes_heatmap.png", help="Output PNG filename")
args <- parser$parse_args()

# -------------------------
# Load CSV
# -------------------------
df <- read_csv(args$file, col_types = cols(), show_col_types = FALSE)

# Rename first column to 'gene' if it's unnamed (starts with ... or is empty)
if(names(df)[1] == "" || grepl("^\\.\\.\\.", names(df)[1])) {
  names(df)[1] <- "gene"
}

# Check required columns
required_cols <- c("gene", "p_val_adj", "avg_log2FC")
missing_cols <- setdiff(required_cols, names(df))
if(length(missing_cols) > 0){
  stop("CSV must have columns: ", paste(required_cols, collapse=", "), 
       "\nFound columns: ", paste(names(df), collapse=", "))
}

# -------------------------
# Select top N genes by avg_log2FC (absolute value)
# -------------------------
top_genes <- df %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n=args$topN)
# -------------------------
# Prepare heatmap matrix - USE RAW avg_log2FC VALUES
# -------------------------
heatmap_mat <- top_genes %>%
  select(gene, avg_log2FC) %>%
  column_to_rownames(var="gene") %>%
  as.matrix()

# -------------------------
# Plot heatmap with raw values
# -------------------------
if(nrow(heatmap_mat) < 1){
  cat("Heatmap skipped: no genes to plot.\n")
} else {
  cat("Plotting heatmap with", nrow(heatmap_mat), "genes...\n")
  png(args$out, width=1200, height=1000, res=150)
  pheatmap(
    heatmap_mat,
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    show_rownames=TRUE,
    show_colnames=TRUE,
    fontsize_row=9,
    color=colorRampPalette(c("blue","white","red"))(100),
    main=paste("Top", args$topN, "Genes by Adjusted P-value"),
    display_numbers = FALSE
  )
  dev.off()
  cat("Heatmap saved to:", args$out, "\n")
}
