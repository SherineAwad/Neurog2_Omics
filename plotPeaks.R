#!/usr/bin/env Rscript

# =========================
# Top N Peaks per Cell Type using annotated CSV (fixed for duplicate gene labels)
# =========================

# -------------------------
# Install packages if needed
# -------------------------
if (!requireNamespace("argparse", quietly = TRUE)) install.packages("argparse", repos = "http://cran.us.r-project.org")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")

# -------------------------
# Load libraries
# -------------------------
library(argparse)
library(dplyr)
library(ggplot2)

# -------------------------
# Command-line parser
# -------------------------
parser <- ArgumentParser(description = "Plot top N differential peaks per cell type using annotated CSV")
parser$add_argument("-f", "--file", required = TRUE, help = "Input annotated CSV file with differential peaks")
parser$add_argument("-n", "--topN", type = "integer", default = 10, help = "Top N peaks per cell type")
parser$add_argument("-o", "--out", default = "top_peaks.png", help = "Output PNG filename")
args <- parser$parse_args()

# -------------------------
# Read input CSV
# -------------------------
dar <- read.csv(args$file, stringsAsFactors = FALSE)

# Use gene names for plotting; fallback to peak coordinates if gene missing
dar$gene_label <- ifelse(!is.na(dar$nearest_gene) & dar$nearest_gene != "", dar$nearest_gene, dar$peak)

# -------------------------
# Get top N peaks per cluster
# -------------------------
top_peaks <- dar %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = args$topN) %>%
  ungroup()

# -------------------------
# Reorder peaks within each cluster safely
# -------------------------
top_peaks <- top_peaks %>%
  group_by(cluster) %>%
  mutate(
    # Make unique labels by appending index for duplicates
    peak_order = make.unique(as.character(gene_label))
  ) %>%
  mutate(peak_order = factor(peak_order, levels = peak_order)) %>%
  ungroup()

# -------------------------
# Plot barplot
# -------------------------
p <- ggplot(top_peaks, aes(x = peak_order, y = avg_log2FC, fill = cluster)) +
  geom_col() +
  facet_wrap(~cluster, scales = "free_y") +
  coord_flip() +
  labs(title = paste("Top", args$topN, "Differential Peaks per Cell Type"),
       x = "Gene / Peak", y = "Average log2 Fold Change") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# -------------------------
# Save plot as PNG (high DPI)
# -------------------------
ggsave(filename = args$out, plot = p, width = 12, height = 8, dpi = 300)

cat("Top peaks plot saved as:", args$out, "\n")

