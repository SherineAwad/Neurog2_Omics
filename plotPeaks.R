#!/usr/bin/env Rscript

# =========================
# Top N Peaks per Cell Type + Peak Category Summary + Cluster-specific Categories
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
parser <- ArgumentParser(description = "Plot top N differential peaks per cell type using annotated CSV and peak category summary")
parser$add_argument("-f", "--file", required = TRUE, help = "Input annotated CSV file with differential peaks")
parser$add_argument("-n", "--topN", type = "integer", default = 10, help = "Top N peaks per cell type")
parser$add_argument("-p", "--peaks", default = "top_peaks.png", help = "Output PNG filename for top peaks plot")
parser$add_argument("-c", "--category", default = "peak_category.png", help = "Output PNG filename for peak category plot")
parser$add_argument("-b", "--bycluster", default = "peak_category_by_cluster.png", help = "Output PNG filename for cluster-specific peak category plot")
args <- parser$parse_args()

# -------------------------
# Read input CSV
# -------------------------
dar <- read.csv(args$file, stringsAsFactors = FALSE)
colnames(dar) <- trimws(colnames(dar))  # trim column names

# Use nearest_gene if exists, fallback to peak
dar$gene_label <- ifelse(!is.na(dar$nearest_gene) & dar$nearest_gene != "", dar$nearest_gene, dar$peak)

# -------------------------
# Get top N peaks per cluster
# -------------------------
top_peaks <- dar %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = args$topN) %>%
  ungroup()

# Reorder peaks safely
top_peaks <- top_peaks %>%
  group_by(cluster) %>%
  mutate(
    peak_order = make.unique(as.character(gene_label))
  ) %>%
  mutate(peak_order = factor(peak_order, levels = peak_order)) %>%
  ungroup()

# -------------------------
# Plot top peaks barplot
# -------------------------
p <- ggplot(top_peaks, aes(x = peak_order, y = avg_log2FC, fill = cluster)) +
  geom_col() +
  facet_wrap(~cluster, scales = "free_y") +
  coord_flip() +
  labs(title = paste("Top", args$topN, "Differential Peaks per Cell Type"),
       x = "Gene / Peak", y = "Average log2 Fold Change") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave(filename = args$peaks, plot = p, width = 12, height = 8, dpi = 300)
cat("Top peaks plot saved as:", args$peaks, "\n")

# -------------------------
# Create Peak Category (Promoter/Enhancer/Distal)
# -------------------------
dar$category <- cut(
  dar$distance,
  breaks = c(-Inf, 2000, 10000, Inf),
  labels = c("Promoter (<2kb)", "Enhancer (2-10kb)", "Distal (>10kb)")
)

# -------------------------
# Plot peak category pie chart (all peaks)
# -------------------------
category_counts <- dar %>%
  group_by(category) %>%
  summarise(count = n()) %>%
  ungroup()

p_cat <- ggplot(category_counts, aes(x = "", y = count, fill = category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  labs(title = "Peak Category Distribution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right")

ggsave(filename = args$category, plot = p_cat, width = 6, height = 6, dpi = 300)
cat("Peak category plot saved as:", args$category, "\n")

# -------------------------
# Cluster-specific peak category summary (stacked bar)
# -------------------------
cluster_summary <- dar %>%
  group_by(cluster, category) %>%
  summarise(count = n()) %>%
  ungroup()

p_cluster <- ggplot(cluster_summary, aes(x = cluster, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(title = "Peak Category Distribution per Cluster",
       x = "Cluster", y = "Peak Count") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

ggsave(filename = args$bycluster, plot = p_cluster, width = 10, height = 6, dpi = 300)
cat("Cluster-specific peak category plot saved as:", args$bycluster, "\n")

