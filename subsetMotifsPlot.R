#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)
library(dplyr)

parser <- ArgumentParser(description = 'Create heatmap of enriched motifs from data.frame')
parser$add_argument('--rds', required = TRUE, help = 'Path to RDS file containing motif data frame')
parser$add_argument('--output', default = 'motif_enrichment_heatmap.png', 
                   help = 'Output file name for the heatmap')
parser$add_argument('--csv', default = 'motif_enrichment_data.csv',
                   help = 'Output file name for the CSV data')
parser$add_argument('--top_n', type = 'integer', default = 30,
                   help = 'Number of top enriched motifs to plot')

args <- parser$parse_args()

# Load the data
motif_data <- readRDS(args$rds)

# Write all motif data to CSV
write.csv(motif_data, file = args$csv, row.names = FALSE)
cat("Full motif data saved to:", args$csv, "\n")

# Create heatmap
heatmap_data <- motif_data %>%
  arrange(desc(fold.enrichment)) %>%
  head(args$top_n) %>%
  mutate(motif.name = factor(motif.name, levels = rev(motif.name)))

p <- ggplot(heatmap_data, aes(x = "Enrichment", y = motif.name)) +
  geom_tile(aes(fill = fold.enrichment), color = "white") +
  geom_text(aes(label = sprintf("%.2f", fold.enrichment)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1,
                      name = "Fold Enrichment") +
  labs(title = paste("Top", args$top_n, "Enriched Motifs"),
       x = "", y = "Transcription Factor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(args$output, p, width = 10, height = 8, dpi = 300)
cat("Heatmap saved to:", args$output, "\n")
