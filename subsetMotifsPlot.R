library(ggplot2)
library(dplyr)
library(argparse)

# Set up argument parser
parser <- ArgumentParser()
parser$add_argument("--rds", required=TRUE, help="Input RDS file with motif data")
parser$add_argument("--output", required=TRUE, help="Output heatmap file")
parser$add_argument("--csv", required=TRUE, help="Output CSV file for all motifs")

args <- parser$parse_args()

cat("Loading:", args$rds, "\n")
motif_data <- readRDS(args$rds)

# Write CSV with motif gene names
write.csv(motif_data, args$csv, row.names = FALSE)
cat("All motifs with gene names written to:", args$csv, "\n")

# Create heatmap of top 20 motifs with motif names
top_motifs <- motif_data %>% 
  arrange(pvalue) %>% 
  head(20)

p <- ggplot(top_motifs, aes(x = "Motifs", y = reorder(motif.name, fold.enrichment))) +
  geom_tile(aes(fill = fold.enrichment), color = "white") +
  geom_text(aes(label = sprintf("%.2f", fold.enrichment)), size = 3) +
  scale_fill_viridis_c(name = "Fold Enrichment") +
  labs(
    title = "Top 20 Enriched Motifs in TH2 Upregulated Peaks",
    x = "",
    y = "Transcription Factor"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(args$output, p, width = 10, height = 8, dpi = 300)
cat("Heatmap saved to:", args$output, "\n")
