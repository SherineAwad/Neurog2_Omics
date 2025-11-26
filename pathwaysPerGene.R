library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

# -----------------------------

# Command-line arguments

# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

# Default values

gene_file <- NULL
prefix <- "GO"

# Parse arguments manually

for (i in seq(1, length(args), by = 2)) {
if (args[i] == "-i") gene_file <- args[i + 1]
if (args[i] == "-o") prefix <- args[i + 1]
}

if (is.null(gene_file)) stop("Usage: Rscript script.R -i gene_file.csv -o prefix")

cat("Loading gene list from:", gene_file, "\n")
cat("Output prefix:", prefix, "\n")

# -----------------------------

# Read gene list

# -----------------------------

genes_df <- read.csv(gene_file, stringsAsFactors = FALSE)
if (!"Gene" %in% colnames(genes_df)) stop("CSV must have a 'Gene' column")

gene_list <- unique(genes_df$Gene)
cat("Number of unique genes:", length(gene_list), "\n")

# -----------------------------

# GO pathway analysis

# -----------------------------

cat("Running GO enrichment analysis...\n")
ego <- enrichGO(
gene = gene_list,
OrgDb = org.Mm.eg.db,
keyType = "SYMBOL",
ont = "BP",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05
)

# -----------------------------

# Save pathway results

# -----------------------------

outfile_pathways <- paste0(prefix, "_pathways.csv")
write.csv(as.data.frame(ego), outfile_pathways, row.names = FALSE)
cat("Pathway results saved to:", outfile_pathways, "\n")

# -----------------------------

# Dotplot

# -----------------------------

outfile_dotplot <- paste0(prefix, "_dotplot.png")
p <- dotplot(ego, showCategory = 20)
ggsave(outfile_dotplot, p, width = 10, height = 8)
cat("Dotplot saved to:", outfile_dotplot, "\n")

cat("Analysis complete!\n")

