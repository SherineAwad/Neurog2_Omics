library(ChIPseeker)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# -----------------------------

# Command-line arguments

# -----------------------------

args <- commandArgs(trailingOnly = TRUE)
rds_file <- args[1]
direction <- as.numeric(args[2])  # 1 or -1

if (!(direction %in% c(1, -1))) stop("Direction must be 1 (UP) or -1 (DOWN)")

suffix <- ifelse(direction == 1, "_UP", "_DOWN")

cat("Loading:", rds_file, "\n")
myObject <- readRDS(rds_file)

# -----------------------------

# Get differential peaks

# -----------------------------

diff_peaks <- myObject@misc$sdiff_peaks
cat("Total differential peaks:", nrow(diff_peaks), "\n")

# Filter peaks based on direction

if (direction == 1) {
selected_peaks <- diff_peaks[diff_peaks$avg_log2FC > 0, ]
} else {
selected_peaks <- diff_peaks[diff_peaks$avg_log2FC < 0, ]
}
cat("Selected peaks:", nrow(selected_peaks), "\n")

# -----------------------------

# Convert peak strings to GRanges

# -----------------------------

peak_strings <- rownames(selected_peaks)
peak_list <- strsplit(peak_strings, "-")
chroms <- sapply(peak_list, `[`, 1)
starts <- as.numeric(sapply(peak_list, `[`, 2))
ends <- as.numeric(sapply(peak_list, `[`, 3))
selected_granges <- GRanges(seqnames = chroms, ranges = IRanges(start = starts, end = ends))

# -----------------------------

# Annotate peaks

# -----------------------------

cat("Annotating peaks with ChIPseeker...\n")
annotated_peaks <- annotatePeak(
selected_granges,
tssRegion = c(-3000, 3000),
TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
annoDb = "org.Mm.eg.db"
)

# -----------------------------

# Extract target genes

# -----------------------------

target_genes <- annotated_peaks@anno$SYMBOL
target_genes <- unique(target_genes[!is.na(target_genes)])
cat("Number of unique target genes near selected peaks:", length(target_genes), "\n")

# -----------------------------

# Save target genes

# -----------------------------

outfile_genes <- paste0("TH2_target_genes", suffix, ".csv")
write.csv(data.frame(Gene = target_genes), outfile_genes, row.names = FALSE)
cat("Target genes saved to:", outfile_genes, "\n")

# -----------------------------

# GO pathway analysis

# -----------------------------

cat("Running GO pathway analysis...\n")
ego <- enrichGO(
gene = target_genes,
OrgDb = org.Mm.eg.db,
keyType = "SYMBOL",
ont = "BP",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05
)

# Save pathway results

outfile_pathways <- paste0("TH2_pathways", suffix, ".csv")
write.csv(as.data.frame(ego), outfile_pathways)
cat("Pathway results saved to:", outfile_pathways, "\n")

# Dotplot

outfile_dotplot <- paste0("TH2_pathways_dotplot", suffix, ".png")
p <- dotplot(ego, showCategory = 20)
ggsave(outfile_dotplot, p, width = 10, height = 8)
cat("Dotplot saved to:", outfile_dotplot, "\n")

cat("Analysis complete!\n")

