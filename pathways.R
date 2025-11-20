library(ChIPseeker)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
rds_file <- args[1]

cat("Loading:", rds_file, "\n")
myObject <- readRDS(rds_file)

# Get differential peaks
diff_peaks <- myObject@misc$sdiff_peaks
cat("Total differential peaks:", nrow(diff_peaks), "\n")

# Filter for upregulated peaks in TH2
up_peaks <- diff_peaks[diff_peaks$avg_log2FC > 0, ]
cat("Upregulated peaks in TH2:", nrow(up_peaks), "\n")

# Convert peak strings to GRanges manually
peak_strings <- rownames(up_peaks)
peak_list <- strsplit(peak_strings, "-")
chroms <- sapply(peak_list, `[`, 1)
starts <- as.numeric(sapply(peak_list, `[`, 2)) 
ends <- as.numeric(sapply(peak_list, `[`, 3))
up_granges <- GRanges(seqnames = chroms, ranges = IRanges(start = starts, end = ends))

# Annotate peaks with ChIPseeker
cat("Annotating peaks with ChIPseeker...\n")
annotated_peaks <- annotatePeak(
  up_granges,
  tssRegion = c(-3000, 3000),
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db"
)

# Get target genes
target_genes <- annotated_peaks@anno$SYMBOL
target_genes <- unique(target_genes[!is.na(target_genes)])

cat("Number of unique target genes near upregulated peaks:", length(target_genes), "\n")

# Save target genes
write.csv(data.frame(Gene = target_genes), "TH2_target_genes.csv", row.names=FALSE)
cat("Target genes saved to: TH2_target_genes.csv\n")

# Pathway analysis
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
write.csv(as.data.frame(ego), "TH2_pathways.csv")
cat("Pathway results saved to: TH2_pathways.csv\n")

# Create dotplot
p <- dotplot(ego, showCategory=20)
ggsave("TH2_pathways_dotplot.png", p, width=10, height=8)
cat("Dotplot saved to: TH2_pathways_dotplot.png\n")

cat("Analysis complete!\n")
