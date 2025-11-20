library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)

set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]
myRDS <- paste0(mysample, "_MG_MGPC_diffPeaks.rds")

cat("Loading:", myRDS, "\n")
myObject <- readRDS(myRDS)

# Get differential peaks results
diff_peaks <- myObject@misc$sdiff_peaks
if (is.null(diff_peaks)) stop("sdiff_peaks missing")

cat("Columns in sdiff_peaks:\n")
print(colnames(diff_peaks))

# Since differential peaks were called as TH2 vs TH1,
# positive logFC means upregulated in TH2
# Look for the logFC column (could be avg_log2FC, logFC, etc.)
if ("avg_log2FC" %in% colnames(diff_peaks)) {
  fc_col <- "avg_log2FC"
} else if ("logFC" %in% colnames(diff_peaks)) {
  fc_col <- "logFC"
} else if ("log2FoldChange" %in% colnames(diff_peaks)) {
  fc_col <- "log2FoldChange"
} else {
  stop("Cannot find fold change column in sdiff_peaks")
}

cat("Using fold change column:", fc_col, "\n")

# Filter for upregulated peaks in TH2 (positive logFC)
up_peaks <- diff_peaks[diff_peaks[[fc_col]] > 0, ]
peak_ids <- rownames(up_peaks)

cat("Number of UPREGULATED peaks in TH2:", length(peak_ids), "\n")
cat("Fold change range in upregulated peaks:", range(up_peaks[[fc_col]]), "\n")

valid_peaks <- intersect(peak_ids, rownames(myObject))
cat("Valid UPREGULATED peaks found in object:", length(valid_peaks), "\n")

if (length(valid_peaks) == 0) stop("No valid UPREGULATED peaks found")

# Add motifs
cat("Adding motifs...\n")
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 10090, all_versions = FALSE)
)

myObject <- AddMotifs(
  object = myObject,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

# Run motif enrichment on UPREGULATED peaks
cat("Running motif enrichment on UPREGULATED peaks...\n")

# Check if GC content is available
gc_content <- myObject@assays$ATAC@meta.features$gc.percent

if (!is.null(gc_content) && any(!is.na(gc_content))) {
  # Use peaks with GC content if available
  peaks_with_gc <- rownames(myObject)[!is.na(gc_content)]
  valid_peaks_with_gc <- intersect(valid_peaks, peaks_with_gc)
  cat("Valid UPREGULATED peaks with GC content:", length(valid_peaks_with_gc), "\n")
  
  if (length(valid_peaks_with_gc) > 0) {
    final_peaks <- valid_peaks_with_gc
    cat("Running motif enrichment WITH GC content correction\n")
  } else {
    cat("Warning: No valid peaks with GC content. Running without GC correction.\n")
    final_peaks <- valid_peaks
  }
} else {
  cat("Warning: No GC content data available. Running motif enrichment without GC correction.\n")
  final_peaks <- valid_peaks
}

# Run motif enrichment
enriched <- FindMotifs(
  object = myObject,
  features = final_peaks
)

# Save output
outfile <- paste0(mysample, "_MG_MGPC_motif_UPREGULATED.rds")
saveRDS(enriched, outfile)
cat("Saved:", outfile, "\n")
cat("Motif enrichment completed successfully!\n")
