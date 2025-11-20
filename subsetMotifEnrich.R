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

# Get peak coordinates from annotated peaks
annotated_peaks <- myObject@misc$sannotated_diff_peaks
if (is.null(annotated_peaks)) stop("sannotated_diff_peaks missing")

peak_ids <- as.character(annotated_peaks$query_region)
cat("Number of differential peaks:", length(peak_ids), "\n")

valid_peaks <- intersect(peak_ids, rownames(myObject))
cat("Valid peaks found in object:", length(valid_peaks), "\n")

if (length(valid_peaks) == 0) stop("No valid peaks found")

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

# Run motif enrichment - FIX: use manual background to avoid GC matching
cat("Running motif enrichment...\n")

# Get all peaks in the object as background
all_peaks <- rownames(myObject)
background_peaks <- setdiff(all_peaks, valid_peaks)

enriched <- FindMotifs(
  object = myObject,
  features = valid_peaks,
  background = background_peaks  # Explicitly set background peaks
)

# Save output
outfile <- paste0(mysample, "_MG_MGPC_motif.rds")
saveRDS(enriched, outfile)
cat("Saved:", outfile, "\n")
