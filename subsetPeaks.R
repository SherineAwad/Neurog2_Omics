library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)

plan(sequential)
options(future.globals.maxSize = 80 * 1024^3)

# ----------------------------
# Input
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

myRDS <- paste0(mysample, "_MG_MGPC_subset.rds")
message("Loading: ", myRDS)
obj <- readRDS(myRDS)

# ----------------------------
# Set default assay
# ----------------------------
DefaultAssay(obj) <- "ATAC"

# ----------------------------
# Set cell identities to Sample
# ----------------------------
Idents(obj) <- obj$Sample

# ----------------------------
# Differential peaks: TH2 vs TH1 samples
# ----------------------------
message("Running differential peak analysis TH2 vs TH1...")

diff_peaks <- FindMarkers(
  object = obj,
  ident.1 = "TH2",
  ident.2 = "TH1",
  assay = "ATAC",
  test.use = "wilcox",
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0.1
)

message("Differential peaks found: ", nrow(diff_peaks))

# ----------------------------
# Save diff peaks
# ----------------------------
out_csv <- paste0(mysample, "_MG_MGPC_diffPeaks.csv")
write.csv(diff_peaks, out_csv, row.names = TRUE)
message("Saved differential peaks: ", out_csv)

# ----------------------------
# Annotate peaks
# ----------------------------
message("Annotating peaks with nearest genes...")

peak_ranges <- StringToGRanges(rownames(diff_peaks), sep = c("-", "-"))
nearest <- ClosestFeature(obj, peak_ranges)

annot_csv <- paste0(mysample, "_MG_MGPC_annotatedDiffPeaks.csv")
write.csv(nearest, annot_csv, row.names = FALSE)
message("Saved annotated peaks: ", annot_csv)

# ----------------------------
# Store in @misc for downstream steps
# ----------------------------
obj@misc$sdiff_peaks <- diff_peaks
obj@misc$sannotated_diff_peaks <- nearest

# Save RDS
out_rds <- paste0(mysample, "_MG_MGPC_diffPeaks.rds")
saveRDS(obj, out_rds)
message("Saved updated object: ", out_rds)

