library(Signac)
library(Seurat)
library(pheatmap)
library(Matrix)
library(dplyr)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

# Load RDS
rds_file <- paste0(mysample, "_MG_MGPC_diffPeaks.rds")
obj <- readRDS(rds_file)

# Extract differential peaks and annotated peaks
diff_peaks <- obj@misc$sdiff_peaks
annot <- obj@misc$sannotated_diff_peaks

# Keep only valid genomic regions
annot <- annot[grepl("^chr[0-9XYM]+-[0-9]+-[0-9]+$", annot$query_region), ]
peak_coords <- annot$query_region

# ATAC matrix
DefaultAssay(obj) <- "ATAC"
mat_all <- GetAssayData(obj, assay = "ATAC", layer = "data")
rownames(mat_all) <- gsub("\\s+", "", rownames(mat_all))

# Match annotation to matrix
matched_peaks <- intersect(peak_coords, rownames(mat_all))
if(length(matched_peaks) == 0) stop("No peaks match between annotation and ATAC matrix! Check formats.")

# Top N peaks
top_n <- 30
if(!is.null(diff_peaks) && "p_val_adj" %in% colnames(diff_peaks)) {
  split_peaks <- do.call(rbind, strsplit(rownames(diff_peaks), "[-:]"))
  rownames(diff_peaks) <- paste0(split_peaks[,1], "-", split_peaks[,2], "-", split_peaks[,3])
  ordered_peaks <- rownames(diff_peaks)[order(diff_peaks$p_val_adj)]
  top_peaks <- head(intersect(ordered_peaks, matched_peaks), top_n)
} else {
  top_peaks <- head(matched_peaks, top_n)
}

# Subset matrix
mat <- mat_all[top_peaks, ]

# Map rownames to gene names
gene_labels <- annot$gene_name
names(gene_labels) <- annot$query_region
rownames(mat) <- gene_labels[top_peaks]

# Convert to dense numeric
mat_dense <- as.matrix(mat)

# Compute row-wise z-score properly with numeric precision
mat_z <- t(apply(mat_dense, 1, function(x) {
  x <- as.numeric(x)
  if(sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)
}))
rownames(mat_z) <- rownames(mat)

# OPTIONAL: scale to visually noticeable range (-3 to 3)
mat_z[mat_z > 3] <- 3
mat_z[mat_z < -3] <- -3

# Plot heatmap without clustering
heatmap_png <- paste0(mysample, "_MG_MGPC_diffPeaksheatmap.png")
png(heatmap_png, width = 1600, height = 1400, res = 150)
pheatmap(
  mat_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10,
  fontsize_col = 8,
  scale = "none",
  color = colorRampPalette(c("blue","white","red"))(100)
)
dev.off()
message("Saved heatmap: ", heatmap_png)
message("DONE.")

