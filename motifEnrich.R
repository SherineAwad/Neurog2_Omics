#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(Signac)
  library(JASPAR2020)
  library(TFBSTools)
  library(Rsamtools)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
})

# ============================================================
# 1. Command line arguments
# ============================================================
parser <- ArgumentParser(description = "Motif enrichment using custom genome + GTF (on differential peaks)")
parser$add_argument("--rds", required=TRUE, help="Input Seurat RDS")
parser$add_argument("--gtf", required=TRUE, help="Custom GTF file")
parser$add_argument("--fa",  required=TRUE, help="Custom genome FASTA")
parser$add_argument("--outdir", default="motif_out", help="Output directory")
parser$add_argument("--top", type="integer", default=20, help="Top motifs to plot")
args <- parser$parse_args()

dir.create(args$outdir, showWarnings=FALSE, recursive=TRUE)

# ============================================================
# 2. Load Seurat object
# ============================================================
message("Loading Seurat object...")
obj <- readRDS(args$rds)

# Detect ATAC assay
message("Detecting ATAC/ChromatinAssay...")
atac_assay <- names(obj@assays)[sapply(obj@assays, inherits, "ChromatinAssay")]
if(length(atac_assay)==0) stop("No ChromatinAssay found in object.")
DefaultAssay(obj) <- atac_assay
message("Using ATAC assay: ", atac_assay)

# ============================================================
# 3. Load custom FASTA genome
# ============================================================
message("Loading custom FASTA genome: ", args$fa)
fa <- FaFile(args$fa)
if(!file.exists(paste0(args$fa,".fai"))) indexFa(fa)

# ============================================================
# 4. Load custom GTF annotation
# ============================================================
message("Loading custom GTF: ", args$gtf)
gtf <- rtracklayer::import(args$gtf)
Annotation(obj) <- gtf
message("GTF annotation loaded.")

# ============================================================
# 5. Add motifs using JASPAR2020 + custom genome
# ============================================================
if(is.null(Motifs(obj))){
  message("Motifs missing. Adding JASPAR2020 motifs with custom genome...")
  pfm <- getMatrixSet(JASPAR2020, opts = list(species=10090, all_versions=FALSE))
  obj <- AddMotifs(obj, pfm=pfm, genome=fa)
} else {
  message("Motifs already exist. Skipping AddMotifs()")
}

# ============================================================
# 6. Select differential peaks
# ============================================================
# Option 1: If you have a column "diff_peaks" in metadata
if("diff_peaks" %in% colnames(obj@meta.data)){
  diff_peaks <- rownames(obj)[obj@meta.data$diff_peaks == TRUE]
} else {
  # Option 2: Use all peaks
  warning("No 'diff_peaks' column found. Using all peaks for motif enrichment.")
  diff_peaks <- rownames(obj)
}
message("Number of peaks for motif enrichment: ", length(diff_peaks))

# ============================================================
# 7. Run motif enrichment
# ============================================================
message("Running motif enrichment on differential peaks...")
motif_res <- FindMotifs(
  object = obj,
  features = diff_peaks,
  background = NULL
)

# ============================================================
# 8. Map motif IDs → TF names correctly
# ============================================================
motif_obj <- Motifs(obj)
motif_res$TF <- sapply(motif_res$motif, function(x){
  if(x %in% names(motif_obj@motif.names)){
    motif_obj@motif.names[[x]]
  } else {
    NA
  }
})

# Use TF column if available
motif_res$motif <- ifelse(is.na(motif_res$TF), motif_res$motif, motif_res$TF)
motif_res$TF <- NULL

# ============================================================
# 9. Save results
# ============================================================
motif_res[] <- lapply(motif_res, function(x) if(is.list(x)) unlist(x) else x)
csv_out <- file.path(args$outdir,"motif_enrichment.csv")
write.csv(motif_res, csv_out, row.names=FALSE)
message("Saved: ", csv_out)

# ============================================================
# 10. Plot top motifs
# ============================================================
top_n <- min(args$top, nrow(motif_res))
df <- motif_res[1:top_n, ]
df$motif <- factor(df$motif, levels=df$motif)

p <- ggplot(df, aes(x=motif, y=fold.enrichment)) +
  geom_col(fill="steelblue") +
  coord_flip() +
  theme_bw() +
  labs(title="Top Enriched Motifs (Differential Peaks)", x="", y="Fold Enrichment")

png_out <- file.path(args$outdir,"motif_enrichment.png")
ggsave(png_out, p, width=7, height=6, dpi=300)
message("Saved: ", png_out)

message("=== Done ===")

