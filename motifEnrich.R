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
parser <- ArgumentParser(description = "Motif enrichment using custom genome + GTF")
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
atac_assay <- NULL
for (nm in names(obj@assays)) {
  if (inherits(obj@assays[[nm]], "ChromatinAssay"))
    atac_assay <- nm
}

if (is.null(atac_assay)) stop("No ChromatinAssay found in object.")
DefaultAssay(obj) <- atac_assay
message("Using ATAC assay: ", atac_assay)

# ============================================================
# 3. Load custom genome FASTA
# ============================================================
message("Loading custom FASTA genome: ", args$fa)
fa <- FaFile(args$fa)
indexFa(fa)   # creates .fai if missing

# ============================================================
# 4. Load custom GTF annotation
# ============================================================
message("Loading custom GTF: ", args$gtf)
gtf <- rtracklayer::import(args$gtf)
Annotation(obj) <- gtf
message("GTF annotation loaded.")

# ============================================================
# 5. Add motifs using JASPAR2020 + custom FASTA genome
# ============================================================
if (is.null(Motifs(obj))) {
  message("Motifs missing. Adding JASPAR 2020 motifs with custom genome...")

  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = 10090, all_versions = FALSE)
  )

  obj <- AddMotifs(
    object = obj,
    pfm = pfm,
    genome = fa    # <<-- CUSTOM FASTA genome
  )

} else {
  message("Motifs already exist — skipping AddMotifs().")
}

# ============================================================
# 6. Motif enrichment
# ============================================================
message("Running motif enrichment...")

motif_res <- FindMotifs(
  object = obj,
  features = rownames(obj),
  background = NULL
)

# ============================================================
# 7. Replace motif IDs with TF names
# ============================================================
motif_obj <- Motifs(obj)
motif_ids <- rownames(motif_obj)

tf_names <- sapply(motif_ids, function(x) {
  name <- motif_obj@motif.names[x]
  if (is.null(name)) x else name
})

motif_res$motif <- tf_names[motif_res$motif]

# ============================================================
# 8. Convert list columns to atomic vectors (CSV-safe)
# ============================================================
motif_res[] <- lapply(motif_res, function(x) if(is.list(x)) unlist(x) else x)

# Save results
csv_out <- file.path(args$outdir, "motif_enrichment.csv")
write.csv(motif_res, csv_out, row.names=FALSE)
message("Saved: ", csv_out)

# ============================================================
# 9. Plot top motifs
# ============================================================
top_n <- min(args$top, nrow(motif_res))
df <- motif_res[1:top_n, ]
df$motif <- factor(df$motif, levels=df$motif)

p <- ggplot(df, aes(x = motif, y = fold.enrichment)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top Enriched Motifs",
    x = "",
    y = "Fold Enrichment"
  )

png_out <- file.path(args$outdir, "motif_enrichment.png")
ggsave(png_out, p, width=7, height=6, dpi=300)
message("Saved: ", png_out)

message("=== Done ===")

