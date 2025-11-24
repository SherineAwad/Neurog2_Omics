library(Seurat)
library(dplyr)
library(ggplot2)

# Load object
mylist <- readRDS("Neurog2_Cells_withDGE.rds")
MyObject <- mylist$object
MyMarkers <- mylist$markers

# Load your gene list
genes_input <- scan("topegenes.txt", what = character(), quiet = TRUE)

# --------------------------
# 1. GET TOP 50 EXACTLY LIKE YOUR PIPELINE
# --------------------------
top50_pval <- MyMarkers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val, n = 50, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, p_val)

genes_top50 <- top50_pval$gene

# --------------------------
# 2. GET SCT MATRIX AND FILTER TO TOP 50
# --------------------------
DefaultAssay(MyObject) <- "SCT"
mat <- GetAssayData(MyObject, assay="SCT", slot="data")

genes_use <- intersect(genes_top50, rownames(mat))
if(length(genes_use) == 0) stop("No overlap between top50 genes and matrix rows.")

mat_sub <- mat[genes_use, , drop = FALSE]

# --------------------------
# 3. ROW-WISE Z-SCORE
# --------------------------
mat_z <- t(scale(t(mat_sub), center=TRUE, scale=TRUE))

MyObject <- SetAssayData(
  MyObject,
  assay = "SCT",
  layer = "scale.data",
  new.data = mat_z
)

# --------------------------
# 4. LABEL ONLY INPUT GENES - FIX OVERLAPPING
# --------------------------
genes_to_label <- intersect(genes_use, genes_input)

# Create labels with empty strings for non-input genes
gene_text <- ifelse(genes_use %in% genes_to_label, genes_use, "")

# --------------------------
# 5. PNG (600 DPI) WITH FIXED GENE LABELS
# --------------------------
png("Top50_pval.rna.Avg_SCT_heatmap_Zscore.png", width=10, height=14, units="in", res=600)

DoHeatmap(
  object = MyObject,
  features = genes_use,
  assay = "SCT",
  slot = "scale.data",
  draw.lines = FALSE,
  size = 5,
  raster = TRUE,
  group.colors = c("#026AB1","#61BFB9","#936DAD","#D11536","#AAA9A9","#EF9000")
) +
  theme(
    axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 1, margin = margin(r = 20)),
    axis.text.x = element_text(angle = 0)
  ) +
  scale_y_discrete(labels = gene_text) +
  scale_fill_gradient2(
    low  = rev(c('#D1E5F0', '#67A9CF', '#2166AC')),
    mid  = "white",
    high = rev(c('#B2182B', '#EF8A62', '#FDDBC7')),
    midpoint = 0,
    guide = "colourbar",
    aesthetics = "fill"
  )

dev.off()

cat("DONE: Top50_pval.rna.Avg_SCT_heatmap_Zscore.png\n")
