library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load object
mylist <- readRDS("Neurog2_Cells_withDGE.rds")
MyObject <- mylist$object
MyMarkers <- mylist$markers

# Load your gene list
genes_input <- scan("topgenes.txt", what = character(), quiet = TRUE)

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
mat <- GetAssayData(MyObject, assay="SCT", layer="data")  # Seurat 5+

genes_use <- genes_top50[genes_top50 %in% rownames(mat)]
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
# 4. LABEL ONLY INPUT GENES OUTSIDE
# --------------------------
gene_text <- sapply(genes_use, function(g) if(g %in% genes_input) g else "")

# --------------------------
# 5. PNG (600 DPI) WITH LABELS OUTSIDE ONLY - FIXED VERSION
# --------------------------
png("Top50_pval.rna.Avg_SCT_heatmap_Zscore.png", width=12, height=18, units="in", res=600)  # Increased height from 14 to 16

# Create base heatmap with NO y-axis labels
p <- DoHeatmap(
  object = MyObject,
  features = genes_use,
  assay = "SCT",
  slot = "scale.data",
  draw.lines = FALSE,
  size = 5,
  raster = TRUE,
  group.colors = c("#026AB1","#61BFB9","#936DAD","#D11536","#AAA9A9","#EF9000")
) +
  scale_y_discrete(labels = rep("", length(genes_use))) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.margin = margin(4, 0.5, 1, 2.8, "cm")  # Increased top margin to 3cm
  ) +
  scale_fill_gradient2(
    low  = rev(c('#D1E5F0', '#67A9CF', '#2166AC')),
    mid  = "white",
    high = rev(c('#B2182B', '#EF8A62', '#FDDBC7')),
    midpoint = 0
  )

# Prepare repel labels - FIXED POSITIONING
y_positions <- rev(seq_along(genes_use))  # REVERSED to match heatmap order
repel_df <- data.frame(
  x = -0.5,  # Position further left to be outside heatmap area
  y = y_positions,
  label = gene_text
)
repel_df <- repel_df[repel_df$label != "", ]

# Add geom_text_repel with adjusted parameters to place labels outside
p + geom_text_repel(
  data = repel_df,
  aes(x = x, y = y, label = label),
  inherit.aes = FALSE,
  direction = "y",
  hjust = 1,
  nudge_x = -0.3,  # Increased nudge to push labels further left
  segment.size = 0,
  segment.color = NA,
  force = 1,
  max.overlaps = Inf,
  xlim = c(-Inf, -0.1),  # Force labels to stay in left region
  box.padding = 0.3,     # Reduced slightly to fit better
  min.segment.length = 0,
  ylim = c(1, length(genes_use))  # Constrain labels to y-axis range
)

dev.off()

cat("DONE: Top50_pval.rna.Avg_SCT_heatmap_Zscore.png\n")
