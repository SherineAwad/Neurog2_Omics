#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(UpSetR)  # for UpSet plot
set.seed(1234)


library(future)

plan(sequential)   # stays the same (or 'multicore' if you want parallel)

# 80 GB limit
options(future.globals.maxSize = 80 * 1024^3)


args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

# --------------------------------------------------------------
# LOAD ANNOTATED OBJECT
# --------------------------------------------------------------
input_rds <- paste0(mysample, "_annotated.rds")
myObject <- readRDS(input_rds)

# --------------------------------------------------------------
# SETUP AND QUALITY CONTROL
# --------------------------------------------------------------
cat("=== DIFFERENTIAL GENE EXPRESSION ANALYSIS ===\n")
cat("Sample:", mysample, "\n")
cat("Number of cells:", ncol(myObject), "\n")
cat("Number of genes:", nrow(myObject), "\n")

celltypes <- levels(Idents(myObject))
cat("Available cell types:", paste(celltypes, collapse = ", "), "\n")

DefaultAssay(myObject) <- "RNA"

if ("SCT" %in% names(myObject@assays)) {
  cat("Using existing SCT assay...\n")
  DefaultAssay(myObject) <- "SCT"
  myObject <- PrepSCTFindMarkers(myObject)
} else {
  cat("Performing SCTransform normalization...\n")
  myObject <- SCTransform(myObject, verbose = FALSE)
  myObject <- PrepSCTFindMarkers(myObject)
}

# --------------------------------------------------------------
# DIFFERENTIAL EXPRESSION BETWEEN CELL TYPES
# --------------------------------------------------------------
cat("\n=== FINDING MARKERS BETWEEN CELL TYPES ===\n")

celltype_markers <- FindAllMarkers(
  myObject,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,      
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

write.csv(celltype_markers, paste0(mysample, "_all_celltype_markers.csv"), row.names = FALSE)

top_markers <- celltype_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

cat("\nTop markers per cell type:\n")
print(top_markers)

# --------------------------------------------------------------
# VISUALIZATION
# --------------------------------------------------------------
cat("\n=== CREATING VISUALIZATIONS ===\n")

top10_markers <- celltype_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

available_genes <- rownames(myObject@assays$SCT@scale.data)
top10_markers_filtered <- top10_markers %>% filter(gene %in% available_genes)

if (nrow(top10_markers_filtered) > 0) {
  # Heatmap
  Idents(myObject) <- factor(Idents(myObject), levels = c('MG', 'MGPC', 'BC', 'AC', 'Rod', 'Cones'))

  heatmap_plot <- DoHeatmap(
    myObject,
    features = top10_markers_filtered$gene,
    assay = "SCT",
    label = TRUE,
    size = 4,
    cells = WhichCells(myObject, idents = levels(Idents(myObject)))
  ) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme(axis.text.y = element_text(size = 10))

  png(paste0(mysample, "_celltype_markers_heatmap.png"), width = 1200, height = 1000)
  print(heatmap_plot)
  dev.off()
} else {
  cat("No genes found in scale.data for heatmap\n")
}

# Dot plot
dot_plot <- DotPlot(
  myObject,
  features = unique(top10_markers$gene),
  assay = "SCT",
  cols = c("blue", "red"),
  dot.scale = 6
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

png(paste0(mysample, "_celltype_markers_dotplot.png"), width = 1200, height = 800)
print(dot_plot)
dev.off()

# Feature plots (robust to 1 gene)
top5_genes <- top_markers$gene[1:min(9, nrow(top_markers))]

feature_plots <- FeaturePlot(
  myObject,
  features = top5_genes,
  reduction = "umap.wnn.harmony",
  ncol = 3,
  pt.size = 0.1,
  order = TRUE,
  combine = FALSE
)

png(paste0(mysample, "_celltype_top_markers_featureplot.png"), width = 1500, height = 1200)
if(length(feature_plots) == 1){
  print(feature_plots[[1]])
} else {
  print(wrap_plots(feature_plots, ncol = 3))
}
dev.off()

# --------------------------------------------------------------
# SAVE FINAL OBJECT
# --------------------------------------------------------------
cat("\n=== SAVING RESULTS ===\n")
myObject@misc$celltype_markers <- celltype_markers
myObject@misc$top_markers <- top_markers
saveRDS(myObject, file = paste0(mysample, "_with_DGE.rds"))

# --------------------------------------------------------------
# FINAL OUTPUT
# --------------------------------------------------------------
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created in current directory:\n")
cat("- ", mysample, "_all_celltype_markers.csv\n")
cat("- ", mysample, "_celltype_markers_heatmap.png\n")
cat("- ", mysample, "_celltype_markers_dotplot.png\n")
cat("- ", mysample, "_celltype_top_markers_featureplot.png\n")
cat("- ", mysample, "_with_DGE.rds\n")
cat("Total significant marker genes (p_adj < 0.05):", sum(celltype_markers$p_val_adj < 0.05), "\n")

