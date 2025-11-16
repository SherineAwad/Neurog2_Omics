#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

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

# The cell types are already in Idents from the annotation script
celltypes <- levels(Idents(myObject))
cat("Available cell types:", paste(celltypes, collapse = ", "), "\n")

# Set default assay to RNA for DGE
DefaultAssay(myObject) <- "RNA"

# Check if SCT assay exists and prepare it properly
if ("SCT" %in% names(myObject@assays)) {
  cat("Using existing SCT assay...\n")
  DefaultAssay(myObject) <- "SCT"
  # Prepare SCT for FindMarkers
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

# Find all markers for each cell type
celltype_markers <- FindAllMarkers(
  myObject,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

# Save all markers
write.csv(celltype_markers, paste0(mysample, "_all_celltype_markers.csv"), row.names = FALSE)

# Get top markers per cell type
top_markers <- celltype_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

cat("\nTop markers per cell type:\n")
print(top_markers)

# --------------------------------------------------------------
# VISUALIZATION
# --------------------------------------------------------------
cat("\n=== CREATING VISUALIZATIONS ===\n")

# 1. Heatmap of top markers
top10_markers <- celltype_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Get genes that exist in scale.data
available_genes <- rownames(myObject@assays$SCT@scale.data)
top10_markers_filtered <- top10_markers %>% filter(gene %in% available_genes)

if (nrow(top10_markers_filtered) > 0) {
  heatmap_plot <- DoHeatmap(
    myObject,
    features = top10_markers_filtered$gene,
    assay = "SCT",
    label = TRUE,
    size = 3
  ) + theme(axis.text.y = element_text(size = 6))
  
  png(paste0(mysample, "_celltype_markers_heatmap.png"), width = 1200, height = 1000)
  print(heatmap_plot)
  dev.off()
} else {
  cat("No genes found in scale.data for heatmap\n")
}

# 2. Dot plot of top markers
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

# 3. Feature plots for top markers on umap.wnn.harmony
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

feature_grid <- wrap_plots(feature_plots, ncol = 3)
png(paste0(mysample, "_celltype_top_markers_featureplot.png"), width = 1500, height = 1200)
print(feature_grid)
dev.off()

# 4. Violin plots for top markers per cell type
for (celltype in unique(top_markers$cluster)) {
  celltype_genes <- top_markers %>% 
    filter(cluster == celltype) %>%
    pull(gene) %>%
    head(3)
  
  if (length(celltype_genes) > 0) {
    vln_plot <- VlnPlot(
      myObject,
      features = celltype_genes,
      idents = celltype,
      ncol = 3,
      pt.size = 0,
      assay = "SCT"
    )
    
    png(paste0(mysample, "_", celltype, "_markers_violin.png"), width = 1200, height = 400)
    print(vln_plot)
    dev.off()
  }
}

# --------------------------------------------------------------
# SUMMARY STATISTICS
# --------------------------------------------------------------
cat("\n=== GENERATING SUMMARY ===\n")

summary_stats <- celltype_markers %>%
  group_by(cluster) %>%
  summarise(
    n_genes = n(),
    n_sig_genes = sum(p_val_adj < 0.05),
    n_sig_up = sum(p_val_adj < 0.05 & avg_log2FC > 0),
    top_gene = first(gene),
    top_log2FC = first(avg_log2FC)
  )

write.csv(summary_stats, paste0(mysample, "_celltype_DGE_summary.csv"))

# --------------------------------------------------------------
# SAVE FINAL OBJECT
# --------------------------------------------------------------
cat("\n=== SAVING RESULTS ===\n")

myObject@misc$celltype_markers <- celltype_markers
myObject@misc$top_markers <- top_markers
myObject@misc$dge_summary <- summary_stats

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
cat("- ", mysample, "_celltype_DGE_summary.csv\n")
cat("- ", mysample, "_with_DGE.rds\n")
cat("Total significant marker genes (p_adj < 0.05):", sum(celltype_markers$p_val_adj < 0.05), "\n")
cat("Average markers per cell type:", round(mean(summary_stats$n_genes), 1), "\n")
