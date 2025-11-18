library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]
annotation_file <- args[2]

# Read annotation file
if (!file.exists(annotation_file)) {
  stop(paste("Annotation file not found:", annotation_file))
}

annotations <- read.csv(annotation_file, stringsAsFactors = FALSE)

# Check if required columns exist
if (!all(c("cluster", "cell_type") %in% colnames(annotations))) {
  stop("Annotation file must contain 'cluster' and 'cell_type' columns")
}

myRDS <- paste(mysample, "_reClustered.rds", sep="")
cat("Processing file:", myRDS, "\n")

if (!file.exists(myRDS)) {
  stop(paste("Input RDS file not found:", myRDS))
}

myObject <- readRDS(myRDS)

# Check available reductions
cat("Available reductions in the object:\n")
print(names(myObject@reductions))

# Create named vector for RenameIdents from annotation file
rename_vector <- setNames(annotations$cell_type, as.character(annotations$cluster))

# Apply the renaming
myObject <- RenameIdents(object = myObject, rename_vector)

cat("Cluster renaming completed. New identities:\n")
print(table(Idents(myObject)))

# Determine which reduction to use for plotting
reduction_to_use <- NULL
if ("umap.wnn.harmony" %in% names(myObject@reductions)) {
  reduction_to_use <- "umap.wnn.harmony"
  cat("Using umap.wnn.harmony reduction for plotting\n")
} else if ("umap" %in% names(myObject@reductions)) {
  reduction_to_use <- "umap"
  cat("Using umap reduction for plotting\n")
} else {
  # Use the first available UMAP reduction
  umap_reductions <- grep("umap", names(myObject@reductions), value = TRUE, ignore.case = TRUE)
  if (length(umap_reductions) > 0) {
    reduction_to_use <- umap_reductions[1]
    cat("Using", reduction_to_use, "reduction for plotting\n")
  } else {
    stop("No UMAP reduction found in the object. Available reductions: ",
         paste(names(myObject@reductions), collapse = ", "))
  }
}

# Generate UMAP plot with new annotations
umap_plot <- DimPlot(myObject,
                     reduction = reduction_to_use,
                     label = TRUE,
                     label.size = 3, pt.size = 1.5,
                     repel = TRUE) +
  ggtitle(paste(mysample, "- Annotated Cell Types")) +
  theme(plot.title = element_text(hjust = 0.5))

# Save UMAP plot to current directory
umap_file <- paste0(mysample, "_annotated_umap.png")
png(umap_file, width = 10, height = 8, units = "in", res = 300)
print(umap_plot)
dev.off()
cat("Saved UMAP plot to:", umap_file, "\n")

# Generate cell type ratio plot by sample
# Assuming the sample information is stored in metadata, adjust the column name as needed
sample_column <- NULL
possible_sample_cols <- c("sample", "orig.ident", "Sample", "sample_id")

for (col in possible_sample_cols) {
  if (col %in% colnames(myObject@meta.data)) {
    sample_column <- col
    break
  }
}

if (is.null(sample_column)) {
  sample_column <- colnames(myObject@meta.data)[1]
  cat("No standard sample column found. Using", sample_column, "as sample identifier\n")
} else {
  cat("Using", sample_column, "as sample identifier\n")
}

# --- FIX: Set CellType as factor with desired order ---
cell_type_ratios <- myObject@meta.data %>%
  mutate(CellType = factor(as.character(Idents(myObject)),
	       levels = c('Cones', 'Rod', 'AC', 'BC', 'MGPC', 'MG')),	
         Sample = .data[[sample_column]]) %>%
  group_by(Sample, CellType) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# Create the ratio plot
ratio_plot <- ggplot(cell_type_ratios, aes(x = Sample, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = paste(mysample, "- Cell Type Ratio by Sample"),
       x = "Sample",
       y = "Percentage of Cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# Save ratio plot to current directory
ratioplot_file <- paste0(mysample, "_celltype_ratio_by_sample.png")
png(ratioplot_file, width = 12, height = 8, units = "in", res = 300)
print(ratio_plot)
dev.off()
cat("Saved cell type ratio plot to:", ratioplot_file, "\n")

# Generate a table summary
summary_table <- cell_type_ratios %>%
  select(Sample, CellType, Count, Percentage) %>%
  arrange(Sample, desc(Percentage))

# Save summary table to current directory
table_file <- paste0(mysample, "_annotation_summary.csv")
write.csv(summary_table, table_file, row.names = FALSE)
cat("Saved annotation summary table to:", table_file, "\n")

# Print summary to console
cat("\n=== ANNOTATION SUMMARY ===\n")
cat("Total cells:", ncol(myObject), "\n")
cat("Number of cell types:", length(unique(Idents(myObject))), "\n")
cat("Samples:", paste(unique(myObject@meta.data[[sample_column]]), collapse = ", "), "\n\n")
print(summary_table)

# Save the renamed object
output_RDS <- paste(mysample, "_annotated.rds", sep="")
saveRDS(myObject, file = output_RDS)
cat("\nSaved renamed object to:", output_RDS, "\n")

# Generate a combined report to current directory
report_file <- paste0(mysample, "_annotation_report.txt")
sink(report_file)
cat("Annotation Report for", mysample, "\n")
cat("Generated on:", date(), "\n\n")
cat("Input files:\n")
cat("- Seurat object:", myRDS, "\n")
cat("- Annotation file:", annotation_file, "\n\n")
cat("Available reductions:", paste(names(myObject@reductions), collapse = ", "), "\n")
cat("Reduction used for plotting:", reduction_to_use, "\n")
cat("Sample column used:", sample_column, "\n\n")
cat("Cell type distribution by sample:\n")
print(summary_table)
cat("\nOutput files generated:\n")
cat("- Annotated Seurat object:", output_RDS, "\n")
cat("- UMAP plot:", umap_file, "\n")
cat("- Cell type ratio plot:", ratioplot_file, "\n")
cat("- Summary table:", table_file, "\n")
sink()

cat("Annotation report saved to:", report_file, "\n")

