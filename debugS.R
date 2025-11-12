# -----------------------------
# Debug script to inspect mitochondrial genes
# -----------------------------
library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]
myh5 <- paste0(mysample, "_filtered_feature_bc_matrix.h5")

cat("Sample:", mysample, "\n")
cat("H5 file:", myh5, "\n\n")

# -----------------------------
# Load RNA counts only
# -----------------------------
counts <- Read10X_h5(myh5)
rna_counts <- counts$`Gene Expression`
cat("RNA matrix dimensions (genes x cells):", dim(rna_counts), "\n")

# -----------------------------
# Inspect first few gene names
# -----------------------------
cat("\nFirst 20 gene names:\n")
print(rownames(rna_counts)[1:20])

# -----------------------------
# Inspect genes that contain "mt" (case-insensitive)
# -----------------------------
mt_candidates <- rownames(rna_counts)[grepl("mt", rownames(rna_counts), ignore.case = TRUE)]
cat("\nNumber of genes containing 'mt' in their name:", length(mt_candidates), "\n")
cat("First 20 'mt' genes:\n")
print(head(mt_candidates, 20))

# -----------------------------
# Inspect RNA counts for the first few cells
# -----------------------------
cat("\nRNA counts for first 5 cells and first 20 genes:\n")
print(as.matrix(rna_counts[1:20, 1:5]))

# -----------------------------
# Optionally, inspect top expressed genes
# -----------------------------
total_counts_per_gene <- rowSums(rna_counts)
top_genes <- sort(total_counts_per_gene, decreasing = TRUE)[1:20]
cat("\nTop 20 expressed genes:\n")
print(top_genes)

# -----------------------------
# Calculate and plot mitochondrial percentage histogram
# -----------------------------
mt_genes <- rownames(rna_counts)[grepl("^mt-", rownames(rna_counts))]
cat("\nMitochondrial genes with '^mt-' pattern:", length(mt_genes), "\n")
cat("Mitochondrial genes found:\n")
print(mt_genes)

if(length(mt_genes) > 0) {
  # Calculate mitochondrial percentage per cell
  mt_counts <- colSums(rna_counts[mt_genes, ])
  total_counts <- colSums(rna_counts)
  percent_mt <- (mt_counts / total_counts) * 100
  
  cat("\nMitochondrial percentage summary:\n")
  print(summary(percent_mt))
  cat("Range:", range(percent_mt), "\n")
  
  # Create histogram
  hist_file <- paste0(mysample, "_mt_percent_histogram.png")
  png(hist_file, width = 800, height = 600)
  hist(percent_mt, 
       main = paste("Mitochondrial Percentage -", mysample),
       xlab = "Percent Mitochondrial Reads", 
       ylab = "Number of Cells",
       col = "lightblue",
       breaks = 50)
  abline(v = median(percent_mt), col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("Median:", round(median(percent_mt), 2), "%"), 
         col = "red", lwd = 2, lty = 2)
  dev.off()
  cat("\nHistogram saved to:", hist_file, "\n")
} else {
  cat("\nNo mitochondrial genes found with pattern '^mt-'\n")
}

cat("\n--- Debugging complete ---\n")
