library(argparse)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(tidyr)

parser <- ArgumentParser()
parser$add_argument("--rds", required=TRUE)
parser$add_argument("-n", type="integer", default=50)
parser$add_argument("-o", required=TRUE)
args <- parser$parse_args()

obj <- readRDS(args$rds)

if (is.null(obj@misc$diff_peaks)) stop("diff_peaks not found.")
if (is.null(obj@misc$peak_annotations)) stop("peak_annotations not found.")

diff <- obj@misc$diff_peaks
ann  <- obj@misc$peak_annotations

# diff$gene is actually the peak ID → rename it
diff <- diff %>% rename(peak_id = gene)

if (!"gene_name" %in% colnames(ann)) stop("peak_annotations must contain gene_name.")
if (!"query_region" %in% colnames(ann)) stop("peak_annotations must contain query_region.")
if (!"avg_log2FC" %in% colnames(diff)) stop("diff_peaks must contain avg_log2FC")
if (!"cluster" %in% colnames(diff)) stop("diff_peaks must contain cluster")

cell_order <- c("MG","MGPC","BC","AC","Rod","Cones")
diff <- diff %>% filter(cluster %in% cell_order)

# Proper merge: peak_id ↔ query_region
merged <- diff %>%
  left_join(ann %>% select(query_region, gene_name),
            by = c("peak_id" = "query_region"))

# --------------------------------------------------------------------
# Build matrix: rows = peaks, columns = clusters
# --------------------------------------------------------------------
mat_all <- merged %>%
  select(peak_id, cluster, avg_log2FC, gene_name) %>%
  distinct() %>%
  pivot_wider(
    names_from = cluster,
    values_from = avg_log2FC,
    values_fill = 0
  )

mat_all <- mat_all[, c("peak_id","gene_name", cell_order)]
peak_ids <- mat_all$peak_id
gene_labels <- mat_all$gene_name

mat_dense <- as.matrix(mat_all[, cell_order])
rownames(mat_dense) <- peak_ids

# --------------------------------------------------------------------
# Compute z-scores
# --------------------------------------------------------------------
mat_z <- t(apply(mat_dense, 1, function(x) {
  x <- as.numeric(x)
  if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)
}))
rownames(mat_z) <- peak_ids
colnames(mat_z) <- cell_order

# --------------------------------------------------------------------
# Rank peaks by max |z| and keep top N
# --------------------------------------------------------------------
peak_rank <- order(apply(mat_z, 1, function(x) max(abs(x))), decreasing=TRUE)
top_peak_ids <- peak_ids[peak_rank][1:args$n]

mat_top <- mat_z[top_peak_ids, , drop=FALSE]

rownames(mat_top) <- gene_labels[match(rownames(mat_top), peak_ids)]

# --------------------------------------------------------------------
# Convert for ggplot
# --------------------------------------------------------------------
gene_scores <- mat_top %>%
  as.data.frame() %>%
  mutate(gene_name = rownames(.)) %>%
  pivot_longer(cols = all_of(cell_order), names_to="cluster", values_to="score")

gene_scores$cluster <- factor(gene_scores$cluster, levels = cell_order)

gene_scores <- gene_scores %>%
  group_by(gene_name) %>%
  mutate(maxpos = which.max(score[match(cell_order, cluster)])) %>%
  ungroup()

gene_scores$gene_name <- factor(
  gene_scores$gene_name,
  levels = unique(gene_scores$gene_name[order(gene_scores$maxpos)])
)

# --------------------------------------------------------------------
# Colors (unchanged)
# --------------------------------------------------------------------
custom_palette <- c(
  "#B5D1E1",
  "#C0DAEA",
  "#FFFFFF",
  "#FDFEFE",
  "#E5A07E",
  "#C94832",
  "#B5332A"
)

value_breaks <- scales::rescale(c(-2, -1, -0.1, 0, 0.3, 1, 2))

# --------------------------------------------------------------------
# Plot heatmap
# --------------------------------------------------------------------
png(args$o, width=1600, height=2000, res=150)
ggplot(gene_scores, aes(x=cluster, y=gene_name, fill=score)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = custom_palette,
    values = value_breaks,
    name = "Z-score"
  ) +
  theme_minimal() +
  labs(y="Nearby Gene", x="Cell Type") +
  theme(
    axis.text.x = element_text(angle=90, hjust=1),
    axis.text.y = element_text(size=10),
    panel.grid = element_blank()
  )
dev.off()

