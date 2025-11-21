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

if (!"gene_name" %in% colnames(ann)) stop("peak_annotations must contain gene_name.")
if (!"avg_log2FC" %in% colnames(diff)) stop("diff_peaks must contain avg_log2FC")
if (!"cluster" %in% colnames(diff)) stop("diff_peaks must contain cluster")

cell_order <- c("MG","MGPC","BC","AC","Rod","Cones")

diff <- diff %>% filter(cluster %in% cell_order)

merged <- cbind(diff, ann["gene_name"])

top_peaks <- merged %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice(1:args$n)

gene_scores <- top_peaks %>%
  group_by(gene_name, cluster) %>%
  summarize(score = mean(avg_log2FC, na.rm=TRUE), .groups="drop")

if (nrow(gene_scores) == 0) stop("No gene scores available.")

gene_scores <- gene_scores %>%
  complete(gene_name, cluster = cell_order, fill = list(score = 0))

gene_scores$cluster <- factor(gene_scores$cluster, levels = cell_order)

gene_scores <- gene_scores %>%
  group_by(gene_name) %>%
  mutate(maxpos = which.max(score[match(cell_order, cluster)])) %>%
  ungroup()

gene_scores$gene_name <- factor(
  gene_scores$gene_name,
  levels = unique(gene_scores$gene_name[order(gene_scores$maxpos)])
)

# --------------------------
# UPDATED COLOR PALETTE
# --------------------------
custom_palette <- c(
  "#B5D1E1",  # light blue (low)
  "#C0DAEA",  # very light blue
  "#FFFFFF",  # white (zero)
  "#FDFEFE",  # near-white
  "#E5A07E",  # light salmon
  "#C94832",  # medium red
  "#B5332A"   # deep red (high)
)

# Define normalized positions for colors
value_breaks <- scales::rescale(c(-2, -1, -0.1, 0, 0.3, 1, 2))

# --------------------------
# FINAL PLOT WITH NEW COLORS
# --------------------------
png(args$o, width=1600, height=2000, res=150)
ggplot(gene_scores, aes(x=cluster, y=gene_name, fill=score)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = custom_palette,
    values = value_breaks,
    name = "Log2FC"
  ) +
  theme_minimal() +
  labs(y="Nearby Gene", x="Cell Type") +
  theme(
    axis.text.x = element_text(angle=90, hjust=1),
    axis.text.y = element_text(size=10),
    panel.grid = element_blank()
  )
dev.off()

