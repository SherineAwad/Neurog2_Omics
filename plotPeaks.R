library(argparse)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ----------------------------

# Arguments

# ----------------------------

parser <- ArgumentParser()
parser$add_argument("--rds", required=TRUE)
parser$add_argument("-n", type="integer", default=50)
parser$add_argument("-o", required=TRUE)
parser$add_argument("-g", required=TRUE, help="Text file with a list of gene names to include")
args <- parser$parse_args()

# ----------------------------

# Load Seurat object and data

# ----------------------------

obj <- readRDS(args$rds)
stopifnot(!is.null(obj@misc$diff_peaks), !is.null(obj@misc$peak_annotations))

diff <- obj@misc$diff_peaks
ann  <- obj@misc$peak_annotations
diff <- diff %>% rename(peak_id = gene)

stopifnot(all(c("gene_name","query_region") %in% colnames(ann)))
stopifnot(all(c("avg_log2FC","cluster") %in% colnames(diff)))

cell_order <- c("MG","MGPC","BC","AC","Rod","Cones")
diff <- diff %>% filter(cluster %in% cell_order)

# ----------------------------

# Merge annotations (allow many-to-many)

# ----------------------------

merged <- diff %>%
left_join(
ann %>% select(query_region, gene_name),
by = c("peak_id" = "query_region"),
relationship = "many-to-many"
)

# ----------------------------

# Build matrix: rows = peaks, columns = clusters

# ----------------------------

mat_all <- merged %>%
select(peak_id, cluster, avg_log2FC, gene_name) %>%
distinct() %>%
pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)

mat_all <- mat_all[, c("peak_id","gene_name", cell_order)]
peak_ids <- mat_all$peak_id
gene_labels <- mat_all$gene_name

mat_dense <- as.matrix(mat_all[, cell_order])
rownames(mat_dense) <- peak_ids

# ----------------------------

# Compute z-scores per peak

# ----------------------------

mat_z <- t(scale(t(mat_dense), center = TRUE, scale = TRUE))
mat_z[is.na(mat_z)] <- 0
rownames(mat_z) <- peak_ids
colnames(mat_z) <- cell_order

# ----------------------------

# Select top N peaks per cluster (both up and down)

# ----------------------------

top_peak_ids <- unique(unlist(lapply(cell_order, function(cl) {
cl_scores <- mat_z[, cl]
top_idx <- order(cl_scores, decreasing = TRUE)[1:min(args$n, length(cl_scores))]
bottom_idx <- order(cl_scores, decreasing = FALSE)[1:min(args$n, length(cl_scores))]
rownames(mat_z)[c(top_idx, bottom_idx)]
})))

mat_top <- mat_z[top_peak_ids, , drop=FALSE]
gene_labels_top <- gene_labels[match(rownames(mat_top), peak_ids)]

# ----------------------------

# Load gene list and label

# ----------------------------

gene_list <- read.table(args$g, stringsAsFactors = FALSE)[,1]
y_labels <- ifelse(gene_labels_top %in% gene_list, gene_labels_top, "")

# ----------------------------

# Convert for ggplot

# ----------------------------

gene_scores <- mat_top %>%
as.data.frame() %>%
mutate(row_id = 1:nrow(mat_top)) %>%
pivot_longer(cols = all_of(cell_order), names_to="cluster", values_to="score")

gene_scores$cluster <- factor(gene_scores$cluster, levels = cell_order)

gene_scores <- gene_scores %>%
group_by(row_id) %>%
mutate(maxpos = which.max(score[match(cell_order, cluster)])) %>%
ungroup()

gene_scores$row_id <- factor(
gene_scores$row_id,
levels = unique(gene_scores$row_id[order(gene_scores$maxpos)])
)

# ----------------------------

# Colors

# ----------------------------

custom_palette <- c(
"#2166AC",
"#67A9CF",
"#D1E5F0",
"white",
"#FDDBC7",
"#EF8A62",
"#B2182B"
)

# ----------------------------

# Plot heatmap

# ----------------------------

png(args$o, width=1600, height=2000, res=150)
ggplot(gene_scores, aes(x=cluster, y=row_id, fill=score)) +
geom_tile() +
scale_fill_gradientn(
colors = custom_palette,
limits = c(min(gene_scores$score), max(gene_scores$score)),
name = "Z-score"
) +
scale_y_discrete(labels = y_labels) +
theme_minimal() +
labs(y="Nearby Gene", x="Cell Type") +
theme(
axis.text.x = element_text(angle=90, hjust=1),
axis.text.y = element_text(size=10),
panel.grid = element_blank()
)
dev.off()

