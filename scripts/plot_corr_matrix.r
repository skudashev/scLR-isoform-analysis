library(data.table)
library(Seurat)
library(pheatmap)

seurat_tx_file <- "sce_seurat_celltypist.rds"
labels_file <- "Developing_Human_Hippocampus_updatedpredicted_labels.csv"
seurat_tx <- readRDS(seurat_tx_file)
labels <- fread(labels_file)
setnames(labels, old = names(labels)[1], new = "cell_id")  # rename first col

# Filter to matching cells
labels <- labels[cell_id %in% colnames(seurat_tx)]
label_vec <- setNames(labels$new_labels, labels$cell_id)
seurat_tx$Developing_Human_Hippocampus_updated_labels <- label_vec[colnames(seurat_tx)]

aggregate_by_subcluster <- function(seurat_obj, min_cells = 100, assay = "originalexp", layer = "counts") {
  # Filter out small subclusters
  cluster_counts <- table(Idents(seurat_obj))
  valid_clusters <- names(cluster_counts[cluster_counts >= min_cells])
  seurat_obj <- subset(seurat_obj, idents = valid_clusters)
  # Re-set identity to valid clusters
  Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = valid_clusters)
  # Aggregate expression by subcluster
  agg_expr <- AggregateExpression(
    seurat_obj,
    group.by = "ident",
    assays = assay,
    layers = layer,
    verbose = FALSE
  )[[1]]

  as.data.table(agg_expr, keep.rownames = "transcript_id")
}


major_groups_to_test <- c("Cajal_Retzius", "Microglia", "Endothelial")
all_heatmaps <- list()

# Define shared colour palette from blue (low) to red (high)
# col_fun <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))

col_fun <- colorRampPalette(rev(c(
  "#004C4C",  # dark teal
  "#006f51",  # teal
  "#009E73",  # Okabe–Ito green
  "#4dbb9d",  # very light teal
  "#b3e2d5"   # near-white teal
)))


breaks <- seq(0.1, 1, length.out = 100)

for (major_group in major_groups_to_test) {
  message("Processing: ", major_group)

  # Subset by major label
  subset_obj <- subset(seurat_tx, subset = Developing_Human_Hippocampus_updated_labels == major_group)
  Idents(subset_obj) <- subset_obj$seurat_clusters

  # Aggregate
  agg_dt <- aggregate_by_subcluster(subset_obj)

  # Format to matrix
  expr_mat <- as.matrix(agg_dt[, -1, with = FALSE])
  rownames(expr_mat) <- agg_dt$transcript_id

  # Correlation matrix (values between 0 and 1)
  cor_mat <- cor(expr_mat, method = "pearson")
  cor_mat[cor_mat < 0] <- 0  

  # Output
  pdf(paste0("subcluster_correlation_", major_group, ".pdf"), width = 7, height = 6)
  pheatmap(
    cor_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    fontsize = 14,
    fontsize_number = 12,
    breaks = breaks,
    number_color = "white",
    angle_col = 0,
    color = col_fun(length(breaks)),
    main = paste(major_group)
  )
  dev.off()
}










