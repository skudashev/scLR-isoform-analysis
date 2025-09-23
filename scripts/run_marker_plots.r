# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

# Define input and output
seurat_file <- "sce_seurat.rds"
out_markers_file <- "cluster_markers_transcript.csv"
out_umap_file <- "umap_plot_clusters.pdf"
out_heatmap_file <- "heatmap_top10_markers.pdf"

okabeito_extended <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F5C710", "#0072B2", 
  "#D55E00", "#CC79A7", "#999999", "#000000", "#A6761D",  "#A52A2A",
  "#1B9E77", "#7570B3","#E7298A","#86b370", "#d95f02"
)

# Load Seurat object
seurat_tx <- readRDS(seurat_file)

# # Check cluster assignments
# print("Cluster levels:")
# print(levels(Idents(seurat_tx)))

# # --- Find All Markers ---
# markers <- FindAllMarkers(
#   seurat_tx,
#   assay = "SCT",
#   only.pos = TRUE,
#   min.pct = 0.1,
#   logfc.threshold = 0.25
# )

# # Save full marker list
# write.csv(markers, out_markers_file, row.names = FALSE)

# # --- UMAP Plot ---
# pdf(out_umap_file, width = 6, height = 5)
# DimPlot(seurat_tx, reduction = "umap", label = TRUE, label.size = 4) +
#   ggtitle("UMAP of Transcript Clusters")
# dev.off()

# # --- Heatmap of Top 10 Markers per Cluster ---
# top10 <- markers %>%
#   group_by(cluster) %>%
#   slice_max(order_by = avg_log2FC, n = 10)

# pdf(out_heatmap_file, width = 10, height = 10)
# DoHeatmap(seurat_tx, features = unique(top10$gene), size = 3) +
#   ggtitle("Top 10 Marker Transcripts per Cluster")
# dev.off()

# ran celltypist; got predicted_labels.csv 


add_celltypist_labels_and_plot <- function(seurat_obj, label_file, colour_vector, reduction = "umap") {
  # Extract prefix from filename
  prefix <- sub("predicted_labels\\.csv$", "", basename(label_file))
  
  # Read in predictions
  preds <- fread(label_file)
  setnames(preds, old = names(preds)[1], new = "cell_id")
  
  # Filter for cells in Seurat object
  preds <- preds[cell_id %in% colnames(seurat_obj)]
  
  # Add labels to metadata
  label_vector <- setNames(preds$majority_voting, preds$cell_id)
  seurat_obj[[paste0(prefix, "_labels")]] <- label_vector[colnames(seurat_obj)]
  
  # Create output UMAP file
  umap_file <- paste0("umap_", prefix, "_labels.pdf")
  
  # Plot UMAP
  pdf(umap_file, width = 7, height = 6)
  DimPlot(seurat_obj, reduction = reduction, group.by = paste0(prefix, "_labels"),
            label = TRUE, label.size = 3, repel = TRUE, alpha = 0.7, cols = colour_vector) +
  ggtitle(paste("UMAP with", prefix, "Cell Labels"))
  dev.off()
  
  return(seurat_obj)
}

# Define colour palette
okabeito_extended <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F5C710", "#0072B2", 
  "#D55E00", "#CC79A7", "#999999", "#000000", "#A6761D", "#A52A2A",
  "#1B9E77", "#7570B3", "#E7298A", "#86b370", "#d95f02"
)

# Load your Seurat object
seurat_tx <- readRDS("sce_seurat.rds")

# Get all *_predicted_labels.csv files in current directory
label_files <- list.files(pattern = "_predicted_labels\\.csv$")

# Apply the function to all label files
for (label_file in label_files) {
  seurat_tx <- add_celltypist_labels_and_plot(seurat_tx, label_file, okabeito_extended)
}

# Save updated Seurat object
saveRDS(seurat_tx, "sce_seurat_celltypist.rds")

# FindMarkers Microglia vs Cajal_Retzius 

microglia_vs_cajal <- FindMarkers(
  seurat_tx,
  ident.1 = "Microglia",
  ident.2 = "Cajal_Retzius",
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)