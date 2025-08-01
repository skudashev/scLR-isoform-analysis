#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(dplyr)
})

# Arguments
args <- commandArgs(trailingOnly = TRUE)
seurat_rds <- args[1]
t2g_file   <- args[2]
output_csv <- args[3]

# Load Seurat object
seurat_gene <- readRDS(seurat_rds)

# Extract counts
gene_counts <- GetAssayData(seurat_gene, layer = "counts")
gene_counts_t <- as.data.frame(Matrix::t(gene_counts))
gene_counts_dt <- as.data.table(gene_counts_t, keep.rownames = "cell_id")

# Keep only ENSG columns + cell_id
keep_cols <- grepl("^ENSG", names(gene_counts_dt)) | names(gene_counts_dt) == "cell_id"
gene_counts_dt <- gene_counts_dt[, ..keep_cols]

# Map ENSG to gene_name
t2g <- fread(t2g_file, header = FALSE, col.names = c("transcript_id", "gene_id", "gene_name"))
id_to_name <- setNames(t2g$gene_name, t2g$gene_id)
gene_cols <- setdiff(names(gene_counts_dt), "cell_id")
renamed_cols <- ifelse(gene_cols %in% names(id_to_name), id_to_name[gene_cols], gene_cols)
setnames(gene_counts_dt, old = gene_cols, new = renamed_cols)

# Keep renamed columns
keep_cols <- !grepl("^ENSG", names(gene_counts_dt)) | names(gene_counts_dt) == "cell_id"
gene_counts_dt <- gene_counts_dt[, ..keep_cols]

# Save
fwrite(gene_counts_dt, output_csv)
