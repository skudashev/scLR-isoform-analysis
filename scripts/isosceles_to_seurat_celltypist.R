#!/usr/bin/env Rscript

suppressPackageStartupMessages({

  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scDblFinder)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tools)
  library(future)
})

args <- commandArgs(trailingOnly = TRUE)

# Simple named parsing: --key value
parse_args <- function(args) {
  out <- list()
  for (i in seq(1, length(args), 2)) {
    key <- gsub("^--", "", args[i])
    val <- args[i+1]
    out[[key]] <- val
  }
  out
}

opts <- parse_args(args)

# set defaults if missing
opts$iso_dir <- opts$iso_dir %||% "."
opts$transcript_rds <- opts$transcript_rds %||% "se_transcript.rds"
opts$gene_rds <- opts$gene_rds %||% "se_gene.rds"
opts$t2g <- opts$t2g %||% "ref_transcript2gene.tsv"
opts$out_dir <- opts$out_dir %||% "."
opts$npcs <- as.integer(opts$npcs %||% 50)


# ---------- Helpers ----------
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse=" ")))

ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

parse_dims <- function(x) {
  # Accept "1:50" or "1,2,3,4"
  if (grepl(":", x)) {
    rng <- strsplit(x, ":", fixed=TRUE)[[1]]
    return(seq(as.integer(rng[1]), as.integer(rng[2])))
  } else {
    return(as.integer(strsplit(x, ",")[[1]]))
  }
}

read_isosceles <- function(dir, fname) {
  f <- file.path(dir, fname)
  if (!file.exists(f)) stop("File not found: ", f)
  obj <- readRDS(f)
  if (!inherits(obj, "SummarizedExperiment")) stop("Object is not a SummarizedExperiment: ", f)
  obj
}

sce_from_se <- function(se, assay_counts) {
  if (!assay_counts %in% assayNames(se)) {
    stop("Assay '", assay_counts, "' not found. Available: ", paste(assayNames(se), collapse=", "))
  }
  as(se, "SingleCellExperiment")
}

run_dbl_finder <- function(sce, seed=42) {
  set.seed(seed)
  sce <- scDblFinder(sce)
  table(sce$scDblFinder.class)
  sce[, sce$scDblFinder.class == "singlet"]
}

seurat_from_sce <- function(sce, assay_counts="counts", min_cells=3, min_features=0) {
  # Create Seurat object from SCE counts
  if (!assay_counts %in% assayNames(sce)) stop("Assay '", assay_counts, "' missing in SCE.")
  counts <- assay(sce, assay_counts)
  so <- CreateSeuratObject(counts = counts, min.cells = min_cells, min.features = min_features)
  # Carry over metadata
  meta <- as.data.frame(colData(sce))
  common <- intersect(colnames(meta), colnames(so[[]]))
  meta_to_add <- meta[, setdiff(colnames(meta), common), drop=FALSE]
  if (ncol(meta_to_add) > 0) so <- AddMetaData(so, metadata = meta_to_add)
  so
}

process_seurat <- function(so, npcs=50, dims=1:50, resolution=0.7, cluster_alg=4, sct_ncells=10000, seed=42) {
  set.seed(seed)
  DefaultAssay(so) <- "RNA"
  plan("multicore", workers = 8)
  so <- SCTransform(so, verbose = TRUE, assay = "RNA", ncells = sct_ncells)
  so <- RunPCA(so, npcs = npcs, assay = "SCT")
  so <- RunUMAP(so, dims = dims, assay = "SCT")
  so <- RunUMAP(so, dims = 1:50, umap.method = "uwot", n.threads = 8)
  so <- FindNeighbors(so, dims = dims, assay = "SCT", n.threads = 8)
  so <- FindClusters(so, algorithm = cluster_alg, resolution = resolution)
  so
}

transfer_clusters_to_gene_se <- function(se_gene, seurat_tx) {
  tx_clusters <- Idents(seurat_tx)
  # Ensure overlapping barcodes
  common_cells <- intersect(colnames(se_gene), names(tx_clusters))
  if (length(common_cells) == 0) stop("No overlapping cell barcodes between gene and transcript objects.")
  colData(se_gene)$tx_cluster <- NA_character_
  colData(se_gene)[common_cells, "tx_cluster"] <- as.character(tx_clusters[common_cells])
  se_gene
}

export_celltypist <- function(seurat_gene, t2g_file, out_dir, prefix="celltypist") {
  # Read mapping: transcript_id, gene_id, gene_name
  if (!file.exists(t2g_file)) stop("Mapping file not found: ", t2g_file)
  t2g <- fread(t2g_file, header = FALSE, col.names = c("transcript_id", "gene_id", "gene_name"))
  if (!all(c("gene_id", "gene_name") %in% names(t2g)))
    stop("Mapping must contain at least gene_id and gene_name columns.")
  id_to_name <- setNames(t2g$gene_name, t2g$gene_id)

  # Use raw counts for CellTypist
  mat <- GetAssayData(seurat_gene, assay = "RNA", slot = "counts")
  # Convert ENSG IDs to symbols where available
  gene_ids <- rownames(mat)
  gene_names <- ifelse(gene_ids %in% names(id_to_name), id_to_name[gene_ids], gene_ids)

  # Collapse duplicates after renaming (sum counts for genes mapping to same symbol)
  # Build a sparse matrix with new row names
  mat_df <- as(mat, "dgTMatrix")
  mat_df@i <- match(gene_names[mat_df@i + 1], unique(gene_names)) - 1
  collapsed <- sparseMatrix(
    i = mat_df@i, j = mat_df@j, x = mat_df@x,
    dims = c(length(unique(gene_names)), ncol(mat)), dimnames = list(unique(gene_names), colnames(mat))
  )

  # Export CSV (cells as rows, genes as columns)
  msg("Writing CellTypist CSV...")
  dt <- as.data.table(Matrix::t(collapsed), keep.rownames = "cell_id")
  fwrite(dt, file.path(out_dir, paste0(prefix, "_counts_cells_as_rows.csv")))

  Matrix::writeMM(collapsed, file.path(out_dir, paste0(prefix, "_counts.mtx")))
  fwrite(data.table(gene=rownames(collapsed)), file.path(out_dir, paste0(prefix, "_genes.tsv")), col.names = FALSE)
  fwrite(data.table(cell=colnames(collapsed)), file.path(out_dir, paste0(prefix, "_barcodes.tsv")), col.names = FALSE)

  invisible(NULL)
}

# ---------- Main ----------
main <- function(opt) {
  set.seed(opt$seed)
  ensure_dir(opt$out_dir)

  msg("Loading Isosceles transcript-level object...")
  se_tx <- read_isosceles(opt$iso_dir, opt$transcript_rds)

  msg("Converting to SCE and removing doublets with scDblFinder...")
  sce_tx <- sce_from_se(se_tx, assay_counts = opt$assay_counts)
  sce_tx <- run_dbl_finder(sce_tx, seed = opt$seed)
  saveRDS(sce_tx, file.path(opt$out_dir, "sce_filtered.rds"))
  msg("Saved filtered SCE: ", file.path(opt$out_dir, "sce_filtered.rds"))

  msg("Creating Seurat object (transcript-level) and processing...")
  seu_tx <- seurat_from_sce(sce_tx, assay_counts = opt$assay_counts,
                            min_cells = opt$min_cells, min_features = opt$min_features)
  seu_tx <- process_seurat(seu_tx,
                           npcs = opt$npcs,
                           dims = parse_dims(opt$dims),
                           resolution = opt$resolution,
                           cluster_alg = opt$cluster_alg,
                           sct_ncells = opt$sct_ncells,
                           seed = opt$seed)
  saveRDS(seu_tx, file.path(opt$out_dir, "sce_seurat.rds"))
  msg("Saved transcript-level Seurat: ", file.path(opt$out_dir, "sce_seurat.rds"))

  # Load gene-level SE and attach transcript clusters
  msg("Loading Isosceles gene-level object and transferring clusters...")
  se_gene <- read_isosceles(opt$iso_dir, opt$gene_rds)
  se_gene <- transfer_clusters_to_gene_se(se_gene, seu_tx)
  saveRDS(se_gene, file.path(opt$out_dir, "se_gene_with_tx_clusters.rds"))

  msg("Creating Seurat object (gene-level)...")
  # Use same assay name for counts if present, otherwise coerce
  if (!opt$assay_counts %in% assayNames(se_gene)) {
    stop("Assay '", opt$assay_counts, "' not present in gene-level SE.")
  }
  seu_gene <- CreateSeuratObject(
    counts = assay(se_gene, opt$assay_counts),
    meta.data = as.data.frame(colData(se_gene)),
    min.cells = opt$min_cells, min.features = opt$min_features
  )
  saveRDS(seu_gene, file.path(opt$out_dir, "sce_gene_seurat.rds"))
  msg("Saved gene-level Seurat: ", file.path(opt$out_dir, "sce_gene_seurat.rds"))

  n_clusters <- length(unique(colData(se_gene)$tx_cluster))
  msg(sprintf("Number of transcript-level clusters: %s", n_clusters))

  msg("Exporting counts for CellTypist...")
  export_celltypist(seu_gene, t2g_file = file.path(opt$iso_dir, opt$t2g),
                    out_dir = opt$out_dir, prefix = opt$export_prefix)

  msg("Done.")
  sessionInfo()
}

if (identical(environment(), globalenv())) {
  main(opt)
}
