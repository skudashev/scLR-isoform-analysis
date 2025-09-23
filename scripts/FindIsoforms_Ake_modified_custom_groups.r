library(Seurat)
library(data.table)
library(dplyr)
library(pbapply)

# -------------------------------
# Inputs you already have
# -------------------------------
seurat_tx <- readRDS("/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/seurat_out/sce_seurat.rds")
t2g <- fread(
  "/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/isosceles_isoquant/sqanti3_filter/transcript2gene.tsv",
  header = FALSE, col.names = c("transcript_id", "gene_name")
)

threshold_abund <- 0.10     # isoform must reach ≥10% in at least one cluster
threshold_var   <- 0.05     # SD of proportions across clusters must be ≥ 0.05
min_cells_per_subcluster <- 100

# Your custom grouping of *Idents(seurat_tx)*
my_groups <- list(
  microglia  = c(9, 16, 10, 13, 8, 3, 6),
  neurons    = c(11, 2, 5, 4, 15),
  astrocytes = c(14, 12, 1, 7)
)

# -------------------------------
# Helpers
# -------------------------------

# Aggregate counts per subcluster (within the current object)
aggregate_by_subcluster <- function(seurat_obj, min_cells = 100) {
  cl_counts <- table(Idents(seurat_obj))
  valid     <- names(cl_counts[cl_counts >= min_cells])
  seurat_obj <- subset(seurat_obj, idents = valid)
  Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = valid)

  agg <- AggregateExpression(seurat_obj, group.by = "ident")[[1]]
  as.data.table(agg, keep.rownames = "transcript_id")
}

# xTest (kept close to your original)
Perform_xTest <- function(tab, threshold_abund = 0.05, threshold_var = 0.05) {
  num_cols <- names(tab)[sapply(tab, is.numeric)]

  # per-isoform totals across groups (not strictly needed below but retained)
  tab <- tab |>
    rowwise() |>
    mutate(row_sum = sum(c_across(all_of(num_cols)))) |>
    ungroup()

  # per-group percentages for plotting/filtering (global per gene set)
  for (col in num_cols) {
    pct_col <- paste0(col, "_pct")
    total_group_expr <- sum(tab[[col]], na.rm = TRUE)
    tab[[pct_col]] <- if (total_group_expr > 0) tab[[col]] / total_group_expr else 0
  }

  # Abundance/variance filters across groups
  # Keep rows (isoforms) that hit the abundance threshold in any group and have variability
  pct_cols <- grep("_pct$", names(tab), value = TRUE)
  tab <- tab |>
    mutate(
      thr.abund = apply(across(all_of(pct_cols)), 1, function(x) max(x, na.rm = TRUE) < threshold_abund),
      thr.var   = apply(across(all_of(num_cols)), 1, function(x) sd(x, na.rm = TRUE) < threshold_var)
    ) |>
    select(-row_sum)

  numTab <- tab |>
    filter(!thr.abund & !thr.var) |>
    select(all_of(num_cols))

  if (nrow(numTab) > 1) {
    chi2 <- suppressWarnings(chisq.test(as.matrix(numTab)))
    if (all(round(chi2$expected, 0) >= 5)) {
      tab <- mutate(tab, p.value = chi2$p.value, statistic = chi2$statistic)
    }
  }
  tab
}

# Adjusted p-values and final filtering
process_results <- function(df, threshold_pval = 0.05) {
  df <- df %>%
    mutate(
      p.value.adjusted = ifelse(
        thr.abund == FALSE & thr.var == FALSE & !is.na(p.value),
        0, NA
      )
    )

  tmpA <- filter(df, p.value.adjusted == 0)
  if (nrow(tmpA) > 0) tmpA$p.value.adjusted <- p.adjust(tmpA$p.value, method = "BH")
  tmpB <- filter(df, is.na(p.value.adjusted))

  df <- bind_rows(tmpA, tmpB) %>%
    arrange(p.value.adjusted, gene_name, transcript_name) %>%
    distinct() %>%
    data.table::data.table()

  df <- df %>%
    filter(p.value.adjusted < threshold_pval) %>%
    select(!c("thr.abund", "thr.var", "transcript_id"))

  df
}

# Subset seurat_tx to a *named group* of cluster IDs and run the DTU test
run_group_switching <- function(seurat_tx, group_name, cluster_ids,
                                min_cells = 100,
                                threshold_abund = 0.10, threshold_var = 0.05) {
  # Ensure identities are characters for matching
  idents_chr <- as.character(Idents(seurat_tx))
  keep_cells <- colnames(seurat_tx)[idents_chr %in% as.character(cluster_ids)]
  sub <- subset(seurat_tx, cells = keep_cells)
  Idents(sub) <- factor(Idents(sub))  # tidy Idents

  # Aggregate per subcluster
  agg_dt <- aggregate_by_subcluster(sub, min_cells = min_cells)

  # Merge gene names and keep numeric columns (subcluster columns)
  agg_dt <- merge(t2g, agg_dt, by = "transcript_id")
  agg_dt[, transcript_name := transcript_id]

  # Keep only numeric columns (clusters) plus IDs
  num_keep <- names(agg_dt)[sapply(agg_dt, is.numeric)]
  agg_dt <- agg_dt[, c("gene_name", "transcript_name", "transcript_id", num_keep), with = FALSE]

  # Multi-isoform genes only
  multi_iso <- agg_dt[, .N, by = gene_name][N > 1, gene_name]
  agg_dt <- agg_dt[gene_name %in% multi_iso]

  # Run xTest per gene
  res <- pblapply(
    split(agg_dt, agg_dt$gene_name),
    Perform_xTest,
    threshold_abund = threshold_abund,
    threshold_var   = threshold_var
  ) |> rbindlist(fill = TRUE)

  res[, major_group := group_name]
  res
}

# -------------------------------
# Run for your custom groups
# -------------------------------
output_dir <- "isoform_switching_results_custom_groups"
dir.create(output_dir, showWarnings = FALSE)

all_group_results <- list()
for (grp in names(my_groups)) {
  message("Processing group: ", grp, " (clusters: ", paste(my_groups[[grp]], collapse = ", "), ")")
  all_group_results[[grp]] <- run_group_switching(
    seurat_tx,
    group_name = grp,
    cluster_ids = my_groups[[grp]],
    min_cells = min_cells_per_subcluster,
    threshold_abund = threshold_abund,
    threshold_var   = threshold_var
  )
}

# Process, save, and report counts
processed_results <- lapply(all_group_results, process_results)
for (grp in names(processed_results)) {
  res <- processed_results[[grp]]
  saveRDS(res, file = file.path(output_dir, paste0(grp, "_isoform_switches.rds")))
  fwrite(res, file = file.path(output_dir, paste0(grp, "_isoform_switches.tsv")),
         sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("Group %-12s: %4d genes with DTU",
                  grp, length(unique(res$gene_name))))
}

# -------------------------------
# Re-run microglia excluding cluster 9
# -------------------------------
if ("microglia" %in% names(my_groups)) {
  mg_no9_ids <- setdiff(my_groups[["microglia"]], 9)
  res_mg_no9 <- run_group_switching(
    seurat_tx,
    group_name = "microglia_no9",
    cluster_ids = mg_no9_ids,
    min_cells = min_cells_per_subcluster,
    threshold_abund = threshold_abund,
    threshold_var   = threshold_var
  )
  res_mg_no9_processed <- process_results(res_mg_no9)

  fwrite(res_mg_no9_processed,
         file = file.path(output_dir, "microglia_no9_isoform_switches.tsv"),
         sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(res_mg_no9_processed,
          file = file.path(output_dir, "microglia_no9_isoform_switches.rds"))
  message(sprintf("Group %-12s: %4d genes with DTU",
                  "microglia_no9", length(unique(res_mg_no9_processed$gene_name))))
}
