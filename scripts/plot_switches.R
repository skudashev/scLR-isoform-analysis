# ---- Libraries ----
library(dplyr)
library(tidyr)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(ggtranscript)
library(patchwork)
library(purrr)
library(data.table)

# ---- Inputs ----
# Set these paths as needed
gtf_path <- "/Volumes/Projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/isosceles_isoquant/transcript_annotation.sorted.gtf.gz"
# switches_path <- "/Volumes/Projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/isosceles_isoquant/isoform_switching_results/Microglia_isoform_switches.tsv"
switches_path <- "/Volumes/Projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/isosceles_isoquant/isoform_switching_results_new/Microglia_isoform_switches.tsv"
capture_isoforms <- "/Volumes/Projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/isosceles_isoquant/sqanti3_filter/Capture_isoforms_with_gene.txt"

# ---- Load data ----
capture_isoforms <- read.delim(capture_isoforms, header = TRUE, stringsAsFactors = FALSE)
gtf <- rtracklayer::import(gtf_path)

# Build transcript metadata map once from the GTF
transcript_metadata <- as.data.frame(gtf) %>%
  dplyr::select(transcript_id, compatible_tx) %>%
  dplyr::distinct()

# read class file
class_file <- "/Volumes/Projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/isosceles_isoquant/sqanti3_output/fsm_classification.txt"
class_df <- fread(class_file, header = TRUE, sep = "\t",
                  select = c("isoform", "structural_category", "subcategory", "associated_transcript"))
# left_join with transcript_metadata and if structural_category == full-splice_match and subcategory == reference_match replace compatible_tx with associated_transcript
transcript_metadata <- transcript_metadata %>%
  dplyr::left_join(class_df, by = c("transcript_id" = "isoform")) %>%
  dplyr::mutate(
    compatible_tx = dplyr::case_when(
      structural_category == "full-splice_match" & subcategory == "reference_match" ~ associated_transcript,
      subcategory %in% c("alternative_3end", "alternative_5end", "alternative_3end5end") ~ paste0(associated_transcript, "altend"),
      TRUE ~ compatible_tx
    )
  ) %>%
  dplyr::select(transcript_id, compatible_tx) %>%
  dplyr::distinct()

# change gtf - replace old compatible_tx with new one
tx_map <- setNames(transcript_metadata$compatible_tx, transcript_metadata$transcript_id)
mcols(gtf)$compatible_tx <- tx_map[mcols(gtf)$transcript_id]


isoform_switches_raw <- read.delim(switches_path, check.names = FALSE)

# Add isoform_id to switches table, preferring compatible_tx when available
isoform_switches <- isoform_switches_raw %>%
  dplyr::left_join(transcript_metadata, by = c("transcript_name" = "transcript_id")) %>%
  dplyr::mutate(
    isoform_id = ifelse(!is.na(compatible_tx) & compatible_tx != "", compatible_tx, transcript_name)
  ) %>%
  dplyr::select(-compatible_tx)

# ---- Helper: colour palette that scales to N isoforms ----
make_palette <- function(n) {
  base <- c("#E69F00", "#56B4E9", "#009E73", "#F5C710", "#0072B2",
            "#D55E00", "#CC79A7", "#999999", "#000000", "#A6761D",
            "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#D95F02",
            "#8C564B", "#17BECF", "#9467BD")
  if (n <= length(base)) {
    base[seq_len(n)]
  } else {
    grDevices::colorRampPalette(base)(n)
  }
}

# ---- Core plotting function ----
# Returns a list with: combined (patchwork), structure_plot, fraction_plot, and data frames used
plot_gene_switch <- function(gtf, switches_table, gene_of_interest,
                             cluster_suffix = "_pct",
                             shorten_intron_min = 300) {
  # 1) Filter switches for the gene and reshape fractions
  sw_gene <- switches_table %>%
    dplyr::filter(gene_name == gene_of_interest)
  
  if (nrow(sw_gene) == 0) {
    stop(sprintf("No rows found for gene '%s' in the switches table.", gene_of_interest))
  }
  
  # Identify fraction columns that end with the chosen suffix
  frac_cols <- grep(paste0(cluster_suffix, "$"), names(sw_gene), value = TRUE)
  if (length(frac_cols) == 0) {
    stop(sprintf("No columns ending with '%s' found. Check your switches table.", cluster_suffix))
  }
  
  frac_long <- sw_gene %>%
    dplyr::select(isoform_id, transcript_name, all_of(frac_cols)) %>%
    tidyr::pivot_longer(cols = all_of(frac_cols), names_to = "cluster", values_to = "fraction") %>%
    dplyr::mutate(
      cluster = sub(paste0(cluster_suffix, "$"), "", cluster),
      # Keep everything after the last dot if clusters are like 'Microglia.7'
      cluster = sub(".*\\.", "", cluster),
      # Order clusters numerically if possible
      cluster_ord = suppressWarnings(as.numeric(cluster)),
      cluster = ifelse(is.na(cluster_ord), cluster, as.character(cluster_ord))
    ) %>%
    dplyr::select(-cluster_ord)
  
  # Define isoform order and colours
  iso_labels <- unique(frac_long$isoform_id)
  pal <- make_palette(length(iso_labels))
  names(pal) <- iso_labels
  
  # 2) Build bar plot of isoform fractions
  p_frac <- ggplot(frac_long, aes(x = cluster, y = fraction, fill = isoform_id)) +
    geom_bar(stat = "identity") +
    labs(x = "Subcluster", y = "Fraction of expression", fill = "Transcript") +
    scale_fill_manual(values = pal) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none")
  
  # 3) Build transcript structure plot from GTF
  # Use the original transcript_name list to pull exons, then relabel to compatible_tx where present
  tx_from_switch <- unique(sw_gene$transcript_name)
  
  exons_df <- as.data.frame(
    gtf[gtf$type == "exon" & gtf$transcript_id %in% tx_from_switch]
  )
  
  if (nrow(exons_df) == 0) {
    stop(sprintf("No exon features in GTF for the transcripts of gene '%s'.", gene_of_interest))
  }
  
  # Prefer compatible_tx when available, so the structure aligns with the isoform_id used in bars
  exons_df <- exons_df %>%
    dplyr::mutate(transcript_plot_id = ifelse(!is.na(compatible_tx) & compatible_tx != "",
                                              compatible_tx, transcript_id))
  
  # Keep only isoforms present in the fraction plot
  exons_df <- exons_df %>%
    dplyr::filter(transcript_plot_id %in% iso_labels)
  
  if (nrow(exons_df) == 0) {
    stop("After mapping to compatible_tx, none of the transcripts overlap isoform_ids used in the fraction plot.")
  }
  
  # Create introns and rescale
  introns_df <- to_intron(exons_df, "transcript_plot_id")
  rescaled <- shorten_gaps(
    exons = exons_df,
    introns = introns_df,
    group_var = "transcript_plot_id"
  )
  
  rescaled_exons  <- rescaled %>% dplyr::filter(.data$type == "exon")
  rescaled_introns <- rescaled %>% dplyr::filter(.data$type == "intron")
  
  # Structure plot
  p_struct <- ggplot(rescaled_exons, aes(
    xstart = start, xend = end, y = transcript_plot_id, fill = transcript_plot_id
  )) +
    geom_range(height = 0.25) +
    geom_intron(
      data = rescaled_introns,
      aes(strand = strand),
      arrow.min.intron.length = shorten_intron_min
    ) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = "Pseudo-genomic position", y = NULL) +
    ggtitle(gene_of_interest) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0, face = "bold"))
  
  combined <- p_struct / p_frac + plot_layout(heights = c(0.5, 1))
  
  list(
    combined = combined,
    structure_plot = p_struct,
    fraction_plot = p_frac,
    data = list(
      fractions_long = frac_long,
      exons_rescaled = rescaled_exons,
      introns_rescaled = rescaled_introns
    )
  )
}

switch_genes_in_capture <- intersect(isoform_switches$gene_name, unique(capture_isoforms$associated_gene))

# plot for each gene in swtich_genes_in_capture
for (gene in switch_genes_in_capture) {
  res <- plot_gene_switch(
    gtf = gtf,
    switches_table = isoform_switches,
    gene_of_interest = gene,          # <- change this to any gene
    cluster_suffix = "_pct",     # adapt if your table uses a different suffix
    shorten_intron_min = 300
  )
  
  # Show the combined plot
  print(res$combined)
  
  # Optionally save
  ggsave(paste0(gene, "_switch_plot.pdf"), res$combined, width = 8, height = 4)
}

capture_res_dir <- "/Volumes/Projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/merged_capture/FINAL_RESULTS_March_2025/merged.secondpass.filt_split/tama_collapse/sqanti3_filter/results"
# list all "topSwitches_" files across subdirectories
files <- list.files(capture_res_dir, pattern = "^topSwitches_.*\\.tsv$", recursive = TRUE, full.names = TRUE)
# read them all and pull out the gene_name column
all_switch_genes <- map(files, ~ read.delim(.x)$gene_name) %>%
  flatten_chr() %>%           # combine into single vector
  unique()

genes_to_plot <- intersect(switch_genes_in_capture, all_switch_genes)


