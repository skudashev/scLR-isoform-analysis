library(docopt)
library(Isosceles)

doc <- 'Usage:
  isosceles.R -b <bam_files> -g <gtf_file> -f <genome_fasta_file> -o <output_dir> [--samples <samples>] [--ncpu <ncpu>]

  Options:
    -b <bam_files>          Comma-separated list of input BAM files from wf-single-cell.
    -g <gtf_file>           Input GTF file.
    -f <genome_fasta_file>  Input genome FASTA file.
    -o <output_dir>         Output directory.
    --samples <samples>     Comma-separated list of sample names (must correspond to the BAM files).
    --ncpu <ncpu>           Number of CPUs to use [default: 1].'

opt <- docopt(doc)

# Split the comma-separated BAM files and sample names
bam_files <- unlist(strsplit(opt$b, ","))
samples <- unlist(strsplit(opt$samples, ","))

# Check if the number of BAM files matches the number of sample names
if (length(bam_files) != length(samples)) {
  stop("The number of BAM files and sample names must match.")
}

# Create a named character vector for BAM files
named_bam_files <- setNames(bam_files, samples)

# Input files
gtf_file <- opt$g
genome_fasta_file <- opt$f
output_dir <- opt$o
ncpu <- as.integer(opt$ncpu)
chunk_size <- 1000000

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Prepare transcript data for the analysis
known_intron_motifs <- c("GT-AG", "GC-AG", "AT-AG", "AT-AC")

if (file.exists(file.path(output_dir, "transcript_data.rds"))) {
  transcript_data <- readRDS(file.path(output_dir, "transcript_data.rds"))
} else {
  # Process each BAM file and prepare transcript data
  bam_parsed <- bam_to_read_structures(bam_files, ncpu = ncpu)
  transcript_data <- prepare_transcripts(
    gtf_file = gtf_file,
    genome_fasta_file = genome_fasta_file,
    known_intron_motifs = known_intron_motifs,
    rescue_annotated_introns = TRUE,
    bam_parsed = bam_parsed,
    min_bam_splice_read_count = 1,
    min_bam_splice_fraction = 0.01)
  saveRDS(transcript_data, file.path(output_dir, "transcript_data.rds"))}

message("Transcript data object saved.")

# Create a TCC (Transcript Compatibility Counts) SummarizedExperiment object:
se_tcc_file <- file.path(output_dir, "se_tcc.rds")
if (file.exists(se_tcc_file)) {
  # If it exists, load the TCC SummarizedExperiment object
  se_tcc <- readRDS(se_tcc_file)
} else {
  # If it doesn't exist, run bam_to_tcc() to generate the object
  se_tcc <- bam_to_tcc(
    bam_files = named_bam_files,
    transcript_data = transcript_data,
    run_mode = "de_novo_loose",
    min_read_count = 1,
    extend_spliced_transcripts = 100,
    min_relative_expression = 0,
    is_single_cell = TRUE,
    barcode_tag = "CB",
    ncpu = ncpu,
    chunk_size = chunk_size
  )
  
  # Save the generated TCC SummarizedExperiment object
  saveRDS(se_tcc, se_tcc_file)
}

message("TCC SummarizedExperiment object created.")

em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- FALSE

# Convert TCC to transcript-level SummarizedExperiment
se_transcript <- tcc_to_transcript(
  se_tcc = se_tcc,
  use_length_normalization = use_length_normalization,
  em.maxiter = em.maxiter,
  em.conv = em.conv,
  ncpu = ncpu
)
saveRDS(se_transcript, file.path(output_dir, "se_transcript.rds"))

# Export transcript annotation as GTF file
# check if the transcript annotation GTF file already exists
if (file.exists(file.path(output_dir, "transcript_annotation.gtf"))) {
  message("Transcript annotation GTF file already exists.")
} else {
  # If it doesn't exist, export the transcript annotation as a GTF file
  export_gtf(se_transcript, file.path(output_dir, "transcript_annotation.gtf")) 
}

# Convert TCC to gene-level SummarizedExperiment
se_gene <- tcc_to_gene(
  se_tcc = se_tcc
)
saveRDS(se_gene, file.path(output_dir, "se_gene.rds"))

se_psi <- transcript_to_psi(
  se_transcript, ncpu = ncpu
)
saveRDS(se_psi, file.path(output_dir, "se_psi.rds"))

message("SummarizedExperiment objects saved.")
