# singularity exec ~/singularity/images/scdblfinder.sif R
library(SingleCellExperiment)
library(scDblFinder)
# Define paths
data_dir <- "/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/isosceles_isoquant"
out_dir <- "/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/isosceles_isoquant"
se_transcript_file <- file.path(data_dir, "se_transcript.rds")

# Load Isosceles SummarizedExperiment objects
se_transcript <- readRDS(se_transcript_file)
sce <- as(se_transcript, "SingleCellExperiment")

sce <- scDblFinder(sce)
table(sce$scDblFinder.class)
sce <- sce[, sce$scDblFinder.class == "singlet"]