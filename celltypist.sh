#!/bin/bash
#SBATCH --job-name=celltypist
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com
#SBATCH -p ei-short
#SBATCH -N 1
#SBATCH -c 5
#SBATCH -t 02:00:00
#SBATCH --mem=15gb
#SBATCH -o logs/celltypist.%N.%j.out
#SBATCH -e logs/celltypist.%N.%j.err

# Usage: sbatch celltypist.sh <model_name.pkl> <Seurat_RDS> <t2g.tsv>
model_file=$1
seurat_rds=$2
t2g_file=$3

# Paths
workdir=$(pwd)
celltypist_img="$HOME/singularity/images/celltypist.sif"
prep_script="prep_celltypist.R"

# Derive prefix from model file (remove extension, path)
prefix=$(basename "$model_file" .pkl)

# Step 1: Prepare input matrix
Rscript $prep_script "$seurat_rds" "$t2g_file" "${prefix}_gene_counts.csv"

# Step 2: Run CellTypist
singularity exec "$celltypist_img" celltypist \
    --model "$model_file" \
    --input "${prefix}_gene_counts.csv" \
    --majority-voting \
    --output "${prefix}_celltypist_results"
