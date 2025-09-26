#!/bin/bash

# -------------------------------------------------------------------------
# Script: match_ujc.sh
# Purpose: Compare unannotated junction chains (UJCs) between two GTF files
# Inputs : gtf1 gtf2 prefix1 prefix2
# Outputs: BED12 files and matched UJC table
# Author : Sofia Kudasheva
# -------------------------------------------------------------------------
set -e
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478 #bedtools
source package 71de5a7a-135b-417a-8de1-ede16dc52660 #UCSC tools

# --- Inputs ---
gtf1="$1"
gtf2="$2"
prefix1="$3"
prefix2="$4"

# --- Scripts and env ---
scripts="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/scripts"
python_script="$scripts/find_matched_UJCs.py"
conda_env="lr_analysis"

# --- UCSC conversion tools (assumed available from loaded module/source) ---
# Tools: gtfToGenePred, genePredToBed

# --- Output files ---
gp1="${prefix1}_1.genePred"
gp2="${prefix2}_2.genePred"
bed1="${prefix1}_1.bed12"
bed2="${prefix2}_2.bed12"
matched="${prefix1}_${prefix2}_matched.txt"

echo "[INFO] Converting GTF to BED12..."
gtfToGenePred "$gtf1" "$gp1"
gtfToGenePred "$gtf2" "$gp2"
genePredToBed "$gp1" "$bed1"
genePredToBed "$gp2" "$bed2"

echo "[INFO] Running junction chain comparison..."
source ~/.bashrc
conda activate "$conda_env"
python "$python_script" "$bed1" "$bed2" "$matched"
