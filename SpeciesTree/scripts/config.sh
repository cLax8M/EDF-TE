#!/usr/bin/env bash
# =============================================================================
# config.sh — Centralized configuration for the phylogenomic pipeline
#
# PURPOSE : Defines all paths, thresholds, and SLURM parameters used
#           throughout the 10-step species tree pipeline.
#
# USAGE   : Source this file from any pipeline script:
#               source "$(dirname "$0")/config.sh"
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# REF     : Liu et al. 2024 (Opisthokonta, 348 taxa, 440 genes)
#           Szánthó et al. 2025 (Fungi timetree, 110 taxa, 225 genes)
#           Merényi et al. 2023 (Fungal genomes, 272 genes)
# =============================================================================

# ─── INPUT DATA ───────────────────────────────────────────────────────────────
# OrthoFinder 3.1.0 results (complete run on 238 fungal+outgroup species)
ORTHO_DIR="/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run"

# Orthogroup gene counts per species (used for fast occupancy + copy-number filtering)
OG_COUNTS_FILE="${ORTHO_DIR}/Orthogroups/Orthogroups.GeneCount.tsv"

# Orthogroups table with gene IDs per species (used for sequence extraction)
OG_TABLE_FILE="${ORTHO_DIR}/Orthogroups/Orthogroups.tsv"

# Per-orthogroup FASTA sequences (unaligned, one file per OG)
OG_SEQ_DIR="${ORTHO_DIR}/Orthogroup_Sequences"

# OrthoFinder species tree (used as starting tree in final IQTree run)
ORTHO_SPECIES_TREE="${ORTHO_DIR}/Species_Tree/SpeciesTree_rooted.txt"

# ─── WORKING DIRECTORY ────────────────────────────────────────────────────────
# All intermediate + final outputs go here
WORK_DIR="/data/carlos/EDF-TE/SpeciesTree"

# Step-by-step output subdirectories (created automatically by each script)
DIR_MARKERS="${WORK_DIR}/01_markers"           # unaligned HOG FASTAs
DIR_ALIGNED="${WORK_DIR}/02_aligned"           # MAFFT alignments
DIR_TRIMMED="${WORK_DIR}/03_trimmed"           # trimAl trimmed alignments
DIR_FILTERED="${WORK_DIR}/04_filtered"         # PhyKIT quality-filtered
DIR_GENE_TREES="${WORK_DIR}/05_gene_trees"     # IQTree gene trees
DIR_FINAL_GENES="${WORK_DIR}/06_final_genes"   # paralog-filtered, 1 seq/species
DIR_CONCAT="${WORK_DIR}/07_concatenated"       # supermatrix + partition file
DIR_GUIDE_TREE="${WORK_DIR}/08_guide_tree"     # LG+F+G4 guide tree
DIR_FINAL_TREE="${WORK_DIR}/09_final_tree"     # LG+C60+F+G4 final species tree

# ─── FILTER THRESHOLDS ────────────────────────────────────────────────────────
# Minimum number of species a HOG must cover (≥40% of 238 species)
# Refs: Liu et al. 2024, Szánthó et al. 2025, Merényi et al. 2023
MIN_SPECIES=95

# Maximum gene copies allowed per species per HOG (allows terminal duplications)
# Paralogs will be removed downstream in step 06
MAX_COPIES=2

# Minimum alignment length (AA columns) after trimming
# Ref: Merényi et al. 2023
MIN_ALN_LEN=60

# Saturation threshold (slope of observed vs. inferred distances)
# HOGs with slope > threshold are excluded due to substitution saturation
# Ref: Szánthó et al. 2025
SAT_THRESHOLD=0.24

# Minimum parsimony-informative sites per alignment
# Below this, the gene contributes little phylogenetic signal
MIN_PIS=25

# Maximum Relative Composition Frequency Variability (compositional bias)
MAX_RCFV=0.1

# Branch length z-score threshold for paralog detection
# Branches with |z| > threshold are flagged as potential ancient paralogs
# Ref: Merényi et al. 2023
PARALOG_ZSCORE=3.5

# ─── COMPUTE RESOURCES ────────────────────────────────────────────────────────
PARTITION="normal.168h"     # use normal.1000h for jobs >7 days           # SLURM partition name
THREADS=8                    # threads for MAFFT, IQTree, trimAl
MAFFT_THREADS=4              # lower for step 2 (many parallel jobs)
MEM_MAFFT="8G"               # memory per MAFFT job
MEM_IQTREE_GENE="16G"        # memory per gene tree IQTree job
MEM_IQTREE_FINAL="128G"      # memory for final species tree (LG+C60)
MEM_IQTREE_GUIDE="64G"       # memory for guide tree
TIME_MAFFT="02:00:00"        # wall time per MAFFT alignment
TIME_TRIMAL="00:30:00"       # wall time per trimAl job
TIME_GENE_TREE="04:00:00"    # wall time per gene tree
TIME_GUIDE_TREE="24:00:00"   # wall time for guide tree
TIME_FINAL_TREE="96:00:00"   # wall time for final tree (LG+C60+F+G4)

# ─── CONDA ENVIRONMENT ────────────────────────────────────────────────────────
CONDA_ENV="esptreeEDFTE"

# ─── TOOL VERSIONS (for reproducibility) ─────────────────────────────────────
# Python 3.10.19, MAFFT v7.526, trimAl v1.4+, IQTree v3.0.1, PhyKIT v2.1.2
# Package manager: micromamba

# ─── EXPECTED GENE COUNTS (for sanity checks) ─────────────────────────────────
# Step 01: ~500–700 HOGs extracted
# Step 04: ~300–400 after PhyKIT filtering
# Step 06: ~225–300 after paralog removal (target per published studies)

echo "[config] Configuration loaded."
echo "  ORTHO_DIR  : ${ORTHO_DIR}"
echo "  WORK_DIR   : ${WORK_DIR}"
echo "  MIN_SPECIES: ${MIN_SPECIES} (≥40% of 238)"
echo "  MAX_COPIES : ${MAX_COPIES}"
echo "  PARTITION  : ${PARTITION}"
