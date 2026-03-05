#!/bin/bash
#SBATCH --job-name=treewas
#SBATCH --array=0-15
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/data/carlos/EDF-TE/EDF-TE/logs/treewas_%A_%a.out
#SBATCH --error=/data/carlos/EDF-TE/EDF-TE/logs/treewas_%A_%a.err
#SBATCH --partition=normal.168h

# ── treeWAS array job: 8 phenotypes × 2 matrices = 16 runs ──────────────────
#
# Run from repo root:
#   sbatch EDF-TE/TreeWAS/scripts/01_submit_treewas.sh
#
# Ensure inputs exist first:
#   python EDF-TE/TreeWAS/scripts/00_prep_data.py
#
# Array index layout:
#   index = phenotype_idx * 2 + matrix_idx
#   phenotype_idx : 0-7  (see PHENOTYPES array below)
#   matrix_idx    : 0=OG, 1=IPR

set -euo pipefail

PHENOTYPES=(
    "log10_Genome_Size"
    "Pct_Intergenic"
    "Avg_Introns_Multiexon"
    "Avg_Exons_Per_Gene"
    "Pct_paralog_genes"
    "Pct_OGs_multicopy"
    "Percent_unknown"
    "mean_r_intra"
)
MATRICES=("OG" "IPR")

PHENO_IDX=$(( SLURM_ARRAY_TASK_ID / 2 ))
MATRIX_IDX=$(( SLURM_ARRAY_TASK_ID % 2 ))

PHENOTYPE="${PHENOTYPES[$PHENO_IDX]}"
MATRIX="${MATRICES[$MATRIX_IDX]}"

echo "================================================"
echo "Job array index : $SLURM_ARRAY_TASK_ID"
echo "Phenotype        : $PHENOTYPE"
echo "Matrix           : $MATRIX"
echo "Node             : $(hostname)"
echo "Started          : $(date)"
echo "================================================"

# ── paths ────────────────────────────────────────────────────────────────────
# SLURM_SUBMIT_DIR = dir from which sbatch was called (repo root when run correctly)
SCRIPT_DIR="${SLURM_SUBMIT_DIR}/TreeWAS/scripts"
RESULTS_DIR="${SLURM_SUBMIT_DIR}/TreeWAS/results"

mkdir -p "$RESULTS_DIR"

# ── activate environment ──────────────────────────────────────────────────────
export MAMBA_EXE="$(which micromamba)"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate treewasEDFTE

# ── run ──────────────────────────────────────────────────────────────────────
Rscript "${SCRIPT_DIR}/01_run_treewas.R" \
    --phenotype "$PHENOTYPE" \
    --matrix    "$MATRIX" \
    --datadir   "$RESULTS_DIR" \
    --outdir    "$RESULTS_DIR" \
    --n_sim     10000

echo "Finished: $(date)"
