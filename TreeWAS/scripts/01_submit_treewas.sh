#!/bin/bash
#SBATCH --job-name=treewas
#SBATCH --array=0-15
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/data/carlos/EDF-TE/EDF-TE/logs/treewas_%A_%a.out
#SBATCH --error=/data/carlos/EDF-TE/EDF-TE/logs/treewas_%A_%a.err
#SBATCH --partition=long

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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/../results"

mkdir -p "$RESULTS_DIR"

# ── activate R environment if needed ─────────────────────────────────────────
# If treeWAS is installed in a conda/mamba env, activate it here:
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate treewas_env
# Or load a module:
# module load R/4.3.0

# ── run ──────────────────────────────────────────────────────────────────────
Rscript "${SCRIPT_DIR}/01_run_treewas.R" \
    --phenotype "$PHENOTYPE" \
    --matrix    "$MATRIX" \
    --datadir   "$RESULTS_DIR" \
    --outdir    "$RESULTS_DIR" \
    --n_sim     10000

echo "Finished: $(date)"
