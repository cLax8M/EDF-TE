#!/bin/bash
#SBATCH --job-name=treewas_annot
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=/data/carlos/EDF-TE/EDF-TE/logs/treewas_annot_%j.out
#SBATCH --error=/data/carlos/EDF-TE/EDF-TE/logs/treewas_annot_%j.err
#SBATCH --partition=short

# Annotate treeWAS results after the array job completes.
#
# Usage (from repo root):
#   sbatch EDF-TE/TreeWAS/scripts/02_submit_annotate.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "================================================"
echo "Job: treewas_annot"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "================================================"

# ── activate environment ──────────────────────────────────────────────────────
export MAMBA_EXE="$(which micromamba)"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate treewasEDFTE

# ── run annotation ────────────────────────────────────────────────────────────
python "${SCRIPT_DIR}/02_annotate_results.py"

echo "Annotation finished: $(date)"
echo "Results: ${SCRIPT_DIR}/../results/treewas_annotated_hits.tsv"
