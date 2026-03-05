#!/bin/bash
#SBATCH --job-name=treewas_prep
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/data/carlos/EDF-TE/EDF-TE/logs/treewas_prep_%j.out
#SBATCH --error=/data/carlos/EDF-TE/EDF-TE/logs/treewas_prep_%j.err
#SBATCH --partition=long

# Prepare treeWAS input matrices.
# Submits 01_submit_treewas.sh as a dependent job after completion.
#
# Usage (from repo root):
#   sbatch EDF-TE/TreeWAS/scripts/00_submit_prep.sh
#
# To run WITHOUT auto-submitting the array:
#   sbatch EDF-TE/TreeWAS/scripts/00_submit_prep.sh --no-submit

set -euo pipefail

AUTO_SUBMIT=true
for arg in "$@"; do
  [[ "$arg" == "--no-submit" ]] && AUTO_SUBMIT=false
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "================================================"
echo "Job: treewas_prep"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "================================================"

# ── activate environment ──────────────────────────────────────────────────────
export MAMBA_EXE="$(which micromamba)"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate treewasEDFTE

# ── run prep ─────────────────────────────────────────────────────────────────
python "${SCRIPT_DIR}/00_prep_data.py"

echo "Prep finished: $(date)"

# ── auto-submit treeWAS array ─────────────────────────────────────────────────
if [[ "$AUTO_SUBMIT" == "true" ]]; then
  echo "Submitting treeWAS array job..."
  JOB_ID=$(sbatch --parsable "${SCRIPT_DIR}/01_submit_treewas.sh")
  echo "Submitted treeWAS array: job ${JOB_ID}"
else
  echo "Skipping auto-submit (--no-submit flag set)"
  echo "To run manually: sbatch TreeWAS/scripts/01_submit_treewas.sh"
fi
