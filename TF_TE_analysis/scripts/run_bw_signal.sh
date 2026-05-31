#!/bin/bash
#SBATCH --job-name=bw_TE_signal
#SBATCH --output=/data/carlos/EDF-TE/logs/bw_TE_signal_%j.out
#SBATCH --error=/data/carlos/EDF-TE/logs/bw_TE_signal_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

PYTHON=/data/carlos/EDF-TE/PanGenome/.venv/bin/python
SCRIPT=/data/carlos/EDF-TE/EDF-TE/TF_TE_analysis/scripts/02_bigwig_signal_TEs.py

echo "Start: $(date)"
echo "Node: $(hostname)"

"${PYTHON}" "${SCRIPT}"

echo "End: $(date)"
