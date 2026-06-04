#!/bin/bash
#SBATCH --job-name=multicopy_te_prox
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/cooption/16_multicopy_te_%j.out
#SBATCH --error=/data/carlos/EDF-TE/logs/cooption/16_multicopy_te_%j.err

set -euo pipefail

mkdir -p /data/carlos/EDF-TE/logs/cooption

# Activate the project venv (has scipy, matplotlib, numpy)
PYTHON=/data/carlos/EDF-TE/.venv/bin/python3

# Ensure bedtools is in PATH via syntenyEDFTE env
eval "$(micromamba shell hook --shell bash)"
micromamba activate syntenyEDFTE

echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Python: $PYTHON"
bedtools --version

# Run for all 238 species (will skip those missing GFF/TE bed)
$PYTHON /data/carlos/EDF-TE/EDF-TE/Cooption/scripts/16_multicopy_te_proximity.py

echo "Done: $(date)"
