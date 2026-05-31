#!/bin/bash
#SBATCH --job-name=gremlin_centromere
#SBATCH --partition=normal.168h
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/gremlin_centromere_%j.out
#SBATCH --error=/data/carlos/EDF-TE/logs/gremlin_centromere_%j.err

echo "Start: $(date)"
echo "Node: $(hostname)"

module load BEDTools/2.31.1-GCC-13.3.0

PYTHON=/data/carlos/EDF-TE/PanGenome/.venv/bin/python
SCRIPT=/data/carlos/EDF-TE/EDF-TE/TF_TE_analysis/scripts/05_gremlin_centromere.py

$PYTHON "$SCRIPT"

echo "End: $(date)"
