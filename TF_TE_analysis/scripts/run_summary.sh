#!/bin/bash
#SBATCH --job-name=TF_TE_summary
#SBATCH --output=/data/carlos/EDF-TE/logs/TF_TE_summary_%j.out
#SBATCH --error=/data/carlos/EDF-TE/logs/TF_TE_summary_%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail

PYTHON=/data/carlos/EDF-TE/PanGenome/.venv/bin/python
SCRIPTS=/data/carlos/EDF-TE/EDF-TE/TF_TE_analysis/scripts
OUTDIR=/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis

echo "Start: $(date)"
echo "Node: $(hostname)"

# Merge FIMO per-class TSVs into final enrichment table + plots
echo "=== Merging FIMO results ==="
"${PYTHON}" "${SCRIPTS}/03_fimo_motif_enrichment.py" \
    --merge-only \
    --outdir "${OUTDIR}/03_fimo_enrichment"

# Generate integrated summary figure
echo "=== Generating summary figure ==="
"${PYTHON}" "${SCRIPTS}/04_summary_figure.py"

echo "End: $(date)"
echo "Outputs at: ${OUTDIR}/figures/TF_TE_summary.pdf"
