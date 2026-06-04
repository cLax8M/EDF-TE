#!/bin/bash
# ============================================================
# run_fimo_cross_species_figures.sh
# Post-array: merge per-species FIMO results and generate figures.
# Run after run_fimo_cross_species.sh array completes.
#
# Usage:
#   sbatch --dependency=afterok:<ARRAY_JOB_ID> run_fimo_cross_species_figures.sh
# ============================================================

#SBATCH --job-name=fimo_figures
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/carlos/EDF-TE/logs/fimo_cross/fimo_figures_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/fimo_cross/fimo_figures_%j.err

set -euo pipefail

SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/TF_TE_analysis/scripts"
RESULTS_DIR="/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/cross_species_fimo"
FIG_DIR="${RESULTS_DIR}/figures"
OUT_TSV="${RESULTS_DIR}/cross_species_fimo_results.tsv"

mkdir -p "${FIG_DIR}"

# Use the PanGenome venv which has pandas, matplotlib, scipy, numpy
VENV="/data/carlos/EDF-TE/PanGenome/.venv/bin/activate"
if [[ -f "${VENV}" ]]; then
    source "${VENV}"
fi

echo "[$(date)] Species in results: $(tail -n +2 ${OUT_TSV} | cut -f1 | sort -u | wc -l)"
echo "[$(date)] Total rows: $(wc -l < ${OUT_TSV})"

echo "[$(date)] Generating figures..."
python3 "${SCRIPTS_DIR}/plot_hmgb3_cross_species.py" \
    --input  "${OUT_TSV}" \
    --outdir "${FIG_DIR}"

echo "[$(date)] DONE — figures in ${FIG_DIR}"
