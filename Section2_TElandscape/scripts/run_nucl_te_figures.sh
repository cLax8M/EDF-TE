#!/bin/bash
# ============================================================
# run_nucl_te_figures.sh
# Post-array: combine per-species F3/F5/F6 outputs and generate
# multi-species figures with F7_plot_nucleosome_analysis.py.
#
# Run after array job completes:
#   sbatch --dependency=afterok:<ARRAY_JOB_ID> run_nucl_te_figures.sh
# ============================================================

#SBATCH --job-name=nucl_te_figures
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/nucl_te/nucl_te_figures_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nucl_te/nucl_te_figures_%j.err

set -euo pipefail

SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts"
OUT_BASE="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes"
FIG_DIR="${OUT_BASE}/figures"
mkdir -p "${FIG_DIR}"

echo "[$(date)] Combining per-species F3 outputs..."
COMBINED_F3="${OUT_BASE}/f3_bias_all_species.tsv"
# Write header from first file, then data from all
head -1 "$(ls ${OUT_BASE}/f3_bias/*_f3_bias.tsv | head -1)" > "${COMBINED_F3}"
for f in "${OUT_BASE}/f3_bias/"*_f3_bias.tsv; do
    tail -n +2 "$f" >> "${COMBINED_F3}"
done
echo " → $(wc -l < "${COMBINED_F3}") rows in combined F3"

echo "[$(date)] Combining per-species F5 outputs..."
COMBINED_F5="${OUT_BASE}/f5_ndr_all_species.tsv"
head -1 "$(ls ${OUT_BASE}/f5_ndr/*_f5_ndr.tsv | head -1)" > "${COMBINED_F5}"
for f in "${OUT_BASE}/f5_ndr/"*_f5_ndr.tsv; do
    tail -n +2 "$f" >> "${COMBINED_F5}"
done
echo " → $(wc -l < "${COMBINED_F5}") rows in combined F5"

echo "[$(date)] Combining per-species F6 meta-profile outputs..."
COMBINED_F6="${OUT_BASE}/f6_meta_all_species.tsv"
head -1 "$(ls ${OUT_BASE}/f6_meta/*_meta_te_profiles.tsv | head -1)" > "${COMBINED_F6}"
for f in "${OUT_BASE}/f6_meta/"*_meta_te_profiles.tsv; do
    tail -n +2 "$f" >> "${COMBINED_F6}"
done
echo " → $(wc -l < "${COMBINED_F6}") rows in combined F6"

# Also update the legacy nucleosome_bias_summary.tsv for compatibility
echo "[$(date)] Updating legacy nucleosome_bias_summary.tsv..."
cp "${COMBINED_F3}" "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosome_bias_summary.tsv"
echo " → nucleosome_bias_summary.tsv updated"

echo "[$(date)] Running F7 figures..."
source /data/carlos/EDF-TE/PanGenome/.venv/bin/activate

python3 "${SCRIPTS_DIR}/F7_plot_nucleosome_analysis.py" \
    --f3 "${COMBINED_F3}" \
    --f5 "${COMBINED_F5}" \
    --f6 "${COMBINED_F6}" \
    --output_dir "${FIG_DIR}"

echo "[$(date)] Figures written to ${FIG_DIR}"
ls -lh "${FIG_DIR}/"
