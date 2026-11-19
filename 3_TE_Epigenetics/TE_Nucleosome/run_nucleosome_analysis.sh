#!/bin/bash
# =============================================================================
# runs nucleosome analysis scripts across 217 fungal species

# executes SLURM array jobs for the EDF-TE nucleosome analysis pipeline
=============================================================================

#SBATCH --job-name=nucl_analysis
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=1-217
#SBATCH --output=/path/to/logs/nucl_analysis/nucl_%A_%a.log
#SBATCH --error=/path/to/logs/nucl_analysis/nucl_%A_%a.err

set -euo pipefail

# ── Paths (edit) ──────────────────────────────────────────────────────────
SCRIPTS_DIR="/path/to/EDF-TE_complete/section3_te_epigenetics/nucleosome_analysis"
WIG_DIR="/path/to/results/Section2_TElandscape/nucleosomes/nuclpred"
TE_DIR="/path/to/results/TE_annotation"
FAI_DIR="/path/to/data/genomes"
OUT_BASE="/path/to/results/Section2_TElandscape/nucleosomes"
SPECIES_LIST="${SCRIPTS_DIR}/species_nucl_te.txt"

mkdir -p "${OUT_BASE}/f1_nuclpred"
mkdir -p "${OUT_BASE}/f3_spacing"
mkdir -p "${OUT_BASE}/f5_meta"

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

WIG="${WIG_DIR}/${SP}_nuclpred.wig"
FAI="${FAI_DIR}/${SP}.fasta.fai"
BED="${TE_DIR}/${SP}_EarlGrey/${SP}_summaryFiles/${SP}.filteredRepeats.bed"

if [[ ! -s "${WIG}" ]] || [[ ! -s "${BED}" ]]; then
    echo "[SKIP] ${SP}: missing WIG or TE BED"
    exit 0
fi

echo "[$(date)] START ${SP}"

# ── 01: NuclPred window metrics + chromatin states + TE bias ─────────────
OUT1="${OUT_BASE}/f1_nuclpred/${SP}"
if [[ ! -s "${OUT1}_te_bias.tsv" ]]; then
    python3 "${SCRIPTS_DIR}/01_nuclpred_analysis.py" \
        --wig "${WIG}" --fai "${FAI}" --te_bed "${BED}" \
        --species "${SP}" --win 500 --out_dir "${OUT1}"
fi

# ── 03: NRL / spacing (10 kb windows) ────────────────────────────────────
OUT3="${OUT_BASE}/f3_spacing/${SP}"
if [[ ! -s "${OUT3}_spacing_summary.tsv" ]]; then
    python3 "${SCRIPTS_DIR}/03_nucleosome_spacing.py" \
        --wig "${WIG}" --fai "${FAI}" --te_bed "${BED}" \
        --species "${SP}" --win 10000 --min_te_frac 0.5 --out_dir "${OUT3}"
fi

# ── 05: meta-TE profiles ─────────────────────────────────────────────────
OUT5="${OUT_BASE}/f5_meta/${SP}"
if [[ ! -s "${OUT5}_meta_te_summary.tsv" ]]; then
    python3 "${SCRIPTS_DIR}/05_mac_te_profiles.py" \
        --wig "${WIG}" --fai "${FAI}" --te_bed "${BED}" \
        --species "${SP}" --window 2000 --min_insertions 10 \
        --max_te_len 5000 --out_dir "${OUT5}"
fi

echo "[$(date)] DONE ${SP}"
