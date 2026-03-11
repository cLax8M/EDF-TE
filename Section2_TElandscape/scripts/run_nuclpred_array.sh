#!/bin/bash
# ============================================================
# run_nuclpred_array.sh
# SLURM array: run NuclPred on all 238 genomes.
#
# Per species:
#   1. Compute DNA shape parameters with DNAshapeR (R)
#   2. Generate input matrix with generate_matrix.py
#   3. Predict nucleosome occupancy with predict_nucleosomes.py
#      → output: {SP}_nuclpred.wig (z-score occupancy)
#
# Requires:
#   conda env nuclpred (from NuclPred_v1/nucleosome_prediction_env.yml)
#   + DNAshapeR: R -e "BiocManager::install('DNAshapeR')"
#   Clone repo: git clone https://github.com/sjmondo/NuclPred_v1 \
#                 /data/carlos/EDF-TE/TE_analysis/NuclPred_v1
#
# Usage:
#   sbatch --array=1-238 run_nuclpred_array.sh
# ============================================================

#SBATCH --job-name=nuclpred
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --array=1-238
#SBATCH --output=/data/carlos/EDF-TE/logs/nuclpred/nuclpred_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nuclpred/nuclpred_%A_%a.err

set -euo pipefail

# ── environment ───────────────────────────────────────────────
conda activate nuclpred

# ── paths ────────────────────────────────────────────────────
GENOMES_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
NUCLPRED_DIR="/data/carlos/EDF-TE/TE_analysis/NuclPred_v1"
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"
SHAPES_BASE="/data/carlos/EDF-TE/TE_analysis/nucleosomes/shapes"
SPECIES_LIST="${SCRIPTS_DIR}/species_nuclpred.txt"

mkdir -p "${OUT_DIR}" "${SHAPES_BASE}"
mkdir -p /data/carlos/EDF-TE/logs/nuclpred

# ── build species list (first run only) ───────────────────────
if [[ ! -f "${SPECIES_LIST}" ]]; then
    ls "${GENOMES_DIR}"/*.fasta 2>/dev/null \
        | xargs -n1 basename \
        | sed 's/\.fasta$//' \
        | sort > "${SPECIES_LIST}"
    echo "[INFO] Species list: ${SPECIES_LIST} ($(wc -l < "${SPECIES_LIST}") species)"
fi

# ── select species for this array task ───────────────────────
SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

GENOME_FA="${GENOMES_DIR}/${SP}.fasta"
OUT_WIG="${OUT_DIR}/${SP}_nuclpred.wig"

# ── skip if already complete ──────────────────────────────────
if [[ -f "${OUT_WIG}" ]]; then
    echo "[SKIP] Already complete: ${SP}"
    exit 0
fi

if [[ ! -f "${GENOME_FA}" ]]; then
    echo "[ERROR] Genome not found: ${GENOME_FA}"
    exit 1
fi

echo "[START] ${SP}"

# ── Step 1: DNA shape parameters ──────────────────────────────
echo "[1/3] Computing DNA shape parameters..."
SHAPES_DIR="${SHAPES_BASE}/${SP}"
mkdir -p "${SHAPES_DIR}"

Rscript "${SCRIPTS_DIR}/run_dnashape.R" "${GENOME_FA}" "${SHAPES_DIR}"

# ── Step 2: generate input matrix ─────────────────────────────
echo "[2/3] Generating NuclPred input matrix..."
MATRIX_FILE="${SHAPES_DIR}/${SP}.matrix"

python3 "${NUCLPRED_DIR}/generate_matrix.py" \
    --assembly "${GENOME_FA}" \
    --data_dir "${SHAPES_DIR}" \
    --organism_scaffold_label "${SP}" \
    > "${MATRIX_FILE}"

echo "[INFO] Matrix rows: $(wc -l < "${MATRIX_FILE}")"

# ── Step 3: predict nucleosome occupancy ──────────────────────
echo "[3/3] Running NuclPred prediction..."
python3 "${NUCLPRED_DIR}/predict_nucleosomes.py" \
    --input "${MATRIX_FILE}" \
    > "${OUT_WIG}"

echo "[DONE] ${SP} → ${OUT_WIG}"

# Sanity check: WIG should have data
WIG_LINES=$(wc -l < "${OUT_WIG}")
if [[ "${WIG_LINES}" -lt 10 ]]; then
    echo "[WARN] WIG file suspiciously small (${WIG_LINES} lines): ${OUT_WIG}"
fi

# Clean up matrix (large intermediate file)
rm -f "${MATRIX_FILE}"
