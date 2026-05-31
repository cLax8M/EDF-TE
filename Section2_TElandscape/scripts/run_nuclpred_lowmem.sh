#!/bin/bash
# ============================================================
# run_nuclpred_lowmem.sh
# Resubmission for NuclPred jobs that OOM-killed with the original
# 80G script. Uses predict_nucleosomes_lowmem.py, which streams
# the matrix one scaffold at a time instead of loading the entire
# genome into RAM at once.
#
# Failed species (indices in species_nuclpred.txt):
#   1  AP_Amoocci  — completely missing (likely never ran)
#  18  BA_Fomme1   — 63M rows OOM-killed
#  25  BA_Rhisol1  — 51M rows OOM-killed
#  61  CH_Neolan1  — 201M rows OOM-killed
#  80  EN_Zoorad1  — 655M rows OOM-killed
#
# Usage:
#   sbatch --array=1,18,25,61,80 run_nuclpred_lowmem.sh
# ============================================================

#SBATCH --job-name=nuclpred_lm
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --array=1,18,25,61,80
#SBATCH --output=/data/carlos/EDF-TE/logs/nuclpred/nuclpred_lm_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nuclpred/nuclpred_lm_%A_%a.err

set -euo pipefail

# ── environment ──────────────────────────────────────────────
eval "$(micromamba shell hook --shell bash)"
micromamba activate nuclpred

# ── paths ────────────────────────────────────────────────────
GENOMES_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
NUCLPRED_DIR="/data/carlos/EDF-TE/TE_analysis/NuclPred_v1"
SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"
SHAPES_BASE="/data/carlos/EDF-TE/TE_analysis/nucleosomes/shapes"
SPECIES_LIST="${SCRIPTS_DIR}/species_nuclpred.txt"

mkdir -p "${OUT_DIR}" "${SHAPES_BASE}"
mkdir -p /data/carlos/EDF-TE/logs/nuclpred

# ── select species for this array task ───────────────────────
SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

GENOME_FA="${GENOMES_DIR}/${SP}.fasta"
OUT_WIG="${OUT_DIR}/${SP}_nuclpred.wig"

# ── skip if already complete (non-empty WIG) ─────────────────
if [[ -s "${OUT_WIG}" ]]; then
    echo "[SKIP] Already complete: ${SP}"
    exit 0
fi

# Remove any partial (0-byte) output left by a previous OOM crash
if [[ -f "${OUT_WIG}" ]]; then
    echo "[WARN] Removing empty/partial WIG from previous run: ${OUT_WIG}"
    rm -f "${OUT_WIG}"
fi

if [[ ! -f "${GENOME_FA}" ]]; then
    echo "[ERROR] Genome not found: ${GENOME_FA}"
    exit 1
fi

echo "[START] ${SP} (array index ${SLURM_ARRAY_TASK_ID})"

# ── Step 1: DNA shape parameters ──────────────────────────────
echo "[1/3] Computing DNA shape parameters..."
SHAPES_DIR="${SHAPES_BASE}/${SP}"
mkdir -p "${SHAPES_DIR}"

# Skip if shapes already computed (re-use from previous run if present)
# DNAshapeR outputs: {SP}.EP, {SP}.HelT, {SP}.MGW, {SP}.ProT, {SP}.Roll
if [[ ! -f "${SHAPES_DIR}/${SP}.MGW" ]]; then
    micromamba run -n treewasEDFTE Rscript "${SCRIPTS_DIR}/run_dnashape.R" \
        "${GENOME_FA}" "${SHAPES_DIR}"
else
    echo "[1/3] Shapes already exist, skipping DNAshapeR."
fi

# ── Step 2: generate input matrix ────────────────────────────
echo "[2/3] Generating NuclPred input matrix..."
MATRIX_FILE="${SHAPES_DIR}/${SP}.matrix"

if [[ ! -f "${MATRIX_FILE}" ]]; then
    python3 "${NUCLPRED_DIR}/generate_matrix.py" \
        --assembly "${GENOME_FA}" \
        --data_dir "${SHAPES_DIR}" \
        --organism_scaffold_label "${SP}" \
        > "${MATRIX_FILE}"
else
    echo "[2/3] Matrix file already exists, re-using: ${MATRIX_FILE}"
fi

echo "[INFO] Matrix rows: $(wc -l < "${MATRIX_FILE}")"

# ── Step 3: predict nucleosome occupancy (LOW-MEMORY VERSION) ──
# predict_nucleosomes_lowmem.py streams the matrix scaffold-by-scaffold
# instead of loading the full genome into RAM at once. Peak memory
# is bounded by the largest single scaffold, not the whole genome.
echo "[3/3] Running low-memory NuclPred prediction..."
(cd "${NUCLPRED_DIR}" && python3 predict_nucleosomes_lowmem.py \
    --input "${MATRIX_FILE}") > "${OUT_WIG}"

echo "[DONE] ${SP} → ${OUT_WIG}"

# Sanity check
WIG_LINES=$(wc -l < "${OUT_WIG}")
if [[ "${WIG_LINES}" -lt 10 ]]; then
    echo "[WARN] WIG file suspiciously small (${WIG_LINES} lines): ${OUT_WIG}"
fi

# Clean up matrix
rm -f "${MATRIX_FILE}"
