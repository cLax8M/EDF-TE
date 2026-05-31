#!/bin/bash
# Same as run_nuclpred_array.sh but SKIPS the DNAshapeR step (step 1)
# when all 5 shape files already exist for a species.
# Used for OOM re-submissions where shapes are already computed.

#SBATCH --job-name=nuclpred
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --array=1-238
#SBATCH --output=/data/carlos/EDF-TE/logs/nuclpred/nuclpred_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nuclpred/nuclpred_%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)" || true
micromamba activate nuclpred

GENOMES_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
NUCLPRED_DIR="/data/carlos/EDF-TE/TE_analysis/NuclPred_v1"
SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"
SHAPES_BASE="/data/carlos/EDF-TE/TE_analysis/nucleosomes/shapes"
SPECIES_LIST="${SCRIPTS_DIR}/species_nuclpred.txt"

mkdir -p "${OUT_DIR}" "${SHAPES_BASE}"
mkdir -p /data/carlos/EDF-TE/logs/nuclpred

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

GENOME_FA="${GENOMES_DIR}/${SP}.fasta"
OUT_WIG="${OUT_DIR}/${SP}_nuclpred.wig"

if [[ -s "${OUT_WIG}" ]]; then
    echo "[SKIP] Already complete: ${SP}"
    exit 0
fi

echo "[START] ${SP}"

SHAPES_DIR="${SHAPES_BASE}/${SP}"
# Check if all 5 shape files exist
SHAPE_COUNT=$(ls "${SHAPES_DIR}"/*.{MGW,EP,Roll,ProT,HelT} 2>/dev/null | wc -l)
if [[ "${SHAPE_COUNT}" -lt 5 ]]; then
    echo "[ERROR] Missing shapes for ${SP} (found ${SHAPE_COUNT}/5). Run full nuclpred script."
    exit 1
fi
echo "[1/3] Shapes already exist (${SHAPE_COUNT}/5), skipping DNAshapeR"

echo "[2/3] Generating NuclPred input matrix..."
MATRIX_FILE="${SHAPES_DIR}/${SP}.matrix"

python3 "${NUCLPRED_DIR}/generate_matrix.py" \
    --assembly "${GENOME_FA}" \
    --data_dir "${SHAPES_DIR}" \
    --organism_scaffold_label "${SP}" \
    > "${MATRIX_FILE}"

echo "[INFO] Matrix rows: $(wc -l < "${MATRIX_FILE}")"

echo "[3/3] Running NuclPred prediction..."
(cd "${NUCLPRED_DIR}" && python3 predict_nucleosomes.py \
    --input "${MATRIX_FILE}") > "${OUT_WIG}"

echo "[DONE] ${SP} → ${OUT_WIG}"

WIG_LINES=$(wc -l < "${OUT_WIG}")
if [[ "${WIG_LINES}" -lt 10 ]]; then
    echo "[WARN] WIG file suspiciously small (${WIG_LINES} lines): ${OUT_WIG}"
fi

rm -f "${MATRIX_FILE}"
