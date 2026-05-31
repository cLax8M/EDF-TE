#!/bin/bash
# ============================================================
# run_nucl_te_interplay_real6.sh
# SLURM array: nucleosome / TE interplay analysis restricted to the
# 6 species with REAL experimental MNase / Nuc-seq data
# (parallel companion to the full 218-species run).
#
# Output is written under nucleosomes_real6/ so it never overwrites
# the main run.
#
# Usage:
#   sbatch run_nucl_te_interplay_real6.sh
# ============================================================

#SBATCH --job-name=nucl_te_real6
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=1-6
#SBATCH --exclude=LEGcompute1
#SBATCH --output=/data/carlos/EDF-TE/logs/nucl_te/nucl_te_real6_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nucl_te/nucl_te_real6_%A_%a.err

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────
SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts"
TE_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_BASE="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes_real6"
SPECIES_LIST="${SCRIPTS_DIR}/species_real6.txt"

mkdir -p /data/carlos/EDF-TE/logs/nucl_te
mkdir -p "${OUT_BASE}/f3_bias"
mkdir -p "${OUT_BASE}/f5_ndr"
mkdir -p "${OUT_BASE}/f6_meta"

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# ── Real WIG dispatch (must match main runner) ────────────────
REAL_WIG_DIR="/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq"
RHIMI_REAL_WIG="/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"

declare -A REAL_WIGS=(
  [BA_Mrafri1]="${REAL_WIG_DIR}/BA_Mrafri1/BA_Mrafri1_danpos/pooled/BA_Mrafri1_dedup.smooth.wig"
  [BA_Mereb1]="${REAL_WIG_DIR}/BA_Mereb1/BA_Mereb1_danpos/pooled/BA_Mereb1_dedup.smooth.wig"
  [AS_Usnflo1]="${REAL_WIG_DIR}/AS_Usnflo1/AS_Usnflo1_danpos/pooled/AS_Usnflo1_dedup.smooth.wig"
  [BL_Catan2]="${REAL_WIG_DIR}/BL_Catan2/BL_Catan2_danpos/pooled/BL_Catan2_dedup.smooth.wig"
  [BA_Rhoto]="${REAL_WIG_DIR}/BA_Rhoto/BA_Rhoto_danpos/pooled/BA_Rhoto_dedup.smooth.wig"
  [MU_Rhimi59]="${RHIMI_REAL_WIG}"
)

if [[ ! -v "REAL_WIGS[${SP}]" ]]; then
    echo "[ERROR] No real WIG mapped for ${SP}"
    exit 1
fi
WIG="${REAL_WIGS[${SP}]}"
if [[ ! -s "${WIG}" ]]; then
    echo "[ERROR] Real WIG missing/empty: ${WIG}"
    exit 1
fi
echo "[INFO] ${SP}: real WIG = ${WIG}"

BED="${TE_DIR}/${SP}_EarlGrey/${SP}_summaryFiles/${SP}.filteredRepeats.bed"
if [[ ! -s "${BED}" ]]; then
    echo "[ERROR] Missing TE BED: ${BED}"
    exit 1
fi

echo "[$(date)] START ${SP}"

# ── F3: TE nucleosome bias ─────────────────────────────────────
F3_OUT="${OUT_BASE}/f3_bias/${SP}_f3_bias.tsv"
if [[ ! -s "${F3_OUT}" ]]; then
    python3 "${SCRIPTS_DIR}/F3_te_nucleosome_bias.py" \
        "${WIG}" "${BED}" "${SP}" "${F3_OUT}"
fi

# ── F5: NDR/NOR TE density enrichment ─────────────────────────
F5_OUT="${OUT_BASE}/f5_ndr/${SP}_f5_ndr.tsv"
if [[ ! -s "${F5_OUT}" ]]; then
    python3 "${SCRIPTS_DIR}/F5_ndr_te_density.py" \
        "${WIG}" "${BED}" "${SP}" "${F5_OUT}" \
        --ndr_pct 25 --nor_pct 75 --min_te_bp 200
fi

# ── F6: meta-TE occupancy profiles ────────────────────────────
F6_OUT_DIR="${OUT_BASE}/f6_meta"
F6_DONE="${F6_OUT_DIR}/${SP}_meta_te_profiles.tsv"
if [[ ! -s "${F6_DONE}" ]]; then
    python3 "${SCRIPTS_DIR}/F6_meta_te_profile.py" \
        "${WIG}" "${BED}" "${SP}" "${F6_OUT_DIR}" \
        --window 2000 --min_insertions 10 --max_te_len 5000
fi

echo "[$(date)] DONE ${SP}"
