#!/bin/bash
# ============================================================
# run_nucl_te_interplay_FL.sh
# SLURM array: nucleosome / TE interplay restricted to FULL-LENGTH
# copies (copy_len >= 0.8 * consensus_len).
#
# Mirrors run_nucl_te_interplay.sh but uses pre-computed FL BEDs
# from /legserv/Temp/Carlos/EDF-TE/data/fl_beds/ and writes to a
# separate output tree (*_FL/) so results are not overwritten.
#
# Usage:
#   sbatch --array=1-218 run_nucl_te_interplay_FL.sh
# ============================================================

#SBATCH --job-name=nucl_te_FL
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=1-218
#SBATCH --exclude=LEGcompute1
#SBATCH --output=/data/carlos/EDF-TE/logs/nucl_te_FL/nucl_te_FL_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nucl_te_FL/nucl_te_FL_%A_%a.err

set -euo pipefail

SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts"
WIG_DIR="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"
FL_BED_DIR="/legserv/Temp/Carlos/EDF-TE/data/fl_beds"
OUT_BASE="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes_FL"
SPECIES_LIST="${SCRIPTS_DIR}/species_nucl_te.txt"

mkdir -p /data/carlos/EDF-TE/logs/nucl_te_FL
mkdir -p "${OUT_BASE}/f3_bias" "${OUT_BASE}/f5_ndr" "${OUT_BASE}/f6_meta"

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
[[ -z "${SP}" ]] && { echo "[SKIP] no species at ${SLURM_ARRAY_TASK_ID}"; exit 0; }

BED="${FL_BED_DIR}/${SP}.filteredRepeats.fl.bed"

# WIG selection — same dispatch table as the main runner
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
NUCLPRED_WIG="${WIG_DIR}/${SP}_nuclpred.wig"

if [[ -v "REAL_WIGS[${SP}]" ]] && [[ -s "${REAL_WIGS[${SP}]}" ]]; then
    WIG="${REAL_WIGS[${SP}]}"; echo "[INFO] real WIG for ${SP}"
elif [[ -s "${NUCLPRED_WIG}" ]]; then
    WIG="${NUCLPRED_WIG}"; echo "[INFO] NuclPred WIG for ${SP}"
else
    echo "[ERROR] no WIG for ${SP}"; exit 1
fi

if [[ ! -s "${BED}" ]]; then
    echo "[ERROR] missing FL BED: ${BED}"; exit 1
fi

echo "[$(date)] START ${SP} (FL-only)"

F3_OUT="${OUT_BASE}/f3_bias/${SP}_f3_bias.tsv"
if [[ ! -s "${F3_OUT}" ]]; then
    python3 "${SCRIPTS_DIR}/F3_te_nucleosome_bias.py" "${WIG}" "${BED}" "${SP}" "${F3_OUT}"
fi

F5_OUT="${OUT_BASE}/f5_ndr/${SP}_f5_ndr.tsv"
if [[ ! -s "${F5_OUT}" ]]; then
    python3 "${SCRIPTS_DIR}/F5_ndr_te_density.py" "${WIG}" "${BED}" "${SP}" "${F5_OUT}" \
        --ndr_pct 25 --nor_pct 75 --min_te_bp 200
fi

F6_OUT_DIR="${OUT_BASE}/f6_meta"
if [[ ! -s "${F6_OUT_DIR}/${SP}_meta_te_profiles.tsv" ]]; then
    python3 "${SCRIPTS_DIR}/F6_meta_te_profile.py" "${WIG}" "${BED}" "${SP}" "${F6_OUT_DIR}" \
        --window 2000 --min_insertions 10 --max_te_len 5000
fi

echo "[$(date)] DONE ${SP} (FL)"
