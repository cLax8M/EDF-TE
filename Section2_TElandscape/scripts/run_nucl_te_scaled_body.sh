#!/bin/bash
# ============================================================
# run_nucl_te_scaled_body.sh
# Per-species scaled-body meta-TE profile for ALL 5 variants:
#   nucleosomes, nucleosomes_recent, nucleosomes_FL,
#   nucleosomes_FL_recent, nucleosomes_real6
#
# Resolves WIG (real Nuc-seq if available, else NuclPred) the same way
# as run_nucl_te_interplay_recent.sh.
# ============================================================
#SBATCH --job-name=nucl_te_scaled
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --array=1-218
#SBATCH --exclude=LEGcompute1
#SBATCH --output=/data/carlos/EDF-TE/logs/nucl_te_scaled/nucl_te_scaled_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nucl_te_scaled/nucl_te_scaled_%A_%a.err

set -euo pipefail
mkdir -p /data/carlos/EDF-TE/logs/nucl_te_scaled

SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts"
PY=/data/carlos/EDF-TE/.venv/bin/python3
WIG_DIR="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"
RESULTS="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape"
SPECIES_LIST="${SCRIPTS_DIR}/species_nucl_te.txt"

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

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
[[ -z "${SP}" ]] && { echo "[SKIP] empty"; exit 0; }

NUCLPRED_WIG="${WIG_DIR}/${SP}_nuclpred.wig"
if [[ -v "REAL_WIGS[${SP}]" ]] && [[ -s "${REAL_WIGS[${SP}]}" ]]; then
    WIG="${REAL_WIGS[${SP}]}";  echo "[INFO] real WIG for ${SP}"
elif [[ -s "${NUCLPRED_WIG}" ]]; then
    WIG="${NUCLPRED_WIG}";      echo "[INFO] NuclPred WIG for ${SP}"
else
    echo "[ERROR] no WIG for ${SP}"; exit 1
fi

# 5 variants: <subdir>:<bed_path>
EG_BED="/data/carlos/EDF-TE/TE_analysis/annotated/${SP}_EarlGrey/${SP}_summaryFiles/${SP}.filteredRepeats.bed"
RECENT_BED="/legserv/Temp/Carlos/EDF-TE/data/recent_beds/${SP}.filteredRepeats.recent.bed"
FL_BED="/legserv/Temp/Carlos/EDF-TE/data/fl_beds/${SP}.filteredRepeats.fl.bed"
FL_REC_BED="/legserv/Temp/Carlos/EDF-TE/data/recent_beds/${SP}.filteredRepeats.fl_recent.bed"

# Real6 only runs the all-TEs view (its own dir).
declare -a VARIANTS=(
  "nucleosomes:${EG_BED}"
  "nucleosomes_recent:${RECENT_BED}"
  "nucleosomes_FL:${FL_BED}"
  "nucleosomes_FL_recent:${FL_REC_BED}"
)
case "${SP}" in
  BA_Mrafri1|BA_Mereb1|AS_Usnflo1|BL_Catan2|BA_Rhoto|MU_Rhimi59)
      VARIANTS+=("nucleosomes_real6:${EG_BED}") ;;
esac

for entry in "${VARIANTS[@]}"; do
    sub="${entry%%:*}"; bed="${entry#*:}"
    if [[ ! -s "${bed}" ]]; then
        echo "[skip ${sub}] no BED ${bed}"; continue
    fi
    OUT_DIR="${RESULTS}/${sub}/scaled_body"
    OUT="${OUT_DIR}/${SP}_scaled.tsv"
    if [[ -s "${OUT}" ]]; then echo "[exists ${sub}] ${OUT}"; continue; fi
    mkdir -p "${OUT_DIR}"
    echo "[$(date)] ${SP} ${sub}"
    ${PY} "${SCRIPTS_DIR}/F9_scaled_body_meta.py" \
        "${WIG}" "${bed}" "${SP}" "${OUT}" \
        --flank 500 --flank_bins 50 --body_bins 100 \
        --min_len 200 --max_len 20000
done
echo "[$(date)] ALL VARIANTS DONE for ${SP}"
