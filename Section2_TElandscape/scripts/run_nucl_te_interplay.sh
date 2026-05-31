#!/bin/bash
# ============================================================
# run_nucl_te_interplay.sh
# SLURM array: nucleosome / TE interplay analysis for 216 species.
#
# Per species (F3 + F5 + F6):
#   F3: mean nucleosome occupancy over TE insertions vs genome background
#   F5: TE density enrichment in NDR / Intermediate / NOR zones
#   F6: meta-TE occupancy profiles (±2 kb around each TE insertion)
#
# After all per-species jobs complete, run F7 for figures:
#   sbatch --dependency=afterok:<ARRAY_JOB_ID> run_nucl_te_figures.sh
#
# Usage:
#   sbatch --array=1-216 run_nucl_te_interplay.sh
# ============================================================

#SBATCH --job-name=nucl_te_interplay
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=1-218
#SBATCH --exclude=LEGcompute1
#SBATCH --output=/data/carlos/EDF-TE/logs/nucl_te/nucl_te_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nucl_te/nucl_te_%A_%a.err

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────
SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts"
WIG_DIR="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"
TE_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_BASE="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes"
SPECIES_LIST="${SCRIPTS_DIR}/species_nucl_te.txt"

mkdir -p /data/carlos/EDF-TE/logs/nucl_te
mkdir -p "${OUT_BASE}/f3_bias"
mkdir -p "${OUT_BASE}/f5_ndr"
mkdir -p "${OUT_BASE}/f6_meta"

# ── Select species for this task ──────────────────────────────
SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

NUCLPRED_WIG="${WIG_DIR}/${SP}_nuclpred.wig"
# Default BED path. Some species have a name mismatch between WIG and EarlGrey dir:
#   AS_AspfumUD1 → WIG genome name, but EarlGrey run as AS_AspfumUD
#   NU_Parvatl_c → WIG genome name, but EarlGrey run as NU_Parvatl
declare -A BED_OVERRIDE=(
  [AS_AspfumUD1]="${TE_DIR}/AS_AspfumUD_EarlGrey/AS_AspfumUD_summaryFiles/AS_AspfumUD.filteredRepeats.bed"
  [NU_Parvatl_c]="${TE_DIR}/NU_Parvatl_EarlGrey/NU_Parvatl_summaryFiles/NU_Parvatl.filteredRepeats.bed"
)

if [[ -v "BED_OVERRIDE[${SP}]" ]]; then
    BED="${BED_OVERRIDE[${SP}]}"
else
    BED="${TE_DIR}/${SP}_EarlGrey/${SP}_summaryFiles/${SP}.filteredRepeats.bed"
fi

# Prefer real experimental nuc-seq WIG over NuclPred prediction when available
# Known species with real DANPOS Nuc-seq WIGs:
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

if [[ -v "REAL_WIGS[${SP}]" ]] && [[ -s "${REAL_WIGS[${SP}]}" ]]; then
    WIG="${REAL_WIGS[${SP}]}"
    echo "[INFO] Using real nuc-seq WIG for ${SP}: ${WIG}"
elif [[ -s "${NUCLPRED_WIG}" ]]; then
    WIG="${NUCLPRED_WIG}"
    echo "[INFO] Using NuclPred WIG for ${SP}"
else
    echo "[ERROR] No WIG found for ${SP}"
    exit 1
fi

# Sanity checks
if [[ ! -s "${BED}" ]]; then
    echo "[ERROR] Missing TE BED: ${BED}"
    exit 1
fi

echo "[$(date)] START ${SP}"

# ── F3: TE nucleosome bias ─────────────────────────────────────
F3_OUT="${OUT_BASE}/f3_bias/${SP}_f3_bias.tsv"
if [[ ! -s "${F3_OUT}" ]]; then
    echo "[$(date)] F3: nucleosome bias..."
    python3 "${SCRIPTS_DIR}/F3_te_nucleosome_bias.py" \
        "${WIG}" "${BED}" "${SP}" "${F3_OUT}"
else
    echo "[$(date)] F3: already done (${F3_OUT}), skipping"
fi

# ── F5: NDR/NOR TE density enrichment ─────────────────────────
F5_OUT="${OUT_BASE}/f5_ndr/${SP}_f5_ndr.tsv"
if [[ ! -s "${F5_OUT}" ]]; then
    echo "[$(date)] F5: NDR/NOR enrichment..."
    python3 "${SCRIPTS_DIR}/F5_ndr_te_density.py" \
        "${WIG}" "${BED}" "${SP}" "${F5_OUT}" \
        --ndr_pct 25 --nor_pct 75 --min_te_bp 200
else
    echo "[$(date)] F5: already done, skipping"
fi

# ── F6: meta-TE occupancy profiles ────────────────────────────
F6_OUT_DIR="${OUT_BASE}/f6_meta"
F6_DONE="${F6_OUT_DIR}/${SP}_meta_te_profiles.tsv"
if [[ ! -s "${F6_DONE}" ]]; then
    echo "[$(date)] F6: meta-TE profiles..."
    python3 "${SCRIPTS_DIR}/F6_meta_te_profile.py" \
        "${WIG}" "${BED}" "${SP}" "${F6_OUT_DIR}" \
        --window 2000 --min_insertions 10 --max_te_len 5000
else
    echo "[$(date)] F6: already done, skipping"
fi

echo "[$(date)] DONE ${SP}"
