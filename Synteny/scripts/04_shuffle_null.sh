#!/usr/bin/env bash
# =============================================================================
# 04_shuffle_null.sh — Compute null distribution of spurious collinear pairs
#                      using 100 rounds of bedtools shuffle per species pair.
#
# For each species pair within a clade:
#   1. Shuffle genes of SpA and SpB independently (bedtools shuffle, preserving
#      scaffold boundaries via -chrom; seed = SEED_BASE + replicate_index)
#   2. Recompute collinear pairs on shuffled positions using 03_collinear_pairs.py
#   3. Record s_ij per replicate → compute mean s_ij over 100 rounds
#
# Output:
#   {SHUFFLE_DIR}/{CLADE}/{SpA}__{SpB}.shuffle_summary.tsv
#   Columns: replicate  c_shuffled
#
#   {SHUFFLE_DIR}/{CLADE}/{CLADE}_null_means.tsv
#   Columns: species_A  species_B  s_ij_mean  s_ij_sd
#
# Requirements: bedtools >= 2.24, GNU ${PARALLEL}, Python 3 (for step 03)
#
# Usage:
#   bash 04_shuffle_null.sh --clade CLADE [--n_shuffles 100] [--threads 8]
#   bash 04_shuffle_null.sh --all [--threads 16]
#
# Author: EDF-TE project
# Last updated: 2026-02-20
# =============================================================================

set -euo pipefail

# Load shared configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ---------------------------------------------------------------------------
# Argument defaults
# ---------------------------------------------------------------------------
CLADES_TO_RUN=()
RUN_ALL=false
N_SHUFFLES=${N_SHUFFLES:-100}
THREADS=8
SEED=${SEED_BASE:-42}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --clade)   shift; CLADES_TO_RUN+=("$1"); shift ;;
        --all)     RUN_ALL=true; shift ;;
        --n_shuffles) shift; N_SHUFFLES="$1"; shift ;;
        --threads) shift; THREADS="$1"; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if $RUN_ALL; then
    CLADES_TO_RUN=("${ALL_CLADES[@]}")
fi

if [[ ${#CLADES_TO_RUN[@]} -eq 0 ]]; then
    echo "ERROR: specify --clade CLADE or --all"
    exit 1
fi

echo "[Step 04] Shuffle null distribution — ${N_SHUFFLES} replicates, ${THREADS} threads"
echo "          Clades: ${CLADES_TO_RUN[*]}"

# ---------------------------------------------------------------------------
# Helper: shuffle one species BED and write a temp BED for replicate i
# ---------------------------------------------------------------------------
shuffle_bed() {
    local sp="$1"
    local rep="$2"
    local out_bed="$3"
    local rep_seed=$(( SEED + rep ))

    ${BEDTOOLS} shuffle \
        -i  "${BED_DIR}/${sp}.bed" \
        -g  "${BED_DIR}/${sp}.genome" \
        -chrom \
        -noOverlapping \
        -seed "${rep_seed}" \
        > "${out_bed}"
}
export -f shuffle_bed
export SEED BED_DIR BEDTOOLS

# ---------------------------------------------------------------------------
# Helper: run step 03 collinear pairs on shuffled BED dir → return c_ij value
# ---------------------------------------------------------------------------
compute_shuffled_collinear() {
    local sp_a="$1"
    local sp_b="$2"
    local rep="$3"
    local clade_code="$4"
    local tmp_bed_dir="$5"
    local pairs_dir="$6"
    local tmp_col_dir="$7"

    # Run step 03 on the shuffled BED directory
    ${PYTHON} "${SCRIPT_DIR}/03_collinear_pairs.py" \
        --bed_dir   "${tmp_bed_dir}" \
        --pairs_dir "${pairs_dir}" \
        --out_dir   "${tmp_col_dir}" \
        --clade     "${clade_code}" 2>/dev/null

    # Extract c_ij from the collinear summary
    local summary="${tmp_col_dir}/${clade_code}/${clade_code}_collinear_summary.tsv"
    if [[ -f "${summary}" ]]; then
        awk -F'\t' -v a="${sp_a}" -v b="${sp_b}" \
            '$1==a && $2==b {print $4}' "${summary}"
    else
        echo "0"
    fi
}
export -f compute_shuffled_collinear

# ---------------------------------------------------------------------------
# Process each clade
# ---------------------------------------------------------------------------
for CLADE in "${CLADES_TO_RUN[@]}"; do
    echo ""
    echo "=== Clade: ${CLADE} ==="

    # Get species list for this clade
    SPECIES_STR="${CLADES[$CLADE]:-}"
    if [[ -z "${SPECIES_STR}" ]]; then
        echo "  Unknown clade ${CLADE} — skipping"
        continue
    fi
    read -ra SPECIES_LIST <<< "${SPECIES_STR}"

    if [[ ${#SPECIES_LIST[@]} -lt 2 ]]; then
        echo "  Only ${#SPECIES_LIST[@]} species — skipping"
        continue
    fi

    CLADE_SHUFFLE_DIR="${SHUFFLE_DIR}/${CLADE}"
    mkdir -p "${CLADE_SHUFFLE_DIR}"

    NULL_MEANS="${CLADE_SHUFFLE_DIR}/${CLADE}_null_means.tsv"
    echo -e "species_A\tspecies_B\ts_ij_mean\ts_ij_sd" > "${NULL_MEANS}"

    # Generate all pairs
    PAIRS=()
    for (( i=0; i<${#SPECIES_LIST[@]}-1; i++ )); do
        for (( j=i+1; j<${#SPECIES_LIST[@]}; j++ )); do
            PAIRS+=("${SPECIES_LIST[$i]}__%_${SPECIES_LIST[$j]}")
        done
    done

    echo "  ${#SPECIES_LIST[@]} species → ${#PAIRS[@]} pairs"

    for PAIR in "${PAIRS[@]}"; do
        SP_A="${PAIR%%__%_*}"
        SP_B="${PAIR##*__%_}"
        echo "  Processing ${SP_A} × ${SP_B} ..."

        PAIR_SHUFFLE_FILE="${CLADE_SHUFFLE_DIR}/${SP_A}__${SP_B}.shuffle_summary.tsv"
        echo -e "replicate\tc_shuffled" > "${PAIR_SHUFFLE_FILE}"

        # Temporary directory for this pair's shuffled BEDs
        TMP_DIR=$(mktemp -d)
        trap "rm -rf ${TMP_DIR}" EXIT

        # Run replicates in parallel using GNU parallel
        run_replicate() {
            local rep="$1"
            local sp_a="${SP_A}"
            local sp_b="${SP_B}"
            local clade="${CLADE}"
            local tmp="${TMP_DIR}/rep_${rep}"
            mkdir -p "${tmp}/beds" "${tmp}/pairs" "${tmp}/collinear"

            # Copy original pairs dir structure for this clade
            cp -r "${ORTHOLOG_PAIRS_DIR}/${clade}" "${tmp}/pairs/${clade}" 2>/dev/null || true

            # Shuffle BEDs for both species
            shuffle_bed "${sp_a}" "${rep}" "${tmp}/beds/${sp_a}.bed"
            shuffle_bed "${sp_b}" "${rep}" "${tmp}/beds/${sp_b}.bed"

            # Copy genome files (unchanged)
            cp "${BED_DIR}/${sp_a}.genome" "${tmp}/beds/" 2>/dev/null || true
            cp "${BED_DIR}/${sp_b}.genome" "${tmp}/beds/" 2>/dev/null || true

            # Run collinear computation on shuffled BEDs
            c_shuffled=$(compute_shuffled_collinear \
                "${sp_a}" "${sp_b}" "${rep}" "${clade}" \
                "${tmp}/beds" "${tmp}/pairs" "${tmp}/collinear")

            echo -e "${rep}\t${c_shuffled:-0}"
            rm -rf "${tmp}"
        }
        export -f run_replicate
        export SP_A SP_B CLADE TMP_DIR ORTHOLOG_PAIRS_DIR BED_DIR SCRIPT_DIR

        seq 1 "${N_SHUFFLES}" | \
            ${PARALLEL} --jobs "${THREADS}" --keep-order run_replicate {} \
            >> "${PAIR_SHUFFLE_FILE}"

        # Compute mean and SD of shuffled c values using awk
        awk -F'\t' 'NR>1 {
            sum += $2; sumsq += $2^2; n++
        } END {
            mean = (n>0) ? sum/n : 0
            sd   = (n>1) ? sqrt((sumsq - n*mean^2)/(n-1)) : 0
            print mean "\t" sd
        }' "${PAIR_SHUFFLE_FILE}" | \
        awk -v a="${SP_A}" -v b="${SP_B}" '{print a "\t" b "\t" $1 "\t" $2}' \
            >> "${NULL_MEANS}"

        echo "    → ${PAIR_SHUFFLE_FILE}"
    done

    echo "  Null means written to ${NULL_MEANS}"
done

echo ""
echo "[Step 04] Complete."
