#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — End-to-end orchestrator for the EDF-TE synteny analysis.
#
# Runs steps 01–06 (and optionally 07) for one or all clades.
# Each step is logged to {OUT_DIR}/logs/{CLADE}_step{N}.log
#
# Usage:
#   bash run_pipeline.sh --clade HO [--threads 8] [--steps 1-6]
#   bash run_pipeline.sh --all [--threads 16] [--skip_shuffle]
#
# Options:
#   --clade CODE     Process a single clade (can repeat)
#   --all            Process all clades
#   --threads N      Number of parallel threads for shuffle step (default: 8)
#   --steps RANGE    Comma-separated steps or range, e.g. 1,3-6 (default: 1-6)
#   --skip_shuffle   Skip step 04 (shuffle null; use if already computed)
#   --circos A B     After step 06, run step 07 for species pair A × B
#   --dry_run        Print commands without executing
#
# Author: EDF-TE project
# Last updated: 2026-02-20
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
CLADES_TO_RUN=()
RUN_ALL=false
THREADS=8
STEPS="1,2,3,4,5,6"
SKIP_SHUFFLE=false
CIRCOS_A=""
CIRCOS_B=""
DRY_RUN=false
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --clade)         shift; CLADES_TO_RUN+=("$1"); shift ;;
        --all)           RUN_ALL=true; shift ;;
        --threads)       shift; THREADS="$1"; shift ;;
        --steps)         shift; STEPS="$1"; shift ;;
        --skip_shuffle)  SKIP_SHUFFLE=true; shift ;;
        --circos)        shift; CIRCOS_A="$1"; shift; CIRCOS_B="$1"; shift ;;
        --dry_run)       DRY_RUN=true; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if $RUN_ALL; then
    CLADES_TO_RUN=("${ALL_CLADES[@]}")
fi

if [[ ${#CLADES_TO_RUN[@]} -eq 0 ]]; then
    echo "ERROR: specify --clade CLADE or --all"; exit 1
fi

# Expand STEPS string into a set
steps_active() {
    local step="$1"
    echo "${STEPS}" | tr ',' '\n' | while read -r s; do
        if [[ "${s}" == *"-"* ]]; then
            lo="${s%-*}"; hi="${s#*-}"
            for i in $(seq "${lo}" "${hi}"); do echo "${i}"; done
        else
            echo "${s}"
        fi
    done | grep -q "^${step}$"
}

run_cmd() {
    local label="$1"; shift
    if $DRY_RUN; then
        echo "[DRY-RUN] ${label}: $*"
    else
        echo "[$(date '+%H:%M:%S')] ${label}"
        "$@"
    fi
}

# ---------------------------------------------------------------------------
# Build clade argument string
# ---------------------------------------------------------------------------
CLADE_ARGS=""
for c in "${CLADES_TO_RUN[@]}"; do
    CLADE_ARGS="${CLADE_ARGS} --clade ${c}"
done
# shellcheck disable=SC2086
CLADE_ARGS_ARRAY=(${CLADE_ARGS})

echo "========================================"
echo " EDF-TE Synteny Pipeline"
echo " Clades : ${CLADES_TO_RUN[*]}"
echo " Steps  : ${STEPS}"
echo " Threads: ${THREADS}"
echo "========================================"

# ---------------------------------------------------------------------------
# Step 01 — GFF → BED
# ---------------------------------------------------------------------------
if steps_active 1; then
    run_cmd "Step 01: GFF → BED" \
        ${PYTHON} "${SCRIPT_DIR}/01_gff_to_bed.py" \
            --gff_dir  "${GFF_DIR}" \
            --faa_dir  "${FAA_DIR}" \
            --out_dir  "${BED_DIR}" \
        2>&1 | tee "${LOG_DIR}/step01.log"
fi

# ---------------------------------------------------------------------------
# Step 02 — Extract 1:1 orthologs
# ---------------------------------------------------------------------------
if steps_active 2; then
    run_cmd "Step 02: Extract ortholog pairs" \
        ${PYTHON} "${SCRIPT_DIR}/02_extract_orthologs.py" \
            --orthologues_dir "${ORTHOLOGUES_DIR}" \
            --out_dir         "${ORTHOLOG_PAIRS_DIR}" \
            "${CLADE_ARGS_ARRAY[@]}" \
        2>&1 | tee "${LOG_DIR}/step02.log"
fi

# ---------------------------------------------------------------------------
# Step 03 — Collinear pairs
# ---------------------------------------------------------------------------
if steps_active 3; then
    run_cmd "Step 03: Collinear pairs" \
        ${PYTHON} "${SCRIPT_DIR}/03_collinear_pairs.py" \
            --bed_dir   "${BED_DIR}" \
            --pairs_dir "${ORTHOLOG_PAIRS_DIR}" \
            --out_dir   "${COLLINEAR_DIR}" \
            "${CLADE_ARGS_ARRAY[@]}" \
        2>&1 | tee "${LOG_DIR}/step03.log"
fi

# ---------------------------------------------------------------------------
# Step 04 — Shuffle null (bedtools + GNU parallel)
# ---------------------------------------------------------------------------
if steps_active 4 && ! $SKIP_SHUFFLE; then
    run_cmd "Step 04: Shuffle null distribution" \
        bash "${SCRIPT_DIR}/04_shuffle_null.sh" \
            "${CLADE_ARGS_ARRAY[@]}" \
            --n_shuffles "${N_SHUFFLES}" \
            --threads    "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/step04.log"
elif $SKIP_SHUFFLE; then
    echo "[Step 04] SKIPPED (--skip_shuffle)"
fi

# ---------------------------------------------------------------------------
# Step 05 — Synteny ratios
# ---------------------------------------------------------------------------
if steps_active 5; then
    run_cmd "Step 05: Compute synteny ratios" \
        ${PYTHON} "${SCRIPT_DIR}/05_synteny_ratio.py" \
            --collinear_dir "${COLLINEAR_DIR}" \
            --shuffle_dir   "${SHUFFLE_DIR}" \
            --out_dir       "${RATIO_DIR}" \
            "${CLADE_ARGS_ARRAY[@]}" \
        2>&1 | tee "${LOG_DIR}/step05.log"
fi

# ---------------------------------------------------------------------------
# Step 06 — Heatmaps
# ---------------------------------------------------------------------------
if steps_active 6; then
    run_cmd "Step 06: Plot heatmaps" \
        ${PYTHON} "${SCRIPT_DIR}/06_plot_synteny.py" \
            --ratio_dir "${RATIO_DIR}" \
            --out_dir   "${FIG_DIR}" \
            "${CLADE_ARGS_ARRAY[@]}" \
        2>&1 | tee "${LOG_DIR}/step06.log"
fi

# ---------------------------------------------------------------------------
# (Optional) Step 07 — Circos
# ---------------------------------------------------------------------------
if [[ -n "${CIRCOS_A}" && -n "${CIRCOS_B}" ]]; then
    # Auto-detect clade from species prefix
    CIRCOS_CLADE="${CIRCOS_A%%_*}"
    CIRCOS_OUT="${FIG_DIR}/circos_${CIRCOS_A}__${CIRCOS_B}"
    run_cmd "Step 07: Circos prep (${CIRCOS_A} × ${CIRCOS_B})" \
        ${PYTHON} "${SCRIPT_DIR}/07_circos_prep.py" \
            --sp_a          "${CIRCOS_A}" \
            --sp_b          "${CIRCOS_B}" \
            --clade         "${CIRCOS_CLADE}" \
            --bed_dir       "${BED_DIR}" \
            --collinear_dir "${COLLINEAR_DIR}" \
            --out_dir       "${CIRCOS_OUT}" \
        2>&1 | tee "${LOG_DIR}/step07_${CIRCOS_A}__${CIRCOS_B}.log"
fi

echo ""
echo "========================================"
echo " Pipeline complete."
echo " Figures : ${FIG_DIR}"
echo " Ratios  : ${RATIO_DIR}"
echo "========================================"
