#!/usr/bin/env bash
# =============================================================================
# 04_synteny_analysis.sh  —  Whole-genome synteny analysis pipeline
#
# Computes gene synteny ratios between all within-clade species pairs using
# OrthoFinder single-copy orthologs and a null model from 100 permutations.
#
# Steps:
#   1. Convert GFF to BED (gene order per scaffold)
#   2. Extract collinear ortholog pairs (OrthoFinder × gene order)
#   3. Count adjacent pair slots (N_ij) and observed collinear pairs (c_ij)
#   4. Generate shuffle-based null model (100 permutations, s_ij = mean)
#   5. Compute normalized synteny ratio r_ij per clade
#   6. Plot synteny ratio per clade (boxplot + heatmap)
#   7. Circos plots for representative species pairs
#   8. Patristic distance correction (partial Mantel, ape ultrametric tree)
#
# Dependencies: Python ≥3.10, bedtools, ape (R), OrthoFinder output
#
# Usage:
#   source ../config.sh
#   bash 04_synteny_analysis.sh [--clades "MU MO GL"] [--dry-run]
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

WORK_DIR="${S1_DIR}/synteny"
LOGS="${WORK_DIR}/logs"
ORTHO_DIR="${RESULTS_DIR}/orthofinder"

mkdir -p "${WORK_DIR}"/{01_bed,02_orthologs,03_collinear,04_shuffle,05_ratios,figures} "${LOGS}"

DRY_RUN=false
CLADES=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --clades)  CLADES="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

log() { echo "[$(date '+%H:%M:%S')] $*"; }
run() { local label="$1"; shift; log "  ${label}"; [[ "${DRY_RUN}" == true ]] && { echo "    DRY: $*"; return; }; "$@"; }

# =============================================================================
# STEP 1 — GFF → BED (gene order per scaffold)
# =============================================================================
log "=== STEP 1: GFF to BED ==="
for gff in "${GFF_DIR}"/*.gff3 "${GFF_DIR}"/*.gff; do
    [[ -f "${gff}" ]] || continue
    sp=$(basename "${gff}" | sed 's/\.\(gff3\|gff\)$//')
    out="${WORK_DIR}/01_bed/${sp}.bed"
    [[ -f "${out}" ]] && continue
    run "gff2bed ${sp}" python3 "${SCRIPT_DIR}/helpers/gff_to_bed.py" \
        --gff "${gff}" --out "${out}"
done

# =============================================================================
# STEP 2 — Extract single-copy orthologs per species pair
# =============================================================================
log "=== STEP 2: Ortholog extraction ==="
run "extract_orthologs" python3 "${SCRIPT_DIR}/helpers/extract_orthologs.py" \
    --ortho_dir "${ORTHO_DIR}" \
    --bed_dir   "${WORK_DIR}/01_bed" \
    --out_dir   "${WORK_DIR}/02_orthologs"

# =============================================================================
# STEP 3 — Collinear pair counting
# =============================================================================
log "=== STEP 3: Collinear pairs ==="
run "collinear_pairs" python3 "${SCRIPT_DIR}/helpers/collinear_pairs.py" \
    --ortho_dir  "${WORK_DIR}/02_orthologs" \
    --out_dir    "${WORK_DIR}/03_collinear"

# =============================================================================
# STEP 4 — Shuffle null model (SLURM array, 100 permutations)
# =============================================================================
log "=== STEP 4: Shuffle null model ==="
if [[ "${DRY_RUN}" == false ]]; then
    sbatch --partition="${SLURM_PARTITION}" \
           --cpus-per-task=4 --mem="8G" \
           --array="1-100" \
           --output="${LOGS}/shuffle_%a.log" \
           --wrap="python3 ${SCRIPT_DIR}/helpers/shuffle_null.py \
               --collinear_dir ${WORK_DIR}/03_collinear \
               --out_dir       ${WORK_DIR}/04_shuffle \
               --seed \${SLURM_ARRAY_TASK_ID}"
fi

# =============================================================================
# STEP 5 — Synteny ratios
# =============================================================================
log "=== STEP 5: Synteny ratios ==="
run "synteny_ratio" python3 "${SCRIPT_DIR}/helpers/synteny_ratio.py" \
    --collinear_dir "${WORK_DIR}/03_collinear" \
    --shuffle_dir   "${WORK_DIR}/04_shuffle" \
    --out_dir       "${WORK_DIR}/05_ratios" \
    --all

# =============================================================================
# STEP 6 — Synteny ratio plots
# =============================================================================
log "=== STEP 6: Synteny plots ==="
run "plot_synteny" python3 "${SCRIPT_DIR}/05_circos_visualization.py" \
    --mode       summary \
    --ratio_dir  "${WORK_DIR}/05_ratios" \
    --out_dir    "${WORK_DIR}/figures"

log "Synteny analysis complete. Outputs in ${WORK_DIR}"
