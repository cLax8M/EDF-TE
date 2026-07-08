#!/usr/bin/env bash
# =============================================================================
# 00_wgd_pipeline.sh  —  WGD pre-processing pipeline (dmd → MCL → KS → syn)
#
# Generates per-species paralog Ks distributions (ks.tsv) and synteny-based
# multiplicon evidence, which feed into 01_wgd_detection.py and downstream
# TE–WGD analyses.
#
# Steps (review each helper before running):
#   1. MCL / OrthoFinder duplications  → paralog pairs per species
#   2. Extract CDS for paralog pairs
#   3. Filter high-confidence pairs
#   4. PAL2NAL + yn00 (PAML) → raw Ks values
#   5. Merge per-species Ks → ks.tsv
#   6. WGD dmd (Ks peak detection)
#   7. MCScanX synteny → collinear blocks / multiplicones
#   8. Integrate synteny evidence with WGD calls
#
# Dependencies: Python ≥3.10, PAML (yn00), PAL2NAL, MCScanX, OrthoFinder output
#
# Usage:
#   source ../config.sh                              # project-wide config
#   source helpers/config_wgd.sh                     # WGD-specific config
#   bash helpers/00_wgd_pipeline.sh [--dry-run]
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"
source "${SCRIPT_DIR}/config_wgd.sh"

DRY_RUN=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

log() { echo "[$(date '+%H:%M:%S')] $*"; }
run() { local label="$1"; shift; log "  ${label}"; [[ "${DRY_RUN}" == true ]] && { echo "    DRY: $*"; return; }; "$@"; }

# ── Step 1: MCL-like paralog extraction from OrthoFinder ──────────────────────
log "=== STEP 1: OrthoFinder duplications → paralog pairs ==="
run "parse_duplications" python3 "${SCRIPT_DIR}/00_parse_orthofinder_duplications.py" \
    --orthofinder-dir "${ORTHOFINDER_DIR}" \
    --out-dir       "${WGD_WORK_DIR}/01_paralog_pairs"

# ── Step 2: Prepare KS pairs ──────────────────────────────────────────────────
log "=== STEP 2: Prepare KS gene pairs per species ==="
run "prepare_ks_pairs" python3 "${SCRIPT_DIR}/00_prepare_ks_pairs.py" \
    --dups_dir    "${WGD_WORK_DIR}/01_paralog_pairs" \
    --out_dir     "${WGD_WORK_DIR}/02_ks_pairs"

# ── Step 3: Extract CDS for paralog pairs ─────────────────────────────────────
log "=== STEP 3: Extract CDS sequences ==="
run "extract_cds" bash "${SCRIPT_DIR}/00_extract_cds_for_ks.sh" \
    "${WGD_WORK_DIR}/02_ks_pairs" \
    "${CDS_DIR}" \
    "${WGD_WORK_DIR}/03_cds"

# ── Step 4: Filter high-confidence pairs ──────────────────────────────────────
log "=== STEP 4: Filter KS pairs ==="
run "filter_pairs" python3 "${SCRIPT_DIR}/00_filter_pairs_for_ks.py" \
    --pairs_dir  "${WGD_WORK_DIR}/02_ks_pairs" \
    --cds_dir    "${WGD_WORK_DIR}/03_cds" \
    --out_dir    "${WGD_WORK_DIR}/04_filtered_pairs"

# ── Step 5: PAL2NAL + yn00 → Ks ───────────────────────────────────────────────
log "=== STEP 5: PAL2NAL + yn00 (SLURM array) ==="
log "  Submit with: sbatch ${SCRIPT_DIR}/slurm_run_paml_ks_array.sh"
log "  Or run locally:  bash ${SCRIPT_DIR}/00_compute_ks_paml.sh"

# ── Step 6: Merge KS results ──────────────────────────────────────────────────
log "=== STEP 6: Merge KS results → ks.tsv per species ==="
run "merge_ks" python3 "${SCRIPT_DIR}/00_merge_ks_results.py" \
    --paml_dir    "${WGD_WORK_DIR}/05_paml_results" \
    --out_dir     "${WGD_WORK_DIR}/06_ks_tsv"

# ── Step 7: Detect WGD peaks (dmd) ────────────────────────────────────────────
log "=== STEP 7: WGD peak detection (dmd) ==="
run "detect_peaks" python3 "${SCRIPT_DIR}/00_detect_wgd_peaks.py" \
    --ks_dir      "${WGD_WORK_DIR}/06_ks_tsv" \
    --out_dir     "${WGD_WORK_DIR}/07_wgd_peaks"

# ── Step 8: MCScanX synteny → multiplicones ───────────────────────────────────
log "=== STEP 8: Synteny (MCScanX) → multiplicons ==="
run "prepare_mcscanx" bash "${SCRIPT_DIR}/00_prepare_mcscanx_input.sh"
run "run_mcscanx"      bash "${SCRIPT_DIR}/00_run_mcscanx.sh"
run "integrate_syn"    python3 "${SCRIPT_DIR}/00_integrate_synteny_wgd.py" \
    --mcscanx_dir  "${WGD_WORK_DIR}/08_mcscanx" \
    --wgd_dir      "${WGD_WORK_DIR}/07_wgd_peaks" \
    --gff_dir      "${GFF_DIR}" \
    --out_dir      "${WGD_WORK_DIR}/09_multiplicons"

log "=== WGD pre-processing complete ==="
log "Outputs:"
log "  KS per species:    ${WGD_WORK_DIR}/06_ks_tsv/"
log "  WGD peaks:         ${WGD_WORK_DIR}/07_wgd_peaks/"
log "  Multiplicons:      ${WGD_WORK_DIR}/09_multiplicons/"
log ""
log "Next: run 01_wgd_detection.py with --ks_dir ${WGD_WORK_DIR}/06_ks_tsv"
