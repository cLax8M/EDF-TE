#!/usr/bin/env bash
# ============================================================
# run_htt_pipeline.sh  — EDF-TE project wrapper
# 
# Submits HTT pipeline steps 02-09 sequentially using SLURM
# dependencies. Paths configured for the EDF-TE project.
# 
# Usage:
#   bash run_htt_pipeline.sh [--start STEP] [--end STEP] [--dry-run]
#   e.g.  bash run_htt_pipeline.sh --start 3 --end 6
#
# Requires: steps 01 already done at HTT_RESULTS/01_phylo_tree/
# ============================================================

set -euo pipefail

# ─── Project paths ────────────────────────────────────────────────────────────
HTT_SCRIPTS=/legserv/Temp/Carlos/EDF-TE/data/Extra/HTT_hpc_scripts
HTT_RESULTS=/legserv/Temp/Carlos/EDF-TE/results/HTT
GENOMES_DIR=/legserv/Temp/Carlos/EDF-TE/data/genomes
TE_ANNOT_BASE=/data/carlos/EDF-TE/TE_analysis/annotated
LOGS_DIR=/data/carlos/EDF-TE/logs/HTT
mkdir -p "$LOGS_DIR"

# ─── Argument parsing ─────────────────────────────────────────────────────────
START_STEP=2
END_STEP=9
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case $1 in
    --start)  START_STEP=$2; shift 2 ;;
    --end)    END_STEP=$2;   shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

echo "HTT pipeline: steps ${START_STEP}–${END_STEP}"
echo "DRY_RUN=${DRY_RUN}"

# ─── helper: patch conda → micromamba in a script ─────────────────────────────
patch_conda_to_micromamba () {
  local SCRIPT="$1"
  # Replace the shell hook initialisation
  sed -i 's|eval "\$(conda shell\.bash hook)"|eval "$(micromamba shell hook --shell bash)"|g' "$SCRIPT"
  # Replace activate / deactivate
  sed -i 's/\bconda activate\b/micromamba activate/g'   "$SCRIPT"
  sed -i 's/\bconda deactivate\b/micromamba deactivate/g' "$SCRIPT"
}

# ─── Step 02: divergence times (PhyloRank + RED values) ───────────────────────
submit_step02 () {
  local TREE="${HTT_RESULTS}/01_phylo_tree/phylogeny_rooted.treefile"
  local RESULTS="${HTT_RESULTS}/02_div_times"
  mkdir -p "$RESULTS"

  # Create a configured version of the script
  local SCRIPT="${LOGS_DIR}/02_div_times_configured.sh"
  sed "s|ROOT=/path/to/HTT|ROOT_UNUSED=1\nTREE_RESULTS=${HTT_RESULTS}/01_phylo_tree\nRESULTS=${RESULTS}\nLOGS=${LOGS_DIR}|" \
    "${HTT_SCRIPTS}/02_div_times.sh" > "$SCRIPT"
  # Patch TREE_FILE variable
  sed -i "s|TREE_FILE=.*|TREE_FILE=${TREE}|" "$SCRIPT"
  # Replace conda with micromamba
  patch_conda_to_micromamba "$SCRIPT"
  chmod +x "$SCRIPT"

  if $DRY_RUN; then
    echo "[DRY] sbatch --partition=normal.168h --cpus-per-task=4 --mem=16G --time=24:00:00 $SCRIPT"
  else
    sbatch --partition=normal.168h --cpus-per-task=4 --mem=16G --time=24:00:00 \
           --job-name=HTT_02_divtimes --output="${LOGS_DIR}/02_%j.out" \
           "$SCRIPT"
  fi
}

# ─── Step 03: Ks divergence ───────────────────────────────────────────────────
submit_step03 () {
  local RESULTS="${HTT_RESULTS}/03_ks_divergence"
  mkdir -p "$RESULTS"

  local SCRIPT="${LOGS_DIR}/03_ks_configured.sh"
  sed "s|ROOT=/path/to/HTT|ROOT_UNUSED=1|" "${HTT_SCRIPTS}/03_ks_divergence.sh" > "$SCRIPT"
  sed -i "s|TREE_RESULTS=\${ROOT}/results/01_results_phylo_tree|TREE_RESULTS=${HTT_RESULTS}/01_phylo_tree|" "$SCRIPT"
  sed -i "s|RESULTS=\${ROOT}/results/03_results_ks_divergence|RESULTS=${RESULTS}|" "$SCRIPT"
  sed -i "s|LOGS=\${ROOT}/logs/03_logs_ks_divergence|LOGS=${LOGS_DIR}|" "$SCRIPT"
  # Replace conda with micromamba
  patch_conda_to_micromamba "$SCRIPT"
  chmod +x "$SCRIPT"

  if $DRY_RUN; then
    echo "[DRY] sbatch --partition=normal.168h --cpus-per-task=8 --mem=32G $SCRIPT"
  else
    sbatch --partition=normal.168h --cpus-per-task=8 --mem=32G --time=48:00:00 \
           --job-name=HTT_03_ks --output="${LOGS_DIR}/03_%j.out" \
           "$SCRIPT"
  fi
}

# ─── Step 04: TE annotation (skip — we already have filteredRepeats) ──────────
setup_step04_symlinks () {
  # Instead of re-running EarlGrey, create symlinks from existing annotations
  local TE_RESULTS="${HTT_RESULTS}/04_te_annotation"
  mkdir -p "$TE_RESULTS"
  
  echo "Creating TE annotation symlinks for all species…"
  for sp_dir in "${TE_ANNOT_BASE}"/*_EarlGrey; do
    sp=$(basename "$sp_dir" | sed 's/_EarlGrey$//')
    bed="${sp_dir}/${sp}_summaryFiles/${sp}.filteredRepeats.bed"
    gff="${sp_dir}/${sp}_summaryFiles/${sp}.filteredRepeats.gff"
    if [[ -f "$bed" && -f "$gff" ]]; then
      ln -sf "$bed" "${TE_RESULTS}/${sp}.filteredRepeats.bed" 2>/dev/null || true
      ln -sf "$gff" "${TE_RESULTS}/${sp}.filteredRepeats.gff" 2>/dev/null || true
    fi
  done
  echo "  Symlinks created in ${TE_RESULTS}"
  ls "${TE_RESULTS}" | wc -l
}

# ─── Step 06: HTT candidates (BLAST-based) ───────────────────────────────────
submit_step06 () {
  local TE_RESULTS="${HTT_RESULTS}/04_te_annotation"
  local KS_RESULTS="${HTT_RESULTS}/03_ks_divergence"
  local RESULTS="${HTT_RESULTS}/06_htt_candidates"
  mkdir -p "$RESULTS"
  
  local SCRIPT="${LOGS_DIR}/06_htt_configured.sh"
  sed "s|ROOT=/path/to/HTT|ROOT_UNUSED=1|" "${HTT_SCRIPTS}/06_htt_candidates.sh" > "$SCRIPT"
  sed -i "s|\${ROOT}/results/04_results_te_annotation|${TE_RESULTS}|g" "$SCRIPT"
  sed -i "s|\${ROOT}/results/03_results_ks_divergence|${KS_RESULTS}|g" "$SCRIPT"
  sed -i "s|\${ROOT}/results/06_results_htt_candidates|${RESULTS}|g"   "$SCRIPT"
  sed -i "s|\${ROOT}/logs|${LOGS_DIR}|g" "$SCRIPT"
  chmod +x "$SCRIPT"
  
  if $DRY_RUN; then
    echo "[DRY] sbatch --partition=normal.168h --cpus-per-task=16 --mem=64G $SCRIPT"
  else
    sbatch --partition=normal.168h --cpus-per-task=16 --mem=64G --time=72:00:00 \
           --job-name=HTT_06_candidates --output="${LOGS_DIR}/06_%j.out" \
           "$SCRIPT"
  fi
}

# ─── DISPATCH ─────────────────────────────────────────────────────────────────
echo ""
echo "=== Pre-flight checks ==="
[[ -f "${HTT_RESULTS}/01_phylo_tree/phylogeny_rooted.treefile" ]] && \
  echo "  [OK] Step 01 tree exists" || echo "  [ERROR] Step 01 tree missing!"

echo ""
for step in $(seq "$START_STEP" "$END_STEP"); do
  echo "=== Submitting step ${step} ==="
  case $step in
    2) submit_step02 ;;
    3) submit_step03 ;;
    4) setup_step04_symlinks ;;
    5) echo "  [SKIP] Step 05 (TE classification) — depends on step 04 output; run after 04" ;;
    6) submit_step06 ;;
    7) echo "  [TODO] Step 07 (flankid) — configure after step 06 completes" ;;
    8) echo "  [TODO] Step 08 (clustering) — configure after step 07 completes" ;;
    9) echo "  [TODO] Step 09 (minimal HTT events) — configure after step 08 completes" ;;
    *) echo "  [SKIP] Unknown step $step" ;;
  esac
done

echo ""
echo "Done. Check jobs with: squeue -u $USER"
