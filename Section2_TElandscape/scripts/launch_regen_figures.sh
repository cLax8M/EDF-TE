#!/usr/bin/env bash
# launch_regen_figures.sh
# Regenerate all tree-dependent and TE-data-dependent figures after:
#   1) Final species tree finalized (final_LG_G_rotated.nwk, May 1 2026)
#   2) TE annotations completed for HO_Homosap and HO_Magisepte_c
#
# Usage:  bash launch_regen_figures.sh
# =============================================================================
set -euo pipefail

MAMBA=/home/carlos/.local/bin/micromamba
ENV_PY=syntenyEDFTE
ENV_R=esptreeEDFTE
SCRIPTS=/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts
SYN_SCRIPTS=/data/carlos/EDF-TE/EDF-TE/Synteny/scripts
LOG=/data/carlos/EDF-TE/logs/regen_figures
mkdir -p "$LOG"

EARLGREY_DIR=/data/carlos/EDF-TE/TE_analysis/annotated
GF_SUMMARY=/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures_V2/ALL_SPECIES_ALL_SUMMARY.txt
D01_OUT=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/D01_te_genomics_merged_provisional.tsv
TREE=/legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/final_LG_G_rotated.nwk
RATIOS=/data/carlos/EDF-TE/EDF-TE/Synteny/data/ratios/all_clades_synteny_ratios.tsv
SYN_CORRECTED=/legserv/Temp/Carlos/EDF-TE/results/Synteny_guide_corrected

BED_DIR=/legserv/Temp/Carlos/EDF-TE/results/Synteny/bed
COLLINEAR_DIR=/data/carlos/EDF-TE/EDF-TE/Synteny/data/collinear
CIRCOS_OUT=/legserv/Temp/Carlos/EDF-TE/results/Synteny/circos_figures
mkdir -p "$CIRCOS_OUT"

# ─── Job 1: Regenerate D01 (TE data for all species incl. Homosap+Magisepte) ─
echo "[SUBMIT] D01 regeneration"
JID1=$(sbatch \
  --job-name="regen_D01" \
  --partition=normal.168h \
  --cpus-per-task=2 \
  --mem=8G \
  --output="${LOG}/regen_D01_%j.out" \
  --error="${LOG}/regen_D01_%j.err" \
  --parsable \
  --wrap="${MAMBA} run -n ${ENV_PY} python3 ${SCRIPTS}/D01_merge_te_genomics.py \
    ${EARLGREY_DIR} \
    ${GF_SUMMARY} \
    ${D01_OUT}")
echo "  D01 job: $JID1"

# ─── Job 2: Update figure_species_data.tsv (needs D01) ───────────────────────
echo "[SUBMIT] prepare_figure_data (after D01)"
JID2=$(sbatch \
  --job-name="regen_figdata" \
  --partition=normal.168h \
  --cpus-per-task=2 \
  --mem=8G \
  --dependency="afterok:${JID1}" \
  --output="${LOG}/regen_figdata_%j.out" \
  --error="${LOG}/regen_figdata_%j.err" \
  --parsable \
  --wrap="${MAMBA} run -n ${ENV_PY} python3 ${SCRIPTS}/prepare_figure_data.py")
echo "  figure_data job: $JID2"

# ─── Job 3: Circular (radial) tree (needs updated figure_species_data) ───────
echo "[SUBMIT] Circular tree (after figure_data)"
JID3=$(sbatch \
  --job-name="regen_circ_tree" \
  --partition=normal.168h \
  --cpus-per-task=4 \
  --mem=16G \
  --dependency="afterok:${JID2}" \
  --output="${LOG}/regen_circ_tree_%j.out" \
  --error="${LOG}/regen_circ_tree_%j.err" \
  --parsable \
  --wrap="${MAMBA} run -n ${ENV_R} Rscript ${SCRIPTS}/plot_circular_tree.R")
echo "  Circular tree job: $JID3"

# ─── Job 4: Rectangular cladogram with TE bars (needs updated figure_data) ───
echo "[SUBMIT] TE cladogram (after figure_data)"
JID4=$(sbatch \
  --job-name="regen_cladogram" \
  --partition=normal.168h \
  --cpus-per-task=2 \
  --mem=8G \
  --dependency="afterok:${JID2}" \
  --output="${LOG}/regen_cladogram_%j.out" \
  --error="${LOG}/regen_cladogram_%j.err" \
  --parsable \
  --wrap="${MAMBA} run -n ${ENV_PY} python3 ${SCRIPTS}/plot_cladogram_te_v2.py")
echo "  Cladogram job: $JID4"

# ─── Job 5: Patristic distance correction for synteny (independent) ──────────
echo "[SUBMIT] Patristic correction (independent)"
JID5=$(sbatch \
  --job-name="regen_patristic" \
  --partition=normal.168h \
  --cpus-per-task=4 \
  --mem=12G \
  --output="${LOG}/regen_patristic_%j.out" \
  --error="${LOG}/regen_patristic_%j.err" \
  --parsable \
  --wrap="${MAMBA} run -n ${ENV_PY} python3 ${SYN_SCRIPTS}/09_patristic_correction.py \
    --ratios ${RATIOS} \
    --tree   ${TREE} \
    --out_dir ${SYN_CORRECTED}")
echo "  Patristic job: $JID5"

# ─── Job 6: Circos (4 representative pairs, needs patristic correction) ──────
# Pairs: MI_Enche1×MI_Encro1, MO_Moramb1×MO_Morwol1, MU_Abscae1×MU_Absrep1, EN_Conth1×EN_Entmai1
echo "[SUBMIT] Circos (after patristic)"
JID6=$(sbatch \
  --job-name="regen_circos" \
  --partition=normal.168h \
  --cpus-per-task=4 \
  --mem=16G \
  --dependency="afterok:${JID5}" \
  --output="${LOG}/regen_circos_%j.out" \
  --error="${LOG}/regen_circos_%j.err" \
  --parsable \
  --wrap="${MAMBA} run -n ${ENV_PY} python3 ${SYN_SCRIPTS}/08_plot_circos.py \
    --batch \
    --bed_dir     ${BED_DIR} \
    --collinear_dir ${COLLINEAR_DIR} \
    --out_dir     ${CIRCOS_OUT} \
    --pairs \
      MI:MI_Enche1:MI_Encro1 \
      MO:MO_Moramb1:MO_Morwol1 \
      MU:MU_Abscae1:MU_Absrep1 \
      EN:EN_Conth1:EN_Entmai1 \
    --simple_style --top_by_size 30")
echo "  Circos job: $JID6"

# ─── Job 7: Synteny boxplot by clade (needs patristic correction) ─────────────
echo "[SUBMIT] Synteny summary figures (after patristic)"
JID7=$(sbatch \
  --job-name="regen_syn_box" \
  --partition=normal.168h \
  --cpus-per-task=2 \
  --mem=8G \
  --dependency="afterok:${JID5}" \
  --output="${LOG}/regen_syn_box_%j.out" \
  --error="${LOG}/regen_syn_box_%j.err" \
  --parsable \
  --wrap="${MAMBA} run -n ${ENV_PY} python3 ${SYN_SCRIPTS}/plot_synteny_summary.py")
echo "  Synteny boxplot job: $JID7"

echo ""
echo "=== Submission summary ==="
echo "  J1 D01 regen:          $JID1"
echo "  J2 figure_data:        $JID2  (after J1)"
echo "  J3 circular tree:      $JID3  (after J2)"
echo "  J4 TE cladogram:       $JID4  (after J2)"
echo "  J5 patristic:          $JID5  (independent)"
echo "  J6 circos 4 pairs:     $JID6  (after J5)"
echo "  J7 synteny boxplot:    $JID7  (after J5)"
echo ""
echo "Logs: $LOG"
echo "Outputs:"
echo "  D01:          $D01_OUT"
echo "  Circular tree: /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures/circular_tree_v4.pdf"
echo "  Cladogram:    /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures/cladogram_te_types_shannon_ri.pdf"
echo "  Patristic:    ${SYN_CORRECTED}/corrected_synteny_ratios.tsv"
echo "  Circos:       ${CIRCOS_OUT}/"
echo "  Syn boxplot:  /legserv/Temp/Carlos/EDF-TE/results/Synteny/figures_summary/"
