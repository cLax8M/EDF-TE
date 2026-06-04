#!/bin/bash
# launch_cooption_figures.sh
# Run all 8 co-option manuscript figures.
# Some can run in parallel; Fig H step 2 (phylo) needs step 1 first.
#
# Usage: bash launch_cooption_figures.sh [figA|figB|figC|figD|figE|figF|figG|figH|all]

set -e

SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Cooption/scripts"
LOG_DIR="/data/carlos/EDF-TE/logs/cooption_figures"
mkdir -p "$LOG_DIR"

ENV_PYTHON="syntenyEDFTE"   # has python3 + bedtools + pandas + matplotlib + networkx
ENV_R="esptreeEDFTE"        # has R + ggtree + ggplot2 + ape + phytools
MAMBA=/home/carlos/.local/bin/micromamba

TARGET="${1:-all}"

run_python() {
  local name="$1"
  local script="$2"
  local mem="${3:-16G}"
  local cpu="${4:-4}"
  echo "[SUBMIT] $name"
  sbatch \
    --job-name="$name" \
    --partition=normal.168h \
    --cpus-per-task="${cpu}" \
    --mem="${mem}" \
    --output="${LOG_DIR}/${name}_%j.out" \
    --error="${LOG_DIR}/${name}_%j.err" \
    --wrap="${MAMBA} run -n ${ENV_PYTHON} python3 ${SCRIPTS_DIR}/${script}"
}

run_r() {
  local name="$1"
  local script="$2"
  local mem="${3:-8G}"
  local cpu="${4:-2}"
  echo "[SUBMIT] $name"
  sbatch \
    --job-name="$name" \
    --partition=normal.168h \
    --cpus-per-task="${cpu}" \
    --mem="${mem}" \
    --output="${LOG_DIR}/${name}_%j.out" \
    --error="${LOG_DIR}/${name}_%j.err" \
    --wrap="${MAMBA} run -n ${ENV_R} Rscript ${SCRIPTS_DIR}/${script}"
}

case "$TARGET" in
  figA|all)
    run_r "figA_tree" "fig_A_tree_cooption.R" "12G" "2"
    ;;&
  figB|all)
    run_r "figB_boxplot" "fig_B_clade_boxplot.R" "8G" "2"
    ;;&
  figC|all)
    run_python "figC_domains" "fig_C_domain_barchart.py" "16G" "4"
    ;;&
  figD|all)
    run_python "figD_multicopy" "fig_D_multicopy_te_figure.py" "8G" "2"
    ;;&
  figE|all)
    run_python "figE_networks" "fig_E_te_centric_networks.py" "24G" "4"
    ;;&
  figF|all)
    run_python "figF_duf" "fig_F_duf_networks.py" "16G" "4"
    ;;&
  figG|all)
    run_python "figG_introns" "fig_G_intron_analysis.py" "16G" "4"
    ;;&
  figH|all)
    # Step 1: find candidates + extract sequences
    JID=$(sbatch \
      --job-name="figH_prep" \
      --partition=normal.168h \
      --cpus-per-task=4 \
      --mem=16G \
      --output="${LOG_DIR}/figH_prep_%j.out" \
      --error="${LOG_DIR}/figH_prep_%j.err" \
      --parsable \
      --wrap="${MAMBA} run -n ${ENV_PYTHON} python3 ${SCRIPTS_DIR}/fig_H_widespread_cooption.py")
    echo "[SUBMIT] figH_prep job=$JID"

    # Step 2: phylogenetic analysis (depends on step 1)
    OUTBASE="/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"
    JID2=$(sbatch \
      --job-name="figH_phylo" \
      --partition=normal.168h \
      --cpus-per-task=8 \
      --mem=16G \
      --dependency="afterok:${JID}" \
      --output="${LOG_DIR}/figH_phylo_%j.out" \
      --error="${LOG_DIR}/figH_phylo_%j.err" \
      --parsable \
      --wrap="${MAMBA} run -n esptreeEDFTE bash -c \"\
        for FAA in ${OUTBASE}/figH_*_seqs.faa; do \
          [ -s \\\"\\\$FAA\\\" ] || continue; \
          PREFIX=\\\"\\\${FAA%.faa}\\\"; \
          mafft --auto --thread 8 \\\"\\\$FAA\\\" > \\\"\\\${PREFIX}_aligned.faa\\\"; \
          iqtree -s \\\"\\\${PREFIX}_aligned.faa\\\" -m TEST -B 1000 -T 8 \
            --prefix \\\"\\\${PREFIX}_iqtree\\\" --seed 42 -redo 2>&1; \
        done\"")
    echo "[SUBMIT] figH_phylo job=$JID2 (after $JID)"

    # Step 3: R plot (depends on step 2)
    sbatch \
      --job-name="figH_plot" \
      --partition=normal.168h \
      --cpus-per-task=2 \
      --mem=8G \
      --dependency="afterok:${JID2}" \
      --output="${LOG_DIR}/figH_plot_%j.out" \
      --error="${LOG_DIR}/figH_plot_%j.err" \
      --wrap="${MAMBA} run -n ${ENV_R} Rscript ${SCRIPTS_DIR}/fig_H_plot_phylo.R"
    ;;
esac

echo ""
echo "Done submitting. Check logs in: $LOG_DIR"
echo "Monitor: squeue -u \$USER"
