#!/bin/bash
# Run Section2 analyses that can be executed with current data.
# Produces provisional outputs while EarlGrey/tree are still incomplete.

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate earlgrey

SCRIPT_DIR="/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts"
ANN_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape"
GF_SUMMARY="/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures_V2/ALL_SPECIES_ALL_SUMMARY.txt"
TREE_CANDIDATES=(
  "/data/carlos/EDF-TE/SpeciesTree/09_final_tree/final.treefile"
  "/data/carlos/EDF-TE/SpeciesTree/09_final_tree/final.nwk"
  "/legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/09_final_tree/final.treefile"
  "/legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/09_final_tree/final.nwk"
)

mkdir -p "${OUT_DIR}"

echo "[RUN] B01 Shannon (provisional)"
python3 "${SCRIPT_DIR}/B01_shannon_diversity.py" \
  "${ANN_DIR}" \
  "${OUT_DIR}/B01_shannon_diversity_provisional.tsv"

echo "[RUN] C01 Kimura (provisional)"
python3 "${SCRIPT_DIR}/C01_kimura_metrics.py" \
  "${ANN_DIR}" \
  "${OUT_DIR}/C01_kimura_species_provisional.tsv" \
  "${OUT_DIR}/C01_kimura_superfamily_provisional.tsv"

echo "[RUN] D01 Merge TE + genomics (provisional)"
python3 "${SCRIPT_DIR}/D01_merge_te_genomics.py" \
  "${ANN_DIR}" \
  "${GF_SUMMARY}" \
  "${OUT_DIR}/D01_te_genomics_merged_provisional.tsv"

TREE_FILE=""
for t in "${TREE_CANDIDATES[@]}"; do
  if [[ -f "${t}" ]]; then
    TREE_FILE="${t}"
    break
  fi
done

if [[ -z "${TREE_FILE}" ]]; then
  echo "[INFO] Final tree not found yet. Skipping B02, D02, A01 for now."
  exit 0
fi

echo "[INFO] Using tree: ${TREE_FILE}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate treewasEDFTE || true

echo "[RUN] B02 Pagel lambda"
Rscript "${SCRIPT_DIR}/B02_pagel_lambda.R" \
  "${OUT_DIR}/B01_shannon_diversity_provisional.tsv" \
  "${TREE_FILE}" \
  "${OUT_DIR}/B02_pagel_lambda_provisional.tsv"

echo "[RUN] D02 PGLS"
Rscript "${SCRIPT_DIR}/D02_pgls_regression.R" \
  "${OUT_DIR}/D01_te_genomics_merged_provisional.tsv" \
  "${TREE_FILE}" \
  "${OUT_DIR}/D02_pgls_provisional"

echo "[RUN] A01 ancestral reconstruction"
Rscript "${SCRIPT_DIR}/A01_ancestral_reconstruction.R" \
  "${ANN_DIR}" \
  "${TREE_FILE}" \
  "${OUT_DIR}/A01_ancestral_provisional"

echo "[DONE] Section2 pending/provisional pipeline finished"
