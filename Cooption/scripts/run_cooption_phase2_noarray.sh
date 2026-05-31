#!/bin/bash
# ============================================================
# run_cooption_phase2_noarray.sh
# Sequential (non-array) execution of co-option phase2 steps 04-06.
# ============================================================

#SBATCH --job-name=cooption_phase2_noarray
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/data/carlos/EDF-TE/logs/cooption_phase2/phase2_noarray_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/cooption_phase2/phase2_noarray_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate earlgrey

SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Cooption/scripts"
COOPTION_DIR="/legserv/Temp/Carlos/EDF-TE/results/Co-option"
FUSIONS_DIR="${COOPTION_DIR}/per_species/fusions"
DOMAINS_DIR="${COOPTION_DIR}/per_species/domains"
CLASS_DIR="${COOPTION_DIR}/per_species/classification"
SPECIES_LIST="${SCRIPTS_DIR}/species_phase2.txt"
mkdir -p "${CLASS_DIR}" "${FUSIONS_DIR}"

if [[ ! -f "${SPECIES_LIST}" ]]; then
    ls "${FUSIONS_DIR}"/*.TE_host_fusions.tsv 2>/dev/null \
        | xargs -n1 basename \
        | sed 's/\.TE_host_fusions\.tsv//' \
        | sort > "${SPECIES_LIST}"
fi

total=$(wc -l < "${SPECIES_LIST}")
done_n=0
fail_n=0

while read -r SP; do
    [[ -z "${SP}" ]] && continue
    OUT_GENE="${CLASS_DIR}/${SP}.gene_classes.tsv"
    OUT_COOPTED="${CLASS_DIR}/${SP}.true_TE_coopted_genes.tsv"
    FUSION_FILE="${FUSIONS_DIR}/${SP}.TE_host_fusions.tsv"
    HOST_DOM_FILE="${FUSIONS_DIR}/${SP}.host_domains.tsv"
    DOMAINS_FILE="${DOMAINS_DIR}/${SP}.domains.tsv"

    if [[ -f "${OUT_COOPTED}" ]]; then
        echo "[SKIP] ${SP} already complete"
        done_n=$((done_n+1))
        continue
    fi
    if [[ ! -f "${FUSION_FILE}" || ! -f "${DOMAINS_FILE}" ]]; then
        echo "[WARN] ${SP} missing inputs"
        fail_n=$((fail_n+1))
        continue
    fi

    echo "[START] ${SP}"
    if python3 "${SCRIPTS_DIR}/04_extract_host_domains.py" "${FUSION_FILE}" "${HOST_DOM_FILE}" \
      && python3 "${SCRIPTS_DIR}/05_classify_genes.py" "${DOMAINS_FILE}" "${OUT_GENE}" \
      && python3 "${SCRIPTS_DIR}/06_filter_true_cooption.py" "${OUT_GENE}" "${OUT_COOPTED}"; then
        echo "[DONE] ${SP}"
        done_n=$((done_n+1))
    else
        echo "[FAIL] ${SP}"
        fail_n=$((fail_n+1))
    fi

done < "${SPECIES_LIST}"

echo "[SUMMARY] total=${total} processed_or_skipped=${done_n} failed=${fail_n}"
