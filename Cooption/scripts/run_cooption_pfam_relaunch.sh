#!/bin/bash
# ============================================================
# run_cooption_pfam_relaunch.sh
# SLURM array: re-run co-option Phase 2 (steps 04-06) for the
# 37 species that previously ran InterProScan without Pfam and
# therefore produced empty true_TE_coopted_genes.tsv files.
#
# PRE-REQUISITE: InterProScan Pfam job (3056288) must have
# finished and its output merged into the per-species .domains.tsv
# files before submitting this script.
#
# What this script does differently from run_cooption_phase2.sh:
#   - Restricted to array indices for the 37 Pfam-relaunch species
#   - DELETES the existing (empty) true_TE_coopted_genes.tsv and
#     gene_classes.tsv for the selected species BEFORE running,
#     bypassing the skip-if-exists guard.
#
# Usage:
#   # after job 3056288 finishes:
#   sbatch run_cooption_pfam_relaunch.sh
# ============================================================

#SBATCH --job-name=cooption_pfam_relaunch
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1,2,31,38,47,53,59,65,70,72,79,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,111,119,172,180,232,233
#SBATCH --output=/data/carlos/EDF-TE/logs/cooption_phase2/pfam_relaunch_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/cooption_phase2/pfam_relaunch_%A_%a.err

set -euo pipefail

# ── environment ───────────────────────────────────────────────
eval "$(micromamba shell hook --shell bash)"
micromamba activate earlgrey

# ── paths ───────────────────────────────────────────────────
SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/Cooption/scripts"
COOPTION_DIR="/legserv/Temp/Carlos/EDF-TE/results/Co-option"
FUSIONS_DIR="${COOPTION_DIR}/per_species/fusions"
DOMAINS_DIR="${COOPTION_DIR}/per_species/domains"
CLASS_DIR="${COOPTION_DIR}/per_species/classification"
SPECIES_LIST="${SCRIPTS_DIR}/species_phase2.txt"

# ── select species for this array task ───────────────────────
SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

echo "[START] ${SP} (Pfam relaunch)"

# ── remove stale output files to bypass skip-if-exists guard ─
OUT_GENE="${CLASS_DIR}/${SP}.gene_classes.tsv"
OUT_COOPTED="${CLASS_DIR}/${SP}.true_TE_coopted_genes.tsv"

[[ -f "${OUT_GENE}" ]]    && rm -f "${OUT_GENE}"    && echo "[RESET] Removed ${OUT_GENE}"
[[ -f "${OUT_COOPTED}" ]] && rm -f "${OUT_COOPTED}" && echo "[RESET] Removed ${OUT_COOPTED}"

# ── validate inputs ───────────────────────────────────────────
FUSION_FILE="${FUSIONS_DIR}/${SP}.TE_host_fusions.tsv"
HOST_DOM_FILE="${FUSIONS_DIR}/${SP}.host_domains.tsv"
DOMAINS_FILE="${DOMAINS_DIR}/${SP}.domains.tsv"

if [[ ! -f "${FUSION_FILE}" ]]; then
    echo "[ERROR] Missing fusion file: ${FUSION_FILE}"
    exit 1
fi
if [[ ! -f "${DOMAINS_FILE}" ]]; then
    echo "[ERROR] Missing domains file (InterProScan Pfam output not yet merged?): ${DOMAINS_FILE}"
    exit 1
fi

# ── step 04: extract host domains ────────────────────────────
echo "[04] Extracting host domains: ${SP}"
python3 "${SCRIPTS_DIR}/04_extract_host_domains.py" \
    "${FUSION_FILE}" \
    "${HOST_DOM_FILE}"

# ── step 05: classify genes ───────────────────────────────────
echo "[05] Classifying genes: ${SP}"
python3 "${SCRIPTS_DIR}/05_classify_genes.py" \
    "${DOMAINS_FILE}" \
    "${OUT_GENE}"

# ── step 06: filter true co-option ───────────────────────────
echo "[06] Filtering true co-option: ${SP}"
python3 "${SCRIPTS_DIR}/06_filter_true_cooption.py" \
    "${OUT_GENE}" \
    "${OUT_COOPTED}"

echo "[DONE] ${SP} — $(wc -l < "${OUT_COOPTED}") co-opted genes (incl. header)"

# ── SETMAR control check (human only) ────────────────────────
if [[ "${SP}" == "HO_Homosap" ]]; then
    echo ""
    echo "=== SETMAR CONTROL CHECK ==="
    if grep -qi "SET\|DDE\|Mariner\|Tc1\|PF01609" "${OUT_COOPTED}" 2>/dev/null; then
        echo "[PASS] SETMAR or related domains found in HO_Homosap — Pfam integration OK"
    else
        echo "[WARN] SETMAR/DDE_Tnp not found in HO_Homosap — check InterProScan output"
    fi
fi
