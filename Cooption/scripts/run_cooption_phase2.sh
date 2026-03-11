#!/bin/bash
# ============================================================
# run_cooption_phase2.sh
# SLURM array: run co-option Phase 2 (steps 04-06) for all
# species that have a TE_host_fusions.tsv file.
#
# Steps:
#   04 - extract host domains (fusion → host_domains.tsv)
#   05 - classify genes        (domains.tsv → gene_classes.tsv)
#   06 - filter true co-option (gene_classes → true_TE_coopted_genes.tsv)
#
# Validation control: after all jobs finish, verify that
#   HO_Homosap.true_TE_coopted_genes.tsv contains SETMAR
#   (gene with DDE_Tnp_IS240/PF01609 + SET domain → canonical
#   human co-option event; if missing, check TE_domains_master.tsv)
#
# Usage:
#   bash run_cooption_phase2.sh [--dry-run]
#   sbatch run_cooption_phase2.sh
# ============================================================

#SBATCH --job-name=cooption_phase2
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-238
#SBATCH --output=/data/carlos/EDF-TE/logs/cooption_phase2/phase2_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/cooption_phase2/phase2_%A_%a.err

set -euo pipefail

# ── paths ───────────────────────────────────────────────────
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COOPTION_DIR="/legserv/Temp/Carlos/EDF-TE/results/Co-option"
FUSIONS_DIR="${COOPTION_DIR}/per_species/fusions"
DOMAINS_DIR="${COOPTION_DIR}/per_species/domains"
CLASS_DIR="${COOPTION_DIR}/per_species/classification"
SPECIES_LIST="${SCRIPTS_DIR}/species_phase2.txt"

# ── dry-run mode ─────────────────────────────────────────────
DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "[DRY-RUN] Listing species that would be processed:"
fi

# ── build species list (first run only) ──────────────────────
if [[ ! -f "${SPECIES_LIST}" ]]; then
    ls "${FUSIONS_DIR}"/*.TE_host_fusions.tsv 2>/dev/null \
        | xargs -n1 basename \
        | sed 's/\.TE_host_fusions\.tsv//' \
        | sort > "${SPECIES_LIST}"
    echo "[INFO] Species list written: ${SPECIES_LIST} ($(wc -l < "${SPECIES_LIST}") species)"
fi

# ── dry-run: just list and exit ───────────────────────────────
if [[ "${DRY_RUN}" == true ]]; then
    echo "[DRY-RUN] Would process $(wc -l < "${SPECIES_LIST}") species:"
    cat "${SPECIES_LIST}"
    already_done=$(comm -12 \
        <(sort "${SPECIES_LIST}") \
        <(ls "${CLASS_DIR}"/*.gene_classes.tsv 2>/dev/null \
            | xargs -n1 basename | sed 's/\.gene_classes\.tsv//' | sort) \
        | wc -l)
    echo "[DRY-RUN] Already done: ${already_done} — would skip these."
    exit 0
fi

# ── select species for this array task ───────────────────────
SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# ── skip if already complete ──────────────────────────────────
OUT_GENE="${CLASS_DIR}/${SP}.gene_classes.tsv"
OUT_COOPTED="${CLASS_DIR}/${SP}.true_TE_coopted_genes.tsv"
if [[ -f "${OUT_COOPTED}" ]]; then
    echo "[SKIP] ${SP} already complete: ${OUT_COOPTED}"
    exit 0
fi

FUSION_FILE="${FUSIONS_DIR}/${SP}.TE_host_fusions.tsv"
HOST_DOM_FILE="${FUSIONS_DIR}/${SP}.host_domains.tsv"
DOMAINS_FILE="${DOMAINS_DIR}/${SP}.domains.tsv"

# ── validate inputs ───────────────────────────────────────────
if [[ ! -f "${FUSION_FILE}" ]]; then
    echo "[ERROR] Missing fusion file: ${FUSION_FILE}"
    exit 1
fi
if [[ ! -f "${DOMAINS_FILE}" ]]; then
    echo "[ERROR] Missing domains file: ${DOMAINS_FILE}"
    exit 1
fi

echo "[START] ${SP}"

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

echo "[DONE] ${SP} — $(wc -l < "${OUT_COOPTED}") co-opted genes"

# ── SETMAR control check (runs once, for human only) ─────────
if [[ "${SP}" == "HO_Homosap" ]]; then
    echo ""
    echo "=== SETMAR CONTROL CHECK ==="
    if grep -qi "SET\|DDE\|Mariner\|Tc1\|PF01609" "${OUT_COOPTED}" 2>/dev/null; then
        echo "[PASS] SETMAR-like entry detected in HO_Homosap co-opted genes"
    else
        echo "[WARN] No SETMAR-like entry in HO_Homosap.true_TE_coopted_genes.tsv"
        echo "[WARN] Check TE_domains_master.tsv for PF01609 (DDE_Tnp_IS240 / Mariner transposase)"
        echo "[WARN] Also verify SETMAR gene ID in HO_Homosap_IPR.tab"
    fi
    echo "============================="
fi
