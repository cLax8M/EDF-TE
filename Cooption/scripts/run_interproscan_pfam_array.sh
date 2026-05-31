#!/bin/bash
# ============================================================
# run_interproscan_pfam_array.sh
# Re-run InterProScan (Pfam only) for species whose existing
# IPR annotation lacks Pfam entries (Gene3D/PRINTS/Coils only).
#
# Output per species:
#   {SP}_pfam.tsv  — raw InterProScan TSV output (Pfam only)
#   {SP}_pfam.domains.tsv — normalized domains.tsv (appended to existing)
#
# Usage:
#   sbatch --array=1-37 run_interproscan_pfam_array.sh
# ============================================================

#SBATCH --job-name=ipr_pfam
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --array=1-37
#SBATCH --output=/data/carlos/EDF-TE/logs/interproscan_pfam/ipr_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/interproscan_pfam/ipr_%A_%a.err

set -euo pipefail

INTERPROSCAN="/opt/easybuild/software/InterProScan/5.77-108.0-GCC-13.3.0/interproscan.sh"
PROTEINS_DIR="/legserv/Temp/Carlos/EDF-TE/data/proteins_annot_faa"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/results/Co-option/per_species/interproscan_pfam"
DOMAINS_DIR="/legserv/Temp/Carlos/EDF-TE/results/Co-option/per_species/domains"
SPECIES_LIST="/data/carlos/EDF-TE/EDF-TE/Cooption/scripts/species_interproscan_pfam.txt"
NORMALIZE_SCRIPT="/data/carlos/EDF-TE/EDF-TE/Cooption/scripts/00_normalize_interpro.py"
TMPDIR_BASE="/data/carlos/EDF-TE/tmp/interproscan_pfam"

mkdir -p "${OUT_DIR}" "${TMPDIR_BASE}"

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

FAA="${PROTEINS_DIR}/${SP}.aa.fasta"
if [[ ! -f "${FAA}" ]]; then
    echo "[ERROR] Protein file not found: ${FAA}"
    exit 1
fi

# Skip if already done
if [[ -f "${OUT_DIR}/${SP}_pfam.tsv" ]]; then
    echo "[SKIP] Already processed: ${SP}"
    exit 0
fi

TMPDIR="${TMPDIR_BASE}/${SP}_$$"
mkdir -p "${TMPDIR}"

echo "[START] ${SP} — $(grep -c '>' ${FAA}) proteins"

# Remove stop-codon asterisks (InterProScan chokes on them)
CLEAN_FAA="${TMPDIR}/${SP}_clean.faa"
sed 's/\*//g' "${FAA}" > "${CLEAN_FAA}"

# Run InterProScan with Pfam only — use custom properties to supply missing release versions
export INTERPROSCAN_CONF="/data/carlos/EDF-TE/EDF-TE/Cooption/scripts/interproscan_pfam.properties"
echo "[IPR] Running InterProScan -appl Pfam for ${SP}"
"${INTERPROSCAN}" \
    --input "${CLEAN_FAA}" \
    --output-file-base "${OUT_DIR}/${SP}_pfam" \
    --formats TSV \
    --applications Pfam \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    --disable-precalc \
    --tempdir "${TMPDIR}" \
    --goterms \
    --iprlookup \
    2>&1

echo "[OK] InterProScan done: ${OUT_DIR}/${SP}_pfam.tsv"

# Normalize the new Pfam TSV and merge into existing domains.tsv
# The new IPR5 TSV has columns:
# protein_accession | md5 | length | analysis | sig_accession | sig_desc |
# start | stop | score | status | date | ipr_accession | ipr_desc | go_terms | pathways
RAW_TSV="${OUT_DIR}/${SP}_pfam.tsv"
DOMAINS_TSV="${DOMAINS_DIR}/${SP}.domains.tsv"

if [[ ! -s "${RAW_TSV}" ]]; then
    echo "[WARN] Empty output for ${SP} — no Pfam hits found"
    rm -rf "${TMPDIR}"
    exit 0
fi

echo "[NORM] Appending Pfam entries to ${DOMAINS_TSV}"
PYTHON3="/data/carlos/EDF-TE/PanGenome/.venv/bin/python3"
"${PYTHON3}" - "${RAW_TSV}" "${DOMAINS_TSV}" <<'PYEOF'
import sys, pandas as pd

raw = sys.argv[1]
out = sys.argv[2]

df = pd.read_csv(raw, sep="\t", header=None, dtype=str).fillna("")
# IPR5 TSV columns (0-indexed):
# 0=protein_accession, 1=md5, 2=length, 3=analysis, 4=sig_accession,
# 5=sig_desc, 6=start, 7=stop, 8=score, 9=status, 10=date,
# 11=ipr_accession, 12=ipr_desc, 13=go_terms, 14=pathways
rows = []
for _, r in df.iterrows():
    gene_id     = r[0]
    source_db   = r[3]   # "Pfam"
    domain_id   = r[4]   # PFxxxxx accession
    domain_name = r[5]   # domain description
    start       = r[6]
    end         = r[7]
    rows.append([gene_id, domain_id, domain_name, source_db, start, end])

new_df = pd.DataFrame(rows, columns=["gene_id","domain_id","domain_name","source_db","start","end"])

import os
if os.path.exists(out) and os.path.getsize(out) > 0:
    existing = pd.read_csv(out, sep="\t", dtype=str).fillna("")
    merged = pd.concat([existing, new_df], ignore_index=True)
else:
    merged = new_df

merged.to_csv(out, sep="\t", index=False)
print(f"[OK] Appended {len(new_df)} Pfam rows to {out}")
PYEOF

rm -rf "${TMPDIR}"
echo "[DONE] ${SP}"
