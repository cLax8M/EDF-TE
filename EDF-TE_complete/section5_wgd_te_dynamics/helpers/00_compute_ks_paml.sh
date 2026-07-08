#!/usr/bin/env bash

# Compute Ks using MAFFT + PAL2NAL + PAML yn00 for one species.
# Usage:
#   03_compute_ks_paml.sh <species>
#
# Expects:
#   - results/ks/paralog_pairs_for_ks.filtered.tsv  (species, gene1, gene2, ...)
#   - results/ks/cds/<species>.cds.fa
#   - proteins under ${PROTEIN_DIR} with usual extensions (.faa/.fa/.fasta/.aa.fasta/.aa.fa)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../core/config_wgd.sh"

if [ "$#" -ne 1 ]; then
  echo "Usage: $(basename "$0") <species>" >&2
  exit 1
fi

SP="$1"

PARALOG_TSV="${RESULTS_DIR}/ks/paralog_pairs_for_ks.filtered.tsv"
CDS_DIR="${RESULTS_DIR}/ks/cds"
OUT_ROOT="${RESULTS_DIR}/ks/paml_yn00"

if [ ! -f "${PARALOG_TSV}" ]; then
  echo "[WGD][KsPAML] ERROR: filtered paralog TSV not found at ${PARALOG_TSV}" >&2
  exit 1
fi

mkdir -p "${OUT_ROOT}"
LOG_DIR="${OUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"
SP_OUT="${OUT_ROOT}/${SP}"
mkdir -p "${SP_OUT}"
PAIRS_DIR="${SP_OUT}/pairs"
mkdir -p "${PAIRS_DIR}"

LOG_FILE="${LOG_DIR}/${SP}.log"
: > "${LOG_FILE}"

# Output Ks table path
KS_TABLE="${SP_OUT}/${SP}.ks.tsv"

# If Ks table already exists with data rows, skip recomputation for this species
if [ -f "${KS_TABLE}" ] && [ "$(wc -l < "${KS_TABLE}")" -gt 1 ]; then
  echo "[WGD][KsPAML] INFO: Ks table already present for ${SP} at ${KS_TABLE}; skipping recomputation." | tee -a "${LOG_FILE}" >&2
  exit 0
fi

# Check tools
for bin in mafft pal2nal.pl yn00; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    echo "[WGD][KsPAML] ERROR: required binary '${bin}' not found in PATH." | tee -a "${LOG_FILE}" >&2
    exit 1
  fi
done

# Resolve protein FASTA
PROT=""
for ext in .faa .fa .fasta .aa.fasta .aa.fa; do
  CAND="${PROTEIN_DIR}/${SP}${ext}"
  if [ -s "${CAND}" ]; then
    PROT="${CAND}"
    break
  fi
done

if [ -z "${PROT}" ]; then
  echo "[WGD][KsPAML] WARNING: protein FASTA not found or empty for ${SP} under ${PROTEIN_DIR}; skipping species." | tee -a "${LOG_FILE}" >&2
  exit 0
fi

CDS_FASTA="${CDS_DIR}/${SP}.cds.fa"
if [ ! -s "${CDS_FASTA}" ]; then
  echo "[WGD][KsPAML] WARNING: CDS FASTA missing or empty for ${SP} at ${CDS_FASTA}; skipping species." | tee -a "${LOG_FILE}" >&2
  exit 0
fi

# Build per-species pairs file (gene1, gene2)
PAIRS_TSV="${PAIRS_DIR}/${SP}.pairs.tsv"
awk -v sp="${SP}" 'NR>1 && $1==sp {print $2"\t"$3}' "${PARALOG_TSV}" > "${PAIRS_TSV}" || true

if [ ! -s "${PAIRS_TSV}" ]; then
  echo "[WGD][KsPAML] WARNING: no paralog pairs for species ${SP}; nothing to do." | tee -a "${LOG_FILE}" >&2
  exit 0
fi

# Diagnose ID coverage (which genes from pairs have both protein and CDS sequences)
python "${SCRIPT_DIR}/diagnose_id_coverage.py" \
  --species "${SP}" \
  --proteins "${PROT}" \
  --cds "${CDS_FASTA}" \
  --pairs "${PAIRS_TSV}" \
  --out "${SP_OUT}/${SP}.id_diagnostics.tsv" >>"${LOG_FILE}" 2>&1 || DIAG_STATUS=$?

if [ "${DIAG_STATUS:-0}" -ne 0 ] && [ "${DIAG_STATUS:-0}" -ne 2 ]; then
  echo "[WGD][KsPAML] ERROR: diagnose_id_coverage.py failed for ${SP} with status ${DIAG_STATUS}; aborting species." | tee -a "${LOG_FILE}" >&2
  exit 1
fi

if [ "${DIAG_STATUS:-0}" -eq 2 ]; then
  echo "[WGD][KsPAML] WARNING: no IDs with both protein and CDS for ${SP}; skipping species before Ks." | tee -a "${LOG_FILE}" >&2
  exit 0
fi

# Prepare per-pair FASTAs using Python helper (will skip only truly missing IDs)
python "${SCRIPT_DIR}/prepare_paml_pair_fastas.py" \
  --proteins "${PROT}" \
  --cds "${CDS_FASTA}" \
  --pairs "${PAIRS_TSV}" \
  --out-dir "${PAIRS_DIR}" >>"${LOG_FILE}" 2>&1

MAPPING_TSV="${PAIRS_DIR}/pairs_index.tsv"
if [ ! -s "${MAPPING_TSV}" ]; then
  echo "[WGD][KsPAML] WARNING: no valid pairs after sequence preparation for ${SP}; nothing to do." | tee -a "${LOG_FILE}" >&2
  exit 0
fi

if [ ! -f "${KS_TABLE}" ]; then
  echo -e "species\tgene1\tgene2\tKs" > "${KS_TABLE}"
fi

# Process each pair sequentially (one species per SLURM job; parallelize across species with arrays).

tail -n +2 "${MAPPING_TSV}" | while IFS=$'\t' read -r PAIR_ID G1 G2; do
  P_BASE="${PAIRS_DIR}/${PAIR_ID}"
  P_PEP="${P_BASE}.pep.fa"
  P_CDS="${P_BASE}.cds.fa"
  if [ ! -s "${P_PEP}" ] || [ ! -s "${P_CDS}" ]; then
    echo "[WGD][KsPAML] WARNING: missing per-pair FASTA for ${SP} ${G1},${G2}; skipping pair." >>"${LOG_FILE}"
    continue
  fi

  P_PEP_ALN="${P_BASE}.pep.aln.fa"
  P_CODON_ALN="${P_BASE}.codon.phy"
  P_CTL="${P_BASE}.yn00.ctl"
  P_OUT="${P_BASE}.yn00.out"

  # 3. Align proteins with MAFFT
  if ! mafft --auto --quiet "${P_PEP}" > "${P_PEP_ALN}" 2>>"${LOG_FILE}"; then
    echo "[WGD][KsPAML] WARNING: MAFFT failed for ${SP} ${G1},${G2}; skipping pair." >>"${LOG_FILE}"
    continue
  fi
  if [ ! -s "${P_PEP_ALN}" ]; then
    echo "[WGD][KsPAML] WARNING: empty protein alignment for ${SP} ${G1},${G2}; skipping pair." >>"${LOG_FILE}"
    continue
  fi

  # 4. Convert to codon alignment with PAL2NAL
  if ! pal2nal.pl "${P_PEP_ALN}" "${P_CDS}" -output paml > "${P_CODON_ALN}" 2>>"${LOG_FILE}"; then
    echo "[WGD][KsPAML] WARNING: PAL2NAL failed for ${SP} ${G1},${G2}; skipping pair." >>"${LOG_FILE}"
    continue
  fi
  if [ ! -s "${P_CODON_ALN}" ]; then
    echo "[WGD][KsPAML] WARNING: empty codon alignment for ${SP} ${G1},${G2}; skipping pair." >>"${LOG_FILE}"
    continue
  fi

  # 5. Run PAML yn00
  cat > "${P_CTL}" <<EOF
seqfile = ${P_CODON_ALN}
outfile = ${P_OUT}
verbose = 0
icode = 0
weighting = 0
commonf3x4 = 0
EOF

  if ! yn00 "${P_CTL}" >>"${LOG_FILE}" 2>&1; then
    echo "[WGD][KsPAML] WARNING: yn00 failed for ${SP} ${G1},${G2}; skipping pair." >>"${LOG_FILE}"
    continue
  fi

  if [ ! -s "${P_OUT}" ]; then
    echo "[WGD][KsPAML] WARNING: yn00 output missing for ${SP} ${G1},${G2}; skipping pair." >>"${LOG_FILE}"
    continue
  fi

  # 6. Parse Ks (dS) from yn00 output
  KS=$(awk '/dS =/ {ks=$3} END { if (ks=="" || ks==0) {print "NA"} else {print ks} }' "${P_OUT}" 2>/dev/null || echo "NA")

  echo -e "${SP}\t${G1}\t${G2}\t${KS}" >> "${KS_TABLE}"

done

echo "[WGD][KsPAML] Finished Ks computation for species ${SP}." | tee -a "${LOG_FILE}" >&2
