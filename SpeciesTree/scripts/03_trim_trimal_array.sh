#!/usr/bin/env bash
# =============================================================================
# 03_trim_trimal_array.sh — Alignment trimming with trimAl (SLURM array)
#
# PURPOSE : Trims poorly aligned columns from each MAFFT alignment using
#           trimAl -gappyout. Then filters alignments shorter than MIN_ALN_LEN
#           amino acid columns. Short alignments are moved to a 'rejected/'
#           subdirectory for inspection.
#
# ALGORITHM: trimAl -gappyout
#   - Removes columns with a high proportion of gaps using a heuristic
#     "gappy" approach. Less aggressive than -strict, better suited for
#     highly divergent fungi where many gap columns are evolutionarily real.
#   - Ref: Capella-Gutiérrez et al. 2009 (Bioinformatics)
#
# FILTER: Min alignment length = 60 AA columns after trimming
#   - Ref: Merényi et al. 2023 (Nat. Ecol. Evol.)
#
# USAGE   : # Generate input list first:
#           ls /data/carlos/EDF-TE/SpeciesTree/02_aligned/*.aln.fa \
#               > /data/carlos/EDF-TE/SpeciesTree/03_trimmed/input_list.txt
#
#           # Then submit:
#           sbatch 03_trim_trimal_array.sh
#
# OUTPUT  : 03_trimmed/<HOG_ID>.trim.fa       — trimmed alignment (passed)
#           03_trimmed/rejected/<HOG_ID>.fa   — too short (failed)
#           03_trimmed/trim_stats.tsv         — lengths before/after trimming
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# =============================================================================

#SBATCH --job-name=trimal_trim
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/data/carlos/EDF-TE/SpeciesTree/logs/03_trimal_%A_%a.out
#SBATCH --error=/data/carlos/EDF-TE/SpeciesTree/logs/03_trimal_%A_%a.err
#SBATCH --array=1-700%100  # Adjust to actual number of aligned files

set -euo pipefail

# ─── Configuration ────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

ALIGNED_DIR="${DIR_ALIGNED}"
OUT_DIR="${DIR_TRIMMED}"
INPUT_LIST="${OUT_DIR}/input_list.txt"
REJECTED_DIR="${OUT_DIR}/rejected"
MIN_LEN=${MIN_ALN_LEN}   # 60 AA columns (from config.sh)

# ─── Setup ────────────────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}" "${REJECTED_DIR}"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"

# ─── Generate input list if needed ────────────────────────────────────────────
if [[ ! -f "${INPUT_LIST}" ]]; then
    ls "${ALIGNED_DIR}"/*.aln.fa | sort > "${INPUT_LIST}"
    N=$(wc -l < "${INPUT_LIST}")
    echo "[03] Found ${N} aligned FASTAs. Adjust --array=1-${N} in #SBATCH header."
fi

# ─── Get input file for this array task ───────────────────────────────────────
INPUT_ALN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INPUT_LIST}")

if [[ -z "${INPUT_ALN}" || ! -f "${INPUT_ALN}" ]]; then
    echo "[ERROR] No file for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

HOG_ID=$(basename "${INPUT_ALN}" .aln.fa)
OUT_TRIM="${OUT_DIR}/${HOG_ID}.trim.fa"
STATS_FILE="${OUT_DIR}/${HOG_ID}.trim_stats.txt"

# Skip if already done
if [[ -f "${OUT_TRIM}" ]] && [[ -s "${OUT_TRIM}" ]]; then
    echo "[03] SKIP: ${HOG_ID}.trim.fa already exists"
    exit 0
fi

echo "[03] Trimming: ${HOG_ID}"
echo "     Input   : ${INPUT_ALN}"
echo "     Started : $(date)"

# Get alignment length before trimming
LEN_BEFORE=$(python3 -c "
from Bio import SeqIO
recs = list(SeqIO.parse('${INPUT_ALN}', 'fasta'))
print(len(recs[0].seq) if recs else 0)
")
N_SEQS=$(grep -c "^>" "${INPUT_ALN}")
echo "     Sequences     : ${N_SEQS}"
echo "     Length before : ${LEN_BEFORE} AA columns"

# ─── Run trimAl ───────────────────────────────────────────────────────────────
TEMP_TRIM="${OUT_DIR}/${HOG_ID}.tmp.fa"

trimal \
    -in "${INPUT_ALN}" \
    -out "${TEMP_TRIM}" \
    -gappyout

# Get alignment length after trimming
LEN_AFTER=$(python3 -c "
from Bio import SeqIO
recs = list(SeqIO.parse('${TEMP_TRIM}', 'fasta'))
print(len(recs[0].seq) if recs else 0)
" 2>/dev/null || echo "0")

echo "     Length after  : ${LEN_AFTER} AA columns"
echo "     Min required  : ${MIN_LEN} AA columns"

# Write stats
echo -e "${HOG_ID}\t${N_SEQS}\t${LEN_BEFORE}\t${LEN_AFTER}" > "${STATS_FILE}"

# ─── Filter by minimum alignment length ───────────────────────────────────────
if [[ ${LEN_AFTER} -lt ${MIN_LEN} ]]; then
    echo "[03] REJECTED: ${HOG_ID} — trimmed length ${LEN_AFTER} < ${MIN_LEN} AA"
    mv "${TEMP_TRIM}" "${REJECTED_DIR}/${HOG_ID}.trim.fa"
    # Write rejection record
    echo -e "${HOG_ID}\tREJECTED_TOO_SHORT\t${LEN_AFTER}" >> "${OUT_DIR}/rejected_list.txt"
else
    mv "${TEMP_TRIM}" "${OUT_TRIM}"
    echo "[03] PASS  : ${HOG_ID} — ${LEN_AFTER} AA columns retained"
fi

echo "[03] Done: ${HOG_ID} — $(date)"

# ─── After all jobs complete: merge stats ─────────────────────────────────────
# Run this manually after array completes:
# cat 03_trimmed/*.trim_stats.txt > 03_trimmed/trim_stats_all.tsv
# echo -e "HOG\tn_seqs\tlen_before\tlen_after" | \
#     cat - 03_trimmed/trim_stats_all.tsv > 03_trimmed/trim_stats.tsv
