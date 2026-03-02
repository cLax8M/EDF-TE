#!/usr/bin/env bash
# =============================================================================
# 02_align_mafft_array.sh — Multiple sequence alignment with MAFFT (SLURM array)
#
# PURPOSE : Aligns each unaligned HOG FASTA from 01_markers/ using MAFFT
#           with --globalpair (L-INS-I algorithm), the most accurate MAFFT
#           strategy. Uses a SLURM job array (one job per FASTA file).
#
# ALGORITHM: L-INS-i (--globalpair --maxiterate 1000)
#   - Best for sequences with conserved domains and variable regions
#   - More accurate than FFT-NS-2 (OrthoFinder default) for diverse fungi
#   - Recommended for phylogenomics with high taxonomic diversity
#   - Ref: Katoh & Standley 2013 (Mol Biol Evol)
#
# USAGE   : # First, generate the input file list:
#           ls /data/carlos/EDF-TE/SpeciesTree/01_markers/*.fa \
#               > /data/carlos/EDF-TE/SpeciesTree/02_aligned/input_list.txt
#
#           # Then submit the array:
#           sbatch 02_align_mafft_array.sh
#
# OUTPUT  : 02_aligned/<HOG_ID>.aln.fa  — aligned FASTA per HOG
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# =============================================================================

#SBATCH --job-name=mafft_align
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/data/carlos/EDF-TE/SpeciesTree/logs/02_mafft_%A_%a.out
#SBATCH --error=/data/carlos/EDF-TE/SpeciesTree/logs/02_mafft_%A_%a.err
#SBATCH --array=1-700%50  # Adjust upper bound to actual number of markers; max 50 concurrent

set -euo pipefail

# ─── Configuration ────────────────────────────────────────────────────────────
SCRIPT_DIR="/data/carlos/EDF-TE/EDF-TE/SpeciesTree/scripts"
source "${SCRIPT_DIR}/config.sh"

MARKERS_DIR="${DIR_MARKERS}"
OUT_DIR="${DIR_ALIGNED}"
INPUT_LIST="${OUT_DIR}/input_list.txt"
THREADS=${MAFFT_THREADS}   # 4 threads (from config.sh)

# ─── Setup ────────────────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "$(dirname "${OUT_DIR}")")/logs"

# Activate micromamba environment
eval "$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"

# ─── Generate input list if it doesn't exist ──────────────────────────────────
# (Run once before submitting array)
if [[ ! -f "${INPUT_LIST}" ]]; then
    echo "[02] Generating input list..."
    ls "${MARKERS_DIR}"/*.fa | sort > "${INPUT_LIST}"
    N=$(wc -l < "${INPUT_LIST}")
    echo "[02] Found ${N} marker FASTAs. Adjust --array=1-${N} in #SBATCH header."
fi

# ─── Get input file for this array task ───────────────────────────────────────
INPUT_FA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INPUT_LIST}")

if [[ -z "${INPUT_FA}" || ! -f "${INPUT_FA}" ]]; then
    echo "[ERROR] No file for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

HOG_ID=$(basename "${INPUT_FA}" .fa)
OUT_FA="${OUT_DIR}/${HOG_ID}.aln.fa"

# Skip if already done
if [[ -f "${OUT_FA}" ]] && [[ -s "${OUT_FA}" ]]; then
    echo "[02] SKIP: ${HOG_ID}.aln.fa already exists"
    exit 0
fi

# ─── Run MAFFT L-INS-i ────────────────────────────────────────────────────────
# --globalpair : L-INS-i algorithm (most accurate, suitable for <~1000 seqs)
# --maxiterate : iterative refinement cycles (1000 = maximum)
# --amino      : force amino acid mode
# --thread     : number of parallel threads
# --quiet      : suppress progress output (goes to SLURM log)

echo "[02] Aligning: ${HOG_ID}"
echo "     Input   : ${INPUT_FA}"
echo "     Output  : ${OUT_FA}"
echo "     Threads : ${THREADS}"
echo "     Started : $(date)"

N_SEQS=$(grep -c "^>" "${INPUT_FA}")
echo "     Sequences: ${N_SEQS}"

# Use L-INS-i for ≤500 seqs, FFT-NS-i for larger (>500 is unusual at max_copies=2)
if [[ ${N_SEQS} -le 500 ]]; then
    mafft \
        --globalpair \
        --maxiterate 1000 \
        --amino \
        --thread "${THREADS}" \
        --quiet \
        "${INPUT_FA}" > "${OUT_FA}"
else
    echo "[WARN] ${N_SEQS} sequences > 500. Using --auto instead of --globalpair."
    mafft \
        --auto \
        --amino \
        --thread "${THREADS}" \
        --quiet \
        "${INPUT_FA}" > "${OUT_FA}"
fi

# Verify output
if [[ ! -s "${OUT_FA}" ]]; then
    echo "[ERROR] MAFFT produced empty output for ${HOG_ID}"
    exit 1
fi

N_OUT=$(grep -c "^>" "${OUT_FA}")
echo "[02] DONE: ${HOG_ID} — ${N_OUT} sequences aligned"
echo "     Finished: $(date)"
