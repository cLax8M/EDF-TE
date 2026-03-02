#!/usr/bin/env bash
# =============================================================================
# 05_gene_trees_array.sh — Individual gene tree inference with IQTree (SLURM array)
#
# PURPOSE : Infers a maximum-likelihood gene tree for each quality-filtered
#           HOG alignment. Gene trees are used to:
#           (a) detect ancient paralogs via anomalous branch lengths (Step 06)
#           (b) provide a per-locus phylogenetic signal assessment
#
# MODEL   : LG+G4
#   - LG: best-fit protein substitution matrix for diverse eukaryotes
#   - G4: discrete Gamma (4 categories) for rate heterogeneity
#   - Simpler than LG+C60 — appropriate for single-gene trees
#   - Ref: Le & Gascuel 2008 (Mol Biol Evol)
#
# BOOTSTRAP: Ultrafast Bootstrap (UFBoot), 1000 replicates
#   - Provides approximate Bayesian-like branch support quickly
#   - Values ≥ 95 considered well-supported
#   - Ref: Minh et al. 2013 (Mol Biol Evol)
#
# USAGE   : # Generate input list:
#           ls /data/carlos/EDF-TE/SpeciesTree/04_filtered/alignments/*.fa \
#               > /data/carlos/EDF-TE/SpeciesTree/05_gene_trees/input_list.txt
#
#           # Submit:
#           sbatch 05_gene_trees_array.sh
#
# OUTPUT  : 05_gene_trees/<HOG>.treefile    — ML tree with branch lengths
#           05_gene_trees/<HOG>.iqtree       — full IQTree log
#           05_gene_trees/<HOG>.log          — run log
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# =============================================================================

#SBATCH --job-name=iqtree_genes
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/carlos/EDF-TE/SpeciesTree/logs/05_genetrees_%A_%a.out
#SBATCH --error=/data/carlos/EDF-TE/SpeciesTree/logs/05_genetrees_%A_%a.err
#SBATCH --array=1-500%30  # Adjust upper bound; max 30 concurrent to avoid overloading

set -euo pipefail

# ─── Configuration ────────────────────────────────────────────────────────────
SCRIPT_DIR="/data/carlos/EDF-TE/EDF-TE/SpeciesTree/scripts"
source "${SCRIPT_DIR}/config.sh"

FILTERED_DIR="${DIR_FILTERED}/alignments"
OUT_DIR="${DIR_GENE_TREES}"
# Allow overriding input list via environment variable (used by retry runs)
INPUT_LIST="${INPUT_LIST_OVERRIDE:-${OUT_DIR}/input_list.txt}"
THREADS=4

# ─── Setup ────────────────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"

# ─── Generate input list if needed ────────────────────────────────────────────
if [[ ! -f "${INPUT_LIST}" ]]; then
    ls "${FILTERED_DIR}"/*.fa | sort > "${INPUT_LIST}"
    N=$(wc -l < "${INPUT_LIST}")
    echo "[05] Found ${N} filtered alignments. Adjust --array=1-${N} in #SBATCH header."
fi

# ─── Get input file for this array task ───────────────────────────────────────
INPUT_ALN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INPUT_LIST}")

if [[ -z "${INPUT_ALN}" || ! -f "${INPUT_ALN}" ]]; then
    echo "[ERROR] No file for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

HOG_ID=$(basename "${INPUT_ALN}" .trim.fa)
PREFIX="${OUT_DIR}/${HOG_ID}"

# Skip if already done
if [[ -f "${PREFIX}.treefile" ]]; then
    echo "[05] SKIP: ${HOG_ID}.treefile already exists"
    exit 0
fi

echo "[05] Gene tree: ${HOG_ID}"
echo "     Input   : ${INPUT_ALN}"
echo "     Prefix  : ${PREFIX}"
echo "     Started : $(date)"

N_SEQS=$(grep -c "^>" "${INPUT_ALN}")
echo "     Sequences: ${N_SEQS}"

# ─── Run IQTree gene tree ──────────────────────────────────────────────────────
# -s  : input alignment
# -m LG+G4 : substitution model (protein, rate heterogeneity)
# -B 1000  : ultrafast bootstrap replicates
# -T 4     : threads
# --prefix : output file prefix
# -redo    : overwrite if partial run exists
# NOTE: use 'iqtree' (v3.0.1 installed), NOT 'iqtree2'

iqtree \
    -s "${INPUT_ALN}" \
    -m LG+G4 \
    -B 1000 \
    -T "${THREADS}" \
    --prefix "${PREFIX}" \
    -redo \
    --quiet

# Verify output
if [[ ! -f "${PREFIX}.treefile" ]]; then
    echo "[ERROR] IQTree failed to produce treefile for ${HOG_ID}"
    exit 1
fi

echo "[05] DONE: ${HOG_ID}"
echo "     Tree    : ${PREFIX}.treefile"
echo "     Finished: $(date)"

# ─── Quick branch length sanity check ─────────────────────────────────────────
# Extract max branch length — very long branches suggest paralogs
# Uses regex on Newick text (no ete3 dependency)
MAX_BL=$(python3 -c "
import re, sys
try:
    nwk = open('${PREFIX}.treefile').read()
    bls = [float(x) for x in re.findall(r':([0-9]+\.?[0-9]*(?:e[+-]?[0-9]+)?)', nwk)]
    print(f'{max(bls):.4f}') if bls else print('0')
except Exception as e:
    print(f'WARN: {e}', file=sys.stderr); print('0')
" 2>/dev/null || echo "0")
echo "     Max branch length: ${MAX_BL}"
