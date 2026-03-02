#!/usr/bin/env bash
# =============================================================================
# 08_species_tree_guide.sh — Guide species tree with LG+F+G4 (SLURM)
#
# PURPOSE : Infers a preliminary species tree under the LG+F+G4 model.
#           This "guide tree" serves two purposes:
#           1. Quality check: compare topology to OrthoFinder's ASTRAL tree
#              and to known fungal relationships
#           2. Starting tree (-t) for the final LG+C60+F+G4 run (Step 09),
#              substantially reducing computation time by ~30-50%
#
# MODEL   : LG+F+G4
#   - LG  : best-fit empirical protein substitution matrix
#   - F   : empirical amino acid frequencies estimated from data
#   - G4  : discrete Gamma rate heterogeneity (4 categories)
#   - Simpler than LG+C60 but sufficient for guide tree topology
#
# EXPECTED TIME : ~6-12 hours for 238 taxa, ~250 genes, ~100k columns
#
# USAGE   : sbatch 08_species_tree_guide.sh
#
# OUTPUT  : 08_guide_tree/guide.treefile    — ML guide tree
#           08_guide_tree/guide.iqtree       — full IQTree log with statistics
#           08_guide_tree/guide.log          — run log
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# REFS    : Liu et al. 2024, Szánthó et al. 2025, Merényi et al. 2023
# =============================================================================

#SBATCH --job-name=iqtree_guide
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=12
#SBATCH --mem=230G
#SBATCH --exclude=LEGcompute1

#SBATCH --output=/data/carlos/EDF-TE/SpeciesTree/logs/08_guide_%j.out
#SBATCH --error=/data/carlos/EDF-TE/SpeciesTree/logs/08_guide_%j.err

set -euo pipefail

# ─── Configuration ────────────────────────────────────────────────────────────
SCRIPT_DIR="/data/carlos/EDF-TE/EDF-TE/SpeciesTree/scripts"
source "${SCRIPT_DIR}/config.sh"

SUPERMATRIX="${DIR_CONCAT}/supermatrix.fa"
OUT_DIR="${DIR_GUIDE_TREE}"
PREFIX="${OUT_DIR}/guide"
THREADS=AUTO

# ─── Setup ────────────────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"

echo "================================================================"
echo " STEP 08: Guide species tree (LG+F+G4)"
echo "================================================================"
echo " Input      : ${SUPERMATRIX}"
echo " Output     : ${OUT_DIR}"
echo " Threads    : ${THREADS}"
echo " Started    : $(date)"
echo " Node       : $(hostname)"

# Verify input
if [[ ! -f "${SUPERMATRIX}" ]]; then
    echo "ERROR: Supermatrix not found: ${SUPERMATRIX}"
    echo "       Run python 07_concatenate.py first."
    exit 1
fi

N_TAXA=$(grep -c "^>" "${SUPERMATRIX}")
ALN_LEN=$(python3 -c "
from Bio import SeqIO
r = next(SeqIO.parse('${SUPERMATRIX}', 'fasta'))
print(len(r.seq))
")
echo " Taxa       : ${N_TAXA}"
echo " Aln length : ${ALN_LEN} AA columns"
echo ""

# ─── Run IQTree — Guide Tree ───────────────────────────────────────────────────
# -s          : supermatrix alignment
# -m LG+F+G4  : model (faster than C60, sufficient for guide tree)
# -B 1000     : ultrafast bootstrap
# -T 16       : CPU threads
# --prefix    : output file prefix
# -redo       : allow re-run if partial output exists
# NOTE: use 'iqtree' (v3.0.1), not 'iqtree2'
# NOTE: -p (partition) intentionally NOT used — all 3 reference papers avoid
#       partitioned models for the final species tree

echo "[08] Running IQTree guide tree (LG+F+G4)..."

iqtree \
    -s "${SUPERMATRIX}" \
    -m LG+F+G4 \
    -T "${THREADS}" \
    --prefix "${PREFIX}" \
    -redo

# Verify output
if [[ ! -f "${PREFIX}.treefile" ]]; then
    echo "ERROR: IQTree did not produce a treefile."
    exit 1
fi

echo ""
echo "[08] Guide tree complete: ${PREFIX}.treefile"
echo " Finished: $(date)"

# ─── Quick topology summary ────────────────────────────────────────────────────
echo ""
echo "--- Tree statistics ---"
python3 -c "
from Bio import Phylo
import numpy as np, io
t = Phylo.read('${PREFIX}.treefile', 'newick')
leaves = t.get_terminals()
print(f'Leaves (taxa)  : {len(leaves)}')
supports = [c.confidence for c in t.find_clades()
            if not c.is_terminal() and c.confidence and c.confidence > 0]
if supports:
    print(f'Bootstrap mean : {np.mean(supports):.1f}')
    print(f'Bootstrap ≥95  : {sum(1 for s in supports if s >= 95)}/{len(supports)} nodes ({100*sum(1 for s in supports if s >= 95)/len(supports):.0f}%)')
" 2>/dev/null || echo "   (tree stats unavailable)"

echo ""
echo "[08] NEXT: Compare guide tree to OrthoFinder tree:"
echo "     diff <(cat ${ORTHO_SPECIES_TREE}) <(cat ${PREFIX}.treefile)"
echo ""
echo "[08] Then submit: sbatch 09_species_tree_final.sh"
