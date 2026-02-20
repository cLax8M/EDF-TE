#!/usr/bin/env bash
# =============================================================================
# 09_species_tree_final.sh — Final species tree with LG+C60+F+G4 (SLURM)
#
# PURPOSE : Infers the definitive maximum-likelihood species tree using the
#           LG+C60+F+G4 site-heterogeneous mixture model. This is the model
#           used in all three reference publications and is the state of the
#           art for phylogenomics with highly divergent eukaryotic taxa.
#
# WHY LG+C60+F+G4?
#   - C60: 60 frequency vector mixture classes capturing among-site variation
#          in amino acid composition (Lartillot & Philippe 2004)
#   - Substantially better fit than LG+G4 for organisms with >500 Ma separation
#   - Reduces long-branch attraction (LBA) artifacts, critical for Microsporidia,
#     Basidiomycetes, and other fast-evolving lineages in this dataset
#   - Used by: Liu et al. 2024 (348 taxa), Szánthó et al. 2025 (110 taxa Fungi),
#             Merényi et al. 2023 (fungal genomics)
#
# STARTING TREE : Guide tree from Step 08 (LG+F+G4) passed via -t flag
#   - Reduces run time by ~30-50% by providing a good topological starting point
#   - Standard practice in large-scale phylogenomics
#
# NO PARTITIONING : All three reference papers use unpartitioned concatenation
#   - Partitioned models require per-gene model selection — computationally
#     prohibitive and not shown to improve accuracy for C60 models
#
# EXPECTED TIME : ~24-96 hours for 238 taxa (depends on cluster speed and
#                 number of genes in supermatrix)
#
# USAGE   : sbatch 09_species_tree_final.sh
#
# OUTPUT  :
#   09_final_tree/final.treefile   — Final ML species tree (annotated with UFBoot)
#   09_final_tree/final.iqtree     — Full IQTree log (model parameters, likelihood)
#   09_final_tree/final.log        — Run log
#
# POST-PROCESSING:
#   - Root tree using outgroup (Holozoa: HO_Pluxyl, Drome, Homosap, etc.)
#   - Visualize with FigTree (local) or iTOL (web)
#   - Compare bootstrap support distribution
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# REFS    : Liu et al. 2024 (Mol Biol Evol)
#           Szánthó et al. 2025 (Nat. Commun.)
#           Merényi et al. 2023 (Nat. Ecol. Evol.)
#           Lartillot & Philippe 2004 (Mol Biol Evol) — C60 model
# =============================================================================

#SBATCH --job-name=iqtree_final
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=96:00:00
#SBATCH --output=/data/carlos/EDF-TE/SpeciesTree/logs/09_final_%j.out
#SBATCH --error=/data/carlos/EDF-TE/SpeciesTree/logs/09_final_%j.err

set -euo pipefail

# ─── Configuration ────────────────────────────────────────────────────────────
SCRIPT_DIR="/data/carlos/EDF-TE/EDF-TE/SpeciesTree/scripts"
source "${SCRIPT_DIR}/config.sh"

SUPERMATRIX="${DIR_CONCAT}/supermatrix.fa"
GUIDE_TREE="${DIR_GUIDE_TREE}/guide.treefile"
OUT_DIR="${DIR_FINAL_TREE}"
PREFIX="${OUT_DIR}/final"
THREADS=32

# ─── Setup ────────────────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"

echo "================================================================"
echo " STEP 09: Final species tree (LG+C60+F+G4)"
echo "================================================================"
echo " Input      : ${SUPERMATRIX}"
echo " Guide tree : ${GUIDE_TREE}"
echo " Output     : ${OUT_DIR}"
echo " Threads    : ${THREADS}"
echo " Started    : $(date)"
echo " Node       : $(hostname)"

# Verify inputs
if [[ ! -f "${SUPERMATRIX}" ]]; then
    echo "ERROR: Supermatrix not found: ${SUPERMATRIX}"
    exit 1
fi

if [[ ! -f "${GUIDE_TREE}" ]]; then
    echo "WARNING: Guide tree not found: ${GUIDE_TREE}"
    echo "         Running without starting tree (-t). This will be slower."
    USE_GUIDE=false
else
    USE_GUIDE=true
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

# ─── Run IQTree — Final Tree ───────────────────────────────────────────────────
# -s              : supermatrix alignment
# -m LG+C60+F+G4  : site-heterogeneous mixture model (state of the art)
# -B 1000         : ultrafast bootstrap replicates
# -T 32           : CPU threads
# -t              : starting tree from Step 08 (reduces computation time)
# --prefix        : output file prefix
# -redo           : allow restart if prior partial run exists
#
# IMPORTANT NOTES:
# - LG+C60 requires substantially more memory than LG+G4 (~4-8x)
# - Run time scales with: n_taxa² × alignment_length × n_C60_classes
# - If memory is insufficient, reduce --mem and --cpus-per-task
# - C60 uses 60 mixture components; C20 is faster if time is critical

echo "[09] Running IQTree final tree (LG+C60+F+G4)..."
echo ""

if [[ "${USE_GUIDE}" == "true" ]]; then
    echo "     Using guide tree as starting topology: ${GUIDE_TREE}"
    iqtree \
        -s "${SUPERMATRIX}" \
        -m LG+C60+F+G4 \
        -B 1000 \
        -T "${THREADS}" \
        -t "${GUIDE_TREE}" \
        --prefix "${PREFIX}" \
        -redo
else
    echo "     No guide tree — running from random starting tree"
    iqtree \
        -s "${SUPERMATRIX}" \
        -m LG+C60+F+G4 \
        -B 1000 \
        -T "${THREADS}" \
        --prefix "${PREFIX}" \
        -redo
fi

# Verify output
if [[ ! -f "${PREFIX}.treefile" ]]; then
    echo "ERROR: IQTree did not produce a treefile."
    exit 1
fi

echo ""
echo "[09] ================================================================"
echo "[09] FINAL SPECIES TREE COMPLETE"
echo "[09] ================================================================"
echo " Tree file  : ${PREFIX}.treefile"
echo " Log file   : ${PREFIX}.iqtree"
echo " Finished   : $(date)"

# ─── Tree statistics ───────────────────────────────────────────────────────────
echo ""
echo "--- Final tree statistics ---"
python3 -c "
from ete3 import Tree
import numpy as np
t = Tree('${PREFIX}.treefile')
leaves = t.get_leaves()
supports = [n.support for n in t.traverse() if not n.is_leaf() and n.support > 0]
print(f'Leaves (taxa)  : {len(leaves)}')
if supports:
    print(f'Bootstrap mean : {np.mean(supports):.1f}')
    print(f'Bootstrap ≥95  : {sum(1 for s in supports if s >= 95)}/{len(supports)} ({100*sum(1 for s in supports if s >= 95)/len(supports):.0f}%)')
    print(f'Bootstrap ≥80  : {sum(1 for s in supports if s >= 80)}/{len(supports)} ({100*sum(1 for s in supports if s >= 80)/len(supports):.0f}%)')
"

# ─── Post-processing guidance ──────────────────────────────────────────────────
echo ""
echo "--- Post-processing steps ---"
echo ""
echo "1. ROOT THE TREE (using Holozoa outgroup):"
echo "   python3 -c \""
echo "   from ete3 import Tree"
echo "   t = Tree('${PREFIX}.treefile')"
echo "   # Outgroup: Holozoa species (HO_ prefix or Choanoflagellates)"
echo "   outgroup = [l.name for l in t.get_leaves() if l.name.startswith('HO_')]"
echo "   t.set_outgroup(t.get_common_ancestor(outgroup))"
echo "   t.write(outfile='${OUT_DIR}/final.rooted.treefile')"
echo "   \""
echo ""
echo "2. VISUALIZE: Open in FigTree or upload to iTOL (itol.embl.de)"
echo ""
echo "3. COPY TO RESULTS:"
echo "   cp -r ${WORK_DIR}/09_final_tree/ /legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/"
echo ""
echo "4. LOG LIKELIHOOD from IQTree log:"
grep "Log-likelihood" "${PREFIX}.iqtree" 2>/dev/null | tail -1 || true
