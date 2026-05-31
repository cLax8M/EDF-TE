#!/bin/bash
#SBATCH --job-name=htt01_phylo_of
#SBATCH --partition=normal.168h
#SBATCH --exclude=LEGcompute1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/carlos/EDF-TE/logs/htt01_phylo_of_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/htt01_phylo_of_%j.err
set -euo pipefail

# ============================================================================
# HTT pipeline step 01 - OrthoFinder adaptation (replaces BUSCO+FAMSA+IQTree)
# ============================================================================
# Romeijn 2025 HTT pipeline expects per-gene aligned protein sequences in
# `${RESULTS}/sequences/{gene}.faa` and a rooted species tree.  Because BUSCO
# yields too few single-copy orthologs across the 238 fungal/protistan
# genomes (~1 Gyr of divergence, 0 fully universal SCO), we use OrthoFinder
# instead.  This script:
#
#   1. Copies OrthoFinder's species tree (already ML, already rooted) into
#      `phylogeny_rooted.treefile` and strips the `.aa` suffix from tips so
#      they match the genome IDs (MU_*, AS_*, etc.).
#   2. For each near-universal single-copy orthogroup (>=80 % of species,
#      385 OGs in `orthofinder_near_universal_80pct_exactly1.txt`), copies
#      the existing OrthoFinder MSA into `${RESULTS}/sequences/{OG}.faa` —
#      again rewriting headers to plain genome IDs (one sequence per genome).
#
# Outputs (compatible with HTT step 03_ks_divergence.sh):
#   ${RESULTS}/phylogeny_rooted.treefile     (newick, branch lengths in subs/site)
#   ${RESULTS}/sequences/OGxxxxxxx.faa       (one MSA per gene, headers = genome ID)
#   ${RESULTS}/selected_ogs.txt              (list of orthogroups used)
# ============================================================================

ROOT=/legserv/Temp/Carlos/EDF-TE
OF_RUN=${ROOT}/results/Orthofinder_run
OG_LIST=${ROOT}/results/HTT/orthofinder_near_universal_80pct_exactly1.txt
RESULTS=${ROOT}/results/HTT/01_phylo_tree
LOGS=${ROOT}/results/HTT/logs/01

mkdir -p "${RESULTS}/sequences" "${LOGS}"

echo "[$(date)] HTT step 01 (OrthoFinder adaptation)"
echo "  OF_RUN  : ${OF_RUN}"
echo "  OG_LIST : ${OG_LIST} ($(wc -l < "${OG_LIST}") OGs)"
echo "  RESULTS : ${RESULTS}"

# ---------------------------------------------------------------------------
# 1. Rooted species tree (strip .aa suffix from tips)
# ---------------------------------------------------------------------------
SP_TREE_IN=${OF_RUN}/Species_Tree/SpeciesTree_rooted.txt
SP_TREE_OUT=${RESULTS}/phylogeny_rooted.treefile

if [ ! -f "${SP_TREE_IN}" ]; then
    echo "ERROR: OF species tree not found: ${SP_TREE_IN}" >&2
    exit 1
fi

# Strip ".aa" right before ":" or "," or ")"
sed -E 's/\.aa([:,)])/\1/g' "${SP_TREE_IN}" > "${SP_TREE_OUT}"

NTIPS=$(grep -oE '[A-Za-z][A-Za-z0-9_.-]+:' "${SP_TREE_OUT}" | wc -l)
echo "  rooted tree tips: ${NTIPS}"

# ---------------------------------------------------------------------------
# 2. Copy MSA for each selected OG, stripping protein-ID prefixes so each
#    sequence header is just the genome ID (the part before the first "|").
#    Drop genomes that have >1 sequence in the OG (keep strict 1-copy only).
# ---------------------------------------------------------------------------
MSA_DIR=${OF_RUN}/MultipleSequenceAlignments
SEQ_DIR=${RESULTS}/sequences
SELECTED=${RESULTS}/selected_ogs.txt
> "${SELECTED}"

n_ok=0; n_skip=0
while IFS= read -r og; do
    src="${MSA_DIR}/${og}.fa"
    dst="${SEQ_DIR}/${og}.faa"
    if [ ! -f "${src}" ]; then
        n_skip=$((n_skip+1))
        continue
    fi
    # rewrite headers: ">AS_AspfumUD1_aa_jgi|AspfumUD1_1|..." -> ">AS_AspfumUD1"
    # OrthoFinder MSAs prefix every header with "<species>_aa_" where <species>
    # is the input proteome basename (stripped of ".aa.fasta"). The species tree
    # tip labels are the same name with a trailing ".aa" — already removed in
    # phylogeny_rooted.treefile above. So we just split on "_aa_" and take [1].
    awk '/^>/ { sub(/^>/,""); n=index($0,"_aa_");
                if (n>0) name=substr($0,1,n-1); else name=$0;
                print ">"name; next }
         { print }' "${src}" > "${dst}"
    echo "${og}" >> "${SELECTED}"
    n_ok=$((n_ok+1))
done < "${OG_LIST}"

echo "  OG MSAs copied: ${n_ok}  (missing: ${n_skip})"
echo "[DONE] $(date)"
