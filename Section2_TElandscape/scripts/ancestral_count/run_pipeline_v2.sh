#!/usr/bin/env bash
#SBATCH --job-name=ancestral_recon_v2
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --exclude=LEGcompute1
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count/logs/pipeline_v2_%j.out
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count/logs/pipeline_v2_%j.err

set -euo pipefail

# Paths
SCRIPTS=/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts/ancestral_count
ANNOT=/data/carlos/EDF-TE/TE_analysis/annotated
RESULTS=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count
TREE=${RESULTS}/05_count/guide_clean.treefile
PY=/data/carlos/EDF-TE/.venv/bin/python3

# Output dirs (v2 suffix to avoid clobbering old results)
FL_DIR=${RESULTS}/01_fl_families_v2
CONS_FA=${RESULTS}/02_fl_consensus_v2.fa
CLUST_PREFIX=${RESULTS}/03_clusters_v2/cluster
PA_DIR=${RESULTS}/04_pa_matrix_v2
COUNT_DIR=${RESULTS}/05_count_v2

mkdir -p "$FL_DIR" "$(dirname "$CLUST_PREFIX")" "$PA_DIR" "$COUNT_DIR" \
         "${RESULTS}/logs"

echo "========== Step 01: Extract FL copies (case-insensitive) =========="
$PY ${SCRIPTS}/01_extract_fl_copies.py \
    --annotated_dir "$ANNOT" \
    --out_dir "$FL_DIR"

# Quick stats
TOTAL_SP=$(ls "$FL_DIR"/*.tsv 2>/dev/null | wc -l)
WITH_FL=$(for f in "$FL_DIR"/*.tsv; do
    awk -F'\t' 'NR>1 && $4>=1 {found=1; exit} END{if(found) print 1}' "$f"
done | wc -l)
echo "  Species: $TOTAL_SP total, $WITH_FL with ≥1 FL family"

echo ""
echo "========== Step 02: Collect FL consensus =========="
$PY ${SCRIPTS}/02_collect_consensus.py \
    --fl_dir "$FL_DIR" \
    --annotated_dir "$ANNOT" \
    --out_fasta "$CONS_FA"
echo "  Sequences: $(grep -c '^>' "$CONS_FA")"

echo ""
echo "========== Step 03: MMseqs2 clustering =========="
bash ${SCRIPTS}/03_cluster_consensus.sh "$CONS_FA" "$CLUST_PREFIX" 0.80 0.80
echo "  Clusters: $(wc -l < "${CLUST_PREFIX}_cluster.tsv")"

echo ""
echo "========== Step 04: Build PA matrix (all 238 tree tips) =========="
$PY ${SCRIPTS}/04_build_pa_matrix.py \
    --cluster_tsv "${CLUST_PREFIX}_cluster.tsv" \
    --fl_dir "$FL_DIR" \
    --tree "$TREE" \
    --out_dir "$PA_DIR"

echo ""
echo "========== Step 05: Count Wagner parsimony =========="
bash ${SCRIPTS}/05_run_count.sh \
    "${PA_DIR}/te_family_pa_matrix.count" \
    "$TREE" \
    "$COUNT_DIR"

echo ""
echo "========== Pipeline complete =========="
echo "FL dir:     $FL_DIR"
echo "Consensus:  $CONS_FA"
echo "Clusters:   ${CLUST_PREFIX}_cluster.tsv"
echo "PA matrix:  ${PA_DIR}/te_family_pa_matrix.count"
echo "Wagner out: ${COUNT_DIR}/wagner_ancestral.txt"
cat "${PA_DIR}/summary.txt"
