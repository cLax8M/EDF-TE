#!/usr/bin/env bash
#SBATCH --job-name=ancestral_v2_resume
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --exclude=LEGcompute1
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count/logs/pipeline_v2_resume_%j.out
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count/logs/pipeline_v2_resume_%j.err

# Resume pipeline v2 from step 03 (steps 01-02 completed in job 4188751)
set -euo pipefail

SCRIPTS=/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts/ancestral_count
RESULTS=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count
TREE=${RESULTS}/05_count/guide_clean.treefile
PY=/data/carlos/EDF-TE/.venv/bin/python3

# Inputs from completed steps
FL_DIR=${RESULTS}/01_fl_families_v2
CONS_FA=${RESULTS}/02_fl_consensus_v2.fa

# Outputs for steps 03-05
CLUST_PREFIX=${RESULTS}/03_clusters_v2/cluster
PA_DIR=${RESULTS}/04_pa_matrix_v2
COUNT_DIR=${RESULTS}/05_count_v2

mkdir -p "$(dirname "$CLUST_PREFIX")" "$PA_DIR" "$COUNT_DIR"

echo "========== Step 03: MMseqs2 clustering (v17) =========="
echo "Input: $CONS_FA ($(grep -c '^>' "$CONS_FA") sequences)"
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
echo "========== Pipeline v2 complete =========="
cat "${PA_DIR}/summary.txt"
echo "Wagner output: ${COUNT_DIR}/wagner_ancestral.txt ($(wc -l < "${COUNT_DIR}/wagner_ancestral.txt") lines)"
