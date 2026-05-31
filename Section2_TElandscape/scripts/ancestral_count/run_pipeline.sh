#!/usr/bin/env bash
# run_pipeline.sh — Orchestrate the FASE 3 ancestral reconstruction pipeline.
#   01: extract full-length families per species (local Python)
#   02: collect FL consensus sequences across all species into one FASTA
#   03: MMseqs2 cluster (SLURM job)
#   04: build presence/absence matrix
#   05: run Count v9.1 AsymmetricWagner
#
# Steps 01-02-04 run locally (lightweight). Step 03 submitted to SLURM.
# Step 05 runs locally after step 04.
#
# Usage: bash run_pipeline.sh [step]
#   step = 01 | 02 | 03 | 04 | 05 | all
# Default: all

set -euo pipefail

# ── paths ────────────────────────────────────────────────────────────────────
ANNOT_DIR=/data/carlos/EDF-TE/TE_analysis/annotated
SCRIPTS_DIR=/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts/ancestral_count
DATA_DIR=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count
TREE=/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count/05_count/guide_clean.treefile
WHITELIST=/data/carlos/EDF-TE/EDF-TE/Synteny/data/qc/species_whitelist.txt

FL_DIR=$DATA_DIR/01_fl_families
FL_FA=$DATA_DIR/02_fl_consensus.fa
CL_PREFIX=$DATA_DIR/03_clusters/cluster
PA_DIR=$DATA_DIR/04_pa_matrix
COUNT_DIR=$DATA_DIR/05_count

mkdir -p "$FL_DIR" "$(dirname "$CL_PREFIX")" "$PA_DIR" "$COUNT_DIR"

STEP="${1:-all}"

run_01() {
    echo "── 01: extract FL families per species ──"
    python3 "$SCRIPTS_DIR/01_extract_fl_copies.py" \
        --annotated_dir "$ANNOT_DIR" \
        --out_dir "$FL_DIR" \
        --fl_frac 0.8
}

run_02() {
    echo "── 02: collect FL consensus across species ──"
    python3 "$SCRIPTS_DIR/02_collect_consensus.py" \
        --fl_dir "$FL_DIR" \
        --annotated_dir "$ANNOT_DIR" \
        --out_fasta "$FL_FA"
}

run_03() {
    echo "── 03: MMseqs2 clustering (SLURM) ──"
    JOB_SCRIPT=$DATA_DIR/03_clusters/submit_mmseqs.sh
    cat > "$JOB_SCRIPT" <<SBATCH
#!/bin/bash
#SBATCH --job-name=ancCount03
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=$DATA_DIR/03_clusters/mmseqs_%j.log
#SBATCH --exclude=LEGcompute1

set -euo pipefail
export PATH=/opt/easybuild/software/MMseqs2/18-8cc5c-gompi-2025a/bin:$PATH

mkdir -p $(dirname $CL_PREFIX)
TMPDIR_MM=$DATA_DIR/03_clusters/tmp_\$SLURM_JOB_ID
mkdir -p \$TMPDIR_MM

mmseqs easy-cluster \\
    $FL_FA \\
    $CL_PREFIX \\
    \$TMPDIR_MM \\
    --min-seq-id 0.80 \\
    -c 0.80 \\
    --cov-mode 0 \\
    --threads 16 \\
    --remove-tmp-files 1

rm -rf \$TMPDIR_MM
NCL=\$(awk '{print \$1}' ${CL_PREFIX}_cluster.tsv | sort -u | wc -l)
NMEM=\$(wc -l < ${CL_PREFIX}_cluster.tsv)
echo "Done: \$NCL clusters | \$NMEM members"
SBATCH
    JID=$(sbatch "$JOB_SCRIPT" | awk '{print $4}')
    echo "Submitted MMseqs cluster job: $JID"
    echo "$JID" > "$DATA_DIR/03_clusters/last_jobid.txt"
    echo "When complete: bash run_pipeline.sh 04"
}

run_04() {
    echo "── 04: build PA matrix ──"
    python3 "$SCRIPTS_DIR/04_build_pa_matrix.py" \
        --cluster_tsv "${CL_PREFIX}_cluster.tsv" \
        --fl_dir "$FL_DIR" \
        --whitelist "$WHITELIST" \
        --out_dir "$PA_DIR"
}

run_05() {
    echo "── 05: run Count Wagner parsimony ──"
    bash "$SCRIPTS_DIR/05_run_count.sh" \
        "$PA_DIR/te_family_pa_matrix.count" \
        "$TREE" \
        "$COUNT_DIR" \
        2.0
}

case "$STEP" in
    01)   run_01 ;;
    02)   run_02 ;;
    03)   run_03 ;;
    04)   run_04 ;;
    05)   run_05 ;;
    all)  run_01; run_02; run_03
          echo "Step 03 submitted asynchronously. Run '04' and '05' once cluster job is done." ;;
    *)    echo "Unknown step: $STEP"; exit 1 ;;
esac
