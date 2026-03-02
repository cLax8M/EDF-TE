#!/bin/bash
# Full TE library curation pipeline — run from Library_curation/
# Usage: bash run_pipeline.sh [--from STEP] [--dry-run]
set -e

SCRIPTS="$(dirname "$0")/scripts"
FROM=1
DRY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

run_step() {
    local n=$1; local cmd=$2
    [[ $n -lt $FROM ]] && return
    echo "════════════════════════════════════════"
    echo "STEP $n: $cmd"
    echo "════════════════════════════════════════"
    [[ $DRY -eq 1 ]] && return
    bash "$SCRIPTS/$cmd" 2>&1 | tee "${SCRIPTS}/../logs/step${n}.log"
}

mkdir -p "$(dirname "$0")/logs"

run_step 1  01_combine_consensi.sh
run_step 2  02_mmseqs_cluster.sh
run_step 3  03_transeq.sh
run_step 4  04_diamond_refseq.sh
run_step 5  05_hmmer_te.sh
run_step 6  06_blast_repeatpeps.sh
run_step 7  07_host_vs_te.py   # called as: python3 scripts/07_host_vs_te.py
run_step 8  08_filter_host_fasta.sh
run_step 9  09_quality_filter.sh
run_step 10 10_evidence_tags.py
run_step 11 11_rename_final.py

echo "Pipeline complete → HoloTE_nonfungal_v1.fa"
