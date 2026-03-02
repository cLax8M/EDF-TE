#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"

mkdir -p "${WORK_DIR}/08_te_only"

module use /data/toby/tools/modulefiles && module load earlgrey

seqkit grep \
    -f "${WORK_DIR}/07_decision/keep_ids.txt" \
    "${WORK_DIR}/02_clustering/consensus_nr.fa" \
    > "${WORK_DIR}/08_te_only/consensus_te_only.fa"

echo "After host filter:"
seqkit stats "${WORK_DIR}/08_te_only/consensus_te_only.fa"
