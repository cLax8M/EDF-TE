#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"

mkdir -p "${WORK_DIR}/02_clustering" "${WORK_DIR}/tmp/mmseqs"

module use /data/toby/tools/modulefiles && module load earlgrey

mmseqs createdb "${WORK_DIR}/01_combined/all_consensus.fa" "${WORK_DIR}/tmp/mmseqs/db"
mmseqs cluster "${WORK_DIR}/tmp/mmseqs/db" "${WORK_DIR}/tmp/mmseqs/clu" "${WORK_DIR}/tmp/mmseqs" \
    --min-seq-id 0.8 -c 0.8 --cov-mode 1 --cluster-reassign \
    --threads ${THREADS}
mmseqs createseqfiledb "${WORK_DIR}/tmp/mmseqs/db" "${WORK_DIR}/tmp/mmseqs/clu" \
    "${WORK_DIR}/tmp/mmseqs/rep"
mmseqs result2flat "${WORK_DIR}/tmp/mmseqs/db" "${WORK_DIR}/tmp/mmseqs/db" \
    "${WORK_DIR}/tmp/mmseqs/rep" "${WORK_DIR}/02_clustering/consensus_nr.fa"

echo "Clustering complete:"
seqkit stats "${WORK_DIR}/02_clustering/consensus_nr.fa"
