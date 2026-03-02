#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"

mkdir -p "${WORK_DIR}/03_translation"

module use /data/toby/tools/modulefiles && module load earlgrey

transeq \
    -sequence "${WORK_DIR}/02_clustering/consensus_nr.fa" \
    -outseq "${WORK_DIR}/03_translation/consensus_nr_6f.faa" \
    -frame 6 \
    -clean \
    -auto
