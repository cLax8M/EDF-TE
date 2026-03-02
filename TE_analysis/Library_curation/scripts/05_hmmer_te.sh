#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"

mkdir -p "${WORK_DIR}/05_te_domains"

module use /data/toby/tools/modulefiles && module load earlgrey

hmmscan \
    --cpu ${THREADS} \
    -E 10 \
    --noali \
    --tblout "${WORK_DIR}/05_te_domains/te_domains.tbl" \
    "${REPET_HMM}" \
    "${WORK_DIR}/03_translation/consensus_nr_6f.faa"

echo "HMMER hits (raw): $(grep -v '^#' "${WORK_DIR}/05_te_domains/te_domains.tbl" | wc -l)"
