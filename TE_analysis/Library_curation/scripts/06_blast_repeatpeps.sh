#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"

mkdir -p "${WORK_DIR}/06_repeatpeps"

module use /data/toby/tools/modulefiles && module load earlgrey

blastp \
    -query "${WORK_DIR}/03_translation/consensus_nr_6f.faa" \
    -db "${REPEATPEPS}" \
    -out "${WORK_DIR}/06_repeatpeps/repeatpeps_hits.tsv" \
    -evalue 1e-3 \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
    -max_target_seqs 1 \
    -num_threads ${THREADS}

echo "RepeatPeps hits: $(wc -l < "${WORK_DIR}/06_repeatpeps/repeatpeps_hits.tsv")"
