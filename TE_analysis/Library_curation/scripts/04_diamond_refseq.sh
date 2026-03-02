#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"

mkdir -p "${WORK_DIR}/04_host_filter"

module use /data/toby/tools/modulefiles && module load earlgrey

diamond blastp \
    -q "${WORK_DIR}/03_translation/consensus_nr_6f.faa" \
    -d "${REFSEQ_DB}" \
    -o "${WORK_DIR}/04_host_filter/refseq_hits.tsv" \
    --sensitive \
    --matrix BLOSUM62 \
    -e 1e-3 \
    --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --threads ${THREADS}

echo "RefSeq hits: $(wc -l < "${WORK_DIR}/04_host_filter/refseq_hits.tsv")"
