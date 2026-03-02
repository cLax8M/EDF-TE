#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"

mkdir -p "${WORK_DIR}/09_quality"

module use /data/toby/tools/modulefiles && module load earlgrey

# Length filter
seqkit seq \
    -m ${MIN_LEN} \
    -M ${MAX_LEN} \
    "${WORK_DIR}/08_te_only/consensus_te_only.fa" \
    > "${WORK_DIR}/09_quality/len_filtered.fa"

echo "After length filter (${MIN_LEN}-${MAX_LEN}bp):"
seqkit stats "${WORK_DIR}/09_quality/len_filtered.fa"

# N-content filter (<=5%)
# seqtk comp gives: name len #A #C #G #T #2 #3 #4 CpG ts tv #ins #del
seqtk comp "${WORK_DIR}/09_quality/len_filtered.fa" | \
    awk '{
        len=$2; a=$3; c=$4; g=$5; t=$6;
        acgt=a+c+g+t;
        n=len-acgt;
        frac=(len>0)?n/len:1;
        if(frac <= '"${MAX_N_FRAC}"') print $1
    }' > "${WORK_DIR}/09_quality/pass_n_filter.txt"

seqkit grep \
    -f "${WORK_DIR}/09_quality/pass_n_filter.txt" \
    "${WORK_DIR}/09_quality/len_filtered.fa" \
    > "${WORK_DIR}/09_quality/consensus_quality_filtered.fa"

echo "After N-content filter (<=5% N):"
seqkit stats "${WORK_DIR}/09_quality/consensus_quality_filtered.fa"
