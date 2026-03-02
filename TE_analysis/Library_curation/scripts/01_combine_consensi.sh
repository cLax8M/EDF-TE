#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"

mkdir -p "${WORK_DIR}/01_combined"

echo "Collecting strained consensi from all EarlGrey runs..."
for f in "${EARlgrey_OUT}"/*_EarlGrey/*_summaryFiles/*.strained; do
    [[ -f "$f" ]] || continue
    sp=$(basename "$(dirname "$(dirname "$f")")" _EarlGrey)
    # Prefix each sequence ID with species name to ensure uniqueness
    awk -v sp="$sp" '/^>/{sub(/^>/,">'"'"'${sp}'"'"'__")} 1' "$f"
done > "${WORK_DIR}/01_combined/all_consensus.fa"

echo "Stats:"
module use /data/toby/tools/modulefiles && module load earlgrey
seqkit stats "${WORK_DIR}/01_combined/all_consensus.fa" | tee "${WORK_DIR}/01_combined/stats.txt"
