#!/usr/bin/env bash
# 03_cluster_consensus.sh — Cross-species TE consensus clustering with MMseqs2.
#
# Usage:
#   bash 03_cluster_consensus.sh <input.fa> <out_prefix> [min_seq_id] [coverage]
#
# Defaults: min_seq_id=0.80, coverage=0.80, cov-mode 0 (bidirectional)
# Produces: <out_prefix>_cluster.tsv  (rep_id<TAB>member_id)
#           <out_prefix>_rep_seq.fasta
# Required: mmseqs in PATH or in /data/toby/tools/

set -euo pipefail

INPUT_FA="${1:?Usage: 03_cluster_consensus.sh <input.fa> <out_prefix> [min_id] [cov]}"
OUT_PREFIX="${2:?out_prefix required}"
MIN_ID="${3:-0.80}"
COVERAGE="${4:-0.80}"

# Locate mmseqs
if ! command -v mmseqs >/dev/null 2>&1; then
    if [[ -x /opt/easybuild/software/MMseqs2/17-b804f-gompi-2024a/bin/mmseqs ]]; then
        export PATH="/opt/easybuild/software/MMseqs2/17-b804f-gompi-2024a/bin:$PATH"
    elif [[ -x /opt/easybuild/software/MMseqs2/15-6f452-gompi-2023b/bin/mmseqs ]]; then
        export PATH="/opt/easybuild/software/MMseqs2/15-6f452-gompi-2023b/bin:$PATH"
    elif [[ -x /data/toby/tools/mmseqs/bin/mmseqs ]]; then
        export PATH="/data/toby/tools/mmseqs/bin:$PATH"
    else
        echo "ERROR: mmseqs not found. Install or set PATH." >&2
        exit 1
    fi
fi

OUT_DIR=$(dirname "$OUT_PREFIX")
TMPDIR_MM="${OUT_DIR}/tmp_mmseqs_$$"
mkdir -p "$OUT_DIR" "$TMPDIR_MM"

echo "Input:    $INPUT_FA  ($(grep -c '^>' "$INPUT_FA") sequences)"
echo "Out:      ${OUT_PREFIX}_cluster.tsv"
echo "Identity: $MIN_ID"
echo "Coverage: $COVERAGE"

# easy-cluster: nucleotide-aware; --cov-mode 0 = bidirectional coverage
mmseqs easy-cluster \
    "$INPUT_FA" \
    "$OUT_PREFIX" \
    "$TMPDIR_MM" \
    --min-seq-id "$MIN_ID" \
    -c "$COVERAGE" \
    --cov-mode 0 \
    --threads 8 \
    --remove-tmp-files 1

# Cleanup
rm -rf "$TMPDIR_MM"

NCL=$(awk '{print $1}' "${OUT_PREFIX}_cluster.tsv" | sort -u | wc -l)
NMEM=$(wc -l < "${OUT_PREFIX}_cluster.tsv")
echo "Done: $NCL clusters | $NMEM members | TSV: ${OUT_PREFIX}_cluster.tsv"
