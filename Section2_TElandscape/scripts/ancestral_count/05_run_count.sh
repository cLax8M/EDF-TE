#!/usr/bin/env bash
# 05_run_count.sh — Run Count v9.1 AsymmetricWagner parsimony for TE family
#                    presence/absence ancestral reconstruction.
#
# Usage: bash 05_run_count.sh <pa_matrix.count> <tree.newick> <out_dir>
#
# Count's command-line tools (in Count.jar):
#   ca.umontreal.iro.evolution.genecontent.AsymmetricWagner
#       Args: <gain_penalty> <tree.newick> <table.txt>
#       Default gain_penalty=2.0 (loss=1) → typical for TE birth-death.
#
# Output:
#   {out_dir}/wagner_ancestral.txt    Family<TAB>ancestral states per node
#   {out_dir}/wagner_log.txt          Stdout + total events
#
# Prerequisites:
#   COUNT_JAR=/data/carlos/EDF-TE/tools/Count/Count/Count.jar
#   PA matrix in Count format (Family<TAB>sp1<TAB>sp2 ... with 0/1 values)
#   Newick tree with tip labels matching matrix column names

set -euo pipefail

PA_FILE="${1:?Usage: 05_run_count.sh <pa_matrix.count> <tree.newick> <out_dir>}"
TREE="${2:?Tree required}"
OUT_DIR="${3:?Out dir required}"
GAIN_PEN="${4:-2.0}"

COUNT_JAR="${COUNT_JAR:-/data/carlos/EDF-TE/tools/Count/Count/Count.jar}"

mkdir -p "$OUT_DIR"

if [[ ! -f "$COUNT_JAR" ]]; then
    echo "ERROR: Count.jar not found at $COUNT_JAR" >&2
    exit 1
fi
if [[ ! -f "$PA_FILE" ]]; then
    echo "ERROR: PA matrix not found: $PA_FILE" >&2
    exit 1
fi
if [[ ! -f "$TREE" ]]; then
    echo "ERROR: Tree not found: $TREE" >&2
    exit 1
fi

echo "Count.jar:    $COUNT_JAR"
echo "PA matrix:    $PA_FILE  ($(($(wc -l < "$PA_FILE") - 1)) families)"
echo "Tree:         $TREE"
echo "Gain penalty: $GAIN_PEN"
echo "Out dir:      $OUT_DIR"

OUT_WAGNER="$OUT_DIR/wagner_ancestral.txt"
LOG="$OUT_DIR/wagner_log.txt"

# Quote tip labels in the tree: Newick standard treats underscores as spaces,
# which breaks Count.jar's matching against the PA matrix columns.
QUOTED_TREE="$OUT_DIR/guide_quoted.treefile"
sed "s/\([A-Za-z][A-Za-z0-9_]*\):/'\1':/g" "$TREE" > "$QUOTED_TREE"
echo "Quoted tree:  $QUOTED_TREE"

# Run Asymmetric Wagner parsimony
java -Xmx32G -cp "$COUNT_JAR" \
    ca.umontreal.iro.evolution.genecontent.AsymmetricWagner \
    -gain "$GAIN_PEN" "$QUOTED_TREE" "$PA_FILE" \
    > "$OUT_WAGNER" 2> "$LOG"

# Quick summary
NL=$(wc -l < "$OUT_WAGNER")
echo "Done. Output: $OUT_WAGNER ($NL lines)"
echo "Log:  $LOG"
echo ""
echo "First 5 lines of output:"
head -n 5 "$OUT_WAGNER"
