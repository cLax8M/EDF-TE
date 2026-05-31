#!/bin/bash
# 05d_merge_genome_results.sh
# ─────────────────────────────────────────────────────────────────────────────
# Run AFTER SLURM array job 05d_identify_gremlins_genome.sh completes.
# Merges per-species *_families.tsv files into gremlin_families_genome.tsv.
#
# Usage:
#   bash 05d_merge_genome_results.sh
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

GENOME_RESULTS="/legserv/Temp/Carlos/EDF-TE/data/gremlin/genome_tblastn_results"
OUT="${GENOME_RESULTS}/gremlin_families_genome.tsv"

echo -e "species\tfamily_name\tnum_hits" > "$OUT"

n_species=0
for F in "$GENOME_RESULTS"/*_families.tsv; do
    [[ -s "$F" ]] || continue
    SPECIES=$(basename "$F" _families.tsv)
    # Each line: "  N  FAMILY_NAME" (uniq -c output)
    awk -v sp="$SPECIES" '{
        hits = $1
        fam  = $2
        # Lower-case normalise (RND-1_FAMILY-N → rnd-1_family-N)
        gsub(/^RND-/, "rnd-", fam)
        # Strip any trailing _N instance index that may have slipped through
        gsub(/_[0-9]+$/, "", fam)
        print sp "\t" fam "\t" hits
    }' "$F" >> "$OUT"
    n_species=$((n_species + 1))
done

echo "Merged $n_species species into $OUT"
echo "Lines: $(wc -l < "$OUT")"

# Quick summary
echo ""
echo "Species with genome tBLASTn Grem-LINE1 hits:"
awk 'NR>1{print $1}' "$OUT" | sort -u

echo ""
echo "Now you can run 05_gremlin_centromere.py — it will automatically pick up"
echo "both ${GENOME_RESULTS}/gremlin_families_genome.tsv and the existing"
echo "tblastn_results/gremlin_families_tblastn_denovo.tsv"
