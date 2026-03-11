#!/bin/bash
# Combine HoloTE_nonfungal + MycoMobilome into HoloTE_combined_v1.fa
# Run this after HoloTE_nonfungal pipeline (job 2540991) completes

HOLOTE="/data/carlos/EDF-TE/TE_analysis/Library_curation/HoloTE_nonfungal_v1.fa"
MYCO="/legserv/Temp/Carlos/EDF-TE/results/TE_pipeline/consensus_myco/MycoMobilome_v1.1-clustered_80_TE_library.fasta"
OUT="/data/carlos/EDF-TE/TE_analysis/Library_curation/HoloTE_combined_v1.fa"

echo "HoloTE_nonfungal: $(grep -c '^>' "$HOLOTE") seqs"
echo "MycoMobilome:     $(grep -c '^>' "$MYCO") seqs"

cat "$HOLOTE" "$MYCO" > "$OUT"

echo "Combined:         $(grep -c '^>' "$OUT") seqs → $OUT"
echo "Size:             $(du -sh "$OUT" | cut -f1)"
