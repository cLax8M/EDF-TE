#!/usr/bin/env python3
"""
03b_fix_alignments.py

Step 3 replacement: write protein alignment files directly from node_gene_pairs.json,
using the pre-aligned sequences from the OrthoFinder MSA.
MAFFT is not needed — the sequences in node_gene_pairs are already aligned slices
from the OrthoFinder MSA.

Removes all-gap columns so PAL2NAL gets clean pairwise protein alignments.
"""

import json
import sys
from pathlib import Path

RESULTS = Path("/legserv/Temp/Carlos/EDF-TE/results/HTT/03_ks_divergence")
PAIRS_FILE = RESULTS / "node_analyses" / "node_gene_pairs.json"
ALIGN_DIR  = RESULTS / "alignments"

ALIGN_DIR.mkdir(parents=True, exist_ok=True)

print(f"Loading pairs from {PAIRS_FILE}")
with open(PAIRS_FILE) as f:
    pairs = json.load(f)
print(f"Total pairs: {len(pairs)}")

written = 0
skipped = 0
for i, pair in enumerate(pairs):
    node_name = pair["node_name"]
    gene      = pair["gene"]
    reps      = pair["representatives"]   # list of [name, aligned_seq]

    if not reps:
        skipped += 1
        continue

    seqs = [(name, seq) for name, seq in reps]
    n_cols = len(seqs[0][1])

    # Verify all seqs have same length (they should, as MSA slices)
    if any(len(s) != n_cols for _, s in seqs):
        # Trim/pad to shortest for safety
        min_len = min(len(s) for _, s in seqs)
        seqs = [(n, s[:min_len]) for n, s in seqs]
        n_cols = min_len

    # Remove all-gap columns (columns where every representative has '-')
    keep = [c for c in range(n_cols) if any(s[c] != "-" for _, s in seqs)]

    if not keep:
        skipped += 1
        continue

    out_file = ALIGN_DIR / f"{node_name}_{gene}_aligned.faa"
    with open(out_file, "w") as fh:
        for name, seq in seqs:
            filtered = "".join(seq[c] for c in keep)
            fh.write(f">{name}\n{filtered}\n")

    written += 1
    if (i + 1) % 5000 == 0:
        print(f"  [{i+1}/{len(pairs)}] written={written} skipped={skipped}")

print(f"Done. Written={written} Skipped={skipped}")
print(f"Alignment files in: {ALIGN_DIR}")
