#!/usr/bin/env python3
"""
hmgb3_jockey_fimo_parser.py
Parse FIMO output and append one row to the summary TSV.
Called by run_hmgb3_jockey_cross_species.sh for each species.

Args: species genome_mb n_elem fimo_tsv summary_tsv
"""
import sys

sp, genome_mb, n_elem, fimo_tsv, summary_tsv = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]

seqs_with_hit = set()
scores = []
qvals = []

try:
    with open(fimo_tsv) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("motif_id"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            seq_name = parts[2]
            try:
                score = float(parts[5])
                qval = float(parts[8]) if len(parts) > 8 else 1.0
            except (ValueError, IndexError):
                continue
            seqs_with_hit.add(seq_name)
            scores.append(score)
            qvals.append(qval)
except FileNotFoundError:
    pass

n_with = len(seqs_with_hit)
frac = n_with / n_elem if n_elem > 0 else 0.0
best = max(scores) if scores else 0.0
mean_s = sum(scores) / len(scores) if scores else 0.0
min_q = min(qvals) if qvals else 1.0

with open(summary_tsv, "a") as out:
    out.write(f"{sp}\t{genome_mb}\t{n_elem}\t{n_with}\t{frac:.4f}\t{best:.3f}\t{mean_s:.3f}\t{min_q:.2e}\n")

print(f"  {sp}: {n_with}/{n_elem} ({frac*100:.1f}%) have HMGb3 site | best_score={best:.1f}")
