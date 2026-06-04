#!/usr/bin/env python3
"""
chimte_merge.py
Lightweight merge of chimTE per-rep TSVs (replaces memory-hungry replicability() in chimTE).

Usage:
    python3 chimte_merge.py <tmp_dir> <out_dir> [--min_reps 2] [--min_coverage 2]

Reads tmp_dir/TE-initiated-repXX.tsv, TE-terminated-repXX.tsv, etc.
Finds TE-gene pairs present in >= min_reps replicates with mean reads >= min_coverage.
Writes final TSVs to out_dir/.
"""

import sys
import os
import argparse
from collections import defaultdict
from pathlib import Path


def load_tsv(filepath):
    """Load chimTE per-rep TSV, return list of row dicts."""
    rows = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            # Standard chimTE format: gene_id, gene_strand, gene_pos, TE_id, TE_strand, TE_pos, chim_reads
            # TE-exonized has 8 cols (extra: exonized_type)
            if len(parts) >= 7:
                rows.append({
                    "gene_id": parts[0],
                    "gene_strand": parts[1],
                    "gene_pos": parts[2],
                    "TE_id": parts[3],
                    "TE_strand": parts[4],
                    "TE_pos": parts[5],
                    "chim_reads": float(parts[6]),
                    "extra": parts[7] if len(parts) > 7 else "",
                })
    return rows


def merge_chimera_type(tmp_dir, chimera_type, min_reps, min_coverage):
    """
    Load all per-rep TSVs for a chimera type, filter by replication & coverage.
    Returns list of merged row dicts.
    """
    # Find all rep TSVs
    pattern = f"{chimera_type}-*.tsv"
    files = sorted(Path(tmp_dir).glob(pattern))
    # Exclude files with 'final' in name
    files = [f for f in files if "final" not in f.name]

    if not files:
        return []

    print(f"  {chimera_type}: found {len(files)} rep TSV(s): {[f.name for f in files]}")

    # Load all rows, track by (gene_id, TE_id, TE_pos) key
    # key → list of (chim_reads, gene_strand, gene_pos, TE_strand, extra) per replicate
    key_data = defaultdict(list)

    for f in files:
        if f.stat().st_size == 0:
            continue
        rows = load_tsv(str(f))
        for row in rows:
            key = (row["gene_id"], row["TE_id"], row["TE_pos"])
            key_data[key].append(row)

    # Filter by min_reps and min_coverage
    results = []
    for key, rep_rows in key_data.items():
        n_reps = len(rep_rows)
        if n_reps < int(min_reps):
            continue
        # Deduplicate by rep (same key in same rep file → keep unique)
        mean_reads = sum(r["chim_reads"] for r in rep_rows) / n_reps
        if mean_reads < float(min_coverage):
            continue
        # Use first row for metadata
        first = rep_rows[0]
        results.append({
            "gene_id": first["gene_id"],
            "gene_strand": first["gene_strand"],
            "gene_pos": first["gene_pos"],
            "TE_id": first["TE_id"],
            "TE_strand": first["TE_strand"],
            "TE_pos": first["TE_pos"],
            "chim_reads": round(mean_reads, 2),
            "extra": first["extra"],
            "n_reps": n_reps,
        })

    return results


def write_final_tsv(results, outpath, include_extra=False):
    """Write merged results to TSV."""
    with open(outpath, "w") as fh:
        if include_extra:
            fh.write("gene_id\tgene_strand\tgene_pos\tTE_id\tTE_strand\tTE_pos\tchim_reads\texonized_type\n")
        else:
            fh.write("gene_id\tgene_strand\tgene_pos\tTE_id\tTE_strand\tTE_pos\tchim_reads\n")
        for r in results:
            row = f"{r['gene_id']}\t{r['gene_strand']}\t{r['gene_pos']}\t{r['TE_id']}\t{r['TE_strand']}\t{r['TE_pos']}\t{r['chim_reads']}"
            if include_extra:
                row += f"\t{r['extra']}"
            fh.write(row + "\n")


def main():
    parser = argparse.ArgumentParser(description="Merge chimTE per-rep TSVs")
    parser.add_argument("tmp_dir", help="Directory containing per-rep TSVs")
    parser.add_argument("out_dir", help="Directory to write final merged TSVs")
    parser.add_argument("--min_reps", type=int, default=2, help="Min replicates (default 2)")
    parser.add_argument("--min_coverage", type=float, default=2.0, help="Min mean reads (default 2)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print(f"\n=== chimte_merge ===")
    print(f"  tmp_dir:     {args.tmp_dir}")
    print(f"  out_dir:     {args.out_dir}")
    print(f"  min_reps:    {args.min_reps}")
    print(f"  min_coverage: {args.min_coverage}")

    chimera_types = [
        ("TE-initiated",            False),
        ("TE-terminated",           False),
        ("TE-exonized_embedded",    True),
        ("TE-exonized_overlapped",  True),
        ("TE-exonized_intronic",    True),
    ]

    total = 0
    for chimera_type, has_extra in chimera_types:
        results = merge_chimera_type(args.tmp_dir, chimera_type, args.min_reps, args.min_coverage)
        outpath = os.path.join(args.out_dir, f"{chimera_type}_final.tsv")
        if results:
            write_final_tsv(results, outpath, include_extra=has_extra)
            print(f"  → {len(results)} {chimera_type} transcripts → {outpath}")
            total += len(results)
        else:
            print(f"  → 0 {chimera_type} transcripts (no output written)")

    print(f"\n=== Total: {total} chimeric TE-gene transcripts ===")


if __name__ == "__main__":
    main()
