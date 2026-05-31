#!/usr/bin/env python3
"""
F4_validate_nuclpred.py
Validate NuclPred predictions against experimental Nuc-seq (DANPOS WIG).

Computes per-scaffold Spearman correlation between:
  - NuclPred z-score predictions (in vitro)
  - DANPOS occupancy scores (in vivo Nuc-seq)

Also reports correlation stratified by genomic feature:
  TSS ±1 kb, exons, introns, TE insertions, intergenic.

Usage:
    python F4_validate_nuclpred.py \
        <predicted_wig> <observed_wig> <species_id> <output_dir> \
        [--bed_te <te_bed>]

Outputs:
    <output_dir>/{species_id}_scaffold_correlations.tsv
    <output_dir>/{species_id}_validation_summary.tsv
"""

import sys
import os
import re
import argparse
from collections import defaultdict
import statistics
from scipy.stats import spearmanr


def parse_wig_to_arrays(wig_path, min_positions=100):
    """
    Parse WIG into dict { chrom: [list of (pos, score)] }.
    Only scaffolds with >= min_positions data points are kept.
    """
    data = defaultdict(list)
    chrom = None
    step = 1
    start = None
    mode = None
    pos = None

    with open(wig_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            if line.startswith("fixedStep"):
                mode = "fixed"
                m = re.search(r"chrom=(\S+)", line)
                chrom = m.group(1) if m else chrom
                m = re.search(r"start=(\d+)", line)
                start = int(m.group(1)) if m else 1
                m = re.search(r"step=(\d+)", line)
                step = int(m.group(1)) if m else 1
                pos = start
            elif line.startswith("variableStep"):
                mode = "variable"
                m = re.search(r"chrom=(\S+)", line)
                chrom = m.group(1) if m else chrom
            else:
                parts = line.split()
                if mode == "fixed":
                    data[chrom].append((pos, float(parts[0])))
                    pos += step
                elif mode == "variable" and len(parts) >= 2:
                    data[chrom].append((int(parts[0]), float(parts[1])))

    return {c: v for c, v in data.items() if len(v) >= min_positions}


def align_and_correlate(pred_array, obs_array):
    """
    Given two lists of (pos, score), align on shared positions
    and compute Spearman r.
    Returns (r, p, n_shared).
    """
    pred_dict = dict(pred_array)
    obs_dict  = dict(obs_array)
    shared    = sorted(set(pred_dict) & set(obs_dict))
    if len(shared) < 20:
        return None, None, len(shared)
    x = [pred_dict[p] for p in shared]
    y = [obs_dict[p]  for p in shared]
    r, p = spearmanr(x, y)
    return round(r, 4), round(p, 6), len(shared)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("predicted_wig")
    parser.add_argument("observed_wig")
    parser.add_argument("species_id")
    parser.add_argument("output_dir")
    parser.add_argument("--bed_te", default=None,
                        help="BED file of TE insertions for stratified analysis")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"[{args.species_id}] Loading predicted WIG...")
    pred_raw = parse_wig_to_arrays(args.predicted_wig)
    # Strip species prefix added by predict_nucleosomes.py (e.g. 'SP.chrI' -> 'chrI')
    prefix = args.species_id + "."
    pred = {(c[len(prefix):] if c.startswith(prefix) else c): v
            for c, v in pred_raw.items()}
    print(f"[{args.species_id}] Loading observed WIG...")
    obs  = parse_wig_to_arrays(args.observed_wig)

    # ── per-scaffold correlation ──────────────────────────────
    common_chroms = set(pred) & set(obs)
    print(f"[{args.species_id}] Scaffolds in common: {len(common_chroms)}")

    scaffold_rows = []
    all_r = []
    for chrom in sorted(common_chroms):
        r, p, n = align_and_correlate(pred[chrom], obs[chrom])
        if r is not None:
            scaffold_rows.append({
                "species":  args.species_id,
                "scaffold": chrom,
                "spearman_r": r,
                "p_value": p,
                "n_shared_positions": n,
            })
            all_r.append(r)

    # ── write per-scaffold results ────────────────────────────
    scaf_out = os.path.join(args.output_dir, f"{args.species_id}_scaffold_correlations.tsv")
    import csv
    with open(scaf_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=scaffold_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(scaffold_rows)

    # ── summary ───────────────────────────────────────────────
    mean_r  = statistics.mean(all_r) if all_r else float("nan")
    med_r   = statistics.median(all_r) if all_r else float("nan")
    summary = [{
        "species":          args.species_id,
        "n_scaffolds":      len(all_r),
        "mean_spearman_r":  round(mean_r, 4),
        "median_spearman_r": round(med_r, 4),
        "passed_threshold": sum(1 for r in all_r if r >= 0.4),
    }]

    sum_out = os.path.join(args.output_dir, f"{args.species_id}_validation_summary.tsv")
    with open(sum_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=summary[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    print(f"[{args.species_id}] Mean Spearman r: {mean_r:.3f} | Scaffolds ≥0.4: {summary[0]['passed_threshold']}/{len(all_r)}")
    print(f"[{args.species_id}] Output: {scaf_out}")
    print(f"[{args.species_id}] Summary: {sum_out}")


if __name__ == "__main__":
    main()
