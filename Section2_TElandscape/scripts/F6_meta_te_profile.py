#!/usr/bin/env python3
"""
F6_meta_te_profile.py
Compute meta-TE nucleosome occupancy profiles: mean occupancy signal
in a window centred on each TE insertion, stratified by superfamily.

Produces averaged occupancy profiles that show nucleosome positioning
patterns relative to TE insertions.

Usage:
    python F6_meta_te_profile.py \
        <wig> <te_bed> <species_id> <output_dir> \
        [--window 2000] [--min_insertions 10] [--max_te_len 5000]

Outputs:
    <output_dir>/{species_id}_meta_te_profiles.tsv
      columns: species, superfamily, relative_pos, mean_occupancy,
               se_occupancy, n_insertions
    <output_dir>/{species_id}_meta_te_summary.tsv
      summary of peak/trough positions per superfamily
"""

import sys
import os
import re
import math
import argparse
from collections import defaultdict
import statistics


# ── WIG parser ─────────────────────────────────────────────────────────────────
def parse_wig(wig_path):
    """Parse WIG → dict {chrom: array list of (pos, score)} keeping ALL values."""
    data = defaultdict(list)
    chrom = None; step = 1; pos = None; mode = None

    with open(wig_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            if line.startswith("fixedStep"):
                mode = "fixed"
                m = re.search(r"chrom=(\S+)", line); chrom = m.group(1) if m else chrom
                m = re.search(r"start=(\d+)", line); pos = int(m.group(1)) if m else 1
                m = re.search(r"step=(\d+)",  line); step = int(m.group(1)) if m else 1
            elif line.startswith("variableStep"):
                mode = "variable"
                m = re.search(r"chrom=(\S+)", line); chrom = m.group(1) if m else chrom
            else:
                parts = line.split()
                if mode == "fixed":
                    data[chrom].append((pos, float(parts[0])))
                    pos += step
                elif mode == "variable" and len(parts) >= 2:
                    data[chrom].append((int(parts[0]), float(parts[1])))
    # Convert to dict for O(1) lookup
    return {c: dict(v) for c, v in data.items()}


def parse_superfamily(name):
    if "#" in name:
        name = name.split("#", 1)[1]
    return name if name else "Unknown"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("wig")
    parser.add_argument("bed")
    parser.add_argument("species_id")
    parser.add_argument("output_dir")
    parser.add_argument("--window",          type=int,   default=2000,
                        help="bp upstream/downstream of TE midpoint (default 2000)")
    parser.add_argument("--min_insertions",  type=int,   default=10,
                        help="Min insertions per superfamily to profile")
    parser.add_argument("--max_te_len",      type=int,   default=5000,
                        help="Skip TEs longer than this (outliers)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    W = args.window

    print(f"[{args.species_id}] Loading WIG: {args.wig}")
    scores_raw = parse_wig(args.wig)

    # Strip species prefix from chrom names
    prefix = args.species_id + "."
    scores = {(c[len(prefix):] if c.startswith(prefix) else c): v
              for c, v in scores_raw.items()}

    # Compute genome-wide mean for normalisation
    all_vals = [v for cv in scores.values() for v in cv.values()]
    global_mean = statistics.mean(all_vals) if all_vals else 1.0
    print(f"[{args.species_id}] Global mean occupancy: {global_mean:.4f}")

    # ── Accumulate profiles per superfamily ───────────────────────────────────
    # sf_profiles[sf][rel_pos] = list of occupancy values
    sf_profiles = defaultdict(lambda: defaultdict(list))

    with open(args.bed) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            te_len = end - start
            if te_len > args.max_te_len or te_len < 50:
                continue
            sf = parse_superfamily(parts[3])
            if chrom not in scores:
                continue

            midpoint  = (start + end) // 2  # BED midpoint (0-based)
            chrom_d   = scores[chrom]
            win_start = midpoint - W
            win_end   = midpoint + W

            # Sample every position in window (1-based WIG)
            for rel in range(-W, W + 1):
                wpos = midpoint + rel + 1  # convert to 1-based
                if wpos < 1:
                    continue
                val = chrom_d.get(wpos, None)
                if val is not None:
                    sf_profiles[sf][rel].append(val)

    # ── Summarise: mean ± SE per position ─────────────────────────────────────
    profile_rows   = []
    summary_rows   = []
    positions      = list(range(-W, W + 1))

    for sf in sorted(sf_profiles):
        # Check minimum coverage at centre
        n_centre = len(sf_profiles[sf].get(0, []))
        if n_centre < args.min_insertions:
            continue

        means = []
        for rel in positions:
            vals = sf_profiles[sf].get(rel, [])
            if vals:
                m  = statistics.mean(vals)
                se = statistics.stdev(vals) / math.sqrt(len(vals)) if len(vals) > 1 else 0.0
            else:
                m = float("nan"); se = float("nan")
            means.append(m)
            profile_rows.append({
                "species":       args.species_id,
                "superfamily":   sf,
                "relative_pos":  rel,
                "mean_occupancy": round(m,  4) if not math.isnan(m)  else "nan",
                "se_occupancy":  round(se, 4) if not math.isnan(se) else "nan",
                "n_insertions":  len(sf_profiles[sf].get(rel, [])),
            })

        # Summary: find trough (NDR) and first flanking peak (NOR)
        valid_means = [(i, v) for i, v in zip(positions, means) if not math.isnan(v)]
        if valid_means:
            trough_pos, trough_val = min(valid_means, key=lambda x: x[1])
            peak_pos,   peak_val   = max(valid_means, key=lambda x: x[1])
            summary_rows.append({
                "species":       args.species_id,
                "superfamily":   sf,
                "n_insertions":  n_centre,
                "global_mean":   round(global_mean, 4),
                "centre_mean":   round(sf_profiles[sf][0][0] if len(sf_profiles[sf][0]) == 1
                                       else statistics.mean(sf_profiles[sf][0]), 4),
                "trough_pos":    trough_pos,
                "trough_val":    round(trough_val, 4),
                "peak_pos":      peak_pos,
                "peak_val":      round(peak_val, 4),
                "delta_trough":  round(trough_val - global_mean, 4),
            })

    # ── Write outputs ──────────────────────────────────────────────────────────
    import csv
    prof_out = os.path.join(args.output_dir, f"{args.species_id}_meta_te_profiles.tsv")
    sum_out  = os.path.join(args.output_dir, f"{args.species_id}_meta_te_summary.tsv")

    if not profile_rows:
        print(f"[{args.species_id}] WARNING: No profiles generated.")
        return

    with open(prof_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=profile_rows[0].keys(), delimiter="\t")
        w.writeheader(); w.writerows(profile_rows)

    with open(sum_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=summary_rows[0].keys(), delimiter="\t")
        w.writeheader(); w.writerows(summary_rows)

    n_sf = len(summary_rows)
    print(f"[{args.species_id}] Profiled {n_sf} superfamilies")
    print(f"[{args.species_id}] Profiles: {prof_out}")
    print(f"[{args.species_id}] Summary:  {sum_out}")


if __name__ == "__main__":
    main()
