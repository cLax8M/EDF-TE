#!/usr/bin/env python3
"""
F5_ndr_te_density.py
Classify genome positions into NDR / NOR / Intermediate based on nucleosome
occupancy, then compute TE density enrichment per superfamily in each zone.

NDR  = Nucleosome Depleted Region  (occupancy < 25th percentile)
NOR  = Nucleosome Occupied Region   (occupancy > 75th percentile)
INT  = Intermediate

Enrichment = (TE_bp_in_zone / total_TE_bp) / (genome_bp_in_zone / total_genome_bp)
Log2 enrichment > 0 → TE prefers that zone.

Usage:
    python F5_ndr_te_density.py \
        <nuclpred_wig> <te_bed> <species_id> <output_tsv> \
        [--ndr_pct 25] [--nor_pct 75] [--min_te_bp 200]

Outputs:
    TSV: species, superfamily, zone, te_bp, genome_bp, enrichment, log2_enrichment, n_insertions
"""

import sys
import re
import os
import math
import argparse
from collections import defaultdict
import statistics


# ── WIG parser ─────────────────────────────────────────────────────────────────
def parse_wig(wig_path):
    """Parse WIG → dict {chrom: {pos: score}}. Only non-zero stored."""
    scores = defaultdict(dict)
    chrom = None; step = 1; span = 1; pos = None; mode = None

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
                m = re.search(r"span=(\d+)",  line); span = int(m.group(1)) if m else 1
            elif line.startswith("variableStep"):
                mode = "variable"
                m = re.search(r"chrom=(\S+)", line); chrom = m.group(1) if m else chrom
                m = re.search(r"span=(\d+)",  line); span = int(m.group(1)) if m else 1
            else:
                parts = line.split()
                if mode == "fixed":
                    val = float(parts[0])
                    if val != 0.0:
                        scores[chrom][pos] = val
                    pos += step
                elif mode == "variable" and len(parts) >= 2:
                    p = int(parts[0]); val = float(parts[1])
                    if val != 0.0:
                        scores[chrom][p] = val
    return scores


def parse_superfamily(name):
    if "#" in name:
        name = name.split("#", 1)[1]
    return name if name else "Unknown"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("wig")
    parser.add_argument("bed")
    parser.add_argument("species_id")
    parser.add_argument("output_tsv")
    parser.add_argument("--ndr_pct",   type=float, default=25.0)
    parser.add_argument("--nor_pct",   type=float, default=75.0)
    parser.add_argument("--min_te_bp", type=int,   default=200,
                        help="Min total TE bp per superfamily to include")
    args = parser.parse_args()

    print(f"[{args.species_id}] Loading WIG: {args.wig}")
    scores_raw = parse_wig(args.wig)

    # Strip species prefix from chrom names (NuclPred adds 'SP.' prefix)
    prefix = args.species_id + "."
    scores = {(c[len(prefix):] if c.startswith(prefix) else c): v
              for c, v in scores_raw.items()}

    # Flatten all values to compute percentile thresholds
    all_vals = [v for cv in scores.values() for v in cv.values()]
    all_vals.sort()
    n = len(all_vals)
    ndr_thresh = all_vals[int(n * args.ndr_pct / 100)]
    nor_thresh = all_vals[int(n * args.nor_pct / 100)]
    print(f"[{args.species_id}] NDR < {ndr_thresh:.3f} | NOR > {nor_thresh:.3f} "
          f"| Total positions: {n:,}")

    # Genome zone counts (in bp, ≈ WIG positions)
    genome_zone = {"NDR": 0, "NOR": 0, "INT": 0}
    for cv in scores.values():
        for v in cv.values():
            if v < ndr_thresh:
                genome_zone["NDR"] += 1
            elif v > nor_thresh:
                genome_zone["NOR"] += 1
            else:
                genome_zone["INT"] += 1
    total_genome_bp = sum(genome_zone.values())

    # ── Process BED: for each TE, compute % of its bp in each zone ────────────
    # Accumulate per superfamily: zone → total bp
    sf_zone_bp = defaultdict(lambda: defaultdict(int))
    sf_n_ins   = defaultdict(int)

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
            sf = parse_superfamily(parts[3])
            if chrom not in scores:
                continue

            chrom_scores = scores[chrom]
            sf_n_ins[sf] += 1
            for p in range(start + 1, end + 1):  # BED 0-based → WIG 1-based
                v = chrom_scores.get(p, 0.0)
                if v < ndr_thresh:
                    sf_zone_bp[sf]["NDR"] += 1
                elif v > nor_thresh:
                    sf_zone_bp[sf]["NOR"] += 1
                else:
                    sf_zone_bp[sf]["INT"] += 1

    # ── Compute enrichment ─────────────────────────────────────────────────────
    rows = []
    for sf, zone_bp in sorted(sf_zone_bp.items()):
        total_te_bp = sum(zone_bp.values())
        if total_te_bp < args.min_te_bp:
            continue
        for zone in ("NDR", "NOR", "INT"):
            te_bp     = zone_bp.get(zone, 0)
            g_bp      = genome_zone.get(zone, 1)
            te_frac   = te_bp / total_te_bp if total_te_bp else 0
            g_frac    = g_bp  / total_genome_bp if total_genome_bp else 0
            enrich    = te_frac / g_frac if g_frac else float("nan")
            log2_e    = math.log2(enrich) if enrich > 0 else float("nan")
            rows.append({
                "species":         args.species_id,
                "superfamily":     sf,
                "zone":            zone,
                "te_bp":           te_bp,
                "genome_bp":       g_bp,
                "te_fraction":     round(te_frac, 4),
                "genome_fraction": round(g_frac,  4),
                "enrichment":      round(enrich,  4) if not math.isnan(enrich) else "nan",
                "log2_enrichment": round(log2_e,  4) if not math.isnan(log2_e) else "nan",
                "n_insertions":    sf_n_ins[sf],
            })

    os.makedirs(os.path.dirname(args.output_tsv) or ".", exist_ok=True)
    import csv
    if not rows:
        print(f"[{args.species_id}] WARNING: No TE data overlapped WIG positions.")
        open(args.output_tsv, "w").close()
        return
    with open(args.output_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_sf = len({r["superfamily"] for r in rows})
    print(f"[{args.species_id}] Superfamilies: {n_sf} | Output: {args.output_tsv}")


if __name__ == "__main__":
    main()
