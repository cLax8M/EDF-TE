#!/usr/bin/env python3
"""
F3_te_nucleosome_bias.py
Compute mean nucleosome occupancy (NuclPred z-scores) over TE insertions
versus genome background per species.

Usage:
    python F3_te_nucleosome_bias.py \
        <nuclpred_wig> <te_bed> <species_id> <output_tsv>

Input:
    nuclpred_wig: WIG file (fixedStep or variableStep) of z-score occupancy
                  from predict_nucleosomes.py
    te_bed:       BED file from EarlGrey ({SP}.filteredRepeats.bed)
                  columns: chrom, start, end, name, score, strand (BED6)
    species_id:   species code for output labelling

Output:
    TSV with per-superfamily TE nucleosome occupancy vs genome background:
    species, superfamily, mean_te_occupancy, median_te_occupancy,
    background_occupancy, delta_occupancy, n_insertions
"""

import sys
import os
import re
import math
from collections import defaultdict
import statistics


def parse_wig(wig_path):
    """
    Parse WIG (fixedStep or variableStep) into dict:
        { chrom: { pos: score } }
    Positions are 1-based as in WIG spec.
    Only stores positions with non-zero scores to save memory.
    """
    scores = defaultdict(dict)
    chrom = None
    step = 1
    span = 1
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
                m = re.search(r"span=(\d+)", line)
                span = int(m.group(1)) if m else 1
                pos = start
            elif line.startswith("variableStep"):
                mode = "variable"
                m = re.search(r"chrom=(\S+)", line)
                chrom = m.group(1) if m else chrom
                m = re.search(r"span=(\d+)", line)
                span = int(m.group(1)) if m else 1
            else:
                parts = line.split()
                if mode == "fixed":
                    val = float(parts[0])
                    if val != 0:
                        scores[chrom][pos] = val
                    pos += step
                elif mode == "variable" and len(parts) >= 2:
                    p = int(parts[0])
                    val = float(parts[1])
                    if val != 0:
                        scores[chrom][p] = val

    return scores


def mean_occupancy_region(scores, chrom, start, end):
    """
    Mean z-score over [start, end) (BED 0-based → WIG 1-based).
    Returns (mean, n_positions_with_data).
    """
    vals = []
    for p in range(start + 1, end + 1):  # convert to 1-based
        if p in scores.get(chrom, {}):
            vals.append(scores[chrom][p])
    if not vals:
        return None, 0
    return statistics.mean(vals), len(vals)


def parse_superfamily_from_bed_name(name):
    """Extract superfamily from EarlGrey BED name field (NAME#Class/Superfamily)."""
    m = re.search(r"#([^/]+)/(.+)$", name)
    if m:
        return m.group(1) + "/" + m.group(2)
    m = re.search(r"#(.+)$", name)
    if m:
        return m.group(1)
    return "Unknown"


def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: F3_te_nucleosome_bias.py <wig> <bed> <species_id> <output_tsv>")

    wig_path = sys.argv[1]
    bed_path = sys.argv[2]
    species  = sys.argv[3]
    out_path = sys.argv[4]

    print(f"[{species}] Loading WIG: {wig_path}")
    scores = parse_wig(wig_path)
    total_positions = sum(len(v) for v in scores.values())
    all_scores_flat = [v for chrom_scores in scores.values() for v in chrom_scores.values()]
    background = statistics.mean(all_scores_flat) if all_scores_flat else 0.0
    print(f"[{species}] WIG positions loaded: {total_positions:,} | background mean: {background:.4f}")

    # ── aggregate TE occupancy per superfamily ────────────────
    sf_occupancy = defaultdict(list)

    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end   = int(parts[2])
            except ValueError:
                continue
            name = parts[3] if len(parts) > 3 else "Unknown"
            sf = parse_superfamily_from_bed_name(name)

            mean_occ, n_pos = mean_occupancy_region(scores, chrom, start, end)
            if mean_occ is not None:
                sf_occupancy[sf].append(mean_occ)

    # ── summarise per superfamily ─────────────────────────────
    rows = []
    for sf, occ_list in sorted(sf_occupancy.items()):
        mean_occ   = statistics.mean(occ_list)
        median_occ = statistics.median(occ_list)
        rows.append({
            "species":              species,
            "superfamily":          sf,
            "mean_te_occupancy":    round(mean_occ, 4),
            "median_te_occupancy":  round(median_occ, 4),
            "background_occupancy": round(background, 4),
            "delta_occupancy":      round(mean_occ - background, 4),
            "n_insertions":         len(occ_list),
        })

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    import csv
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[{species}] Superfamilies analysed: {len(rows)}")
    print(f"[{species}] Output: {out_path}")


if __name__ == "__main__":
    main()
