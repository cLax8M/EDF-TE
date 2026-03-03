#!/usr/bin/env python3
"""
08_perspecies_synteny.py — Per-species mean intra-clade synteny summary.

Reads the combined pairwise synteny ratios (all_clades_synteny_ratios.tsv from
step 05) and computes for each species:
  - mean_r_intra   : mean r_ij_norm with all clade partners (excluding self and
                     fragmented-assembly species)
  - median_r_intra : median r_ij_norm with clade partners
  - max_r_intra    : maximum r_ij_norm observed
  - n_partners     : number of valid clade partners used
  - fragmented     : True if species itself was auto-detected as fragmented

Assembly-fragmentation filter (same logic as 06_plot_synteny.py):
  Species with max(r_ij_norm)==0 AND mean(n_orthologs)>--min_orthologs across
  ALL their pairs are flagged as fragmented and excluded from partner counts so
  they do not bias the mean_r of their neighbours.

Usage:
  python 08_perspecies_synteny.py --ratios all_clades_synteny_ratios.tsv
                                  --out    perspecies_synteny.tsv
                                  [--min_orthologs 2000]
                                  [--exclude_species SP1 SP2 ...]

Output TSV columns:
  species | clade | mean_r_intra | median_r_intra | max_r_intra | n_partners | fragmented

Author: EDF-TE project
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ratios", required=True,
                   help="all_clades_synteny_ratios.tsv (step 05 combined output)")
    p.add_argument("--out", required=True,
                   help="Output TSV path")
    p.add_argument("--min_orthologs", type=int, default=2000,
                   help="Ortholog threshold for fragmentation detection (default 2000)")
    p.add_argument("--exclude_species", nargs="*", default=[],
                   help="Additional species to exclude from partner counts")
    return p.parse_args()


def detect_fragmented(df, min_orthologs=2000):
    """
    Return set of species where r_ij_norm == 0 in ALL pairs AND
    mean n_orthologs > min_orthologs.
    """
    fragmented = set()
    all_sp = set(df["species_A"]) | set(df["species_B"])
    for sp in all_sp:
        mask = (df["species_A"] == sp) | (df["species_B"] == sp)
        sub = df[mask]
        if sub.empty:
            continue
        if sub["r_ij_norm"].max() == 0 and sub["n_orthologs"].mean() > min_orthologs:
            fragmented.add(sp)
    return fragmented


def main():
    args = parse_args()

    df = pd.read_csv(args.ratios, sep="\t")

    # Detect fragmented assemblies
    fragmented = detect_fragmented(df, min_orthologs=args.min_orthologs)
    exclude_set = fragmented | set(args.exclude_species)

    if fragmented:
        print(f"Fragmented-assembly species detected ({len(fragmented)}): "
              + ", ".join(sorted(fragmented)), file=sys.stderr)

    # Build long-form table: one row per (species, partner) with clade and r
    rows_a = df[["clade", "species_A", "species_B", "n_orthologs", "r_ij_norm"]].copy()
    rows_a.columns = ["clade", "species", "partner", "n_orthologs", "r"]
    rows_b = df[["clade", "species_B", "species_A", "n_orthologs", "r_ij_norm"]].copy()
    rows_b.columns = ["clade", "species", "partner", "n_orthologs", "r"]
    long = pd.concat([rows_a, rows_b], ignore_index=True)

    # Remove rows where PARTNER is fragmented/excluded (don't let them bias mean)
    long = long[~long["partner"].isin(exclude_set)]

    # Compute per-species stats
    records = []
    all_species = sorted(set(df["species_A"]) | set(df["species_B"]))

    for sp in all_species:
        sub = long[long["species"] == sp]
        if sub.empty:
            clade = "?"
            mean_r = np.nan
            median_r = np.nan
            max_r = np.nan
            n = 0
        else:
            clade = sub["clade"].iloc[0]
            mean_r   = sub["r"].mean()
            median_r = sub["r"].median()
            max_r    = sub["r"].max()
            n        = len(sub)

        records.append({
            "species":        sp,
            "clade":          clade,
            "mean_r_intra":   round(mean_r,   6) if not np.isnan(mean_r)   else np.nan,
            "median_r_intra": round(median_r, 6) if not np.isnan(median_r) else np.nan,
            "max_r_intra":    round(max_r,    6) if not np.isnan(max_r)    else np.nan,
            "n_partners":     n,
            "fragmented":     sp in fragmented,
        })

    out_df = pd.DataFrame(records).sort_values(["clade", "mean_r_intra"],
                                                ascending=[True, False])

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"Written: {args.out}  ({len(out_df)} species)")

    # Print quick per-clade summary to stderr
    print("\n--- Per-clade mean_r summary ---", file=sys.stderr)
    summary = (out_df[~out_df["fragmented"]]
               .groupby("clade")["mean_r_intra"]
               .agg(["mean", "min", "max", "count"])
               .round(3))
    print(summary.to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
