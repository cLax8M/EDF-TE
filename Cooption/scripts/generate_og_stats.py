#!/usr/bin/env python3
"""
Generate species_og_cooption_stats.tsv from step8 and OrthoFinder gene counts.

Outputs:
  species_og_cooption_stats.tsv  (species, n_og_with_coopt, n_og_total, pct_og_coopted)

Usage:
  python3 generate_og_stats.py [--step8 PATH] [--genecount PATH] [--out PATH]
"""

import argparse
import os
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--step8",
        default="/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures/step8_coopted_genes_with_dup_class.tsv",
        help="step8_coopted_genes_with_dup_class.tsv from 09_duplication_cooption_orthofinder.py",
    )
    p.add_argument(
        "--genecount",
        default="/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups/Orthogroups.GeneCount.tsv",
        help="Orthogroups.GeneCount.tsv from OrthoFinder",
    )
    p.add_argument(
        "--out",
        default="/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures/species_og_cooption_stats.tsv",
        help="Output TSV path",
    )
    args = p.parse_args()

    # --- load step8: per-species/per-OG coopted gene table ---
    step8 = pd.read_csv(args.step8, sep="\t")
    n_coopt = (
        step8.groupby("species")["orthogroup"]
        .nunique()
        .rename("n_og_with_coopt")
    )

    # --- load OrthoFinder gene counts; strip .aa suffix from column names ---
    gc = pd.read_csv(args.genecount, sep="\t", index_col=0)
    gc.columns = [c.replace(".aa", "") for c in gc.columns]

    # n_og_total = number of orthogroups in which the species has ≥1 gene
    n_total = (gc > 0).sum().rename("n_og_total")

    # --- merge and compute percentage ---
    df = pd.DataFrame({"n_og_with_coopt": n_coopt, "n_og_total": n_total})
    df = df.dropna(subset=["n_og_total"])
    df["n_og_with_coopt"] = df["n_og_with_coopt"].fillna(0).astype(int)
    df["n_og_total"] = df["n_og_total"].astype(int)
    df["pct_og_coopted"] = (df["n_og_with_coopt"] / df["n_og_total"] * 100).round(4)
    df.index.name = "species"
    df = df.sort_index().reset_index()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] Written {len(df)} species to: {args.out}")


if __name__ == "__main__":
    main()
