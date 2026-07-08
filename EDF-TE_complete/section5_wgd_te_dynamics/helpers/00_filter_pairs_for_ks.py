#!/usr/bin/env python3

"""Filter intra-species paralog pairs before Ks calculation.

Uses:
  - Orthogroups/Orthogroups.GeneCount.tsv from OrthoFinder
  - results/ks/paralog_pairs_for_ks.tsv (output of 02_prepare_ks_pairs.py)

Strategy:
  - Remove orthogroups where ANY species has more than a given copy number
    (e.g. >5–10 genes). This prunes very multi-copy families that create
    huge, noisy Ks clouds.
  - Keep table format identical, but reduced in size.

Output:
  - results/ks/paralog_pairs_for_ks.filtered.tsv
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Filter intra-species paralog pairs before Ks."
    )
    p.add_argument(
        "--orthofinder-dir",
        required=True,
        help=(
            "Root directory of OrthoFinder run "
            "(contains Orthogroups/Orthogroups.GeneCount.tsv)"
        ),
    )
    p.add_argument(
        "--pairs-in",
        required=True,
        help="Input TSV with paralog pairs (paralog_pairs_for_ks.tsv)",
    )
    p.add_argument(
        "--out-dir",
        required=True,
        help="Output directory (filtered TSV will be written here)",
    )
    p.add_argument(
        "--max-copies-per-species",
        type=int,
        default=10,
        help=(
            "Maximum allowed copy number in ANY species for an orthogroup to be kept "
            "(default: 10). Orthogroups with any species > this value are discarded."
        ),
    )
    return p.parse_args()


def load_gene_counts(orthofinder_dir: Path) -> pd.DataFrame:
    path = orthofinder_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"Orthogroups.GeneCount.tsv not found at {path}")

    df = pd.read_csv(path, sep="\t")
    if "Orthogroup" not in df.columns:
        raise ValueError("Missing 'Orthogroup' column in Orthogroups.GeneCount.tsv")

    return df


def load_pairs(pairs_path: Path) -> pd.DataFrame:
    if not pairs_path.is_file():
        raise FileNotFoundError(f"Pairs file not found at {pairs_path}")
    return pd.read_csv(pairs_path, sep="\t")


def filter_orthogroups_by_copy_number(
    df_counts: pd.DataFrame, max_copies_per_species: int
) -> pd.Series:
    """Return Series of orthogroup IDs that pass the copy-number filter."""

    species_cols = [
        c
        for c in df_counts.columns
        if c not in ("Orthogroup", "Total") and not c.startswith("#")
    ]
    if not species_cols:
        raise ValueError("No species columns found in Orthogroups.GeneCount.tsv")

    numeric = df_counts[species_cols].fillna(0)
    too_many = (numeric > max_copies_per_species).any(axis=1)

    keep_mask = ~too_many
    kept_ogs = df_counts.loc[keep_mask, "Orthogroup"]

    return kept_ogs


def main() -> None:
    args = parse_args()

    orthofinder_dir = Path(args.orthofinder_dir).resolve()
    pairs_in = Path(args.pairs_in).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    df_counts = load_gene_counts(orthofinder_dir)
    df_pairs = load_pairs(pairs_in)

    required_cols = {"species", "gene_id_1", "gene_id_2", "orthogroup"}
    missing = required_cols - set(df_pairs.columns)
    if missing:
        raise ValueError(
            f"Input pairs file {pairs_in} is missing required columns: {missing}"
        )

    kept_ogs = filter_orthogroups_by_copy_number(
        df_counts, args.max_copies_per_species
    )
    kept_og_set = set(kept_ogs)

    before_n_pairs = len(df_pairs)
    before_n_ogs = df_pairs["orthogroup"].nunique()

    df_filtered = df_pairs[df_pairs["orthogroup"].isin(kept_og_set)].copy()

    after_n_pairs = len(df_filtered)
    after_n_ogs = df_filtered["orthogroup"].nunique()

    out_path = out_dir / "paralog_pairs_for_ks.filtered.tsv"
    df_filtered.to_csv(out_path, sep="\t", index=False)

    import sys

    sys.stderr.write(
        "[WGD][Ks] Filtered paralog pairs by orthogroup copy number: "
        f"{before_n_pairs} pairs in {before_n_ogs} OGs -> "
        f"{after_n_pairs} pairs in {after_n_ogs} OGs; "
        f"max_copies_per_species={args.max_copies_per_species}; "
        f"output={out_path}\n"
    )


if __name__ == "__main__":
    main()
