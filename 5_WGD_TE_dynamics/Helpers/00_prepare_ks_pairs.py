#!/usr/bin/env python3

#Prepare intra-species paralog pairs for Ks calculation.

import argparse
import itertools
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare intra-species paralog pairs for Ks using Orthogroups.tsv",
    )
    parser.add_argument(
        "--orthofinder-dir",
        required=True,
        help="Root directory of OrthoFinder run (contains Orthogroups/Orthogroups.tsv)",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for Ks pair table",
    )
    return parser.parse_args()


def load_orthogroups(orthofinder_dir: Path) -> pd.DataFrame:
    og_path = orthofinder_dir / "Orthogroups" / "Orthogroups.tsv"
    if not og_path.is_file():
        raise FileNotFoundError(f"Orthogroups.tsv not found at {og_path}")

    df = pd.read_csv(og_path, sep="\t")
    if "Orthogroup" not in df.columns:
        raise ValueError("Orthogroups.tsv is missing the 'Orthogroup' column")
    return df


def build_paralog_pairs(df_og: pd.DataFrame) -> pd.DataFrame:
    """Generate all intra-species paralog pairs from Orthogroups.tsv.

    For each orthogroup and species column, we parse the comma-separated
    list of gene IDs and generate all unordered pairs.
    """

    species_cols = [c for c in df_og.columns if c != "Orthogroup"]
    if not species_cols:
        raise ValueError("No species columns found in Orthogroups.tsv")

    records = []

    for _, row in df_og.iterrows():
        og_id = row["Orthogroup"]
        for sp in species_cols:
            val = row.get(sp)
            if isinstance(val, float) and pd.isna(val):
                continue
            if not isinstance(val, str) or not val.strip():
                continue

            # Genes are stored as "gene1, gene2, gene3" in OrthoFinder output
            genes = [g.strip() for g in val.split(",") if g.strip()]
            if len(genes) < 2:
                continue

            for g1, g2 in itertools.combinations(sorted(genes), 2):
                records.append(
                    {
                        "species": sp,
                        "gene_id_1": g1,
                        "gene_id_2": g2,
                        "orthogroup": og_id,
                        "dup_node_id": pd.NA,
                    }
                )

    if not records:
        return pd.DataFrame(
            columns=["species", "gene_id_1", "gene_id_2", "orthogroup", "dup_node_id"],
        )

    df_pairs = pd.DataFrame.from_records(records)
    # Ordenar para hacer la tabla un poco más legible
    df_pairs = df_pairs.sort_values(["species", "orthogroup", "gene_id_1", "gene_id_2"])
    return df_pairs


def main() -> None:
    args = parse_args()
    orthofinder_dir = Path(args.orthofinder_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        df_og = load_orthogroups(orthofinder_dir)
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"[WGD][Ks] Error loading Orthogroups.tsv: {exc}\n")
        sys.exit(1)

    df_pairs = build_paralog_pairs(df_og)
    out_path = out_dir / "paralog_pairs_for_ks.tsv"
    df_pairs.to_csv(out_path, sep="\t", index=False)

    n_pairs = len(df_pairs)
    n_species = df_pairs["species"].nunique() if n_pairs else 0
    sys.stderr.write(
        f"[WGD][Ks] Wrote {n_pairs} paralog pairs across {n_species} species to {out_path}\n",
    )


if __name__ == "__main__":
    main()
