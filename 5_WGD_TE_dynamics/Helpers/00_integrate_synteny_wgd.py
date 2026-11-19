#!/usr/bin/env python3

"""Skeleton: integrate synteny (MCScanX) with WGD candidates.

Planned behavior:
  - parse MCScanX .collinearity outputs
  - overlap with Ks-based WGD candidates
  - summarize paralogons supporting each putative WGD event

Current behavior: just writes an empty summary table.
"""

import argparse
from pathlib import Path
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Integrate synteny and WGD (skeleton)")
    p.add_argument("--out-dir", required=True, help="Output directory for synteny summaries")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cols = ["species", "block_id", "n_genes", "n_pairs_wgd_like"]
    df_empty = pd.DataFrame(columns=cols)
    out_path = out_dir / "paralogons_per_species.tsv"
    df_empty.to_csv(out_path, sep="\t", index=False)
    print(f"[WGD][Syn] Wrote empty skeleton synteny summary to {out_path}")


if __name__ == "__main__":
    main()
