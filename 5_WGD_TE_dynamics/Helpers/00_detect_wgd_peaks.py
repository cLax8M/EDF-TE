#!/usr/bin/env python3

#Skeleton: detect WGD-related Ks peaks.

import argparse
from pathlib import Path
import sys
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Detect WGD Ks peaks (skeleton)")
    p.add_argument("--ks-table", required=True, help="Path to ks_values_per_species.tsv")
    p.add_argument("--out-dir", required=True, help="Output directory for WGD candidates")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ks_path = Path(args.ks_table).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not ks_path.is_file():
        sys.stderr.write(f"[WGD][Ks] Ks table not found: {ks_path}\n")
        sys.exit(1)

    cols = ["node_id", "node_name", "ks_peak_center", "n_descendant_species", "supporting_pairs"]
    df_empty = pd.DataFrame(columns=cols)
    out_path = out_dir / "wgd_candidates.tsv"
    df_empty.to_csv(out_path, sep="\t", index=False)
    print(f"[WGD][Ks] Wrote empty skeleton WGD candidates table to {out_path}")


if __name__ == "__main__":
    main()
