#!/usr/bin/env python3
"""
05_synteny_ratio.py — Compute the gene synteny ratio r_ij for all within-clade
                      species pairs and normalize within each clade.

Formula (Sebé-Pedrós et al. 2013):

         (c_ij - s_ij)   /   N_ij
  r_ij = ──────────────────────────────────
         (c_max - s_max) /   N_max

where:
  c_ij   = observed collinear ortholog pairs (from step 03)
  s_ij   = mean spurious pairs over 100 shuffles (from step 04)
  N_ij   = number of adjacent ortholog-pair slots (from step 03)
  max    = the pair with the highest (c-s)/N value within the clade
           (used as normalization reference, analogous to Sphaeroforma-Creolimax
            in the original paper)

Negative raw ratios (c_ij < s_ij, i.e. less synteny than random) are floored to 0.

Output per clade:
  {RATIO_DIR}/{CLADE}/{CLADE}_synteny_ratios.tsv
  Columns: species_A  species_B  c_ij  s_ij  N_ij  raw_score  r_ij_norm

Summary across all clades:
  {RATIO_DIR}/all_clades_synteny_ratios.tsv

Usage:
  python 05_synteny_ratio.py --clade CLADE [--clade ...] | --all
                             --collinear_dir COLLINEAR_DIR
                             --shuffle_dir SHUFFLE_DIR
                             --out_dir RATIO_DIR

Author: EDF-TE project
Last updated: 2026-02-20
"""

import argparse
import sys
from pathlib import Path


CLADES_ALL = ["AP", "AS", "BA", "BL", "BS", "CH", "CR", "EN",
              "GL", "HO", "KX", "MI", "MO", "MU", "NU", "OL", "ZO"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--collinear_dir", required=True,
                   help="Directory with per-clade collinear summary TSVs (step 03)")
    p.add_argument("--shuffle_dir", required=True,
                   help="Directory with per-clade null mean TSVs (step 04)")
    p.add_argument("--out_dir", required=True, help="Output directory")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--clade", action="append", choices=CLADES_ALL)
    group.add_argument("--all", action="store_true")
    return p.parse_args()


def load_collinear_summary(collinear_dir, clade):
    """Load step 03 summary: returns dict (sp_a,sp_b) -> {c_ij, N_ij, n_orthologs}."""
    summary_file = Path(collinear_dir) / clade / f"{clade}_collinear_summary.tsv"
    data = {}
    if not summary_file.exists():
        print(f"  WARNING: collinear summary not found: {summary_file}", file=sys.stderr)
        return data
    with open(summary_file) as fh:
        next(fh)  # header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            sp_a, sp_b = cols[0], cols[1]
            data[(sp_a, sp_b)] = {
                "n_orthologs": int(cols[2]),
                "c_ij": int(cols[3]),
                "N_ij": int(cols[4]),
            }
    return data


def load_null_means(shuffle_dir, clade):
    """Load step 04 null means: returns dict (sp_a,sp_b) -> {s_ij_mean, s_ij_sd}."""
    null_file = Path(shuffle_dir) / clade / f"{clade}_null_means.tsv"
    data = {}
    if not null_file.exists():
        print(f"  WARNING: null means not found: {null_file}", file=sys.stderr)
        return data
    with open(null_file) as fh:
        next(fh)  # header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            sp_a, sp_b = cols[0], cols[1]
            try:
                data[(sp_a, sp_b)] = {
                    "s_ij_mean": float(cols[2]),
                    "s_ij_sd":   float(cols[3]),
                }
            except ValueError:
                pass
    return data


def compute_ratios(clade, collinear_dir, shuffle_dir, out_dir):
    """Compute normalized synteny ratios for one clade."""
    print(f"  [{clade}]", end=" ")

    col_data  = load_collinear_summary(collinear_dir, clade)
    null_data = load_null_means(shuffle_dir, clade)

    if not col_data:
        print("no collinear data — skipping")
        return []

    # Merge collinear + null data
    rows = []
    for (sp_a, sp_b), col in col_data.items():
        null = null_data.get((sp_a, sp_b), {"s_ij_mean": 0.0, "s_ij_sd": 0.0})
        c   = col["c_ij"]
        s   = null["s_ij_mean"]
        N   = col["N_ij"]
        n_orthologs = col["n_orthologs"]

        raw = (c - s) / N if N > 0 else 0.0
        raw = max(raw, 0.0)  # floor at 0

        rows.append({
            "species_A": sp_a,
            "species_B": sp_b,
            "n_orthologs": n_orthologs,
            "c_ij": c,
            "s_ij": s,
            "N_ij": N,
            "raw_score": raw,
        })

    # Normalize within clade: divide by maximum raw_score
    max_raw = max((r["raw_score"] for r in rows), default=0.0)
    ref_pair = None
    for r in rows:
        if r["raw_score"] == max_raw:
            ref_pair = f"{r['species_A']} × {r['species_B']}"
            break

    for r in rows:
        r["r_ij_norm"] = r["raw_score"] / max_raw if max_raw > 0 else 0.0

    print(f"{len(rows)} pairs | max_raw={max_raw:.4f} (ref: {ref_pair})")

    # Write per-clade output
    clade_out = Path(out_dir) / clade
    clade_out.mkdir(parents=True, exist_ok=True)
    out_file = clade_out / f"{clade}_synteny_ratios.tsv"

    with open(out_file, "w") as fh:
        fh.write("species_A\tspecies_B\tn_orthologs\tc_ij\ts_ij\tN_ij\t"
                 "raw_score\tr_ij_norm\n")
        for r in sorted(rows, key=lambda x: -x["r_ij_norm"]):
            fh.write(
                f"{r['species_A']}\t{r['species_B']}\t{r['n_orthologs']}\t"
                f"{r['c_ij']}\t{r['s_ij']:.2f}\t{r['N_ij']}\t"
                f"{r['raw_score']:.6f}\t{r['r_ij_norm']:.6f}\n"
            )

    return rows


def main():
    args = parse_args()
    target_clades = CLADES_ALL if args.all else args.clade

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for clade in target_clades:
        rows = compute_ratios(clade, args.collinear_dir, args.shuffle_dir, out_dir)
        for r in rows:
            r["clade"] = clade
        all_rows.extend(rows)

    # Write global summary
    global_out = out_dir / "all_clades_synteny_ratios.tsv"
    with open(global_out, "w") as fh:
        fh.write("clade\tspecies_A\tspecies_B\tn_orthologs\tc_ij\ts_ij\tN_ij\t"
                 "raw_score\tr_ij_norm\n")
        for r in all_rows:
            fh.write(
                f"{r['clade']}\t{r['species_A']}\t{r['species_B']}\t"
                f"{r['n_orthologs']}\t{r['c_ij']}\t{r['s_ij']:.2f}\t"
                f"{r['N_ij']}\t{r['raw_score']:.6f}\t{r['r_ij_norm']:.6f}\n"
            )

    print(f"\nGlobal summary: {global_out}")
    print("Step 05 complete.")


if __name__ == "__main__":
    main()
