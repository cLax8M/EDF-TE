#!/usr/bin/env python3
"""Compute normalized synteny ratios from collinear summary and null means."""

import argparse
import sys
from pathlib import Path

CLADES_ALL = ["AP", "AS", "BA", "BL", "BS", "CH", "CR", "EN", "GL", "HO", "KX", "MI", "MO", "MU", "NU", "OL", "ZO"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--collinear_dir", required=True)
    p.add_argument("--shuffle_dir", required=True)
    p.add_argument("--out_dir", required=True)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--clade", action="append", choices=CLADES_ALL)
    group.add_argument("--all", action="store_true")
    return p.parse_args()


def load_collinear_summary(collinear_dir, clade):
    path = Path(collinear_dir) / clade / f"{clade}_collinear_summary.tsv"
    rows = []
    if not path.exists():
        return rows
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 5:
                rows.append((cols[0], cols[1], int(cols[2]), int(cols[3]), int(cols[4])))
    return rows


def load_null_means(shuffle_dir, clade):
    path = Path(shuffle_dir) / clade / f"{clade}_null_means.tsv"
    data = {}
    if not path.exists():
        return data
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 4:
                data[(cols[0], cols[1])] = float(cols[2])
    return data


def main():
    args = parse_args()
    target_clades = CLADES_ALL if args.all else args.clade
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for clade in target_clades:
        collinear_rows = load_collinear_summary(args.collinear_dir, clade)
        nulls = load_null_means(args.shuffle_dir, clade)
        if not collinear_rows:
            continue
        out_sub = out_dir / clade
        out_sub.mkdir(parents=True, exist_ok=True)
        out_path = out_sub / f"{clade}_synteny_ratios.tsv"
        rows = []
        for sp_a, sp_b, n_orthologs, c_ij, n_adj in collinear_rows:
            s_ij = nulls.get((sp_a, sp_b), 0.0)
            raw = max((c_ij - s_ij) / n_adj, 0.0) if n_adj > 0 else 0.0
            rows.append((sp_a, sp_b, n_orthologs, c_ij, s_ij, n_adj, raw))
        if not rows:
            continue
        max_raw = max(r[-1] for r in rows) if rows else 0.0
        with open(out_path, "w") as fh:
            fh.write("species_A\tspecies_B\tn_orthologs\tc_ij\ts_ij\tN_ij\traw_score\tr_ij_norm\n")
            for sp_a, sp_b, n_orthologs, c_ij, s_ij, n_adj, raw in sorted(rows, key=lambda x: -x[-1]):
                norm = raw / max_raw if max_raw > 0 else 0.0
                fh.write(f"{sp_a}\t{sp_b}\t{n_orthologs}\t{c_ij}\t{s_ij:.4f}\t{n_adj}\t{raw:.6f}\t{norm:.6f}\n")
                all_rows.append((clade, sp_a, sp_b, n_orthologs, c_ij, s_ij, n_adj, raw, norm))
        print(f"Wrote {out_path}")
    global_out = out_dir / "all_clades_synteny_ratios.tsv"
    with open(global_out, "w") as fh:
        fh.write("clade\tspecies_A\tspecies_B\tn_orthologs\tc_ij\ts_ij\tN_ij\traw_score\tr_ij_norm\n")
        for row in all_rows:
            fh.write("\t".join(map(str, row)) + "\n")
    print(f"Wrote {global_out}")


if __name__ == "__main__":
    main()
