#!/usr/bin/env python3

# compute meta-TE nucleosome occupancy profiles by TE class

# normalize each TE insertion and calculate average profiles

# compare occupancy between the 5' flank and TE body

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, kruskal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
from nucleosome_common import (  # noqa: E402
    load_wig,
    load_chrom_sizes,
    load_te_bed,
    TE_CLASS_COLORS,
    TE_CLASSES,
    style_ax,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--wig", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--te_bed", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--window", type=int, default=2000)
    p.add_argument("--min_insertions", type=int, default=10)
    p.add_argument("--max_te_len", type=int, default=5000)
    p.add_argument("--n_bins", type=int, default=200)
    p.add_argument("--normalize", choices=["zscore", "minmax"], default="zscore")
    p.add_argument("--out_dir", default="mac_te_out")
    return p.parse_args()


def extract_profile(arr, mid, window, n_bins):
    """Extract window around mid, rebin to n_bins; None if out of range."""
    s = mid - window
    e = mid + window
    if s < 0 or e >= len(arr):
        return None
    seg = arr[s:e].astype(float)
    if seg.std() < 1e-9:
        return None
    x = np.linspace(0, len(seg) - 1, n_bins)
    rebin = np.interp(x, np.arange(len(seg)), seg)
    return rebin


def normalize(seg, mode):
    if mode == "minmax":
        lo, hi = seg.min(), seg.max()
        return (seg - lo) / (hi - lo + 1e-9)
    return (seg - seg.mean()) / (seg.std() + 1e-9)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    chrom_sizes = load_chrom_sizes(args.fai)
    occ = load_wig(args.wig, chrom_sizes)
    te = load_te_bed(args.te_bed)
    te = te[te["length"] <= args.max_te_len].copy()

    # restrict classes to those present
    present = [c for c in TE_CLASSES if (te["te_class"] == c).sum() >= args.min_insertions]
    profiles = {c: [] for c in present}

    for _, row in te.iterrows():
        cls = row["te_class"]
        if cls not in profiles:
            continue
        arr = occ.get(row["chrom"])
        if arr is None:
            continue
        mid = (row["start"] + row["end"]) // 2
        seg = extract_profile(arr, mid, args.window, args.n_bins)
        if seg is None:
            continue
        profiles[cls].append(normalize(seg, args.normalize))

    if not any(profiles.values()):
        sys.exit("No valid profiles computed")

    mean_p = {}
    se_p = {}
    for c, stack in profiles.items():
        if not stack:
            continue
        mat = np.vstack(stack)
        mean_p[c] = mat.mean(axis=0)
        se_p[c] = mat.std(axis=0) / np.sqrt(mat.shape[0])

    x = np.linspace(-args.window, args.window, args.n_bins)

    # ---- output averaged profile -------------------------------------
    rows = [{"species": args.species, "te_class": c, "rel_pos_bp": float(xi),
             "mean": float(m), "se": float(s)}
            for c in mean_p
            for xi, m, s in zip(x, mean_p[c], se_p[c])]
    pd.DataFrame(rows).to_csv(out_dir / f"{args.species}_meta_te_profiles.tsv",
                              sep="\t", index=False)

    # ---- flank vs body stats -----------------------------------------
    n_flank = args.n_bins // 4          # outermost 25 % = flank
    n_body = slice(n_flank, -n_flank)   # central 50 % = body
    summ_rows = []
    body_vals_all = []
    for c in mean_p:
        mat = np.vstack(profiles[c])
        body = mat[:, n_body].mean(axis=1)
        flank5 = mat[:, :n_flank].mean(axis=1)
        body_vals_all.append(body)
        if len(body) >= 5:
            try:
                _, p = wilcoxon(body, flank5)
            except ValueError:
                p = np.nan
        else:
            p = np.nan
        summ_rows.append({
            "species": args.species,
            "te_class": c,
            "n_insertions": mat.shape[0],
            "mean_flank": float(flank5.mean()),
            "mean_body": float(body.mean()),
            "flank_minus_body": float(flank5.mean() - body.mean()),
            "wilcoxon_body_vs_flank_p": p,
        })
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(out_dir / f"{args.species}_meta_te_summary.tsv", sep="\t", index=False)

    # cross-class Kruskal-Wallis on body occupancy
    if len(body_vals_all) >= 3:
        _, kp = kruskal(*body_vals_all)
        with open(out_dir / "cross_class_test.txt", "w") as fh:
            fh.write(f"kruskal_wallis_body_occupancy_p={kp:.4g}\n")

    # ---- figures ------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for c in mean_p:
        ax.plot(x / 1000, mean_p[c], color=TE_CLASS_COLORS.get(c, "#999"),
                lw=1.5, label=c)
        ax.fill_between(x / 1000, mean_p[c] - se_p[c], mean_p[c] + se_p[c],
                        color=TE_CLASS_COLORS.get(c, "#999"), alpha=0.15)
    ax.axvline(0, color="grey", ls="--", lw=0.8)
    ax.axvspan(-args.window / 2000, 0, color="grey", alpha=0.05)
    ax.set_xlabel("Distance from TE midpoint (kb)", fontsize=9)
    ax.set_ylabel(f"Mean normalised occupancy ({args.normalize})", fontsize=9)
    ax.set_title(f"{args.species} — meta-TE nucleosome profiles", fontsize=10)
    ax.legend(fontsize=8, frameon=False, ncol=2)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.species}_meta_te_profiles.pdf")
    fig.savefig(fig_dir / f"{args.species}_meta_te_profiles.png", dpi=200)
    plt.close(fig)

    # per-class body vs flank summary bar
    if len(summ):
        fig, ax = plt.subplots(figsize=(max(5, 0.7 * len(summ) + 1), 4))
        x = np.arange(len(summ))
        w = 0.38
        ax.bar(x - w / 2, summ["mean_body"], w, label="TE body", color="#d62728", alpha=0.85)
        ax.bar(x + w / 2, summ["mean_flank"], w, label="5' flank", color="#1f77b4", alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(summ["te_class"], rotation=30, ha="right", fontsize=8)
        ax.set_ylabel(f"Mean normalised occupancy ({args.normalize})", fontsize=9)
        ax.legend(fontsize=8, frameon=False)
        style_ax(ax)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{args.species}_body_vs_flank.pdf")
        plt.close(fig)

    print(f"[{args.species}] Done. profiles for {len(mean_p)} classes.")


if __name__ == "__main__":
    main()
