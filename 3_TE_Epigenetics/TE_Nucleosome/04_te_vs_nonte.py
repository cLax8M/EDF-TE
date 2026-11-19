#!/usr/bin/env python3
"""
04_te_vs_nonte.py  —  TE vs non-TE (gene / background) occupancy comparison.

Compares the nucleosome occupancy signal between TE-overlapping regions and
non-TE regions at the 500 bp window level, for one or more species / tracks:

  A. Mean occupancy of TE+ windows (≥50 % EarlGrey overlap) vs TE− windows
     (no TE overlap).
  B. Same comparison restricted to TE+ vs *gene* windows and TE+ vs
     *background* (non-TE, non-genic) windows.
  C. Per-superfamily delta-occupancy (TE − background).

An optional ``--paired`` flag (default on when >1 species passed) averages
per-species deltas and tests the cross-species TE+/TE− difference with the
Wilcoxon signed-rank test.

Statistics: Mann–Whitney U (per species) + BH-FDR; Wilcoxon signed-rank
(cross-species).

Inputs:
  --wig          occupancy WIG (repeatable) or --wig_list file
  --fai          genome .fasta.fai (repeatable, matching --wig order)
  --te_bed       EarlGrey filteredRepeats.bed (repeatable)
  --gff          gene GFF (repeatable; optional)
  --species      species label (repeatable; default s1, s2, ...)
  --win          window size (default 500)
  --min_te_frac  TE-overlap threshold (default 0.5)
  --out_dir      output directory

Outputs:
  te_vs_nonte_summary.tsv
  per_species/<species>_windows.tsv.gz
  figures/<species>_te_vs_nonte_boxplot.pdf
  figures/te_vs_nonte_deltas.pdf
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
from nucleosome_common import (  # noqa: E402
    load_wig,
    load_chrom_sizes,
    load_te_bed,
    te_intervals_by_chrom,
    load_gene_intervals,
    window_metrics,
    bh_fdr,
    cliffs_delta,
    style_ax,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--wig", action="append", default=[])
    p.add_argument("--wig_list", default=None,
                   help="file with one WIG path per line (alternative to --wig)")
    p.add_argument("--fai", action="append", default=[])
    p.add_argument("--te_bed", action="append", default=[])
    p.add_argument("--gff", action="append", default=[])
    p.add_argument("--species", action="append", default=[])
    p.add_argument("--win", type=int, default=500)
    p.add_argument("--min_te_frac", type=float, default=0.5)
    p.add_argument("--out_dir", default="te_vs_nonte_out")
    return p.parse_args()


def window_overlap(chrom_sizes, win, ivs_dict):
    """Per-window fractional coverage of intervals (returns DataFrame)."""
    rows = []
    for chrom, size in chrom_sizes.items():
        ivs = ivs_dict.get(chrom, np.zeros((0, 2), dtype=np.int32))
        s0 = ivs[:, 0]
        e0 = ivs[:, 1]
        n = size // win
        for i in range(n):
            s = i * win
            e = min((i + 1) * win, size)
            lo = np.searchsorted(s0, s, side="right") - 1
            hi = np.searchsorted(e0, e, side="left")
            lo = max(lo, 0)
            hi = min(hi, len(s0) - 1)
            cov = 0.0
            if lo <= hi and len(s0) > 0:
                seg = np.clip(np.minimum(e0[lo:hi + 1], e) -
                              np.maximum(s0[lo:hi + 1], s), 0, None)
                cov = seg.sum()
            rows.append((chrom, s, e, cov / (e - s)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "frac"])


def process_species(wig, fai, te_bed, gff, sp, win, min_te_frac):
    chrom_sizes = load_chrom_sizes(fai)
    occ = load_wig(wig, chrom_sizes)
    te = load_te_bed(te_bed)
    te_by_chrom = te_intervals_by_chrom(te)
    gene_by_chrom = load_gene_intervals(gff) if gff else None

    te_win = window_overlap(chrom_sizes, win, te_by_chrom)
    gene_win = window_overlap(chrom_sizes, win, gene_by_chrom) if gene_by_chrom else None

    m_list = []
    for chrom, size in chrom_sizes.items():
        arr = occ.get(chrom)
        if arr is None:
            continue
        m = window_metrics(arr, win)
        n = len(m["mean_occ"])
        for i in range(n):
            m_list.append((m["mean_occ"][i], m["cv"][i],
                           m["ac_val"][i], bool(m["valid"][i])))
    mdf = pd.DataFrame(m_list, columns=["mean_occ", "cv", "ac_val", "valid"])
    mdf = mdf.reindex(range(len(te_win)))
    wdf = pd.concat([te_win, mdf[["mean_occ", "cv", "ac_val", "valid"]]], axis=1)
    wdf = wdf[wdf["valid"]].copy()

    te_flag = wdf["frac"] >= min_te_frac
    wdf["group"] = np.where(te_flag, "TE+", "TE-")
    if gene_win is not None:
        gene_flag = gene_win["frac"].reindex(wdf.index, fill_value=0) >= 0.5
        bg = (~te_flag) & (~gene_flag)
        wdf["context"] = np.where(te_flag, "TE",
                                  np.where(gene_flag, "Gene", "Background"))
    else:
        wdf["context"] = np.where(te_flag, "TE", "Background")
    wdf["species"] = sp
    return wdf, te


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    per_dir = out_dir / "per_species"
    fig_dir = out_dir / "figures"
    for d in (out_dir, per_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    wigs = list(args.wig)
    if args.wig_list:
        wigs += [l.strip() for l in open(args.wig_list) if l.strip()]

    n = len(wigs)
    if n == 0:
        sys.exit("No WIG provided (use --wig or --wig_list)")
    fai = list(args.fai)
    teb = list(args.te_bed)
    gff = list(args.gff) or [None] * n
    sps = list(args.species) or [f"s{i+1}" for i in range(n)]
    if len(fai) != n or len(teb) != n or len(sps) != n:
        sys.exit("--wig, --fai, --te_bed and --species must have the same count")

    summary_rows = []
    deltas = []
    all_wdf = []

    for i in range(n):
        print(f"[{sps[i]}] Processing ...", flush=True)
        wdf, te = process_species(wigs[i], fai[i], teb[i],
                                  gff[i], sps[i], args.win, args.min_te_frac)
        wdf.to_csv(per_dir / f"{sps[i]}_windows.tsv.gz", sep="\t",
                   index=False, compression="gzip")

        te_vals = wdf.loc[wdf["group"] == "TE+", "mean_occ"].to_numpy(float)
        nt_vals = wdf.loc[wdf["group"] == "TE-", "mean_occ"].to_numpy(float)
        bg_vals = wdf.loc[wdf["context"] == "Background", "mean_occ"].to_numpy(float)
        gene_vals = wdf.loc[wdf["context"] == "Gene", "mean_occ"].to_numpy(float)

        def _mw(x, y):
            if len(x) >= 5 and len(y) >= 5:
                _, p = mannwhitneyu(x, y, alternative="two-sided")
            else:
                p = np.nan
            return p

        row = {
            "species": sps[i],
            "n_te_plus": len(te_vals),
            "n_te_minus": len(nt_vals),
            "mean_occ_TE+": te_vals.mean() if len(te_vals) else np.nan,
            "mean_occ_TE-": nt_vals.mean() if len(nt_vals) else np.nan,
            "delta_TE_vs_TEminus": (te_vals.mean() - nt_vals.mean()) if len(te_vals) and len(nt_vals) else np.nan,
            "mw_TE_vs_TEminus_p": _mw(te_vals, nt_vals),
            "cliffs_delta_TE_vs_TEminus": cliffs_delta(te_vals, nt_vals),
        }
        if len(bg_vals):
            row["mean_occ_background"] = bg_vals.mean()
            row["delta_TE_vs_bg"] = te_vals.mean() - bg_vals.mean()
            row["mw_TE_vs_bg_p"] = _mw(te_vals, bg_vals)
        if len(gene_vals):
            row["mean_occ_gene"] = gene_vals.mean()
            row["delta_TE_vs_gene"] = te_vals.mean() - gene_vals.mean()
            row["mw_TE_vs_gene_p"] = _mw(te_vals, gene_vals)
        summary_rows.append(row)
        deltas.append(row)
        all_wdf.append(wdf)

        # per-species boxplot
        fig, ax = plt.subplots(figsize=(5, 4))
        data = [wdf.loc[wdf["context"] == c, "mean_occ"] for c in ("TE", "Gene", "Background")]
        data = [d for d in data if len(d)]
        labels = [c for c, d in zip(("TE", "Gene", "Background"), data) if len(d)]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False,
                        medianprops=dict(color="black", lw=1.2))
        colors = {"TE": "#d62728", "Gene": "#1f77b4", "Background": "#999999"}
        for patch, lab in zip(bp["boxes"], labels):
            patch.set_facecolor(colors.get(lab, "#999"))
            patch.set_alpha(0.8)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("Mean nucleosome occupancy", fontsize=9)
        ax.set_title(f"{sps[i]} — TE vs non-TE", fontsize=10)
        style_ax(ax)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{sps[i]}_te_vs_nonte_boxplot.pdf")
        plt.close(fig)

    summ = pd.DataFrame(summary_rows)
    # BH-FDR over per-species Mann-Whitney p-values
    mwcols = [c for c in summ.columns if c.endswith("_p")]
    pvals = summ[mwcols].to_numpy(float).ravel()
    mask = np.isfinite(pvals)
    if mask.any():
        adj = bh_fdr(pvals[mask])
        full = np.full_like(pvals, np.nan)
        full[mask] = adj
        for j, c in enumerate(mwcols):
            summ[f"{c}_bh"] = full[j::len(mwcols)]

    # cross-species signed-rank on delta_TE_vs_bg
    d = pd.DataFrame(deltas)["delta_TE_vs_bg"].dropna()
    if len(d) >= 3:
        try:
            _, wp = wilcoxon(d)
        except ValueError:
            wp = np.nan
        summ.attrs["cross_species_wilcoxon_p"] = wp
        with open(out_dir / "cross_species_test.txt", "w") as fh:
            fh.write(f"delta_TE_vs_bg  n={len(d)}  "
                     f"median={d.median():.4f}  wilcoxon_p={wp:.4g}\n")

    summ.to_csv(out_dir / "te_vs_nonte_summary.tsv", sep="\t", index=False)

    # delta plot across species
    if len(summ) > 1:
        fig, ax = plt.subplots(figsize=(max(5, 0.5 * len(summ) + 2), 4))
        x = np.arange(len(summ))
        cols = ["#d62728" if v > 0 else "#1f77b4"
                for v in summ["delta_TE_vs_bg"].fillna(0)]
        ax.bar(x, summ["delta_TE_vs_bg"], color=cols)
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(summ["species"], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Δ occupancy (TE − background)", fontsize=9)
        ax.set_title("TE-vs-background nucleosome occupancy", fontsize=10)
        style_ax(ax)
        fig.tight_layout()
        fig.savefig(fig_dir / "te_vs_nonte_deltas.pdf")
        plt.close(fig)

    print("Done.")


if __name__ == "__main__":
    main()
