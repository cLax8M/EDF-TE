#!/usr/bin/env python3

# Nucleosome repeat length (NRL) and spacing analysis.

# Estimates NRL from occupancy tracks using autocorrelation and evaluates

# differences in nucleosome spacing across TE-density groups.

# Analyses:

#   - TE-rich vs TE-poor genomic windows.

#   - NRL distribution across TE-density quartiles.

# imports
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
from nucleosome_common import (  # noqa: E402
    load_wig,
    load_chrom_sizes,
    load_te_bed,
    te_intervals_by_chrom,
    window_metrics,
    bh_fdr,
    cliffs_delta,
    style_ax,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--wig",         required=True)
    p.add_argument("--fai",         required=True)
    p.add_argument("--te_bed",      required=True)
    p.add_argument("--species",     required=True)
    p.add_argument("--win",         type=int, default=10000)
    p.add_argument("--min_te_frac", type=float, default=0.5)
    p.add_argument("--min_te_bp",   type=int, default=0)
    p.add_argument("--out_dir",     default="spacing_out")
    return p.parse_args()


def te_fraction_per_window(chrom_sizes, win, te_by_chrom):
    """Fraction of each tiled window covered by TE intervals."""
    rows = []
    for chrom, size in chrom_sizes.items():
        n = size // win
        ivs = te_by_chrom.get(chrom)
        if ivs is None:
            ivs = np.zeros((0, 2), dtype=np.int32)
        s0 = ivs[:, 0]
        e0 = ivs[:, 1]
        for i in range(n):
            s = i * win
            e = min((i + 1) * win, size)
            if e <= s:
                continue
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
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "te_frac"])


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    chrom_sizes = load_chrom_sizes(args.fai)
    print(f"[{args.species}] Loading WIG ...", flush=True)
    occ = load_wig(args.wig, chrom_sizes)
    te_by_chrom = te_intervals_by_chrom(load_te_bed(args.te_bed))

    # per-window TE fraction + spacing metrics
    wdf = te_fraction_per_window(chrom_sizes, args.win, te_by_chrom)
    m_list = []
    for chrom, size in chrom_sizes.items():
        arr = occ.get(chrom)
        if arr is None:
            continue
        m = window_metrics(arr, args.win)
        n = len(m["mean_occ"])
        for i in range(n):
            m_list.append((m["ac_lag"][i], m["ac_val"][i],
                           m["mean_occ"][i], bool(m["valid"][i])))
    mdf = pd.DataFrame(m_list, columns=["ac_lag", "ac_val", "mean_occ", "valid"])
    if len(mdf) != len(wdf):
        # chromosomes may be missing from WIG; keep alignment by index
        mdf = mdf.reindex(range(len(wdf)))
    wdf = pd.concat([wdf, mdf[["ac_lag", "ac_val", "mean_occ", "valid"]]], axis=1)
    wdf = wdf[wdf["valid"]].copy()
    wdf = wdf[wdf["ac_lag"].notna()].copy()

    # TE-rich / TE-poor classification
    te_bp = (wdf["te_frac"] * args.win).astype(int)
    wdf["group"] = np.where(te_bp >= max(args.min_te_bp, args.min_te_frac * args.win),
                            "TE_rich", "TE_poor")

    wdf.to_csv(out_dir / f"{args.species}_window_spacing.tsv", sep="\t", index=False)

    # ---- comparison --------------------------------------------------
    rich = wdf.loc[wdf["group"] == "TE_rich", "ac_lag"].to_numpy(float)
    poor = wdf.loc[wdf["group"] == "TE_poor", "ac_lag"].to_numpy(float)
    if len(rich) >= 5 and len(poor) >= 5:
        u, p = mannwhitneyu(rich, poor, alternative="two-sided")
    else:
        u, p = np.nan, np.nan
    summ = pd.DataFrame([{
        "species": args.species,
        "win_bp": args.win,
        "n_te_rich": len(rich),
        "n_te_poor": len(poor),
        "median_nrl_te_rich": np.median(rich) if len(rich) else np.nan,
        "median_nrl_te_poor": np.median(poor) if len(poor) else np.nan,
        "delta_nrl_bp": (np.median(rich) - np.median(poor)) if len(rich) and len(poor) else np.nan,
        "mannwhitney_p": p,
        "bh_fdr": p,
        "cliffs_delta": cliffs_delta(rich, poor),
    }])
    summ.to_csv(out_dir / f"{args.species}_spacing_summary.tsv", sep="\t", index=False)

    # ---- figures -----------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.hist(poor, bins=60, range=(120, 260), alpha=0.6, color="#1f77b4",
            label=f"TE-poor (n={len(poor)})")
    ax.hist(rich, bins=60, range=(120, 260), alpha=0.6, color="#d62728",
            label=f"TE-rich (n={len(rich)})")
    ax.axvline(np.median(poor) if len(poor) else np.nan, color="#1f77b4", ls="--", lw=1)
    ax.axvline(np.median(rich) if len(rich) else np.nan, color="#d62728", ls="--", lw=1)
    ax.set_xlabel("NRL (bp, autocorrelation lag 140–240)", fontsize=9)
    ax.set_ylabel("windows", fontsize=9)
    ax.set_title(f"{args.species} — NRL distribution", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.species}_nrl_distribution.pdf")
    plt.close(fig)

    # NRL vs TE-density quartiles
    q = pd.qcut(wdf["te_frac"].rank(method="first"), 4,
                labels=["Q1", "Q2", "Q3", "Q4"])
    wdf["q"] = q
    meds = wdf.groupby("q", observed=True)["ac_lag"].median()
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(meds.index.astype(str), meds.values, marker="o", color="#333333")
    ax.set_xlabel("TE-density quartile", fontsize=9)
    ax.set_ylabel("Median NRL (bp)", fontsize=9)
    ax.set_title(f"{args.species} — NRL by TE density", fontsize=10)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.species}_nrl_by_te_density.pdf")
    plt.close(fig)

    med_r = np.median(rich) if len(rich) else np.nan
    med_p = np.median(poor) if len(poor) else np.nan
    print(f"[{args.species}] Done. median NRL TE-rich={med_r:.1f} "
          f"TE-poor={med_p:.1f} p={p:.4g}")


if __name__ == "__main__":
    main()
