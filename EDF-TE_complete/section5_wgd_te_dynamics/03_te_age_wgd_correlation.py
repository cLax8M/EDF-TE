#!/usr/bin/env python3
"""
03_te_age_wgd_correlation.py  —  TE burst timing relative to WGD events.

Tests whether TE expansions (measured by Kimura divergence peaks) are
temporally associated with inferred WGD events (measured by Ks peaks).

Analysis:
  1. For each species with a WGD call: extract WGD Ks peak and translate to
     absolute time using a species-specific substitution rate (estimated from
     calibrated nodes in the species tree)
  2. For each TE superfamily: identify the Kimura divergence peak (mode of
     distribution using KDE)
  3. Translate Kimura divergence to absolute time using the same substitution
     rate (assuming clock-like evolution of TE copies after insertion)
  4. Correlate WGD time with TE burst time: scatter plot per TE class
  5. Classify TEs as: pre-WGD, co-WGD (within 2× TE rate CI), post-WGD

Inputs:
  --wgd_summary   wgd_summary.tsv from 01_wgd_detection.py
  --kimura_tsv    kimura_species.tsv from Section2/04_kimura_divergence.py
  --rate_tsv      TSV with per-species substitution rate estimates
                  (species, rate_per_site_per_myr)
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory

Outputs:
  te_wgd_timing.tsv
  figures/te_age_wgd_scatter.pdf
  figures/te_class_wgd_timing_bar.pdf
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#984ea3",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}

CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--wgd_summary",  required=True)
    p.add_argument("--kimura_tsv",   required=True)
    p.add_argument("--rate_tsv",     required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="te_wgd_timing_out")
    p.add_argument("--co_wgd_fold",  type=float, default=2.0,
                   help="Fold range around WGD time considered 'co-WGD' (default: 2)")
    return p.parse_args()


def ks_to_mya(ks: float, rate_per_myr: float) -> float:
    """Convert Ks to million years ago using synonymous substitution rate."""
    if rate_per_myr <= 0:
        return np.nan
    return ks / (2 * rate_per_myr)


def kimura_to_mya(kimura: float, rate_per_myr: float) -> float:
    """Convert Kimura divergence to million years ago (insertion date)."""
    if rate_per_myr <= 0:
        return np.nan
    return kimura / (2 * rate_per_myr)


def classify_te_timing(te_mya: float, wgd_mya: float, fold: float) -> str:
    if np.isnan(te_mya) or np.isnan(wgd_mya):
        return "unknown"
    lo = wgd_mya / fold
    hi = wgd_mya * fold
    if lo <= te_mya <= hi:
        return "co-WGD"
    elif te_mya > wgd_mya * fold:
        return "pre-WGD"
    else:
        return "post-WGD"


# =============================================================================
# Figures
# =============================================================================

def plot_scatter(df: pd.DataFrame, clade_map: dict, out_dir: Path) -> None:
    """Scatter: WGD Mya vs TE peak Mya, coloured by TE class."""
    valid = df.dropna(subset=["wgd_mya", "te_peak_mya"])
    if valid.empty:
        print("No valid WGD × TE timing data for scatter.")
        return

    classes = valid["te_class"].unique()
    fig, ax = plt.subplots(figsize=(7, 6))

    for cls in classes:
        sub = valid[valid["te_class"] == cls]
        color = TE_CLASS_COLORS.get(cls, "#999")
        ax.scatter(sub["wgd_mya"], sub["te_peak_mya"],
                   color=color, s=40, alpha=0.7, label=cls, zorder=3)

    max_val = max(valid["wgd_mya"].max(), valid["te_peak_mya"].max()) * 1.1
    ax.plot([0, max_val], [0, max_val], "--", color="grey", linewidth=0.8, zorder=1)

    r, p = spearmanr(valid["wgd_mya"], valid["te_peak_mya"])
    ax.text(0.05, 0.97, f"Spearman r={r:.2f}, p={p:.2g}",
            transform=ax.transAxes, fontsize=9, va="top")
    ax.set_xlabel("WGD time (Mya)", fontsize=10)
    ax.set_ylabel("TE burst peak time (Mya)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "te_age_wgd_scatter.pdf")
    fig.savefig(out_dir / "te_age_wgd_scatter.png", dpi=200)
    plt.close(fig)
    print("Saved: te_age_wgd_scatter.pdf")


def plot_timing_bar(df: pd.DataFrame, out_dir: Path) -> None:
    """Stacked bar: fraction pre/co/post WGD per TE class."""
    agg = df.groupby(["te_class", "wgd_timing"]).size().unstack(
        "wgd_timing", fill_value=0)
    for cat in ["pre-WGD", "co-WGD", "post-WGD"]:
        if cat not in agg.columns:
            agg[cat] = 0
    agg_pct = agg[["pre-WGD", "co-WGD", "post-WGD"]].div(
        agg.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(8, 5))
    bottom = np.zeros(len(agg_pct))
    timing_colors = {"pre-WGD": "#a6cee3", "co-WGD": "#e41a1c", "post-WGD": "#b2df8a"}
    for cat, col in timing_colors.items():
        ax.bar(range(len(agg_pct)), agg_pct[cat], bottom=bottom,
               color=col, label=cat, width=0.8)
        bottom += agg_pct[cat].values
    ax.set_xticks(range(len(agg_pct)))
    ax.set_xticklabels(agg_pct.index, rotation=30, fontsize=9)
    ax.set_ylabel("% TE families", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "te_class_wgd_timing_bar.pdf")
    plt.close(fig)
    print("Saved: te_class_wgd_timing_bar.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))

    wgd     = pd.read_csv(args.wgd_summary, sep="\t")
    kimura  = pd.read_csv(args.kimura_tsv,  sep="\t")
    rates   = pd.read_csv(args.rate_tsv,    sep="\t").set_index("species")

    # Only species with WGD calls
    wgd_sp = wgd[wgd["wgd_call"] == True].copy()
    wgd_sp["rate"] = wgd_sp["species"].map(
        rates["rate_per_site_per_myr"] if "rate_per_site_per_myr" in rates.columns
        else rates.iloc[:, 0])

    wgd_sp["wgd_mya"] = wgd_sp.apply(
        lambda r: ks_to_mya(float(r["peak_ks"].split(";")[0])
                             if r["peak_ks"] else np.nan,
                             r.get("rate", np.nan)), axis=1)

    rows = []
    for _, wgd_row in wgd_sp.iterrows():
        sp       = wgd_row["species"]
        wgd_mya  = wgd_row["wgd_mya"]
        rate     = wgd_row.get("rate", np.nan)

        sp_kim = kimura[kimura["species"] == sp] if "species" in kimura.columns else kimura
        if "mean_kimura" in sp_kim.columns:
            # Use per-class mean from summary
            for cls, grp in sp_kim.groupby("te_class") if "te_class" in sp_kim.columns else []:
                te_mya = kimura_to_mya(float(grp["mean_kimura"].mean()), rate)
                timing = classify_te_timing(te_mya, wgd_mya, args.co_wgd_fold)
                rows.append({"species": sp, "clade": clade_map.get(sp, "?"),
                             "te_class": cls, "te_peak_mya": te_mya,
                             "wgd_mya": wgd_mya, "wgd_timing": timing})

    result = pd.DataFrame(rows)
    if result.empty:
        print("No timing data generated — check column names in inputs.", file=sys.stderr)
        sys.exit(0)

    result.to_csv(out_dir / "te_wgd_timing.tsv", sep="\t", index=False)
    plot_scatter(result, clade_map, fig_dir)
    plot_timing_bar(result, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
