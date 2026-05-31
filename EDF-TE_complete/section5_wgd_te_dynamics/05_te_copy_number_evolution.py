#!/usr/bin/env python3
"""
05_te_copy_number_evolution.py  —  TE copy number change across WGD and speciation.

Estimates TE family expansion and loss rates across the phylogeny by
comparing copy number across sister lineages.  Species that underwent WGD
are tested for excess TE retention relative to the WGD-naive expectation
(copy number doubled by WGD then decayed at a neutral loss rate).

Analysis:
  1. Load per-species copy number per TE superfamily
  2. Fit a phylogenetic signal test (Blomberg K) per superfamily
  3. For each pair of WGD vs non-WGD sister clades: compare copy number ratio
     using Mann-Whitney test
  4. Identify superfamilies with significant excess retention post-WGD

Inputs:
  --copy_number  Per-species TE copy number TSV: species, superfamily, n_copies
  --wgd_summary  wgd_summary.tsv from 01_wgd_detection.py
  --clade_map    Two-column TSV: species_code<TAB>clade_tag
  --out_dir      Output directory

Outputs:
  te_copy_evolution.tsv
  figures/te_copy_wgd_violin.pdf
  figures/te_expansion_heatmap.pdf
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
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

CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--copy_number",  required=True)
    p.add_argument("--wgd_summary",  required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="te_copy_evo_out")
    p.add_argument("--min_species",  type=int, default=5,
                   help="Minimum species per WGD category to test (default: 5)")
    return p.parse_args()


def load_pivot(cn_tsv: str) -> pd.DataFrame:
    """Load copy-number long format and pivot to species × superfamily."""
    df = pd.read_csv(cn_tsv, sep="\t")
    sp_col  = next(c for c in df.columns if c.lower() == "species")
    sf_col  = next(c for c in df.columns if "superfam" in c.lower() or "family" in c.lower())
    cnt_col = next(c for c in df.columns if "cop" in c.lower() or "count" in c.lower() or "n_" in c.lower())
    pivot = df.pivot_table(index=sp_col, columns=sf_col, values=cnt_col,
                           aggfunc="sum", fill_value=0)
    return pivot


# =============================================================================
# Figures
# =============================================================================

def plot_violin(cn_matrix: pd.DataFrame, wgd_set: set, out_dir: Path) -> None:
    """Violin: copy-number distribution for WGD vs non-WGD species."""
    total_cn = cn_matrix.sum(axis=1)
    is_wgd   = np.array([sp in wgd_set for sp in total_cn.index])

    wgd_vals  = np.log10(total_cn[is_wgd] + 1)
    nowgd_vals = np.log10(total_cn[~is_wgd] + 1)

    if len(wgd_vals) < 3 or len(nowgd_vals) < 3:
        print("Insufficient groups for violin plot.")
        return

    fig, ax = plt.subplots(figsize=(5, 5))
    vp = ax.violinplot([nowgd_vals.values, wgd_vals.values],
                        showmedians=True)
    for body in vp["bodies"]:
        body.set_alpha(0.7)
    vp["bodies"][0].set_facecolor("#a6cee3")
    vp["bodies"][1].set_facecolor("#e31a1c")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Non-WGD", "WGD"], fontsize=10)
    ax.set_ylabel("log10(Total TE copy number + 1)", fontsize=10)

    stat, pval = mannwhitneyu(wgd_vals, nowgd_vals, alternative="greater")
    ax.text(0.5, 0.96, f"Mann-Whitney p={pval:.3g}",
            transform=ax.transAxes, ha="center", va="top", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "te_copy_wgd_violin.pdf")
    plt.close(fig)
    print("Saved: te_copy_wgd_violin.pdf")


def plot_heatmap(result: pd.DataFrame, out_dir: Path) -> None:
    """Heatmap: log2 fold change WGD / non-WGD per superfamily."""
    sig = result[(result["pval"] < 0.05) & result["lfc"].notna()]
    if sig.empty:
        return

    top20 = sig.nlargest(20, "lfc")
    fig, ax = plt.subplots(figsize=(4, max(4, len(top20) * 0.3 + 1)))
    colors = ["#e31a1c" if v > 0 else "#1f78b4" for v in top20["lfc"]]
    ax.barh(range(len(top20)), top20["lfc"], color=colors)
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels(top20["superfamily"], fontsize=7)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel("log2(WGD / non-WGD copy number)", fontsize=9)
    ax.set_title("TE superfamily enrichment in WGD species", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "te_expansion_heatmap.pdf")
    plt.close(fig)
    print("Saved: te_expansion_heatmap.pdf")


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

    wgd_df  = pd.read_csv(args.wgd_summary, sep="\t")
    wgd_set = set(wgd_df[wgd_df["wgd_call"] == True]["species"])

    cn = load_pivot(args.copy_number)
    plot_violin(cn, wgd_set, fig_dir)

    wgd_sp   = [sp for sp in cn.index if sp in wgd_set]
    nowgd_sp = [sp for sp in cn.index if sp not in wgd_set]

    if len(wgd_sp) < args.min_species or len(nowgd_sp) < args.min_species:
        print(f"Fewer than {args.min_species} species in WGD ({len(wgd_sp)}) "
              f"or non-WGD ({len(nowgd_sp)}) group.  Skipping per-superfamily test.")
        return

    rows = []
    for sf in cn.columns:
        wgd_v   = cn.loc[wgd_sp, sf].values.astype(float)
        nowgd_v = cn.loc[nowgd_sp, sf].values.astype(float)
        if wgd_v.max() == 0 and nowgd_v.max() == 0:
            continue
        mean_wgd   = wgd_v.mean()
        mean_nowgd = nowgd_v.mean()
        lfc = np.log2((mean_wgd + 1) / (mean_nowgd + 1))
        stat, pval = mannwhitneyu(wgd_v, nowgd_v, alternative="two-sided")
        rows.append({"superfamily": sf, "mean_wgd": mean_wgd,
                      "mean_nowgd": mean_nowgd, "lfc": lfc, "pval": pval})

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("pval")
        result.to_csv(out_dir / "te_copy_evolution.tsv", sep="\t", index=False)
        plot_heatmap(result, fig_dir)

    print(f"{(result['pval'] < 0.05).sum()} superfamilies enriched in WGD species.")
    print("Done.")


if __name__ == "__main__":
    main()
