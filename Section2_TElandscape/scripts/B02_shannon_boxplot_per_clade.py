#!/usr/bin/env python3
"""B02_shannon_boxplot_per_clade.py — Shannon H' diversity boxplot by clade.

Uses the same clade color scheme as the cladogram tree figure.
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SHANNON = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/B01_shannon_diversity.tsv"
OUTDIR  = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures"

# Canonical clade palette (single source of truth)
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from clade_palette import CLADE_COLORS, CLADE_NAMES, CLADE_ORDER, get_clade  # noqa: E402


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    df = pd.read_csv(SHANNON, sep="\t")
    df["clade"] = df["species"].apply(get_clade)

    present = [c for c in CLADE_ORDER if c in df["clade"].values]

    fig, ax = plt.subplots(figsize=(12, 5))

    data = []
    colors = []
    labels = []
    for c in present:
        vals = df[df["clade"] == c]["H_prime"].dropna().values
        data.append(vals)
        colors.append(CLADE_COLORS.get(c, "#cccccc"))
        n = len(vals)
        labels.append(f"{c}\n(n={n})")

    bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6,
                    medianprops=dict(color="black", linewidth=1.5))

    rng = np.random.default_rng(42)
    for i, (patch, c) in enumerate(zip(bp["boxes"], present)):
        patch.set_facecolor(colors[i])
        patch.set_alpha(0.7)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.5)
        vals = data[i]
        jitter = rng.normal(0, 0.1, len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                   s=12, alpha=0.5, c="k", zorder=3, linewidths=0)

    ax.set_xticks(range(1, len(present) + 1))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Shannon H' (TE diversity)", fontsize=11)
    ax.set_title("Shannon diversity index per clade", fontsize=13, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "shannon_boxplot_per_clade.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print("[DONE] shannon_boxplot_per_clade.pdf")


if __name__ == "__main__":
    main()
