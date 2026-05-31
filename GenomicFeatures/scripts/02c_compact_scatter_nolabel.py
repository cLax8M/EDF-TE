#!/usr/bin/env python3
"""02c_compact_scatter_nolabel.py — Compact no-labels exon-intron scatter.

Identical to the no_labels mode of 02_plot_exon_intron_scatter.py but at
reduced heights (×0.85, ×0.75, ×0.65). NO highlights, NO black circles,
NO size changes — pure clade-colored scatter only.
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SUMMARY = "/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY.txt"
OUTDIR  = "/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures"

CLADE_ORDER = ['CR','OL','ZO','GL','BL','EN','KX','NU',
               'BS','BA','AS','CH','MI','MO','MU','HO','AP']

CLADE_COLORS = {
    'CR': '#a6cee3', 'OL': '#1f78b4', 'ZO': '#b2df8a', 'GL': '#33a02c',
    'BL': '#fb9a99', 'EN': '#e31a1c', 'KX': '#fdbf6f', 'NU': '#ff7f00',
    'BS': '#cab2d6', 'BA': '#6a3d9a', 'AS': '#ffff99', 'CH': '#b15928',
    'MI': '#8dd3c7', 'MO': '#bebada', 'MU': '#fb8072', 'HO': '#80b1d3',
    'AP': '#fdb462',
}

CLADE_NAMES = {
    'CR': 'Cryptomycota', 'OL': 'Olpidiomycota', 'ZO': 'Zoopagomycota',
    'GL': 'Glomeromycota', 'BL': 'Blastocladiomycota', 'EN': 'Entomophthoromycota',
    'KX': 'Kickxellomycota', 'NU': 'Nucleariida', 'BS': 'Basidiobolomycota',
    'BA': 'Basidiomycota', 'AS': 'Ascomycota', 'CH': 'Chytridiomycota',
    'MI': 'Microsporidia', 'MO': 'Mortierellomycota', 'MU': 'Mucoromycota',
    'HO': 'Holozoa', 'AP': 'Aphelidiomycota',
}


def load_data(path):
    df = pd.read_csv(path, sep="\t")
    df["Species"] = df["Species"].str.replace(r"\.aa$", "", regex=True)
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["Clade"] = df["Species"].str[:2]
    df = df[df["Clade"].isin(CLADE_ORDER)].copy()
    df["Clade"] = pd.Categorical(df["Clade"], categories=CLADE_ORDER, ordered=True)
    return df.sort_values("Clade").reset_index(drop=True)


def plot(df, out_dir, height_scale, suffix):
    fig, ax = plt.subplots(figsize=(9, 9 * height_scale))

    # Identity line
    ax.plot([0, 90], [0, 90], "k:", linewidth=1.2, alpha=0.4, zorder=0)
    ax.text(2, 5, "identity", rotation=45, fontsize=8, color="#888888",
            fontstyle="italic", ha="left", va="bottom")

    # Plot each clade — identical to original no_labels version
    for clade in CLADE_ORDER:
        sub = df[df["Clade"] == clade]
        if sub.empty:
            continue
        color = CLADE_COLORS.get(clade, "#aaaaaa")
        ax.scatter(sub["Pct_Exonic"], sub["Pct_Intronic"],
                   c=color, s=45, alpha=0.75, edgecolors="white",
                   linewidths=0.3, zorder=2, label=clade)

    ax.set_xlabel("% exon in genome", fontsize=12, fontweight="bold")
    ax.set_ylabel("% intron in genome", fontsize=12, fontweight="bold")
    ax.set_title("B. % exon and intron sequence in genome",
                 fontsize=14, fontweight="bold", loc="left")
    ax.set_xlim(-1, 95)
    ax.set_ylim(-1, 50)
    ax.set_xticks(range(0, 100, 10))
    ax.set_yticks(range(0, 50, 5))
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)

    present = [c for c in CLADE_ORDER if c in df["Clade"].values]
    handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#aaa"),
                              label=f"{c} ({CLADE_NAMES.get(c, c)})")
               for c in present]
    ax.legend(handles=handles, fontsize=6.5, ncol=1,
              loc="upper left", bbox_to_anchor=(1.02, 1.0),
              borderaxespad=0.0, framealpha=0.8,
              title="Clade", title_fontsize=7)

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, f"exon_vs_intron_scatter_compact_{suffix}.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  OK  compact_{suffix}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    df = load_data(SUMMARY)
    print(f"Loaded {len(df)} species")
    for scale, suffix in [(0.85, "h85"), (0.75, "h75"), (0.65, "h65")]:
        plot(df, OUTDIR, scale, suffix)
    print("Done")


if __name__ == "__main__":
    main()
