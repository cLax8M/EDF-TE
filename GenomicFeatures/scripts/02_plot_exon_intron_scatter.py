#!/usr/bin/env python3
"""
02_plot_exon_intron_scatter.py — Scatter plot of % exon vs % intron in genome.

Recreates the classic exon-intron-intergenic scatter (Sebé-Pedrós et al. style)
from the ALL_SPECIES_ALL_SUMMARY.txt. Points colored by clade, key species labeled.

Usage:
  python 02_plot_exon_intron_scatter.py <ALL_SPECIES_ALL_SUMMARY.txt> <out_dir>
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

CLADE_ORDER = ['CR', 'OL', 'ZO', 'GL', 'BL', 'EN', 'KX', 'NU',
               'BS', 'BA', 'AS', 'CH', 'MI', 'MO', 'MU', 'HO', 'AP']

CLADE_COLORS = {
    'CR': '#a6cee3', 'OL': '#1f78b4', 'ZO': '#b2df8a', 'GL': '#33a02c',
    'BL': '#fb9a99', 'EN': '#e31a1c', 'KX': '#fdbf6f', 'NU': '#ff7f00',
    'BS': '#cab2d6', 'BA': '#6a3d9a', 'AS': '#ffff99', 'CH': '#b15928',
    'MI': '#8dd3c7', 'MO': '#bebada', 'MU': '#fb8072', 'HO': '#80b1d3',
    'AP': '#fdb462'
}

CLADE_NAMES = {
    'CR': 'Cryptomycota', 'OL': 'Olpidiomycota', 'ZO': 'Zoopagomycota',
    'GL': 'Glomeromycota', 'BL': 'Blastocladiomycota', 'EN': 'Entomophthoromycota',
    'KX': 'Kickxellomycota', 'NU': 'Nucleariida', 'BS': 'Basidiobolomycota',
    'BA': 'Basidiomycota', 'AS': 'Ascomycota', 'CH': 'Chytridiomycota',
    'MI': 'Microsporidia', 'MO': 'Mortierellomycota', 'MU': 'Mucoromycota',
    'HO': 'Holozoa', 'AP': 'Aphelidiomycota'
}

# Species to label on the plot (short label → full species code)
HIGHLIGHT = {
    "HO_Homosap":    ("Hsap",   "#E15759"),
    "HO_Drome":      ("Dmel",   "#E15759"),
    "HO_Salpirose":  ("Sros",   "#F28E2B"),
    "HO_Monbrevi":   ("Mbre",   "#F28E2B"),
    "HO_Capowcza":   ("Cowc",   "#F28E2B"),
    "HO_Sphaearti":  ("Sarc",   "#59A14F"),
    "HO_Trichadha":  ("Tadh",   "#59A14F"),
    "AS_Sacce1":     ("Scer",   "#B07AA1"),
    "AS_Schpo1":     ("Spom",   "#B07AA1"),
    "AS_Ncrass2":    ("Ncra",   "#B07AA1"),
    "MI_Enccu1":     ("Ecun",   "#8dd3c7"),
    "MI_Edhae1":     ("Ehae",   "#8dd3c7"),
    "MI_Nosbo1":     ("Nbo",    "#8dd3c7"),
    "MU_Rhiarbi":    ("Rarb",   "#fb8072"),
    "MU_Mucci2":     ("Mcir",   "#fb8072"),
    "MO_Morvert2":   ("Mver",   "#bebada"),
    "CH_Batra1":     ("Bden",   "#b15928"),
    "HO_Naegrub":    ("Ngru",   "#80b1d3"),
    "HO_Acancaste":  ("Acas",   "#80b1d3"),
    "HO_Dicdisco":   ("Ddis",   "#80b1d3"),
    "AP_Amoprot":    ("Apro",   "#fdb462"),
    "HO_Hydravulg":  ("Hvul",   "#E15759"),
}


def load_data(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t")
    df["Species"] = df["Species"].str.replace(r"\.aa$", "", regex=True)
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["Clade"] = df["Species"].str[:2]
    df = df[df["Clade"].isin(CLADE_ORDER)].copy()
    df["Clade"] = pd.Categorical(df["Clade"], categories=CLADE_ORDER, ordered=True)
    return df.sort_values("Clade").reset_index(drop=True)


def plot_exon_intron(df, out_dir, label_mode="highlight"):
    """% exon vs % intron scatter with identity line and configurable labels.

    label_mode:
      - "highlight": only a curated subset is labeled
      - "all": all species codes are labeled
      - "none": no labels
    """
    fig, ax = plt.subplots(figsize=(9, 9))

    # Identity line
    ax.plot([0, 90], [0, 90], "k:", linewidth=1.2, alpha=0.4, zorder=0)
    ax.text(2, 5, "identity", rotation=45, fontsize=8, color="#888888",
            fontstyle="italic", ha="left", va="bottom")

    # Plot each clade
    for clade in CLADE_ORDER:
        sub = df[df["Clade"] == clade]
        if sub.empty:
            continue
        color = CLADE_COLORS.get(clade, "#aaaaaa")
        ax.scatter(sub["Pct_Exonic"], sub["Pct_Intronic"],
                   c=color, s=45, alpha=0.75, edgecolors="white",
                   linewidths=0.3, zorder=2, label=clade)

    # Labels / highlights according to mode
    if label_mode == "highlight":
        for sp, (label, color) in HIGHLIGHT.items():
            row = df[df["Species"] == sp]
            if row.empty:
                continue
            x = row["Pct_Exonic"].values[0]
            y = row["Pct_Intronic"].values[0]
            ax.scatter(x, y, c=color, s=70, edgecolors="black",
                       linewidths=0.8, zorder=3)
            offset_x, offset_y = 1.2, 0.8
            if x > 70:
                offset_x = -1.5
            if y > 35:
                offset_y = -1.2
            ax.annotate(label, (x, y),
                        xytext=(x + offset_x, y + offset_y),
                        fontsize=7.5, fontweight="bold", color=color,
                        zorder=4)
    elif label_mode == "all":
        for _, row in df.iterrows():
            x = row["Pct_Exonic"]
            y = row["Pct_Intronic"]
            sp = row["Species"]
            ax.annotate(sp, (x, y),
                        xytext=(x + 0.4, y + 0.2),
                        fontsize=4.5, color="#2b2b2b",
                        zorder=4)

    # Axes
    ax.set_xlabel("% exon in genome", fontsize=12, fontweight="bold")
    ax.set_ylabel("% intron in genome", fontsize=12, fontweight="bold")
    ax.set_title("B. % exon and intron sequence in genome",
                 fontsize=14, fontweight="bold", loc="left")

    # Grid-like appearance
    ax.set_xlim(-1, 95)
    ax.set_ylim(-1, 50)
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45])
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)

    # Legend by clade (outside panel)
    clades_present = [c for c in CLADE_ORDER if c in df["Clade"].values]
    handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#aaa"),
                              label=f"{c} ({CLADE_NAMES.get(c, c)})")
               for c in clades_present]
    ax.legend(handles=handles, fontsize=6.5, ncol=1,
              loc="upper left", bbox_to_anchor=(1.02, 1.0),
              borderaxespad=0.0, framealpha=0.8,
              title="Clade", title_fontsize=7)

    plt.tight_layout()

    if label_mode == "highlight":
        out_base = "exon_vs_intron_scatter"
    elif label_mode == "all":
        out_base = "exon_vs_intron_scatter_all_species_labels"
    else:
        out_base = "exon_vs_intron_scatter_no_labels"

    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(out_dir, f"{out_base}.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  OK  {out_base}")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    tsv_path, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    print(f"Loading {tsv_path}")
    df = load_data(tsv_path)
    print(f"  {len(df)} species | {df['Clade'].nunique()} clades")

    # Summary for Microsporidia
    mi = df[df["Clade"] == "MI"]
    print(f"\nMicrosporidia summary (n={len(mi)}):")
    print(f"  Genome size: {mi['Genome_Size_bp'].min()/1e6:.1f} - {mi['Genome_Size_bp'].max()/1e6:.1f} Mb")
    print(f"  Pct_Exonic:     {mi['Pct_Exonic'].min():.1f} - {mi['Pct_Exonic'].max():.1f}%")
    print(f"  Pct_Intronic:   {mi['Pct_Intronic'].min():.1f} - {mi['Pct_Intronic'].max():.1f}%")
    print(f"  Pct_Intergenic: {mi['Pct_Intergenic'].min():.1f} - {mi['Pct_Intergenic'].max():.1f}%")

    print("\nPlotting exon vs intron scatter (highlight labels)...")
    plot_exon_intron(df, out_dir, label_mode="highlight")
    print("Plotting exon vs intron scatter (all species labels)...")
    plot_exon_intron(df, out_dir, label_mode="all")
    print("Plotting exon vs intron scatter (no labels)...")
    plot_exon_intron(df, out_dir, label_mode="none")

    print(f"\nDone — plots in {out_dir}")


if __name__ == "__main__":
    main()
