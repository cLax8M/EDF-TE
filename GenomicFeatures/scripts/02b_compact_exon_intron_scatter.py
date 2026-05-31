#!/usr/bin/env python3
"""
02b_compact_exon_intron_scatter.py — Generate compact versions of the exon-intron scatter.

Reads the same ALL_SPECIES_ALL_SUMMARY.txt used by 02_plot_exon_intron_scatter.py
and produces three additional versions with reduced heights (×0.85, ×0.75, ×0.65).

Usage:
  python3 02b_compact_exon_intron_scatter.py
"""
import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────
SUMMARY = "/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY.txt"
OUTDIR  = "/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures"

# ── Config (same as 02_plot_exon_intron_scatter.py) ────────────────────────
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


def plot_compact(df, out_dir, height_scale, suffix):
    """Plot exon-intron scatter at a given height scale."""
    w = 9
    h = 9 * height_scale

    fig, ax = plt.subplots(figsize=(w, h))

    # Identity line
    ax.plot([0, 90], [0, 90], "k:", linewidth=1.2, alpha=0.4, zorder=0)

    # Plot each clade
    for clade in CLADE_ORDER:
        sub = df[df["Clade"] == clade]
        if sub.empty:
            continue
        color = CLADE_COLORS.get(clade, "#aaaaaa")
        ax.scatter(sub["Pct_Exonic"], sub["Pct_Intronic"],
                   c=color, s=40, alpha=0.75, edgecolors="white",
                   linewidths=0.3, zorder=2, label=clade)

    # Highlights
    for sp, (label, color) in HIGHLIGHT.items():
        row = df[df["Species"] == sp]
        if row.empty:
            continue
        x = row["Pct_Exonic"].values[0]
        y = row["Pct_Intronic"].values[0]
        ax.scatter(x, y, c=color, s=60, edgecolors="black",
                   linewidths=0.8, zorder=3)
        offset_x, offset_y = 1.0, 0.6
        if x > 70:
            offset_x = -1.5
        if y > 35:
            offset_y = -1.0
        ax.annotate(label, (x, y),
                    xytext=(x + offset_x, y + offset_y),
                    fontsize=6.5 * height_scale, fontweight="bold", color=color,
                    zorder=4)

    ax.set_xlabel("% exon in genome", fontsize=11)
    ax.set_ylabel("% intron in genome", fontsize=11)
    ax.set_title("B. % exon and intron sequence in genome",
                 fontsize=12, fontweight="bold", loc="left")
    ax.set_xlim(-1, 95)
    ax.set_ylim(-1, 50)
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45])
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)

    # Legend
    clades_present = [c for c in CLADE_ORDER if c in df["Clade"].values]
    handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#aaa"),
               label=f"{c} ({CLADE_NAMES.get(c, c)})")
               for c in clades_present]
    ax.legend(handles=handles, fontsize=5.5, ncol=1,
              loc="upper left", bbox_to_anchor=(1.02, 1.0),
              borderaxespad=0.0, framealpha=0.8,
              title="Clade", title_fontsize=6)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"exon_vs_intron_scatter_compact_{suffix}.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] compact_{suffix} (height={h:.1f})")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    df = load_data(SUMMARY)
    print(f"Loaded {len(df)} species")

    for scale, suffix in [(0.85, "h85"), (0.75, "h75"), (0.65, "h65")]:
        plot_compact(df, OUTDIR, scale, suffix)

    print(f"\nDone — compact scatter versions in {OUTDIR}")


if __name__ == "__main__":
    main()
