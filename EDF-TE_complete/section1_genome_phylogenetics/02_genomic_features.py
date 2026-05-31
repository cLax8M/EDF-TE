#!/usr/bin/env python3
"""
02_genomic_features.py  —  Genome-wide feature analysis across 238 fungal species.

Computes and visualises six genomic metrics per species and clade:
  - Total gene count
  - Average exons per gene
  - Average introns per multi-exon gene
  - Average isoforms per gene
  - Percentage of genome that is exonic
  - Percentage of genome that is intergenic

Inputs (set via command-line or config.sh-derived environment variables):
  SUMMARY_TSV   : tab-separated table with one row per species and columns
                  Species, Total_Genes, Avg_Exons_Per_Gene,
                  Avg_Introns_Multiexon, Avg_Isoforms_Per_Gene,
                  Pct_Exonic, Pct_Intergenic
  CLADE_MAP     : two-column TSV mapping species codes to two-letter clade tags
  OUT_DIR       : directory for output figures

Outputs (PDF + PNG):
  - boxplot_<metric>.pdf         per-metric boxplot per clade
  - stacked_genome_composition   stacked bar: exonic / intronic / intergenic %
  - heatmap_genomic_features     z-score heatmap, species × 6 metrics

Usage:
  python3 02_genomic_features.py --summary SUMMARY.tsv --clade_map CLADE.tsv
                                 --out_dir figures/
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats

# ── Phylogenetic order (earliest-diverging first) ────────────────────────────
CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

# ── Colour palette (consistent with all downstream figures) ──────────────────
CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

CLADE_NAMES = {
    "HO": "Holozoa", "NU": "Nucleariida", "AP": "Archaeplastida",
    "CR": "Cryptomycetes", "MI": "Microsporidia", "BL": "Blastocladiomycota",
    "OL": "Olpidiales", "CH": "Chytridiomycetes", "NE": "Neocallimastigomycota",
    "BS": "Bacillariophyta", "EN": "Entomophthorales", "ZO": "Zoopagomycota",
    "KX": "Kickxellomycotina", "GL": "Glomerales", "MO": "Mortierellomycota",
    "MU": "Mucoromycotina", "BA": "Basidiomycota", "AS": "Ascomycota",
}

METRICS = [
    ("Total_Genes",           "Total genes",                    True),
    ("Avg_Exons_Per_Gene",    "Avg exons per gene",             False),
    ("Avg_Introns_Multiexon", "Avg introns/gene (multi-exon)",  False),
    ("Avg_Isoforms_Per_Gene", "Avg isoforms per gene",          False),
    ("Pct_Exonic",            "Exonic (%)",                     False),
    ("Pct_Intergenic",        "Intergenic (%)",                 False),
]


# =============================================================================
# I/O helpers
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary",   required=True,
                   help="Genomic features summary TSV (one row per species).")
    p.add_argument("--clade_map", required=True,
                   help="Two-column TSV: species_code<TAB>clade_tag.")
    p.add_argument("--out_dir",   default="figures",
                   help="Output directory for PDF/PNG figures.")
    return p.parse_args()


def load_data(summary_path: str, clade_path: str) -> pd.DataFrame:
    df = pd.read_csv(summary_path, sep="\t")
    df["Species"] = df["Species"].str.replace(r"\.aa$", "", regex=True)

    clade_map = pd.read_csv(clade_path, sep="\t", header=None,
                            names=["Species", "Clade"])
    df = df.merge(clade_map, on="Species", how="left")

    numeric_cols = [col for col, *_ in METRICS]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[df["Clade"].isin(CLADE_ORDER)].copy()
    df["Clade"] = pd.Categorical(df["Clade"], categories=CLADE_ORDER, ordered=True)
    return df.sort_values("Clade").reset_index(drop=True)


# =============================================================================
# Figures
# =============================================================================

def boxplots(df: pd.DataFrame, out_dir: Path) -> None:
    """One boxplot per metric: clades on x-axis, points jittered."""
    for col, label, log_scale in METRICS:
        fig, ax = plt.subplots(figsize=(14, 4))
        clades_present = [c for c in CLADE_ORDER if c in df["Clade"].values]
        data_by_clade = [df.loc[df["Clade"] == c, col].dropna().values
                         for c in clades_present]
        colors = [CLADE_COLORS.get(c, "#999999") for c in clades_present]

        bp = ax.boxplot(data_by_clade, patch_artist=True,
                        medianprops=dict(color="black", linewidth=1.5),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8),
                        flierprops=dict(marker=".", markersize=3, alpha=0.4))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)

        # Jitter
        rng = np.random.default_rng(42)
        for i, (vals, color) in enumerate(zip(data_by_clade, colors), start=1):
            jx = rng.uniform(-0.25, 0.25, size=len(vals)) + i
            ax.scatter(jx, vals, color=color, s=6, alpha=0.5, zorder=3)

        ax.set_xticks(range(1, len(clades_present) + 1))
        ax.set_xticklabels(clades_present, fontsize=8)
        ax.set_ylabel(label, fontsize=10)
        if log_scale:
            ax.set_yscale("log")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()

        stem = col.lower()
        fig.savefig(out_dir / f"boxplot_{stem}.pdf")
        fig.savefig(out_dir / f"boxplot_{stem}.png", dpi=200)
        plt.close(fig)


def stacked_composition(df: pd.DataFrame, out_dir: Path) -> None:
    """Stacked bar: mean exonic / intronic / intergenic % per clade."""
    clades = [c for c in CLADE_ORDER if c in df["Clade"].values]
    pct_exon = [df.loc[df["Clade"] == c, "Pct_Exonic"].mean() for c in clades]
    pct_inter = [df.loc[df["Clade"] == c, "Pct_Intergenic"].mean() for c in clades]
    pct_intron = [100 - e - i for e, i in zip(pct_exon, pct_inter)]

    x = np.arange(len(clades))
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.bar(x, pct_exon,   label="Exonic",     color="#4daf4a", alpha=0.85)
    ax.bar(x, pct_intron, bottom=pct_exon,    label="Intronic",    color="#377eb8", alpha=0.85)
    ax.bar(x, pct_inter,  bottom=[e + n for e, n in zip(pct_exon, pct_intron)],
           label="Intergenic", color="#e41a1c", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(clades, fontsize=8)
    ax.set_ylabel("Genome composition (%)", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "stacked_genome_composition.pdf")
    fig.savefig(out_dir / "stacked_genome_composition.png", dpi=200)
    plt.close(fig)


def heatmap(df: pd.DataFrame, out_dir: Path) -> None:
    """Z-score heatmap: species × 6 metrics, ordered by clade."""
    cols = [c for c, *_ in METRICS]
    labels = [lab for _, lab, *_ in METRICS]
    mat = df[cols].apply(stats.zscore, nan_policy="omit")
    mat.index = df["Species"]

    clade_colors = df["Clade"].map(CLADE_COLORS)
    g = sns.clustermap(mat, col_cluster=False, row_cluster=False,
                       row_colors=clade_colors.values,
                       cmap="RdBu_r", vmin=-3, vmax=3,
                       xticklabels=labels, yticklabels=False,
                       figsize=(8, 14), dendrogram_ratio=(0, 0))
    g.savefig(out_dir / "heatmap_genomic_features.pdf")
    g.savefig(out_dir / "heatmap_genomic_features.png", dpi=200)
    plt.close()


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data from {args.summary}")
    df = load_data(args.summary, args.clade_map)
    print(f"  {len(df)} species across {df['Clade'].nunique()} clades")

    print("Generating boxplots ...")
    boxplots(df, out_dir)

    print("Generating stacked composition bar ...")
    stacked_composition(df, out_dir)

    print("Generating heatmap ...")
    heatmap(df, out_dir)

    print(f"Figures written to: {out_dir}")


if __name__ == "__main__":
    main()
