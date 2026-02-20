#!/usr/bin/env python3
"""
03_plot_per_group.py
--------------------
Generate one figure per taxonomic group, each with two subpanels:
  - Left : pangenome accumulation (total gene families vs. N genomes)
  - Right: core genome decay      (core gene families vs. N genomes)

Both subpanels show mean ± 1 SD bands from 100 random subsamples.

Outputs:
  figures/per_group/{CODE}_pangenome.pdf  -- individual vector PDF per group
  figures/per_group/{CODE}_pangenome.png  -- raster preview
  figures/per_group_all.pdf               -- all groups in a single multipage PDF

Usage:
    python 03_plot_per_group.py
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages

# ── Paths ──────────────────────────────────────────────────────────────────────
RESULTS_DIR  = Path("/data/carlos/EDF-TE/PanGenome/results")
FIGURES_DIR  = Path("/data/carlos/EDF-TE/PanGenome/figures/per_group")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ── Group metadata ─────────────────────────────────────────────────────────────
GROUP_ORDER = ["BL", "CH", "EN", "GL", "KX", "MI", "MO", "MU", "ZO", "CR", "BS"]
GROUP_NAMES = {
    "BL": "Blastocladiomycota",
    "CH": "Chytridiomycota",
    "EN": "Entomophthoromycota",
    "GL": "Glomeromycota",
    "KX": "Kickxellomycota",
    "MI": "Microsporidia",
    "MO": "Mortierellomycota",
    "MU": "Mucoromycota",
    "ZO": "Zoopagomycota",
    "CR": "Cryptomycota",
    "BS": "Basidiobolomycota",
}

# One color per group (same palette as script 02)
PALETTE = {
    "BL": "#1f77b4",
    "CH": "#ff7f0e",
    "EN": "#2ca02c",
    "GL": "#d62728",
    "KX": "#9467bd",
    "MI": "#8c564b",
    "MO": "#e377c2",
    "MU": "#7f7f7f",
    "ZO": "#bcbd22",
    "CR": "#17becf",
    "BS": "#aec7e8",
}

ALPHA_BAND  = 0.18
LINEWIDTH   = 1.6
FONTSIZE_LABEL  = 9
FONTSIZE_TICK   = 8
FONTSIZE_TITLE  = 10
FONTSIZE_ANNOT  = 7.5

plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "font.size":          FONTSIZE_LABEL,
    "axes.labelsize":     FONTSIZE_LABEL,
    "xtick.labelsize":    FONTSIZE_TICK,
    "ytick.labelsize":    FONTSIZE_TICK,
    "axes.linewidth":     0.8,
    "xtick.major.width":  0.6,
    "ytick.major.width":  0.6,
    "xtick.major.size":   3,
    "ytick.major.size":   3,
    "pdf.fonttype":       42,   # TrueType → editable text in Illustrator
    "ps.fonttype":        42,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
})


def fmt_k(v, _):
    """Format axis tick as integer with thousands separator."""
    return f"{int(v):,}"


def make_group_figure(code: str, df: pd.DataFrame) -> plt.Figure:
    """Create a 2-panel figure for one group."""
    name  = GROUP_NAMES[code]
    color = PALETTE[code]
    n_sp  = len(df)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.8))
    fig.suptitle(f"{code} – {name}  (n = {n_sp} genomes)",
                 fontsize=FONTSIZE_TITLE, fontweight="bold", y=1.01)

    x = df["n"].values

    # ── Colours for the four curves ───────────────────────────────────────────
    c_pan    = color          # group colour  → pangenome
    c_core   = "#222222"      # near-black    → core
    c_shell  = "#e8a020"      # amber         → shell
    c_unique = "#cc3333"      # red           → unique

    # ── Panel a: pangenome + shell + unique ───────────────────────────────────
    curves_a = [
        ("pangenome_mean", "pangenome_sd", c_pan,    "Pangenome", "-"),
        ("shell_mean",     "shell_sd",     c_shell,  "Shell",     "--"),
        ("unique_mean",    "unique_sd",    c_unique, "Unique",    ":"),
    ]
    for y_col, sd_col, c, label, ls in curves_a:
        y  = df[y_col].values
        sd = df[sd_col].values
        ax1.plot(x, y, color=c, linewidth=LINEWIDTH, linestyle=ls,
                 label=label, zorder=3)
        ax1.fill_between(x, np.maximum(0, y - sd), y + sd,
                         color=c, alpha=ALPHA_BAND, linewidth=0)
        # Annotate final value
        ax1.annotate(f"{int(round(y[-1])):,}",
                     xy=(x[-1], y[-1]), xytext=(5, 0),
                     textcoords="offset points",
                     fontsize=FONTSIZE_ANNOT, color=c,
                     va="center", ha="left")

    ax1.legend(fontsize=FONTSIZE_ANNOT, frameon=False, loc="upper left")
    ax1.set_xlabel("Number of genomes", fontsize=FONTSIZE_LABEL)
    ax1.set_ylabel("Gene families", fontsize=FONTSIZE_LABEL)
    ax1.set_title("Pangenome, Shell & Unique", fontsize=FONTSIZE_LABEL,
                  pad=4, color="#333333")
    ax1.yaxis.set_major_formatter(mticker.FuncFormatter(fmt_k))
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=6))
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.text(-0.12, 1.06, "a", transform=ax1.transAxes,
             fontsize=11, fontweight="bold", va="top", ha="left")

    # Annotate per-genome gene count at n=1
    y1  = df["pangenome_mean"].values[0]
    sd1 = df["pangenome_sd"].values[0]
    ax1.annotate(f"n=1: {int(round(y1)):,} ± {int(round(sd1)):,}",
                 xy=(1, y1), xytext=(5, -12), textcoords="offset points",
                 fontsize=FONTSIZE_ANNOT - 0.5, color="#555555",
                 va="top", ha="left")

    # ── Panel b: core ─────────────────────────────────────────────────────────
    y_c  = df["core_mean"].values
    sd_c = df["core_sd"].values
    ax2.plot(x, y_c, color=c_core, linewidth=LINEWIDTH, label="Core", zorder=3)
    ax2.fill_between(x, np.maximum(0, y_c - sd_c), y_c + sd_c,
                     color=c_core, alpha=ALPHA_BAND, linewidth=0)
    ax2.annotate(f"{int(round(y_c[-1])):,}",
                 xy=(x[-1], y_c[-1]), xytext=(5, 0),
                 textcoords="offset points",
                 fontsize=FONTSIZE_ANNOT, color=c_core,
                 va="center", ha="left")

    ax2.set_xlabel("Number of genomes", fontsize=FONTSIZE_LABEL)
    ax2.set_ylabel("Gene families", fontsize=FONTSIZE_LABEL)
    ax2.set_title("Core genome decay", fontsize=FONTSIZE_LABEL,
                  pad=4, color="#333333")
    ax2.yaxis.set_major_formatter(mticker.FuncFormatter(fmt_k))
    ax2.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=6))
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.text(-0.12, 1.06, "b", transform=ax2.transAxes,
             fontsize=11, fontweight="bold", va="top", ha="left")

    fig.tight_layout()
    return fig


def main():
    print("Loading per-group TSVs ...")

    # Load summary for final counts
    summary_path = RESULTS_DIR / "summary_final_counts.tsv"
    summary = pd.read_csv(summary_path, sep="\t").set_index("group_code")

    # Multipage PDF
    multipage_path = Path("/data/carlos/EDF-TE/PanGenome/figures/per_group_all.pdf")
    with PdfPages(multipage_path) as pdf:
        for code in GROUP_ORDER:
            tsv_path = RESULTS_DIR / "per_group" / f"{code}_pangenome.tsv"
            if not tsv_path.exists():
                print(f"  [WARNING] Missing: {tsv_path}")
                continue

            df   = pd.read_csv(tsv_path, sep="\t")
            name = GROUP_NAMES[code]
            n_sp = len(df)
            print(f"\n  {code} ({name}): {n_sp} genomes, "
                  f"pangenome={summary.loc[code,'pangenome']:,}, "
                  f"core={summary.loc[code,'core']:,}, "
                  f"shell={summary.loc[code,'shell']:,}, "
                  f"unique={summary.loc[code,'unique']:,}")

            fig = make_group_figure(code, df)

            # Individual files
            pdf_path = FIGURES_DIR / f"{code}_pangenome.pdf"
            png_path = FIGURES_DIR / f"{code}_pangenome.png"
            fig.savefig(pdf_path, format="pdf")
            fig.savefig(png_path, format="png", dpi=300)
            print(f"    Saved → {pdf_path}")
            print(f"    Saved → {png_path}")

            # Add to multipage PDF
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"\n  Multipage PDF → {multipage_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
