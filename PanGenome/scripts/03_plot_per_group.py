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

    for ax, y_col, sd_col, ylabel, sublabel in [
        (ax1, "pangenome_mean", "pangenome_sd",
         "Pangenome size\n(gene families)", "a"),
        (ax2, "core_mean",      "core_sd",
         "Core genome size\n(gene families)", "b"),
    ]:
        y  = df[y_col].values
        sd = df[sd_col].values

        ax.plot(x, y, color=color, linewidth=LINEWIDTH, zorder=3)
        ax.fill_between(x, np.maximum(0, y - sd), y + sd,
                        color=color, alpha=ALPHA_BAND, linewidth=0)

        # Annotate final value
        ax.annotate(f"{int(round(y[-1])):,}",
                    xy=(x[-1], y[-1]),
                    xytext=(5, 0), textcoords="offset points",
                    fontsize=FONTSIZE_ANNOT, color=color,
                    va="center", ha="left")

        ax.set_xlabel("Number of genomes", fontsize=FONTSIZE_LABEL)
        ax.set_ylabel(ylabel, fontsize=FONTSIZE_LABEL)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(fmt_k))
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=6))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.12, 1.06, sublabel, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top", ha="left")

    # Annotate first-genome mean ± SD on pangenome panel (= genome size proxy)
    y1  = df["pangenome_mean"].values[0]
    sd1 = df["pangenome_sd"].values[0]
    ax1.annotate(f"n=1: {int(round(y1)):,} ± {int(round(sd1)):,}",
                 xy=(1, y1),
                 xytext=(5, -12), textcoords="offset points",
                 fontsize=FONTSIZE_ANNOT - 0.5, color="#555555",
                 va="top", ha="left")

    # Stack bar inset (pie-like summary) in top-right of ax2
    last = df.iloc[-1]
    core_f  = last["core_mean"]
    shell_f = last["shell_mean"]
    uniq_f  = last["unique_mean"]
    total   = core_f + shell_f + uniq_f
    if total > 0:
        pct_c = core_f  / total * 100
        pct_s = shell_f / total * 100
        pct_u = uniq_f  / total * 100
        txt = (f"Core:   {int(round(core_f)):>6,}  ({pct_c:.0f}%)\n"
               f"Shell:  {int(round(shell_f)):>6,}  ({pct_s:.0f}%)\n"
               f"Unique: {int(round(uniq_f)):>6,}  ({pct_u:.0f}%)")
        ax2.text(0.97, 0.97, txt, transform=ax2.transAxes,
                 fontsize=FONTSIZE_ANNOT, va="top", ha="right",
                 family="monospace",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                           edgecolor="#cccccc", linewidth=0.6))

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
