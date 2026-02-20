#!/usr/bin/env python3
"""
02_plot_pangenome.py
--------------------
Generate publication-quality pangenome accumulation figures from
the TSV outputs of 01_pangenome_compute.py.

Produces:
  figures/pangenome_curves.pdf / .png   — pangenome size vs N genomes (Fig. 2a style)
  figures/core_curves.pdf    / .png     — core genome size vs N genomes (Fig. 2b style)

PDFs are fully vector (matplotlib PDF backend) with editable text — suitable
for import and editing in Adobe Illustrator or Inkscape.

Usage:
    python 02_plot_pangenome.py
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("pdf")          # vector PDF backend
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages

# ── Paths ──────────────────────────────────────────────────────────────────────
RESULTS_DIR  = Path("/data/carlos/EDF-TE/PanGenome/results")
FIGURES_DIR  = Path("/data/carlos/EDF-TE/PanGenome/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ── Group metadata (same order as compute script) ─────────────────────────────
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

# 12-color palette: 11 groups + 1 for combined
# Designed to be distinguishable and Illustrator-friendly (not perceptually uniform)
PALETTE = [
    "#1f77b4",  # BL  – steel blue
    "#ff7f0e",  # CH  – orange
    "#2ca02c",  # EN  – green
    "#d62728",  # GL  – red
    "#9467bd",  # KX  – purple
    "#8c564b",  # MI  – brown
    "#e377c2",  # MO  – pink
    "#7f7f7f",  # MU  – grey
    "#bcbd22",  # ZO  – olive
    "#17becf",  # CR  – cyan
    "#aec7e8",  # BS  – light blue
    "#000000",  # ALL – black (combined)
]

# ── Styling ────────────────────────────────────────────────────────────────────
ALPHA_BAND  = 0.15    # transparency of SD error bands
LINEWIDTH   = 1.4
FONT_FAMILY = "Arial"
FONTSIZE_LABEL  = 9
FONTSIZE_TICK   = 8
FONTSIZE_LEGEND = 7
FONTSIZE_TITLE  = 10

plt.rcParams.update({
    "font.family":            FONT_FAMILY,
    "font.size":              FONTSIZE_LABEL,
    "axes.labelsize":         FONTSIZE_LABEL,
    "xtick.labelsize":        FONTSIZE_TICK,
    "ytick.labelsize":        FONTSIZE_TICK,
    "legend.fontsize":        FONTSIZE_LEGEND,
    "axes.linewidth":         0.8,
    "xtick.major.width":      0.6,
    "ytick.major.width":      0.6,
    "xtick.major.size":       3,
    "ytick.major.size":       3,
    "pdf.use14corefonts":     False,   # embed full font, not just Type1 core
    "pdf.fonttype":           42,      # TrueType → editable text in Illustrator
    "ps.fonttype":            42,
    "savefig.dpi":            300,
    "savefig.bbox":           "tight",
})


# ── Data loading ───────────────────────────────────────────────────────────────
def load_all_data():
    """Load per-group and combined TSVs. Return dict label→DataFrame."""
    data = {}
    for code in GROUP_ORDER:
        path = RESULTS_DIR / "per_group" / f"{code}_pangenome.tsv"
        if path.exists():
            data[code] = pd.read_csv(path, sep="\t")
        else:
            print(f"  [WARNING] Missing: {path}")

    combined_path = RESULTS_DIR / "combined" / "combined_pangenome.tsv"
    if combined_path.exists():
        data["ALL"] = pd.read_csv(combined_path, sep="\t")
    else:
        print(f"  [WARNING] Missing: {combined_path}")

    return data


# ── Generic plot function ──────────────────────────────────────────────────────
def make_accumulation_figure(
    data: dict,
    y_mean_col: str,
    y_sd_col:   str,
    ylabel:     str,
    title:      str,
    combined_on_separate_axis: bool = True,
) -> plt.Figure:
    """
    Draw accumulation curves for all groups + combined.

    If combined_on_separate_axis=True, the combined curve is plotted on a
    second y-axis (right) because it spans far more genomes than any single group.
    """
    keys_groups  = [k for k in data.keys() if k != "ALL"]
    has_combined = "ALL" in data

    if combined_on_separate_axis and has_combined:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax2 = ax.twinx()
        ax2.set_ylabel("Combined pangenome (all groups)", fontsize=FONTSIZE_LABEL - 1,
                        color="#000000")
        ax2.tick_params(axis="y", labelcolor="#000000", labelsize=FONTSIZE_TICK)
    else:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax2 = None

    for i, code in enumerate(keys_groups):
        if code not in data:
            continue
        df    = data[code]
        color = PALETTE[GROUP_ORDER.index(code)]
        label = f"{code} – {GROUP_NAMES[code]} (n={len(df)})"
        x     = df["n"].values
        y     = df[y_mean_col].values
        sd    = df[y_sd_col].values

        ax.plot(x, y, color=color, linewidth=LINEWIDTH, label=label, zorder=3)
        ax.fill_between(x, y - sd, y + sd, color=color, alpha=ALPHA_BAND, linewidth=0)

    # Combined on second axis or same axis
    if has_combined:
        df_c  = data["ALL"]
        xc    = df_c["n"].values
        yc    = df_c[y_mean_col].values
        sdc   = df_c[y_sd_col].values
        _ax   = ax2 if ax2 is not None else ax
        _ax.plot(xc, yc, color="#000000", linewidth=LINEWIDTH + 0.4,
                 linestyle="--", label="Combined (all groups)", zorder=4)
        _ax.fill_between(xc, yc - sdc, yc + sdc, color="#000000",
                         alpha=ALPHA_BAND * 0.7, linewidth=0)

    ax.set_xlabel("Number of genomes", fontsize=FONTSIZE_LABEL)
    ax.set_ylabel(ylabel, fontsize=FONTSIZE_LABEL)
    ax.set_title(title, fontsize=FONTSIZE_TITLE, pad=8)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    if ax2 is not None:
        ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))

    # Legend: group lines on ax, combined on ax2
    handles1, labels1 = ax.get_legend_handles_labels()
    if ax2 is not None:
        handles2, labels2 = ax2.get_legend_handles_labels()
        all_handles = handles1 + handles2
        all_labels  = labels1 + labels2
    else:
        all_handles, all_labels = handles1, labels1

    ax.legend(all_handles, all_labels,
              loc="upper left", frameon=True, framealpha=0.85,
              edgecolor="#cccccc", ncol=2)

    ax.spines["top"].set_visible(False)
    if ax2 is None:
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    return fig


# ── Save helpers ───────────────────────────────────────────────────────────────
def save_fig(fig: plt.Figure, stem: str):
    """Save figure as both PDF (vector) and PNG (raster preview)."""
    pdf_path = FIGURES_DIR / f"{stem}.pdf"
    png_path = FIGURES_DIR / f"{stem}.png"
    fig.savefig(pdf_path, format="pdf")
    fig.savefig(png_path, format="png", dpi=300)
    print(f"  Saved → {pdf_path}")
    print(f"  Saved → {png_path}")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print("Loading results TSVs ...")
    data = load_all_data()
    print(f"  Loaded data for: {list(data.keys())}")

    # ── Figure 1: Pangenome accumulation ──────────────────────────────────────
    print("\nGenerating pangenome accumulation figure ...")
    fig_pan = make_accumulation_figure(
        data,
        y_mean_col = "pangenome_mean",
        y_sd_col   = "pangenome_sd",
        ylabel     = "Pangenome size (gene families)",
        title      = "Pangenome accumulation curves",
        combined_on_separate_axis = True,
    )
    save_fig(fig_pan, "pangenome_curves")
    plt.close(fig_pan)

    # ── Figure 2: Core genome decay ───────────────────────────────────────────
    print("\nGenerating core genome decay figure ...")
    fig_core = make_accumulation_figure(
        data,
        y_mean_col = "core_mean",
        y_sd_col   = "core_sd",
        ylabel     = "Core genome size (gene families)",
        title      = "Core genome decay curves",
        combined_on_separate_axis = True,
    )
    save_fig(fig_core, "core_curves")
    plt.close(fig_core)

    # ── Figure 3: Multi-panel (both on one page, publication-style) ───────────
    print("\nGenerating combined multi-panel figure ...")
    matplotlib.use("pdf")  # ensure PDF backend
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax_idx, (y_mean_col, y_sd_col, ylabel, subtitle) in enumerate([
        ("pangenome_mean", "pangenome_sd", "Pangenome size (gene families)", "a"),
        ("core_mean",      "core_sd",     "Core genome size (gene families)", "b"),
    ]):
        ax = axes[ax_idx]
        for i, code in enumerate(GROUP_ORDER):
            if code not in data:
                continue
            df    = data[code]
            color = PALETTE[i]
            label = f"{code} – {GROUP_NAMES[code]}"
            x, y, sd = df["n"].values, df[y_mean_col].values, df[y_sd_col].values
            ax.plot(x, y, color=color, linewidth=LINEWIDTH, label=label, zorder=3)
            ax.fill_between(x, y - sd, y + sd, color=color, alpha=ALPHA_BAND, linewidth=0)

        ax.set_xlabel("Number of genomes", fontsize=FONTSIZE_LABEL)
        ax.set_ylabel(ylabel, fontsize=FONTSIZE_LABEL)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.10, 1.05, subtitle, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top", ha="left")

    # Single legend on the right panel
    handles, labels = axes[0].get_legend_handles_labels()
    axes[1].legend(handles, labels, loc="upper right", frameon=True,
                   framealpha=0.85, edgecolor="#cccccc",
                   fontsize=FONTSIZE_LEGEND, ncol=1)

    fig.tight_layout()
    save_fig(fig, "pangenome_multipanel")
    plt.close(fig)

    print("\nAll figures saved.")


if __name__ == "__main__":
    main()
