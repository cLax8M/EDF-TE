#!/usr/bin/env python3
"""plot_synteny_summary.py — Missing synteny deliverables.

Produces:
  1. Boxplot: synteny ratio by clade
  2. Scatter: synteny ratio vs TE content (colored by clade)

Data sources:
  - perspecies_synteny.tsv    (mean_r_intra per species)
  - figure_species_data.tsv   (pct_TE_total per species)
"""
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────
SYNTENY_TSV  = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny/perspecies_synteny.tsv")
SPECIES_DAT  = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figure_species_data.tsv")
OUTDIR       = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny/figures_summary")

CLADE_COLORS = {
    "MU": "#E07B39", "BL": "#9C5FA8", "AP": "#59A14F",
    "MI": "#76B7B2", "CH": "#F1CE63", "BA": "#4E79A7",
    "AS": "#E15759", "NU": "#B07AA1", "MO": "#FF9DA7",
    "GL": "#9C755F", "HO": "#80b1d3", "CR": "#a6cee3",
    "OL": "#1f78b4", "ZO": "#b2df8a", "EN": "#e31a1c",
    "KX": "#fdbf6f", "BS": "#cab2d6",
}
CLADE_NAMES = {
    "MU": "Mucoromycota", "BL": "Blastocladiomycota",
    "AP": "Aphelidiomycota", "MI": "Microsporidia",
    "CH": "Chytridiomycota", "BA": "Basidiomycota",
    "AS": "Ascomycota", "NU": "Neocallimastigomycota",
    "MO": "Mortierellomycota", "GL": "Glomeromycotina",
    "HO": "Holozoa", "CR": "Cryptomycota",
    "OL": "Olpidiomycota", "ZO": "Zoopagomycota",
    "EN": "Entomophthoromycota", "KX": "Kickxellomycota",
    "BS": "Basidiobolomycota",
}

# Clade ordering for boxplot (phylogenetic order)
CLADE_ORDER = ["HO", "CR", "MI", "AP", "OL", "BL", "CH", "NU", "BS",
               "EN", "ZO", "KX", "MO", "GL", "MU", "BA", "AS"]


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────
    df_syn = pd.read_csv(SYNTENY_TSV, sep="\t")
    df_sp = pd.read_csv(SPECIES_DAT, sep="\t")

    # Merge
    df_syn["clade"] = df_syn["species"].str[:2]
    te_map = dict(zip(df_sp["taxon"], df_sp["pct_TE_total"]))
    df_syn["pct_TE"] = df_syn["species"].map(te_map)

    # Filter to species with data (non-trivial synteny = has >1 partner)
    df_plot = df_syn[df_syn["n_partners"] > 1].copy()
    print(f"[INFO] {len(df_plot)} species with synteny data (>1 partner)")

    # ── Figure 1: Boxplot synteny ratio by clade ─────────────────────────
    present_clades = [c for c in CLADE_ORDER if c in df_plot["clade"].values]

    fig, ax = plt.subplots(figsize=(12, 5))
    box_data = [df_plot[df_plot["clade"] == c]["mean_r_intra"].dropna().values
                for c in present_clades]
    box_labels = [CLADE_NAMES.get(c, c) for c in present_clades]
    box_colors = [CLADE_COLORS.get(c, "#ddd") for c in present_clades]

    bp = ax.boxplot(box_data, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black", linewidth=1.5),
                    flierprops=dict(marker="o", markersize=3, alpha=0.5))
    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Overlay individual points
    for i, (clade, data) in enumerate(zip(present_clades, box_data)):
        jitter = np.random.normal(0, 0.08, len(data))
        ax.scatter(np.full(len(data), i + 1) + jitter, data,
                  c=CLADE_COLORS.get(clade, "#ddd"), s=12, alpha=0.6,
                  edgecolor="k", linewidth=0.3, zorder=3)

    # Add n= labels
    for i, data in enumerate(box_data):
        ax.text(i + 1, ax.get_ylim()[0] - 0.02, f"n={len(data)}",
                ha="center", va="top", fontsize=6, color="#555")

    ax.set_xticklabels(box_labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean intra-clade synteny ratio", fontsize=10)
    ax.set_title("TE synteny conservation by fungal clade", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out1 = OUTDIR / "synteny_ratio_by_clade_box.pdf"
    fig.savefig(out1, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "synteny_ratio_by_clade_box.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] {out1}")

    # ── Figure 2: Synteny vs TE content scatter ──────────────────────────
    df_scatter = df_plot.dropna(subset=["pct_TE", "mean_r_intra"])
    print(f"[INFO] {len(df_scatter)} species for scatter plot")

    fig, ax = plt.subplots(figsize=(8, 6))

    for clade in present_clades:
        sub = df_scatter[df_scatter["clade"] == clade]
        if len(sub) == 0:
            continue
        ax.scatter(sub["pct_TE"], sub["mean_r_intra"],
                  c=CLADE_COLORS.get(clade, "#ddd"),
                  label=CLADE_NAMES.get(clade, clade),
                  s=25, alpha=0.7, edgecolor="k", linewidth=0.3)

    # Regression line
    x = df_scatter["pct_TE"].values
    y = df_scatter["mean_r_intra"].values
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() > 10:
        slope, intercept, r, p, se = stats.linregress(x[mask], y[mask])
        x_line = np.linspace(x[mask].min(), x[mask].max(), 100)
        ax.plot(x_line, slope * x_line + intercept, "k--", lw=1, alpha=0.7)
        ax.text(0.02, 0.98,
                f"$r$ = {r:.3f}, $p$ = {p:.2e}\n$y$ = {slope:.4f}$x$ + {intercept:.3f}",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                         edgecolor="gray", alpha=0.8))

    # Spearman as well
    rho, p_sp = stats.spearmanr(x[mask], y[mask])
    ax.text(0.02, 0.85,
            f"Spearman $\\rho$ = {rho:.3f}, $p$ = {p_sp:.2e}",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                     edgecolor="gray", alpha=0.8))

    ax.set_xlabel("TE content (% genome)", fontsize=10)
    ax.set_ylabel("Mean intra-clade synteny ratio", fontsize=10)
    ax.set_title("TE synteny conservation vs TE content", fontsize=12)
    ax.legend(fontsize=6, loc="lower right", ncol=2, framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out2 = OUTDIR / "synteny_vs_te_content.pdf"
    fig.savefig(out2, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "synteny_vs_te_content.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED] {out2}")

    # ── Summary stats ─────────────────────────────────────────────────────
    print(f"\n{'='*50}")
    print("Synteny summary by clade:")
    for c in present_clades:
        sub = df_plot[df_plot["clade"] == c]["mean_r_intra"]
        print(f"  {CLADE_NAMES.get(c,c):25s}  n={len(sub):3d}  "
              f"median={sub.median():.3f}  mean={sub.mean():.3f}")


if __name__ == "__main__":
    main()
