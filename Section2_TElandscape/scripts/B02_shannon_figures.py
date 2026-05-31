#!/usr/bin/env python3
"""
B02_shannon_figures.py
======================
Generate Shannon diversity figures from B01 output table.

Figures produced (all PDF):
  1. shannon_h_by_clade_boxplot.pdf
        — Boxplot of H' per clade, coloured, points overlaid.
  2. shannon_j_by_clade_boxplot.pdf
        — Boxplot of Pielou Evenness J per clade.
  3. shannon_h_vs_total_te.pdf
        — Scatter H' vs total TE bp (log scale), coloured by clade.
  4. shannon_dominant_sf_barplot.pdf
        — Stacked barplot of dominant superfamily frequency per clade.
  5. shannon_h_vs_n_superfamilies.pdf
        — Scatter H' vs n_superfamilies, coloured by clade.

Usage:
    python B02_shannon_figures.py [FIGURES_DIR]
"""

import sys, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

# ── paths ─────────────────────────────────────────────────────────────────────
RESULTS = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape"
FIG_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RESULTS, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

SHANNON = os.path.join(RESULTS, "B01_shannon_diversity.tsv")

# ── colour palette per clade ──────────────────────────────────────────────────
CLADE_COLOURS = {
    "AP": "#e6194b", "AS": "#3cb44b", "BA": "#4363d8", "CH": "#f58231",
    "EF": "#911eb4", "GI": "#42d4f4", "HO": "#f032e6", "KI": "#bfef45",
    "MI": "#fabebe", "MO": "#469990", "MU": "#e6beff", "NU": "#9A6324",
    "OP": "#fffac8", "OR": "#800000", "ZO": "#aaffc3", "ZY": "#808000",
    "PE": "#ffd8b1",
}

def clade_of(sp):
    return sp[:2]

def load_data():
    df = pd.read_csv(SHANNON, sep="\t")
    df["clade"] = df["species"].apply(clade_of)
    return df


# ── Figure 1: H' by clade ───────────────────────────────────────────────────
def fig_h_by_clade(df):
    clades_sorted = df.groupby("clade")["H_prime"].median().sort_values(ascending=False).index.tolist()
    fig, ax = plt.subplots(figsize=(12, 5))
    data = [df.loc[df.clade == c, "H_prime"].values for c in clades_sorted]
    bp = ax.boxplot(data, labels=clades_sorted, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", linewidth=1.5))
    for i, c in enumerate(clades_sorted):
        col = CLADE_COLOURS.get(c, "#888888")
        bp["boxes"][i].set_facecolor(col)
        bp["boxes"][i].set_alpha(0.5)
        vals = df.loc[df.clade == c, "H_prime"].values
        jitter = np.random.normal(0, 0.06, len(vals))
        ax.scatter(np.full_like(vals, i+1) + jitter, vals, c=col, s=18, alpha=0.7, edgecolors="none", zorder=3)
    ax.set_ylabel("Shannon Diversity (H')")
    ax.set_xlabel("Clade")
    ax.set_title("Shannon Diversity H' by Clade")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "shannon_h_by_clade_boxplot.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 2: J by clade ────────────────────────────────────────────────────
def fig_j_by_clade(df):
    clades_sorted = df.groupby("clade")["J_evenness"].median().sort_values(ascending=False).index.tolist()
    fig, ax = plt.subplots(figsize=(12, 5))
    data = [df.loc[df.clade == c, "J_evenness"].values for c in clades_sorted]
    bp = ax.boxplot(data, labels=clades_sorted, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", linewidth=1.5))
    for i, c in enumerate(clades_sorted):
        col = CLADE_COLOURS.get(c, "#888888")
        bp["boxes"][i].set_facecolor(col)
        bp["boxes"][i].set_alpha(0.5)
        vals = df.loc[df.clade == c, "J_evenness"].values
        jitter = np.random.normal(0, 0.06, len(vals))
        ax.scatter(np.full_like(vals, i+1) + jitter, vals, c=col, s=18, alpha=0.7, edgecolors="none", zorder=3)
    ax.set_ylabel("Pielou Evenness (J)")
    ax.set_xlabel("Clade")
    ax.set_title("Pielou Evenness J by Clade")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "shannon_j_by_clade_boxplot.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 3: H' vs total TE bp ─────────────────────────────────────────────
def fig_h_vs_total_te(df):
    fig, ax = plt.subplots(figsize=(8, 6))
    for c in sorted(df.clade.unique()):
        sub = df[df.clade == c]
        ax.scatter(sub.total_te_bp, sub.H_prime, label=c,
                   c=CLADE_COLOURS.get(c, "#888"), s=25, alpha=0.7, edgecolors="none")
    ax.set_xscale("log")
    ax.set_xlabel("Total TE content (bp, log scale)")
    ax.set_ylabel("Shannon Diversity (H')")
    ax.set_title("Shannon H' vs Total TE Content")
    from scipy.stats import spearmanr
    r, p = spearmanr(df.total_te_bp, df.H_prime)
    ax.text(0.02, 0.98, f"Spearman r={r:.3f}, p={p:.1e}", transform=ax.transAxes,
            va="top", fontsize=8, bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    ax.legend(fontsize=5, ncol=3, loc="lower right", framealpha=0.5)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "shannon_h_vs_total_te.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 4: Dominant superfamily frequency per clade ───────────────────────
def fig_dominant_sf_barplot(df):
    # Top superfamilies + "Other"
    top_sf = df["dominant_superfamily"].value_counts().head(8).index.tolist()
    df2 = df.copy()
    df2["dom_group"] = df2["dominant_superfamily"].where(df2["dominant_superfamily"].isin(top_sf), "Other")
    ct = pd.crosstab(df2["clade"], df2["dom_group"], normalize="index")

    SF_COLOURS = {
        "LTR/Gypsy": "#4363d8", "LTR/Copia": "#3cb44b", "Unknown": "#aaaaaa",
        "DNA/TcMar-Fot1": "#e6194b", "LINE/Tad1": "#f58231", "LINE/I": "#911eb4",
        "DNA/hAT-Ac": "#42d4f4", "LTR/DIRS": "#f032e6", "Other": "#cccccc",
    }

    clades_sorted = ct.index.tolist()
    fig, ax = plt.subplots(figsize=(12, 5))
    bottom = np.zeros(len(clades_sorted))
    for col in ct.columns:
        vals = ct[col].values
        ax.bar(clades_sorted, vals, bottom=bottom, label=col,
               color=SF_COLOURS.get(col, "#999999"), edgecolor="white", linewidth=0.3)
        bottom += vals

    ax.set_ylabel("Proportion")
    ax.set_xlabel("Clade")
    ax.set_title("Dominant Superfamily Composition per Clade")
    ax.legend(fontsize=6, ncol=2, loc="upper right", framealpha=0.5)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "shannon_dominant_sf_barplot.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 5: H' vs n_superfamilies ─────────────────────────────────────────
def fig_h_vs_nsf(df):
    fig, ax = plt.subplots(figsize=(8, 6))
    for c in sorted(df.clade.unique()):
        sub = df[df.clade == c]
        ax.scatter(sub.n_superfamilies, sub.H_prime, label=c,
                   c=CLADE_COLOURS.get(c, "#888"), s=25, alpha=0.7, edgecolors="none")
    ax.set_xlabel("Number of Superfamilies")
    ax.set_ylabel("Shannon Diversity (H')")
    ax.set_title("Shannon H' vs Superfamily Richness")
    from scipy.stats import spearmanr
    r, p = spearmanr(df.n_superfamilies, df.H_prime)
    ax.text(0.02, 0.98, f"Spearman r={r:.3f}, p={p:.1e}", transform=ax.transAxes,
            va="top", fontsize=8, bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    ax.legend(fontsize=5, ncol=3, loc="lower right", framealpha=0.5)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "shannon_h_vs_n_superfamilies.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("[B02] Generating Shannon figures …")
    df = load_data()
    print(f"  Loaded {len(df)} species")

    fig_h_by_clade(df)
    fig_j_by_clade(df)
    fig_h_vs_total_te(df)
    fig_dominant_sf_barplot(df)
    fig_h_vs_nsf(df)

    print("[B02] Done — 5 figures generated.")
