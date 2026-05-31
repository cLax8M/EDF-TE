#!/usr/bin/env python3
"""
report_findings_v2.py
====================
Comprehensive findings report — all figures generated from data with matplotlib.
Sections:
  1. Kimura distance & Recency Index (conceptual diagram + data)
  2. Shannon TE diversity
  3. Microsynteny & TE influence (diagram + scatter + heatmap)
  4. Nucleosome insertion bias (computed from WIG + BED)
"""

import os, re, sys, io, random, gzip, math
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from scipy.stats import spearmanr
from pypdf import PdfWriter, PdfReader

# ─── Paths ─────────────────────────────────────────────────────────────────────
ROOT      = Path("/legserv/Temp/Carlos/EDF-TE")
SEC2      = ROOT / "results/Section2_TElandscape"
SYN_DIR   = ROOT / "results/Synteny"
TE_ANN    = Path("/data/carlos/EDF-TE/TE_analysis/annotated")
WIG_DIR   = SEC2 / "nucleosomes/nuclpred"
OUTDIR    = SEC2
OUT_PDF   = OUTDIR / "TE_findings_report_v2.pdf"

# ─── Clade colours ─────────────────────────────────────────────────────────────
CLADE_COLORS = {
    "AP": "#8B5E3C", "AS": "#E07B39", "BA": "#B5541E",
    "BL": "#6B9E3E", "BS": "#4B7A2E", "CH": "#3B7FBF",
    "CR": "#7B5EA7", "EN": "#C45B9E", "GL": "#3BA07F",
    "KX": "#E0A020", "MI": "#2B6B9F", "MO": "#9B6B3F",
    "MU": "#D04040", "NU": "#7B7B7B", "IC": "#5B9BA0",
}
CLADE_NAMES = {
    "AP": "Aphelida", "AS": "Ascomycota", "BA": "Basidiomycota",
    "BL": "Blastocladiomycota", "BS": "Basidiobolomycota",
    "CH": "Chytridiomycota", "CR": "Cryptomycota",
    "EN": "Entomophthoromycota", "GL": "Glomeromycota",
    "KX": "Kickxellomycota", "MI": "Microsporida",
    "MO": "Mortierellomycota", "MU": "Mucoromycota",
    "NU": "Nucleariida", "IC": "Ichthyosporia",
}

# ─── Utility: save matplotlib fig to bytes ─────────────────────────────────────
def fig_to_bytes(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.read()

def append_bytes(writer, pdf_bytes):
    reader = PdfReader(io.BytesIO(pdf_bytes))
    for page in reader.pages:
        writer.add_page(page)

def append_file(writer, path):
    if Path(path).exists():
        reader = PdfReader(str(path))
        for page in reader.pages:
            writer.add_page(page)
        print(f"  [+] {Path(path).name}")
    else:
        print(f"  [skip] {Path(path).name}")

# ─── Load data ──────────────────────────────────────────────────────────────────
print("Loading data...")
kimura = pd.read_csv(SEC2/"C01_kimura_species.tsv", sep="\t")
kimura["clade"] = kimura["species"].str.split("_").str[0]

shannon = pd.read_csv(SEC2/"B01_shannon_diversity.tsv", sep="\t")
shannon["clade"] = shannon["species"].str.split("_").str[0]

outliers = pd.read_csv(SEC2/"clade_outliers/clade_outliers.tsv", sep="\t")
synteny  = pd.read_csv(SYN_DIR/"perspecies_synteny.tsv", sep="\t")

# Merge all species data
merged = outliers.merge(kimura[["species","recency_index","burst_score","mean_div_weighted","n_superfamilies"]], on="species", how="left")
merged = merged.merge(shannon[["species","H_prime","J_evenness"]], on="species", how="left")
merged = merged.merge(synteny[["species","mean_r_intra"]], on="species", how="left")
print(f"  Merged: {len(merged)} species")

# Major clades for display
MAJOR_CLADES = ["MU","CH","GL","EN","AS","BA","BL","KX","MI","AP","BS","CR"]
clade_order  = [c for c in MAJOR_CLADES if c in merged["clade"].unique()]

# ══════════════════════════════════════════════════════════════════════════════
# COVER PAGE
# ══════════════════════════════════════════════════════════════════════════════
def make_cover():
    fig = plt.figure(figsize=(8.27, 11.69))   # A4
    fig.patch.set_facecolor("#1B2A3B")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # Top colour band
    ax.add_patch(plt.Rectangle((0, 0.55), 1, 0.45, color="#2E6DA4", zorder=1))
    ax.add_patch(plt.Rectangle((0, 0.545), 1, 0.01, color="#E86E2D", zorder=2))

    # Title
    ax.text(0.5, 0.82, "TE Landscape Analysis",
            ha="center", va="center", color="white",
            fontsize=32, fontweight="bold", zorder=3)
    ax.text(0.5, 0.73, "Key Findings Report",
            ha="center", va="center", color="#D6E8F7",
            fontsize=18, zorder=3)
    ax.text(0.5, 0.67, "Early Diverging Fungi & Holozoa — 238 species",
            ha="center", va="center", color="#A0BDD8",
            fontsize=13, style="italic", zorder=3)

    # Sections box
    sections = [
        ("1", "#E86E2D", "Kimura Distance & TE Activity",    "Recency Index, burst detection, activity age"),
        ("2", "#2E8B57", "Shannon TE Diversity",             "Superfamily richness and evenness across clades"),
        ("3", "#7B2D8B", "Microsynteny & TE Influence",      "Gene neighbourhood conservation vs TE load"),
        ("4", "#1A6B9F", "Nucleosome Insertion Bias",        "NuclPred occupancy at TE insertion sites"),
    ]
    y0 = 0.46
    for num, color, title, sub in sections:
        ax.add_patch(plt.Rectangle((0.06, y0-0.005), 0.88, 0.085, color="#203245", zorder=3, clip_on=False))
        ax.add_patch(plt.Rectangle((0.06, y0-0.005), 0.07, 0.085, color=color, zorder=4, clip_on=False))
        ax.text(0.095, y0+0.035, num, ha="center", va="center", color="white",
                fontsize=20, fontweight="bold", zorder=5)
        ax.text(0.155, y0+0.052, title, ha="left", va="center", color="white",
                fontsize=12, fontweight="bold", zorder=5)
        ax.text(0.155, y0+0.022, sub, ha="left", va="center", color="#A0BDD8",
                fontsize=9, style="italic", zorder=5)
        y0 -= 0.105

    ax.text(0.5, 0.05, "April 2026",
            ha="center", va="center", color="#A0BDD8", fontsize=12, zorder=3)
    return fig_to_bytes(fig)

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — KIMURA / RECENCY INDEX
# ══════════════════════════════════════════════════════════════════════════════
def make_kimura_concept():
    """Conceptual figure explaining Kimura distance, RI and 'burst'."""
    fig = plt.figure(figsize=(11, 7))
    fig.patch.set_facecolor("white")

    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.12, left=0.04, right=0.97,
                           top=0.88, bottom=0.1)

    # ── Panel A: Kimura histogram concept ────────────────────────────────────
    ax = fig.add_subplot(gs[0])
    ax.set_facecolor("#F8F9FA")

    # Simulate a "bursty" species and an "old" species
    np.random.seed(42)
    x = np.linspace(0, 45, 400)
    # Recent burst (many copies at low div)
    burst = (0.55 * np.exp(-((x-2)**2)/(2*1.5**2)) +
             0.15 * np.exp(-((x-12)**2)/(2*4**2)) +
             0.03 * np.exp(-((x-30)**2)/(2*5**2)))
    # Old/silenced species (spread distribution)
    old   = (0.08 * np.exp(-((x-18)**2)/(2*7**2)) +
             0.05 * np.exp(-((x-30)**2)/(2*5**2)) +
             0.025)

    ax.fill_between(x, burst, alpha=0.45, color="#D04040", label="TE-active species (high RI)")
    ax.plot(x, burst, color="#D04040", lw=2)
    ax.fill_between(x, old, alpha=0.35, color="#3B7FBF", label="TE-silent species (low RI)")
    ax.plot(x, old, color="#3B7FBF", lw=2)

    # RI threshold line
    ax.axvline(5, color="#E86E2D", lw=2, ls="--", zorder=5)
    ax.annotate("Kimura < 5\n(recent copies)\n→ counts for RI",
                xy=(5, 0.42), xytext=(12, 0.52),
                fontsize=8.5, color="#E86E2D",
                arrowprops=dict(arrowstyle="->", color="#E86E2D", lw=1.5))

    # Burst arrow
    ax.annotate("", xy=(1.5, 0.48), xytext=(1.5, 0.02),
                arrowprops=dict(arrowstyle="->", color="#D04040", lw=2))
    ax.text(2.2, 0.25, "Recent\nburst", color="#D04040", fontsize=8.5, fontweight="bold")

    ax.set_xlabel("Kimura substitution distance (CpG-corrected)", fontsize=11)
    ax.set_ylabel("Proportion of TE copies", fontsize=11)
    ax.set_title("A   Kimura Distance: TE Age Distribution", fontsize=12, fontweight="bold", loc="left")
    ax.legend(fontsize=9, framealpha=0.8, loc="upper right")
    ax.set_xlim(0, 45)
    ax.set_ylim(0, 0.68)
    ax.spines[["top","right"]].set_visible(False)

    # RI formula box
    from matplotlib.patches import FancyBboxPatch
    box = FancyBboxPatch((22, 0.48), 20, 0.16,
                         boxstyle="round,pad=0.5", linewidth=1.5,
                         edgecolor="#E86E2D", facecolor="#FFF3EC")
    ax.add_patch(box)
    ax.text(32, 0.60, "Recency Index (RI)", ha="center", va="center",
            fontsize=9, fontweight="bold", color="#E86E2D")
    ax.text(32, 0.535,
            r"$RI = \frac{\sum bp_{K<5}}{\sum bp_{total}}$",
            ha="center", va="center", fontsize=10, color="#333333")

    # ── Panel B: Timeline diagram ──────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[1])
    ax2.set_facecolor("white")
    ax2.set_xlim(0, 10); ax2.set_ylim(0, 10); ax2.axis("off")
    ax2.set_title("B   From Kimura Distance to Evolutionary Time", fontsize=12,
                  fontweight="bold", loc="left")

    # Timeline arrow
    ax2.annotate("", xy=(9.5, 5), xytext=(0.5, 5),
                 arrowprops=dict(arrowstyle="-|>", color="#333333", lw=2))
    ax2.text(5, 4.2, "Evolutionary time  →", ha="center", fontsize=10, color="#555555")
    ax2.text(0.5, 4.55, "Present", fontsize=9, color="#333333", ha="left")
    ax2.text(9.2,4.55, "Ancient", fontsize=9, color="#333333", ha="right")

    # Kimura axis below
    for k, x in zip(["K=0","K=5","K=15","K=30","K=45"], [0.8, 2.0, 4.5, 7.0, 9.2]):
        ax2.plot(x, 4.95, "v", color="#888888", ms=6)
        ax2.text(x, 4.6, k, ha="center", fontsize=7.5, color="#888888")

    # RI zone highlight
    ax2.add_patch(plt.Rectangle((0.5, 5.1), 1.6, 1.0, color="#FFDDC0", alpha=0.7, zorder=2))
    ax2.text(1.3, 5.85, "RI zone\n(K < 5)", ha="center", fontsize=8, color="#E86E2D", fontweight="bold")

    # Three scenarios
    scenarios = [
        (2.5, 7.5, "#D04040", "Recent\nburst", "Many copies at\nlow K → high RI"),
        (5.0, 7.5, "#888888", "Continuous\ninsertion", "Spread distribution\n→ moderate RI"),
        (7.5, 7.5, "#3B7FBF", "Ancient /\nsilenced", "Copies at high K\n→ low RI"),
    ]
    for x, y, col, label, desc in scenarios:
        ax2.add_patch(FancyBboxPatch((x-1.0, y-0.6), 2.0, 0.9,
                                    boxstyle="round,pad=0.1",
                                    facecolor=col, alpha=0.18,
                                    edgecolor=col, lw=1.5))
        ax2.text(x, y, label, ha="center", va="center",
                 fontsize=8.5, fontweight="bold", color=col)
        ax2.text(x, y-1.1, desc, ha="center", va="center",
                 fontsize=7.5, color="#555555", style="italic")
        # Arrow down to timeline
        ax2.annotate("", xy=(x, 5.1), xytext=(x, y-0.65),
                     arrowprops=dict(arrowstyle="-|>", color=col,
                                    lw=1.2, connectionstyle="arc3,rad=0"))

    # TE class burst examples
    ax2.text(5, 1.5,
             "Burst score = peak height at K < 5 relative to background\n"
             "(higher burst score = stronger evidence of recent expansion)",
             ha="center", fontsize=9, color="#444444",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#EEF4FF",
                       edgecolor="#3B7FBF", lw=1.2))

    fig.suptitle("Section 1 — Kimura Distance, Recency Index & TE Bursts",
                 fontsize=14, fontweight="bold", y=0.97, color="#1B2A3B")
    return fig_to_bytes(fig)


def make_kimura_results():
    """RI per clade boxplot + scatter RI vs TE%."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.patch.set_facecolor("white")

    # ── Panel A: RI per clade ─────────────────────────────────────────────────
    ax = axes[0]
    data_by_clade = [merged.loc[merged["clade"]==c, "recency_index"].dropna().values
                     for c in clade_order]
    colors = [CLADE_COLORS.get(c, "#888888") for c in clade_order]
    labels = [CLADE_NAMES.get(c, c) for c in clade_order]

    bp = ax.boxplot(data_by_clade, patch_artist=True, notch=False,
                    medianprops=dict(color="black", lw=2),
                    whiskerprops=dict(lw=1.2),
                    capprops=dict(lw=1.2),
                    flierprops=dict(marker="o", ms=3, alpha=0.5))
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col); patch.set_alpha(0.7)

    # Overlay dots
    for i, (vals, col) in enumerate(zip(data_by_clade, colors)):
        jitter = np.random.RandomState(i).uniform(-0.22, 0.22, len(vals))
        ax.scatter(np.full(len(vals), i+1) + jitter, vals,
                   color=col, alpha=0.5, s=18, zorder=3, edgecolors="none")

    ax.set_xticks(range(1, len(clade_order)+1))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=9)
    ax.set_ylabel("Recency Index (RI)", fontsize=11)
    ax.set_title("A   RI Distribution per Clade", fontsize=12,
                 fontweight="bold", loc="left")
    ax.spines[["top","right"]].set_visible(False)
    ax.set_facecolor("#F8F9FA")

    # Annotation
    ax.text(0.98, 0.97, "Higher RI → more\nrecent TE activity",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color="#E86E2D",
            bbox=dict(facecolor="#FFF3EC", edgecolor="#E86E2D",
                      boxstyle="round,pad=0.3"))

    # ── Panel B: RI vs total TE% ──────────────────────────────────────────────
    ax2 = axes[1]
    for clade in clade_order:
        sub = merged[merged["clade"] == clade]
        ax2.scatter(sub["total_TE_pct"], sub["recency_index"],
                    color=CLADE_COLORS.get(clade, "#888"), alpha=0.7,
                    s=40, label=CLADE_NAMES.get(clade, clade),
                    edgecolors="white", lw=0.4)

    # Trend line
    xy = merged[["total_TE_pct","recency_index"]].dropna()
    if len(xy) > 10:
        z = np.polyfit(xy["total_TE_pct"], xy["recency_index"], 1)
        xfit = np.linspace(xy["total_TE_pct"].min(), xy["total_TE_pct"].max(), 100)
        ax2.plot(xfit, np.poly1d(z)(xfit), "k--", lw=1.5, alpha=0.6)
        r, p = spearmanr(xy["total_TE_pct"], xy["recency_index"])
        ax2.text(0.97, 0.05, f"ρ = {r:.2f}, p = {p:.1e}",
                 transform=ax2.transAxes, ha="right", fontsize=9,
                 bbox=dict(facecolor="white", edgecolor="#CCCCCC",
                           boxstyle="round,pad=0.3"))

    ax2.set_xlabel("Total TE content (%)", fontsize=11)
    ax2.set_ylabel("Recency Index (RI)", fontsize=11)
    ax2.set_title("B   RI vs Total TE Content", fontsize=12,
                  fontweight="bold", loc="left")
    ax2.legend(fontsize=7.5, ncol=2, loc="upper left",
               framealpha=0.8, markerscale=1.2)
    ax2.spines[["top","right"]].set_visible(False)
    ax2.set_facecolor("#F8F9FA")

    fig.suptitle("Section 1 — Kimura Distance & TE Activity", fontsize=13,
                 fontweight="bold", y=1.01, color="#1B2A3B")
    fig.tight_layout()
    return fig_to_bytes(fig)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — SHANNON DIVERSITY
# ══════════════════════════════════════════════════════════════════════════════
def make_shannon_results():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.patch.set_facecolor("white")

    # ── Panel A: H' per clade ────────────────────────────────────────────────
    ax = axes[0]
    data_by_clade = [merged.loc[merged["clade"]==c, "H_prime"].dropna().values
                     for c in clade_order]
    colors = [CLADE_COLORS.get(c, "#888888") for c in clade_order]
    bp = ax.boxplot(data_by_clade, patch_artist=True,
                    medianprops=dict(color="black", lw=2),
                    whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2),
                    flierprops=dict(marker="o", ms=3, alpha=0.5))
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col); patch.set_alpha(0.7)
    for i, (vals, col) in enumerate(zip(data_by_clade, colors)):
        jitter = np.random.RandomState(i+100).uniform(-0.2, 0.2, len(vals))
        ax.scatter(np.full(len(vals), i+1)+jitter, vals,
                   color=col, alpha=0.5, s=18, zorder=3, edgecolors="none")
    ax.set_xticks(range(1, len(clade_order)+1))
    ax.set_xticklabels([CLADE_NAMES.get(c, c) for c in clade_order],
                       rotation=40, ha="right", fontsize=9)
    ax.set_ylabel("Shannon H' (TE superfamily diversity)", fontsize=11)
    ax.set_title("A   TE Diversity (H') per Clade", fontsize=12,
                 fontweight="bold", loc="left")
    ax.spines[["top","right"]].set_visible(False)
    ax.set_facecolor("#F8F9FA")

    box_text = ("High H' = many superfamilies\n"
                "with similar abundance\n"
                "(diverse TE landscape)")
    ax.text(0.98, 0.97, box_text, transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, color="#1A6B45",
            bbox=dict(facecolor="#EEFFF4", edgecolor="#1A6B45",
                      boxstyle="round,pad=0.3"))

    # ── Panel B: H' vs TE% ───────────────────────────────────────────────────
    ax2 = axes[1]
    for clade in clade_order:
        sub = merged[merged["clade"] == clade]
        ax2.scatter(sub["total_TE_pct"], sub["H_prime"],
                    color=CLADE_COLORS.get(clade, "#888"), alpha=0.7,
                    s=40, label=CLADE_NAMES.get(clade, clade),
                    edgecolors="white", lw=0.4)

    xy = merged[["total_TE_pct","H_prime"]].dropna()
    if len(xy) > 10:
        z = np.polyfit(xy["total_TE_pct"], xy["H_prime"], 1)
        xfit = np.linspace(0, xy["total_TE_pct"].max(), 100)
        ax2.plot(xfit, np.poly1d(z)(xfit), "k--", lw=1.5, alpha=0.6)
        r, p = spearmanr(xy["total_TE_pct"], xy["H_prime"])
        ax2.text(0.97, 0.05, f"ρ = {r:.2f}, p = {p:.1e}",
                 transform=ax2.transAxes, ha="right", fontsize=9,
                 bbox=dict(facecolor="white", edgecolor="#CCCCCC",
                           boxstyle="round,pad=0.3"))

    ax2.set_xlabel("Total TE content (%)", fontsize=11)
    ax2.set_ylabel("Shannon H'", fontsize=11)
    ax2.set_title("B   H' vs Total TE Content", fontsize=12,
                  fontweight="bold", loc="left")
    ax2.legend(fontsize=7.5, ncol=2, loc="upper left", framealpha=0.8)
    ax2.spines[["top","right"]].set_visible(False)
    ax2.set_facecolor("#F8F9FA")

    fig.suptitle("Section 2 — Shannon Diversity of TE Superfamilies",
                 fontsize=13, fontweight="bold", y=1.01, color="#1B2A3B")
    fig.tight_layout()
    return fig_to_bytes(fig)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — MICROSYNTENY CONCEPT + RESULTS
# ══════════════════════════════════════════════════════════════════════════════
def make_synteny_concept():
    """Diagram explaining how microsynteny is measured."""
    fig = plt.figure(figsize=(11, 6.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 11); ax.set_ylim(0, 7); ax.axis("off")

    fig.suptitle("Section 3 — Microsynteny: How It Is Measured",
                 fontsize=14, fontweight="bold", y=0.97, color="#1B2A3B")

    # ── Species A (conserved) ─────────────────────────────────────────────────
    y1 = 5.2
    ax.text(0.3, y1+0.35, "Species A", fontsize=11, fontweight="bold", color="#1B2A3B")
    ax.text(0.3, y1+0.05, "(reference genome)", fontsize=8.5, color="#666")

    genes_A = [(1.2, "G1"), (2.2, "G2"), (3.2, "G3"), (4.2, "G4"), (5.2, "G5"),
               (6.2, "G6"), (7.2, "G7"), (8.2, "G8"), (9.2, "G9")]
    # Chromosome line
    ax.plot([0.9, 10.1], [y1, y1], color="#888888", lw=3, solid_capstyle="round")
    for x, name in genes_A:
        ax.add_patch(plt.Rectangle((x, y1-0.22), 0.65, 0.44,
                                   color="#3B7FBF", alpha=0.85, zorder=3))
        ax.text(x+0.325, y1, name, ha="center", va="center",
                fontsize=8, color="white", fontweight="bold", zorder=4)

    # Sliding window bracket
    wx = 2.1; ww = 3.5
    ax.add_patch(mpatches.FancyBboxPatch((wx, y1-0.55), ww, 1.1,
                                         boxstyle="round,pad=0.1",
                                         edgecolor="#E86E2D", facecolor="#FFF3EC",
                                         lw=2, alpha=0.5, zorder=2))
    ax.text(wx + ww/2, y1-0.72, "sliding window (k genes)", ha="center",
            fontsize=8.5, color="#E86E2D", fontweight="bold")

    # ── Species B (conserved neighbour) ──────────────────────────────────────
    y2 = 3.4
    ax.text(0.3, y2+0.35, "Species B", fontsize=11, fontweight="bold", color="#1B2A3B")
    ax.text(0.3, y2+0.05, "(query genome — conserved)", fontsize=8.5, color="#1A6B45")
    ax.plot([0.9, 10.1], [y2, y2], color="#888888", lw=3, solid_capstyle="round")

    # Same order, synteny conserved
    genes_B = [(1.2, "G1"), (2.2, "G2"), (3.2, "G3"), (4.2, "G4"), (5.2, "G5"),
               (6.2, "G6"), (7.2, "G7"), (8.2, "G8"), (9.2, "G9")]
    for x, name in genes_B:
        ax.add_patch(plt.Rectangle((x, y2-0.22), 0.65, 0.44,
                                   color="#1A6B45", alpha=0.85, zorder=3))
        ax.text(x+0.325, y2, name, ha="center", va="center",
                fontsize=8, color="white", fontweight="bold", zorder=4)

    # Match arrows (conserved)
    for xA, xB in [(2.525, 2.525), (3.525, 3.525), (4.525, 4.525), (5.525, 5.525)]:
        ax.annotate("", xy=(xB, y2+0.22), xytext=(xA, y1-0.22),
                    arrowprops=dict(arrowstyle="-", color="#1A6B45",
                                   lw=1, alpha=0.4))
    ax.add_patch(mpatches.FancyBboxPatch((wx, y2-0.55), ww, 1.1,
                                         boxstyle="round,pad=0.1",
                                         edgecolor="#1A6B45", facecolor="#EEFFF4",
                                         lw=2, alpha=0.5, zorder=2))
    ax.text(4.0, y2-0.72, "✓ collinear pairs found → r = 1.0", ha="center",
            fontsize=8.5, color="#1A6B45", fontweight="bold")

    # ── Species C (TE disrupted) ──────────────────────────────────────────────
    y3 = 1.5
    ax.text(0.3, y3+0.35, "Species C", fontsize=11, fontweight="bold", color="#1B2A3B")
    ax.text(0.3, y3+0.05, "(query genome — TE disrupted)", fontsize=8.5, color="#D04040")
    ax.plot([0.9, 10.1], [y3, y3], color="#888888", lw=3, solid_capstyle="round")

    # Rearranged + TE insertions
    genes_C = [(1.2, "G1"), (2.5, "TE"), (3.5, "TE"), (4.5, "G5"), (5.5, "TE"),
               (6.5, "G3"), (7.5, "G9"), (8.5, "TE"), (9.5, "G2")]
    for x, name in genes_C:
        col = "#E86E2D" if name == "TE" else "#D04040"
        ax.add_patch(plt.Rectangle((x, y3-0.22), 0.6 if name!="TE" else 0.8, 0.44,
                                   color=col, alpha=0.85, zorder=3))
        ax.text(x+0.3, y3, name, ha="center", va="center",
                fontsize=7.5, color="white", fontweight="bold", zorder=4)

    ax.add_patch(mpatches.FancyBboxPatch((wx, y3-0.55), ww, 1.1,
                                         boxstyle="round,pad=0.1",
                                         edgecolor="#D04040", facecolor="#FFF0EE",
                                         lw=2, alpha=0.5, zorder=2))
    ax.text(4.0, y3-0.72, "✗ few collinear pairs → r ≈ 0.0", ha="center",
            fontsize=8.5, color="#D04040", fontweight="bold")

    # ── Formula box ────────────────────────────────────────────────────────────
    ax.add_patch(FancyBboxPatch((7.0, 0.4), 3.6, 2.2,
                                boxstyle="round,pad=0.2",
                                edgecolor="#7B2D8B", facecolor="#F5F0FF", lw=2))
    ax.text(8.8, 2.35, "Microsynteny Score (r)", ha="center",
            fontsize=10, fontweight="bold", color="#7B2D8B")
    ax.text(8.8, 1.9,
            r"$r_{AB} = \frac{|\text{collinear pairs}_{AB}|}{k}$",
            ha="center", fontsize=11, color="#333333")
    ax.text(8.8, 1.45, "k = window size", ha="center", fontsize=8.5, color="#666")
    ax.text(8.8, 1.1, "mean_r_intra = mean across\nall within-clade pairs",
            ha="center", va="center", fontsize=8.5, color="#555")
    ax.text(8.8, 0.65, "r = 1.0 → perfect conservation\nr = 0.0 → complete disruption",
            ha="center", va="center", fontsize=8.5, color="#555", style="italic")

    return fig_to_bytes(fig)


def make_synteny_results():
    """Scatter TE% vs synteny score, coloured by clade."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.patch.set_facecolor("white")

    # ── Panel A: scatter TE% vs mean_r_intra ──────────────────────────────────
    ax = axes[0]
    ax.set_facecolor("#F8F9FA")

    # Merge in n_partners
    syn_cols = synteny[["species","mean_r_intra","n_partners"]]
    merged_syn = merged.merge(syn_cols, on="species", how="left",
                              suffixes=("","_syn"))
    # Use mean_r_intra from synteny (may already be in merged, but ensure n_partners)
    if "mean_r_intra_syn" in merged_syn.columns:
        merged_syn["mean_r_intra"] = merged_syn["mean_r_intra_syn"]
    if "n_partners" not in merged_syn.columns:
        merged_syn["n_partners"] = 99

    # Species with few partners → artefact (r=1.0 trivially)
    MIN_PARTNERS = 5
    artefact = merged_syn[merged_syn["n_partners"] < MIN_PARTNERS]
    valid    = merged_syn[merged_syn["n_partners"] >= MIN_PARTNERS]

    # Draw artefact species as faded grey crosses first
    ax.scatter(artefact["total_TE_pct"], artefact["mean_r_intra"],
               color="#BBBBBB", alpha=0.4, s=35, marker="x", lw=1,
               zorder=2, label=f"< {MIN_PARTNERS} within-clade partners\n(r trivially = 1.0)")

    for clade in clade_order:
        sub = valid[(valid["clade"] == clade) & valid["mean_r_intra"].notna()
                    & valid["total_TE_pct"].notna()]
        ax.scatter(sub["total_TE_pct"], sub["mean_r_intra"],
                   color=CLADE_COLORS.get(clade, "#888"), alpha=0.8,
                   s=55, label=CLADE_NAMES.get(clade, clade),
                   edgecolors="white", lw=0.5, zorder=3)

    # Trend on valid data only
    xy_v = valid[["total_TE_pct","mean_r_intra"]].dropna()
    if len(xy_v) > 10:
        z = np.polyfit(xy_v["total_TE_pct"], xy_v["mean_r_intra"], 1)
        xfit = np.linspace(0, xy_v["total_TE_pct"].max(), 100)
        ax.plot(xfit, np.poly1d(z)(xfit), "k--", lw=2, alpha=0.7, zorder=4)
        r_val, p_val = spearmanr(xy_v["total_TE_pct"], xy_v["mean_r_intra"])
        sig = "**" if p_val < 0.01 else "*" if p_val < 0.05 else "(n.s.)"
        ax.text(0.97, 0.97,
                f"Spearman ρ = {r_val:.2f}\np = {p_val:.2e}  {sig}\n"
                f"(n = {len(xy_v)}, n_partners ≥ {MIN_PARTNERS})",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                color="#7B2D8B", fontweight="bold",
                bbox=dict(facecolor="#F5F0FF", edgecolor="#7B2D8B",
                          boxstyle="round,pad=0.4"))

    ax.set_xlabel("Total TE content (%)", fontsize=12)
    ax.set_ylabel("Mean microsynteny score (r_intra)", fontsize=12)
    ax.set_title("A   TE Content vs Microsynteny Conservation\n"
                 "(grey × = excluded: trivial r=1.0 in small clades)",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=7.5, ncol=2, loc="lower left", framealpha=0.85)
    ax.spines[["top","right"]].set_visible(False)

    # ── Panel B: mean_r_intra per clade boxplot ───────────────────────────────
    ax2 = axes[1]
    ax2.set_facecolor("#F8F9FA")
    clades_syn = [c for c in clade_order
                  if merged.loc[merged["clade"]==c, "mean_r_intra"].dropna().shape[0] > 1]
    data_syn = [merged.loc[merged["clade"]==c, "mean_r_intra"].dropna().values
                for c in clades_syn]
    colors2 = [CLADE_COLORS.get(c, "#888") for c in clades_syn]

    bp = ax2.boxplot(data_syn, patch_artist=True,
                     medianprops=dict(color="black", lw=2),
                     whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2),
                     flierprops=dict(marker="o", ms=3, alpha=0.5))
    for patch, col in zip(bp["boxes"], colors2):
        patch.set_facecolor(col); patch.set_alpha(0.7)
    for i, (vals, col) in enumerate(zip(data_syn, colors2)):
        jitter = np.random.RandomState(i+200).uniform(-0.2, 0.2, len(vals))
        ax2.scatter(np.full(len(vals), i+1)+jitter, vals,
                    color=col, alpha=0.5, s=18, zorder=3, edgecolors="none")

    ax2.set_xticks(range(1, len(clades_syn)+1))
    ax2.set_xticklabels([CLADE_NAMES.get(c, c) for c in clades_syn],
                        rotation=40, ha="right", fontsize=9)
    ax2.set_ylabel("Mean microsynteny score", fontsize=11)
    ax2.set_title("B   Microsynteny Conservation per Clade",
                  fontsize=12, fontweight="bold", loc="left")
    ax2.spines[["top","right"]].set_visible(False)

    fig.suptitle("Section 3 — Microsynteny & TE Influence",
                 fontsize=13, fontweight="bold", y=1.01, color="#1B2A3B")
    fig.tight_layout()
    return fig_to_bytes(fig)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — NUCLEOSOME INSERTION BIAS
# ══════════════════════════════════════════════════════════════════════════════

def compute_nuc_bias_streaming(species_list, max_lines=3_000_000, seed=42):
    """
    Streaming approach: read up to max_lines of each WIG file, collect all
    (position, score) pairs, then for every TE whose midpoint falls within
    the read region compute mean score ±500 bp and compare to background.

    max_lines=3M covers ~3 Mb at step=1, enough for a complete yeast genome
    or the first few scaffolds of larger fungi.
    """
    rng = random.Random(seed)
    records = []

    for sp in species_list:
        wig_path = WIG_DIR / f"{sp}_nuclpred.wig"
        bed_path = TE_ANN / f"{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
        if not wig_path.exists() or not bed_path.exists():
            continue
        clade = sp.split("_")[0]
        print(f"  {sp} ...", flush=True, end=" ")

        # ── Read WIG (streaming, stop at max_lines) ──────────────────────────
        # Store as {chrom: (start_pos, np.array of scores)}
        wig_store = {}
        chrom = None; wig_start_pos = 1; step = 1
        cur_vals = []; cur_start = 1; lines_read = 0

        def store_current():
            if chrom and cur_vals:
                wig_store[chrom] = (cur_start, np.array(cur_vals, dtype=np.float32))

        with open(wig_path) as fh:
            for raw in fh:
                if lines_read >= max_lines:
                    break
                raw = raw.rstrip()
                if not raw or raw.startswith("track") or raw.startswith("#"):
                    continue
                if raw.startswith("fixedStep"):
                    store_current(); cur_vals = []
                    m = re.search(r"chrom=(\S+)", raw)
                    raw_chrom = m.group(1) if m else "?"
                    # Strip species prefix: "AS_Sacce1.chrI" → "chrI"
                    chrom = raw_chrom.split(".", 1)[1] if "." in raw_chrom else raw_chrom
                    m = re.search(r"start=(\d+)", raw); cur_start = int(m.group(1)) if m else 1
                    m = re.search(r"step=(\d+)", raw); step = int(m.group(1)) if m else 1
                else:
                    try:
                        cur_vals.append(float(raw))
                        lines_read += 1
                    except ValueError:
                        pass
        store_current()

        if not wig_store:
            print("no WIG data"); continue

        # ── Match TEs to read regions ─────────────────────────────────────────
        tes = []
        with open(bed_path) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 4: continue
                ch, s, e, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                if ch not in wig_store: continue
                top = name.split("/")[0].split("#")[-1].split("_")[0]
                tes.append((ch, s, e, top))

        if not tes:
            print("no TEs in WIG region"); continue

        sampled = rng.sample(tes, min(150, len(tes)))

        def get_mean(ch, pos, half=300):
            if ch not in wig_store: return np.nan
            w_start, arr = wig_store[ch]
            i0 = max(0, (pos - half - w_start) // step)
            i1 = min(len(arr)-1, (pos + half - w_start) // step)
            if i1 < i0: return np.nan
            v = arr[i0:i1+1]
            return float(np.mean(v)) if len(v) > 3 else np.nan

        class_te = defaultdict(list)
        class_bg = defaultdict(list)

        for ch, s, e, cls in sampled:
            mid = (s + e) // 2
            vte = get_mean(ch, mid)
            shift = rng.choice([-1,1]) * rng.randint(8000, 40000)
            vbg = get_mean(ch, max(1, mid + shift))
            if not np.isnan(vte): class_te[cls].append(vte)
            if not np.isnan(vbg): class_bg[cls].append(vbg)

        all_te, all_bg = [], []
        for cls in class_te:
            te_v = class_te[cls]; bg_v = class_bg.get(cls, [])
            if len(te_v) < 3: continue
            mt = float(np.mean(te_v)); mb = float(np.mean(bg_v)) if bg_v else np.nan
            delta = mt - mb if not np.isnan(mb) else np.nan
            records.append({"species": sp, "clade": clade,
                            "te_class": cls, "mean_te": mt,
                            "mean_bg": mb, "delta": delta, "n": len(te_v)})
            all_te.extend(te_v); all_bg.extend(bg_v)

        if all_te and all_bg:
            records.append({"species": sp, "clade": clade, "te_class": "ALL",
                            "mean_te": float(np.mean(all_te)),
                            "mean_bg": float(np.mean(all_bg)),
                            "delta": float(np.mean(all_te)) - float(np.mean(all_bg)),
                            "n": len(all_te)})
            print(f"OK ({len(sampled)} TEs, Δ={float(np.mean(all_te))-float(np.mean(all_bg)):+.3f})")
        else:
            print("insufficient data")

    return pd.DataFrame(records)


def make_nucleosome_figure(nuc_df):
    """Summary figure: NuclPred score at TE vs background per clade and TE class."""
    if nuc_df.empty:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    fig.patch.set_facecolor("white")

    # ── Panel A: Δ occupancy (TE - bg) per clade ────────────────────────────
    ax = axes[0]
    ax.set_facecolor("#F8F9FA")

    clade_delta = (nuc_df[nuc_df["te_class"]=="ALL"]
                   .groupby("clade")["delta"].agg(["mean","sem","size"])
                   .reset_index()
                   .rename(columns={"mean":"delta_mean","sem":"delta_sem","size":"n"}))
    clade_delta = clade_delta[clade_delta["n"] >= 2].sort_values("delta_mean")

    y_pos = range(len(clade_delta))
    colors_bar = [CLADE_COLORS.get(c, "#888") for c in clade_delta["clade"]]
    ax.barh(list(y_pos), clade_delta["delta_mean"],
            xerr=clade_delta["delta_sem"],
            color=colors_bar, alpha=0.8, edgecolor="white",
            error_kw=dict(ecolor="black", lw=1.5, capsize=4))
    ax.axvline(0, color="black", lw=1.5, ls="--")
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels([CLADE_NAMES.get(c, c) for c in clade_delta["clade"]],
                       fontsize=9)
    ax.set_xlabel("Δ NuclPred score\n(TE insertion sites − genome background)", fontsize=10.5)
    ax.set_title("A   Nucleosome Occupancy at TE Insertions\nvs Genome Background",
                 fontsize=11, fontweight="bold", loc="left")
    ax.spines[["top","right"]].set_visible(False)

    # Annotation
    ax.text(0.97, 0.04,
            "Δ < 0 → TEs insert into\nnucleosome-depleted regions (NDRs)\n"
            "Δ > 0 → TEs insert into\nnucleosome-enriched regions",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color="#555",
            bbox=dict(facecolor="#EEF4FF", edgecolor="#3B7FBF",
                      boxstyle="round,pad=0.3"))

    # ── Panel B: Δ per TE class (all species averaged) ────────────────────────
    ax2 = axes[1]
    ax2.set_facecolor("#F8F9FA")

    class_delta = (nuc_df[nuc_df["te_class"]!="ALL"]
                   .groupby("te_class")["delta"].agg(["mean","sem","size"])
                   .reset_index()
                   .rename(columns={"mean":"delta_mean","sem":"delta_sem","size":"n"}))
    class_delta = class_delta[class_delta["n"] >= 3].sort_values("delta_mean")

    # Colours per class
    class_colors = {
        "LTR": "#D04040", "LINE": "#E07B39", "SINE": "#E0A020",
        "DNA": "#3B7FBF", "RC": "#1A6B45", "Rolling_Circle": "#1A6B45",
        "Penelope": "#7B2D8B", "Unknown": "#888888",
        "Unclassified": "#AAAAAA", "MITE": "#4B7A2E",
    }
    y2 = range(len(class_delta))
    cb2 = [class_colors.get(c, "#888888") for c in class_delta["te_class"]]
    ax2.barh(list(y2), class_delta["delta_mean"],
             xerr=class_delta["delta_sem"],
             color=cb2, alpha=0.8, edgecolor="white",
             error_kw=dict(ecolor="black", lw=1.5, capsize=4))
    ax2.axvline(0, color="black", lw=1.5, ls="--")
    ax2.set_yticks(list(y2))
    ax2.set_yticklabels([c.replace("_"," ") for c in class_delta["te_class"]],
                        fontsize=9)
    ax2.set_xlabel("Δ NuclPred score (TE − background)", fontsize=10.5)
    ax2.set_title("B   Nucleosome Bias by TE Class\n(averaged across all species)",
                  fontsize=11, fontweight="bold", loc="left")
    ax2.spines[["top","right"]].set_visible(False)

    n_sp = nuc_df["species"].nunique()
    fig.suptitle(f"Section 4 — Nucleosome Occupancy at TE Insertion Sites  (n = {n_sp} species)",
                 fontsize=13, fontweight="bold", y=1.01, color="#1B2A3B")
    fig.tight_layout()
    return fig_to_bytes(fig)


def make_nucleosome_concept():
    """Conceptual figure explaining nucleosome insertion bias."""
    fig = plt.figure(figsize=(11, 5.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.03, 0.05, 0.94, 0.85])
    ax.set_xlim(0, 11); ax.set_ylim(0, 6); ax.axis("off")

    fig.suptitle("Section 4 — Nucleosome Insertion Bias: Concept",
                 fontsize=14, fontweight="bold", y=0.99, color="#1B2A3B")

    # DNA strand
    y_dna = 3.8
    ax.plot([0.3, 10.7], [y_dna, y_dna], color="#BBBBBB", lw=4,
            solid_capstyle="round", zorder=1)
    ax.text(0.1, y_dna, "Genome", fontsize=9, color="#888", va="center", rotation=90)

    # Nucleosomes (circles on DNA)
    nuc_positions = [1.0, 2.2, 3.4, 5.8, 7.0, 8.2, 9.4]
    for x in nuc_positions:
        circ = plt.Circle((x, y_dna), 0.38, color="#4472C4", alpha=0.55, zorder=2)
        ax.add_patch(circ)
        ax.text(x, y_dna, "NUC", ha="center", va="center",
                fontsize=6, color="white", fontweight="bold", zorder=3)

    # NDR (nucleosome-depleted region)
    ndr_x = 4.6
    ax.add_patch(plt.Rectangle((4.0, y_dna-0.5), 1.4, 1.0,
                                color="#E86E2D", alpha=0.12, zorder=2))
    ax.text(ndr_x, y_dna+0.7, "NDR", ha="center", fontsize=10,
            fontweight="bold", color="#E86E2D")
    ax.text(ndr_x, y_dna+0.42, "(Nucleosome-\nDepleted Region)", ha="center",
            fontsize=7.5, color="#E86E2D")

    # NuclPred score profile
    ax_inset = fig.add_axes([0.1, 0.06, 0.55, 0.35])
    x_prof = np.linspace(0, 10, 500)
    score = np.zeros(500)
    # Nucleosome peaks
    for xn in nuc_positions:
        score += 0.6 * np.exp(-((x_prof - xn)**2) / (2 * 0.25**2))
    # NDR dip
    score -= 0.3 * np.exp(-((x_prof - 4.6)**2) / (2 * 0.3**2))
    score = np.clip(score, 0, None)
    ax_inset.fill_between(x_prof, score, alpha=0.4, color="#4472C4")
    ax_inset.plot(x_prof, score, color="#4472C4", lw=1.5)
    ax_inset.axvspan(4.0, 5.2, alpha=0.15, color="#E86E2D")
    ax_inset.set_xlabel("Genome position", fontsize=9)
    ax_inset.set_ylabel("NuclPred score", fontsize=9)
    ax_inset.set_title("NuclPred occupancy profile", fontsize=9, fontweight="bold")
    ax_inset.spines[["top","right"]].set_visible(False)

    # TE insertion arrows
    # Into NDR
    ax.annotate("", xy=(4.6, y_dna+0.02), xytext=(4.6, y_dna+1.7),
                arrowprops=dict(arrowstyle="-|>", color="#D04040", lw=2.5))
    ax.text(4.6, y_dna+1.95, "TE inserting\ninto NDR\n(Δ < 0)", ha="center",
            fontsize=8.5, color="#D04040", fontweight="bold")

    # Into nucleosome
    ax.annotate("", xy=(1.0, y_dna+0.4), xytext=(2.5, y_dna+1.7),
                arrowprops=dict(arrowstyle="-|>", color="#1A6B45", lw=2.5))
    ax.text(2.0, y_dna+1.95, "TE inserting\ninto NUC\n(Δ > 0)", ha="center",
            fontsize=8.5, color="#1A6B45", fontweight="bold")

    # Method box
    ax.add_patch(FancyBboxPatch((7.2, 0.3), 3.5, 3.2,
                                boxstyle="round,pad=0.2",
                                edgecolor="#1A6B9F", facecolor="#EEF4FF", lw=1.5))
    ax.text(8.95, 3.2, "Method", ha="center", fontsize=10,
            fontweight="bold", color="#1A6B9F")
    steps = [
        "1. For each TE: extract NuclPred",
        "   scores at insertion midpoint ± w",
        "2. For each TE: sample random",
        "   background region (same chrom)",
        "3. Δ = mean(TE scores) − mean(BG)",
        "4. Negative Δ → preference for NDRs",
        "5. Summarise by TE class & clade",
    ]
    for i, s in enumerate(steps):
        ax.text(7.4, 2.85 - i*0.38, s, fontsize=8, color="#333")

    return fig_to_bytes(fig)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    writer = PdfWriter()

    print("\n[Cover]")
    append_bytes(writer, make_cover())

    print("\n[Section 1 — Kimura concept]")
    append_bytes(writer, make_kimura_concept())

    print("[Section 1 — Kimura results]")
    append_bytes(writer, make_kimura_results())

    print("\n[Section 2 — Shannon]")
    append_bytes(writer, make_shannon_results())

    print("\n[Section 3 — Microsynteny concept]")
    append_bytes(writer, make_synteny_concept())

    print("[Section 3 — Microsynteny results]")
    append_bytes(writer, make_synteny_results())

    print("[Section 3 — MU synteny heatmap]")
    append_file(writer, SYN_DIR/"figures_filtered_test/MU_synteny_heatmap.pdf")

    # ── Section 4: Nucleosomes ─────────────────────────────────────────────────
    print("\n[Section 4 — Nucleosome concept]")
    append_bytes(writer, make_nucleosome_concept())

    # Pick representative species (max 3 per clade, prefer smaller genomes = faster WIG)
    wig_species = {f.replace("_nuclpred.wig","")
                   for f in os.listdir(WIG_DIR) if f.endswith(".wig")}
    bed_species = {p.parent.parent.name.replace("_EarlGrey","")
                   for p in Path(TE_ANN).glob("*_EarlGrey/*_summaryFiles/*.filteredRepeats.bed")}
    available = sorted(wig_species & bed_species)

    # Sort by genome size ASC (smaller = faster WIG), pick top 2 per clade
    gsize = outliers.set_index("species")["genome_size_mb"].to_dict()
    te_pct = outliers.set_index("species")["total_TE_pct"].to_dict()
    by_clade = defaultdict(list)
    for sp in available:
        clade = sp.split("_")[0]
        by_clade[clade].append((gsize.get(sp, 9999), sp))
    selected = []
    for clade, sps in sorted(by_clade.items()):
        sps.sort()           # ascending genome size
        selected.extend([sp for _, sp in sps[:2]])

    print(f"  Running nucleosome bias for {len(selected)} species:")
    for sp in selected:
        print(f"    {sp}: {gsize.get(sp,'?'):.1f} Mb, {te_pct.get(sp,'?'):.1f}% TE")

    nuc_df = compute_nuc_bias_streaming(selected, max_lines=2_000_000)
    print(f"  Result: {len(nuc_df)} records, {nuc_df['species'].nunique() if not nuc_df.empty else 0} species")

    if not nuc_df.empty:
        nuc_df.to_csv(SEC2/"nucleosome_bias_summary.tsv", sep="\t", index=False)
        print("[Section 4 — Nucleosome results]")
        fig_bytes = make_nucleosome_figure(nuc_df)
        if fig_bytes:
            append_bytes(writer, fig_bytes)
    else:
        print("  [!] No nucleosome data — skipping results figure")

    print(f"\nWriting {OUT_PDF} ...")
    with open(OUT_PDF, "wb") as f:
        writer.write(f)
    print(f"Done — {OUT_PDF.stat().st_size/1024:.0f} KB")
    print(f"Pages: {len(writer.pages)}")


if __name__ == "__main__":
    main()
