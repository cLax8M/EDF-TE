#!/usr/bin/env python3
"""
C02_kimura_figures.py
=====================
Generate Kimura divergence figures from C01 output tables.

Figures produced (all PDF):
  1. kimura_ri_by_clade_boxplot.pdf
        — Boxplot of Recency Index per clade, coloured by clade, points overlaid.
  2. kimura_ri_vs_burst.pdf
        — Scatter RI vs burst_score, coloured by clade.
  3. kimura_ri_vs_mean_div.pdf
        — Scatter RI vs mean_div_weighted, coloured by clade.
  4. kimura_mean_div_by_clade.pdf
        — Boxplot of mean_div_weighted per clade.
  5. kimura_clade_heatmap.pdf
        — Heatmap: clade × metric (RI, burst, mean_div) with means.
  6. kimura_landscape_stacked.pdf
        — Representative Kimura landscape (stacked area by TE class) for
          one species per clade (species with median TE content in clade).
  7. kimura_ri_vs_shannon.pdf
        — Scatter RI (Kimura activity) vs H' (Shannon diversity).

Usage:
    python C02_kimura_figures.py [FIGURES_DIR]
"""

import sys, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

warnings.filterwarnings("ignore")

# ── paths ─────────────────────────────────────────────────────────────────────
RESULTS = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape"
FIG_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RESULTS, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

KIMURA_SP  = os.path.join(RESULTS, "C01_kimura_species.tsv")
KIMURA_SF  = os.path.join(RESULTS, "C01_kimura_superfamily.tsv")
SHANNON    = os.path.join(RESULTS, "B01_shannon_diversity.tsv")
EARLGREY   = "/data/carlos/EDF-TE/TE_analysis/annotated"

# ── colour palette per clade (canonical, single source) ─────────────────────
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from clade_palette import CLADE_COLORS as CLADE_COLOURS, get_clade  # noqa: E402

def clade_of(sp):
    return get_clade(sp)


def load_data():
    sp = pd.read_csv(KIMURA_SP, sep="\t")
    sp["clade"] = sp["species"].apply(clade_of)
    sf = pd.read_csv(KIMURA_SF, sep="\t")
    sf["clade"] = sf["species"].apply(clade_of)
    sh = pd.read_csv(SHANNON, sep="\t")
    sh["clade"] = sh["species"].apply(clade_of)
    return sp, sf, sh


# ── Figure 1: RI by clade boxplot ────────────────────────────────────────────
def fig_ri_by_clade(sp):
    from clade_palette import CLADE_ORDER as _PHYLO_ORDER, CLADE_NAMES
    present = set(sp["clade"].unique())
    clades_sorted = [c for c in _PHYLO_ORDER if c in present] + \
                    sorted(present - set(_PHYLO_ORDER))
    fig, ax = plt.subplots(figsize=(12, 5))
    data = [sp.loc[sp.clade == c, "recency_index"].values for c in clades_sorted]
    bp = ax.boxplot(data, labels=clades_sorted, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", linewidth=1.5))
    for i, c in enumerate(clades_sorted):
        col = CLADE_COLOURS.get(c, "#888888")
        bp["boxes"][i].set_facecolor(col)
        bp["boxes"][i].set_alpha(0.5)
        vals = sp.loc[sp.clade == c, "recency_index"].values
        jitter = np.random.normal(0, 0.06, len(vals))
        ax.scatter(np.full_like(vals, i+1) + jitter, vals, c=col, s=18, alpha=0.7, edgecolors="none", zorder=3)
    ax.set_ylabel("Recency Index (RI)")
    ax.set_xlabel("Clade")
    ax.set_title("Kimura Recency Index by Clade")
    ax.axhline(0.1, ls="--", c="grey", alpha=0.5, lw=0.8)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "kimura_ri_by_clade_boxplot.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 2: RI vs burst_score scatter ──────────────────────────────────────
def fig_ri_vs_burst(sp):
    fig, ax = plt.subplots(figsize=(8, 6))
    for c in sorted(sp.clade.unique()):
        sub = sp[sp.clade == c]
        ax.scatter(sub.recency_index, sub.burst_score, label=c,
                   c=CLADE_COLOURS.get(c, "#888888"), s=25, alpha=0.7, edgecolors="none")
    ax.set_xlabel("Recency Index (RI)")
    ax.set_ylabel("Burst Score")
    ax.set_title("Kimura: Recency Index vs Burst Score")
    ax.legend(fontsize=6, ncol=3, loc="best", framealpha=0.5)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "kimura_ri_vs_burst.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 3: RI vs mean_div scatter ─────────────────────────────────────────
def fig_ri_vs_meandiv(sp):
    fig, ax = plt.subplots(figsize=(8, 6))
    for c in sorted(sp.clade.unique()):
        sub = sp[sp.clade == c]
        ax.scatter(sub.recency_index, sub.mean_div_weighted, label=c,
                   c=CLADE_COLOURS.get(c, "#888888"), s=25, alpha=0.7, edgecolors="none")
    ax.set_xlabel("Recency Index (RI)")
    ax.set_ylabel("Mean Divergence (weighted)")
    ax.set_title("Kimura: Recency Index vs Mean Divergence")
    ax.legend(fontsize=6, ncol=3, loc="best", framealpha=0.5)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "kimura_ri_vs_mean_div.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 4: mean_div by clade boxplot ──────────────────────────────────────
def fig_meandiv_by_clade(sp):
    clades_sorted = sp.groupby("clade")["mean_div_weighted"].median().sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(12, 5))
    data = [sp.loc[sp.clade == c, "mean_div_weighted"].values for c in clades_sorted]
    bp = ax.boxplot(data, labels=clades_sorted, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", linewidth=1.5))
    for i, c in enumerate(clades_sorted):
        col = CLADE_COLOURS.get(c, "#888888")
        bp["boxes"][i].set_facecolor(col)
        bp["boxes"][i].set_alpha(0.5)
        vals = sp.loc[sp.clade == c, "mean_div_weighted"].values
        jitter = np.random.normal(0, 0.06, len(vals))
        ax.scatter(np.full_like(vals, i+1) + jitter, vals, c=col, s=18, alpha=0.7, edgecolors="none", zorder=3)
    ax.set_ylabel("Mean Kimura Divergence (weighted)")
    ax.set_xlabel("Clade")
    ax.set_title("Mean TE Divergence by Clade")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "kimura_mean_div_by_clade.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 5: Clade heatmap ──────────────────────────────────────────────────
def fig_clade_heatmap(sp):
    metrics = ["recency_index", "burst_score", "mean_div_weighted"]
    labels  = ["Recency Index", "Burst Score", "Mean Divergence"]
    agg = sp.groupby("clade")[metrics].mean()
    agg = agg.sort_values("recency_index", ascending=False)

    fig, ax = plt.subplots(figsize=(7, max(4, len(agg)*0.4)))
    data = agg.values
    # Normalise each column to [0,1] for colour mapping
    norm_data = (data - data.min(axis=0)) / (data.max(axis=0) - data.min(axis=0) + 1e-12)
    im = ax.imshow(norm_data, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yticks(range(len(agg)))
    ax.set_yticklabels(agg.index)
    # Annotate with actual values
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, f"{data[i,j]:.3f}", ha="center", va="center", fontsize=7,
                    color="white" if norm_data[i, j] > 0.6 else "black")
    ax.set_title("Mean Kimura Metrics per Clade")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "kimura_clade_heatmap.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 6: Representative Kimura landscapes ──────────────────────────────
NON_TE_SF = frozenset({
    "Simple_repeat", "Low_complexity", "Satellite",
    "rRNA", "tRNA", "snRNA", "scRNA", "sno", "srpRNA",
    "ARTEFACT", "acro", "centr", "subtelo", "DNA_virus", "RNA", "Y",
})

CLASS_COLOURS = {
    "DNA": "#e6194b",
    "LINE": "#3cb44b",
    "LTR": "#4363d8",
    "SINE": "#f58231",
    "Rolling Circle": "#911eb4",
    "Penelope": "#42d4f4",
    "Unclassified": "#aaaaaa",
    "Other": "#cccccc",
}

def load_divergence_raw(sp_name):
    """Load raw divergence_summary_table and return per-family rows."""
    pattern = os.path.join(EARLGREY, f"{sp_name}_EarlGrey",
                           f"{sp_name}_summaryFiles",
                           f"{sp_name}_divergence_summary_table.tsv")
    if not os.path.isfile(pattern):
        return None
    df = pd.read_csv(pattern, sep="\t")
    df = df.dropna(subset=["mean_div", "total_bp"])
    df["mean_div"] = pd.to_numeric(df["mean_div"], errors="coerce")
    df["total_bp"] = pd.to_numeric(df["total_bp"], errors="coerce")
    df = df.dropna(subset=["mean_div", "total_bp"])
    non_te = (df["subclass"] == "Other") & (df["superfamily"].isin(NON_TE_SF))
    df = df[~non_te].copy()
    return df

def fig_kimura_landscapes(sp_df):
    """Pick one representative species per clade (median TE content) and plot stacked landscapes."""
    # Pick representative species per clade
    sp_df = sp_df.copy()
    sp_df["te_pct_rank"] = sp_df.groupby("clade")["total_te_bp"].rank(method="first")
    sp_df["clade_n"] = sp_df.groupby("clade")["species"].transform("count")
    sp_df["median_rank"] = (sp_df["clade_n"] + 1) / 2
    sp_df["rank_dist"] = (sp_df["te_pct_rank"] - sp_df["median_rank"]).abs()
    reps = sp_df.loc[sp_df.groupby("clade")["rank_dist"].idxmin()]
    reps = reps.sort_values("clade")

    # Limit to 12 representative species for readability
    if len(reps) > 12:
        reps = reps.head(12)

    n = len(reps)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3.5*nrows), squeeze=False)

    bins = np.arange(0, 0.55, 0.01)  # 0 to 0.50 in steps of 0.01
    bin_centres = (bins[:-1] + bins[1:]) / 2

    classes_order = ["LTR", "DNA", "LINE", "SINE", "Rolling Circle", "Penelope", "Unclassified"]

    for idx, (_, row) in enumerate(reps.iterrows()):
        ax = axes[idx // ncols][idx % ncols]
        sp_name = row["species"]
        raw = load_divergence_raw(sp_name)
        if raw is None:
            ax.set_title(f"{sp_name}\n(no data)", fontsize=8)
            ax.set_visible(False)
            continue

        # Build histogram per class
        stacks = {}
        for cls in classes_order:
            mask = raw["subclass"].str.contains(cls, case=False, na=False) if cls != "Unclassified" else (raw["subclass"] == "Unclassified")
            if cls == "Rolling Circle":
                mask = raw["subclass"].str.contains("Rolling", case=False, na=False)
            subset = raw[mask]
            hist = np.zeros(len(bin_centres))
            for _, r in subset.iterrows():
                bi = np.searchsorted(bins, r["mean_div"]) - 1
                if 0 <= bi < len(hist):
                    hist[bi] += r["total_bp"]
            stacks[cls] = hist

        # Stacked area plot
        bottom = np.zeros(len(bin_centres))
        for cls in classes_order:
            vals = stacks.get(cls, np.zeros(len(bin_centres)))
            if vals.sum() > 0:
                ax.fill_between(bin_centres, bottom, bottom + vals,
                               label=cls, color=CLASS_COLOURS.get(cls, "#999"),
                               alpha=0.8, linewidth=0.3, edgecolor="white")
                bottom += vals

        ax.set_title(f"{sp_name} ({row['clade']})", fontsize=8)
        ax.set_xlim(0, 0.50)
        ax.set_xlabel("Kimura divergence", fontsize=7)
        ax.set_ylabel("bp", fontsize=7)
        ax.tick_params(labelsize=6)

    # Hide empty subplots
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Shared legend
    handles = [mpatches.Patch(color=CLASS_COLOURS[c], label=c) for c in classes_order]
    fig.legend(handles=handles, loc="lower center", ncol=len(classes_order), fontsize=7,
               framealpha=0.5, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Representative Kimura Landscapes (one per clade)", fontsize=11, y=1.01)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "kimura_landscape_stacked.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


# ── Figure 7: RI vs Shannon H' ──────────────────────────────────────────────
def fig_ri_vs_shannon(sp, sh):
    merged = sp.merge(sh[["species", "H_prime", "J_evenness"]], on="species", how="inner")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: RI vs H'
    ax = axes[0]
    for c in sorted(merged.clade.unique()):
        sub = merged[merged.clade == c]
        ax.scatter(sub.recency_index, sub.H_prime, label=c,
                   c=CLADE_COLOURS.get(c, "#888"), s=25, alpha=0.7, edgecolors="none")
    ax.set_xlabel("Recency Index (RI)")
    ax.set_ylabel("Shannon Diversity (H')")
    ax.set_title("A) RI vs Shannon H'")
    # Add correlation
    from scipy.stats import spearmanr
    r, p = spearmanr(merged.recency_index, merged.H_prime)
    ax.text(0.02, 0.98, f"Spearman r={r:.3f}, p={p:.1e}", transform=ax.transAxes,
            va="top", fontsize=8, bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    # Panel B: RI vs J
    ax = axes[1]
    for c in sorted(merged.clade.unique()):
        sub = merged[merged.clade == c]
        ax.scatter(sub.recency_index, sub.J_evenness, label=c,
                   c=CLADE_COLOURS.get(c, "#888"), s=25, alpha=0.7, edgecolors="none")
    ax.set_xlabel("Recency Index (RI)")
    ax.set_ylabel("Pielou Evenness (J)")
    ax.set_title("B) RI vs Evenness J")
    r, p = spearmanr(merged.recency_index, merged.J_evenness)
    ax.text(0.02, 0.98, f"Spearman r={r:.3f}, p={p:.1e}", transform=ax.transAxes,
            va="top", fontsize=8, bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    ax.legend(fontsize=5, ncol=3, loc="lower right", framealpha=0.5)

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "kimura_ri_vs_shannon.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# ── main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("[C02] Generating Kimura figures …")
    sp, sf, sh = load_data()
    print(f"  Loaded {len(sp)} species, {len(sf)} superfamily rows, {len(sh)} Shannon rows")

    fig_ri_by_clade(sp)
    fig_ri_vs_burst(sp)
    fig_ri_vs_meandiv(sp)
    fig_meandiv_by_clade(sp)
    fig_clade_heatmap(sp)
    fig_kimura_landscapes(sp)
    fig_ri_vs_shannon(sp, sh)

    print("[C02] Done — 7 figures generated.")
