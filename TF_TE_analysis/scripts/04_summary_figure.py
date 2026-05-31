#!/usr/bin/env python3
"""
04_summary_figure.py
────────────────────
Integrate the three analyses (peak overlap, BigWig signal, FIMO enrichment)
into a single multi-panel summary figure for publication / report.

Panel layout:
  A  Peak-TE overlap heatmap (fraction of peaks on each TE class)
  B  BigWig signal heatmap (z-scored log2 mean signal, TF × TE class)
  C  FIMO enrichment heatmap (median log2OR per motif family × TE class)
  D  Scatter: correlation between peak-overlap fraction and BigWig signal (per TF × TE class)

Outputs
  <outdir>/TF_TE_summary.pdf

Usage
  python 04_summary_figure.py

Requirements
  pip install pandas scipy matplotlib
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# ── CONFIG ──────────────────────────────────────────────────────────────────
RES_ROOT = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis")
OUTDIR   = RES_ROOT / "figures"

# Input files (produced by scripts 01, 02, 03)
OVERLAP_MATRIX = RES_ROOT / "01_peak_overlap" / "peak_TE_overlap_matrix.tsv"
OVERLAP_STATS  = RES_ROOT / "01_peak_overlap" / "peak_TE_overlap_stats.tsv"
SIGNAL_CLASS   = RES_ROOT / "02_bigwig_signal" / "bw_TE_signal_class_mean.tsv"
FIMO_ENRICH    = RES_ROOT / "03_fimo_enrichment" / "fimo_enrichment.tsv"

SKIP_TE_CLASSES = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA",
                   "snRNA", "Unknown", "Unspecified", "Other"}
# TE prefixes to skip (tRNA-derived SINEs and RNA genes)
SKIP_TE_PREFIXES = ("tRNA", "SINE/tRNA", "Retroposon/sno", "RNA", "rRNA",
                    "snRNA", "srpRNA")

TF_FAMILY_COLORS = {
    "HMGb":        "#E07B39",
    "HOX":         "#1F77B4",
    "HOX_KN":      "#AEC7E8",
    "HOX_Myb":     "#9AC7E8",
    "SANT_Myb":    "#2CA02C",
    "bHLH_Myc":    "#D62728",
    "bZIP":        "#9467BD",
    "ForkHead":    "#8C564B",
    "GATA_Znf":    "#E377C2",
    "HSF":         "#7F7F7F",
    "MADS":        "#BCBD22",
    "Znf_C2H2":    "#17BECF",
    "CP2":         "#98DF8A",
    "CopperFist":  "#FFBB78",
}


def _skip_te(cls_name):
    """Return True if this TE class should be excluded from analysis."""
    if cls_name in SKIP_TE_CLASSES:
        return True
    for pref in SKIP_TE_PREFIXES:
        if cls_name.startswith(pref):
            return True
    return False


def _cluster_order(mat, axis=0, method="ward", metric="euclidean"):
    """Return hierarchical-clustering leaf order for rows (axis=0) or cols (axis=1)."""
    data = mat if axis == 0 else mat.T
    if data.shape[0] < 2:
        return list(range(data.shape[0]))
    try:
        dist = pdist(data, metric=metric)
        Z = linkage(dist, method=method)
        return list(leaves_list(Z))
    except Exception:
        return list(range(data.shape[0]))

# ────────────────────────────────────────────────────────────────────────────


def load_data():
    """Load all analysis outputs; return None for missing datasets."""
    data = {}
    for key, path in [("overlap_mat", OVERLAP_MATRIX), ("overlap_stats", OVERLAP_STATS),
                      ("signal", SIGNAL_CLASS), ("fimo", FIMO_ENRICH)]:
        if path.exists():
            data[key] = pd.read_csv(path, sep="\t", index_col=0 if "matrix" in str(path) or "class" in str(path) else None)
        else:
            print(f"  WARNING: {path} not found — panel will be skipped")
            data[key] = None
    return data


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    print("Loading data …")
    d = load_data()

    # ── Figure layout ─────────────────────────────────────────────────────
    fig = plt.figure(figsize=(22, 18))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.45)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_labels = {"A": ax_a, "B": ax_b, "C": ax_c, "D": ax_d}
    for label, ax in panel_labels.items():
        ax.text(-0.08, 1.05, label, transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top")

    # ── Panel A: Peak-overlap fraction heatmap ───────────────────────────
    _panel_overlap_heatmap(ax_a, d["overlap_mat"], d["overlap_stats"])

    # ── Panel B: BigWig signal heatmap ───────────────────────────────────
    _panel_signal_heatmap(ax_b, d["signal"], d["overlap_stats"])

    # ── Panel C: FIMO enrichment heatmap ─────────────────────────────────
    _panel_fimo_heatmap(ax_c, d["fimo"])

    # ── Panel D: Peak-overlap vs signal correlation scatter ──────────────
    _panel_scatter(ax_d, d["overlap_mat"], d["overlap_stats"], d["signal"])

    fig.suptitle(
        "TF binding to transposable elements in Rhizopus microsporus ATCC11559\n"
        "DAP-seq, bigwig signal, and FIMO motif enrichment",
        fontsize=13, y=0.99
    )
    out = OUTDIR / "TF_TE_summary.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def _panel_overlap_heatmap(ax, mat_df, stats_df):
    """Panel A: fraction of peaks overlapping each TE class (clustered)."""
    ax.set_title("A  Peak overlap with TE classes (hierarchical clustering)", fontsize=10, loc="left")
    if mat_df is None or stats_df is None:
        ax.text(0.5, 0.5, "Data not available", ha="center", va="center",
                transform=ax.transAxes)
        return

    # Normalise by total peaks per TF
    totals = stats_df.set_index("tf_label")["n_peaks"]
    # Filter TE columns: remove tRNA, snRNA, simple/low complexity
    te_cols = [c for c in mat_df.columns
               if not _skip_te(c) and mat_df[c].max() > 0]
    frac = mat_df[te_cols].div(totals, axis=0).fillna(0)
    # Keep only TE classes with meaningful signal in at least 1 TF
    frac = frac.loc[:, frac.max() > 0.01]

    if frac.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes)
        return

    # Hierarchical clustering of rows (TFs by signal profile)
    row_idx = _cluster_order(frac.values, axis=0)
    col_idx = _cluster_order(frac.values, axis=1)
    frac = frac.iloc[row_idx, col_idx]

    # Build TF → family map
    fam_map = stats_df.set_index("tf_label")["tf_family"].to_dict()

    fig = ax.get_figure()
    # Reserve space for family color bar: narrow axes left of heatmap
    pos = ax.get_position()
    cbar_width = 0.012
    ax_cbar = fig.add_axes([pos.x0 - cbar_width - 0.005, pos.y0,
                            cbar_width, pos.height])

    families = [fam_map.get(tf, "Other") for tf in frac.index]
    fam_set   = sorted(set(families))
    fam_colors = [TF_FAMILY_COLORS.get(f, "#CCCCCC") for f in families]
    for i, col in enumerate(fam_colors):
        ax_cbar.add_patch(
            plt.Rectangle((0, i / len(frac)), 1, 1 / len(frac),
                           color=col, transform=ax_cbar.transAxes,
                           clip_on=False))
    ax_cbar.set_xlim(0, 1); ax_cbar.set_ylim(0, 1)
    ax_cbar.axis("off")
    # Legend patches
    legend_h = [mpatches.Patch(color=TF_FAMILY_COLORS.get(f, "#CCCCCC"), label=f)
                for f in fam_set if f in families]
    ax.legend(handles=legend_h, fontsize=4.5, loc="upper right",
              bbox_to_anchor=(1.0, 1.04), ncol=2, framealpha=0.7,
              title="TF family", title_fontsize=5)

    im = ax.imshow(frac.values, aspect="auto", cmap="YlOrRd",
                   vmin=0, vmax=frac.values.max())
    ax.set_xticks(range(len(frac.columns)))
    ax.set_xticklabels(frac.columns, rotation=40, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(frac.index)))
    tf_labels = [f"{tf}" for tf in frac.index]
    ax.set_yticklabels(tf_labels, fontsize=5.5)
    # Color y-tick labels by TF family
    for ytick, tf in zip(ax.get_yticklabels(), frac.index):
        fam = fam_map.get(tf, "Other")
        ytick.set_color(TF_FAMILY_COLORS.get(fam, "black"))
    ax.set_xlabel("TE class", fontsize=8)
    ax.set_ylabel("TF", fontsize=8)
    plt.colorbar(im, ax=ax, label="Fraction of peaks on TE", fraction=0.035, pad=0.02)


def _panel_signal_heatmap(ax, signal_df, stats_df=None):
    """Panel B: z-scored BigWig signal, TF×TE (clustered), TF labels colored by family."""
    ax.set_title("B  DAP-seq signal over TE classes — z-score (clustered)", fontsize=10, loc="left")
    if signal_df is None:
        ax.text(0.5, 0.5, "Data not available", ha="center", va="center",
                transform=ax.transAxes)
        return

    # Rows = TE classes, cols = TFs — filter out tRNA and other skipped classes
    sig = signal_df.loc[[i for i in signal_df.index if not _skip_te(str(i))]]
    if sig.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes)
        return

    from scipy.stats import zscore as _zscore
    z_arr = _zscore(sig.values.astype(float), axis=1, ddof=1, nan_policy='omit')
    z_arr = np.where(np.isnan(z_arr), 0.0, z_arr)   # constant rows → 0
    z = pd.DataFrame(z_arr, index=sig.index, columns=sig.columns)

    # Hierarchical clustering of both axes
    row_idx = _cluster_order(z.values, axis=0)
    col_idx = _cluster_order(z.values, axis=1)
    z = z.iloc[row_idx, col_idx]
    vmax = float(np.nanpercentile(np.abs(z.values), 95))

    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(z.columns)))
    ax.set_xticklabels(z.columns, rotation=40, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index, fontsize=5.5)
    ax.set_xlabel("TF", fontsize=8)
    ax.set_ylabel("TE class", fontsize=8)
    plt.colorbar(im, ax=ax, label="z-score (signal)", fraction=0.035, pad=0.02)

    # Color x-tick labels by TF family (if stats_df supplied)
    if stats_df is not None:
        fam_map = stats_df.set_index("tf_label")["tf_family"].to_dict()
        for xtick in ax.get_xticklabels():
            fam = fam_map.get(xtick.get_text(), "Other")
            xtick.set_color(TF_FAMILY_COLORS.get(fam, "black"))


def _panel_fimo_heatmap(ax, fimo_df):
    """Panel C: FIMO motif enrichment (log2OR) per TF family × TE class."""
    ax.set_title("C  Motif enrichment in TE sequences (FIMO)", fontsize=10, loc="left")
    if fimo_df is None:
        ax.text(0.5, 0.5, "Data not available", ha="center", va="center",
                transform=ax.transAxes)
        return

    sig = fimo_df.copy()
    sig.loc[sig["fisher_pval"] > 0.05, "log2_odds_ratio"] = 0
    pivot = sig.pivot_table(
        index="motif_family", columns="te_class",
        values="log2_odds_ratio", aggfunc="median"
    ).fillna(0)
    pivot = pivot.loc[:, [c for c in pivot.columns if not _skip_te(str(c))]]
    pivot = pivot.loc[(pivot.abs() > 0).any(axis=1)]

    if pivot.empty:
        ax.text(0.5, 0.5, "No significant enrichments\n(run script 03 first)",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        return

    vmax = max(1.0, np.nanpercentile(np.abs(pivot.values), 95))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=40, ha="right", fontsize=6)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_xlabel("TE class", fontsize=8)
    ax.set_ylabel("TF motif family", fontsize=8)
    plt.colorbar(im, ax=ax, label="Median log₂OR (p<0.05)", fraction=0.04, pad=0.02)


def _panel_scatter(ax, mat_df, stats_df, signal_df):
    """Panel D: per-TF × TE-class: peak-overlap fraction vs mean BigWig signal."""
    ax.set_title("D  Peak overlap fraction vs signal intensity", fontsize=10, loc="left")
    if mat_df is None or stats_df is None or signal_df is None:
        ax.text(0.5, 0.5, "Data not available\n(requires all 3 analyses)",
                ha="center", va="center", transform=ax.transAxes)
        return

    # Build common overlap-fraction and signal matrices with matching TF and TE class coords
    totals = stats_df.set_index("tf_label")["n_peaks"]
    te_cols = [c for c in mat_df.columns if not _skip_te(str(c))]
    frac = mat_df[te_cols].div(totals, axis=0).fillna(0)

    # Align signal_df: rows=TE classes, cols=TFs
    sig = signal_df.loc[
        signal_df.index.intersection(frac.columns), frac.index.intersection(signal_df.columns)
    ] if signal_df is not None else None

    if sig is None or sig.empty:
        ax.text(0.5, 0.5, "Insufficient shared TF/TE data", ha="center", va="center",
                transform=ax.transAxes)
        return

    # Flatten both to 1D arrays for scatter
    common_tes = frac.columns.intersection(signal_df.index)
    common_tfs = frac.index.intersection(signal_df.columns)
    x_vals, y_vals = [], []
    for tf in common_tfs:
        for te in common_tes:
            x_vals.append(frac.loc[tf, te])
            y_vals.append(signal_df.loc[te, tf])

    x_arr = np.array(x_vals)
    y_arr = np.array(y_vals)

    ax.scatter(x_arr, y_arr, alpha=0.4, s=10, color="steelblue", edgecolors="none")

    # Spearman correlation
    if len(x_arr) > 5:
        rho, pval = spearmanr(x_arr, y_arr)
        ax.text(0.05, 0.95,
                f"ρ = {rho:.3f}\np = {pval:.2e}\nn = {len(x_arr)}",
                transform=ax.transAxes, fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))

    ax.set_xlabel("Fraction of peaks on TE class", fontsize=8)
    ax.set_ylabel("log₂ mean BigWig signal (TE class)", fontsize=8)


if __name__ == "__main__":
    main()
