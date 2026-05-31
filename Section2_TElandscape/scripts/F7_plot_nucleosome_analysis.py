#!/usr/bin/env python3
"""
F7_plot_nucleosome_analysis.py
Generate all nucleosome-TE figures from F3/F4/F5/F6 outputs.

Figures produced:
  Fig A — TE nucleosome bias (delta_occupancy) by superfamily [F3 output]
  Fig B — NuclPred validation: Spearman r per scaffold violin [F4 output]
  Fig C — NDR/NOR TE enrichment heatmap by superfamily [F5 output]
  Fig D — Meta-TE occupancy profiles by superfamily (line plot) [F6 output]

Can run on single-species test data or multi-species combined TSVs.

Usage:
    python F7_plot_nucleosome_analysis.py \
        --f3 <f3_tsv>          (or glob pattern for multi-species) \
        --f4 <f4_scaffold_tsv> \
        --f5 <f5_tsv>          \
        --f6 <f6_profiles_tsv> \
        --output_dir <dir>     \
        [--species <SP>]       (for single-species title)

Each --fN argument can be a single TSV or a glob (quoted).
"""

import sys
import os
import glob
import argparse
import math
import warnings

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore")

# ── Colour palettes ────────────────────────────────────────────────────────────
TE_CLASS_COLORS = {
    "LTR":        "#e63946",
    "LTR/Gypsy":  "#e63946",
    "LTR/Copia":  "#f4a261",
    "LTR/ERV1":   "#e76f51",
    "LINE":       "#2a9d8f",
    "LINE/L1":    "#2a9d8f",
    "LINE/R1":    "#52b788",
    "LINE/R2":    "#40916c",
    "SINE":       "#457b9d",
    "DNA":        "#8338ec",
    "DNA/hAT":    "#8338ec",
    "DNA/hAT-Ac": "#8338ec",
    "DNA/Tc1":    "#6d28d9",
    "DNA/Mariner":"#7c3aed",
    "DNA/Maverick":"#4c1d95",
    "DNA/Dada":   "#9d4edd",
    "Helitron":   "#fb8500",
    "RC/Helitron":"#fb8500",
    "TIR":        "#f72585",
    "Unknown":    "#adb5bd",
    "Simple_repeat": "#dee2e6",
    "Low_complexity": "#dee2e6",
}

ZONE_COLORS = {"NDR": "#2196F3", "INT": "#9E9E9E", "NOR": "#F44336"}


def get_te_color(sf):
    if sf in TE_CLASS_COLORS:
        return TE_CLASS_COLORS[sf]
    broad = sf.split("/")[0]
    return TE_CLASS_COLORS.get(broad, "#adb5bd")


def load_tsvs(pattern):
    files = sorted(glob.glob(pattern)) if "*" in pattern else [pattern]
    frames = []
    for f in files:
        if os.path.isfile(f) and os.path.getsize(f) > 0:
            frames.append(pd.read_csv(f, sep="\t"))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ── Figure A: delta_occupancy by superfamily (F3) ─────────────────────────────
def plot_f3_bias(df, out_path, title=""):
    if df.empty:
        print("[F3] No data, skipping Fig A")
        return

    # Filter out repeat masker artefacts
    te_classes = ["Simple_repeat", "Low_complexity", "Satellite", "rRNA"]
    df = df[~df["superfamily"].isin(te_classes)].copy()
    if df.empty:
        print("[F3] No TE data after filtering, skipping Fig A")
        return

    # Sort by median delta
    order = (df.groupby("superfamily")["delta_occupancy"]
               .median().sort_values().index.tolist())
    df["superfamily"] = pd.Categorical(df["superfamily"], categories=order, ordered=True)
    df = df.sort_values("superfamily")

    n_sp = df["species"].nunique()
    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.7 + 2), 5))

    if n_sp == 1:
        # Single species: bar chart with value
        sp = df["species"].iloc[0]
        colors = [get_te_color(sf) for sf in order]
        medians = df.groupby("superfamily")["delta_occupancy"].median().reindex(order)
        counts  = df.groupby("superfamily")["n_insertions"].sum().reindex(order)
        bars = ax.bar(range(len(order)), medians.values, color=colors, edgecolor="white", width=0.7)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
        for i, (v, n) in enumerate(zip(medians.values, counts.values)):
            ax.text(i, v + (0.002 if v >= 0 else -0.008), f"n={int(n)}", ha="center",
                    va="bottom" if v >= 0 else "top", fontsize=6, color="#555")
    else:
        # Multi-species: boxplot
        groups = [df[df["superfamily"] == sf]["delta_occupancy"].values for sf in order]
        bp = ax.boxplot(groups, patch_artist=True, widths=0.6,
                        medianprops=dict(color="black", lw=1.5),
                        flierprops=dict(marker=".", markersize=3, alpha=0.4))
        for patch, sf in zip(bp["boxes"], order):
            patch.set_facecolor(get_te_color(sf))
            patch.set_alpha(0.8)
        ax.set_xticks(range(1, len(order) + 1))
        ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)

    ax.axhline(0, color="black", lw=0.8, linestyle="--", alpha=0.5)
    ax.set_ylabel("Δ Nucleosome occupancy\n(TE mean − genome background)", fontsize=9)
    ax.set_xlabel("TE superfamily", fontsize=9)
    ax.set_title(title or "TE nucleosome occupancy bias", fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    # Annotation: negative = TE avoids nucleosomes
    ax.text(0.02, 0.97, "← TE in NDR (depleted)", transform=ax.transAxes,
            fontsize=7, color="#2196F3", va="top")
    ax.text(0.02, 0.02, "→ TE in NOR (occupied)", transform=ax.transAxes,
            fontsize=7, color="#F44336", va="bottom")

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"[F3] Saved: {out_path}")


# ── Figure B: NuclPred validation Spearman r (F4) ─────────────────────────────
def plot_f4_validation(df, out_path, title=""):
    if df.empty:
        print("[F4] No data, skipping Fig B")
        return

    n_sp = df["species"].nunique()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Left: violin/strip of per-scaffold r per species
    ax = axes[0]
    species_list = sorted(df["species"].unique())
    palette = plt.cm.tab10(np.linspace(0, 0.9, len(species_list)))

    for i, sp in enumerate(species_list):
        sub = df[df["species"] == sp]["spearman_r"].dropna().values
        if len(sub) == 0:
            continue
        parts = ax.violinplot([sub], positions=[i], widths=0.7, showmedians=True,
                               showextrema=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(palette[i])
            pc.set_alpha(0.7)
        # strip jitter
        jitter = np.random.uniform(-0.15, 0.15, len(sub))
        ax.scatter(i + jitter, sub, color=palette[i], s=12, alpha=0.6, zorder=3)

    ax.axhline(0.4, color="red", lw=1, linestyle="--", alpha=0.7, label="r=0.4 threshold")
    ax.set_xticks(range(len(species_list)))
    ax.set_xticklabels(species_list, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Spearman r (per scaffold)", fontsize=9)
    ax.set_title("NuclPred vs Nuc-seq validation", fontsize=10)
    ax.legend(fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(-0.1, 1.0)

    # Right: summary bar with mean r ± SD per species
    ax2 = axes[1]
    summary = df.groupby("species")["spearman_r"].agg(["mean", "std", "count"]).reset_index()
    summary = summary.sort_values("mean", ascending=True)
    bars = ax2.barh(range(len(summary)), summary["mean"].values,
                    xerr=summary["std"].values, color=palette[:len(summary)],
                    capsize=3, error_kw=dict(elinewidth=1), alpha=0.8)
    ax2.set_yticks(range(len(summary)))
    ax2.set_yticklabels(summary["species"].values, fontsize=8)
    ax2.axvline(0.4, color="red", lw=1, linestyle="--", alpha=0.7)
    ax2.set_xlabel("Mean Spearman r", fontsize=9)
    ax2.set_title("Per-species mean ± SD", fontsize=10)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.set_xlim(0, 1)

    # Add n= annotations
    for i, row in enumerate(summary.itertuples()):
        pct = df[(df["species"] == row.species) & (df["spearman_r"] >= 0.4)].shape[0]
        tot = row.count
        ax2.text(row.mean + (row.std or 0) + 0.02, i,
                 f"n={pct}/{tot}", va="center", fontsize=7, color="#555")

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"[F4] Saved: {out_path}")


# ── Figure C: NDR/NOR enrichment heatmap (F5) ──────────────────────────────────
def plot_f5_ndr(df, out_path, title=""):
    if df.empty:
        print("[F5] No data, skipping Fig C")
        return

    # Filter non-TE rows
    mask = ~df["superfamily"].isin(["Simple_repeat", "Low_complexity", "Satellite",
                                     "rRNA", "Other"])
    df = df[mask].copy()
    df["log2_enrichment"] = pd.to_numeric(df["log2_enrichment"], errors="coerce")
    df = df.dropna(subset=["log2_enrichment"])
    if df.empty:
        return

    # Pivot: superfamily × zone, value = log2_enrichment
    pivot = df.pivot_table(index="superfamily", columns="zone",
                           values="log2_enrichment", aggfunc="mean")
    # Order: NDR, INT, NOR
    for col in ["NDR", "INT", "NOR"]:
        if col not in pivot.columns:
            pivot[col] = 0
    pivot = pivot[["NDR", "INT", "NOR"]]
    # Sort by NDR enrichment
    pivot = pivot.sort_values("NDR", ascending=False)

    fig, ax = plt.subplots(figsize=(5, max(4, len(pivot) * 0.35 + 1.5)))
    vmax = max(abs(pivot.values[~np.isnan(pivot.values)].max()),
               abs(pivot.values[~np.isnan(pivot.values)].min()), 1)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", norm=norm,
                   interpolation="nearest")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["NDR", "Intermediate", "NOR"], fontsize=9)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title(title or "TE enrichment in chromatin zones", fontsize=11, fontweight="bold")

    # Text annotations
    for i in range(len(pivot)):
        for j in range(3):
            v = pivot.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                        fontsize=6.5, color="white" if abs(v) > vmax * 0.5 else "black")

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("log₂ enrichment", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # Color y-labels by TE class
    for label in ax.get_yticklabels():
        label.set_color(get_te_color(label.get_text()))

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"[F5] Saved: {out_path}")


# ── Figure D: Meta-TE profiles (F6) ────────────────────────────────────────────
def plot_f6_meta(df, out_path, title="", max_superfamilies=8, smooth_sigma=0):
    if df.empty:
        print("[F6] No data, skipping Fig D")
        return

    # Filter non-TE
    mask = ~df["superfamily"].isin(["Simple_repeat", "Low_complexity", "Satellite", "rRNA"])
    df = df[mask].copy()
    df["mean_occupancy"] = pd.to_numeric(df["mean_occupancy"], errors="coerce")
    df["se_occupancy"]   = pd.to_numeric(df["se_occupancy"],   errors="coerce")
    df = df.dropna(subset=["mean_occupancy"])

    # Pick top superfamilies by number of insertions (most represented)
    top_sf = (df.groupby("superfamily")["n_insertions"].max()
                .sort_values(ascending=False)
                .head(max_superfamilies).index.tolist())
    df = df[df["superfamily"].isin(top_sf)]
    if df.empty:
        return

    n_sf = len(top_sf)
    ncols = min(3, n_sf)
    nrows = math.ceil(n_sf / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 4.5, nrows * 3),
                             sharey=False, sharex=True)
    if n_sf == 1:
        axes = [[axes]]
    elif nrows == 1:
        axes = [axes]

    for idx, sf in enumerate(top_sf):
        row_i, col_i = divmod(idx, ncols)
        ax = axes[row_i][col_i]
        sub = df[df["superfamily"] == sf].sort_values("relative_pos")
        x   = sub["relative_pos"].values
        y   = sub["mean_occupancy"].values
        se  = sub["se_occupancy"].fillna(0).values
        n   = sub["n_insertions"].max()
        color = get_te_color(sf)

        if smooth_sigma > 0:
            y_plot  = gaussian_filter1d(y,  sigma=smooth_sigma)
            se_plot = gaussian_filter1d(se, sigma=smooth_sigma)
        else:
            y_plot, se_plot = y, se

        ax.fill_between(x, y_plot - se_plot, y_plot + se_plot, alpha=0.25, color=color)
        ax.plot(x, y_plot, color=color, lw=1.5)
        ax.axvline(0, color="black", lw=0.8, linestyle="--", alpha=0.5, label="TE centre")
        ax.axhline(y_plot.mean(), color="grey", lw=0.6, linestyle=":", alpha=0.6)
        ax.set_title(f"{sf}\n(n={n})", fontsize=8, color=color, fontweight="bold")
        ax.set_xlabel("Position relative to TE centre (bp)", fontsize=7)
        ax.set_ylabel("Nucleosome occupancy", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.spines[["top", "right"]].set_visible(False)

    # Hide empty subplots
    for idx in range(n_sf, nrows * ncols):
        row_i, col_i = divmod(idx, ncols)
        axes[row_i][col_i].set_visible(False)

    sm_label = f" (Gaussian σ={smooth_sigma})" if smooth_sigma > 0 else ""
    fig.suptitle((title or "Meta-TE nucleosome profiles") + sm_label,
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"[F6] Saved: {out_path}")


# ── Figure E: combined 4-panel summary ─────────────────────────────────────────
def plot_combined(f3_df, f4_df, f5_df, f6_df, out_path, species="", smooth_sigma=2):
    """Compact 2×2 summary figure."""
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

    # ── Panel A: delta_occupancy bar ───────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])
    if not f3_df.empty:
        mask = ~f3_df["superfamily"].isin(["Simple_repeat", "Low_complexity", "Satellite"])
        sub  = f3_df[mask].copy()
        if not sub.empty:
            order = sub.groupby("superfamily")["delta_occupancy"].median().sort_values().index
            medians = sub.groupby("superfamily")["delta_occupancy"].median().reindex(order)
            colors  = [get_te_color(sf) for sf in order]
            ax_a.bar(range(len(order)), medians.values, color=colors, edgecolor="white")
            ax_a.set_xticks(range(len(order)))
            ax_a.set_xticklabels(order, rotation=40, ha="right", fontsize=7)
            ax_a.axhline(0, color="black", lw=0.8, linestyle="--", alpha=0.5)
            ax_a.set_ylabel("Δ Occupancy", fontsize=8)
    ax_a.set_title("A. TE nucleosome bias", fontsize=9, fontweight="bold")
    ax_a.spines[["top", "right"]].set_visible(False)

    # ── Panel B: Spearman r scatter ────────────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])
    if not f4_df.empty:
        species_list = sorted(f4_df["species"].unique())
        palette = plt.cm.tab10(np.linspace(0, 0.9, max(len(species_list), 1)))
        for i, sp in enumerate(species_list):
            sub = f4_df[f4_df["species"] == sp]["spearman_r"].dropna().values
            if len(sub) == 0:
                continue
            jitter = np.random.uniform(-0.2, 0.2, len(sub))
            ax_b.scatter(i + jitter, sub, color=palette[i], s=15, alpha=0.6)
            ax_b.scatter([i], [sub.mean()], color=palette[i], s=60, marker="D",
                         edgecolors="black", linewidths=0.5, zorder=5)
        ax_b.axhline(0.4, color="red", lw=1, linestyle="--", alpha=0.7)
        ax_b.set_xticks(range(len(species_list)))
        ax_b.set_xticklabels(species_list, rotation=30, ha="right", fontsize=7)
        ax_b.set_ylim(-0.1, 1.0)
        ax_b.set_ylabel("Spearman r", fontsize=8)
    ax_b.set_title("B. NuclPred validation (vs Nuc-seq)", fontsize=9, fontweight="bold")
    ax_b.spines[["top", "right"]].set_visible(False)

    # ── Panel C: NDR/NOR enrichment ────────────────────────────────────────────
    ax_c = fig.add_subplot(gs[1, 0])
    if not f5_df.empty:
        mask  = ~f5_df["superfamily"].isin(["Simple_repeat", "Low_complexity", "Satellite"])
        sub   = f5_df[mask].copy()
        sub["log2_enrichment"] = pd.to_numeric(sub["log2_enrichment"], errors="coerce")
        sub   = sub.dropna(subset=["log2_enrichment"])
        if not sub.empty:
            pivot = sub.pivot_table(index="superfamily", columns="zone",
                                    values="log2_enrichment", aggfunc="mean")
            for col in ["NDR", "INT", "NOR"]:
                if col not in pivot.columns:
                    pivot[col] = 0
            pivot = pivot[["NDR", "INT", "NOR"]].sort_values("NDR", ascending=False)
            vmax  = max(abs(pivot.values[~np.isnan(pivot.values)]).max(), 1)
            norm  = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
            ax_c.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                        norm=norm, interpolation="nearest")
            ax_c.set_xticks([0, 1, 2])
            ax_c.set_xticklabels(["NDR", "INT", "NOR"], fontsize=7)
            ax_c.set_yticks(range(len(pivot)))
            ax_c.set_yticklabels(pivot.index, fontsize=7)
    ax_c.set_title("C. TE enrichment NDR/NOR", fontsize=9, fontweight="bold")

    # ── Panel D: meta-TE profile (top 3 superfamilies) ────────────────────────
    ax_d = fig.add_subplot(gs[1, 1])
    if not f6_df.empty:
        mask  = ~f6_df["superfamily"].isin(["Simple_repeat", "Low_complexity"])
        sub   = f6_df[mask].copy()
        sub["mean_occupancy"] = pd.to_numeric(sub["mean_occupancy"], errors="coerce")
        sub   = sub.dropna(subset=["mean_occupancy"])
        if not sub.empty:
            top_sf = (sub.groupby("superfamily")["n_insertions"]
                        .max().sort_values(ascending=False).head(4).index)
            for sf in top_sf:
                s = sub[sub["superfamily"] == sf].sort_values("relative_pos")
                y = s["mean_occupancy"].values
                if smooth_sigma > 0:
                    y = gaussian_filter1d(y, sigma=smooth_sigma)
                ax_d.plot(s["relative_pos"], y,
                          label=sf, color=get_te_color(sf), lw=1.5)
            ax_d.axvline(0, color="black", lw=0.8, linestyle="--", alpha=0.5)
            ax_d.legend(fontsize=6, loc="upper right")
            ax_d.set_xlabel("Position (bp)", fontsize=7)
            ax_d.set_ylabel("Occupancy", fontsize=7)
    ax_d.set_title("D. Meta-TE profiles (top superfamilies)", fontsize=9, fontweight="bold")
    ax_d.spines[["top", "right"]].set_visible(False)

    sp_title = f" — {species}" if species else ""
    fig.suptitle(f"Nucleosome–TE landscape{sp_title}", fontsize=13, fontweight="bold")

    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"[COMBINED] Saved: {out_path}")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--f3",         default=None, help="F3 TSV (or glob)")
    parser.add_argument("--f4",         default=None, help="F4 scaffold correlations TSV")
    parser.add_argument("--f5",         default=None, help="F5 NDR/NOR enrichment TSV")
    parser.add_argument("--f6",         default=None, help="F6 meta-TE profiles TSV")
    parser.add_argument("--output_dir",    required=True)
    parser.add_argument("--species",        default="", help="Species label for titles")
    parser.add_argument("--smooth_sigma",   type=float, default=2,
                        help="Gaussian sigma for meta-TE profile smoothing (0=off, default=2)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    f3_df = load_tsvs(args.f3) if args.f3 else pd.DataFrame()
    f4_df = load_tsvs(args.f4) if args.f4 else pd.DataFrame()
    f5_df = load_tsvs(args.f5) if args.f5 else pd.DataFrame()
    f6_df = load_tsvs(args.f6) if args.f6 else pd.DataFrame()

    sp    = args.species or (f3_df["species"].iloc[0] if not f3_df.empty else "")

    plot_f3_bias(f3_df, os.path.join(args.output_dir, f"{sp}_figA_te_bias.pdf"),
                 title=f"TE nucleosome bias — {sp}")
    plot_f4_validation(f4_df, os.path.join(args.output_dir, f"{sp}_figB_validation.pdf"),
                       title=f"NuclPred validation — {sp}")
    plot_f5_ndr(f5_df, os.path.join(args.output_dir, f"{sp}_figC_ndr_enrichment.pdf"),
                title=f"TE enrichment in NDR/NOR — {sp}")
    plot_f6_meta(f6_df, os.path.join(args.output_dir, f"{sp}_figD_meta_te.pdf"),
                 title=f"Meta-TE profiles — {sp}", smooth_sigma=args.smooth_sigma)
    plot_combined(f3_df, f4_df, f5_df, f6_df,
                  os.path.join(args.output_dir, f"{sp}_fig_combined.pdf"),
                  species=sp, smooth_sigma=args.smooth_sigma)

    print(f"\n[DONE] All figures in: {args.output_dir}")


if __name__ == "__main__":
    main()
