#!/usr/bin/env python3
"""plot_kimura_compact.py — Compact 2-panel Kimura summary for main figure.

Panel A: RI boxplot by clade with individual points
Panel B: Heatmap of weighted-mean Kimura divergence per superfamily per clade,
         with optional row+column hierarchical clustering. Colour scale is
         REVERSED (dark = recent/low divergence, light = old/high divergence)
         and uses a strong sequential palette ("inferno_r") for high contrast.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

KIMURA_SP   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "C01_kimura_species.tsv")
KIMURA_SF   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "C01_kimura_superfamily.tsv")
OUTDIR      = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")

# Canonical clade palette (single source of truth)
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from clade_palette import CLADE_COLORS, CLADE_NAMES, CLADE_ORDER, get_clade  # noqa: E402


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df_sp = pd.read_csv(KIMURA_SP, sep="\t")
    df_sf = pd.read_csv(KIMURA_SF, sep="\t")
    # Drop control yeasts (VL clade) globally — user request 2026-04-21
    df_sp = df_sp[~df_sp["species"].str.startswith("VL_")].reset_index(drop=True)
    df_sf = df_sf[~df_sf["species"].str.startswith("VL_")].reset_index(drop=True)
    df_sp["clade"] = df_sp["species"].apply(get_clade)
    df_sf["clade"] = df_sf["species"].apply(get_clade)

    present = [c for c in CLADE_ORDER if c in df_sp["clade"].values]

    fig, (ax_box, ax_heat) = plt.subplots(1, 2, figsize=(11.5, 5.5),
                                           gridspec_kw={"width_ratios": [0.85, 1]})

    # ── Panel A: RI boxplot by clade (vertical, compact spacing) ─────────
    box_data = [df_sp[df_sp["clade"] == c]["recency_index"].dropna().values
                for c in present]
    box_colors = [CLADE_COLORS.get(c, "#ddd") for c in present]

    # Tight positions: gap of 0.72 between centres → boxes almost touching
    positions = np.arange(1, len(present) + 1) * 0.72
    bp = ax_box.boxplot(box_data, patch_artist=True, widths=0.55, vert=True,
                        positions=positions,
                        medianprops=dict(color="black", linewidth=1.5),
                        flierprops=dict(marker="o", markersize=2, alpha=0.3),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8))
    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Overlay jittered points (narrower jitter to match tighter boxes)
    for i, (clade, data) in enumerate(zip(present, box_data)):
        jitter = np.random.normal(0, 0.06, len(data))
        ax_box.scatter(np.full(len(data), positions[i]) + jitter, data,
                      c=CLADE_COLORS.get(clade, "#ddd"), s=10, alpha=0.5,
                      edgecolor="k", linewidth=0.2, zorder=3)

    short_labels = [CLADE_NAMES.get(c, c)[:12] for c in present]
    ax_box.set_xticks(positions)
    ax_box.set_xticklabels(short_labels, rotation=45, ha="right", fontsize=7)
    ax_box.set_ylabel("Recency Index (RI)\nfraction TE bp with div < 0.1", fontsize=9)
    ax_box.set_title("A. TE Recency Index by clade", fontsize=10, fontweight="bold",
                     loc="left")
    ax_box.spines["top"].set_visible(False)
    ax_box.spines["right"].set_visible(False)
    ax_box.set_ylim(-0.02, 1.02)

    # ── Panel B: Heatmap weighted mean Kimura per superfamily per clade ──
    # Keep only TE superfamilies (exclude Simple_repeat, rRNA, etc.)
    exclude = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA", "tRNA",
               "snRNA", "scRNA", "5S", "7SL", "Unknown"}
    df_te = df_sf[~df_sf["superfamily"].isin(exclude)].copy()

    # Top 15 superfamilies by total coverage across all species
    top_sf = (df_te.groupby("superfamily")["total_bp"].sum()
              .sort_values(ascending=False).head(15).index.tolist())
    df_te = df_te[df_te["superfamily"].isin(top_sf)]

    # Weighted mean divergence per clade × superfamily
    pivot = df_te.groupby(["clade", "superfamily"]).apply(
        lambda g: np.average(g["mean_div_weighted"], weights=g["total_bp"])
        if g["total_bp"].sum() > 0 else np.nan
    ).unstack(level="superfamily")
    pivot = pivot.reindex(index=present, columns=top_sf)

    # Friendlier superfamily labels (LTR/Gypsy etc.)
    SF_LABEL = {
        # LINE
        "L1":            "LINE/L1",
        "L2":            "LINE/L2",
        "R1":            "LINE/R1",
        "R2":            "LINE/R2",
        "RTE":           "LINE/RTE",
        "RTE-X":         "LINE/RTE-X",
        "RTE-BovB":      "LINE/RTE-BovB",
        "Jockey":        "LINE/Jockey",
        "I":             "LINE/I",
        "CR1":           "LINE/CR1",
        "Tx1":           "LINE/Tx1",
        "Penelope":      "PLE/Penelope",
        # LTR
        "Gypsy":         "LTR/Gypsy",
        "Copia":         "LTR/Copia",
        "Bel-Pao":       "LTR/Bel-Pao",
        "Pao":           "LTR/Pao",
        "ERV":           "LTR/ERV",
        "ERV1":          "LTR/ERV1",
        "ERVK":          "LTR/ERVK",
        "DIRS":          "LTR/DIRS",
        "Ngaro":         "LTR/Ngaro",
        # DNA
        "Tc1-Mariner":   "DNA/Tc1-Mariner",
        "hAT":           "DNA/hAT",
        "PIF-Harbinger": "DNA/PIF-Harbinger",
        "MULE":          "DNA/MULE",
        "PiggyBac":      "DNA/PiggyBac",
        "Crypton":       "DNA/Crypton",
        "Maverick":      "DNA/Maverick",
        "Polinton":      "DNA/Polinton",
        "TIR":           "DNA/TIR",
        "Sola":          "DNA/Sola",
        "Kolobok":       "DNA/Kolobok",
        "Academ":        "DNA/Academ",
        "Zator":         "DNA/Zator",
        # Rolling Circle
        "Helitron":      "RC/Helitron",
        # SINE
        "tRNA":          "SINE/tRNA",
        "5S":            "SINE/5S",
        "7SL":           "SINE/7SL",
        # Other
        "PLE":           "PLE",
    }
    sf_labels = [SF_LABEL.get(sf, f"?/{sf}") for sf in top_sf]

    # Fixed colour scale 0.15–0.30 (user request 2026-04-21):
    # values <0.15 saturate to darkest red, >0.30 saturate to white.
    vmin = 0.15
    vmax = 0.30
    # Red sequential colour scale: dark red = low divergence (RECENT bursts),
    # white = high divergence (old / decayed). Same colour family the user
    # originally requested, but with intensity reversed.
    im = ax_heat.imshow(pivot.values, aspect="auto", cmap="Reds_r",
                        vmin=vmin, vmax=vmax)
    ax_heat.set_xticks(range(len(top_sf)))
    ax_heat.set_xticklabels(sf_labels, rotation=45, ha="right", fontsize=7)
    ax_heat.set_yticks(range(len(present)))
    ax_heat.set_yticklabels([CLADE_NAMES.get(c, c)[:14] for c in present], fontsize=7)
    ax_heat.set_title("B. Weighted-mean Kimura divergence by superfamily\n"
                      "(dark red = recent / low divergence)",
                      fontsize=10, fontweight="bold", loc="left")

    cbar = fig.colorbar(im, ax=ax_heat, shrink=0.6, pad=0.02)
    cbar.set_label("Weighted mean Kimura\n(per-bp average over families)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.tight_layout()
    out = OUTDIR / "kimura_compact_summary.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "kimura_compact_summary.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[DONE] {out}")

    # ── Bonus: clustered version with explicit dendrograms ───────────────
    # Show the *raw* weighted-mean Kimura divergence (single-gradient
    # Reds_r colour map: dark red = recent/young TEs, white = old/divergent).
    # Cluster rows AND columns with **Ward linkage** (Euclidean).
    from scipy.cluster.hierarchy import dendrogram
    mat_clust = pivot.copy()
    mat_clust = mat_clust.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if mat_clust.shape[0] >= 3 and mat_clust.shape[1] >= 3:
        # Impute remaining NaNs with column mean for clustering only
        mat_imp = mat_clust.copy()
        for c in mat_imp.columns:
            mat_imp[c] = mat_imp[c].fillna(mat_imp[c].mean())
        # Ward linkage requires Euclidean distance on the raw values
        row_link = linkage(mat_imp.values,   method="ward", metric="euclidean")
        col_link = linkage(mat_imp.values.T, method="ward", metric="euclidean")
        row_order = leaves_list(row_link)
        col_order = leaves_list(col_link)

        clust_mat = mat_clust.iloc[row_order, col_order]

        # Two clustered heatmap variants requested by user (2026-04-22):
        #   v1 (default): vmax=0.25 — saturate above 0.25 to white (high contrast for the dark range)
        #   v2 (_full):   vmax=1.00 — full Kimura range
        for _suffix, _vmin, _vmax in [("",      0.15, 0.25),
                                       ("_full", 0.15, 1.00)]:
            fig2 = plt.figure(figsize=(11, 6.8))
            gs = fig2.add_gridspec(
                2, 3,
                width_ratios=[0.18, 1, 0.05],
                height_ratios=[0.18, 1],
                wspace=0.04, hspace=0.04,
            )
            ax_top  = fig2.add_subplot(gs[0, 1])
            ax_left = fig2.add_subplot(gs[1, 0])
            ax2     = fig2.add_subplot(gs[1, 1])
            ax_cb   = fig2.add_subplot(gs[1, 2])

            dendrogram(col_link, ax=ax_top, no_labels=True,
                       color_threshold=0, above_threshold_color="#444")
            ax_top.set_xticks([]); ax_top.set_yticks([])
            for s in ("top", "right", "bottom", "left"):
                ax_top.spines[s].set_visible(False)

            dendrogram(row_link, ax=ax_left, orientation="left", no_labels=True,
                       color_threshold=0, above_threshold_color="#444")
            ax_left.set_xticks([]); ax_left.set_yticks([])
            for s in ("top", "right", "bottom", "left"):
                ax_left.spines[s].set_visible(False)

            plot_mat = clust_mat.iloc[::-1]
            im2 = ax2.imshow(plot_mat.values, aspect="auto", cmap="Reds_r",
                             vmin=_vmin, vmax=_vmax)
            ax2.set_xticks(range(plot_mat.shape[1]))
            ax2.set_xticklabels(
                [SF_LABEL.get(sf, f"?/{sf}") for sf in plot_mat.columns],
                rotation=45, ha="right", fontsize=9)
            ax2.set_yticks(range(plot_mat.shape[0]))
            ax2.set_yticklabels(
                [CLADE_NAMES.get(c, c) for c in plot_mat.index], fontsize=9)
            ax2.tick_params(axis="x", which="both", top=False)

            for i in range(plot_mat.shape[0]):
                for j in range(plot_mat.shape[1]):
                    v = plot_mat.iloc[i, j]
                    if np.isfinite(v):
                        txt_col = "white" if v < (_vmin + (_vmax - _vmin) * 0.4) else "black"
                        ax2.text(j, i, f"{v*100:.0f}", ha="center", va="center",
                                 fontsize=6.5, color=txt_col)

            cb2 = fig2.colorbar(im2, cax=ax_cb)
            cb2.set_label("Weighted mean Kimura divergence\n"
                          "(dark red = recent / low div., white = old / high div.)",
                          fontsize=8)
            cb2.ax.tick_params(labelsize=7)

            fig2.suptitle(f"Weighted-mean Kimura divergence \u2014 Ward clustering "
                          f"(scale {_vmin:.2f}\u2013{_vmax:.2f})\n"
                          "(rows: clades, cols: TE superfamilies; "
                          "Ward linkage, Euclidean distance; numbers = mean div \u00d7 100)",
                          fontsize=11, fontweight="bold", y=0.995)
            out2 = OUTDIR / f"kimura_compact_summary_clustered{_suffix}.pdf"
            fig2.savefig(out2, dpi=300, bbox_inches="tight")
            fig2.savefig(out2.with_suffix(".png"), dpi=150, bbox_inches="tight")
            plt.close(fig2)
            print(f"[DONE] {out2}")


if __name__ == "__main__":
    main()
