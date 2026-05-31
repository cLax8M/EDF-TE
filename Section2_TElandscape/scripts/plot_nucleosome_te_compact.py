#!/usr/bin/env python3
"""plot_nucleosome_te_compact.py — Compact nucleosome-TE interaction figure.

2-panel figure for main text:
  A. Delta occupancy (TE insertion vs background) by TE class, across clades
  B. NDR/NOR enrichment heatmap: TE enrichment in nucleosome-depleted vs
     occupied regions, per TE class per clade

Data: F3_bias and F5_ndr per-species TSVs (205 species).
"""
from pathlib import Path
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

F3_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
              "nucleosomes/f3_bias")
F5_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
              "nucleosomes/f5_ndr")
OUTDIR = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")

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
    "AS": "Ascomycota", "NU": "Nucleariida",
    "MO": "Mortierellomycota", "GL": "Glomeromycotina",
    "HO": "Holozoa", "CR": "Cryptomycota",
    "OL": "Olpidiomycota", "ZO": "Zoopagomycota",
    "EN": "Entomophthoromycota", "KX": "Kickxellomycota",
    "BS": "Basidiobolomycota",
}

CLADE_ORDER = ["HO", "CR", "MI", "AP", "OL", "BL", "CH", "NU", "BS",
               "EN", "ZO", "KX", "MO", "GL", "MU", "BA", "AS"]

# Broad TE classes to show
TE_BROAD = ["LTR", "DNA", "LINE", "SINE", "Unknown"]
TE_BROAD_COLORS = {"LTR": "#e41a1c", "DNA": "#377eb8", "LINE": "#4daf4a",
                   "SINE": "#f1c232", "Unknown": "#bdbdbd"}


def load_multi_tsv(directory, suffix):
    frames = []
    for f in sorted(directory.glob(f"*{suffix}")):
        if f.stat().st_size > 0:
            try:
                frames.append(pd.read_csv(f, sep="\t"))
            except Exception:
                pass
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def to_broad_class(sf):
    if "/" in sf:
        sf = sf.split("/")[0]
    if sf in TE_BROAD:
        return sf
    return None  # skip non-TE superfamilies


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── Load F3 (bias) data ───────────────────────────────────────────────
    df_f3 = load_multi_tsv(F3_DIR, "_f3_bias.tsv")
    print(f"[INFO] F3 bias: {len(df_f3)} rows, {df_f3['species'].nunique()} species")

    df_f3["clade"] = df_f3["species"].str[:2]
    df_f3["broad_class"] = df_f3["superfamily"].apply(to_broad_class)
    df_f3 = df_f3.dropna(subset=["broad_class"])
    # Filter to insertions with at least 5 copies
    df_f3 = df_f3[df_f3["n_insertions"] >= 5]

    # Per-species per-broad-class: weighted mean delta
    f3_agg = (df_f3.groupby(["species", "clade", "broad_class"])
              .apply(lambda g: pd.Series({
                  "delta": np.average(g["delta_occupancy"], weights=g["n_insertions"]),
                  "n": g["n_insertions"].sum()
              }))
              .reset_index())

    # ── Load F5 (NDR/NOR) data ────────────────────────────────────────────
    df_f5 = load_multi_tsv(F5_DIR, "_f5_ndr.tsv")
    print(f"[INFO] F5 NDR: {len(df_f5)} rows, {df_f5['species'].nunique()} species")

    df_f5["clade"] = df_f5["species"].str[:2]
    df_f5["broad_class"] = df_f5["superfamily"].apply(to_broad_class)
    df_f5 = df_f5.dropna(subset=["broad_class"])
    df_f5 = df_f5[df_f5["n_insertions"] >= 5]

    # ── Figure ────────────────────────────────────────────────────────────
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6),
                                      gridspec_kw={"width_ratios": [1.2, 1]})

    # ── Panel A: Delta occupancy boxplot by TE class ──────────────────────
    # Group by broad class, show distribution across species
    positions = []
    labels = []
    pos = 0
    group_centers = []
    for cls in TE_BROAD:
        sub = f3_agg[f3_agg["broad_class"] == cls]["delta"].dropna().values
        if len(sub) == 0:
            continue
        # Clip extreme outliers for visualization
        sub_clip = np.clip(sub, -0.2, 0.2)
        pos += 1
        bp = ax_a.boxplot([sub_clip], positions=[pos], widths=0.6,
                          patch_artist=True,
                          medianprops=dict(color="black", linewidth=1.5),
                          flierprops=dict(marker=".", markersize=2, alpha=0.3))
        bp["boxes"][0].set_facecolor(TE_BROAD_COLORS.get(cls, "#aaa"))
        bp["boxes"][0].set_alpha(0.7)
        # Jittered points colored by clade
        jitter = np.random.normal(0, 0.08, len(sub_clip))
        clades_sub = f3_agg[f3_agg["broad_class"] == cls]["clade"].values
        colors = [CLADE_COLORS.get(c, "#aaa") for c in clades_sub]
        ax_a.scatter(np.full(len(sub_clip), pos) + jitter, sub_clip,
                    c=colors, s=8, alpha=0.5, edgecolor="k", linewidth=0.15, zorder=3)
        labels.append(cls)
        group_centers.append(pos)

    ax_a.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax_a.set_xticks(group_centers)
    ax_a.set_xticklabels(labels, fontsize=9)
    ax_a.set_ylabel("Δ nucleosome occupancy\n(TE insertion − background)", fontsize=9)
    ax_a.set_ylim(-0.22, 0.22)
    ax_a.set_title("A. TE nucleosome insertion bias", fontsize=10,
                   fontweight="bold", loc="left")
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    # ── Panel B: NDR/NOR enrichment heatmap ───────────────────────────────
    # For each clade × TE class, compute mean log2 enrichment in NDR vs NOR
    df_ndr = df_f5[df_f5["zone"] == "NDR"].copy()
    df_nor = df_f5[df_f5["zone"] == "NOR"].copy()

    # Mean log2_enrichment per clade × broad_class for NDR
    ndr_agg = (df_ndr.groupby(["clade", "broad_class"])
               .apply(lambda g: np.average(g["log2_enrichment"],
                                           weights=g["n_insertions"])
                      if g["n_insertions"].sum() > 0 else np.nan)
               .unstack(level="broad_class"))

    nor_agg = (df_nor.groupby(["clade", "broad_class"])
               .apply(lambda g: np.average(g["log2_enrichment"],
                                           weights=g["n_insertions"])
                      if g["n_insertions"].sum() > 0 else np.nan)
               .unstack(level="broad_class"))

    # Combine: NDR enrichment - NOR enrichment = preference for NDR
    diff = ndr_agg.sub(nor_agg)
    present_clades = [c for c in CLADE_ORDER if c in diff.index]
    present_classes = [c for c in TE_BROAD if c in diff.columns]
    diff = diff.reindex(index=present_clades, columns=present_classes)

    norm = TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.5)
    im = ax_b.imshow(diff.values, aspect="auto", cmap="RdBu_r", norm=norm)
    ax_b.set_xticks(range(len(present_classes)))
    ax_b.set_xticklabels(present_classes, fontsize=9, rotation=45, ha="right")
    ax_b.set_yticks(range(len(present_clades)))
    ax_b.set_yticklabels([CLADE_NAMES.get(c, c)[:14] for c in present_clades],
                         fontsize=7)
    ax_b.set_title("B. NDR vs NOR preference\n(log₂ NDR enrichment − log₂ NOR enrichment)",
                   fontsize=10, fontweight="bold", loc="left")

    # Annotate cells
    for i in range(len(present_clades)):
        for j in range(len(present_classes)):
            val = diff.values[i, j]
            if not np.isnan(val):
                color = "white" if abs(val) > 0.3 else "black"
                ax_b.text(j, i, f"{val:.2f}", ha="center", va="center",
                         fontsize=6, color=color)

    cbar = fig.colorbar(im, ax=ax_b, shrink=0.7, pad=0.02)
    cbar.set_label("NDR − NOR preference", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    fig.tight_layout()
    out = OUTDIR / "nucleosome_te_compact.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "nucleosome_te_compact.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[DONE] {out}")

    # ── Summary statistics ────────────────────────────────────────────────
    print(f"\n{'='*50}")
    print("Delta occupancy (TE vs background) by class:")
    for cls in TE_BROAD:
        sub = f3_agg[f3_agg["broad_class"] == cls]["delta"]
        if len(sub) > 0:
            print(f"  {cls:12s}  median={sub.median():+.4f}  "
                  f"mean={sub.mean():+.4f}  n_species={len(sub)}")

    print("\nNDR vs NOR preference by clade (positive = prefer NDR):")
    for c in present_clades:
        row = diff.loc[c]
        vals = row.dropna()
        if len(vals) > 0:
            print(f"  {CLADE_NAMES.get(c,c):20s}  mean={vals.mean():+.3f}")


if __name__ == "__main__":
    main()
