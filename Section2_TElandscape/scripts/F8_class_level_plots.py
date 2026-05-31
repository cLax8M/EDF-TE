#!/usr/bin/env python3
"""
F8_class_level_plots.py
Aggregate F3 (Δ-occupancy bias), F5 (NDR/NOR/INT enrichment) and F6
(meta-TE midpoint profile) to TE-class level (DNA, LTR, LINE, SINE,
Rolling Circle, PLE) and produce per-class summary figures.

Usage:
    python F8_class_level_plots.py \\
        --f3 <combined_f3.tsv> --f5 <combined_f5.tsv> \\
        --f6 <combined_f6.tsv> --output_dir <dir> --label <variant>
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clade_palette import TE_CLASS_COLORS  # noqa: E402

CLASSES = ["DNA", "LTR", "LINE", "SINE", "Rolling Circle", "PLE"]


def superfamily_to_class(sf: str) -> str:
    if not isinstance(sf, str):
        return "Other"
    head = sf.split("/", 1)[0]
    head = head.rstrip("?")
    mapping = {
        "DNA": "DNA",
        "LTR": "LTR",
        "LINE": "LINE",
        "SINE": "SINE",
        "RC": "Rolling Circle",
        "Helitron": "Rolling Circle",
        "PLE": "PLE",
        "Penelope": "PLE",
        "Retroposon": "LINE",
    }
    return mapping.get(head, "Other")


# ────────────────────────────────────────────────────────────────────────
# Figure A: Δ-occupancy boxplot per class across species
# ────────────────────────────────────────────────────────────────────────
def plot_class_bias(f3_df, out_path, label):
    df = f3_df.copy()
    df["te_class"] = df["superfamily"].map(superfamily_to_class)
    df = df[df["te_class"].isin(CLASSES)]
    df["delta_occupancy"] = pd.to_numeric(df["delta_occupancy"], errors="coerce")
    df = df.dropna(subset=["delta_occupancy"])
    if df.empty:
        print("[A] no class-level data, skip"); return
    # weighted mean per (species, class)
    df["w"] = pd.to_numeric(df["n_insertions"], errors="coerce").fillna(0)
    g = (df.groupby(["species", "te_class"])
           .apply(lambda x: np.average(x["delta_occupancy"],
                                       weights=x["w"]) if x["w"].sum() > 0
                  else x["delta_occupancy"].mean())
           .rename("delta").reset_index())

    fig, ax = plt.subplots(figsize=(8, 5))
    data, labels, colors = [], [], []
    for c in CLASSES:
        sub = g[g["te_class"] == c]["delta"].dropna().values
        if len(sub) == 0:
            continue
        data.append(sub); labels.append(f"{c}\n(n={len(sub)})"); colors.append(TE_CLASS_COLORS[c])
    bp = ax.boxplot(data, patch_artist=True, widths=0.65,
                    medianprops=dict(color="black", lw=1.5),
                    flierprops=dict(marker=".", markersize=3, alpha=0.4))
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col); patch.set_alpha(0.85)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=9)
    ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.set_ylabel("Δ Nucleosome occupancy\n(TE − genome background)", fontsize=10)
    ax.set_title(f"TE-class nucleosome bias — {label}\n"
                 "Each point = one species (insertion-weighted)", fontsize=11)
    ax.text(0.02, 0.97, "← TE in NDR (depleted)", transform=ax.transAxes,
            fontsize=8, color="#2196F3", va="top")
    ax.text(0.02, 0.02, "→ TE in NOR (occupied)", transform=ax.transAxes,
            fontsize=8, color="#F44336", va="bottom")
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"[A] saved {out_path}")


# ────────────────────────────────────────────────────────────────────────
# Figure B: NDR/NOR/INT log2-enrichment per class
# ────────────────────────────────────────────────────────────────────────
def plot_class_ndr(f5_df, out_path, label):
    df = f5_df.copy()
    df["te_class"] = df["superfamily"].map(superfamily_to_class)
    df = df[df["te_class"].isin(CLASSES)]
    df["log2_enrichment"] = pd.to_numeric(df["log2_enrichment"], errors="coerce")
    df["n_insertions"] = pd.to_numeric(df["n_insertions"], errors="coerce").fillna(0)
    df = df.dropna(subset=["log2_enrichment"])
    df = df[np.isfinite(df["log2_enrichment"])]
    if df.empty:
        print("[B] no data"); return

    # weighted mean log2 by zone × class across species
    g = (df.groupby(["zone", "te_class"])
           .apply(lambda x: np.average(x["log2_enrichment"], weights=x["n_insertions"])
                  if x["n_insertions"].sum() > 0 else x["log2_enrichment"].mean())
           .unstack("te_class")
           .reindex(["NDR", "INT", "NOR"])[CLASSES])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                             gridspec_kw={"width_ratios": [1.2, 1]})

    # --- left: heatmap mean log2 enrichment ---
    ax = axes[0]
    arr = g.values
    vmax = max(0.5, np.nanmax(np.abs(arr)))
    im = ax.imshow(arr, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(CLASSES)))
    ax.set_xticklabels(CLASSES, rotation=30, ha="right")
    ax.set_yticks(range(len(g.index)))
    ax.set_yticklabels(g.index)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                        fontsize=9,
                        color="white" if abs(v) > vmax * 0.55 else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Mean log2(observed/expected)\nweighted by n_insertions",
                   fontsize=9)
    ax.set_title(f"TE-class enrichment in NDR/INT/NOR — {label}", fontsize=11)

    # --- right: per-species boxplot for NDR only ---
    ax2 = axes[1]
    ndr = df[df["zone"] == "NDR"]
    data, labels, colors = [], [], []
    for c in CLASSES:
        sub = ndr[ndr["te_class"] == c]
        # weighted mean per species
        gs = (sub.groupby("species")
                 .apply(lambda x: np.average(x["log2_enrichment"], weights=x["n_insertions"])
                        if x["n_insertions"].sum() > 0 else np.nan)
                 .dropna().values)
        if len(gs) == 0: continue
        data.append(gs); labels.append(f"{c}\n(n={len(gs)})")
        colors.append(TE_CLASS_COLORS[c])
    bp = ax2.boxplot(data, patch_artist=True, widths=0.65,
                     medianprops=dict(color="black", lw=1.5),
                     flierprops=dict(marker=".", markersize=3, alpha=0.4))
    for p, col in zip(bp["boxes"], colors):
        p.set_facecolor(col); p.set_alpha(0.85)
    ax2.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax2.set_xticks(range(1, len(labels) + 1))
    ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("log2(observed / expected) in NDR", fontsize=10)
    ax2.set_title("Per-species NDR preference by TE class", fontsize=11)
    ax2.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"[B] saved {out_path}")


# ────────────────────────────────────────────────────────────────────────
# Figure C: meta-TE midpoint profile per class (collapsing subfamilies)
# ────────────────────────────────────────────────────────────────────────
def plot_class_meta(f6_df, out_path, label):
    df = f6_df.copy()
    df["te_class"] = df["superfamily"].map(superfamily_to_class)
    df = df[df["te_class"].isin(CLASSES)]
    df["mean_occupancy"] = pd.to_numeric(df["mean_occupancy"], errors="coerce")
    df["n_insertions"]   = pd.to_numeric(df["n_insertions"], errors="coerce").fillna(0)
    df["relative_pos"]   = pd.to_numeric(df["relative_pos"], errors="coerce")
    df = df.dropna(subset=["mean_occupancy", "relative_pos"])
    if df.empty:
        print("[C] no data"); return

    # within each species: normalise occupancy by species mean (to remove
    # global sequencing-depth differences)
    df["norm_occ"] = df.groupby("species")["mean_occupancy"].transform(
        lambda x: x / x.mean() if x.mean() else x)

    # collapse to class × position with insertion-weighted mean across (species, sf)
    g = (df.groupby(["te_class", "relative_pos"])
           .apply(lambda x: np.average(x["norm_occ"], weights=x["n_insertions"])
                  if x["n_insertions"].sum() > 0 else x["norm_occ"].mean())
           .rename("occ").reset_index())

    fig, ax = plt.subplots(figsize=(9, 5))
    for c in CLASSES:
        sub = g[g["te_class"] == c].sort_values("relative_pos")
        if sub.empty: continue
        # smooth (rolling 31bp)
        y = sub["occ"].rolling(31, center=True, min_periods=1).mean().values
        ax.plot(sub["relative_pos"], y, color=TE_CLASS_COLORS[c],
                lw=1.8, label=c)
    ax.axvline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axhline(1.0, color="grey", lw=0.6, ls=":", alpha=0.6)
    ax.set_xlabel("Position relative to TE midpoint (bp)", fontsize=10)
    ax.set_ylabel("Normalised nucleosome occupancy\n(species mean = 1)", fontsize=10)
    ax.set_title(f"Meta-TE midpoint profile by class — {label}\n"
                 "Insertion-weighted across species; 31 bp rolling mean",
                 fontsize=11)
    ax.legend(loc="best", fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"[C] saved {out_path}")


# ────────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--f3", required=True)
    p.add_argument("--f5", required=True)
    p.add_argument("--f6", required=True)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--label", required=True)
    a = p.parse_args()
    os.makedirs(a.output_dir, exist_ok=True)

    print(f"[load] F3 = {a.f3}")
    f3 = pd.read_csv(a.f3, sep="\t")
    print(f"[load] F5 = {a.f5}")
    f5 = pd.read_csv(a.f5, sep="\t")
    print(f"[load] F6 = {a.f6}")
    f6 = pd.read_csv(a.f6, sep="\t")
    print(f"  rows F3={len(f3)} F5={len(f5)} F6={len(f6)}")

    plot_class_bias(f3, os.path.join(a.output_dir, f"{a.label}_classA_bias.pdf"), a.label)
    plot_class_ndr(f5,  os.path.join(a.output_dir, f"{a.label}_classB_ndr.pdf"),  a.label)
    plot_class_meta(f6, os.path.join(a.output_dir, f"{a.label}_classC_metaTE.pdf"), a.label)
    print(f"[DONE] {a.output_dir}")


if __name__ == "__main__":
    main()
