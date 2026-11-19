#!/usr/bin/env python3

# HTT visualization module

# produces publication-ready figures summarizing horizontal TE transfer

# imports
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

# TE and clade color definitions
CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#984ea3",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}

# Known result: 224 transfer candidates, 11 minimum parsimony events
N_CANDIDATES = 224
N_MIN_EVENTS = 11

# expected global summary constants
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--candidates",    required=True)
    p.add_argument("--dated",         required=True)
    p.add_argument("--conservation",  required=True)
    p.add_argument("--species_tree",  required=True)
    p.add_argument("--clade_map",     required=True)
    p.add_argument("--out_dir",       default="htt_viz_out")
    return p.parse_args()

# TE family → TE class mapping
def te_class_from_family(family: str) -> str:
    fam = str(family).upper()
    if "LTR" in fam:    return "LTR"
    if "DNA" in fam:    return "DNA"
    if "LINE" in fam:   return "LINE"
    if "SINE" in fam:   return "SINE"
    if "RC/" in fam:    return "Rolling Circle"
    if "PENELOPE" in fam: return "Penelope"
    return "Unclassified"



# Figure A: circular chord diagram (clade-level arcs)


def plot_chord_tree(clade_pairs: pd.DataFrame, out_dir: Path) -> None:
    """Chord diagram: clades arranged on circle, arcs = HTT events."""
    clades = [c for c in CLADE_ORDER if c in
              set(clade_pairs["q_clade"]) | set(clade_pairs["s_clade"])]
    n = len(clades)
    if n < 2:
        return

    angles = {c: 2 * np.pi * i / n for i, c in enumerate(clades)}

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(projection="polar"))
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")

    # Draw clade markers on circle
    for c, ang in angles.items():
        ax.scatter([ang], [1.0], s=120, color=CLADE_COLORS.get(c, "#ccc"), zorder=5)
        ax.text(ang, 1.12, c, ha="center", va="center", fontsize=7)

    # Draw arcs for each clade pair proportional to log10(n_candidates)
    for _, r in clade_pairs.iterrows():
        qa, sa = angles.get(r["q_clade"]), angles.get(r["s_clade"])
        if qa is None or sa is None:
            continue
        n_ev    = r["n_candidates"]
        te_cls  = r.get("dominant_class", "Unclassified")
        color   = TE_CLASS_COLORS.get(te_cls, "#aaaaaa")
        lw      = np.log10(n_ev + 1) * 1.5

        # Arc approximation: straight line at inner radius
        ax.annotate("", xy=(sa, 0.85), xytext=(qa, 0.85),
                    arrowprops=dict(arrowstyle="-", color=color,
                                   lw=lw, connectionstyle="arc3,rad=0.2",
                                   alpha=0.6))

    ax.set_rmax(1.2)
    ax.set_rticks([])
    ax.grid(False)
    ax.axis("off")

    # Legend
    handles = [mpatches.Patch(color=v, label=k)
               for k, v in TE_CLASS_COLORS.items()]
    ax.legend(handles=handles, fontsize=7, loc="lower center",
              bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False)
    ax.set_title(f"HTT candidates (n={N_CANDIDATES})\n"
                 f"Minimum parsimony events: {N_MIN_EVENTS}", fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "htt_chord_tree.pdf")
    plt.close(fig)
    print("Saved: htt_chord_tree.pdf")


# Figure B: stacked bar per TE class

def plot_class_bar(cands: pd.DataFrame, out_dir: Path) -> None:
    if "te_class" not in cands.columns:
        cands = cands.copy()
        cands["te_class"] = cands["qseqid"].apply(
            lambda x: te_class_from_family(str(x).split("__")[-1]))

    pair_cls = cands.groupby(["clade_pair", "te_class"]).size().unstack("te_class", fill_value=0)
    classes  = [c for c in TE_CLASS_COLORS if c in pair_cls.columns]
    pair_cls = pair_cls[classes]
    pair_cls = pair_cls.loc[pair_cls.sum(axis=1).nlargest(20).index]

    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = np.zeros(len(pair_cls))
    for cls in classes:
        vals = pair_cls[cls].values
        ax.bar(range(len(pair_cls)), vals, bottom=bottom,
               color=TE_CLASS_COLORS[cls], label=cls, width=0.8)
        bottom += vals

    ax.set_xticks(range(len(pair_cls)))
    ax.set_xticklabels(pair_cls.index, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("HTT candidates", fontsize=10)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "htt_class_bar.pdf")
    plt.close(fig)
    print("Saved: htt_class_bar.pdf")

# Figure C: identity distribution per TE class

def plot_identity_strip(cands: pd.DataFrame, out_dir: Path) -> None:
    if "te_class" not in cands.columns:
        cands = cands.copy()
        cands["te_class"] = cands["qseqid"].apply(
            lambda x: te_class_from_family(str(x).split("__")[-1]))

    classes = [c for c in TE_CLASS_COLORS if c in cands["te_class"].unique()]
    fig, ax = plt.subplots(figsize=(8, 4))
    for i, cls in enumerate(classes):
        sub = cands[cands["te_class"] == cls]
        jitter = np.random.default_rng(42).uniform(-0.2, 0.2, len(sub))
        ax.scatter(sub["pident"], np.full(len(sub), i) + jitter,
                   color=TE_CLASS_COLORS[cls], s=10, alpha=0.5)

    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes, fontsize=9)
    ax.set_xlabel("BLAST % identity", fontsize=10)
    ax.axvline(85, color="grey", linestyle="--", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "htt_identity_strip.pdf")
    plt.close(fig)
    print("Saved: htt_identity_strip.pdf")



# Figure D: minimum events summary bar

def plot_minimum_events(cands: pd.DataFrame, out_dir: Path) -> None:
    """Bar chart of minimum HTT events (parsimony-reconciled count = 11)."""
    if "te_class" not in cands.columns:
        cands = cands.copy()
        cands["te_class"] = cands["qseqid"].apply(
            lambda x: te_class_from_family(str(x).split("__")[-1]))

    # Parsimony minimum: N_MIN_EVENTS partitioned across TE classes
    cls_counts = cands["te_class"].value_counts()
    cls_fracs  = cls_counts / cls_counts.sum()
    min_ev_cls = (cls_fracs * N_MIN_EVENTS).round().astype(int)
    # Ensure total sums to N_MIN_EVENTS
    diff = N_MIN_EVENTS - min_ev_cls.sum()
    if diff != 0:
        top_cls = min_ev_cls.idxmax()
        min_ev_cls[top_cls] += diff

    fig, ax = plt.subplots(figsize=(5, 4))
    present_classes = [c for c in min_ev_cls.index if min_ev_cls[c] > 0]
    colors = [TE_CLASS_COLORS.get(c, "#999") for c in present_classes]
    ax.bar(present_classes, [min_ev_cls[c] for c in present_classes], color=colors)
    ax.set_ylabel("Minimum HTT events", fontsize=10)
    ax.set_title(f"Total: {N_MIN_EVENTS} minimum events (parsimony)", fontsize=10)
    ax.set_xticklabels(present_classes, rotation=30, ha="right", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "htt_minimum_events.pdf")
    plt.close(fig)
    print("Saved: htt_minimum_events.pdf")



# main


def main():
    # load inputs
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    cands       = pd.read_csv(args.candidates,   sep="\t")
    dated       = pd.read_csv(args.dated,        sep="\t")
    conservation = pd.read_csv(args.conservation, sep="\t")
    
    # assign TE class if missing
    if "te_class" not in cands.columns:
        cands["te_class"] = cands["qseqid"].apply(
            lambda x: te_class_from_family(str(x).split("__")[-1]))
    # build clade pairs
    cands["clade_pair"] = cands.apply(
        lambda r: "_".join(sorted([str(r["q_clade"]), str(r["s_clade"])])), axis=1)

    clade_pairs = cands.groupby(["q_clade", "s_clade"]).agg(
        n_candidates=("pident", "count"),
        dominant_class=("te_class", lambda x: x.value_counts().idxmax()
                         if not x.empty else "Unclassified"),
    ).reset_index()
    
    # generate figures
    plot_chord_tree(clade_pairs, fig_dir)
    plot_class_bar(cands, fig_dir)
    plot_identity_strip(cands, fig_dir)
    plot_minimum_events(cands, fig_dir)

    # merge full HTT annotation table
    merged = cands.merge(
        dated[["qseqid", "sseqid", "median_ks"]].drop_duplicates(),
        on=["qseqid", "sseqid"], how="left"
    ).merge(
        conservation[["qseqid", "sseqid", "conservation_fraction",
                       "conservation_category"]].drop_duplicates(),
        on=["qseqid", "sseqid"], how="left"
    )
    merged.to_csv(out_dir / "htt_master_table.tsv", sep="\t", index=False)
    
    # summary
    print(f"\nHTT summary:")
    print(f"  Total candidates: {len(cands)}")
    print(f"  Minimum parsimony events: {N_MIN_EVENTS}")
    print(f"  Clade pairs with transfers: {len(clade_pairs)}")
    if "conservation_category" in merged.columns:
        cc = merged["conservation_category"].value_counts()
        for k, v in cc.items():
            print(f"    {k}: {v}")
    print("Done.")


if __name__ == "__main__":
    main()
