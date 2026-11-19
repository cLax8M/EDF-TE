#!/usr/bin/env python3

# builds TE composition cladogram with Shannon diversity and Recency Index per species

# imports
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

try:
    from ete3 import Tree
except ImportError:
    print("ete3 not installed: pip install ete3")
    sys.exit(1)


# define clade metadata (color + naming scheme)
CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]


# define TE class color scheme
CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

CLADE_NAMES = {
    "HO": "Holozoa", "NU": "Nucleariida", "AP": "Archaeplastida",
    "CR": "Cryptomycetes", "MI": "Microsporidia", "BL": "Blastocladiomycota",
    "OL": "Olpidiales", "CH": "Chytridiomycetes", "NE": "Neocallimastigomycota",
    "BS": "Bacillariophyta", "EN": "Entomophthorales", "ZO": "Zoopagomycota",
    "KX": "Kickxellomycotina", "GL": "Glomerales", "MO": "Mortierellomycota",
    "MU": "Mucoromycotina", "BA": "Basidiomycota", "AS": "Ascomycota",
}

# Neocallimastigomycota species are annotated with the CH_ prefix but belong
# to a separate clade
NE_SPECIES = frozenset({
    "CH_Anasp1", "CH_Neocon1", "CH_Neolan1",
    "CH_Orpsp1", "CH_Pecora1", "CH_Caecom1", "CH_Pirfi3",
})

TE_CLASS_COLORS = {
    "LTR":             "#e41a1c",
    "DNA":             "#377eb8",
    "LINE":            "#4daf4a",
    "SINE":            "#984ea3",
    "Rolling Circle":  "#ff7f00",
    "Penelope":        "#a65628",
    "Unclassified":    "#999999",
}

# define TE class order for plotting
TE_CLASS_ORDER = ["LTR", "DNA", "LINE", "SINE", "Rolling Circle",
                  "Penelope", "Unclassified"]


# parse arguments
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tree",       required=True)
    p.add_argument("--te_matrix",  required=True)
    p.add_argument("--shannon_tsv",required=True)
    p.add_argument("--kimura_tsv", required=True)
    p.add_argument("--clade_map",  required=True)
    p.add_argument("--out_dir",    default="figures")
    return p.parse_args()





# define fallback clade assignment if not in map
def get_clade(species: str) -> str:
    if species in NE_SPECIES:
        return "NE"
    return species[:2] if len(species) >= 2 else "??"

# load input tables and merge metadata
def load_data(args):
    clade_df = pd.read_csv(args.clade_map, sep="\t", header=None,
                           names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))

    te = pd.read_csv(args.te_matrix, sep="\t")
    te["species"] = te["species"].str.replace(r"\.aa$", "", regex=True)
    te["clade"] = te["species"].apply(
        lambda s: clade_map.get(s, get_clade(s)))

    shannon = pd.read_csv(args.shannon_tsv, sep="\t")
    shannon["species"] = shannon["species"].str.replace(r"\.aa$", "", regex=True)

    kimura = pd.read_csv(args.kimura_tsv, sep="\t")
    kimura["species"] = kimura["species"].str.replace(r"\.aa$", "", regex=True)

    df = te.merge(shannon[["species", "shannon_H"]], on="species", how="left")
    df = df.merge(kimura[["species", "pct_recent"]], on="species", how="left")

    # Name remapping (tree tips may differ from data labels)
    name_map = {}
    return df, name_map



# extract species order from phylogenetic tree
def tree_species_order(tree_file: str) -> list[str]:
    t = Tree(tree_file, format=1)
    return [leaf.name for leaf in t.iter_leaves()]


# Figure

def draw_figure(df: pd.DataFrame, tree_order: list[str], out_dir: Path) -> None:
    # Reindex by tree order
    df = df.set_index("species")
    species_in_plot = [s for s in tree_order if s in df.index]
    df = df.loc[species_in_plot]
    n = len(species_in_plot)

    # create figure layout
    fig = plt.figure(figsize=(18, max(12, n * 0.22)))
    gs = GridSpec(1, 4, figure=fig, width_ratios=[0.5, 3.5, 0.9, 0.9],
                  wspace=0.05)

    ax_clade = fig.add_subplot(gs[0])  # clade color strip
    ax_bars  = fig.add_subplot(gs[1])  # TE class stacked bars
    ax_shan  = fig.add_subplot(gs[2])  # Shannon H
    ax_ri    = fig.add_subplot(gs[3])  # Recency Index

    y_pos = np.arange(n)

    # ── Clade colour strip ────────────────────────────────────────────────────
    for i, sp in enumerate(species_in_plot):
        clade = df.at[sp, "clade"] if "clade" in df.columns else get_clade(sp)
        color = CLADE_COLORS.get(clade, "#cccccc")
        ax_clade.barh(i, 1, left=0, height=0.85, color=color)

    ax_clade.set_xlim(0, 1)
    ax_clade.set_yticks(y_pos)
    ax_clade.set_yticklabels(species_in_plot, fontsize=5)
    ax_clade.set_xticks([])
    ax_clade.spines[:].set_visible(False)

    # draw TE composition stacked bars
    lefts = np.zeros(n)
    for te_class in TE_CLASS_ORDER:
        col = te_class if te_class in df.columns else None
        if col is None:
            continue
        vals = df[col].fillna(0).values
        ax_bars.barh(y_pos, vals, left=lefts, height=0.75,
                     color=TE_CLASS_COLORS.get(te_class, "#bbbbbb"),
                     label=te_class)
        lefts += vals

    ax_bars.set_xlabel("TE content (% genome)", fontsize=9)
    ax_bars.set_yticks([])
    ax_bars.spines[["top", "right", "left"]].set_visible(False)

    # draw Shannon diversity
    clades_arr = df["clade"].values if "clade" in df.columns else \
                 [get_clade(s) for s in species_in_plot]
    colors_arr = [CLADE_COLORS.get(c, "#999") for c in clades_arr]
    shan_vals = df["shannon_H"].fillna(0).values

    ax_shan.barh(y_pos, shan_vals, height=0.75, color=colors_arr, alpha=0.8)
    ax_shan.set_xlabel("Shannon H", fontsize=9)
    ax_shan.set_yticks([])
    ax_shan.spines[["top", "right", "left"]].set_visible(False)

    # draw Recency Index 
    ri_vals = df["pct_recent"].fillna(0).values
    ax_ri.barh(y_pos, ri_vals, height=0.75, color=colors_arr, alpha=0.8)
    ax_ri.set_xlabel("Recency Index", fontsize=9)
    ax_ri.set_yticks([])
    ax_ri.spines[["top", "right", "left"]].set_visible(False)

    # TE class legend 
    patches = [mpatches.Patch(color=TE_CLASS_COLORS[c], label=c)
               for c in TE_CLASS_ORDER if c in df.columns]
    ax_bars.legend(handles=patches, loc="lower right", frameon=False,
                   fontsize=7, ncol=1)

    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "cladogram_te_composition.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "cladogram_te_composition.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: cladogram_te_composition.pdf ({n} species)")



# main

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)

    print("Loading data ...")
    df, _ = load_data(args)
    print(f"  {len(df)} species with TE data")

    print(f"Reading tree: {args.tree}")
    tree_order = tree_species_order(args.tree)
    print(f"  {len(tree_order)} tips in tree")

    print("Drawing figure ...")
    draw_figure(df, tree_order, out_dir)


if __name__ == "__main__":
    main()
