#!/usr/bin/env python3

# generates a summary figure with phylogeny, WGD calls,
# TE content and synteny multiplicons

# imports
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

# parse arguments
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species_tree",  required=True)
    p.add_argument("--wgd_summary",   required=True)
    p.add_argument("--te_content",    required=True)
    p.add_argument("--multiplicons",  required=True)
    p.add_argument("--clade_map",     required=True)
    p.add_argument("--out_dir",       default="wgd_viz_out")
    return p.parse_args()

# get species order from the phylogeny
def get_tree_leaf_order(nwk: str, species_list: set[str]) -> list[str]:
    try:
        from ete3 import Tree
        t = Tree(nwk, format=1)
        t.prune([l for l in t.get_leaves() if l.name in species_list],
                preserve_branch_length=True)
        return [l.name for l in t.get_leaves()]
    except Exception:
        return sorted(species_list)

# main
def main():
# parse args
    args = parse_args()
    
    # define output dirs
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # load metadata
    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))

    # load input tables
    wgd   = pd.read_csv(args.wgd_summary, sep="\t").set_index("species")
    te_c  = pd.read_csv(args.te_content, sep="\t").set_index("species")
    multi = pd.read_csv(args.multiplicons, sep="\t").set_index("species")

    # define species order from phylogeny
    all_sp = set(wgd.index) & set(te_c.index)
    nwk    = open(args.species_tree).read()
    sp_order = get_tree_leaf_order(nwk, all_sp)

    # align data to leaf order (species order)
    te_vals    = [te_c.loc[sp, "te_pct_genome"] if sp in te_c.index else 0.0
                  for sp in sp_order]
    wgd_calls  = [bool(wgd.loc[sp, "wgd_call"])  if sp in wgd.index else False
                  for sp in sp_order]
    multi_vals = [int(multi.loc[sp, "n_multiplicons"]) if sp in multi.index else 0
                  for sp in sp_order]
    
    # define species colors
    clades = [clade_map.get(sp, "MU") for sp in sp_order]
    colors = [CLADE_COLORS.get(c, "#cccccc") for c in clades]

    # inicialize figure
    n = len(sp_order)
    fig = plt.figure(figsize=(16, max(10, n * 0.25 + 3)))
    gs  = gridspec.GridSpec(1, 3, width_ratios=[3, 2, 1.5], figure=fig)
    ax_tree = fig.add_subplot(gs[0])
    ax_te   = fig.add_subplot(gs[1])
    ax_mult = fig.add_subplot(gs[2])

    # Panel A: phylogeny and WGD calls
    y_positions = {sp: i for i, sp in enumerate(sp_order)}
    for sp, y in y_positions.items():
        clade  = clade_map.get(sp, "MU")
        color  = CLADE_COLORS.get(clade, "#cccccc")
        ax_tree.barh(y, 1, left=0, height=0.8, color=color, alpha=0.5)
        ax_tree.text(-0.05, y, sp.replace("_", " "), ha="right", va="center",
                     fontsize=5, color=color)
        if wgd_calls[list(sp_order).index(sp)]:
            ax_tree.scatter([1.05], [y], marker="*", color="gold",
                            s=60, zorder=5)
    ax_tree.set_xlim(-0.2, 1.3)
    ax_tree.axis("off")
    ax_tree.set_title("Phylogeny", fontsize=10)

    # Panel B: TE content
    ax_te.barh(range(n), te_vals, color=colors, alpha=0.75)
    ax_te.set_yticks(range(n))
    ax_te.set_yticklabels([])
    ax_te.set_xlabel("TE content (% genome)", fontsize=9)
    ax_te.set_title("TE content", fontsize=10)
    ax_te.spines[["top", "right"]].set_visible(False)
    ax_te.set_ylim(-0.5, n - 0.5)

    # Panel C: synteny multiplicons 
    ax_mult.scatter(multi_vals, range(n), color=colors, s=20, alpha=0.75)
    ax_mult.set_yticks([])
    ax_mult.set_xlabel("Multiplicons", fontsize=9)
    ax_mult.set_title("Synteny blocks", fontsize=10)
    ax_mult.spines[["top", "right"]].set_visible(False)
    ax_mult.set_ylim(-0.5, n - 0.5)

    # Generate legend
    seen = {}
    for c in set(clades):
        if c not in seen:
            seen[c] = mpatches.Patch(color=CLADE_COLORS.get(c, "#ccc"), label=c)
    seen["WGD"] = mpatches.Patch(facecolor="none", edgecolor="none",
                                  label="")  # placeholder
    wgd_marker = plt.scatter([], [], marker="*", color="gold", s=80, label="WGD")
    fig.legend(handles=list(seen.values()) + [wgd_marker],
               fontsize=7, ncol=4, loc="lower center", frameon=False)

    fig.tight_layout(rect=[0, 0.05, 1, 1])
    
    # save figure
    fig.savefig(fig_dir / "wgd_te_tree_figure.pdf")
    fig.savefig(fig_dir / "wgd_te_tree_figure.png", dpi=200)
    plt.close(fig)
    print("Saved: wgd_te_tree_figure.pdf")
    print("Done.")


if __name__ == "__main__":
    main()
