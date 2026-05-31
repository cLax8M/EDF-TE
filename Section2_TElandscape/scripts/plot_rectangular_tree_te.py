#!/usr/bin/env python3
"""
plot_rectangular_tree_te.py — Rectangular phylogeny + TE content bars + Shannon diversity.

Produces a publication figure with:
  Panel A: Rectangular phylogenetic tree (guide_clean.treefile)
  Panel B: Stacked horizontal bars of TE content (ClassI, ClassII, Unknown)
  Panel C: Shannon diversity bubbles (H') and evenness (J')
  Clade color strip on the left

Usage:
  python3 plot_rectangular_tree_te.py
"""
from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from io import StringIO

# ── Paths ──────────────────────────────────────────────────────────────────
TREE_FILE  = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count/05_count/guide_clean.treefile")
SPECIES_DATA = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figure_species_data.tsv")
SHANNON_TSV  = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/B01_shannon_diversity.tsv")
OUTDIR       = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")

# ── Clade config ───────────────────────────────────────────────────────────
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

# ── Simple recursive-descent Newick parser ─────────────────────────────────
class Node:
    __slots__ = ["name", "length", "children", "x", "y"]
    def __init__(self, name="", length=0.0, children=None):
        self.name = name
        self.length = length
        self.children = children or []
        self.x = 0.0
        self.y = 0.0

def parse_newick(s):
    s = s.strip().rstrip(";").strip()
    pos = [0]

    def parse_node():
        node = Node()
        if s[pos[0]] == "(":
            pos[0] += 1  # skip (
            children = []
            children.append(parse_node())
            while s[pos[0]] == ",":
                pos[0] += 1
                children.append(parse_node())
            if s[pos[0]] == ")":
                pos[0] += 1
            node.children = children
        # Parse name
        name_chars = []
        while pos[0] < len(s) and s[pos[0]] not in ":,();":
            name_chars.append(s[pos[0]])
            pos[0] += 1
        node.name = "".join(name_chars)
        # Parse branch length
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            bl_chars = []
            while pos[0] < len(s) and s[pos[0]] not in ",();":
                bl_chars.append(s[pos[0]])
                pos[0] += 1
            try:
                node.length = float("".join(bl_chars))
            except ValueError:
                node.length = 0.0
        return node

    return parse_node()


def get_leaves(node):
    if not node.children:
        return [node]
    leaves = []
    for c in node.children:
        leaves.extend(get_leaves(c))
    return leaves


def assign_positions(node, x_offset=0.0, y_counter=None):
    """Assign x (depth) and y (ladder order) to each node."""
    if y_counter is None:
        y_counter = [0]
    node.x = x_offset + node.length
    if not node.children:
        node.y = y_counter[0]
        y_counter[0] += 1
    else:
        for c in node.children:
            assign_positions(c, node.x, y_counter)
        node.y = np.mean([c.y for c in node.children])


def draw_tree(ax, node, tip_names_order):
    """Draw rectangular tree on ax. Returns list of (tip_name, y_position)."""
    lines_h = []  # horizontal (branch)
    lines_v = []  # vertical (connecting children)

    def traverse(n, parent_x=None):
        if parent_x is not None:
            lines_h.append((parent_x, n.x, n.y))
        if n.children:
            ys = [c.y for c in n.children]
            lines_v.append((n.x, min(ys), max(ys)))
            for c in n.children:
                traverse(c, n.x)

    traverse(node)

    for x0, x1, y in lines_h:
        ax.plot([x0, x1], [y, y], "k-", linewidth=0.4, solid_capstyle="round")
    for x, y0, y1 in lines_v:
        ax.plot([x, x], [y0, y1], "k-", linewidth=0.4, solid_capstyle="round")

    tip_positions = {}
    for leaf in get_leaves(node):
        tip_positions[leaf.name] = leaf.y
    return tip_positions


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── Load tree ────────────────────────────────────────────────────────
    tree_str = TREE_FILE.read_text().strip()
    root = parse_newick(tree_str)
    assign_positions(root, x_offset=0.0)
    leaves = get_leaves(root)
    n_tips = len(leaves)
    print(f"[INFO] Tree has {n_tips} tips")

    # ── Load species data ────────────────────────────────────────────────
    df_sp = pd.read_csv(SPECIES_DATA, sep="\t")
    # tree_label has .aa suffix, taxon does not
    sp_data = {}
    for _, row in df_sp.iterrows():
        sp_data[row["taxon"]] = row

    # ── Load Shannon ─────────────────────────────────────────────────────
    df_sh = pd.read_csv(SHANNON_TSV, sep="\t")
    shannon_map = dict(zip(df_sh["species"], df_sh["H_prime"]))
    evenness_map = dict(zip(df_sh["species"], df_sh["J_evenness"]))
    nsf_map = dict(zip(df_sh["species"], df_sh["n_superfamilies"]))

    # ── Set up figure ────────────────────────────────────────────────────
    fig_height = max(12, n_tips * 0.065)
    fig, axes = plt.subplots(1, 4, figsize=(18, fig_height),
                              gridspec_kw={"width_ratios": [0.3, 3, 3, 1.5],
                                           "wspace": 0.02})
    ax_clade = axes[0]  # clade color strip
    ax_tree = axes[1]   # phylogeny
    ax_bars = axes[2]   # TE content bars
    ax_shan = axes[3]   # Shannon diversity

    # ── Draw tree ────────────────────────────────────────────────────────
    tip_pos = draw_tree(ax_tree, root, [l.name for l in leaves])
    max_x = max(l.x for l in leaves)

    ax_tree.set_xlim(-0.05, max_x * 1.05)
    ax_tree.set_ylim(-0.5, n_tips - 0.5)
    ax_tree.invert_yaxis()
    ax_tree.set_xlabel("Substitutions/site", fontsize=8)
    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)
    ax_tree.spines["left"].set_visible(False)
    ax_tree.set_yticks([])

    # ── Clade color strip ────────────────────────────────────────────────
    ax_clade.set_xlim(0, 1)
    ax_clade.set_ylim(-0.5, n_tips - 0.5)
    ax_clade.invert_yaxis()
    ax_clade.axis("off")

    for leaf in leaves:
        y = tip_pos[leaf.name]
        clade = leaf.name[:2]
        color = CLADE_COLORS.get(clade, "#ddd")
        ax_clade.barh(y, 1, height=0.9, color=color, edgecolor="none")

    # ── TE content stacked bars ──────────────────────────────────────────
    ax_bars.set_xlim(0, 50)
    ax_bars.set_ylim(-0.5, n_tips - 0.5)
    ax_bars.invert_yaxis()
    ax_bars.set_xlabel("% genome", fontsize=8)
    ax_bars.set_title("TE content", fontsize=9)
    ax_bars.spines["top"].set_visible(False)
    ax_bars.spines["right"].set_visible(False)
    ax_bars.spines["left"].set_visible(False)
    ax_bars.set_yticks([])

    for leaf in leaves:
        y = tip_pos[leaf.name]
        sp = leaf.name
        if sp in sp_data:
            row = sp_data[sp]
            classI = float(row.get("pct_ClassI", 0) or 0)
            classII = float(row.get("pct_ClassII", 0) or 0)
            unknown = float(row.get("pct_TE_unknown", 0) or 0)
            # Stack: ClassI → ClassII → Unknown
            ax_bars.barh(y, classI, height=0.7, color="#4575b4", edgecolor="none")
            ax_bars.barh(y, classII, left=classI, height=0.7, color="#d73027", edgecolor="none")
            ax_bars.barh(y, unknown, left=classI + classII, height=0.7,
                        color="#999999", edgecolor="none")

    # ── Shannon diversity bubbles ────────────────────────────────────────
    ax_shan.set_ylim(-0.5, n_tips - 0.5)
    ax_shan.invert_yaxis()
    ax_shan.set_xlabel("Shannon H'", fontsize=8)
    ax_shan.set_title("Diversity", fontsize=9)
    ax_shan.spines["top"].set_visible(False)
    ax_shan.spines["right"].set_visible(False)
    ax_shan.spines["left"].set_visible(False)
    ax_shan.set_yticks([])

    for leaf in leaves:
        y = tip_pos[leaf.name]
        sp = leaf.name
        # Shannon species names may differ slightly (e.g. AS_AspfumUD vs AS_AspfumUD1)
        h = shannon_map.get(sp, np.nan)
        if pd.isna(h):
            # Try without trailing digit
            sp_short = re.sub(r'\d+$', '', sp)
            h = shannon_map.get(sp_short, np.nan)
        if not pd.isna(h):
            j = evenness_map.get(sp, evenness_map.get(re.sub(r'\d+$', '', sp), 0.3))
            clade = sp[:2]
            ax_shan.scatter(h, y, s=max(10, j * 80), alpha=0.7,
                          c=CLADE_COLORS.get(clade, "#aaa"),
                          edgecolor="k", linewidth=0.3, zorder=3)

    ax_shan.set_xlim(0, max(shannon_map.values()) * 1.1 if shannon_map else 3)

    # ── Legend ────────────────────────────────────────────────────────────
    # TE content legend
    te_handles = [
        mpatches.Patch(color="#4575b4", label="Class I (retro)"),
        mpatches.Patch(color="#d73027", label="Class II (DNA)"),
        mpatches.Patch(color="#999999", label="Unknown"),
    ]
    ax_bars.legend(handles=te_handles, fontsize=6, loc="lower right", framealpha=0.9)

    # Clade legend (on separate text)
    present_clades = sorted(set(l.name[:2] for l in leaves))
    clade_handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#ddd"),
                     label=CLADE_NAMES.get(c, c))
                     for c in present_clades if c in CLADE_COLORS]
    fig.legend(handles=clade_handles, fontsize=5, loc="upper left",
               bbox_to_anchor=(0.01, 0.99), ncol=2, framealpha=0.9,
               title="Clades", title_fontsize=6)

    fig.suptitle("Rectangular phylogeny with TE landscape and Shannon diversity",
                 fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(OUTDIR / "rectangular_tree_te_shannon.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "rectangular_tree_te_shannon.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[DONE] rectangular_tree_te_shannon.pdf ({n_tips} species)")


if __name__ == "__main__":
    main()
