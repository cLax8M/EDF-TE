#!/usr/bin/env python3
"""plot_cladogram_te_types.py — Cladogram + stacked TE-type bars + Shannon diversity.

Produces a publication figure with:
  - Clade color strip
  - Cladogram (equal branch lengths / all tips aligned)
  - Stacked horizontal bars: LTR, DNA, LINE, SINE, Rolling Circle, Penelope, Unclassified
  - Shannon H' bubbles

Data sources:
  - guide_clean.treefile  (species phylogeny)
  - class_level_matrix.tsv (TE composition per species, proportions within TE)
  - figure_species_data.tsv (pct_TE_total per species for genome %)
  - B01_shannon_diversity.tsv (H' per species)
"""
from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────
TREE_FILE   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/ancestral_count/05_count/guide_clean.treefile")
CLASS_MAT   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures/class_level_matrix.tsv")
SPECIES_DAT = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figure_species_data.tsv")
SHANNON_TSV = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/B01_shannon_diversity.tsv")
OUTDIR      = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")

# ── Clade config (same colours throughout project) ─────────────────────────
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

# ── TE type config (matching user's reference figure) ──────────────────────
TE_ORDER = ["LTR", "DNA", "LINE", "SINE", "Rolling Circle", "Penelope", "Unclassified"]
TE_COLORS = {
    "LTR":            "#e41a1c",  # red
    "DNA":            "#377eb8",  # blue
    "LINE":           "#4daf4a",  # green
    "SINE":           "#f1c232",  # yellow-gold
    "Rolling Circle": "#ff7f00",  # orange
    "Penelope":       "#66c2a5",  # teal
    "Unclassified":   "#bdbdbd",  # grey
}

# ── Newick parser ──────────────────────────────────────────────────────────
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

    def _node():
        node = Node()
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            children = [_node()]
            while pos[0] < len(s) and s[pos[0]] == ",":
                pos[0] += 1
                children.append(_node())
            if pos[0] < len(s) and s[pos[0]] == ")":
                pos[0] += 1
            node.children = children
        name_chars = []
        while pos[0] < len(s) and s[pos[0]] not in ":,();":
            name_chars.append(s[pos[0]])
            pos[0] += 1
        node.name = "".join(name_chars)
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            bl = []
            while pos[0] < len(s) and s[pos[0]] not in ",();":
                bl.append(s[pos[0]])
                pos[0] += 1
            try:
                node.length = float("".join(bl))
            except ValueError:
                node.length = 0.0
        return node

    return _node()


def get_leaves(node):
    if not node.children:
        return [node]
    out = []
    for c in node.children:
        out.extend(get_leaves(c))
    return out


def get_max_depth(node, depth=0):
    if not node.children:
        return depth
    return max(get_max_depth(c, depth + 1) for c in node.children)


def assign_cladogram_positions(node, y_counter, depth=0, max_depth=None):
    """All tips aligned at x = max_depth, internal nodes at x = their depth."""
    if not node.children:
        node.x = max_depth
        node.y = y_counter[0]
        y_counter[0] += 1
    else:
        node.x = depth
        for c in node.children:
            assign_cladogram_positions(c, y_counter, depth + 1, max_depth)
        node.y = np.mean([c.y for c in node.children])


def draw_cladogram(ax, node):
    """Draw rectangular cladogram. Returns {tip_name: y_pos}."""
    def traverse(n, parent_x=None):
        if parent_x is not None:
            # horizontal line from parent to this node
            ax.plot([parent_x, n.x], [n.y, n.y], "k-", lw=0.3, solid_capstyle="round")
        if n.children:
            ys = [c.y for c in n.children]
            # vertical line connecting children
            ax.plot([n.x, n.x], [min(ys), max(ys)], "k-", lw=0.3, solid_capstyle="round")
            for c in n.children:
                traverse(c, n.x)
    traverse(node)
    return {leaf.name: leaf.y for leaf in get_leaves(node)}


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Load tree
    tree_str = TREE_FILE.read_text().strip()
    root = parse_newick(tree_str)
    max_d = get_max_depth(root)
    assign_cladogram_positions(root, [0], depth=0, max_depth=max_d)
    leaves = get_leaves(root)
    n_tips = len(leaves)
    print(f"[INFO] Tree has {n_tips} tips, max depth {max_d}")

    # Load TE class composition
    df_class = pd.read_csv(CLASS_MAT, sep="\t")
    class_map = {}
    for _, row in df_class.iterrows():
        class_map[row["species"]] = {t: float(row.get(t, 0) or 0) for t in TE_ORDER}

    # Load pct_TE_total for genome % conversion
    df_sp = pd.read_csv(SPECIES_DAT, sep="\t")
    te_total_map = dict(zip(df_sp["taxon"], df_sp["pct_TE_total"]))

    # Load Shannon
    df_sh = pd.read_csv(SHANNON_TSV, sep="\t")
    shannon_map = dict(zip(df_sh["species"], df_sh["H_prime"]))

    # ── Figure layout ────────────────────────────────────────────────────
    fig_height = max(14, n_tips * 0.065)
    fig, axes = plt.subplots(1, 4, figsize=(18, fig_height),
                              gridspec_kw={"width_ratios": [0.3, 3, 3.5, 1.2],
                                           "wspace": 0.02})
    ax_clade, ax_tree, ax_bars, ax_shan = axes

    # ── Draw cladogram ───────────────────────────────────────────────────
    tip_pos = draw_cladogram(ax_tree, root)

    ax_tree.set_xlim(-0.5, max_d + 0.5)
    ax_tree.set_ylim(-0.5, n_tips - 0.5)
    ax_tree.invert_yaxis()
    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)
    ax_tree.spines["left"].set_visible(False)
    ax_tree.spines["bottom"].set_visible(False)
    ax_tree.set_xticks([])
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

    # ── Stacked TE bars (genome %) ───────────────────────────────────────
    max_te = 0
    bar_data = {}
    for leaf in leaves:
        sp = leaf.name
        te_total = te_total_map.get(sp, 0)
        if sp in class_map and te_total > 0:
            vals = {}
            for t in TE_ORDER:
                vals[t] = class_map[sp].get(t, 0) * te_total
            bar_data[sp] = vals
            total = sum(vals.values())
            if total > max_te:
                max_te = total

    ax_bars.set_xlim(0, min(max_te * 1.1, 60))
    ax_bars.set_ylim(-0.5, n_tips - 0.5)
    ax_bars.invert_yaxis()
    ax_bars.set_xlabel("% genome", fontsize=8)
    ax_bars.spines["top"].set_visible(False)
    ax_bars.spines["right"].set_visible(False)
    ax_bars.spines["left"].set_visible(False)
    ax_bars.set_yticks([])

    for leaf in leaves:
        y = tip_pos[leaf.name]
        sp = leaf.name
        if sp not in bar_data:
            continue
        left = 0
        for t in TE_ORDER:
            val = bar_data[sp][t]
            if val > 0:
                ax_bars.barh(y, val, left=left, height=0.7,
                            color=TE_COLORS[t], edgecolor="none")
                left += val

    # ── Shannon bubbles ──────────────────────────────────────────────────
    ax_shan.set_ylim(-0.5, n_tips - 0.5)
    ax_shan.invert_yaxis()
    ax_shan.set_xlabel("Shannon H'", fontsize=8)
    ax_shan.spines["top"].set_visible(False)
    ax_shan.spines["right"].set_visible(False)
    ax_shan.spines["left"].set_visible(False)
    ax_shan.set_yticks([])

    h_vals = []
    for leaf in leaves:
        sp = leaf.name
        h = shannon_map.get(sp, np.nan)
        if not pd.isna(h):
            h_vals.append(h)
    max_h = max(h_vals) if h_vals else 3

    for leaf in leaves:
        y = tip_pos[leaf.name]
        sp = leaf.name
        h = shannon_map.get(sp, np.nan)
        if pd.isna(h):
            continue
        clade = sp[:2]
        ax_shan.scatter(h, y, s=max(8, h * 25), alpha=0.7,
                       c=CLADE_COLORS.get(clade, "#aaa"),
                       edgecolor="k", linewidth=0.2, zorder=3)

    ax_shan.set_xlim(0, max_h * 1.1)

    # ── Legends ──────────────────────────────────────────────────────────
    # TE type legend
    te_handles = [mpatches.Patch(color=TE_COLORS[t], label=t) for t in TE_ORDER]
    ax_bars.legend(handles=te_handles, fontsize=5, loc="lower right",
                   framealpha=0.9, title="TE type", title_fontsize=6)

    # Clade legend
    present_clades = sorted(set(l.name[:2] for l in leaves))
    clade_handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#ddd"),
                     label=CLADE_NAMES.get(c, c))
                     for c in present_clades if c in CLADE_COLORS]
    fig.legend(handles=clade_handles, fontsize=5, loc="upper left",
               bbox_to_anchor=(0.01, 0.99), ncol=2, framealpha=0.9,
               title="Clades", title_fontsize=6)

    fig.suptitle("Cladogram with TE composition and Shannon diversity",
                 fontsize=12, y=0.998)
    fig.tight_layout(rect=[0, 0, 1, 0.995])
    out = OUTDIR / "cladogram_te_types_shannon.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "cladogram_te_types_shannon.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[DONE] {out.name} ({n_tips} species)")


if __name__ == "__main__":
    main()
