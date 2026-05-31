#!/usr/bin/env python3
"""plot_ancestral_gains_losses.py — Phylogeny with TE family gains/losses bubbles.

Produces a publication figure:
  - Clade color strip
  - Phylogram (branch lengths proportional)
  - Gain bubbles (green) and loss bubbles (red) at each node
  - Present-family count as tip bar chart
  - Internal-node annotations for key nodes

Data sources:
  - guide_clean.treefile         (species phylogeny, 238 tips)
  - wagner_ancestral.txt         (Count AsymmetricWagner output)
  - gains_losses_summary.tsv     (pre-parsed gains/losses/present per node)
"""
from pathlib import Path
import re, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ── Paths ──────────────────────────────────────────────────────────────────
TREE_FILE   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "ancestral_count/05_count/guide_clean.treefile")
WAGNER_FILE = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "ancestral_count/05_count_v2b/wagner_ancestral.txt")
SUMMARY_TSV = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "ancestral_count/05_count_v2b/gains_losses_summary.tsv")
OUTDIR      = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")

# ── Clade colours ──────────────────────────────────────────────────────────
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

# ── Newick parser ──────────────────────────────────────────────────────────
class Node:
    __slots__ = ["name", "length", "children", "x", "y", "count_id",
                 "gains", "losses", "present"]
    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.x = 0.0
        self.y = 0.0
        self.count_id = ""
        self.gains = 0
        self.losses = 0
        self.present = 0


def parse_newick(s):
    s = s.strip().rstrip(";").strip()
    pos = [0]
    def _node():
        node = Node()
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            node.children.append(_node())
            while pos[0] < len(s) and s[pos[0]] == ",":
                pos[0] += 1
                node.children.append(_node())
            if pos[0] < len(s) and s[pos[0]] == ")":
                pos[0] += 1
        name_chars = []
        while pos[0] < len(s) and s[pos[0]] not in ":,();":
            name_chars.append(s[pos[0]])
            pos[0] += 1
        node.name = "".join(name_chars).strip("'")
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


def get_all_nodes(node):
    """Return all nodes in pre-order traversal."""
    out = [node]
    for c in node.children:
        out.extend(get_all_nodes(c))
    return out


# ── Assign Count internal node IDs ────────────────────────────────────────
def assign_count_ids(root, header_nodes):
    """Count.jar numbers nodes in a specific order: leaves first (in tree
    traversal order), then internal nodes in post-order.
    
    We match by parsing the header from the Wagner output and mapping to our
    tree via the same traversal order.
    """
    # Leaves in header order should match tree-traversal-order leaves
    leaves = get_leaves(root)
    leaf_names = [l.name for l in leaves]
    
    # The header order: first come 238 tips (in tree pre-order), then 237 internal
    # nodes (in post-order). Count assigns $N_noname to internal nodes starting
    # from N = n_tips (0-indexed).
    n_tips = len(leaves)
    
    # Map tip names
    for leaf in leaves:
        leaf.count_id = leaf.name
    
    # Map internal nodes: Count does post-order traversal (left-right-root)
    internal_counter = [n_tips]
    def _postorder_assign(node):
        if not node.children:
            return
        for c in node.children:
            _postorder_assign(c)
        node.count_id = f"${internal_counter[0]}_noname"
        internal_counter[0] += 1
    _postorder_assign(root)
    
    # Verify against header
    header_internal = [h for h in header_nodes if h.startswith("$")]
    internal_in_tree = []
    def _collect_postorder(node):
        if not node.children:
            return
        for c in node.children:
            _collect_postorder(c)
        internal_in_tree.append(node.count_id)
    _collect_postorder(root)
    
    matched = sum(1 for a, b in zip(internal_in_tree, header_internal) if a == b)
    print(f"[INFO] Internal node matching: {matched}/{len(header_internal)}")
    if matched < len(header_internal) * 0.9:
        print("[WARN] Poor matching — node IDs may not align correctly")


# ── Phylogram layout ──────────────────────────────────────────────────────
def assign_phylogram_positions(node, y_counter):
    """Tips evenly spaced; x = cumulative branch length from root."""
    if not node.children:
        node.y = y_counter[0]
        y_counter[0] += 1
    else:
        for c in node.children:
            c.x = node.x + c.length
            assign_phylogram_positions(c, y_counter)
        node.y = np.mean([c.y for c in node.children])


# ── Cladogram layout (for cleaner look) ──────────────────────────────────
def assign_cladogram_positions(node, y_counter, depth=0, max_depth=None):
    if not node.children:
        node.x = max_depth
        node.y = y_counter[0]
        y_counter[0] += 1
    else:
        node.x = depth
        for c in node.children:
            assign_cladogram_positions(c, y_counter, depth + 1, max_depth)
        node.y = np.mean([c.y for c in node.children])


def get_max_depth(node, depth=0):
    if not node.children:
        return depth
    return max(get_max_depth(c, depth + 1) for c in node.children)


# ── Draw tree ─────────────────────────────────────────────────────────────
def draw_tree(ax, node, cladogram=True):
    def traverse(n, parent_x=None):
        if parent_x is not None:
            ax.plot([parent_x, n.x], [n.y, n.y], "k-", lw=0.3,
                    solid_capstyle="round")
        if n.children:
            ys = [c.y for c in n.children]
            ax.plot([n.x, n.x], [min(ys), max(ys)], "k-", lw=0.3,
                    solid_capstyle="round")
            for c in n.children:
                traverse(c, n.x)
    traverse(node)


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── Parse tree ────────────────────────────────────────────────────────
    tree_str = TREE_FILE.read_text().strip()
    root = parse_newick(tree_str)
    root.x = 0.0

    leaves = get_leaves(root)
    n_tips = len(leaves)

    # ── Load Wagner header to get node order ──────────────────────────────
    with open(WAGNER_FILE) as f:
        for line in f:
            if line.startswith("# FAMILY\tname"):
                cols = line.rstrip("\n").split("\t")
                header_nodes = cols[2:-4]  # skip Family, name, and last 4 summary cols
                break

    # ── Assign Count IDs ──────────────────────────────────────────────────
    assign_count_ids(root, header_nodes)

    # ── Load gains/losses data ────────────────────────────────────────────
    df = pd.read_csv(SUMMARY_TSV, sep="\t")
    gains_map = dict(zip(df["node"], df["family_gain"]))
    losses_map = dict(zip(df["node"], df["family_loss"]))
    present_map = dict(zip(df["node"], df["family_present"]))

    # Also parse PRESENT section from Wagner output for root node
    # (root has no CHANGE line, only PRESENT)
    with open(WAGNER_FILE) as f:
        for line in f:
            if line.startswith("# PRESENT"):
                parts = line.strip().split("\t")
                if len(parts) >= 4 and parts[1] not in present_map:
                    if parts[3] == "families": continue
                    present_map[parts[1]] = int(parts[3])  # families column

    # Assign to tree nodes
    all_nodes = get_all_nodes(root)
    mapped = 0
    for node in all_nodes:
        cid = node.count_id
        if cid in gains_map or cid in present_map:
            node.gains = gains_map.get(cid, 0)
            node.losses = losses_map.get(cid, 0)
            node.present = present_map.get(cid, 0)
            mapped += 1
    print(f"[INFO] Mapped gains/losses to {mapped}/{len(all_nodes)} nodes")

    # ── Layout ────────────────────────────────────────────────────────────
    max_d = get_max_depth(root)
    assign_cladogram_positions(root, [0], depth=0, max_depth=max_d)

    # ── Figure ────────────────────────────────────────────────────────────
    fig_height = max(16, n_tips * 0.065)
    fig, axes = plt.subplots(1, 5, figsize=(24, fig_height),
                              gridspec_kw={"width_ratios": [0.3, 4, 2.5, 2.5, 2],
                                           "wspace": 0.03})
    ax_clade, ax_tree, ax_gains, ax_losses, ax_present = axes

    # ── Clade strip ───────────────────────────────────────────────────────
    ax_clade.set_xlim(0, 1)
    ax_clade.set_ylim(-0.5, n_tips - 0.5)
    ax_clade.invert_yaxis()
    ax_clade.axis("off")
    for leaf in leaves:
        y = leaf.y
        clade = leaf.name[:2]
        color = CLADE_COLORS.get(clade, "#ddd")
        ax_clade.barh(y, 1, height=0.9, color=color, edgecolor="none")

    # ── Draw tree ─────────────────────────────────────────────────────────
    draw_tree(ax_tree, root)
    ax_tree.set_xlim(-0.5, max_d + 0.5)
    ax_tree.set_ylim(-0.5, n_tips - 0.5)
    ax_tree.invert_yaxis()
    for sp in ax_tree.spines.values():
        sp.set_visible(False)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])

    # ── Bubble-size scaling ───────────────────────────────────────────────
    # Gains and losses use sqrt-proportional bubble sizes
    all_gains = [n.gains for n in all_nodes if n.gains > 0]
    all_losses = [n.losses for n in all_nodes if n.losses > 0]
    max_gain = max(all_gains) if all_gains else 1
    max_loss = max(all_losses) if all_losses else 1

    def gain_size(g):
        if g <= 0:
            return 0
        return max(4, np.sqrt(g / max_gain) * 200)

    def loss_size(l):
        if l <= 0:
            return 0
        return max(4, np.sqrt(l / max_loss) * 200)

    # ── Gain bubbles on tree + dedicated axis ─────────────────────────────
    # Draw gain/loss bubbles directly on the tree for internal nodes
    for node in all_nodes:
        if node.gains > 0:
            ax_tree.scatter(node.x, node.y, s=gain_size(node.gains),
                          c="#2ca02c", alpha=0.5, edgecolor="#1a6f1a",
                          linewidth=0.3, zorder=5)
        if node.losses > 0:
            ax_tree.scatter(node.x, node.y, s=loss_size(node.losses),
                          c="#d62728", alpha=0.5, edgecolor="#8b1a1a",
                          linewidth=0.3, zorder=5)

    # ── Gains bar chart (tips only) ───────────────────────────────────────
    ax_gains.set_ylim(-0.5, n_tips - 0.5)
    ax_gains.invert_yaxis()
    ax_gains.set_xlabel("Family gains", fontsize=8)
    ax_gains.spines["top"].set_visible(False)
    ax_gains.spines["right"].set_visible(False)
    ax_gains.spines["left"].set_visible(False)
    ax_gains.set_yticks([])
    ax_gains.tick_params(axis='x', labelsize=6)

    max_tip_gain = max((l.gains for l in leaves), default=1)
    for leaf in leaves:
        if leaf.gains > 0:
            clade = leaf.name[:2]
            ax_gains.barh(leaf.y, leaf.gains, height=0.7,
                         color=CLADE_COLORS.get(clade, "#aaa"),
                         edgecolor="none", alpha=0.85)
    ax_gains.set_xlim(0, max_tip_gain * 1.05)

    # ── Losses bar chart (tips only) ──────────────────────────────────────
    ax_losses.set_ylim(-0.5, n_tips - 0.5)
    ax_losses.invert_yaxis()
    ax_losses.set_xlabel("Family losses", fontsize=8)
    ax_losses.spines["top"].set_visible(False)
    ax_losses.spines["right"].set_visible(False)
    ax_losses.spines["left"].set_visible(False)
    ax_losses.set_yticks([])
    ax_losses.tick_params(axis='x', labelsize=6)

    max_tip_loss = max((l.losses for l in leaves), default=1)
    for leaf in leaves:
        if leaf.losses > 0:
            ax_losses.barh(leaf.y, leaf.losses, height=0.7,
                          color="#d62728", edgecolor="none", alpha=0.7)
    ax_losses.set_xlim(0, max(max_tip_loss * 1.05, 5))

    # ── Present families bar chart ────────────────────────────────────────
    ax_present.set_ylim(-0.5, n_tips - 0.5)
    ax_present.invert_yaxis()
    ax_present.set_xlabel("FL families present", fontsize=8)
    ax_present.spines["top"].set_visible(False)
    ax_present.spines["right"].set_visible(False)
    ax_present.spines["left"].set_visible(False)
    ax_present.set_yticks([])
    ax_present.tick_params(axis='x', labelsize=6)

    max_present = max((l.present for l in leaves), default=1)
    for leaf in leaves:
        if leaf.present > 0:
            clade = leaf.name[:2]
            ax_present.barh(leaf.y, leaf.present, height=0.7,
                           color=CLADE_COLORS.get(clade, "#aaa"),
                           edgecolor="none", alpha=0.6)
    ax_present.set_xlim(0, max_present * 1.05)

    # ── Annotate key internal nodes ───────────────────────────────────────
    # Find the root and nodes with highest gains
    root_node = root
    if root_node.present > 0:
        ax_tree.annotate(f"Root\n{root_node.present} fam",
                        xy=(root_node.x, root_node.y),
                        fontsize=5, ha="right", va="center",
                        color="#333", style="italic")

    # Annotate top 10 internal nodes by gains (excluding tips)
    internal = [n for n in all_nodes if n.children and n.gains >= 50]
    internal.sort(key=lambda n: -n.gains)
    for i, node in enumerate(internal[:10]):
        ax_tree.annotate(f"+{node.gains}",
                        xy=(node.x, node.y), xytext=(node.x - 0.3, node.y),
                        fontsize=4, ha="right", va="center",
                        color="#1a6f1a", fontweight="bold",
                        arrowprops=dict(arrowstyle="-", lw=0.3, color="#1a6f1a")
                        if i < 5 else None)

    # ── Legends ───────────────────────────────────────────────────────────
    # Size legend for bubbles
    ref_sizes = [10, 100, 500, 2000]
    gain_handles = []
    for rs in ref_sizes:
        gain_handles.append(Line2D([0], [0], marker='o', color='w',
                            markerfacecolor='#2ca02c', markersize=np.sqrt(gain_size(rs))/2,
                            alpha=0.5, markeredgecolor='#1a6f1a',
                            markeredgewidth=0.3, label=f'{rs}'))
    loss_handles = []
    for rs in ref_sizes[:3]:  # fewer loss reference
        loss_handles.append(Line2D([0], [0], marker='o', color='w',
                            markerfacecolor='#d62728', markersize=np.sqrt(loss_size(rs))/2,
                            alpha=0.5, markeredgecolor='#8b1a1a',
                            markeredgewidth=0.3, label=f'{rs}'))

    leg1 = ax_tree.legend(handles=gain_handles, fontsize=5, loc="lower left",
                          title="Gains", title_fontsize=6, framealpha=0.9,
                          labelspacing=1.2)
    ax_tree.add_artist(leg1)
    ax_tree.legend(handles=loss_handles, fontsize=5, loc="lower right",
                   title="Losses", title_fontsize=6, framealpha=0.9,
                   labelspacing=1.2)

    # Clade legend
    present_clades = sorted(set(l.name[:2] for l in leaves))
    clade_handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#ddd"),
                     label=CLADE_NAMES.get(c, c))
                     for c in present_clades if c in CLADE_COLORS]
    fig.legend(handles=clade_handles, fontsize=5, loc="upper left",
               bbox_to_anchor=(0.01, 0.99), ncol=2, framealpha=0.9,
               title="Clades", title_fontsize=6)

    # ── Title & save ──────────────────────────────────────────────────────
    fig.suptitle("Ancestral reconstruction of TE family gains and losses\n"
                 f"(AsymmetricWagner, gain penalty=2.0 | {n_tips} species, "
                 f"45,840 families | Root: {root_node.present} families)",
                 fontsize=11, y=0.998)
    fig.tight_layout(rect=[0, 0, 1, 0.993])

    out_pdf = OUTDIR / "ancestral_gains_losses_tree.pdf"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "ancestral_gains_losses_tree.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[DONE] {out_pdf}")

    # ── Summary stats ─────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"Total gains:  {sum(n.gains for n in all_nodes):,}")
    print(f"Total losses: {sum(n.losses for n in all_nodes):,}")
    print(f"Root families: {root_node.present}")
    print(f"Tip with most gains: "
          f"{max(leaves, key=lambda l: l.gains).name} "
          f"({max(leaves, key=lambda l: l.gains).gains})")
    print(f"Tip with most present: "
          f"{max(leaves, key=lambda l: l.present).name} "
          f"({max(leaves, key=lambda l: l.present).present})")


if __name__ == "__main__":
    main()
