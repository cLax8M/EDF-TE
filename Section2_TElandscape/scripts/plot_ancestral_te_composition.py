#!/usr/bin/env python3
"""plot_ancestral_te_composition.py — Phylogeny with TE class composition at ancestral nodes.

Shows a cladogram where each internal (ancestral) node has a small stacked
horizontal bar indicating the proportion of DNA, LTR, LINE, SINE, RC, and
Unknown TE families inferred to be present at that ancestor by
AsymmetricWagner parsimony.

Data sources:
  - guide_clean.treefile                   (species phylogeny, 238 tips)
  - wagner_ancestral.txt                   (Count AsymmetricWagner full output)
  - combined_library.fasta per species     (family → TE class mapping)
"""
from pathlib import Path
import re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter, defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────
TREE_FILE   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "ancestral_count/05_count/guide_clean.treefile")
WAGNER_FILE = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "ancestral_count/05_count_v2b/wagner_ancestral.txt")
LIB_DIR     = Path("/data/carlos/EDF-TE/TE_analysis/annotated")
OUTDIR      = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")

# ── Clade definitions ─────────────────────────────────────────────────────
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

# TE classes
TE_CLASSES = ["LTR", "DNA", "LINE", "SINE", "RC", "Unknown"]
TE_CLASS_COLORS = {
    "LTR": "#e41a1c", "DNA": "#377eb8", "LINE": "#4daf4a",
    "SINE": "#f1c232", "RC": "#984ea3", "Unknown": "#bdbdbd",
}

# ── Newick parser ──────────────────────────────────────────────────────────
class Node:
    __slots__ = ["name", "length", "children", "x", "y", "count_id"]
    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.x = 0.0
        self.y = 0.0
        self.count_id = ""


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
    out = [node]
    for c in node.children:
        out.extend(get_all_nodes(c))
    return out


# ── Assign Count IDs ──────────────────────────────────────────────────────
def assign_count_ids(root):
    """Count.jar numbers internal nodes in post-order, starting from n_tips."""
    leaves = get_leaves(root)
    n_tips = len(leaves)
    for leaf in leaves:
        leaf.count_id = leaf.name
    counter = [n_tips]
    def _postorder(node):
        if not node.children:
            return
        for c in node.children:
            _postorder(c)
        node.count_id = f"${counter[0]}_noname"
        counter[0] += 1
    _postorder(root)


# ── Cladogram layout ──────────────────────────────────────────────────────
def get_max_depth(node, depth=0):
    if not node.children:
        return depth
    return max(get_max_depth(c, depth + 1) for c in node.children)


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


# ── Draw tree ─────────────────────────────────────────────────────────────
def draw_tree(ax, node):
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


# ── Build family → class mapping ─────────────────────────────────────────
def build_family_class_map():
    """Parse all combined libraries to map family_name → TE class."""
    family_class = {}
    libs = sorted(LIB_DIR.glob("*_EarlGrey/*_Curated_Library/*_combined_library.fasta"))
    for lib in libs:
        species = lib.parent.name.replace("_Curated_Library", "")
        with open(lib) as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    if "#" in header:
                        name, cls_sf = header.split("#", 1)
                        broad = cls_sf.split("/")[0]
                        if broad not in TE_CLASSES:
                            broad = "Unknown"
                        family_class[(species, name.upper())] = broad
    return family_class


# ── Parse Wagner matrix ──────────────────────────────────────────────────
def parse_wagner_matrix(family_class_map):
    """Parse the full Wagner output matrix. Return:
       - node_composition: dict  node_id → Counter({class: count})
       for all nodes (tips + internal).
    """
    node_composition = defaultdict(Counter)
    header_nodes = None

    with open(WAGNER_FILE, encoding="latin-1") as f:
        for line in f:
            if not line.startswith("# FAMILY"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[1] == "name":
                # Header row
                header_nodes = parts[2:-4]  # skip FAMILY, name, and last 4 columns
                continue

            # Data row: # FAMILY <family_name> <0/1 values...> <gains> <losses> <exp> <red>
            family_name = parts[1]
            values = parts[2:-4]

            # Determine TE class for this family
            if "|" in family_name:
                sp, local = family_name.split("|", 1)
                cls = family_class_map.get((sp, local.upper()), "Unknown")
            else:
                cls = "Unknown"

            # Mark which nodes have this family present
            for i, v in enumerate(values):
                if v == "1":
                    node_composition[header_nodes[i]][cls] += 1

    return node_composition, header_nodes


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── Parse tree ────────────────────────────────────────────────────────
    print("[1/4] Parsing tree...")
    tree_str = TREE_FILE.read_text().strip()
    root = parse_newick(tree_str)
    root.x = 0.0
    assign_count_ids(root)

    leaves = get_leaves(root)
    n_tips = len(leaves)
    all_nodes = get_all_nodes(root)
    print(f"       {n_tips} tips, {len(all_nodes) - n_tips} internal nodes")

    # ── Build family→class mapping ────────────────────────────────────────
    print("[2/4] Building family→class mapping from combined libraries...")
    family_class_map = build_family_class_map()
    print(f"       {len(family_class_map)} entries")

    # ── Parse Wagner matrix ───────────────────────────────────────────────
    print("[3/4] Parsing Wagner matrix (45,840 families × 475 nodes)...")
    node_composition, header_nodes = parse_wagner_matrix(family_class_map)
    print(f"       Composition data for {len(node_composition)} nodes")

    # Map composition to tree nodes
    node_comp_by_id = {}
    for node in all_nodes:
        cid = node.count_id
        if cid in node_composition:
            node_comp_by_id[cid] = node_composition[cid]

    # ── Layout ────────────────────────────────────────────────────────────
    print("[4/4] Generating figure...")
    max_d = get_max_depth(root)
    assign_cladogram_positions(root, [0], depth=0, max_depth=max_d)

    # ── Figure (2 panels: proportional composition + absolute count) ─────
    fig_height = max(20, n_tips * 0.065)
    fig, axes = plt.subplots(1, 5, figsize=(24, fig_height),
                              gridspec_kw={"width_ratios": [0.3, 5, 2.5, 1.5, 0.5],
                                           "wspace": 0.02})
    ax_clade, ax_tree, ax_prop, ax_abs, ax_total = axes

    # ── Clade strip ───────────────────────────────────────────────────────
    ax_clade.set_xlim(0, 1)
    ax_clade.set_ylim(-0.5, n_tips - 0.5)
    ax_clade.invert_yaxis()
    ax_clade.axis("off")
    for leaf in leaves:
        clade = leaf.name[:2]
        color = CLADE_COLORS.get(clade, "#ddd")
        ax_clade.barh(leaf.y, 1, height=0.9, color=color, edgecolor="none")

    # ── Draw tree on ax_tree ──────────────────────────────────────────────
    draw_tree(ax_tree, root)
    ax_tree.set_xlim(-0.5, max_d + 0.5)
    ax_tree.set_ylim(-0.5, n_tips - 0.5)
    ax_tree.invert_yaxis()
    for sp in ax_tree.spines.values():
        sp.set_visible(False)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])

    # ── Proportional stacked bars (% composition) ─────────────────────────
    bar_height = 0.8
    internal_bar_height = 0.6

    ax_prop.set_ylim(-0.5, n_tips - 0.5)
    ax_prop.invert_yaxis()
    ax_prop.set_xlabel("TE class composition (%)", fontsize=9)
    ax_prop.spines["top"].set_visible(False)
    ax_prop.spines["right"].set_visible(False)
    ax_prop.spines["left"].set_visible(False)
    ax_prop.set_yticks([])
    ax_prop.tick_params(axis='x', labelsize=6)

    for node in all_nodes:
        comp = node_comp_by_id.get(node.count_id)
        if not comp:
            continue
        total = sum(comp.values())
        if total == 0:
            continue

        is_leaf = len(node.children) == 0
        h = bar_height if is_leaf else internal_bar_height
        alpha = 0.9 if is_leaf else 0.7
        edge = "none" if is_leaf else "#333"
        lw = 0 if is_leaf else 0.15

        left = 0.0
        for cls in TE_CLASSES:
            cnt = comp.get(cls, 0)
            if cnt > 0:
                pct = cnt / total * 100.0
                ax_prop.barh(node.y, pct, left=left, height=h,
                            color=TE_CLASS_COLORS[cls], edgecolor=edge,
                            linewidth=lw, alpha=alpha)
                left += pct

    ax_prop.set_xlim(0, 100)

    # ── Absolute count bars (log-scaled) ──────────────────────────────────
    ax_abs.set_ylim(-0.5, n_tips - 0.5)
    ax_abs.invert_yaxis()
    ax_abs.set_xlabel("Total families", fontsize=9)
    ax_abs.spines["top"].set_visible(False)
    ax_abs.spines["right"].set_visible(False)
    ax_abs.spines["left"].set_visible(False)
    ax_abs.set_yticks([])
    ax_abs.tick_params(axis='x', labelsize=6)
    ax_abs.set_xscale("symlog", linthresh=10)

    for node in all_nodes:
        comp = node_comp_by_id.get(node.count_id)
        if not comp:
            continue
        total = sum(comp.values())
        if total == 0:
            continue

        is_leaf = len(node.children) == 0
        h = bar_height if is_leaf else internal_bar_height
        clade = node.name[:2] if is_leaf else "#888"
        color = CLADE_COLORS.get(clade, "#888") if is_leaf else "#aaa"
        alpha = 0.7 if is_leaf else 0.5

        ax_abs.barh(node.y, total, height=h, color=color,
                   edgecolor="none", alpha=alpha)

    # ── Total families labels ─────────────────────────────────────────────
    ax_total.set_ylim(-0.5, n_tips - 0.5)
    ax_total.invert_yaxis()
    ax_total.axis("off")

    for node in all_nodes:
        if node.children:
            comp = node_comp_by_id.get(node.count_id)
            if comp:
                total = sum(comp.values())
                if total >= 20:
                    ax_total.text(0.5, node.y, str(total), fontsize=3,
                                 ha="center", va="center", color="#555")

    # ── Legends ───────────────────────────────────────────────────────────
    class_handles = [mpatches.Patch(color=TE_CLASS_COLORS[c], label=c)
                     for c in TE_CLASSES]
    clade_handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#ddd"),
                     label=CLADE_NAMES.get(c, c))
                     for c in CLADE_ORDER if c in CLADE_COLORS]

    leg1 = fig.legend(handles=class_handles, fontsize=6, loc="upper right",
                      bbox_to_anchor=(0.98, 0.99), title="TE class",
                      title_fontsize=7, framealpha=0.95)
    fig.legend(handles=clade_handles, fontsize=5, loc="upper left",
               bbox_to_anchor=(0.01, 0.99), ncol=2, framealpha=0.95,
               title="Clades", title_fontsize=6)

    # ── Annotate root ─────────────────────────────────────────────────────
    root_comp = node_comp_by_id.get(root.count_id, Counter())
    root_total = sum(root_comp.values())
    ax_tree.annotate(f"Root: {root_total} families",
                    xy=(root.x, root.y), xytext=(root.x - 1.5, root.y - 3),
                    fontsize=5, ha="center", va="center", color="#333",
                    style="italic",
                    arrowprops=dict(arrowstyle="->", lw=0.5, color="#333"))

    # ── Title & save ──────────────────────────────────────────────────────
    fig.suptitle("Ancestral TE class composition inferred by parsimony\n"
                 f"(AsymmetricWagner, gain=2.0 | {n_tips} species, "
                 f"45,840 families | Root: {root_total} families)",
                 fontsize=11, y=0.998)
    fig.tight_layout(rect=[0, 0, 1, 0.993])

    out_pdf = OUTDIR / "ancestral_te_composition.pdf"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "ancestral_te_composition.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[DONE] {out_pdf}")

    # ── Print some stats ──────────────────────────────────────────────────
    print(f"\nRoot composition ({root_total} families):")
    for cls in TE_CLASSES:
        print(f"  {cls}: {root_comp.get(cls, 0)}")

    # Some key internal nodes
    internal_nodes = [n for n in all_nodes if n.children]
    internal_nodes.sort(key=lambda n: -sum(node_comp_by_id.get(n.count_id, Counter()).values()))
    print(f"\nTop 10 internal nodes by total families:")
    for n in internal_nodes[:10]:
        comp = node_comp_by_id.get(n.count_id, Counter())
        total = sum(comp.values())
        print(f"  {n.count_id}: {total} families — "
              + ", ".join(f"{c}:{comp.get(c,0)}" for c in TE_CLASSES if comp.get(c,0)))


if __name__ == "__main__":
    main()
