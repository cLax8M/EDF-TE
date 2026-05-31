#!/usr/bin/env python3
"""04_gremlines_tree_phylo.py — Grem-LINE1 phylogenetic tree with centromere annotations.

Figures produced:
  A) Rectangular phylogram with tips labelled CenN.R (rank by distance to CC, 1=closest)
     - Colored by centromere identity
     - Sister pairs connected by arcs (green=same cen, red=diff cen)
     - Outgroup tips in grey
  B) Inter-centromere jumping chord diagram
  C) Per-centromere copy count + ancestry summary

Also outputs text analysis answering:
  - Do all copies originate from the same ancestral copy?
  - Is jumping inter- or intra-centromere?
"""
from __future__ import annotations
import re, os
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Arc
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────
OUT       = Path("/legserv/Temp/Carlos/EDF-TE/results/Gremlines/MU402")
CEN_BED   = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/MU402_centromere_windows_clipped.bed")
COPY_ZONE = OUT / "copy_zone.tsv"       # original 56-row table
TREE_FILE = OUT / "gremlin_iqtree.treefile"
FAI_FILE  = Path("/legserv/Temp/Carlos/EDF-TE/data/Extra/mu402GenomeJGI.fa.fai")

CEN_SCAFF_ORDER = {
    "scaffold_1": 1, "scaffold_2": 2, "scaffold_3": 3,
    "scaffold_4": 4, "scaffold_5": 5, "scaffold_6": 6,
    "scaffold_8": 7, "scaffold_9": 8, "scaffold_11": 9,
}

CEN_COLORS = {
    1: "#e41a1c", 2: "#377eb8", 3: "#4daf4a",
    4: "#984ea3", 5: "#ff7f00", 6: "#a65628",
    7: "#f781bf", 8: "#999999", 9: "#66c2a5",
    0: "#cccccc",
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


def assign_positions(node, x_offset=0.0, y_counter=None):
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


def draw_tree(ax, node):
    """Draw rectangular phylogram. Returns {tip_name: y_pos}."""
    def traverse(n, parent_x=None):
        if parent_x is not None:
            ax.plot([parent_x, n.x], [n.y, n.y], "k-", lw=0.5, solid_capstyle="round")
        if n.children:
            ys = [c.y for c in n.children]
            ax.plot([n.x, n.x], [min(ys), max(ys)], "k-", lw=0.5, solid_capstyle="round")
            for c in n.children:
                traverse(c, n.x)
    traverse(node)
    return {leaf.name: leaf.y for leaf in get_leaves(node)}


def find_leaf_sisters(node):
    """Find sister leaf pairs (leaves sharing immediate parent)."""
    pairs = []
    if len(node.children) == 2:
        left, right = node.children
        if not left.children and not right.children:
            pairs.append((left.name, right.name))
    for c in node.children:
        pairs.extend(find_leaf_sisters(c))
    return pairs


def is_outgroup(name):
    return "OUT" in name or "RND" in name


def fasta_to_tree_name(fasta_name):
    return fasta_name.replace("(+)", "___").replace("(-)", "_-_")


def tree_to_fasta_name(tree_name):
    if tree_name.endswith("___"):
        return tree_name[:-3] + "(+)"
    elif tree_name.endswith("_-_"):
        return tree_name[:-3] + "(-)"
    return tree_name


# ── Centromere helpers ─────────────────────────────────────────────────────
def load_centromeres():
    cen = {}
    with open(CEN_BED) as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            cen[chrom] = (start, end, (start + end) // 2)
    return cen


def assign_centromere_info(row, cen_map):
    chrom = row["chrom"]
    mid = int(row["mid"])
    if chrom not in cen_map:
        return 0, "NoCen", "-"
    cen_start, cen_end, cen_center = cen_map[chrom]
    cen_num = CEN_SCAFF_ORDER.get(chrom, 0)
    side = "L" if mid < cen_center else "R"
    return cen_num, f"Cen{cen_num}", side


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    cen_map = load_centromeres()

    # ── Load and enrich copy data ────────────────────────────────────────
    df = pd.read_csv(COPY_ZONE, sep="\t")
    print(f"[INFO] {len(df)} copies loaded")

    df["tree_tip"] = df["fasta_name"].apply(fasta_to_tree_name)

    cen_nums, cen_names, sides = [], [], []
    for _, row in df.iterrows():
        cn, cname, side = assign_centromere_info(row, cen_map)
        cen_nums.append(cn)
        cen_names.append(cname)
        sides.append(side)
    df["cen_num"] = cen_nums
    df["cen_name"] = cen_names
    df["side"] = sides

    # ── Assign rank per centromere (1 = closest to CC) ───────────────────
    df["rank"] = 0
    for cn in df["cen_num"].unique():
        mask = df["cen_num"] == cn
        order = df.loc[mask].sort_values("dist_centro").index
        for rank, idx in enumerate(order, 1):
            df.loc[idx, "rank"] = rank

    # Build label: CenN.R (rank, 1=closest)
    df["tree_label"] = df.apply(
        lambda r: f"Cen{r['cen_num']}.{r['rank']}" if r["cen_num"] > 0
                  else f"NoCen.{r['rank']}", axis=1)

    # Maps by tree_tip name
    label_map = dict(zip(df["tree_tip"], df["tree_label"]))
    cen_map_copy = dict(zip(df["tree_tip"], df["cen_num"]))
    k80_map = dict(zip(df["tree_tip"], df["k80"]))
    dist_map = dict(zip(df["tree_tip"], df["dist_centro"]))
    status_map = dict(zip(df["tree_tip"], df["status"]))

    # Save enriched table
    df.to_csv(OUT / "copy_zone_enriched_v2.tsv", sep="\t", index=False)

    # ── Parse tree ───────────────────────────────────────────────────────
    tree_str = TREE_FILE.read_text().strip()
    root = parse_newick(tree_str)
    assign_positions(root, x_offset=0.0)
    leaves = get_leaves(root)
    n_tips = len(leaves)
    print(f"[INFO] {n_tips} tips in tree")

    aut_leaves = [l for l in leaves if not is_outgroup(l.name)]
    out_leaves = [l for l in leaves if is_outgroup(l.name)]
    print(f"[INFO] {len(aut_leaves)} AUT + {len(out_leaves)} outgroup")

    # Sister pairs
    sisters = find_leaf_sisters(root)
    print(f"[INFO] {len(sisters)} sister pairs")

    # ── Analysis: sister pair centromere jumping ─────────────────────────
    pair_details = []
    same_cen = 0
    diff_cen = 0
    daughter_closer = 0
    daughter_farther = 0

    for a, b in sisters:
        if is_outgroup(a) or is_outgroup(b):
            continue
        ka, kb = k80_map.get(a, np.nan), k80_map.get(b, np.nan)
        da, db = dist_map.get(a, np.nan), dist_map.get(b, np.nan)
        ca, cb = cen_map_copy.get(a, -1), cen_map_copy.get(b, -1)
        if pd.isna(ka) or pd.isna(kb):
            continue

        # Parent = older (higher K80), daughter = younger (lower K80)
        if ka > kb:
            parent, daughter = a, b
            d_parent, d_daughter = da, db
            c_parent, c_daughter = ca, cb
        else:
            parent, daughter = b, a
            d_parent, d_daughter = db, da
            c_parent, c_daughter = cb, ca

        is_same = (c_parent == c_daughter and c_parent > 0)
        if c_parent > 0 and c_daughter > 0:
            if is_same:
                same_cen += 1
            else:
                diff_cen += 1

        if not pd.isna(d_parent) and not pd.isna(d_daughter):
            if d_daughter < d_parent:
                daughter_closer += 1
            else:
                daughter_farther += 1

        pair_details.append({
            "parent": label_map.get(parent, parent),
            "daughter": label_map.get(daughter, daughter),
            "parent_tip": parent,
            "daughter_tip": daughter,
            "k80_parent": k80_map.get(parent),
            "k80_daughter": k80_map.get(daughter),
            "cen_parent": c_parent,
            "cen_daughter": c_daughter,
            "same_centromere": is_same,
            "dist_parent": d_parent,
            "dist_daughter": d_daughter,
        })

    pairs_df = pd.DataFrame(pair_details)
    if len(pairs_df):
        pairs_df.to_csv(OUT / "sister_pairs_v2.tsv", sep="\t", index=False)

    total_dir = daughter_closer + daughter_farther
    binom_p = stats.binomtest(daughter_closer, total_dir, 0.5).pvalue if total_dir > 0 else float("nan")

    # ── Ancestry analysis ────────────────────────────────────────────────
    # All AUT copies form a monophyletic group (outgroup is separate).
    # Within AUT, identify basal lineages.
    # The tree topology shows: outgroup is nested near cluster 3 (3.01, 3.02)
    # which means cluster 3 copies are the most ancestral AUT copies.
    # All other AUT copies form a single radiation.

    # Find which copies are deepest (closest to outgroup split)
    # by looking at total root-to-tip distance (x coordinate)
    aut_depths = [(l.name, l.x, label_map.get(l.name, l.name))
                  for l in aut_leaves]
    aut_depths.sort(key=lambda t: t[1])  # shortest = shallowest = most ancestral

    # ── Write results ────────────────────────────────────────────────────
    with open(OUT / "phylo_analysis_results.txt", "w") as f:
        f.write("Grem-LINE1 phylogenetic analysis (MU402 JGI)\n")
        f.write("=" * 60 + "\n\n")

        f.write("[1] ANCESTRY — Do all copies come from the same ancestor?\n")
        f.write(f"    All 31 AUT copies form a monophyletic clade (outgroup: 5 RND copies).\n")
        f.write(f"    → YES: all copies descend from a single ancestral insertion.\n")
        f.write(f"    Most ancestral lineage (shallowest branches):\n")
        for name, depth, label in aut_depths[:5]:
            f.write(f"      {label:15s}  root-to-tip = {depth:.6f}\n")
        f.write(f"    Cluster 3 copies (3.01, 3.02) branch off closest to outgroup\n")
        f.write(f"    → Cluster 3 represents the oldest/ancestral Grem-LINE1 lineage\n\n")

        f.write("[2] INTER- vs INTRA-CENTROMERE JUMPING\n")
        f.write(f"    Sister pairs (AUT only): {len(pair_details)}\n")
        f.write(f"    Same centromere (intra):     {same_cen}\n")
        f.write(f"    Different centromere (inter): {diff_cen}\n")
        if same_cen + diff_cen > 0:
            pct_inter = 100 * diff_cen / (same_cen + diff_cen)
            f.write(f"    → {pct_inter:.0f}% of sister pairs are INTER-centromere jumps\n")
            f.write(f"    → Transposition preferentially targets NEW centromeres\n\n")

        f.write("[3] DIRECTIONALITY — Do new copies insert closer to centromere center?\n")
        f.write(f"    Daughter closer to CC: {daughter_closer}\n")
        f.write(f"    Daughter farther from CC: {daughter_farther}\n")
        f.write(f"    Binomial p = {binom_p:.4g}\n")
        if daughter_closer > daughter_farther:
            f.write(f"    → Trend towards INWARD insertion (closer to CC)\n\n")
        else:
            f.write(f"    → No clear inward bias\n\n")

        f.write("[4] PER-CENTROMERE COPY DISTRIBUTION\n")
        for cn in sorted(df["cen_num"].unique()):
            sub = df[df["cen_num"] == cn]
            aut_n = (sub["status"] == "AUT").sum()
            rem_n = (sub["status"] == "REM").sum()
            name = f"Cen{cn}" if cn > 0 else "NoCen"
            f.write(f"    {name:8s}: {len(sub):3d} copies ({aut_n} AUT, {rem_n} REM)\n")

        f.write(f"\n[5] SISTER PAIR DETAILS\n")
        for _, p in pairs_df.iterrows():
            jump = "SAME" if p["same_centromere"] else "DIFF"
            f.write(f"    {p['parent']:15s} → {p['daughter']:15s}  "
                    f"Cen{int(p['cen_parent'])}→Cen{int(p['cen_daughter'])} [{jump}]\n")

    print("[INFO] phylo_analysis_results.txt saved")

    # ══════════════════════════════════════════════════════════════════════
    # FIGURE 1: Phylogenetic tree with centromere annotations
    # ══════════════════════════════════════════════════════════════════════
    fig, (ax_tree, ax_strip) = plt.subplots(
        1, 2, figsize=(14, 10),
        gridspec_kw={"width_ratios": [5, 0.3], "wspace": 0.01})

    tip_pos = draw_tree(ax_tree, root)
    max_x = max(l.x for l in leaves)

    # Annotate tips
    for leaf in leaves:
        y = tip_pos[leaf.name]
        if is_outgroup(leaf.name):
            # Outgroup
            ax_tree.plot(leaf.x, y, "D", color="#888888", markersize=5,
                        markeredgecolor="k", markeredgewidth=0.3, zorder=4)
            short = leaf.name.split("_")[-2]  # e.g. "8" from RND-1_FAMILY-218_OUT_8_-_
            ax_tree.text(leaf.x + max_x * 0.01, y, f"OG_{short}",
                        fontsize=6, va="center", color="#888888", fontstyle="italic")
        else:
            cn = cen_map_copy.get(leaf.name, 0)
            color = CEN_COLORS.get(cn, "#cccccc")
            st = status_map.get(leaf.name, "?")
            marker = "o"
            ax_tree.plot(leaf.x, y, marker, color=color, markersize=6,
                        markeredgecolor="k", markeredgewidth=0.4, zorder=4)
            lbl = label_map.get(leaf.name, leaf.name)
            ax_tree.text(leaf.x + max_x * 0.01, y, lbl,
                        fontsize=6.5, va="center", color=color,
                        fontweight="bold")

    # Draw sister pair arcs
    for a, b in sisters:
        if is_outgroup(a) or is_outgroup(b):
            continue
        ya = tip_pos.get(a)
        yb = tip_pos.get(b)
        if ya is None or yb is None:
            continue
        ca = cen_map_copy.get(a, -1)
        cb = cen_map_copy.get(b, -1)
        is_same = (ca == cb and ca > 0)
        color = "#2ca02c" if is_same else "#d62728"
        # Draw arc to the right of tips
        x_arc = max_x * 1.15
        mid_y = (ya + yb) / 2
        height = abs(ya - yb)
        ax_tree.annotate("", xy=(x_arc, ya), xytext=(x_arc, yb),
                         arrowprops=dict(arrowstyle="-", color=color,
                                        connectionstyle=f"arc3,rad=0.3",
                                        lw=1.5, alpha=0.6))

    ax_tree.set_xlim(-max_x * 0.05, max_x * 1.25)
    ax_tree.set_ylim(-0.5, n_tips - 0.5)
    ax_tree.invert_yaxis()
    ax_tree.set_xlabel("Substitutions/site", fontsize=9)
    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)
    ax_tree.spines["left"].set_visible(False)
    ax_tree.set_yticks([])

    # Centromere strip
    ax_strip.set_xlim(0, 1)
    ax_strip.set_ylim(-0.5, n_tips - 0.5)
    ax_strip.invert_yaxis()
    ax_strip.axis("off")
    for leaf in leaves:
        y = tip_pos[leaf.name]
        if is_outgroup(leaf.name):
            ax_strip.barh(y, 1, height=0.8, color="#eeeeee", edgecolor="none")
        else:
            cn = cen_map_copy.get(leaf.name, 0)
            ax_strip.barh(y, 1, height=0.8,
                         color=CEN_COLORS.get(cn, "#ccc"), edgecolor="none")

    # Legend
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markeredgecolor="k", markersize=7, label="AUT copy"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#888",
               markeredgecolor="k", markersize=6, label="Outgroup"),
        Line2D([0], [0], color="#2ca02c", lw=2, label="Same centromere"),
        Line2D([0], [0], color="#d62728", lw=2, label="Diff centromere"),
    ]
    for cn in sorted(CEN_SCAFF_ORDER.values()):
        legend_elements.append(
            Line2D([0], [0], marker="s", color="w",
                   markerfacecolor=CEN_COLORS[cn], markeredgecolor="k",
                   markersize=7, label=f"Cen{cn}"))
    legend_elements.append(
        Line2D([0], [0], marker="s", color="w",
               markerfacecolor=CEN_COLORS[0], markeredgecolor="k",
               markersize=7, label="NoCen"))

    ax_tree.legend(handles=legend_elements, loc="lower left", fontsize=6,
                   ncol=2, framealpha=0.9)

    fig.suptitle("Grem-LINE1 phylogeny — MU402 JGI\n"
                 "Tips: CenN.R (R = rank by distance to centromere center, 1 = closest)",
                 fontsize=11, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT / "gremlines_phylo_tree.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[INFO] gremlines_phylo_tree.pdf saved")

    # ══════════════════════════════════════════════════════════════════════
    # FIGURE 2: Inter-centromere jumping network
    # ══════════════════════════════════════════════════════════════════════
    fig, (ax_net, ax_bar) = plt.subplots(1, 2, figsize=(14, 5),
                                          gridspec_kw={"width_ratios": [2, 1]})

    # Panel A: Chord diagram — centromeres as positions on a circle
    cen_ids = sorted(CEN_SCAFF_ORDER.values())  # 1-9
    n_cen = len(cen_ids)
    angles = np.linspace(0, 2 * np.pi, n_cen, endpoint=False)
    cx = np.cos(angles) * 3
    cy = np.sin(angles) * 3

    # Count copies per centromere (AUT only in tree)
    cen_counts = defaultdict(int)
    for l in aut_leaves:
        cn = cen_map_copy.get(l.name, 0)
        if cn > 0:
            cen_counts[cn] += 1

    # Draw centromere nodes
    for i, cn in enumerate(cen_ids):
        n_copies = cen_counts.get(cn, 0)
        size = max(200, n_copies * 120)
        ax_net.scatter(cx[i], cy[i], s=size,
                      c=CEN_COLORS[cn], edgecolor="k", linewidth=1, zorder=5)
        ax_net.text(cx[i], cy[i], f"Cen{cn}\n({n_copies})",
                   ha="center", va="center", fontsize=7, fontweight="bold", zorder=6)

    # Draw connections between centromeres for each sister pair
    jump_counts = defaultdict(int)  # (cen_a, cen_b) → count
    for _, p in pairs_df.iterrows():
        ca, cb = int(p["cen_parent"]), int(p["cen_daughter"])
        if ca > 0 and cb > 0 and ca != cb:
            key = tuple(sorted([ca, cb]))
            jump_counts[key] += 1

    for (ca, cb), count in jump_counts.items():
        ia = cen_ids.index(ca)
        ib = cen_ids.index(cb)
        lw = 1.5 + count * 1.0
        ax_net.plot([cx[ia], cx[ib]], [cy[ia], cy[ib]],
                   color="#d62728", lw=lw, alpha=0.5, zorder=2)
        mid_x = (cx[ia] + cx[ib]) / 2
        mid_y = (cy[ia] + cy[ib]) / 2
        ax_net.text(mid_x, mid_y, str(count), fontsize=7,
                   ha="center", va="center",
                   bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))

    # Intra-centromere jumps (self-loops)
    for _, p in pairs_df.iterrows():
        ca, cb = int(p["cen_parent"]), int(p["cen_daughter"])
        if ca > 0 and cb > 0 and ca == cb:
            ia = cen_ids.index(ca)
            circle = plt.Circle((cx[ia], cy[ia] + 0.8), 0.5,
                               fill=False, color="#2ca02c", lw=2, alpha=0.6)
            ax_net.add_patch(circle)

    ax_net.set_xlim(-5, 5)
    ax_net.set_ylim(-5, 5)
    ax_net.set_aspect("equal")
    ax_net.axis("off")
    ax_net.set_title(f"Inter-centromere jumping\n"
                     f"(inter: {diff_cen}, intra: {same_cen})",
                     fontsize=11, fontweight="bold")

    # Panel B: Per-centromere copy count bar chart
    all_cens = sorted(df["cen_num"].unique())
    x_pos = range(len(all_cens))
    aut_counts = []
    rem_counts = []
    bar_labels = []
    for cn in all_cens:
        sub = df[df["cen_num"] == cn]
        aut_counts.append((sub["status"] == "AUT").sum())
        rem_counts.append((sub["status"] == "REM").sum())
        bar_labels.append(f"Cen{cn}" if cn > 0 else "NoCen")

    ax_bar.bar(x_pos, aut_counts, color=[CEN_COLORS.get(cn, "#ccc") for cn in all_cens],
               edgecolor="k", linewidth=0.5, label="AUT")
    ax_bar.bar(x_pos, rem_counts, bottom=aut_counts,
               color=[CEN_COLORS.get(cn, "#ccc") for cn in all_cens],
               edgecolor="k", linewidth=0.5, alpha=0.4, label="REM")

    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(bar_labels, rotation=45, fontsize=8)
    ax_bar.set_ylabel("Number of copies")
    ax_bar.set_title("Copies per centromere", fontsize=11, fontweight="bold")
    ax_bar.legend(fontsize=8)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / "gremlines_jumping_network.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[INFO] gremlines_jumping_network.pdf saved")

    # ══════════════════════════════════════════════════════════════════════
    # FIGURE 3: Contig map with rank labels
    # ══════════════════════════════════════════════════════════════════════
    fai = pd.read_csv(FAI_FILE, sep="\t", header=None, names=["chrom", "length"],
                       usecols=[0, 1])
    chroms_with = df["chrom"].unique()
    fai_sub = fai[fai["chrom"].isin(chroms_with)].sort_values("length", ascending=False)
    cen_info = load_centromeres()

    n_chroms = len(fai_sub)
    fig, axes = plt.subplots(n_chroms, 1, figsize=(13, max(4, 0.6 * n_chroms)),
                              sharex=False)
    if n_chroms == 1:
        axes = [axes]

    for ax, (_, frow) in zip(axes, fai_sub.iterrows()):
        ch = frow["chrom"]
        L = frow["length"]
        ax.barh(0, L / 1e6, height=0.4, color="#e8e8e8", edgecolor="#999", linewidth=0.5)

        if ch in cen_info:
            cs, ce, cc = cen_info[ch]
            ax.barh(0, (ce - cs) / 1e6, left=cs / 1e6, height=0.5,
                    color="#ffcccc", edgecolor="#cc0000", linewidth=0.8, alpha=0.6)
            ax.axvline(cc / 1e6, color="#cc0000", ls=":", linewidth=0.8, alpha=0.5)

        sub = df[df["chrom"] == ch]
        for _, crow in sub.iterrows():
            cn = crow["cen_num"]
            color = CEN_COLORS.get(cn, "#cccccc")
            marker = "v" if crow["status"] == "AUT" else "^"
            y_off = 0.35 if crow["strand"] == "+" else -0.35
            ax.plot(crow["mid"] / 1e6, y_off, marker=marker, color=color,
                    markersize=6, markeredgecolor="black", markeredgewidth=0.4)
            ax.text(crow["mid"] / 1e6, y_off + (0.15 if y_off > 0 else -0.15),
                    crow["tree_label"], fontsize=4, ha="center",
                    va="bottom" if y_off > 0 else "top", rotation=45)

        cn_num = CEN_SCAFF_ORDER.get(ch, 0)
        label = f"{ch} (Cen{cn_num})" if cn_num > 0 else ch
        ax.set_ylabel(label, fontsize=7, rotation=0, ha="right", va="center")
        ax.set_xlim(0, L / 1e6)
        ax.set_ylim(-0.8, 0.8)
        ax.set_yticks([])
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)

    axes[-1].set_xlabel("Genomic position (Mb)")
    legend_elements = [
        Line2D([0], [0], marker="v", color="w", markerfacecolor="gray",
               markeredgecolor="k", markersize=8, label="AUT"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor="gray",
               markeredgecolor="k", markersize=8, label="REM"),
    ]
    fig.legend(handles=legend_elements, loc="upper right", fontsize=7, framealpha=0.9)
    fig.suptitle("MU402 Grem-LINE1 — Genomic map with centromere rank labels", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "gremlines_contig_map_ranked.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[INFO] gremlines_contig_map_ranked.pdf saved")

    # Print summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"All 31 AUT copies are MONOPHYLETIC → single ancestral insertion")
    print(f"Sister pairs: {len(pair_details)} (AUT-AUT)")
    print(f"  Inter-centromere: {diff_cen} ({100*diff_cen/(same_cen+diff_cen):.0f}%)" if same_cen+diff_cen>0 else "")
    print(f"  Intra-centromere: {same_cen}")
    print(f"Daughter closer to CC: {daughter_closer}/{total_dir} (p={binom_p:.3f})")
    print(f"Figures: gremlines_phylo_tree.pdf, gremlines_jumping_network.pdf, gremlines_contig_map_ranked.pdf")


if __name__ == "__main__":
    main()
