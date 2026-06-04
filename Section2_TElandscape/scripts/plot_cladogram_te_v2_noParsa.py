#!/usr/bin/env python3
"""plot_cladogram_te_v2.py — Cladogram + TE bars + Shannon + Recency Index.

Improvements over v1:
  - Added Recency Index (RI) track (fraction TE bp with div < 0.1)
  - Narrower Shannon and RI panels
  - Ensured ALL present clades appear in legend (including NU)
  - Output: cladogram_te_types_shannon_ri.pdf
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
TREE_FILE   = Path("/legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/final_LG_G_rotated.nwk")
CLASS_MAT   = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "figures/class_level_matrix.tsv")
SPECIES_DAT = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "figure_species_data.tsv")
SHANNON_TSV = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "B01_shannon_diversity.tsv")
KIMURA_TSV  = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "C01_kimura_species.tsv")
OUTDIR      = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")

# ── Clade config (canonical palette) ─────────────────────────────
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from clade_palette import CLADE_COLORS, CLADE_NAMES  # noqa: E402

# Species whose prefix is CH_ but belong to Neocallimastigomycota
NE_SPECIES = {
    "CH_Anasp1", "CH_Neocon1", "CH_Neolan1",
    "CH_Orpsp1", "CH_Pecora1", "CH_Caecom1", "CH_Pirfi3",
}

# Name mapping: tree tip name → data file name (when different)
NAME_MAP = {
    "NU_Parvatl_c": "NU_Parvatl",
    "AS_AspfumUD1": "AS_AspfumUD",
    "HO_Magisepte_c": "HO_Magisepte",
}

TE_ORDER = ["LTR", "DNA", "LINE", "SINE", "Rolling Circle", "Penelope", "Unclassified"]

# ── Manual two-letter code overrides ─────────────────────────────────────
MANUAL_CODE_OVERRIDES = {
    "CR_Rozal1":    "Ra",  # Rozella allomyces [algorithm gives Rl]
    "HO_Hydravulg": "Hv",  # Hydra vulgaris    [algorithm gives Hr]
}


def _raw_sp_code(name):
    import re as _re
    if name in MANUAL_CODE_OVERRIDES:
        return MANUAL_CODE_OVERRIDES[name]
    base = name.split("_", 1)[1] if "_" in name else name
    base = _re.sub(r'[\d_].*', '', base)
    if len(base) < 2:
        return (base[0].upper() + "?") if base else "??"
    if len(base) < 4:
        return base[0].upper() + base[-1].lower()
    c1 = base[0].upper()
    vowels = set('aeiou')
    c2 = base[3].lower()
    if c2 in vowels and len(base) > 4:
        c2 = base[4].lower()
    return f"{c1}{c2}"


def compute_tip_codes(ordered_names):
    import re as _re
    used = set()
    final = {}
    for name in ordered_names:
        raw = _raw_sp_code(name)
        if raw not in used:
            final[name] = raw
            used.add(raw)
        else:
            base = name.split("_", 1)[1] if "_" in name else name
            base = _re.sub(r'[\d_].*', '', base)
            c1 = base[0].upper() if base else "?"
            assigned = False
            for i in range(1, len(base)):
                candidate = f"{c1}{base[i].lower()}"
                if candidate not in used:
                    final[name] = candidate
                    used.add(candidate)
                    assigned = True
                    break
            if not assigned:
                full = (name.split("_", 1)[1] if "_" in name else name).lower()
                for i in range(1, len(full)):
                    candidate = f"{c1}{full[i]}"
                    if candidate not in used:
                        final[name] = candidate
                        used.add(candidate)
                        assigned = True
                        break
            if not assigned:
                final[name] = raw
    return final


def get_clade(species_name):
    """Return clade code for a species, handling NE overrides."""
    if species_name in NE_SPECIES:
        return "NE"
    return species_name[:2]


TE_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#f1c232",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#66c2a5",
    "Unclassified":   "#bdbdbd",
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


def drop_tip(node, tip_name):
    """Remove a tip by name; prune resulting single-child internal nodes."""
    if not node.children:
        return None if node.name == tip_name else node
    node.children = [c for c in [drop_tip(c, tip_name) for c in node.children] if c is not None]
    if len(node.children) == 0:
        return None
    if len(node.children) == 1 and node.name == "":
        return node.children[0]
    return node


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
    if not node.children:
        node.x = max_depth
        node.y = y_counter[0]
        y_counter[0] += 1
    else:
        node.x = depth
        for c in node.children:
            assign_cladogram_positions(c, y_counter, depth + 1, max_depth)
        node.y = np.mean([c.y for c in node.children])


def assign_phylogram_positions(node, y_counter, x=0.0):
    """Assign x = cumulative branch length (root → tip in subs/site).
    y is an in-order tip index (children get integers; internal nodes the mean).
    """
    node.x = x
    if not node.children:
        node.y = y_counter[0]
        y_counter[0] += 1
    else:
        for c in node.children:
            assign_phylogram_positions(c, y_counter, x + max(c.length, 0.0))
        node.y = float(np.mean([c.y for c in node.children]))


def get_max_x(node):
    if not node.children:
        return node.x
    return max(get_max_x(c) for c in node.children)


def draw_cladogram(ax, node):
    def traverse(n, parent_x=None):
        if parent_x is not None:
            ax.plot([parent_x, n.x], [n.y, n.y], "k-", lw=0.3, solid_capstyle="round")
        if n.children:
            ys = [c.y for c in n.children]
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
    root = drop_tip(root, "CR_Parsa1")
    print("[INFO] Dropped CR_Parsa1 (Paramicrosporidium) from tree")

    # ── Rotate nodes to match desired clade display order ────────────────────
    # bottom-to-top: index 0 = bottom (first leaf, lowest y); last = top
    _DESIRED_CLADE_ORDER = ["HO","NU","AP","CR","MI","BL","NE","CH","OL","BS","ZO","EN","KX","MO","GL","MU","AS","BA"]
    _clade_rank = {cl: i for i, cl in enumerate(_DESIRED_CLADE_ORDER)}
    def _rotate_to_order(node):
        if not node.children:
            return
        for c in node.children:
            _rotate_to_order(c)
        def _min_idx(n):
            lv = get_leaves(n)
            idxs = [_clade_rank[l.name[:2]] for l in lv if l.name[:2] in _clade_rank]
            return min(idxs) if idxs else len(_DESIRED_CLADE_ORDER)
        node.children.sort(key=_min_idx)
    _rotate_to_order(root)
    print("[INFO] Tree rotated to match desired clade order")

    # Cladogram: all tips at the same x position (equal branch lengths)
    max_d = get_max_depth(root)
    assign_cladogram_positions(root, [0], depth=0, max_depth=max_d)
    leaves = get_leaves(root)
    n_tips = len(leaves)
    print(f"[INFO] Tree has {n_tips} tips, max depth {max_d} levels")

    # Helper: look up data using tree tip name or mapped name
    def data_name(tree_tip):
        return NAME_MAP.get(tree_tip, tree_tip)

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

    # Load Recency Index
    df_ki = pd.read_csv(KIMURA_TSV, sep="\t")
    ri_map = dict(zip(df_ki["species"], df_ki["recency_index"]))
    print(f"[INFO] Recency index data for {len(ri_map)} species")

    # Check NU coverage
    nu_tips = [l.name for l in get_leaves(root) if l.name.startswith("NU_")]
    for nt in nu_tips:
        dn = data_name(nt)
        has_class = dn in class_map
        has_te = dn in te_total_map or nt in te_total_map
        has_sh = dn in shannon_map or nt in shannon_map
        has_ri = dn in ri_map or nt in ri_map
        print(f"[INFO] {nt} → {dn}: class={has_class} te_total={has_te} shannon={has_sh} ri={has_ri}")

    # ── Figure layout: 5 panels ──────────────────────────────────────────
    # clade | tree | TE bars | Shannon | Recency
    fig_height = max(14, n_tips * 0.065)
    fig, axes = plt.subplots(1, 5, figsize=(20, fig_height),
                              gridspec_kw={
                                  "width_ratios": [0.25, 3, 3.5, 0.8, 0.8],
                                  "wspace": 0.02
                              })
    ax_clade, ax_tree, ax_bars, ax_shan, ax_ri = axes

    # ── Draw cladogram ───────────────────────────────────────────────────
    tip_pos = draw_cladogram(ax_tree, root)
    ax_tree.set_xlim(-max_d * 0.02, max_d * 1.12)
    ax_tree.set_ylim(-0.5, n_tips - 0.5)
    for sp in ax_tree.spines.values():
        sp.set_visible(False)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])

    # ── 2-letter italic species code at each tip (e.g. Sc = Saccharomyces cerevisiae) ──
    # Pre-compute all codes with deduplication (tree order: bottom→top)
    tip_codes = compute_tip_codes([leaf.name for leaf in leaves])

    for leaf in leaves:
        y = tip_pos[leaf.name]
        code = tip_codes[leaf.name]
        ax_tree.text(max_d * 1.005, y, code, fontsize=4.5,
                     ha="left", va="center", style="italic",
                     fontweight="bold", color="#222")

    # ── Scale bar removed (cladogram has no branch length units) ────────────────
    ax_tree.set_ylim(-0.5, n_tips - 0.5)

    # ── Clade color strip ────────────────────────────────────────────────
    ax_clade.set_xlim(0, 1)
    ax_clade.set_ylim(-0.5, n_tips - 0.5)
    ax_clade.axis("off")
    for leaf in leaves:
        y = tip_pos[leaf.name]
        clade = get_clade(leaf.name)
        color = CLADE_COLORS.get(clade, "#ddd")
        ax_clade.barh(y, 1, height=0.9, color=color, edgecolor="none")

    # ── Stacked TE bars (genome %) ───────────────────────────────────────
    max_te = 0
    bar_data = {}
    for leaf in leaves:
        sp = leaf.name
        dn = data_name(sp)
        te_total = te_total_map.get(sp, te_total_map.get(dn, 0))
        cm = class_map.get(sp, class_map.get(dn, None))
        if cm is not None and te_total > 0:
            vals = {}
            for t in TE_ORDER:
                vals[t] = cm.get(t, 0) * te_total
            bar_data[sp] = vals
            total = sum(vals.values())
            # Exclude HO/NU outgroups from axis-scale calculation
            # (animal genomes have 40-80% TE and would dwarf fungal bars)
            if get_clade(sp) not in ("HO", "NU") and total > max_te:
                max_te = total

    ax_bars.set_xlim(0, min(max_te * 1.1, 60))
    ax_bars.set_ylim(-0.5, n_tips - 0.5)
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

    # ── Horizontal zebra stripes for species tracking ────────────────────
    for i in range(n_tips):
        if i % 2 == 0:
            for ax in [ax_bars, ax_shan, ax_ri]:
                ax.axhspan(i - 0.5, i + 0.5, color="#f5f5f5", zorder=0)

    # ── Shannon bubbles (narrower panel) ─────────────────────────────────
    ax_shan.set_ylim(-0.5, n_tips - 0.5)
    ax_shan.set_xlabel("Shannon H'", fontsize=7)
    ax_shan.spines["top"].set_visible(False)
    ax_shan.spines["right"].set_visible(False)
    ax_shan.spines["left"].set_visible(False)
    ax_shan.set_yticks([])
    ax_shan.tick_params(axis='x', labelsize=6)

    h_vals = [shannon_map.get(l.name, shannon_map.get(data_name(l.name), np.nan)) for l in leaves]
    h_vals_clean = [h for h in h_vals if not pd.isna(h)]
    max_h = max(h_vals_clean) if h_vals_clean else 3

    for leaf in leaves:
        y = tip_pos[leaf.name]
        sp = leaf.name
        dn = data_name(sp)
        h = shannon_map.get(sp, shannon_map.get(dn, np.nan))
        if pd.isna(h):
            continue
        clade = get_clade(sp)
        ax_shan.scatter(h, y, s=max(6, h * 18), alpha=0.7,
                       c=CLADE_COLORS.get(clade, "#aaa"),
                       edgecolor="k", linewidth=0.15, zorder=3)
    ax_shan.set_xlim(0, max_h * 1.1)
    ax_shan.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax_shan.grid(axis="x", color="#ddd", linewidth=0.4, zorder=1)

    # ── Recency Index bubbles ────────────────────────────────────────────
    ax_ri.set_ylim(-0.5, n_tips - 0.5)
    ax_ri.set_xlabel("Recency (RI)", fontsize=7)
    ax_ri.spines["top"].set_visible(False)
    ax_ri.spines["right"].set_visible(False)
    ax_ri.spines["left"].set_visible(False)
    ax_ri.set_yticks([])
    ax_ri.tick_params(axis='x', labelsize=6)

    ri_vals_clean = [ri_map.get(l.name, ri_map.get(data_name(l.name), np.nan)) for l in leaves]
    ri_vals_clean = [r for r in ri_vals_clean if not pd.isna(r)]
    max_ri = max(ri_vals_clean) if ri_vals_clean else 1

    for leaf in leaves:
        y = tip_pos[leaf.name]
        sp = leaf.name
        dn = data_name(sp)
        ri = ri_map.get(sp, ri_map.get(dn, np.nan))
        if pd.isna(ri):
            continue
        clade = get_clade(sp)
        ax_ri.scatter(ri, y, s=max(6, ri * 60), alpha=0.7,
                     c=CLADE_COLORS.get(clade, "#aaa"),
                     edgecolor="k", linewidth=0.15, zorder=3)
    ax_ri.set_xlim(-0.02, max_ri * 1.1)  # small left pad to avoid clipping
    ax_ri.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax_ri.grid(axis="x", color="#ddd", linewidth=0.4, zorder=1)

    # Light box borders for Shannon and RI panels
    for ax in [ax_shan, ax_ri]:
        for side in ["top", "right", "left"]:
            ax.spines[side].set_visible(True)
            ax.spines[side].set_color("#ccc")
            ax.spines[side].set_linewidth(0.5)

    # ── Legends ──────────────────────────────────────────────────────────
    # TE type legend
    te_handles = [mpatches.Patch(color=TE_COLORS[t], label=t) for t in TE_ORDER]
    ax_bars.legend(handles=te_handles, fontsize=5, loc="lower right",
                   framealpha=0.9, title="TE type", title_fontsize=6)

    # Clade legend — explicitly include ALL clades present in tree
    present_clades = sorted(set(get_clade(l.name) for l in leaves))
    clade_handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#ddd"),
                     label=CLADE_NAMES.get(c, c))
                     for c in present_clades if c in CLADE_COLORS]
    n_clades = len(clade_handles)
    print(f"[INFO] {n_clades} clades in legend: {present_clades}")
    fig.legend(handles=clade_handles, fontsize=5, loc="upper left",
               bbox_to_anchor=(0.01, 0.99), ncol=2, framealpha=0.9,
               title="Clades", title_fontsize=6)

    fig.suptitle("Cladogram with TE composition, Shannon diversity, and Recency Index",
                 fontsize=11, y=0.998)
    fig.tight_layout(rect=[0, 0, 1, 0.995])

    out = OUTDIR / "cladogram_te_types_shannon_ri_noParsa.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "cladogram_te_types_shannon_ri_noParsa.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[DONE] {out.name} ({n_tips} species)")


if __name__ == "__main__":
    main()
