#!/usr/bin/env python3
"""
F02_plot_dnmt_tree.py
=====================
Plot a cladogram of Mucoromycota (MU) + Glomeromycota (GL) species with:
  - Corrected synteny ratio (patristic-corrected, from perspecies_corrected_synteny.tsv)
  - C-5 cytosine DNA methyltransferase (DNMT) orthogroup presence/absence
    detected via IPR001525 / IPR018117 signatures

Output:
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures/
      F02_dnmt_mucorales_tree.pdf
      F02_dnmt_mucorales_tree.png

Usage:
    python F02_plot_dnmt_tree.py
"""

import os, re, glob
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable

# ── paths ──────────────────────────────────────────────────────────────────────
_DATA    = "/legserv/Temp/Carlos/EDF-TE"
TREE_F   = f"{_DATA}/results/SpeciesTree/final_LG_G_rotated.nwk"
SYNTENY  = f"{_DATA}/results/Synteny_guide_corrected/perspecies_corrected_synteny.tsv"
IPR_DIR  = f"{_DATA}/data/functional_interpro"
OG_FILE  = f"{_DATA}/results/Orthofinder_run/Orthogroups/Orthogroups.tsv"
OUT_DIR  = f"{_DATA}/results/Section2_TElandscape/figures"

TARGET_IPRS = {"IPR001525", "IPR018117"}   # C-5 cytosine methyltransferase

# ── 1. Collect C-5 methyltransferase proteins per species (MU+GL only) ────────
print("[1/6] Collecting methyltransferase proteins ...", flush=True)
species_meth: dict = defaultdict(set)   # full_sp -> set(protein_id)
for fpath in sorted(glob.glob(os.path.join(IPR_DIR, "*_IPR.tab"))):
    sp = os.path.basename(fpath).replace("_IPR.tab", "")
    if not (sp.startswith("MU_") or sp.startswith("GL_")):
        continue
    with open(fpath) as fh:
        hdr = fh.readline().strip().split("\t")
        hdr = [c.lstrip("#") for c in hdr]
        if "iprId" not in hdr or "proteinId" not in hdr:
            continue
        ii, pi = hdr.index("iprId"), hdr.index("proteinId")
        for line in fh:
            fs = line.split("\t")
            if len(fs) > max(ii, pi) and fs[ii] in TARGET_IPRS:
                species_meth[sp].add(fs[pi])

mu_gl_sps = sorted(species_meth)
print(f"  {len(mu_gl_sps)} species with DNMT hits, "
      f"{sum(len(v) for v in species_meth.values())} total proteins")

# lookup set for fast OG scanning
sp_pid_lookup = set()
sp_short_map: dict = {}   # sp_short -> full_sp
for sp, pids in species_meth.items():
    sp_short = sp.split("_", 1)[1] if "_" in sp else sp
    sp_short_map[sp_short] = sp
    for p in pids:
        sp_pid_lookup.add((sp_short, p))

# ── 2. Scan Orthogroups.tsv ────────────────────────────────────────────────────
print("[2/6] Scanning Orthogroups.tsv ...", flush=True)
og_hits: dict = defaultdict(lambda: defaultdict(set))   # OG -> full_sp -> prot IDs
with open(OG_FILE) as fh:
    fh.readline()
    for ln, line in enumerate(fh, 2):
        if ln % 30000 == 0:
            print(f"  ... {ln:,} OGs", flush=True)
        for m in re.finditer(r'jgi\|(\w+)\|(\d+)\|', line):
            sp_short, pid = m.group(1), m.group(2)
            if (sp_short, pid) in sp_pid_lookup:
                og = line.split('\t', 1)[0]
                og_hits[og][sp_short_map[sp_short]].add(pid)
print(f"  Found {len(og_hits)} DNMT-containing OGs")

# ── 3. Build presence/absence matrix ──────────────────────────────────────────
print("[3/6] Building presence/absence matrix ...", flush=True)
# Keep OGs present in at least 2 MU/GL species
all_tree_sps = set()
with open(TREE_F) as f:
    tree_str = f.read()
tips = [x for x in re.findall(r'([A-Za-z0-9_.-]+):', tree_str)]
mu_gl_tree = sorted(t for t in tips if t.startswith("MU_") or t.startswith("GL_"))

og_df = pd.DataFrame(0, index=mu_gl_tree,
                     columns=sorted(og_hits.keys()))
for og, sp_dict in og_hits.items():
    for sp in sp_dict:
        if sp in og_df.index:
            og_df.loc[sp, og] = 1

# Filter to OGs with presence in ≥2 species
og_df = og_df.loc[:, og_df.sum(axis=0) >= 2]
# Sort columns by prevalence (most prevalent first)
og_df = og_df[og_df.sum().sort_values(ascending=False).index]
print(f"  {og_df.shape[1]} OGs present in ≥2 MU/GL species")

# ── 4. Load synteny data ───────────────────────────────────────────────────────
print("[4/6] Loading corrected synteny ratios ...", flush=True)
syn_df = pd.read_csv(SYNTENY, sep="\t")
syn_map = dict(zip(syn_df["species"], syn_df["median_r_corrected"].astype(float)))

# ── 5. Parse + prune tree to MU+GL ────────────────────────────────────────────
print("[5/6] Parsing and pruning tree ...", flush=True)

# Minimal Newick parser (no ete3/Bio dependency here)
class Node:
    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children: list = []
        self.x = 0.0
        self.y = 0.0

def parse_newick(s: str) -> Node:
    """Iterative stack-based Newick parser (safe for deep trees)."""
    s = s.strip().rstrip(";").rstrip()
    stack = []
    current = Node()
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c == "(":
            new_node = Node()
            current.children.append(new_node)
            stack.append(current)
            current = new_node
            i += 1
        elif c == ",":
            # finish current leaf/internal, start sibling
            parent = stack[-1]
            new_node = Node()
            parent.children.append(new_node)
            current = new_node
            i += 1
        elif c == ")":
            # go back to parent
            current = stack.pop()
            i += 1
            # read optional label after ")"
            start = i
            while i < n and s[i] not in "(),;":
                i += 1
            label = s[start:i]
            if label:
                if ":" in label:
                    parts = label.rsplit(":", 1)
                    current.name = parts[0]
                    try:
                        current.length = float(parts[1])
                    except ValueError:
                        pass
                else:
                    current.name = label
        elif c == ";":
            break
        else:
            # leaf label
            start = i
            while i < n and s[i] not in "(),;":
                i += 1
            label = s[start:i]
            if ":" in label:
                parts = label.rsplit(":", 1)
                current.name = parts[0]
                try:
                    current.length = float(parts[1])
                except ValueError:
                    pass
            else:
                current.name = label
    return current

def get_leaves(node):
    if not node.children:
        return [node]
    out = []
    for c in node.children:
        out.extend(get_leaves(c))
    return out

def prune_to(node, keep_set):
    """Remove all leaves not in keep_set; also collapse internal nodes."""
    if not node.children:
        return node if node.name in keep_set else None
    new_children = []
    for c in node.children:
        pruned = prune_to(c, keep_set)
        if pruned is not None:
            new_children.append(pruned)
    node.children = new_children
    if not new_children:
        return None
    if len(new_children) == 1:
        # collapse singleton internal node
        return new_children[0]
    return node

def get_max_depth(node, depth=0):
    if not node.children:
        return depth
    return max(get_max_depth(c, depth + 1) for c in node.children)

def assign_cladogram_pos(node, y_counter, depth=0, max_depth=None):
    if not node.children:
        node.x = max_depth
        node.y = y_counter[0]
        y_counter[0] += 1
    else:
        node.x = depth
        for c in node.children:
            assign_cladogram_pos(c, y_counter, depth + 1, max_depth)
        node.y = np.mean([c.y for c in node.children])

def draw_cladogram(ax, node):
    tip_pos = {}
    def _draw(n):
        if not n.children:
            tip_pos[n.name] = n.y
        for c in n.children:
            ax.plot([n.x, c.x], [n.y, n.y], "k-", lw=0.6)
            ax.plot([c.x, c.x], [n.y, c.y], "k-", lw=0.6)
            _draw(c)
    _draw(node)
    return tip_pos

root = parse_newick(open(TREE_F).read())
keep = set(mu_gl_tree)
pruned = prune_to(root, keep)
md = get_max_depth(pruned)
assign_cladogram_pos(pruned, [0], max_depth=md)
leaves = get_leaves(pruned)
n_tips = len(leaves)

# ── 6. Plot ────────────────────────────────────────────────────────────────────
print("[6/6] Plotting ...", flush=True)
os.makedirs(OUT_DIR, exist_ok=True)

n_og = og_df.shape[1]
fig_h = max(14, n_tips * 0.15)

# Column widths: clade label | tree | synteny bar | DNMT columns | colorbar
ncols_og = min(n_og, 30)   # show at most 30 OGs
width_ratios = [0.25, 2.5, 0.35] + [0.12] * ncols_og + [0.3]
fig, axes = plt.subplots(1, len(width_ratios),
                         figsize=(len(width_ratios) * 0.55 + 4, fig_h),
                         gridspec_kw={"width_ratios": width_ratios, "wspace": 0.01})

ax_clade, ax_tree, ax_syn = axes[0], axes[1], axes[2]
ax_og = axes[3:3 + ncols_og]
ax_cbar = axes[-1]

# ─ cladogram ─────────────────────────────────────────────────────────────────
tip_pos = draw_cladogram(ax_tree, pruned)
ax_tree.set_xlim(-0.5, md + 0.5)
ax_tree.set_ylim(-0.5, n_tips - 0.5)
ax_tree.invert_yaxis()
ax_tree.axis("off")

# ─ clade labels ──────────────────────────────────────────────────────────────
CLADE_COLORS = {"MU": "#5d8aa8", "GL": "#8cb369"}
for leaf in leaves:
    y = tip_pos[leaf.name]
    clade = leaf.name.split("_")[0]
    ax_clade.add_patch(mpatches.Rectangle((0, y - 0.45), 1, 0.9,
                                          fc=CLADE_COLORS.get(clade, "#aaa"),
                                          ec="none"))
    ax_clade.text(0.5, y, clade, ha="center", va="center",
                  fontsize=4.5, color="white", fontweight="bold")

ax_clade.set_xlim(0, 1)
ax_clade.set_ylim(-0.5, n_tips - 0.5)
ax_clade.invert_yaxis()
ax_clade.axis("off")

# ─ synteny bar ───────────────────────────────────────────────────────────────
SYN_CMAP = plt.cm.RdYlGn
SYN_NORM = mcolors.Normalize(vmin=0, vmax=1)
for leaf in leaves:
    y = tip_pos[leaf.name]
    val = syn_map.get(leaf.name, np.nan)
    if not np.isnan(val):
        color = SYN_CMAP(SYN_NORM(val))
        ax_syn.add_patch(mpatches.Rectangle((0, y - 0.45), val, 0.9,
                                            fc=color, ec="none"))
ax_syn.set_xlim(0, 1)
ax_syn.set_ylim(-0.5, n_tips - 0.5)
ax_syn.invert_yaxis()
ax_syn.set_xticks([0, 0.5, 1])
ax_syn.set_xticklabels(["0", "0.5", "1"], fontsize=5)
ax_syn.tick_params(axis="y", left=False, labelleft=False)
ax_syn.set_xlabel("Corrected\nsynteny", fontsize=6, labelpad=2)
ax_syn.spines["top"].set_visible(False)
ax_syn.spines["right"].set_visible(False)
ax_syn.spines["left"].set_visible(False)

# tip labels on ax_syn
for leaf in leaves:
    y = tip_pos[leaf.name]
    sp_label = leaf.name.replace("MU_", "").replace("GL_", "")
    ax_tree.text(md + 0.05, y, leaf.name, ha="left", va="center",
                 fontsize=3.8, color="#222")

# ─ DNMT heatmap ──────────────────────────────────────────────────────────────
og_cols = og_df.columns[:ncols_og].tolist()
DNMT_COLORS = {0: "#f5f5f5", 1: "#c0392b"}
for col_i, og in enumerate(og_cols):
    ax = ax_og[col_i]
    for leaf in leaves:
        y = tip_pos[leaf.name]
        present = og_df.loc[leaf.name, og] if leaf.name in og_df.index else 0
        ax.add_patch(mpatches.Rectangle((0, y - 0.45), 1, 0.9,
                                        fc=DNMT_COLORS[int(present)], ec="none"))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n_tips - 0.5)
    ax.invert_yaxis()
    ax.axis("off")
    # OG label on top (rotated)
    ax.text(0.5, -0.8, og.replace("OG", ""), ha="center", va="top",
            fontsize=4, rotation=90, transform=ax.transData)

# ─ colorbar for synteny ──────────────────────────────────────────────────────
cb = ColorbarBase(ax_cbar, cmap=SYN_CMAP, norm=SYN_NORM,
                  orientation="vertical")
cb.set_label("Corrected\nsynteny ratio", fontsize=6)
cb.ax.tick_params(labelsize=5)

# ─ legend for DNMT ───────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(fc="#c0392b", ec="none", label="DNMT present"),
    mpatches.Patch(fc="#f5f5f5", ec="#ccc", lw=0.5, label="DNMT absent"),
    mpatches.Patch(fc=CLADE_COLORS["MU"], ec="none", label="Mucoromycota (MU)"),
    mpatches.Patch(fc=CLADE_COLORS["GL"], ec="none", label="Glomeromycota (GL)"),
]
fig.legend(handles=legend_patches, loc="lower center",
           ncol=4, fontsize=6, frameon=False,
           bbox_to_anchor=(0.5, -0.01))

fig.suptitle(
    "C-5 cytosine DNA methyltransferase (DNMT) orthogroup distribution\n"
    "Mucoromycota + Glomeromycota — corrected synteny ratio",
    fontsize=9, y=1.005)

for ext in ("pdf", "png"):
    out = os.path.join(OUT_DIR, f"F02_dnmt_mucorales_tree.{ext}")
    fig.savefig(out, bbox_inches="tight", dpi=200)
    print(f"  → {out}")

plt.close(fig)
print("Done.")
