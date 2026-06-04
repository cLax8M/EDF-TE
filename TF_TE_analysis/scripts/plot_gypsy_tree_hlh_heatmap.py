#!/usr/bin/env python3
"""
Figure: Gypsy RT-domain phylogenetic tree (IQ-TREE) coloured by TEsorter clade,
with an adjacent heatmap of DAP-seq peak overlap and HLH_Myc motif hits per element.

Panel A: Cladogram (circular or rectangular) — nodes coloured by Gypsy RT clade
Panel B: Strip heatmap adjacent to tree leaves:
  Rows = 103 elements with RT domain (same order as tree leaves)
  Columns = 3 bHLH_Myc TFs (peak overlap: yes/no) + 3 HLH_Myc motifs (hit count)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from Bio import Phylo
from io import StringIO
from pathlib import Path

# ─── Paths ────────────────────────────────────────────────────────────────────
TREE_FILE   = Path("/tmp/gypsy_RT_tree.treefile")
CLADES_TSV  = Path("/tmp/gypsy_rt_clades.tsv")
FIMO_DIR    = Path("/tmp")
PEAK_DIR    = Path("/tmp")
OUT_DIR     = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# TFs and their motifs
TFS = ["bHLH_Myc_6330", "bHLH_Myc_4449", "bHLH_Myc_8146"]
# Top HLH_Myc motifs by enrichment in LTR/Gypsy
MOTIFS = ["HLH_Myc14", "HLH_Myc5", "HLH_Myc15"]

# RT clade colours
CLADE_COLORS = {
    "non-chromo-outgroup": "#4878CF",   # blue
    "Chlamyvir":           "#6ACC65",   # green
    "chromo-outgroup":     "#D65F5F",   # red
    "Reina":               "#B47CC7",   # purple
    "unclassified":        "#888888",   # grey
}

# ─── Load clade annotations ──────────────────────────────────────────────────
clades_df = pd.read_csv(CLADES_TSV, sep="\t")
clades_df["safe_id"] = clades_df["safe_id"].str.strip()
clade_map = dict(zip(clades_df["safe_id"], clades_df["rt_clade"]))
rt_clade_map = dict(zip(clades_df["safe_id"], clades_df["rt_clade"]))

# ─── Load FIMO hits: motif hits per element ───────────────────────────────────
def load_fimo_hits(motif):
    fimo_tsv = FIMO_DIR / f"fimo_gypsy_{motif}" / "fimo.tsv"
    df = pd.read_csv(fimo_tsv, sep="\t", comment="#")
    # sequence_name col
    seqcol = "sequence_name"
    if seqcol not in df.columns:
        seqcol = df.columns[2]
    # Convert to safe_id
    df["safe_id"] = (df[seqcol]
                     .str.split("::").str[0]
                     .str.replace(":", "_")
                     .str.replace("(","").str.replace(")",""))
    return df.groupby("safe_id").size().rename(motif)

motif_data = {}
for motif in MOTIFS:
    s = load_fimo_hits(motif)
    motif_data[motif] = s

# ─── Load DAP-seq peak overlap ────────────────────────────────────────────────
def load_peaks(tf):
    peak_file = PEAK_DIR / f"peak_overlap_{tf}.txt"
    ids = set()
    with open(peak_file) as fh:
        for line in fh:
            ids.add(line.strip())
    return ids

peak_sets = {tf: load_peaks(tf) for tf in TFS}

# ─── Load and root tree ───────────────────────────────────────────────────────
tree = Phylo.read(str(TREE_FILE), "newick")
tree.root_at_midpoint()
tree.ladderize()

# Get leaf order from tree
leaves = [c.name for c in tree.get_terminals()]
# Sanitize leaf names (IQ-TREE may quote or modify them)
def sanitize_leaf(name):
    if not name:
        return name
    clean = name.strip("'").replace("::","__").replace(":","_").replace("(","").replace(")","").rstrip("-")
    # IQ-TREE doubles the ID in leaf name (e.g. scaffold_1_78095-82595__scaffold_1_78095-82595)
    # Extract just the first part (before __)
    if "__" in clean:
        clean = clean.split("__")[0]
    return clean

leaves_safe = [sanitize_leaf(l) for l in leaves]
n_leaves = len(leaves_safe)

print(f"Tree leaves: {n_leaves}")
print(f"Example leaf: {leaves_safe[0]}")
print(f"Example clade key: {list(clade_map.keys())[0]}")

# Build matrix: rows = leaves, cols = TFs + motifs
col_labels = TFS + MOTIFS
data_matrix = np.zeros((n_leaves, len(col_labels)))

for i, sid in enumerate(leaves_safe):
    for j, tf in enumerate(TFS):
        matched = sid in peak_sets[tf]
        data_matrix[i, j] = 1.0 if matched else 0.0
    for j, motif in enumerate(MOTIFS, start=len(TFS)):
        hits = motif_data[motif].get(sid, 0)
        data_matrix[i, j] = min(hits, 5)   # cap at 5

# ─── Plot ─────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
gs = GridSpec(1, 2, figure=fig, width_ratios=[3, 1.2], wspace=0.05)

ax_tree = fig.add_subplot(gs[0])
ax_heat = fig.add_subplot(gs[1])

# --- Draw tree manually -------------------------------------------------------
# Use BioPython draw_ascii approach but we need matplotlib coords
# Draw horizontal rectangular cladogram
def get_leaf_positions(tree):
    """Assign y-positions to leaves (0..n-1 from bottom to top)."""
    leaves_list = tree.get_terminals()
    n = len(leaves_list)
    pos = {}
    for i, leaf in enumerate(leaves_list):
        pos[id(leaf)] = i
    return pos, leaves_list

def assign_y(clade, pos, internal_y):
    """Recursively assign y positions for internal nodes (mean of children)."""
    if clade.is_terminal():
        return pos[id(clade)]
    child_ys = [assign_y(c, pos, internal_y) for c in clade.clades]
    y = np.mean(child_ys)
    internal_y[id(clade)] = y
    return y

def get_x_from_root(tree):
    """Assign x positions (cumulative branch lengths)."""
    x_pos = {}
    def traverse(clade, x):
        x_pos[id(clade)] = x
        for child in clade.clades:
            bl = child.branch_length if child.branch_length else 0
            traverse(child, x + bl)
    traverse(tree.root, 0)
    return x_pos

leaf_pos, leaves_ordered = get_leaf_positions(tree)
internal_y = {}
assign_y(tree.root, leaf_pos, internal_y)
x_pos = get_x_from_root(tree)
# Assign y to terminals too
for leaf in leaves_ordered:
    internal_y[id(leaf)] = leaf_pos[id(leaf)]

# Max x for scaling
max_x = max(x_pos.values())

def draw_clade(ax, clade, leaf_pos, internal_y, x_pos, clade_map):
    """Draw branches recursively."""
    y_self = internal_y[id(clade)]
    x_self = x_pos[id(clade)]
    color = "#444444"
    if clade.is_terminal():
        sid = sanitize_leaf(clade.name)
        rt_c = clade_map.get(sid, "unclassified")
        color = CLADE_COLORS.get(rt_c, "#888888")
        # Leaf dot
        ax.plot(x_self, y_self, "o", ms=4, color=color, zorder=3)
    for child in clade.clades:
        y_child = internal_y[id(child)]
        x_child = x_pos[id(child)]
        # Vertical connector
        ax.plot([x_self, x_self], [y_self, y_child], "-", color="#aaaaaa", lw=0.6, zorder=1)
        # Horizontal branch
        c_sid = sanitize_leaf(child.name) if child.is_terminal() else None
        c_rt = clade_map.get(c_sid, "unclassified") if c_sid else "unclassified"
        c_color = CLADE_COLORS.get(c_rt, "#888888") if child.is_terminal() else "#aaaaaa"
        ax.plot([x_self, x_child], [y_child, y_child], "-", color=c_color, lw=0.8, zorder=2)
        draw_clade(ax, child, leaf_pos, internal_y, x_pos, clade_map)

draw_clade(ax_tree, tree.root, leaf_pos, internal_y, x_pos, clade_map)

ax_tree.set_xlim(-0.01 * max_x, max_x * 1.05)
ax_tree.set_ylim(-1, n_leaves)
ax_tree.axis("off")
ax_tree.set_title("LTR/Gypsy RT-domain phylogeny\n(Rhizopus microsporus ATCC11559, ≥3 kb elements)",
                  fontsize=11, pad=8)

# Legend for clades
clade_patches = [mpatches.Patch(color=v, label=k) for k, v in CLADE_COLORS.items()]
ax_tree.legend(handles=clade_patches, title="Gypsy RT clade", loc="lower left",
               fontsize=7, title_fontsize=8, framealpha=0.8)

# --- Draw heatmap adjacent to tree -------------------------------------------
# Reorder data_matrix to match tree leaf order
leaves_safe_ordered = [sanitize_leaf(l.name) for l in leaves_ordered]
idx_map = {sid: i for i, sid in enumerate(leaves_safe)}
reordered = np.zeros_like(data_matrix)
for tree_row, sid in enumerate(leaves_safe_ordered):
    orig_row = idx_map.get(sid)
    if orig_row is not None:
        reordered[tree_row] = data_matrix[orig_row]

# TF columns: binary (0/1)
tf_matrix  = reordered[:, :len(TFS)]
# Motif columns: 0-5
mot_matrix = reordered[:, len(TFS):]

# Use imshow with discrete colourmaps
# TF: white/darkblue
# Motif: white/darkred gradient

combined = np.hstack([tf_matrix, mot_matrix / 5.0])

cmap_tf  = mcolors.LinearSegmentedColormap.from_list("tf",  ["#f5f5f5", "#2166ac"])
cmap_mot = mcolors.LinearSegmentedColormap.from_list("mot", ["#f5f5f5", "#b2182b"])

# Draw two imshow panels side by side using axes_divider alternative: just plot on same ax
x_tf_start  = 0
x_mot_start = len(TFS) + 0.2

ax_heat.set_xlim(-0.3, len(col_labels) + 0.5)
ax_heat.set_ylim(-1, n_leaves)
ax_heat.axis("off")

cell_w = 0.85
cell_h = 0.85

# Peak (TF) columns
for j, tf in enumerate(TFS):
    for i in range(n_leaves):
        val = tf_matrix[i, j]
        color = cmap_tf(val)
        rect = plt.Rectangle([j, i - 0.4], cell_w, cell_h, color=color, ec="none")
        ax_heat.add_patch(rect)
    # Column label
    ax_heat.text(j + cell_w/2, n_leaves + 0.2, tf.replace("bHLH_Myc_","TF\n"),
                 ha="center", va="bottom", fontsize=6.5, rotation=45)

# Motif columns
for j, motif in enumerate(MOTIFS):
    x = len(TFS) + 0.3 + j
    for i in range(n_leaves):
        val = mot_matrix[i, j] / 5.0
        color = cmap_mot(val)
        rect = plt.Rectangle([x, i - 0.4], cell_w, cell_h, color=color, ec="none")
        ax_heat.add_patch(rect)
    # Column label
    ax_heat.text(x + cell_w/2, n_leaves + 0.2, motif,
                 ha="center", va="bottom", fontsize=6.5, rotation=45)

# Column separator line
ax_heat.axvline(len(TFS) + 0.15, color="#666666", lw=0.8, ls="--", ymin=0, ymax=1)

# Colourbar legend patches
tf_patch  = mpatches.Patch(color=cmap_tf(1.0),  label="DAP-seq peak overlap")
mot_patch = mpatches.Patch(color=cmap_mot(1.0), label="HLH_Myc motif hits (max=5)")
ax_heat.legend(handles=[tf_patch, mot_patch], loc="lower right", fontsize=7, framealpha=0.8)

ax_heat.set_title("Peak overlap\n& motif hits", fontsize=9)

plt.tight_layout()
for ext in ("pdf", "png"):
    out = OUT_DIR / f"gypsy_tree_hlh_heatmap.{ext}"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved: {out}")
