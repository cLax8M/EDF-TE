#!/usr/bin/env python3
"""
fig_htt_chord.py  –  Chord diagram of HTT candidates on clade wheel
Clades are arranged in phylogenetic order; arcs show n HTT candidates
between clade pairs (width ∝ log(n)).  Colours match plot_circular_tree_v43.R.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patheffects as pe
from pathlib import Path as P

# ── Output ────────────────────────────────────────────────────────────────────
OUT_DIR = P("/legserv/Temp/Carlos/EDF-TE/results/HTT/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load HTT data ─────────────────────────────────────────────────────────────
HTT_TSV = "/legserv/Temp/Carlos/EDF-TE/results/HTT/06_htt_candidates/htt_candidates/htt_candidates.tsv"
htt = pd.read_csv(HTT_TSV, sep="\t")
htt["ca"] = htt[["qclade", "sclade"]].min(axis=1)
htt["cb"] = htt[["qclade", "sclade"]].max(axis=1)
pair_n = htt.groupby(["ca", "cb"]).size().reset_index(name="n")

# ── Clade definitions ─────────────────────────────────────────────────────────
# Phylogenetic order (clockwise from top), colours from plot_circular_tree_v43.R
CLADES = [
    ("HO", "Holozoa",             "#b4d0e8"),
    ("AP", "Aphelida",            "#5c9e3a"),
    ("NU", "Nucleariida",         "#1a2850"),
    ("CH", "Chytridiomycota",     "#e0d040"),
    ("BL", "Blastocladiomycota",  "#7030a0"),
    ("ZO", "Zoopagomycota",       "#c8dc70"),
    ("EN", "Entomophthoromycota", "#e01010"),
    ("BA", "Basidiobolales",      "#d0b4e0"),
    ("KX", "Kickxellomycotina",   "#f0a820"),
    ("GL", "Glomeromycota",       "#8b5e48"),
    ("MO", "Mortierellomycota",   "#f4a8b0"),
    ("MU", "Mucoromycota",        "#e07840"),
    ("BS", "Basidiomycota",       "#4472c4"),
    ("AS", "Ascomycota",          "#e8887a"),
]
abbr2idx  = {c[0]: i for i, c in enumerate(CLADES)}
abbr2col  = {c[0]: c[2] for c in CLADES}
abbr2full = {c[0]: c[1] for c in CLADES}
N = len(CLADES)

# ── Angular positions ──────────────────────────────────────────────────────────
# Allocate arc width proportional to species count (proxy: n HTT involvements)
species_count = {
    "HO": 21, "AP": 2, "NU": 3, "CH": 32, "BL": 5, "ZO": 4,
    "EN": 5, "BA": 14, "KX": 13, "GL": 7, "MO": 26, "MU": 64,
    "BS": 2, "AS": 14,
}
total = sum(species_count.values())
GAP_DEG = 2.0          # gap between sectors (degrees)
gap_total = GAP_DEG * N
arc_total = 360.0 - gap_total

# Start at top (90°), go clockwise (negative direction in matplotlib)
starts = {}
ends   = {}
mids   = {}
angle = 90.0
for abbr, _, _ in CLADES:
    w = (species_count[abbr] / total) * arc_total
    starts[abbr] = angle
    ends[abbr]   = angle - w
    mids[abbr]   = angle - w / 2
    angle -= w + GAP_DEG

def deg2xy(deg, r=1.0):
    rad = np.radians(deg)
    return r * np.cos(rad), r * np.sin(rad)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect("equal")
ax.axis("off")

R_OUTER = 1.00   # outer radius of clade arcs
R_INNER = 0.85   # inner radius
R_LABEL = 1.12   # text label radius

# ── Draw clade arcs ────────────────────────────────────────────────────────────
theta_res = 200
for abbr, fullname, col in CLADES:
    s, e = starts[abbr], ends[abbr]
    thetas = np.linspace(np.radians(s), np.radians(e), theta_res)
    xo = np.cos(thetas) * R_OUTER
    yo = np.sin(thetas) * R_OUTER
    xi = np.cos(thetas[::-1]) * R_INNER
    yi = np.sin(thetas[::-1]) * R_INNER
    xs = np.concatenate([xo, xi, [xo[0]]])
    ys = np.concatenate([yo, yi, [yo[0]]])
    ax.fill(xs, ys, color=col, alpha=0.92, zorder=3)
    ax.plot(np.cos(thetas) * R_OUTER, np.sin(thetas) * R_OUTER,
            color="white", lw=0.5, zorder=4)
    # Label
    mid = mids[abbr]
    lx, ly = deg2xy(mid, R_LABEL)
    ha = "left" if lx > 0.05 else ("right" if lx < -0.05 else "center")
    ax.text(lx, ly, abbr, fontsize=8.5, ha=ha, va="center",
            fontweight="bold", color=col, zorder=5)

# ── Draw HTT chords ────────────────────────────────────────────────────────────
max_n = pair_n["n"].max()
for _, row in pair_n.iterrows():
    ca, cb, n = row["ca"], row["cb"], row["n"]
    if ca not in abbr2idx or cb not in abbr2idx:
        continue
    # Arc anchor points: midpoint of each sector's inner edge
    a1 = mids[ca]
    a2 = mids[cb]
    x1, y1 = deg2xy(a1, R_INNER)
    x2, y2 = deg2xy(a2, R_INNER)
    # Bezier control point at origin (gives nice circular bow)
    depth = 0.0  # 0 = straight through center
    cx, cy = depth * (x1 + x2) / 2, depth * (y1 + y2) / 2

    alpha = 0.15 + 0.55 * (np.log1p(n) / np.log1p(max_n))
    lw    = 0.5  + 3.5  * (np.log1p(n) / np.log1p(max_n))
    # Blend colours of the two clades
    def hex2rgb(h):
        h = h.lstrip("#")
        return tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))
    c1, c2 = hex2rgb(abbr2col[ca]), hex2rgb(abbr2col[cb])
    col_blend = tuple((a + b) / 2 for a, b in zip(c1, c2))

    verts = [(x1, y1), (cx, cy), (x2, y2)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    path  = Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor="none",
                               edgecolor=col_blend,
                               linewidth=lw, alpha=alpha, zorder=2)
    ax.add_patch(patch)
    # Annotate if n ≥ 10
    if n >= 10:
        mx = (x1 + x2 + cx) / 3
        my = (y1 + y2 + cy) / 3
        ax.text(mx, my, str(n), fontsize=5.5, ha="center", va="center",
                color=col_blend, alpha=min(1.0, alpha + 0.4), zorder=6,
                fontweight="bold")

# ── Legend for clade names ─────────────────────────────────────────────────────
legend_patches = [mpatches.Patch(color=c[2], label=f"{c[0]}  {c[1]}")
                  for c in CLADES]
ax.legend(handles=legend_patches, loc="lower center",
          bbox_to_anchor=(0.5, -0.12), ncol=4,
          fontsize=7.5, frameon=False, title="Clades", title_fontsize=8)

# ── Title ─────────────────────────────────────────────────────────────────────
ax.set_title("Horizontal TE Transfers (HTT) between fungal clades\n"
             "Arc width ∝ log(n candidates)  |  223 candidates  |  14 clades",
             fontsize=11, fontweight="bold", pad=10)

ax.set_xlim(-1.35, 1.35)
ax.set_ylim(-1.35, 1.35)

fig.savefig(str(OUT_DIR / "figHTT_chord.pdf"), dpi=200, bbox_inches="tight")
fig.savefig(str(OUT_DIR / "figHTT_chord.png"), dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_DIR}/figHTT_chord.pdf")
print(f"Saved: {OUT_DIR}/figHTT_chord.png")
