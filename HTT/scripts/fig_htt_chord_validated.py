#!/usr/bin/env python3
"""
fig_htt_chord_validated.py  –  Same chord design as fig_htt_chord.py
but restricted to the 15 community-validated transfer events (08_htt_clustering).
Arc width ∝ log(n pairs per community cluster).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from pathlib import Path as P

# ── Output ────────────────────────────────────────────────────────────────────
OUT_DIR = P("/legserv/Temp/Carlos/EDF-TE/results/HTT/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load VALIDATED community data ─────────────────────────────────────────────
COMM_TSV = "/legserv/Temp/Carlos/EDF-TE/results/HTT/08_htt_clustering/summary/htt_pairs_with_community.tsv"
comm = pd.read_csv(COMM_TSV, sep="\t")
comm["ca"] = comm[["cladeA", "cladeB"]].min(axis=1)
comm["cb"] = comm[["cladeA", "cladeB"]].max(axis=1)

# n = total pairs per (clade_pair, community); for arc width use n pairs per clade pair
pair_n = comm.groupby(["ca", "cb"]).agg(
    n=("pair_id", "count"),
    n_comm=("community_id", "nunique")
).reset_index()
print(f"Validated clade pairs: {len(pair_n)}  |  Total pairs: {pair_n['n'].sum()}")

# ── Clade definitions — ORDER and SIZES match the phylogenetic tree in htt_chord_v2.pdf
# Order: bottom→top in fan = clockwise in Circos from gap at top
# Colors: identical to CLADE_COL in plot_htt_chord.R
CLADES = [
    # (abbr,  full name,               R-script colour)
    ("HO", "Holozoa",               "#7B241C"),
    ("NU", "Nuclearida",            "#1C2833"),
    ("AP", "Aphelida",              "#27AE60"),
    ("CR", "Rozellomycota",         "#76D7C4"),
    ("MI", "Microsporidia",         "#5DADE2"),
    ("BL", "Blastocladiomycota",   "#E67E22"),
    ("CH", "Chytridiomycota",       "#E74C3C"),
    ("OL", "Olpidiomycota",         "#F7DC6F"),
    ("BS", "Basidiomycota",         "#1A5276"),
    ("ZO", "Zoopagomycota",         "#85C1E9"),
    ("EN", "Entomophthoromycota",  "#D35400"),
    ("KX", "Kickxellomycotina",     "#F39C12"),
    ("MO", "Mortierellomycota",     "#7D3C98"),
    ("GL", "Glomeromycotina",       "#D4AC0D"),
    ("MU", "Mucoromycotina",        "#BB8FCE"),
    ("AS", "Ascomycota",            "#1E8449"),
    ("BA", "Basidiobolales",        "#2471A3"),
]
abbr2col = {c[0]: c[2] for c in CLADES}
abbr2full = {c[0]: c[1] for c in CLADES}
abbr2idx  = {c[0]: i for i, c in enumerate(CLADES)}
N = len(CLADES)

# ── Angular positions ──────────────────────────────────────────────────────────
# Species counts from the phylogenetic tree (same as htt_chord_v2.pdf)
species_count = {
    "HO": 20, "NU": 2,  "AP": 2,  "CR": 2,  "MI": 22,
    "BL": 5,  "CH": 37, "OL": 1,  "BS": 2,  "ZO": 4,
    "EN": 5,  "KX": 13, "MO": 26, "GL": 7,  "MU": 63,
    "AS": 13, "BA": 14,
}
total      = sum(species_count.values())
GAP_DEG    = 2.0
arc_total  = 360.0 - GAP_DEG * N

starts, ends, mids = {}, {}, {}
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

def hex2rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) / 255 for i in (0, 2, 4))

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect("equal")
ax.axis("off")

R_OUTER = 1.00
R_INNER = 0.85
R_LABEL = 1.12
theta_res = 200

# ── Clade arcs ────────────────────────────────────────────────────────────────
for abbr, fullname, col in CLADES:
    s, e = starts[abbr], ends[abbr]
    thetas = np.linspace(np.radians(s), np.radians(e), theta_res)
    xo = np.cos(thetas) * R_OUTER;  yo = np.sin(thetas) * R_OUTER
    xi = np.cos(thetas[::-1]) * R_INNER;  yi = np.sin(thetas[::-1]) * R_INNER
    xs = np.concatenate([xo, xi, [xo[0]]])
    ys = np.concatenate([yo, yi, [yo[0]]])
    ax.fill(xs, ys, color=col, alpha=0.92, zorder=3)
    ax.plot(np.cos(thetas) * R_OUTER, np.sin(thetas) * R_OUTER,
            color="white", lw=0.5, zorder=4)
    mid = mids[abbr]
    lx, ly = deg2xy(mid, R_LABEL)
    ha = "left" if lx > 0.05 else ("right" if lx < -0.05 else "center")
    ax.text(lx, ly, abbr, fontsize=8.5, ha=ha, va="center",
            fontweight="bold", color=col, zorder=5)

# ── HTT chords (validated only) ───────────────────────────────────────────────
max_n = pair_n["n"].max()
for _, row in pair_n.iterrows():
    ca, cb, n = row["ca"], row["cb"], int(row["n"])
    if ca not in abbr2idx or cb not in abbr2idx:
        continue
    a1 = mids[ca];  a2 = mids[cb]
    x1, y1 = deg2xy(a1, R_INNER)
    x2, y2 = deg2xy(a2, R_INNER)
    cx, cy = 0.0, 0.0   # straight through centre

    alpha = 0.25 + 0.55 * (np.log1p(n) / np.log1p(max_n))
    lw    = 1.0  + 4.0  * (np.log1p(n) / np.log1p(max_n))
    c1, c2 = hex2rgb(abbr2col[ca]), hex2rgb(abbr2col[cb])
    col_blend = tuple((a + b) / 2 for a, b in zip(c1, c2))

    verts  = [(x1, y1), (cx, cy), (x2, y2)]
    codes  = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    patch  = mpatches.PathPatch(Path(verts, codes),
                                facecolor="none", edgecolor=col_blend,
                                linewidth=lw, alpha=alpha, zorder=2)
    ax.add_patch(patch)
    # Label n pairs
    mx = (x1 + x2 + cx) / 3
    my = (y1 + y2 + cy) / 3
    ax.text(mx, my, str(n), fontsize=6, ha="center", va="center",
            color=col_blend, alpha=min(1.0, alpha + 0.4), zorder=6,
            fontweight="bold")

# ── Legend ────────────────────────────────────────────────────────────────────
legend_patches = [mpatches.Patch(color=c[2], label=f"{c[0]}  {c[1]}")
                  for c in CLADES]
ax.legend(handles=legend_patches, loc="lower center",
          bbox_to_anchor=(0.5, -0.12), ncol=4,
          fontsize=7.5, frameon=False, title="Clades", title_fontsize=8)

ax.set_title(
    "Horizontal TE Transfers (HTT) — validated events\n"
    f"Arc width ∝ log(n pairs)  |  {pair_n['n'].sum()} pairs  |  {len(pair_n)} clade pairs  |  15 multi-genome communities",
    fontsize=10, fontweight="bold", pad=10)

ax.set_xlim(-1.35, 1.35)
ax.set_ylim(-1.35, 1.35)

out_pdf = OUT_DIR / "figHTT_chord_validated.pdf"
out_png = OUT_DIR / "figHTT_chord_validated.png"
fig.savefig(str(out_pdf), dpi=200, bbox_inches="tight")
fig.savefig(str(out_png), dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
