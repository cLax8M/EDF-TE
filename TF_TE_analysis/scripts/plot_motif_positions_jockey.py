#!/usr/bin/env python3
"""
Motif position figure for LINE/I-Jockey sequences (revised):
  Panel A: 4 full-length copies in Rhimi59 — only HMGb3 hits as tick marks
           (each tick = one FIMO hit; clustered ticks reflect overlapping windows
            inside a (CAA)n / (TGT)n tandem repeat)
  Panel B: 44 GremLINE1 centromeric copies from Mucor circinelloides MU_Mucci3
           — horizontal sequence lines with vertical tick marks at HMGb3 hit positions
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────
FIMO_RHIMI_HMG = Path("/tmp/fimo_rhimi59_jockey/fimo.tsv")
FIMO_GREM_HMG  = Path("/tmp/fimo_gremline1_hmgb3_nopgc/fimo.tsv")
FASTA_GREM     = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_centromeric.fasta")

OUT_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Sequence metadata ──────────────────────────────────────────────────────────
RHIMI_SEQS = [
    ("Copy 1 (448 bp, scf8)",  448,  "LINE/I-Jockey_448bp"),
    ("Copy 2 (331 bp, scf7)",  331,  "LINE/I-Jockey_331bp"),
    ("Copy 3 (280 bp, scf3)",  280,  "LINE/I-Jockey_280bp"),
    ("Copy 4 (244 bp, scf10)", 244,  "LINE/I-Jockey_244bp"),
]

# ── Load FIMO results ──────────────────────────────────────────────────────────
def load_fimo(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "motif_id":
                continue
            rows.append(parts)
    if not rows:
        return pd.DataFrame(columns=["motif_id","alt","seq_name","start","stop","strand","score","pval","qval","matched"])
    cols = ["motif_id","alt","seq_name","start","stop","strand","score","pval","qval","matched"]
    df = pd.DataFrame(rows, columns=cols[:len(rows[0])])
    df["start"] = pd.to_numeric(df["start"])
    df["stop"]  = pd.to_numeric(df["stop"])
    df["pval"]  = pd.to_numeric(df["pval"], errors="coerce")
    return df

# ── Load data ──────────────────────────────────────────────────────────────────
rhimi_hmg = load_fimo(FIMO_RHIMI_HMG)

# GremLINE1 FIMO + sequence lengths from FASTA
grem_hmg = load_fimo(FIMO_GREM_HMG)
grem_lens = {}   # header (full, with ::coords) -> length
with open(FASTA_GREM) as fh:
    name, seq = None, []
    for line in fh:
        if line.startswith(">"):
            if name:
                grem_lens[name] = len("".join(seq))
            name = line[1:].strip()
            seq  = []
        else:
            seq.append(line.strip())
    if name:
        grem_lens[name] = len("".join(seq))

# Map short seq_name (FIMO uses full header) -> length
# FIMO --no-pgc: seq_name in output = full FASTA header
print(f"GremLINE1 sequences: {len(grem_lens)}")
print(f"HMGb3 hits in GremLINE1: {len(grem_hmg)}")
print(f"HMGb3 hits in Rhimi59:   {len(rhimi_hmg)}")

# ── Figure layout ──────────────────────────────────────────────────────────────
# Panel A: 4 rows (one per copy)   — short, square-ish
# Panel B: 44 rows (one per copy)  — tall
N_GREM     = len(grem_lens)
ROW_H_A    = 0.55   # inches per row in panel A
ROW_H_B    = 0.20   # inches per row in panel B
fig_w      = 12.0
fig_h_a    = 4 * ROW_H_A + 0.7   # 4 copies + labels
fig_h_b    = N_GREM * ROW_H_B + 0.9

fig = plt.figure(figsize=(fig_w, fig_h_a + fig_h_b + 0.5))
gs  = fig.add_gridspec(2, 1,
                       height_ratios=[fig_h_a, fig_h_b],
                       hspace=0.45)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])

COL_HMG  = "#CC2222"   # dark red
TICK_H   = 0.14        # same height as backbone bar (y ± 0.07)
BACK_Y   = 0.5         # y-position of backbone within each row (0–1 normalised)

# ──────────────────────────────────────────────────────────────────────────────
# PANEL A: 4 Rhimi59 LINE/I-Jockey copies — HMGb3 tick marks only
# ──────────────────────────────────────────────────────────────────────────────
max_len_a = max(l for _, l, _ in RHIMI_SEQS)
n_a       = len(RHIMI_SEQS)
gap_a     = 1.0   # distance between row centres in data units

for i, (label, seqlen, seq_prefix) in enumerate(RHIMI_SEQS):
    y = (n_a - 1 - i) * gap_a

    # Backbone line (grey rectangle)
    ax_a.add_patch(plt.Rectangle((0, y - 0.07), seqlen, 0.14,
                                  color="#BBBBBB", zorder=1))
    # Right arrowhead for directionality
    ax_a.annotate("", xy=(seqlen + max_len_a*0.005, y),
                  xytext=(seqlen, y),
                  arrowprops=dict(arrowstyle="-|>", color="#888888",
                                  lw=1.2, mutation_scale=10),
                  zorder=2)

    # HMGb3 hit tick marks (thin vertical lines at start position)
    sub = rhimi_hmg[rhimi_hmg["seq_name"].str.startswith(seq_prefix)]
    for _, row in sub.iterrows():
        pos = row["start"]    # 1-based start
        ax_a.plot([pos, pos], [y - TICK_H/2, y + TICK_H/2],
                  color=COL_HMG, lw=1.8, alpha=0.9, zorder=3)

    # Label: left of backbone
    ax_a.text(-max_len_a * 0.015, y, label,
              ha="right", va="center", fontsize=8)
    # Length label: right of backbone
    ax_a.text(seqlen + max_len_a * 0.025, y, f"{seqlen} bp",
              ha="left", va="center", fontsize=7.5, color="#555555")

ax_a.set_xlim(-max_len_a * 0.12, max_len_a * 1.09)
ax_a.set_ylim(-gap_a * 0.6, (n_a - 1) * gap_a + gap_a * 0.6)
ax_a.set_xlabel("Position within element (bp)", fontsize=9)
ax_a.set_yticks([])
for sp in ax_a.spines.values():
    sp.set_visible(False)
ax_a.spines["bottom"].set_visible(True)
ax_a.tick_params(axis="x", labelsize=8)

legend_a = [mpatches.Patch(color=COL_HMG, label="HMGb3 motif hit (p < 0.01)")]
ax_a.legend(handles=legend_a, loc="upper right", fontsize=8, framealpha=0.9)
ax_a.set_title(
    "A   HMGb3 motif positions in 4 LINE/I-Jockey copies"
    "  (Rhizopus microsporus ATCC11559)",
    fontsize=9, loc="left", fontweight="bold", pad=4)

# ──────────────────────────────────────────────────────────────────────────────
# PANEL B: 44 GremLINE1 copies — individual linear maps with HMGb3 tick marks
# ──────────────────────────────────────────────────────────────────────────────
# Sort by family number, then by copy number within family
def sort_key(header):
    short = header.split("::")[0]   # e.g. GremLINE1_1.03
    parts = short.replace("GremLINE1_", "").split(".")
    fam  = int(parts[0]) if parts[0].isdigit() else 99
    copy = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
    return (fam, copy)

grem_sorted = sorted(grem_lens.keys(), key=sort_key)
max_len_b   = max(grem_lens.values())
n_b         = len(grem_sorted)
gap_b       = 1.0

TICK_H_B    = 0.12        # same height as backbone
BACK_H_B    = 0.12

for i, header in enumerate(grem_sorted):
    seqlen = grem_lens[header]
    y      = (n_b - 1 - i) * gap_b
    short  = header.split("::")[0]   # GremLINE1_1.03

    # Backbone (scaled to actual length on x-axis 0..max_len_b)
    ax_b.add_patch(plt.Rectangle((0, y - BACK_H_B/2), seqlen, BACK_H_B,
                                  color="#BBBBBB", zorder=1))

    # HMGb3 tick marks — match by full header (FIMO seq_name = full header)
    sub = grem_hmg[grem_hmg["seq_name"] == header]
    for _, row in sub.iterrows():
        pos = row["start"]
        ax_b.plot([pos, pos], [y - TICK_H_B/2, y + TICK_H_B/2],
                  color=COL_HMG, lw=1.5, alpha=0.9, zorder=3)

    # Label on right
    ax_b.text(max_len_b * 1.005, y, short,
              ha="left", va="center", fontsize=4.5, color="#333333")

ax_b.set_xlim(-max_len_b * 0.01, max_len_b * 1.20)
ax_b.set_ylim(-gap_b * 0.5, (n_b - 1) * gap_b + gap_b * 0.5)
ax_b.set_xlabel("Position within element (bp)", fontsize=9)
ax_b.set_yticks([])
for sp in ax_b.spines.values():
    sp.set_visible(False)
ax_b.spines["bottom"].set_visible(True)
ax_b.tick_params(axis="x", labelsize=8)

legend_b = [mpatches.Patch(color=COL_HMG, label="HMGb3 motif hit (p < 0.0001)")]
ax_b.legend(handles=legend_b, loc="upper right", fontsize=8, framealpha=0.9)
ax_b.set_title(
    "B   HMGb3 motif positions in 44 GremLINE1 (LINE/I-Jockey) copies"
    "  (Mucor circinelloides CBS 402.73 / MU_Mucci3)",
    fontsize=9, loc="left", fontweight="bold", pad=4)

# ── Save ───────────────────────────────────────────────────────────────────────
out_pdf = OUT_DIR / "motif_positions_jockey_gremlin.pdf"
out_png = OUT_DIR / "motif_positions_jockey_gremlin.png"
fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
fig.savefig(out_png, bbox_inches="tight", dpi=200)
plt.close()
print(f"\nSaved: {out_pdf}")
print(f"Saved: {out_png}")
