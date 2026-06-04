#!/usr/bin/env python3
"""
fig_htt_te_per_clade.py  –  TE-class proportions per clade (HTT candidates)
For each clade: proportion of LINE / LTR / RC / SINE / DNA among all HTT
candidate pairs that involve that clade (as query or subject).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────────
BASE    = Path("/legserv/Temp/Carlos/EDF-TE/results/HTT")
HTT_TSV = BASE / "06_htt_candidates/htt_candidates/htt_candidates.tsv"
TE_INFO = BASE / "05_te_classification/all_tes_info.tsv"
OUT_DIR = BASE / "figures"
OUT_PDF = OUT_DIR / "figHTT_te_per_clade.pdf"
OUT_PNG = OUT_DIR / "figHTT_te_per_clade.png"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── phylogenetic order & clade colours ────────────────────────────────────────
TREE_CLADE_ORDER = ["HO","NU","AP","CR","MI","BL","CH","OL","BS","ZO",
                    "EN","KX","MO","GL","MU","AS","BA"]
CLADE_COLORS = {
    "HO":"#7B241C","NU":"#1C2833","AP":"#27AE60","CR":"#76D7C4","MI":"#5DADE2",
    "BL":"#E67E22","CH":"#E74C3C","OL":"#F7DC6F","BS":"#1A5276","ZO":"#85C1E9",
    "EN":"#D35400","KX":"#F39C12","MO":"#7D3C98","GL":"#D4AC0D","MU":"#BB8FCE",
    "AS":"#1E8449","BA":"#2471A3",
}
CLADE_NAMES = {
    "HO":"Holozoa","NU":"Nuclearida","AP":"Aphelida","CR":"Rozellomycota",
    "MI":"Microsporidia","BL":"Blastocladiomycota","CH":"Chytridiomycota",
    "OL":"Olpidiomycota","BS":"Basidiomycota","ZO":"Zoopagomycota",
    "EN":"Entomophthoromycota","KX":"Kickxellomycotina","MO":"Mortierellomycota",
    "GL":"Glomeromycotina","MU":"Mucoromycotina","AS":"Ascomycota",
    "BA":"Basidiobolales",
}

# ── TE class colours (same as enrichment figure) ───────────────────────────────
TE_CLASS_ORDER  = ["LINE", "LTR", "RC", "SINE", "DNA", "Other"]
TE_CLASS_COLORS = {
    "LINE":  "#c0392b",
    "LTR":   "#2980b9",
    "RC":    "#27ae60",
    "SINE":  "#f39c12",
    "DNA":   "#8e44ad",
    "Other": "#aaaaaa",
}

# ── load data ──────────────────────────────────────────────────────────────────
print("Loading HTT candidates…")
htt = pd.read_csv(HTT_TSV, sep="\t")

print("Loading TE classification (chunked)…")
chunks = pd.read_csv(TE_INFO, sep="\t", low_memory=False,
                     usecols=["te_id","te_class"],
                     chunksize=500_000)
te_class_map = {}
for chunk in chunks:
    for _, r in chunk.iterrows():
        te_class_map[r["te_id"]] = r["te_class"]

# ── merge TE class ─────────────────────────────────────────────────────────────
htt["te_class_raw"] = htt["qseqid"].map(te_class_map).fillna("Unknown")
MAIN = {"LINE", "LTR", "RC", "SINE", "DNA"}
htt["te_class"] = htt["te_class_raw"].apply(lambda x: x if x in MAIN else "Other")

# ── build per-clade counts ─────────────────────────────────────────────────────
# Each pair contributes to BOTH clades involved
rows = []
for _, r in htt.iterrows():
    rows.append({"clade": r["qclade"], "te_class": r["te_class"]})
    if r["sclade"] != r["qclade"]:          # avoid double-counting within-clade
        rows.append({"clade": r["sclade"], "te_class": r["te_class"]})

df_long = pd.DataFrame(rows)

pivot = (df_long.groupby(["clade","te_class"])
                .size()
                .unstack(fill_value=0)
                .reindex(columns=TE_CLASS_ORDER, fill_value=0))

# totals per clade
totals = pivot.sum(axis=1)

# proportions
props = pivot.div(totals, axis=0)

# ── order clades phylogenetically (HO at top) ──────────────────────────────────
present_clades = [c for c in TREE_CLADE_ORDER if c in props.index]
# For horizontal barh: index 0 → bottom, index n-1 → top
# So reverse the list → BA first (bottom), HO last (top)
ordered = present_clades[::-1]
props   = props.loc[ordered]
totals  = totals.loc[ordered]

# ── plot ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

left   = np.zeros(len(props))
n_rows = len(props)

for cls in TE_CLASS_ORDER:
    if cls not in props.columns:
        continue
    vals = props[cls].values
    bars = ax.barh(range(n_rows), vals, left=left,
                   color=TE_CLASS_COLORS[cls],
                   label=cls, edgecolor="none", height=0.7)
    # Annotate bars ≥ 5 %
    for i, (v, l) in enumerate(zip(vals, left)):
        if v >= 0.05:
            ax.text(l + v / 2, i, f"{v*100:.0f}%",
                    ha="center", va="center", fontsize=6.5,
                    color="white", fontweight="bold")
    left += vals

# ── axes labels (clade abbreviation, coloured) ─────────────────────────────────
ytick_labels = [f"{c}  (n={int(totals[c])})" for c in ordered]
ax.set_yticks(range(n_rows))
ax.set_yticklabels(ytick_labels, fontsize=8.5)
for tick, cl in zip(ax.get_yticklabels(), ordered):
    tick.set_color(CLADE_COLORS.get(cl, "black"))
    tick.set_fontweight("bold")

# Add clade full names as secondary y-axis ticks (right side)
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())
ax2.set_yticks(range(n_rows))
ax2.set_yticklabels([CLADE_NAMES.get(c, c) for c in ordered], fontsize=7.5)
ax2.tick_params(right=False)
for tick, cl in zip(ax2.get_yticklabels(), ordered):
    tick.set_color(CLADE_COLORS.get(cl, "black"))

ax.set_xlabel("Proportion of HTT candidates", fontsize=10)
ax.set_xlim(0, 1)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
ax.set_title("TE class composition of HTT candidates by clade",
             fontsize=11, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Legend
legend_patches = [mpatches.Patch(color=TE_CLASS_COLORS[c], label=c)
                  for c in TE_CLASS_ORDER if c in props.columns]
ax.legend(handles=legend_patches, loc="upper center",
          bbox_to_anchor=(0.5, -0.12), ncol=len(TE_CLASS_ORDER),
          frameon=False, fontsize=9, title="TE class",
          title_fontsize=9)

plt.tight_layout()
fig.savefig(OUT_PDF, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"Saved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")

# ── text summary ───────────────────────────────────────────────────────────────
print("\n=== TE class proportions per clade (% LINE | LTR | RC | SINE | DNA) ===")
for cl in present_clades:  # top-to-bottom order
    p = props.loc[cl]
    n = int(totals[cl])
    print(f"  {cl:3s} (n={n:3d}): "
          f"LINE={p.get('LINE',0)*100:4.1f}% "
          f"LTR={p.get('LTR',0)*100:4.1f}% "
          f"RC={p.get('RC',0)*100:4.1f}% "
          f"SINE={p.get('SINE',0)*100:4.1f}% "
          f"DNA={p.get('DNA',0)*100:4.1f}%")
