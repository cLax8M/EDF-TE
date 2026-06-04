#!/usr/bin/env python3
"""
fig_htt_summary.py  –  HTT summary figure
3 panels:
  A – Ks distribution of 223 candidates, stacked by TE class
  B – Heatmap: n HTT candidates per clade pair (14×14 triangular)
  C – Minimal HTT events per node (bubble plot: node depth vs events, size=n_pairs)
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
BASE    = Path("/legserv/Temp/Carlos/EDF-TE/results/HTT")
HTT_TSV = BASE / "06_htt_candidates/htt_candidates/htt_candidates.tsv"
TE_INFO = BASE / "05_te_classification/all_tes_info.tsv"
EVENTS  = BASE / "09_minimal_htt_events/summary/htt_minimal_events_summary.txt"
OUT_DIR = BASE / "figures"
OUT_PDF = OUT_DIR / "figHTT_summary.pdf"
OUT_PNG = OUT_DIR / "figHTT_summary.png"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── load data ─────────────────────────────────────────────────────────────────
htt = pd.read_csv(HTT_TSV, sep="\t")
te_info = pd.read_csv(TE_INFO, sep="\t", low_memory=False)

# Merge TE class
merged = htt.merge(te_info[["te_id", "te_class", "family"]],
                   left_on="qseqid", right_on="te_id", how="left")
merged["te_class"] = merged["te_class"].fillna("Unknown")

# Node depth (busco_ks as proxy)
node_depth = merged.groupby("mrca_node")["busco_ks_threshold"].first()

# Minimal events per node (from summary)
events_per_node = pd.DataFrame({
    "node": ["node_32", "node_38", "node_25", "node_60", "node_37", "node_42", "node_61", "node_2"],
    "min_events": [6, 3, 1, 1, 1, 1, 1, 1]
})
events_per_node["busco_ks"] = events_per_node["node"].map(node_depth)
events_per_node["n_pairs"]  = events_per_node["node"].map(
    merged.groupby("mrca_node").size())
# Node labels: two main clades separated by "→"
node_labels = {
    "node_32": "multi-phylum\n(Ks≈1.43)",
    "node_38": "BL/BS/EN/KX\nvs MO/MU/ZO\n(Ks≈1.75)",
    "node_25": "BS/HO↔MO\n(Ks≈1.03)",
    "node_60": "MO↔MU\n(Ks≈1.11)",
    "node_37": "AS↔CH\n(Ks≈1.23)",
    "node_42": "CH↔NU\n(Ks≈1.55)",
    "node_61": "CH↔HO\n(Ks=n/a)",
    "node_2":  "AP↔fungi\n(Ks≈2.07)",
}
events_per_node["label"] = events_per_node["node"].map(node_labels)
events_per_node = events_per_node.sort_values("busco_ks", ascending=True)

# Clade pair counts
merged["clade_a"] = merged[["qclade", "sclade"]].min(axis=1)
merged["clade_b"] = merged[["qclade", "sclade"]].max(axis=1)
pair_counts = merged.groupby(["clade_a", "clade_b"]).size().reset_index(name="n")

# ── Phylogenetic clade order & colours (matches htt_chord_v2.pdf / R tree) ──────
TREE_CLADE_ORDER = ["HO","NU","AP","CR","MI","BL","CH","OL","BS","ZO","EN","KX","MO","GL","MU","AS","BA"]
CLADE_COLORS = {
    "HO":"#7B241C","NU":"#1C2833","AP":"#27AE60","CR":"#76D7C4","MI":"#5DADE2",
    "BL":"#E67E22","CH":"#E74C3C","OL":"#F7DC6F","BS":"#1A5276","ZO":"#85C1E9",
    "EN":"#D35400","KX":"#F39C12","MO":"#7D3C98","GL":"#D4AC0D","MU":"#BB8FCE",
    "AS":"#1E8449","BA":"#2471A3"
}

present_clades = set(merged["qclade"].tolist() + merged["sclade"].tolist())
all_clades = [c for c in TREE_CLADE_ORDER if c in present_clades]

hm_data = pd.DataFrame(0, index=all_clades, columns=all_clades, dtype=float)
for _, row in pair_counts.iterrows():
    if row["clade_a"] in all_clades and row["clade_b"] in all_clades:
        hm_data.loc[row["clade_a"], row["clade_b"]] = row["n"]
        hm_data.loc[row["clade_b"], row["clade_a"]] = row["n"]
arr = hm_data.to_numpy(copy=True, dtype=float)
np.fill_diagonal(arr, np.nan)
hm_data = pd.DataFrame(arr, index=hm_data.index, columns=hm_data.columns)

# ── TE class colours ──────────────────────────────────────────────────────────
TE_COLOURS = {
    "LINE":    "#c0392b",
    "LTR":     "#2980b9",
    "RC":      "#27ae60",
    "SINE":    "#f39c12",
    "DNA":     "#8e44ad",
    "Unknown": "#aaaaaa",
}
CLASS_ORDER = ["LINE", "LTR", "RC", "SINE", "DNA", "Unknown"]

# ── figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(20, 7.5))
gs  = fig.add_gridspec(1, 3, width_ratios=[1, 1.9, 1], wspace=0.45)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])
axC = fig.add_subplot(gs[2])

# ─── Panel A: Ks histogram ───────────────────────────────────────────────────
bins   = np.linspace(0, 1.8, 30)
bottom = np.zeros(len(bins) - 1)
for cls in CLASS_ORDER:
    vals = merged.loc[merged["te_class"] == cls, "te_ks"].dropna().values
    counts, _ = np.histogram(vals, bins=bins)
    axA.bar(bins[:-1], counts, width=np.diff(bins), bottom=bottom,
            color=TE_COLOURS[cls], label=cls, align="edge",
            edgecolor="white", linewidth=0.3)
    bottom += counts

axA.set_xlabel("TE Ks (synonymous substitutions / site)", fontsize=9)
axA.set_ylabel("HTT candidates (n)", fontsize=9)
axA.set_title("A", loc="left", fontweight="bold", fontsize=12)
axA.set_xlim(0, 1.8)
axA.tick_params(labelsize=8)
axA.legend(fontsize=7.5, frameon=False, title="TE class", title_fontsize=8)
# vertical median line
med_ks = merged["te_ks"].median()
axA.axvline(med_ks, color="#333333", lw=1.2, linestyle="--", alpha=0.8)
axA.text(med_ks + 0.03, axA.get_ylim()[1] * 0.92,
         f"median\n{med_ks:.2f}", fontsize=7, va="top", color="#333333")

# ─── Panel B: triangular interaction matrix ──────────────────────────────────
n_cl = len(all_clades)
# lower triangle only (mask upper triangle + diagonal)
mask = np.triu(np.ones((n_cl, n_cl), dtype=bool), k=0)
hm_plot = hm_data.values.copy().astype(float)
hm_plot[mask] = np.nan

max_val = np.nanmax(hm_plot)
im = axB.imshow(hm_plot, aspect="equal",
                cmap="YlOrRd",
                norm=LogNorm(vmin=1, vmax=max_val))

# White background for upper triangle (cleaner than default grey)
axB.set_facecolor("white")

# Diagonal line separating triangles
axB.plot([-0.5, n_cl - 0.5], [-0.5, n_cl - 0.5],
         color="grey", lw=0.8, ls="-", zorder=5)

axB.set_xticks(range(n_cl))
axB.set_xticklabels(all_clades, rotation=45, ha="right", fontsize=7.5)
axB.set_yticks(range(n_cl))
axB.set_yticklabels(all_clades, fontsize=7.5)
axB.set_title("B", loc="left", fontweight="bold", fontsize=12)

# Colour tick labels to match clade colours
for tick, cl in zip(axB.get_xticklabels(), all_clades):
    tick.set_color(CLADE_COLORS.get(cl, "black"))
    tick.set_fontweight("bold")
for tick, cl in zip(axB.get_yticklabels(), all_clades):
    tick.set_color(CLADE_COLORS.get(cl, "black"))
    tick.set_fontweight("bold")

# Annotate cells
for i in range(n_cl):
    for j in range(n_cl):
        if not np.isnan(hm_plot[i, j]) and hm_plot[i, j] > 0:
            v = int(hm_plot[i, j])
            c = "white" if v >= 15 else "black"
            axB.text(j, i, str(v), ha="center", va="center",
                     fontsize=6.5, color=c, fontweight="bold")

plt.colorbar(im, ax=axB, shrink=0.6, label="n HTT candidates", pad=0.02)

# ─── Panel C: minimal events bubble chart ────────────────────────────────────
ev = events_per_node.copy().reset_index(drop=True)
# Clip extreme busco_ks for display (node_61 = 9.99 is a huge outlier for CH-HO)
ev = ev[ev["busco_ks"] <= 2.0]
y_pos = range(len(ev))

sizes = ev["n_pairs"].values * 12   # scale bubble size
sc = axC.scatter(ev["busco_ks"].values, list(y_pos),
                 s=sizes,
                 c=ev["min_events"].values,
                 cmap="Blues", vmin=0.5, vmax=ev["min_events"].max() + 0.5,
                 edgecolors="#333333", linewidths=0.7, zorder=3)

axC.set_yticks(list(y_pos))
axC.set_yticklabels(ev["label"].tolist(), fontsize=7.5)
axC.set_xlabel("Node depth (BUSCO Ks threshold)", fontsize=9)
axC.set_title("C", loc="left", fontweight="bold", fontsize=12)
axC.grid(axis="x", linestyle="--", alpha=0.4)
axC.set_xlim(0.8, 2.0)
axC.invert_yaxis()

# colour-bar for number of events
cb = plt.colorbar(sc, ax=axC, shrink=0.55, label="Min. HTT events", pad=0.02)
cb.set_ticks([1, 2, 3, 4, 5, 6])

# bubble size legend
for ns, label in [(5, "5 pairs"), (20, "20 pairs"), (45, "45 pairs")]:
    axC.scatter([], [], s=ns * 12, c="#4a90d9",
                edgecolors="#333333", linewidths=0.7, label=label)
axC.legend(title="n candidates", fontsize=7.5, title_fontsize=7.5,
           loc="lower right", frameon=True, framealpha=0.8)

# ── save ──────────────────────────────────────────────────────────────────────
fig.suptitle(f"Horizontal Transposon Transfers (HTT): 223 candidates, "
             f"15 minimal events, 14 fungal clades",
             fontsize=10, fontweight="bold", y=1.01)
fig.savefig(str(OUT_PDF), dpi=200, bbox_inches="tight")
fig.savefig(str(OUT_PNG), dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")

# ── text summary ──────────────────────────────────────────────────────────────
print("\n=== HTT Summary ===")
print(f"Total candidates: {len(htt)}")
print(f"Ks range: {htt['te_ks'].min():.3f} – {htt['te_ks'].max():.3f}  "
      f"(median {htt['te_ks'].median():.3f})")
print("\nTE class:")
for cls in CLASS_ORDER:
    n = (merged["te_class"] == cls).sum()
    print(f"  {cls:8s}: {n:3d}  ({100*n/len(merged):.1f}%)")
print("\nTop clade pairs (n HTT):")
top = pair_counts.sort_values("n", ascending=False).head(10)
for _, r in top.iterrows():
    print(f"  {r['clade_a']}-{r['clade_b']:5s}: {int(r['n'])}")
print("\nMinimal events per node:")
for _, r in events_per_node.iterrows():
    print(f"  {r['node']}: {int(r['min_events'])} events  "
          f"({int(r['n_pairs'])} pairs, Ks≈{r['busco_ks']:.2f})")
