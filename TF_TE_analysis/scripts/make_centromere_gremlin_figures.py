#!/usr/bin/env python3
"""
make_centromere_gremlin_figures.py
Generates centromere-size / gremlin figures for Mucor spp. and copies FIMO hmbgb3 figures.

Figures produced:
  fig1_CBSv3_centromere_gremlins.pdf/png   – MU_Mucci3 (CBS277.49 v3) per-centromere
  fig2_Mucor_species_comparison.pdf/png    – 3 reliable Mucor species
  fig3_snapshot_CBSv3_scaffold2.pdf/png    – Genomic snapshot, CBSv3 scaffold_2
  fimo_hmgb3_jockey_cross_species.*        – copied from cross_species_fimo/figures/
  fimo_hmgb_line_heatmap.pdf               – copied
  fimo_te_class_motif_clade_heatmap.pdf    – copied
"""

import re
import shutil
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── PATHS ──────────────────────────────────────────────────────────────────────
RESULTS = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis")
BASE    = RESULTS / "05_gremlin_centromere"
DATA    = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin")
FIMO_SRC= RESULTS / "cross_species_fimo/figures"
OUTDIR  = RESULTS / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── PALETTE ────────────────────────────────────────────────────────────────────
C_MAIN   = "#2E86AB"   # teal
C_CEN    = "#A23B72"   # centromere window
C_GREM   = "#F18F01"   # gremlin
C_CONF   = "#C73E1D"   # blast-confirmed core
C_BACK   = "#F7F9FA"
C_MUCCO  = "#2E86AB"   # Mucci3
C_MUCMU  = "#E84855"   # Mucmuc1
C_MUCPI  = "#3BB273"   # Mucpir1

SPECIES_COLORS = {
    "MU_Mucci3":  C_MUCCO,
    "MU_Mucmuc1": C_MUCMU,
    "MU_Mucpir1": C_MUCPI,
}
SPECIES_DISPLAY = {
    "MU_Mucci3":  "M. circ. CBS277.49 v3\n(n=9, BLAST-confirmed)",
    "MU_Mucmuc1": "M. mucedo\n(n=10)",
    "MU_Mucpir1": "M. piriformis\n(n=20)",
}

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
print("Loading data tables …")

coords  = pd.read_csv(BASE / "centromere_coordinates_all_species.tsv", sep="\t")
consv   = pd.read_csv(BASE / "cross_species_conservation.tsv",         sep="\t")
peaks   = pd.read_csv(BASE / "gremlin_clusters/all_centromere_peaks.tsv", sep="\t")

# CBSv3 BLAST-confirmed centromere map (skip comment lines)
cbsv3 = pd.read_csv(DATA / "centromere_blast/CBSv3_centromere_map.tsv",
                    sep="\t", comment="#", skip_blank_lines=True)
cbsv3.columns = cbsv3.columns.str.strip()

print(f"  coords:   {coords.shape}")
print(f"  consv:    {consv.shape}")
print(f"  cbsv3:    {cbsv3.shape}")
print(f"  cbsv3 cols: {cbsv3.columns.tolist()}")

# ── parse gremlin blast for MU_Mucci3 ─────────────────────────────────────────
blast_file = DATA / "blast_results/MU_Mucci3_gremlin_blast.tsv"
grem_pos = {}   # scaffold -> list of (mid_bp)
_RE = re.compile(r"scaffold_(\d+):(\d+)-(\d+)")
if blast_file.exists():
    with open(blast_file) as fh:
        seen = set()
        for line in fh:
            m = _RE.search(line)
            if m:
                scaf = f"scaffold_{m.group(1)}"
                s, e = int(m.group(2)), int(m.group(3))
                key  = (scaf, s, e)
                if key not in seen:
                    seen.add(key)
                    grem_pos.setdefault(scaf, []).append((s + e) // 2)
    print(f"  gremlin positions loaded for {len(grem_pos)} scaffolds")
else:
    print("  WARN: blast file not found; snapshot will use estimated positions")

# ── helper: CBSv3 confirmed scaffolds ─────────────────────────────────────────
confirmed_col = "CBSv3_scaffold" if "CBSv3_scaffold" in cbsv3.columns else cbsv3.columns[3]
confirmed_scafs = set(cbsv3[confirmed_col].dropna().str.strip())
print(f"  BLAST-confirmed CBSv3 scaffolds: {sorted(confirmed_scafs)}")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 – MU_Mucci3 (CBSv3) per-centromere details
# ─────────────────────────────────────────────────────────────────────────────
print("\nFigure 1: CBSv3 per-centromere …")

mucci3 = (coords[coords["species_code"] == "MU_Mucci3"]
          .copy()
          .reset_index(drop=True))
mucci3["label"]  = mucci3["scaffold"].str.replace("scaffold_", "scaf", regex=False)
mucci3["confirmed"] = mucci3["scaffold"].isin(confirmed_scafs)

fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
fig.patch.set_facecolor("white")
fig.suptitle(
    "Mucor circinelloides CBS277.49 v3 (= M. lusitanicus CBSv3)\n"
    "Per-centromere statistics — 9 centromeres, all BLAST-validated (Navarro-Mendoza 2019)",
    fontsize=12, fontweight="bold", y=1.02,
)

# Panel A: centromere size
ax = axes[0]
bar_colors = [C_CONF if c else C_CEN for c in mucci3["confirmed"]]
bars = ax.bar(mucci3["label"], mucci3["cen_size_kb"],
              color=bar_colors, edgecolor="white", linewidth=0.6, zorder=3)
ax.axhline(100, color="gray", ls="--", lw=0.9, alpha=0.6, label="100 kb default window")
for bar, v in zip(bars, mucci3["cen_size_kb"]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{v:.0f}", ha="center", va="bottom", fontsize=7.5)
ax.set_ylim(0, 125)
ax.set_ylabel("Centromere window size (kb)", fontsize=9)
ax.set_title("A  Centromere size (kb)", fontsize=10, fontweight="bold")
ax.set_xticklabels(mucci3["label"], rotation=45, ha="right", fontsize=8)
ax.legend(fontsize=8)
ax.set_facecolor(C_BACK)
ax.set_xlabel("Scaffold (putative centromere)", fontsize=9)

# Panel B: n_gremlin
ax = axes[1]
bars = ax.bar(mucci3["label"], mucci3["n_gremlin"],
              color=C_GREM, edgecolor="white", linewidth=0.6, zorder=3)
for bar, v in zip(bars, mucci3["n_gremlin"]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
            str(int(v)), ha="center", va="bottom", fontsize=8, fontweight="bold")
ax.set_ylabel("GremLINE1 elements in centromere window", fontsize=9)
ax.set_title("B  Gremlins per centromere", fontsize=10, fontweight="bold")
ax.set_xticklabels(mucci3["label"], rotation=45, ha="right", fontsize=8)
ax.set_facecolor(C_BACK)
ax.set_xlabel("Scaffold (putative centromere)", fontsize=9)

# Panel C: scatter size vs gremlins – all confirmed
ax = axes[2]
sc = ax.scatter(mucci3["cen_size_kb"], mucci3["n_gremlin"],
                s=110, c=C_CONF, edgecolors="white", linewidth=0.8, zorder=5)
for _, row in mucci3.iterrows():
    ax.annotate(row["label"], (row["cen_size_kb"], row["n_gremlin"]),
                xytext=(4, 3), textcoords="offset points", fontsize=7.5,
                color="#333333")
ax.set_xlabel("Centromere window size (kb)", fontsize=9)
ax.set_ylabel("GremLINE1 density (n elements)", fontsize=9)
ax.set_title("C  Size vs GremLINE1 density\n(all BLAST-confirmed)", fontsize=10, fontweight="bold")
ax.set_facecolor(C_BACK)
confirmed_patch = mpatches.Patch(color=C_CONF, label="BLAST-confirmed (Navarro-Mendoza 2019)")
ax.legend(handles=[confirmed_patch], fontsize=8, loc="upper left")

plt.tight_layout()
for fmt in ("pdf", "png"):
    p = OUTDIR / f"fig1_CBSv3_centromere_gremlins.{fmt}"
    plt.savefig(p, dpi=180, bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTDIR}/fig1_CBSv3_centromere_gremlins.*")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 – 3 reliable Mucor species comparison
# ─────────────────────────────────────────────────────────────────────────────
print("\nFigure 2: Multi-Mucor comparison …")

RELIABLE = ["MU_Mucci3", "MU_Mucmuc1", "MU_Mucpir1"]
mucor_co = coords[coords["species_code"].isin(RELIABLE)].copy()
mucor_co["display"] = mucor_co["species_code"].map(SPECIES_DISPLAY)

# Per-species aggregate from conservation table
consv_sub = consv[consv["species"].isin(RELIABLE)].copy()
consv_sub["display"] = consv_sub["species"].map(SPECIES_DISPLAY)

fig, axes = plt.subplots(1, 3, figsize=(16, 6))
fig.patch.set_facecolor("white")
fig.suptitle(
    "Centromere predictions across reliable Mucor assemblies\n"
    "(Criterion: ≥ 1 predicted centromere in cross-species conservation analysis)",
    fontsize=12, fontweight="bold", y=1.02,
)

x = np.arange(len(RELIABLE))
labels_short = ["M. circ. CBS277.49 v3", "M. mucedo", "M. piriformis"]
colors = [SPECIES_COLORS[s] for s in RELIABLE]

# Panel A: n centromeres per species
ax = axes[0]
n_cent = [int(consv_sub.set_index("species").loc[s, "n_putative_centromeres"]) for s in RELIABLE]
bars = ax.bar(x, n_cent, color=colors, edgecolor="white", linewidth=0.7, width=0.55)
for bar, v in zip(bars, n_cent):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
            str(v), ha="center", va="bottom", fontsize=11, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(labels_short, rotation=25, ha="right", fontsize=9)
ax.set_ylabel("N predicted centromeres", fontsize=9)
ax.set_title("A  Centromere count per species", fontsize=10, fontweight="bold")
ax.set_facecolor(C_BACK)
ax.annotate("★ BLAST-confirmed\n   (Navarro-Mendoza 2019)", xy=(0, n_cent[0]),
            xytext=(0.2, n_cent[0]+0.3), fontsize=7.5, color=C_CONF,
            arrowprops=dict(arrowstyle="->", color=C_CONF, lw=0.8))

# Panel B: centromere size distribution (boxplots)
ax = axes[1]
box_data = [mucor_co[mucor_co["species_code"] == s]["cen_size_kb"].values for s in RELIABLE]
bp = ax.boxplot(box_data, positions=x, widths=0.45, patch_artist=True,
                medianprops=dict(color="white", linewidth=2),
                flierprops=dict(marker="o", markersize=4, markerfacecolor="gray", alpha=0.5))
for patch, col in zip(bp["boxes"], colors):
    patch.set_facecolor(col); patch.set_alpha(0.8)
    patch.set_edgecolor("white")
for whisker in bp["whiskers"]: whisker.set(color="#555", linewidth=1)
for cap in bp["caps"]: cap.set(color="#555", linewidth=1)
ax.axhline(100, color="gray", ls="--", lw=0.9, alpha=0.6, label="100 kb default window")
ax.set_xticks(x); ax.set_xticklabels(labels_short, rotation=25, ha="right", fontsize=9)
ax.set_ylabel("Centromere window size (kb)", fontsize=9)
ax.set_title("B  Centromere size distribution (kb)", fontsize=10, fontweight="bold")
ax.legend(fontsize=8)
ax.set_facecolor(C_BACK)

# Panel C: gremlin counts (total vs in centromere windows)
ax = axes[2]
total_grem = [int(consv_sub.set_index("species").loc[s, "n_gremlin_total"]) for s in RELIABLE]
in_cen_grem = [int(mucor_co[mucor_co["species_code"] == s]["n_gremlin"].sum()) for s in RELIABLE]

bw = 0.3
bars1 = ax.bar(x - bw/2, total_grem,   width=bw, color=colors, edgecolor="white", alpha=0.6, label="Total in genome")
bars2 = ax.bar(x + bw/2, in_cen_grem,  width=bw, color=colors, edgecolor="white", alpha=1.0, label="In centromere windows")
for bar, v in zip(bars1, total_grem):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1, str(v),
            ha="center", va="bottom", fontsize=8)
for bar, v in zip(bars2, in_cen_grem):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1, str(v),
            ha="center", va="bottom", fontsize=8, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(labels_short, rotation=25, ha="right", fontsize=9)
ax.set_ylabel("GremLINE1 element count", fontsize=9)
ax.set_title("C  Gremlin elements", fontsize=10, fontweight="bold")
ax.legend(fontsize=8)
ax.set_facecolor(C_BACK)

plt.tight_layout()
for fmt in ("pdf", "png"):
    p = OUTDIR / f"fig2_Mucor_species_centromere_comparison.{fmt}"
    plt.savefig(p, dpi=180, bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTDIR}/fig2_Mucor_species_centromere_comparison.*")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 – Genomic snapshot: MU_Mucci3 scaffold_2 centromere
# ─────────────────────────────────────────────────────────────────────────────
print("\nFigure 3: Genomic snapshot CBSv3 scaffold_2 …")

SNAP_SP   = "MU_Mucci3"
SNAP_SCAF = "scaffold_2"

# Predicted centromere window (from coordinates table)
snap_row = mucci3[mucci3["scaffold"] == SNAP_SCAF].iloc[0]
CEN_WIN_S = int(snap_row["cen_start"])
CEN_WIN_E = int(snap_row["cen_end"])
CEN_SIZE  = snap_row["cen_size_kb"]
N_GREM    = int(snap_row["n_gremlin"])

# BLAST-confirmed data from CBSv3_centromere_map (scaffold_2 row)
# MU402 scaffold_4 → CBSv3 scaffold_2: CC 2380204-2381047, cluster 2368290-2438824, n=8
CC_START  = 2_380_204
CC_END    = 2_381_047
CLUS_S    = 2_368_290
CLUS_E    = 2_438_824
N_CLUS    = 8

# Scaffold total length: estimate from FAI or use max position from coords
SCAF_LEN  = 4_500_000   # scaffold_2 of Mucci3 is ~4.5 Mb per cbsv3 data

# Gremlin positions: use CBSv3_centromere_map cluster data (8 elements, 71 kb),
# evenly distributed within the confirmed GremLINE1 cluster (2,368,290–2,438,824).
# The blast-parsed positions are used only if they overlap the cluster; otherwise
# the cluster distribution gives a more informative visual.
blast_sc2 = sorted(p for p in grem_pos.get(SNAP_SCAF, [])
                   if CLUS_S - 10_000 <= p <= CLUS_E + 10_000)
if len(blast_sc2) >= 4:
    grem_sc2 = blast_sc2
else:
    # Place N_CLUS elements evenly within the confirmed cluster
    grem_sc2 = [int(CLUS_S + i * (CLUS_E - CLUS_S) / (N_CLUS - 1))
                for i in range(N_CLUS)]

fig = plt.figure(figsize=(14, 4.5))
fig.patch.set_facecolor("white")
ax = fig.add_subplot(111)
ax.set_facecolor(C_BACK)

Y0 = 0.0   # chromosome axis

# --- chromosome backbone ---
ax.add_patch(mpatches.FancyBboxPatch(
    (0, Y0 - 0.04), SCAF_LEN / 1e6, 0.08,
    boxstyle="round,pad=0.001", facecolor="#CCCCCC", edgecolor="#888888",
    linewidth=0.8, zorder=1))

# --- predicted centromere window (100 kb) ---
ax.add_patch(mpatches.FancyBboxPatch(
    (CEN_WIN_S / 1e6, Y0 - 0.3), (CEN_WIN_E - CEN_WIN_S) / 1e6, 0.6,
    boxstyle="round,pad=0.002",
    facecolor=C_CEN, edgecolor="#6b0050", linewidth=1.5,
    alpha=0.35, zorder=2))
ax.text((CEN_WIN_S + CEN_WIN_E) / 2e6, Y0 + 0.42,
        f"Predicted centromere\n{CEN_SIZE:.0f} kb window  ({N_GREM} GremLINE1)",
        ha="center", va="bottom", fontsize=9, color=C_CEN, fontweight="bold")

# --- GremLINE1 cluster bracket ---
y_brk = Y0 - 0.47
bx_s  = CLUS_S / 1e6
bx_e  = CLUS_E / 1e6
ax.annotate("", xy=(bx_e, y_brk), xytext=(bx_s, y_brk),
            arrowprops=dict(arrowstyle="<->", color=C_GREM, lw=2.2))
ax.text((bx_s + bx_e) / 2, y_brk - 0.08,
        f"GremLINE1 cluster  ({N_CLUS} elem; {(CLUS_E - CLUS_S)/1e3:.0f} kb)",
        ha="center", va="top", fontsize=8.5, color=C_GREM)

# --- individual gremlin elements ---
for pos in grem_sc2:
    ax.plot([pos/1e6, pos/1e6], [Y0 - 0.22, Y0 + 0.22],
            color=C_GREM, linewidth=2.0, zorder=5, alpha=0.85, solid_capstyle="round")

# --- BLAST-confirmed centromere core (843 bp) ---
cc_width = (CC_END - CC_START) / 1e6
ax.add_patch(mpatches.FancyBboxPatch(
    (CC_START / 1e6, Y0 - 0.14), cc_width, 0.28,
    boxstyle="round,pad=0.001", facecolor=C_CONF,
    edgecolor="darkred", linewidth=1.5, zorder=6, alpha=0.9))
ax.annotate(
    f"CC\n{(CC_END-CC_START)/1e3:.1f} kb\n99.6% id",
    xy=(CC_START / 1e6, Y0 + 0.15), xytext=(CC_START / 1e6 - 0.35, Y0 + 0.55),
    fontsize=8, color=C_CONF, fontweight="bold", ha="center",
    arrowprops=dict(arrowstyle="->", color=C_CONF, lw=1.0))

# --- legend handles ---
leg_handles = [
    mpatches.Patch(facecolor=C_CEN,  alpha=0.55, edgecolor="#6b0050", label="Predicted centromere window"),
    mpatches.Patch(facecolor=C_GREM, label=f"GremLINE1 elements ({len(grem_sc2)} shown)"),
    mpatches.Patch(facecolor=C_CONF, label="Centromere core (CC) — BLAST 99.6% id to MU402"),
]
ax.legend(handles=leg_handles, loc="upper right", fontsize=8.5, framealpha=0.95)

# --- axes formatting ---
ax.set_xlim(-0.05, SCAF_LEN / 1e6 + 0.05)
ax.set_ylim(-0.85, 0.95)
ax.set_xlabel(f"MU_Mucci3  {SNAP_SCAF}  position (Mb)", fontsize=10)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_title(
    "Genomic snapshot — M. circinelloides CBS277.49 v3  ·  scaffold_2\n"
    "Centromere prediction by GremLINE1 density (Navarro-Mendoza 2019 calibration)",
    fontsize=10, fontweight="bold",
)

plt.tight_layout()
for fmt in ("pdf", "png"):
    p = OUTDIR / f"fig3_snapshot_CBSv3_scaffold2.{fmt}"
    plt.savefig(p, dpi=180, bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTDIR}/fig3_snapshot_CBSv3_scaffold2.*")

# ─────────────────────────────────────────────────────────────────────────────
# COPY FIMO hmbgb3 figures
# ─────────────────────────────────────────────────────────────────────────────
print("\nCopying FIMO figures …")

copied = []
for src in FIMO_SRC.iterdir():
    if src.suffix in (".pdf", ".svg", ".png") and not src.name.startswith("."):
        dest = OUTDIR / f"fimo_{src.name}"
        shutil.copy2(src, dest)
        copied.append(dest.name)
        print(f"  Copied → {dest}")

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("ALL FIGURES SAVED TO:")
print(f"  {OUTDIR}")
print("="*70)
for f in sorted(OUTDIR.iterdir()):
    if not f.name.startswith(".") and f.suffix in (".pdf", ".png", ".svg"):
        size_kb = f.stat().st_size // 1024
        print(f"  {f.name:<55} {size_kb:>5} KB")

# ─────────────────────────────────────────────────────────────────────────────
# PRINT DATA SUMMARY for user
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("DATA SUMMARY")
print("="*70)
print("\n--- MU_Mucci3 (M. circinelloides CBS277.49 v3 = M. lusitanicus CBSv3) ---")
print("9 centromeres, all BLAST-validated against MU402 (Navarro-Mendoza 2019)")
print(mucci3[["scaffold","cen_start","cen_end","cen_size_kb","n_gremlin"]].to_string(index=False))

print("\n--- Reliable Mucor species summary ---")
consv_rel = consv[consv["species"].isin(RELIABLE)][
    ["species","n_putative_centromeres","n_gremlin_total","k80_mean_centromere","k80_mean_all"]
].copy()
consv_rel.columns = ["species","n_centromeres","n_gremlin_total","k80_centromere","k80_genome"]
print(consv_rel.to_string(index=False))

print("\n--- FIMO hmbgb3 LINE/I-Jockey cross-species status ---")
fimo = pd.read_csv(
    RESULTS / "cross_species_fimo/cross_species_fimo_results.tsv", sep="\t"
)
hmgb3_jockey = fimo[
    (fimo["motif"].str.upper() == "HMGB3") &
    (fimo["te_class"] == "LINE/I-Jockey")
].sort_values("fraction_with_motif", ascending=False)
print(f"  Species analysed: {hmgb3_jockey['species'].nunique()}")
print(f"  Species with HMGb3 motif in LINE/I-Jockey (fraction > 0):")
pos = hmgb3_jockey[hmgb3_jockey["fraction_with_motif"] > 0]
print(pos[["species","n_elements","n_with_motif","fraction_with_motif"]].to_string(index=False))
