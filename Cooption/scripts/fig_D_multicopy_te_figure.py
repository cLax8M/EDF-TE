#!/usr/bin/env python3
"""
fig_D_multicopy_te_figure.py
Improved visualization of TE proximity around multicopy vs single-copy genes.

Uses:
  - 16_multicopy_te_summary.tsv: per-species stats
  - clade_assignment.tsv: clade info
  - per_species/classification/*.true_TE_coopted_genes.tsv: co-option status

Figures produced:
  figD_multicopy_ratio_by_clade.pdf/.png
    Violin/boxplot per clade: ratio (median_multi / median_single) of TE bp
    with dots = species, colored by significance

  figD_mule_enrichment.pdf/.png
    % genes with MULE nearby: multi vs single, across clades
"""
import sys, os, glob, collections
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu

# ─── Paths ──────────────────────────────────────────────────────────────────
SUMMARY_FILE  = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures/16_multicopy_te_summary.tsv"
CLADE_FILE    = "/legserv/Temp/Carlos/EDF-TE/results/HTT/clade_assignment.tsv"
OUT_DIR       = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"

CLADE_LABELS = {
    "AP":"Aphelida","AS":"Ascomycota","BA":"Basidiomycota",
    "BL":"Blastocladiomycota","BS":"Basidiobolomycota","CH":"Chytridiomycota",
    "CR":"Cryptomycota","EN":"Entomophthoromycota","GL":"Glomeromycota",
    "HO":"Holozoa","KX":"Kickxellomycotina","MI":"Microsporidia",
    "MO":"Mortierellomycota","MU":"Mucoromycotina","NU":"Nuclearida",
    "OL":"Olpidiomycota","ZO":"Zoopagomycota",
}

CLADE_COLORS = {
    "HO":"#e07b54","MI":"#b0b0b0","NU":"#d4b8e0","AP":"#5ba3d9",
    "GL":"#2e7d32","MO":"#81c784","MU":"#c8e6c9","ZO":"#ff8a65",
    "EN":"#ffb74d","BS":"#f06292","KX":"#ba68c8","BL":"#64b5f6",
    "CH":"#4dd0e1","CR":"#fff176","OL":"#a5d6a7","AS":"#ff7043",
    "BA":"#8d6e63",
}

# ─── Load data ───────────────────────────────────────────────────────────────
df = pd.read_csv(SUMMARY_FILE, sep="\t")
clades = pd.read_csv(CLADE_FILE, sep="\t", header=None, names=["species","clade"])
df = df.merge(clades, on="species", how="left")
df["clade"] = df["clade"].fillna(df["species"].str.split("_").str[0])
df["clade_label"] = df["clade"].map(CLADE_LABELS).fillna(df["clade"])

# Ratio: median TE bp in multicopy / single-copy flanking regions
# Avoid division by zero
df["ratio"] = np.where(
    df["median_single"] > 0,
    df["median_multi"] / df["median_single"],
    np.nan
)
df["log2_ratio"] = np.log2(df["ratio"].clip(lower=1e-3))

print(f"[INFO] Loaded {len(df)} species with TE proximity data")
print(f"[INFO] Significant (p<0.05) multicopy enrichment: "
      f"{df['signif'].sum()} of {len(df)} species")

# Order clades by median ratio, filter ≥3 species
clade_n = df["clade"].value_counts()
df = df[df["clade"].isin(clade_n[clade_n >= 3].index)].copy()

clade_meds  = df.groupby("clade_label")["log2_ratio"].median().sort_values(ascending=False)
clade_order = clade_meds.index.tolist()

# ─── Figure 1: log2 ratio per clade ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 6))

for i, cl in enumerate(clade_order):
    sub = df[df["clade_label"] == cl]
    code = sub["clade"].iloc[0] if len(sub) > 0 else "AS"
    color = CLADE_COLORS.get(code, "#888888")

    # Violin (if ≥5 points) else boxplot
    vals = sub["log2_ratio"].dropna().values
    if len(vals) >= 4:
        vp = ax.violinplot([vals], positions=[i], widths=0.7, showmedians=True,
                           showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor(color)
            body.set_alpha(0.6)
            body.set_edgecolor("grey")
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.5)

    # Individual points (colored by significance)
    sig  = sub["signif"].fillna(False)
    ax.scatter([i - 0.04]*len(sub[sig]),   sub.loc[sig, "log2_ratio"],
               color="black", s=12, zorder=5, alpha=0.8)
    ax.scatter([i + 0.04]*len(sub[~sig]),  sub.loc[~sig, "log2_ratio"],
               color="grey",  s=12, zorder=5, alpha=0.5)

ax.axhline(0, color="black", lw=0.8, ls="--", label="Ratio = 1 (no difference)")
ax.set_xticks(range(len(clade_order)))
ax.set_xticklabels(clade_order, rotation=40, ha="right", fontsize=9)
ax.set_ylabel("log₂(median TE bp multicopy / single-copy) ±5 kb", fontsize=10)
ax.set_title("TE content near multicopy vs single-copy gene flanks, by clade\n"
             "(colored = p<0.05 Mann-Whitney; grey = not significant)", fontsize=11)

sig_patch = mpatches.Patch(color="black", label="p < 0.05")
ns_patch  = mpatches.Patch(color="grey",  label="n.s.")
ax.legend(handles=[sig_patch, ns_patch], fontsize=9)
plt.tight_layout()

fig.savefig(f"{OUT_DIR}/figD_multicopy_ratio_by_clade.pdf", dpi=150, bbox_inches="tight")
fig.savefig(f"{OUT_DIR}/figD_multicopy_ratio_by_clade.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"[OK] figD_multicopy_ratio_by_clade.pdf/.png")

# ─── Figure 2: MULE enrichment multi vs single by clade ──────────────────────
df["mule_diff"] = df["pct_mule_multi"] - df["pct_mule_single"]
clade_mule_med = df.groupby("clade_label")["mule_diff"].median().sort_values(ascending=False)
clade_order2 = clade_mule_med.index.tolist()

fig2, ax2 = plt.subplots(figsize=(12, 5))
colors2 = []
for cl in clade_order2:
    code = df[df["clade_label"] == cl]["clade"].iloc[0]
    colors2.append(CLADE_COLORS.get(code, "#888888"))

medians2 = [clade_mule_med[cl] for cl in clade_order2]
bars = ax2.bar(range(len(clade_order2)), medians2, color=colors2, edgecolor="grey", lw=0.5)
ax2.axhline(0, color="black", lw=0.8, ls="--")
ax2.set_xticks(range(len(clade_order2)))
ax2.set_xticklabels(clade_order2, rotation=40, ha="right", fontsize=9)
ax2.set_ylabel("Median % MULE-associated (multicopy − single-copy)", fontsize=10)
ax2.set_title("MULE proximity enrichment: multicopy vs single-copy genes\n"
              "(positive = more MULE near multicopy)", fontsize=11)
plt.tight_layout()
fig2.savefig(f"{OUT_DIR}/figD_mule_enrichment.pdf", dpi=150, bbox_inches="tight")
fig2.savefig(f"{OUT_DIR}/figD_mule_enrichment.png", dpi=150, bbox_inches="tight")
plt.close(fig2)
print(f"[OK] figD_mule_enrichment.pdf/.png")

# ─── Figure 3: summary scatter — ratio vs %co-option per species ─────────────
TRAITS_FILE = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures/14_cooption_traits.tsv"
try:
    traits = pd.read_csv(TRAITS_FILE, sep="\t")
    df_merge = df.merge(traits[["species","pct_coopted"]], on="species", how="left")

    fig3, ax3 = plt.subplots(figsize=(7, 6))
    for code, grp in df_merge.groupby("clade"):
        color = CLADE_COLORS.get(code, "#888888")
        label = CLADE_LABELS.get(code, code)
        ax3.scatter(grp["pct_coopted"], grp["log2_ratio"],
                    color=color, s=20, alpha=0.75, label=label)
    ax3.axhline(0, color="grey", lw=0.8, ls="--")
    ax3.set_xlabel("% co-opted genes per species", fontsize=10)
    ax3.set_ylabel("log₂ TE bp ratio (multicopy / single-copy)", fontsize=10)
    ax3.set_title("TE proximity ratio vs co-option level", fontsize=11)
    ax3.legend(fontsize=7, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.tight_layout()
    fig3.savefig(f"{OUT_DIR}/figD_ratio_vs_pct_coopted.pdf", dpi=150, bbox_inches="tight")
    fig3.savefig(f"{OUT_DIR}/figD_ratio_vs_pct_coopted.png", dpi=150, bbox_inches="tight")
    plt.close(fig3)
    print(f"[OK] figD_ratio_vs_pct_coopted.pdf/.png")
except Exception as e:
    print(f"[SKIP] scatter vs cooption: {e}")

print("[DONE] Fig D complete.")
