#!/usr/bin/env python3
"""
Single-row heatmap: proportion of LINE/I-Jockey copies with HMGb3 motif,
per fungal clade, across the cross-species FIMO dataset.
Columns = clades (phylogenetic order), color = mean fraction (0-1).
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
CS_FIMO = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/cross_species_fimo/cross_species_fimo_results.tsv")
OUT_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────────────────────
MIN_N     = 5        # minimum elements per species to include
TE_CLASS  = "LINE/I-Jockey"
MOTIF     = "HMGb3"

# Phylogenetic order and display labels
# (excludes HO = non-fungal outgroup)
CLADE_ORDER = [
    "NU",   # Nucleariida (outgroup)
    "AP",   # Aphelida
    "CR",   # Cryptomycota
    "BL",   # Blastocladiomycota
    "CH",   # Chytridiomycota
    "OL",   # Olpidiomycota
    "EN",   # Entomophthoromycotina (Zoopagomycota)
    "ZO",   # Zoopagomycota
    "KX",   # Kickxellomycotina (Zoopagomycota)
    "GL",   # Glomeromycota (Mucoromycota)
    "MO",   # Mortierellomycotina
    "MU",   # Mucoromycotina
    "BA",   # Basidiomycota
    "AS",   # Ascomycota
]

CLADE_LABELS = {
    "NU":  "Nucleariida",
    "AP":  "Aphelida",
    "CR":  "Cryptomycota",
    "BL":  "Blastocladio.",
    "CH":  "Chytrid.",
    "OL":  "Olpidio.",
    "EN":  "Entomophtho.",
    "ZO":  "Zoopago.",
    "KX":  "Kickxello.",
    "GL":  "Glomero.",
    "MO":  "Mortierello.",
    "MU":  "Mucoro.",
    "BA":  "Basidio.",
    "AS":  "Ascomycota",
}

# ── Load and filter data ───────────────────────────────────────────────────────
cs = pd.read_csv(CS_FIMO, sep="\t")
cs["clade"] = cs["species"].str.split("_").str[0]

# Filter to target TE × motif, minimum n per species
lij = cs[
    (cs["te_class"] == TE_CLASS) &
    (cs["motif"]    == MOTIF) &
    (cs["n_elements"] >= MIN_N)
].copy()

# Per-clade aggregation
clade_stats = (
    lij.groupby("clade")
    .agg(
        n_species           = ("species",              "nunique"),
        mean_frac           = ("fraction_with_motif",  "mean"),
        median_frac         = ("fraction_with_motif",  "median"),
        se_frac             = ("fraction_with_motif",  lambda x: x.std(ddof=1) / np.sqrt(len(x))),
        n_elements_total    = ("n_elements",           "sum"),
    )
    .reset_index()
)

# Keep only clades present in data and in our order list
present = [c for c in CLADE_ORDER if c in clade_stats["clade"].values]
clade_stats = clade_stats.set_index("clade").loc[present].reset_index()

print("Clades included in heatmap:")
print(clade_stats.to_string(index=False))

# ── Plot ───────────────────────────────────────────────────────────────────────
n_clades = len(clade_stats)
fig_width = max(10, n_clades * 0.85)

fig, ax = plt.subplots(figsize=(fig_width, 2.8))

cmap = plt.cm.YlOrRd
norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

fracs   = clade_stats["mean_frac"].values
labels  = [CLADE_LABELS.get(c, c) for c in clade_stats["clade"]]
n_sp    = clade_stats["n_species"].values

# Draw rectangles (single row)
rect_h = 1.0
for i, (frac, n) in enumerate(zip(fracs, n_sp)):
    color = cmap(norm(frac))
    rect  = plt.Rectangle([i, 0], 1, rect_h, color=color, lw=0.5, ec="white")
    ax.add_patch(rect)

    # Annotate: fraction value
    ax.text(i + 0.5, 0.62, f"{frac:.2f}",
            ha="center", va="center", fontsize=8.5, fontweight="bold",
            color="white" if frac > 0.55 else "black")
    # n species below
    ax.text(i + 0.5, 0.28, f"n={n}",
            ha="center", va="center", fontsize=7.0,
            color="white" if frac > 0.55 else "#444444")

ax.set_xlim(0, n_clades)
ax.set_ylim(0, rect_h)

# X-axis labels (clade names)
ax.set_xticks([i + 0.5 for i in range(n_clades)])
ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=9)
ax.xaxis.tick_top()
ax.xaxis.set_label_position("top")
ax.set_xticklabels(labels, rotation=40, ha="left", fontsize=9)

# Remove y axis
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation="horizontal",
                    fraction=0.04, pad=0.22, aspect=30, shrink=0.5)
cbar.set_label("Mean fraction of LINE/I-Jockey copies\nwith HMGb3 motif", fontsize=9)
cbar.ax.tick_params(labelsize=8)

ax.set_title(f"Conservation of HMGb3 binding motif in LINE/I-Jockey (Grem-LINE1) across fungi\n"
             f"(species with ≥{MIN_N} full-length copies)",
             fontsize=10, pad=18)

plt.tight_layout()
out_pdf = OUT_DIR / "hmgb3_jockey_clade_heatmap.pdf"
out_png = OUT_DIR / "hmgb3_jockey_clade_heatmap.png"
fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
fig.savefig(out_png, bbox_inches="tight", dpi=200)
plt.close()
print(f"\nSaved: {out_pdf}")
print(f"Saved: {out_png}")
