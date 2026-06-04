#!/usr/bin/env python3
"""
Three horizontal stacked-proportion bars:
  Bar A – co-opted gene count by TE Class (Class I Retrotransposons vs Class II DNA transposons)
  Bar B – co-opted gene count by TE domain (host_domain), all, top 20 + Other
  Bar C – non-helicase TE domains only, re-proportioned among themselves

Input : figC_domain_tepair_summary.tsv  (host_domain, te_superfam, te_category, n_genes, n_species)
Output: figC_TE_class_proportions.pdf / .png
"""

import os, re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"
OUT_DIR  = DATA_DIR
TSV      = os.path.join(DATA_DIR, "figC_domain_tepair_summary.tsv")

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_csv(TSV, sep="\t")

# ── Aggregate ─────────────────────────────────────────────────────────────────
# Bar A: by TE class
by_class = df.groupby("te_category")["n_genes"].sum()
class_order = ["DNA transposon", "Retrotransposon"]
by_class = by_class.reindex(class_order)

# Bar B+C: ALL TE structural domains (te_superfam) — 8 types, no top-N needed
# Rename te_superfam values to TE-domain style labels
TE_DOMAIN_RENAME = {
    "Helitron":              "Helicase/Helitron",
    "Transposase":           "Transposase (DDE)",
    "LINE":                  "APE endonuclease (LINE)",
    "LINE/LTR":              "Zn-finger/misc. (LINE/LTR)",
    "Reverse_transcriptase": "RVT (reverse transcriptase)",
    "Integrase":             "Integrase (RVE)",
    "Protease":              "Aspartyl protease",
    "Gag":                   "Gag / capsid",
}
df["te_domain"] = df["te_superfam"].map(TE_DOMAIN_RENAME).fillna(df["te_superfam"])

# Bar B: all TE domains (sorted by n_genes)
by_te_b = df.groupby("te_domain")["n_genes"].sum().sort_values(ascending=False)

# Bar C: non-helicase TE domains (exclude Helitron/Helicase), re-proportioned
NON_HEL_NAMES = [v for k,v in TE_DOMAIN_RENAME.items() if k != "Helitron"]
by_te_c = by_te_b.reindex([d for d in by_te_b.index if d in NON_HEL_NAMES]).dropna()
total_nh  = by_te_c.sum()
total_all = by_te_b.sum()
nh_pct    = total_nh / total_all * 100

# ── Colours ───────────────────────────────────────────────────────────────────
CLASS_COLORS = {
    "DNA transposon":  "#2980b9",
    "Retrotransposon": "#e67e22",
}

# TE domain colours: blues/teals for Class II (DNA), warm for Class I (Retro)
TE_DOM_COLORS = {
    "Helicase/Helitron":             "#1a5fa8",  # deep blue   (Class II)
    "Transposase (DDE)":             "#5dade2",  # light blue  (Class II)
    "APE endonuclease (LINE)":       "#e67e22",  # orange      (Class I)
    "Zn-finger/misc. (LINE/LTR)":   "#f39c12",  # amber       (Class I)
    "RVT (reverse transcriptase)":   "#e74c3c",  # red         (Class I)
    "Integrase (RVE)":               "#c0392b",  # dark red    (Class I)
    "Aspartyl protease":             "#8e44ad",  # purple      (Class I)
    "Gag / capsid":                  "#d35400",  # dark orange (Class I)
}

# ── Plot helper ────────────────────────────────────────────────────────────────
def _stacked_hbar(ax, counts, colors, title, label_min_pct=2.0):
    total = counts.sum()
    left  = 0.0
    for key, val in counts.items():
        pct = val / total * 100
        color = colors.get(key, "#aaaaaa")
        ax.barh(0, pct, left=left, height=0.55,
                color=color, edgecolor="white", linewidth=0.8)
        if pct >= label_min_pct:
            ax.text(left + pct / 2, 0, f"{pct:.1f}%",
                    ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        left += pct
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.6, 0.6)
    ax.set_yticks([])
    ax.set_xlabel("Proportion (%)", fontsize=9)
    ax.set_title(title, fontsize=9.5, fontweight="bold", loc="left", pad=4)
    ax.spines[["top", "left", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(11, 7),
                         gridspec_kw={"height_ratios": [1, 1.2, 1.2], "hspace": 1.6})

# Bar A — TE Class
_stacked_hbar(axes[0], by_class, CLASS_COLORS,
    "A   Co-opted genes by TE class", label_min_pct=1.5)
legend_a = [mpatches.Patch(
    facecolor=CLASS_COLORS[c],
    label=f"Class II  ({c})" if "DNA" in c else f"Class I  ({c})")
    for c in class_order]
axes[0].legend(handles=legend_a, loc="upper left",
               bbox_to_anchor=(0.0, -0.55), fontsize=8, framealpha=0.9, ncol=2)

# Bar B — ALL TE structural domains
_stacked_hbar(axes[1], by_te_b, TE_DOM_COLORS,
    "B   Co-opted genes by TE structural domain (all)", label_min_pct=1.0)
legend_b = [mpatches.Patch(
    facecolor=TE_DOM_COLORS.get(d, "#aaaaaa"),
    label=f"{d}  ({int(by_te_b[d]):,})") for d in by_te_b.index]
axes[1].legend(handles=legend_b, loc="upper left",
               bbox_to_anchor=(0.0, -0.55), fontsize=7.5, framealpha=0.9, ncol=2)

# Bar C — non-helicase TE structural domains, re-proportioned
title_c = (f"C   Non-helicase TE domains — re-proportioned "
           f"({nh_pct:.1f}% of all co-opted genes, n\u2009=\u2009{total_nh:,})")
_stacked_hbar(axes[2], by_te_c, TE_DOM_COLORS, title_c, label_min_pct=2.0)
legend_c = [mpatches.Patch(
    facecolor=TE_DOM_COLORS.get(d, "#aaaaaa"),
    label=f"{d}  ({int(by_te_c[d]):,})") for d in by_te_c.index]
axes[2].legend(handles=legend_c, loc="upper left",
               bbox_to_anchor=(0.0, -0.55), fontsize=7.5, framealpha=0.9, ncol=2)

fig.suptitle("TE origin of co-opted protein domains across Fungi",
             fontsize=11, fontweight="bold", y=1.01)

plt.tight_layout()
for ext in ("pdf", "png"):
    out = os.path.join(OUT_DIR, f"figC_TE_class_proportions.{ext}")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"[OK] Saved {out}")
plt.close()

