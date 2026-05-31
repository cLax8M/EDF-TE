#!/usr/bin/env python3
"""
15_clade_cooption_analysis.py
Clade-level comparison of TE co-option burden.

For each fungal/host clade (defined by 2-letter species prefix) computes:
  - Distribution of pct_coopted (median, IQR, min/max)
  - Kruskal-Wallis H-test across all clades (non-parametric ANOVA)
  - Pairwise Dunn post-hoc tests with Benjamini-Hochberg FDR correction
  - TE superfamily prevalence per clade (dominant_te_superfam frequencies)
  - Top host-domain usage per clade (dominant_host_domain frequencies)

Inputs:
  - figures/14_cooption_traits.tsv          (per-species co-option table)

Outputs:
  - figures/15_clade_stats.tsv             (per-clade summary statistics)
  - figures/15_dunn_pairwise.tsv           (pairwise Dunn test results)
  - figures/15_clade_boxplot.pdf           (pct_coopted boxplots per clade)
  - figures/15_te_superfam_heatmap.pdf     (TE superfamily × clade heatmap)
  - figures/15_host_domain_heatmap.pdf     (host domain × clade heatmap)

Usage:
    python3 15_clade_cooption_analysis.py
"""

import sys
import os
import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LogNorm
import scipy.stats as ss
from itertools import combinations
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = "/legserv/Temp/Carlos/EDF-TE/results/Co-option"
TRAITS_FILE = f"{BASE}/figures/14_cooption_traits.tsv"
OUT_DIR = f"{BASE}/figures"
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Full clade names for display
# ---------------------------------------------------------------------------
CLADE_LABELS = {
    "AP": "Aphelida",
    "AS": "Ascomycota",
    "BA": "Basidiomycota",
    "BL": "Blastocladiomycota",
    "BS": "Basidiobolomycota",
    "CH": "Chytridiomycota",
    "CR": "Cryptomycota",
    "EN": "Entomophthoromycota",
    "GL": "Glomeromycota",
    "HO": "Holozoa",
    "KX": "Kickxellomycotina",
    "MI": "Microsporidia",
    "MO": "Mortierellomycota",
    "MU": "Mucoromycotina",
    "NU": "Nuclearida",
    "OL": "Olpidiomycota",
    "ZO": "Zoopagomycota",
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df = pd.read_csv(TRAITS_FILE, sep="\t")
df["clade"] = df["species"].str.split("_").str[0]
df["clade_label"] = df["clade"].map(CLADE_LABELS).fillna(df["clade"])

# drop clades with <3 species (insufficient for stats)
clade_sizes = df["clade"].value_counts()
valid_clades = clade_sizes[clade_sizes >= 3].index.tolist()
df_valid = df[df["clade"].isin(valid_clades)].copy()

print(f"[INFO] Species: {len(df)} total, {len(df_valid)} in clades with ≥3 species")
print(f"[INFO] Clades analysed: {sorted(valid_clades)}")

# ---------------------------------------------------------------------------
# 1. Per-clade summary statistics
# ---------------------------------------------------------------------------
stats_rows = []
for clade, grp in df_valid.groupby("clade"):
    vals = grp["pct_coopted"].dropna()
    stats_rows.append({
        "clade":       clade,
        "clade_label": CLADE_LABELS.get(clade, clade),
        "n_species":   len(vals),
        "mean_pct":    round(vals.mean(), 4),
        "median_pct":  round(vals.median(), 4),
        "sd_pct":      round(vals.std(), 4),
        "q25_pct":     round(vals.quantile(0.25), 4),
        "q75_pct":     round(vals.quantile(0.75), 4),
        "min_pct":     round(vals.min(), 4),
        "max_pct":     round(vals.max(), 4),
        "mean_n_coopted":   round(grp["n_coopted"].mean(), 1),
        "median_n_coopted": round(grp["n_coopted"].median(), 1),
    })

stats_df = pd.DataFrame(stats_rows).sort_values("median_pct", ascending=False)
stats_df.to_csv(f"{OUT_DIR}/15_clade_stats.tsv", sep="\t", index=False)
print(f"[OK] Clade stats → {OUT_DIR}/15_clade_stats.tsv")
print(stats_df[["clade_label", "n_species", "median_pct", "sd_pct"]].to_string(index=False))

# ---------------------------------------------------------------------------
# 2. Kruskal-Wallis test (non-parametric, no normality assumption)
# ---------------------------------------------------------------------------
groups = [grp["pct_coopted"].dropna().values
          for _, grp in df_valid.groupby("clade")]
kw_stat, kw_p = ss.kruskal(*groups)
print(f"\n[STATS] Kruskal-Wallis H={kw_stat:.3f}, p={kw_p:.2e}, df={len(groups)-1}")

# ---------------------------------------------------------------------------
# 3. Pairwise Dunn post-hoc test with BH FDR correction
# ---------------------------------------------------------------------------
def dunn_posthoc(df_in, group_col, val_col):
    """Dunn's test: z-score-based pairwise comparison after KW."""
    groups_ord = sorted(df_in[group_col].unique())
    all_vals = df_in[val_col].dropna().values
    n_total = len(all_vals)
    # ranks on the combined sample
    ranked = ss.rankdata(df_in[val_col].values)
    df_in = df_in.copy()
    df_in["rank"] = ranked

    rows = []
    for g1, g2 in combinations(groups_ord, 2):
        r1 = df_in.loc[df_in[group_col] == g1, "rank"].values
        r2 = df_in.loc[df_in[group_col] == g2, "rank"].values
        n1, n2 = len(r1), len(r2)
        mean_r1 = r1.mean()
        mean_r2 = r2.mean()
        # standard error
        se = np.sqrt((n_total * (n_total + 1) / 12.0) * (1/n1 + 1/n2))
        z = (mean_r1 - mean_r2) / se
        p = 2 * ss.norm.sf(abs(z))
        rows.append({"clade1": g1, "clade2": g2, "z": round(z, 3), "p_raw": p})

    result = pd.DataFrame(rows)
    _, bh_p, _, _ = multipletests(result["p_raw"].values, method="fdr_bh")
    result["p_adj_BH"] = bh_p
    result["significant"] = result["p_adj_BH"] < 0.05
    return result.sort_values("p_adj_BH")

dunn_df = dunn_posthoc(df_valid, "clade", "pct_coopted")
# Add full labels
dunn_df["label1"] = dunn_df["clade1"].map(CLADE_LABELS).fillna(dunn_df["clade1"])
dunn_df["label2"] = dunn_df["clade2"].map(CLADE_LABELS).fillna(dunn_df["clade2"])
dunn_df.to_csv(f"{OUT_DIR}/15_dunn_pairwise.tsv", sep="\t", index=False)
n_sig = dunn_df["significant"].sum()
print(f"[OK] Dunn pairwise → {OUT_DIR}/15_dunn_pairwise.tsv ({n_sig} significant pairs)")
print(dunn_df[dunn_df["significant"]][["label1", "label2", "z", "p_adj_BH"]].to_string(index=False))

# ---------------------------------------------------------------------------
# 4. Boxplot — pct_coopted per clade (sorted by median)
# ---------------------------------------------------------------------------
clade_order = stats_df["clade"].tolist()   # already sorted by median desc
clade_labels_ord = [CLADE_LABELS.get(c, c) for c in clade_order]
plot_data = [df_valid.loc[df_valid["clade"] == c, "pct_coopted"].dropna().values
             for c in clade_order]

fig, ax = plt.subplots(figsize=(14, 6))
bp = ax.boxplot(plot_data, patch_artist=True, notch=False,
                medianprops=dict(color="black", linewidth=2),
                flierprops=dict(marker="o", markersize=3, alpha=0.5))

# Colour by median value (warm → cool gradient)
medians = [np.median(d) if len(d) > 0 else 0 for d in plot_data]
cmap = plt.cm.RdYlBu_r
norm = plt.Normalize(vmin=min(medians), vmax=max(medians))
for patch, m in zip(bp["boxes"], medians):
    patch.set_facecolor(cmap(norm(m)))
    patch.set_alpha(0.8)

# Sample size annotation
for i, (d, c) in enumerate(zip(plot_data, clade_order), start=1):
    ax.text(i, ax.get_ylim()[0] - 0.15, f"n={len(d)}", ha="center", va="top",
            fontsize=7)

# KW annotation
ax.text(0.98, 0.97,
        f"Kruskal-Wallis\nH={kw_stat:.1f}, p={kw_p:.1e}, df={len(groups)-1}\n{n_sig} sig. pairs (Dunn BH)",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=8, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

ax.set_xticks(range(1, len(clade_order) + 1))
ax.set_xticklabels(clade_labels_ord, rotation=35, ha="right", fontsize=9)
ax.set_ylabel("% TE-coopted genes", fontsize=11)
ax.set_title("TE co-option burden by fungal/host clade", fontsize=13)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
plt.tight_layout()
fig.savefig(f"{OUT_DIR}/15_clade_boxplot.pdf", dpi=150)
plt.close(fig)
print(f"[OK] Boxplot → {OUT_DIR}/15_clade_boxplot.pdf")

# ---------------------------------------------------------------------------
# 5. TE superfamily heatmap — proportion of species per clade where
#    each superfamily is dominant
# ---------------------------------------------------------------------------
def build_freq_matrix(df_in, col, clade_col="clade", min_freq=2):
    """
    For each (clade, category) pair: count species where category == dominant.
    Returns a proportion matrix (clades × categories), keeping only categories
    that appear in ≥ min_freq species globally.
    """
    # Expand multi-value cells (some dominant fields contain '+'-joined labels)
    exploded = []
    for _, row in df_in.iterrows():
        val = str(row[col]) if pd.notna(row[col]) else "Unknown"
        val = val.strip()
        if val in ("\\N", "nan", "", "None"):
            val = "Unknown"
        exploded.append({"clade": row[clade_col], "category": val})
    exp_df = pd.DataFrame(exploded)

    # Raw counts
    ct = exp_df.groupby(["clade", "category"]).size().unstack(fill_value=0)

    # Filter low-frequency categories
    total_counts = ct.sum(axis=0)
    ct = ct[total_counts[total_counts >= min_freq].index]

    # Normalise per clade (proportion of species)
    clade_totals = df_in.groupby(clade_col).size()
    ct = ct.div(clade_totals, axis=0)
    return ct

def plot_heatmap(matrix, title, outpath, figsize=None, cbar_label="Proportion of species"):
    if matrix.empty or matrix.shape[1] == 0:
        print(f"[WARN] Skipping {outpath}: empty matrix")
        return
    # Sort clades by descending median co-option (use stats_df order)
    row_order = [c for c in clade_order if c in matrix.index]
    matrix = matrix.loc[row_order]
    row_labels = [CLADE_LABELS.get(c, c) for c in row_order]

    if figsize is None:
        figsize = (max(10, matrix.shape[1] * 0.55), max(5, len(row_order) * 0.4 + 2))

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix.values, aspect="auto", cmap="YlOrRd",
                   vmin=0, vmax=matrix.values.max())
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9)

    # Annotate cells with value if >0.05
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = matrix.values[i, j]
            if v > 0.05:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=6, color="black" if v < 0.6 else "white")

    plt.colorbar(im, ax=ax, label=cbar_label, shrink=0.6)
    ax.set_title(title, fontsize=12)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[OK] Heatmap → {outpath}")

# TE superfamily heatmap
te_sf_matrix = build_freq_matrix(df_valid, "dominant_te_superfam", min_freq=2)
plot_heatmap(te_sf_matrix,
             "Dominant TE superfamily per clade (proportion of species)",
             f"{OUT_DIR}/15_te_superfam_heatmap.pdf")

# Host domain heatmap
host_dom_matrix = build_freq_matrix(df_valid, "dominant_host_domain", min_freq=2)
plot_heatmap(host_dom_matrix,
             "Dominant host domain in co-opted genes per clade",
             f"{OUT_DIR}/15_host_domain_heatmap.pdf")

# ---------------------------------------------------------------------------
# 6. Scatter: pct_coopted vs pct_genome_TE (colour = clade)
# ---------------------------------------------------------------------------
df_scatter = df_valid[df_valid["pct_genome_TE"].notna() &
                      (df_valid["pct_genome_TE"] > 0)].copy()

if len(df_scatter) > 10:
    fig, ax = plt.subplots(figsize=(9, 7))
    clades_uniq = sorted(df_scatter["clade"].unique())
    cmap2 = plt.cm.get_cmap("tab20", len(clades_uniq))
    color_map = {c: cmap2(i) for i, c in enumerate(clades_uniq)}

    for clade, grp in df_scatter.groupby("clade"):
        ax.scatter(grp["pct_genome_TE"], grp["pct_coopted"],
                   label=CLADE_LABELS.get(clade, clade),
                   color=color_map[clade], s=40, alpha=0.75, edgecolors="none")

    # Spearman correlation
    rho, p_spear = ss.spearmanr(df_scatter["pct_genome_TE"], df_scatter["pct_coopted"])
    ax.text(0.02, 0.97,
            f"Spearman rho={rho:.3f}, p={p_spear:.2e}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

    ax.set_xlabel("% genome covered by TEs", fontsize=11)
    ax.set_ylabel("% TE-coopted genes", fontsize=11)
    ax.set_title("Co-option rate vs. TE genome content", fontsize=13)
    ax.legend(fontsize=7, ncol=2, loc="upper right",
              bbox_to_anchor=(1.0, 0.95), framealpha=0.7)
    plt.tight_layout()
    fig.savefig(f"{OUT_DIR}/15_scatter_te_vs_cooption.pdf", dpi=150)
    plt.close(fig)
    print(f"[OK] Scatter → {OUT_DIR}/15_scatter_te_vs_cooption.pdf")
else:
    print("[WARN] Not enough species with pct_genome_TE for scatter plot")

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
print("\n=== Summary ===")
print(f"Kruskal-Wallis: H={kw_stat:.3f}, p={kw_p:.2e} → ", end="")
if kw_p < 0.001:
    print("SIGNIFICANT (p<0.001) — clades differ in co-option rate")
elif kw_p < 0.05:
    print("SIGNIFICANT (p<0.05)")
else:
    print("NOT significant")
print(f"Pairwise significant pairs (Dunn BH<0.05): {n_sig}")
print("\nTop 5 clades by median co-option %:")
print(stats_df[["clade_label", "n_species", "median_pct", "sd_pct"]].head(5).to_string(index=False))
print("\nBottom 5 clades:")
print(stats_df[["clade_label", "n_species", "median_pct", "sd_pct"]].tail(5).to_string(index=False))
