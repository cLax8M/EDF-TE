#!/usr/bin/env python3
"""
fig_htt_te_enrichment.py
Enrichment analysis: are certain TE classes over/under-represented among
HTT candidates compared to their genomic abundance (all_tes_info.tsv)?

Obs  = count of HTT candidates per TE class (via qseqid)
Exp  = (class_loci / total_loci) × total_HTT_candidates
Test = chi-square goodness-of-fit
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import chisquare
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

BASE    = Path("/legserv/Temp/Carlos/EDF-TE/results/HTT")
HTT_TSV = BASE / "06_htt_candidates/htt_candidates/htt_candidates.tsv"
TE_INFO = BASE / "05_te_classification/all_tes_info.tsv"
OUT_DIR = BASE / "figures"
OUT_PDF = OUT_DIR / "figHTT_te_enrichment.pdf"
OUT_PNG = OUT_DIR / "figHTT_te_enrichment.png"
OUT_TSV = OUT_DIR / "figHTT_te_enrichment.tsv"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 1. Background: TE loci abundance per class ────────────────────────────────
print("Computing background TE class abundance (chunked read)…", flush=True)
bg_counts: dict[str, int] = {}
for chunk in pd.read_csv(TE_INFO, sep="\t", usecols=["te_class"],
                          chunksize=1_000_000, low_memory=False):
    for cls, n in chunk["te_class"].value_counts().items():
        bg_counts[cls] = bg_counts.get(cls, 0) + n
bg = pd.Series(bg_counts, name="bg_loci").rename_axis("te_class")
print(f"  Total TE loci: {bg.sum():,}")
print(f"  Classes: {bg.to_dict()}")

# ── 2. Observed: HTT candidates per class ────────────────────────────────────
print("Loading HTT candidates…", flush=True)
htt = pd.read_csv(HTT_TSV, sep="\t")

# Map qseqid → te_class using a fast chunked merge
print("Mapping qseqid → te_class…", flush=True)
te_map: dict[str, str] = {}
for chunk in pd.read_csv(TE_INFO, sep="\t", usecols=["te_id", "te_class"],
                          chunksize=1_000_000, low_memory=False):
    hits = chunk[chunk["te_id"].isin(htt["qseqid"])]
    for _, row in hits.iterrows():
        te_map[row["te_id"]] = row["te_class"]

htt["te_class"] = htt["qseqid"].map(te_map).fillna("Unknown")
obs = htt["te_class"].value_counts().rename("obs")
print(f"  HTT per class:\n{obs.to_string()}")

# ── 3. Build enrichment table ─────────────────────────────────────────────────
EXCLUDE = {"Unknown"}   # too noisy; exclude from test and plot
df = (pd.DataFrame({"obs": obs, "bg": bg})
        .fillna(0)
        .reset_index()
        .rename(columns={"index": "te_class"}))
df = df[~df["te_class"].isin(EXCLUDE)].copy()
df = df[df["bg"] > 0].copy()

total_obs = df["obs"].sum()
total_bg  = df["bg"].sum()
df["expected"]    = (df["bg"] / total_bg) * total_obs
df["obs_exp"]     = df["obs"] / df["expected"].replace(0, np.nan)
df["log2_obs_exp"]= np.log2(df["obs_exp"])
df["bg_frac"]     = df["bg"] / total_bg
df["obs_frac"]    = df["obs"] / total_obs
df = df.sort_values("log2_obs_exp", ascending=False)

# For the plot only include classes with obs > 0 (avoids -inf in bars)
df_plot = df[df["obs"] > 0].copy().sort_values("log2_obs_exp", ascending=True)

# ── 4. Chi-square goodness-of-fit ─────────────────────────────────────────────
stat, pval = chisquare(df["obs"], df["expected"])
print(f"\nChi-square: χ²={stat:.2f}  df={len(df)-1}  p={pval:.2e}")

# Per-class residuals
df["residual"] = (df["obs"] - df["expected"]) / np.sqrt(df["expected"])

# ── 5. Save table ─────────────────────────────────────────────────────────────
df.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4f")
print(f"\nTable saved: {OUT_TSV}")
print(df[["te_class","obs","expected","obs_exp","log2_obs_exp","bg_frac","residual"]].to_string(index=False))

# ── 6. Plot ───────────────────────────────────────────────────────────────────
TE_COLORS = {
    "DNA":         "#e07030",
    "LINE":        "#4080c0",
    "LTR":         "#50a050",
    "SINE":        "#c05090",
    "RC":          "#a0a030",
    "PLE":         "#806020",
    "Retroposon":  "#30b0b0",
    "SINE?":       "#c090c0",
    "LINE?":       "#90b0e0",
}
DEFAULT_COL = "#888888"

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# ── Panel A: log2(obs/exp) bar chart ─────────────────────────────────────────
ax = axes[0]
colors = [TE_COLORS.get(c, DEFAULT_COL) for c in df_plot["te_class"]]
bars = ax.barh(df_plot["te_class"], df_plot["log2_obs_exp"], color=colors,
               edgecolor="white", linewidth=0.4)
ax.axvline(0, color="black", lw=0.8, ls="--")
ax.set_xlabel("log₂(Observed / Expected)", fontsize=10)
ax.set_title(f"HTT enrichment by TE class\nχ²={stat:.1f}, p={pval:.2e}",
             fontsize=10, fontweight="bold")
ax.set_ylabel("TE class")

# Add n_obs labels outside each bar
for bar, (_, row) in zip(bars, df_plot.iterrows()):
    w = bar.get_width()
    pad = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.015
    x_off = w + pad if w >= 0 else w - pad
    ha    = "left" if w >= 0 else "right"
    ax.text(x_off, bar.get_y() + bar.get_height() / 2,
            f"n={int(row['obs'])}", va="center", ha=ha, fontsize=7.5, color="black")

# ── Panel B: obs vs exp bubble/dot plot ──────────────────────────────────────
ax2 = axes[1]
xs = df["bg_frac"] * 100   # genome % (background)
ys = df["obs_frac"] * 100  # HTT %

for i, (_, row) in enumerate(df_plot.iterrows()):
    col = TE_COLORS.get(row["te_class"], DEFAULT_COL)
    ax2.scatter(row["bg_frac"] * 100, row["obs_frac"] * 100,
                s=max(20, row["obs"] * 4), color=col, edgecolors="white",
                linewidth=0.5, zorder=3, label=row["te_class"])
    ax2.annotate(row["te_class"],
                 (row["bg_frac"] * 100, row["obs_frac"] * 100),
                 textcoords="offset points", xytext=(5, 3),
                 fontsize=7.5, color=col)

# Diagonal = expected
max_val = max(xs.max(), ys.max()) * 1.1
ax2.plot([0, max_val], [0, max_val], ls="--", color="grey", lw=0.9, zorder=1)
ax2.set_xlabel("Genomic abundance (% of all TE loci)", fontsize=10)
ax2.set_ylabel("Fraction of HTT candidates (%)", fontsize=10)
ax2.set_title("Observed vs expected\n(dashed = expected if proportional to abundance)",
              fontsize=10, fontweight="bold")
ax2.set_xlim(left=0)
ax2.set_ylim(bottom=0)

plt.tight_layout(pad=2.0)
fig.savefig(str(OUT_PDF), dpi=200, bbox_inches="tight")
fig.savefig(str(OUT_PNG), dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\nFigure saved: {OUT_PDF}")
