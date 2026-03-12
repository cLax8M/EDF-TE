#!/usr/bin/env python3
"""
E02_clade_composition_viz.py
============================
Visualise TE class composition outliers within each clade.

Outputs:
  1. clade_stacked_bars.pdf  — one subplot per clade, each species = stacked bar
                               coloured by TE class; outlier species labelled in red
  2. clade_deviation_heatmap.pdf — heatmap of (species_proportion - clade_mean)
                               per TE class; clades separated by white lines;
                               outlier species marked with ★

Reads:
  - EarlGrey highLevelCount.txt files  (same as E01)
  - clade_outliers.tsv produced by E01 (to get outlier flags)

Usage:
  python3 E02_clade_composition_viz.py [EARLGREY_DIR] [OUTLIERS_TSV] [OUT_DIR]
"""

import sys, os, glob, warnings
import numpy as np
import pandas as pd
from scipy.stats import kruskal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns

# ── config ────────────────────────────────────────────────────────────────────
EARLGREY_DIR  = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/carlos/EDF-TE/TE_analysis/annotated"
OUTLIERS_TSV  = sys.argv[2] if len(sys.argv) > 2 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/clade_outliers/clade_outliers.tsv"
OUT_DIR       = sys.argv[3] if len(sys.argv) > 3 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/clade_outliers"

os.makedirs(OUT_DIR, exist_ok=True)

CLASS_MAP = {
    "DNA":                   "DNA",
    "DNA-nested":            "DNA",
    "Rolling Circle":        "Rolling Circle",
    "Rolling Circle-nested": "Rolling Circle",
    "Penelope":              "Penelope",
    "Penelope-nested":       "Penelope",
    "LINE":                  "LINE",
    "LINE-nested":           "LINE",
    "SINE":                  "SINE",
    "SINE-nested":           "SINE",
    "LTR":                   "LTR",
    "LTR-nested":            "LTR",
    "Unclassified":          "Unclassified",
    "Unclassified-nested":   "Unclassified",
}
TE_CLASSES = ["LTR", "DNA", "LINE", "SINE", "Rolling Circle", "Penelope", "Unclassified"]

# Colour palette — one per TE class, fixed
CLASS_COLORS = {
    "LTR":           "#e63946",
    "DNA":           "#457b9d",
    "LINE":          "#2a9d8f",
    "SINE":          "#e9c46a",
    "Rolling Circle":"#f4a261",
    "Penelope":      "#a8dadc",
    "Unclassified":  "#adb5bd",
}

# ── 1. parse EarlGrey files ───────────────────────────────────────────────────
records = []
for fp in sorted(glob.glob(os.path.join(EARLGREY_DIR,
                 "*_EarlGrey", "*_summaryFiles", "*.highLevelCount.txt"))):
    sp = os.path.basename(fp).replace(".highLevelCount.txt", "")
    clade = sp[:2]
    bp_by_class = {c: 0.0 for c in TE_CLASSES}
    genome_size = None
    total_bp = 0.0
    try:
        df = pd.read_csv(fp, sep="\t")
        df.columns = [c.strip() for c in df.columns]
        class_col = [c for c in df.columns if "Classification" in c][0]
        bp_col    = [c for c in df.columns if "Coverage (bp)"   in c][0]
        size_col  = [c for c in df.columns if "Genome Size"      in c][0]
        for _, r in df.iterrows():
            cls = str(r[class_col]).strip()
            bp  = float(r[bp_col]) if str(r[bp_col]) not in ("NA","nan","") else 0.0
            if cls in CLASS_MAP:
                bp_by_class[CLASS_MAP[cls]] += bp
            if "Total Interspersed" in cls:
                total_bp = bp
            if genome_size is None:
                try: genome_size = float(r[size_col])
                except: pass
        row = {"species": sp, "clade": clade,
               "total_TE_pct": total_bp / genome_size * 100 if genome_size else np.nan}
        for c in TE_CLASSES:
            row[c] = bp_by_class[c]
        records.append(row)
    except Exception as e:
        print(f"[WARN] {sp}: {e}")

df_all = pd.DataFrame(records)

# Compute proportions (fraction of total TE bp)
total_te_bp = df_all[TE_CLASSES].sum(axis=1).replace(0, np.nan)
for c in TE_CLASSES:
    df_all[f"prop_{c}"] = df_all[c] / total_te_bp
prop_cols = [f"prop_{c}" for c in TE_CLASSES]

print(f"[INFO] {len(df_all)} species, {df_all['clade'].nunique()} clades")

# ── 2. load outlier flags ─────────────────────────────────────────────────────
if os.path.exists(OUTLIERS_TSV):
    df_flags = pd.read_csv(OUTLIERS_TSV, sep="\t")[
        ["species", "outlier_univariate", "outlier_composition"]]
    df_all = df_all.merge(df_flags, on="species", how="left")
    df_all["outlier_univariate"]  = df_all["outlier_univariate"].fillna(False)
    df_all["outlier_composition"] = df_all["outlier_composition"].fillna(False)
else:
    df_all["outlier_univariate"]  = False
    df_all["outlier_composition"] = False

# ── 3. PLOT 1: stacked bars per clade ─────────────────────────────────────────
clades = sorted(df_all["clade"].unique())
n_clades = len(clades)

# Figure layout: grid of subplots, ~4 per row
ncols = 4
nrows = int(np.ceil(n_clades / ncols))
fig, axes = plt.subplots(nrows, ncols,
                          figsize=(ncols * 5, nrows * 4.5),
                          squeeze=False)
fig.suptitle("TE class composition by clade\n"
             "(red labels = composition outliers  |  ★ = total-TE outlier)",
             fontsize=14, y=1.01)

for ax_idx, clade in enumerate(clades):
    row_i, col_i = divmod(ax_idx, ncols)
    ax = axes[row_i][col_i]

    sub = df_all[df_all["clade"] == clade].copy()
    # Sort species by total TE% so bars go low→high
    sub = sub.sort_values("total_TE_pct", ascending=True).reset_index(drop=True)
    n = len(sub)

    x = np.arange(n)
    bottom = np.zeros(n)
    for cls in TE_CLASSES:
        vals = sub[f"prop_{cls}"].fillna(0).values
        ax.bar(x, vals, bottom=bottom, color=CLASS_COLORS[cls],
               label=cls, width=0.85, edgecolor="none")
        bottom += vals

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, 1)
    ax.set_title(f"{clade}  (n={n})", fontsize=10, fontweight="bold")
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=7)
    ax.set_xticks(x)

    # X-tick labels: plain text (coloured below via annotations)
    tick_labels = []
    for _, sp_row in sub.iterrows():
        sp = sp_row["species"]
        short = sp[3:] if sp[:2] == clade else sp
        tick_labels.append(short)

    ax.set_xticklabels(tick_labels, rotation=90, fontsize=5.5)

    # Colour labels and add markers AFTER drawing ticks
    ax.figure.canvas.draw()  # needed to position tick labels
    for tick, (_, sp_row) in zip(ax.get_xticklabels(), sub.iterrows()):
        is_comp = sp_row["outlier_composition"]
        is_univ = sp_row["outlier_univariate"]
        if is_comp or is_univ:
            tick.set_color("red")
            tick.set_fontweight("bold")

    # Add markers above bars: ● (coloured dot = driving TE class) and/or *
    # Compute which TE class deviates most for composition outliers
    if os.path.exists(OUTLIERS_TSV):
        df_dev_clade = df_all[df_all["clade"] == clade].copy()
        clade_mean = df_dev_clade[prop_cols].mean()
        for i, (_, sp_row) in enumerate(sub.iterrows()):
            is_comp = sp_row["outlier_composition"]
            is_univ = sp_row["outlier_univariate"]
            y_top = 1.02
            offset = 0
            if is_univ:
                ax.text(i, y_top + offset, "*", ha="center", va="bottom",
                        fontsize=8, color="black", fontweight="bold")
                offset += 0.06
            if is_comp:
                # Find which class drives the outlier (max |deviation|)
                devs = {c: abs(sp_row[f"prop_{c}"] - clade_mean[f"prop_{c}"])
                        for c in TE_CLASSES}
                driver = max(devs, key=devs.get)
                ax.plot(i, y_top + offset, "o",
                        color=CLASS_COLORS[driver], markersize=5,
                        markeredgecolor="black", markeredgewidth=0.4,
                        transform=ax.get_xaxis_transform(),
                        clip_on=False)

    ax.grid(axis="y", alpha=0.25, lw=0.5)
    ax.spines[["top","right"]].set_visible(False)

# Hide empty subplots
for ax_idx in range(n_clades, nrows * ncols):
    row_i, col_i = divmod(ax_idx, ncols)
    axes[row_i][col_i].set_visible(False)

# Shared legend
legend_patches = [mpatches.Patch(color=CLASS_COLORS[c], label=c) for c in TE_CLASSES]
fig.legend(handles=legend_patches, loc="lower center", ncol=len(TE_CLASSES),
           fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))

plt.tight_layout(rect=[0, 0.04, 1, 1])
out1 = os.path.join(OUT_DIR, "clade_stacked_bars.pdf")
fig.savefig(out1, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Stacked bars → {out1}")

# ── 3b. PLOT 1b: EDF-only stacked bars (≥3 species) ─────────────────────────
EDF_CLADES = {"AP","BL","BS","CH","CR","EN","GL","KX","MI","MO","MU","OL","ZO"}
edf_clades_ok = [cl for cl in sorted(df_all["clade"].unique())
                 if cl in EDF_CLADES and (df_all["clade"] == cl).sum() >= 3]
print(f"[INFO] EDF clades with ≥3 species: {edf_clades_ok}")

n_edf = len(edf_clades_ok)
ncols_e = 4
nrows_e = int(np.ceil(n_edf / ncols_e))
fig_e, axes_e = plt.subplots(nrows_e, ncols_e,
                               figsize=(ncols_e * 5, nrows_e * 4.5),
                               squeeze=False)
fig_e.suptitle("TE class composition — EDF clades (≥3 genomes)\n"
               "(red labels = composition outliers  |  ★ = total-TE outlier)",
               fontsize=14, y=1.01)

for ax_idx, clade in enumerate(edf_clades_ok):
    row_i, col_i = divmod(ax_idx, ncols_e)
    ax = axes_e[row_i][col_i]

    sub = df_all[df_all["clade"] == clade].copy()
    sub = sub.sort_values("total_TE_pct", ascending=True).reset_index(drop=True)
    n = len(sub)
    x = np.arange(n)
    bottom = np.zeros(n)
    for cls in TE_CLASSES:
        vals = sub[f"prop_{cls}"].fillna(0).values
        ax.bar(x, vals, bottom=bottom, color=CLASS_COLORS[cls],
               label=cls, width=0.85, edgecolor="none")
        bottom += vals

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, 1)
    ax.set_title(f"{clade}  (n={n})", fontsize=10, fontweight="bold")
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=7)
    ax.set_xticks(x)

    tick_labels = []
    for _, sp_row in sub.iterrows():
        sp = sp_row["species"]
        short = sp[3:] if sp[:2] == clade else sp
        tick_labels.append(short)
    ax.set_xticklabels(tick_labels, rotation=90, fontsize=5.5)

    ax.figure.canvas.draw()
    for tick, (_, sp_row) in zip(ax.get_xticklabels(), sub.iterrows()):
        if sp_row["outlier_composition"] or sp_row["outlier_univariate"]:
            tick.set_color("red")
            tick.set_fontweight("bold")

    df_dev_clade = df_all[df_all["clade"] == clade].copy()
    clade_mean = df_dev_clade[prop_cols].mean()
    for i, (_, sp_row) in enumerate(sub.iterrows()):
        is_comp = sp_row["outlier_composition"]
        is_univ = sp_row["outlier_univariate"]
        y_top = 1.02
        offset = 0
        if is_univ:
            ax.text(i, y_top + offset, "*", ha="center", va="bottom",
                    fontsize=8, color="black", fontweight="bold")
            offset += 0.06
        if is_comp:
            devs = {c: abs(sp_row[f"prop_{c}"] - clade_mean[f"prop_{c}"])
                    for c in TE_CLASSES}
            driver = max(devs, key=devs.get)
            ax.plot(i, y_top + offset, "o",
                    color=CLASS_COLORS[driver], markersize=5,
                    markeredgecolor="black", markeredgewidth=0.4,
                    transform=ax.get_xaxis_transform(),
                    clip_on=False)
    ax.grid(axis="y", alpha=0.25, lw=0.5)
    ax.spines[["top","right"]].set_visible(False)

for ax_idx in range(n_edf, nrows_e * ncols_e):
    row_i, col_i = divmod(ax_idx, ncols_e)
    axes_e[row_i][col_i].set_visible(False)

fig_e.legend(handles=legend_patches, loc="lower center", ncol=len(TE_CLASSES),
             fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.04, 1, 1])
out1b = os.path.join(OUT_DIR, "clade_stacked_bars_EDF_only.pdf")
fig_e.savefig(out1b, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] EDF-only stacked bars → {out1b}")

# ── 4. PLOT 2: deviation-from-clade-mean heatmap ──────────────────────────────
# For each species compute: prop_class - clade_mean_prop_class
records_dev = []
for clade, grp in df_all.groupby("clade"):
    clade_mean = grp[prop_cols].mean()
    for _, r in grp.iterrows():
        dev = {c.replace("prop_",""): r[c] - clade_mean[c] for c in prop_cols}
        dev["species"] = r["species"]
        dev["clade"]   = clade
        dev["outlier_composition"] = r["outlier_composition"]
        dev["outlier_univariate"]  = r["outlier_univariate"]
        dev["total_TE_pct"]        = r["total_TE_pct"]
        records_dev.append(dev)

df_dev = pd.DataFrame(records_dev)

# Sort: by clade, then by total_TE_pct within clade
df_dev = df_dev.sort_values(["clade", "total_TE_pct"])

# ── significance: permutation test per species ────────────────────────────────
# For each species compute the total absolute deviation from clade mean.
# Then permute (reshuffle rows within clade) 999× to get null distribution.
# p-value = fraction of permutations where max_abs_dev ≥ observed.
np.random.seed(42)
N_PERM = 999

def _compute_abs_dev(sub_props, clade_mean):
    return np.abs(sub_props - clade_mean).sum()

sig_neg_log10_p = []
for _, row in df_dev.iterrows():
    clade = row["clade"]
    grp = df_all[df_all["clade"] == clade][prop_cols].values   # shape (n, 7)
    n_grp = len(grp)
    sp_vals = np.array([row[c.replace("prop_","")] + 
                        df_all[df_all["clade"]==clade][c].mean() - 
                        df_all[df_all["clade"]==clade][c].mean()
                        for c in prop_cols])
    # easier: just recompute from df_all
    sp_vals = df_all[df_all["species"]==row["species"]][prop_cols].values[0]
    clade_mean = grp.mean(axis=0)
    obs = np.abs(sp_vals - clade_mean).sum()
    if n_grp < 3:
        sig_neg_log10_p.append(np.nan)
        continue
    # null: for each permutation, pick a random row from same clade and compute its deviation
    null_obs = np.array([np.abs(grp[np.random.randint(n_grp)] - clade_mean).sum()
                         for _ in range(N_PERM)])
    p = (null_obs >= obs).sum() / N_PERM
    p = max(p, 1 / N_PERM)   # floor at min p
    sig_neg_log10_p.append(-np.log10(p))

df_dev["neg_log10_p"] = sig_neg_log10_p

# Build matrix
heat_mat = df_dev[TE_CLASSES].values          # shape (n_sp, n_classes)
species_labels = df_dev["species"].values
clade_labels   = df_dev["clade"].values

# Row colours for clade annotation
clade_palette = sns.color_palette("tab20", n_colors=len(clades))
clade_color_map = dict(zip(clades, clade_palette))
row_colors = pd.Series(clade_labels, name="Clade").map(clade_color_map)

# Build short labels for y-axis
short_labels = []
for sp, cl in zip(species_labels, clade_labels):
    short = sp[3:] if sp.startswith(cl) else sp
    short_labels.append(short)

fig_h = max(14, len(species_labels) * 0.22)
fig, axes_row = plt.subplots(
    1, 3, figsize=(12, fig_h),
    gridspec_kw={"width_ratios": [1, 0.06, 0.04]})
ax_heat, ax_sig, ax_cbar = axes_row

# Draw heatmap
vmax = max(abs(heat_mat.max()), abs(heat_mat.min()))
vmax = max(vmax, 0.1)
im = ax_heat.imshow(heat_mat, aspect="auto", cmap="RdBu_r",
                    vmin=-vmax, vmax=vmax, interpolation="nearest")

# Y-axis
ax_heat.set_yticks(range(len(species_labels)))
ax_heat.set_yticklabels(short_labels, fontsize=5)

# Colour outlier y-labels red
for i, (sp, comp_out, univ_out) in enumerate(
        zip(species_labels, df_dev["outlier_composition"], df_dev["outlier_univariate"])):
    if comp_out or univ_out:
        ax_heat.get_yticklabels()[i].set_color("red")
        ax_heat.get_yticklabels()[i].set_fontweight("bold")

# Clade annotation strip on the left
ax_heat.set_yticks(range(len(species_labels)))  # ensure ticks set
clade_strip_x = -0.8  # in data coords; we'll use a twin axis instead

# Draw coloured rectangles along y-axis for clade annotation
# Use inset axes as colour strip
ax_strip = ax_heat.inset_axes([-0.12, 0, 0.10, 1], transform=ax_heat.transAxes)
ax_strip.set_xlim(0, 1)
ax_strip.set_ylim(-0.5, len(species_labels) - 0.5)
ax_strip.invert_yaxis()
ax_strip.axis("off")

prev_clade = None
block_start = 0
for i, cl in enumerate(clade_labels):
    if cl != prev_clade:
        if prev_clade is not None:
            ax_strip.barh(y=(block_start + i - 1) / 2,
                          width=1, height=i - block_start,
                          left=0, color=clade_color_map[prev_clade], alpha=0.9)
            ax_strip.text(0.5, (block_start + i - 1) / 2, prev_clade,
                          ha="center", va="center", fontsize=6,
                          color="white", fontweight="bold")
        block_start = i
        prev_clade = cl
# last block
ax_strip.barh(y=(block_start + len(clade_labels) - 1) / 2,
              width=1, height=len(clade_labels) - block_start,
              left=0, color=clade_color_map[prev_clade], alpha=0.9)
ax_strip.text(0.5, (block_start + len(clade_labels) - 1) / 2, prev_clade,
              ha="center", va="center", fontsize=6,
              color="white", fontweight="bold")

# Draw horizontal lines between clades
clade_breaks = np.where(np.array(clade_labels[:-1]) != np.array(clade_labels[1:]))[0]
for b in clade_breaks:
    ax_heat.axhline(b + 0.5, color="white", lw=1.5)

# Mark outlier species with ★ on the right
for i, (comp_out, univ_out) in enumerate(
        zip(df_dev["outlier_composition"], df_dev["outlier_univariate"])):
    if comp_out or univ_out:
        ax_heat.text(len(TE_CLASSES) - 0.4, i, "★",
                     ha="left", va="center", fontsize=7, color="red")

# X-axis
ax_heat.set_xticks(range(len(TE_CLASSES)))
ax_heat.set_xticklabels(TE_CLASSES, rotation=40, ha="right", fontsize=8)
ax_heat.set_title("TE class composition deviation from clade mean\n"
                  "(red = outlier species;  ★ = significant;  blue = less, red = more)",
                  fontsize=10, pad=10)

# Colorbar
plt.colorbar(im, cax=ax_cbar, label="Δ proportion (species − clade mean)")

# Significance strip: -log10(p) from permutation test
sig_vals = df_dev["neg_log10_p"].values
sig_mat = sig_vals[:, np.newaxis]
sig_max = np.nanmax(sig_vals) if not np.all(np.isnan(sig_vals)) else 1.0
ax_sig.imshow(sig_mat, aspect="auto", cmap="Purples",
              vmin=0, vmax=sig_max, interpolation="nearest")
ax_sig.set_xticks([0])
ax_sig.set_xticklabels(["-log₁₀(p)"], rotation=40, ha="right", fontsize=7)
ax_sig.set_yticks([])
# Add asterisk annotations for significance thresholds
for i, v in enumerate(sig_vals):
    if np.isnan(v):
        continue
    if v >= -np.log10(0.01):   # p < 0.01
        ax_sig.text(0, i, "**", ha="center", va="center", fontsize=5,
                    color="white", fontweight="bold")
    elif v >= -np.log10(0.05): # p < 0.05
        ax_sig.text(0, i, "*", ha="center", va="center", fontsize=5,
                    color="white", fontweight="bold")
ax_sig.spines[:].set_visible(False)

plt.tight_layout()
out2 = os.path.join(OUT_DIR, "clade_deviation_heatmap.pdf")
fig.savefig(out2, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Heatmap → {out2}")
