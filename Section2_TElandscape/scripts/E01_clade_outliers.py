#!/usr/bin/env python3
"""
E01_clade_outliers.py
=====================
Detect species that are outliers within their phylogenetic clade in terms of:
  1. Total TE content (% genome)
  2. TE class composition (proportion vector)

Outlier criteria:
  - Univariate (total TE%): Z-score > Z_THRESH within clade
  - Multivariate (composition): Mahalanobis distance > chi2 p<0.05 threshold,
    fallback to Euclidean Z-score when n < MIN_N_MAHAL

Inputs (read directly from EarlGrey summaryFiles):
  EARLGREY_DIR/{SP}_EarlGrey/{SP}_summaryFiles/{SP}.highLevelCount.txt

Outputs:
  OUT_DIR/clade_outliers.tsv          — per-species outlier table
  OUT_DIR/clade_total_te_boxplot.pdf  — boxplot total TE% per clade, outliers labelled
  OUT_DIR/clade_composition_pca.pdf   — PCA of TE composition, coloured by clade,
                                        outliers marked with a star

Usage:
  python3 E01_clade_outliers.py [EARLGREY_DIR] [OUT_DIR]

Defaults:
  EARLGREY_DIR = /data/carlos/EDF-TE/TE_analysis/annotated
  OUT_DIR      = /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/clade_outliers
"""

import sys, os, glob, warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import mahalanobis
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ── config ────────────────────────────────────────────────────────────────────
EARLGREY_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/clade_outliers"

os.makedirs(OUT_DIR, exist_ok=True)

# Z-score threshold for univariate outlier
Z_THRESH = 2.5
# Minimum clade size for Mahalanobis (needs n > n_features)
MIN_N_MAHAL = 8

# TE classes to keep (merge nested into parent, drop "Other" and "Non-Repeat")
CLASS_MAP = {
    "DNA":                    "DNA",
    "DNA-nested":             "DNA",
    "Rolling Circle":         "Rolling_Circle",
    "Rolling Circle-nested":  "Rolling_Circle",
    "Penelope":               "Penelope",
    "Penelope-nested":        "Penelope",
    "LINE":                   "LINE",
    "LINE-nested":            "LINE",
    "SINE":                   "SINE",
    "SINE-nested":            "SINE",
    "LTR":                    "LTR",
    "LTR-nested":             "LTR",
    "Unclassified":           "Unclassified",
    "Unclassified-nested":    "Unclassified",
}
TE_CLASSES = sorted(set(CLASS_MAP.values()))

# ── 1. parse EarlGrey highLevelCount files ────────────────────────────────────
records = []
pattern = os.path.join(EARLGREY_DIR, "*_EarlGrey", "*_summaryFiles", "*.highLevelCount.txt")
files = glob.glob(pattern)
if not files:
    sys.exit(f"[ERROR] No highLevelCount.txt files found in {EARLGREY_DIR}")

print(f"[INFO] Found {len(files)} highLevelCount.txt files")

for fp in sorted(files):
    sp = os.path.basename(fp).replace(".highLevelCount.txt", "")
    clade = sp[:2]  # first two chars = phylogenetic group code
    row = {"species": sp, "clade": clade}

    bp_by_class = {c: 0.0 for c in TE_CLASSES}
    genome_size = None
    total_bp = 0.0

    try:
        df = pd.read_csv(fp, sep="\t")
        # Normalise column names (strip whitespace)
        df.columns = [c.strip() for c in df.columns]
        # Find relevant columns
        class_col = [c for c in df.columns if "Classification" in c][0]
        bp_col    = [c for c in df.columns if "Coverage (bp)" in c][0]
        pct_col   = [c for c in df.columns if "% Genome" in c][0]
        size_col  = [c for c in df.columns if "Genome Size" in c][0]

        for _, r in df.iterrows():
            cls = str(r[class_col]).strip()
            bp  = float(r[bp_col]) if str(r[bp_col]) not in ("NA", "nan", "") else 0.0
            if cls in CLASS_MAP:
                bp_by_class[CLASS_MAP[cls]] += bp
            if "Total Interspersed" in cls:
                total_bp = bp
            if genome_size is None:
                try:
                    genome_size = float(r[size_col])
                except (ValueError, TypeError):
                    pass

        if genome_size and genome_size > 0:
            row["total_TE_pct"]  = total_bp / genome_size * 100
            row["genome_size_mb"] = genome_size / 1e6
        else:
            row["total_TE_pct"]  = np.nan
            row["genome_size_mb"] = np.nan

        for c in TE_CLASSES:
            row[c] = bp_by_class[c]

    except Exception as e:
        print(f"[WARN] Skipping {sp}: {e}")
        continue

    records.append(row)

df_all = pd.DataFrame(records)
print(f"[INFO] Parsed {len(df_all)} species, {df_all['clade'].nunique()} clades")

# ── 2. compute TE class proportions (relative to total TE bp) ─────────────────
for c in TE_CLASSES:
    df_all[f"prop_{c}"] = df_all[c] / df_all[TE_CLASSES].sum(axis=1).replace(0, np.nan)
prop_cols = [f"prop_{c}" for c in TE_CLASSES]

# ── 3. per-clade outlier detection ───────────────────────────────────────────
outlier_rows = []

for clade, grp in df_all.groupby("clade"):
    n = len(grp)
    if n < 3:
        # too small to define outliers — flag all as "small clade"
        for sp in grp["species"]:
            outlier_rows.append({"species": sp, "clade": clade, "n_clade": n,
                                  "total_TE_pct": grp.loc[grp.species==sp,"total_TE_pct"].values[0],
                                  "z_total_TE": np.nan, "outlier_univariate": False,
                                  "mahal_dist": np.nan, "outlier_composition": False,
                                  "note": "small_clade"})
        continue

    # — univariate: Z-score on total TE% ——————————————————————
    vals = grp["total_TE_pct"].values.astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        z = stats.zscore(vals, nan_policy="omit")

    # — multivariate: Mahalanobis or Euclidean Z on TE composition ——————————
    comp = grp[prop_cols].fillna(0).values.astype(float)
    mahal_dists = np.full(n, np.nan)

    if n >= MIN_N_MAHAL:
        try:
            cov = np.cov(comp.T)
            # regularise if near-singular
            cov += np.eye(cov.shape[0]) * 1e-6
            VI = np.linalg.inv(cov)
            centroid = comp.mean(axis=0)
            mahal_dists = np.array([mahalanobis(row, centroid, VI) for row in comp])
            # chi2 threshold: df = n_features, p=0.05
            from scipy.stats import chi2
            chi2_thresh = chi2.ppf(0.95, df=comp.shape[1])
            outlier_comp = mahal_dists**2 > chi2_thresh
        except np.linalg.LinAlgError:
            # fallback: Euclidean distance Z-score
            centroid = comp.mean(axis=0)
            euc = np.linalg.norm(comp - centroid, axis=1)
            euc_z = stats.zscore(euc, nan_policy="omit")
            mahal_dists = euc
            outlier_comp = np.abs(euc_z) > Z_THRESH
    else:
        # Euclidean distance Z-score for small clades
        comp_z = np.abs(stats.zscore(comp, axis=0, nan_policy="omit"))
        # outlier if any single proportion is extreme
        outlier_comp = comp_z.max(axis=1) > Z_THRESH

    for i, (idx, row) in enumerate(grp.iterrows()):
        outlier_rows.append({
            "species":             row["species"],
            "clade":               clade,
            "n_clade":             n,
            "genome_size_mb":      row["genome_size_mb"],
            "total_TE_pct":        row["total_TE_pct"],
            "z_total_TE":          round(z[i], 3) if not np.isnan(z[i]) else np.nan,
            "outlier_univariate":  bool(abs(z[i]) > Z_THRESH) if not np.isnan(z[i]) else False,
            "mahal_dist":          round(float(mahal_dists[i]), 3) if not np.isnan(mahal_dists[i]) else np.nan,
            "outlier_composition": bool(outlier_comp[i]),
            "note":                ""
        })

df_out = pd.DataFrame(outlier_rows)

# add dominant TE class per species
dom = df_all[prop_cols].idxmax(axis=1).str.replace("prop_", "")
df_all["dominant_class"] = dom
df_out = df_out.merge(df_all[["species","dominant_class"] + prop_cols], on="species", how="left")

out_tsv = os.path.join(OUT_DIR, "clade_outliers.tsv")
df_out.to_csv(out_tsv, sep="\t", index=False)
print(f"[OK] Outlier table → {out_tsv}")

n_univ = df_out["outlier_univariate"].sum()
n_comp = df_out["outlier_composition"].sum()
n_both = (df_out["outlier_univariate"] & df_out["outlier_composition"]).sum()
print(f"[INFO] Outliers — univariate: {n_univ}  composition: {n_comp}  both: {n_both}")
print(df_out[df_out["outlier_univariate"] | df_out["outlier_composition"]]
      [["species","clade","total_TE_pct","z_total_TE","dominant_class","outlier_univariate","outlier_composition"]]
      .to_string(index=False))

# ── 4. PLOT 1: boxplot total TE% per clade ────────────────────────────────────
clade_order = sorted(df_all["clade"].unique())
clade_colors = plt.cm.tab20(np.linspace(0, 1, len(clade_order)))
color_map = dict(zip(clade_order, clade_colors))

fig, ax = plt.subplots(figsize=(max(10, len(clade_order) * 0.8), 6))
bp_data = [df_all[df_all["clade"] == cl]["total_TE_pct"].dropna().values for cl in clade_order]
bp = ax.boxplot(bp_data, patch_artist=True, notch=False, showfliers=False,
                medianprops=dict(color="black", lw=2))

for patch, cl in zip(bp["boxes"], clade_order):
    patch.set_facecolor(color_map[cl])
    patch.set_alpha(0.7)

# overlay individual points + label outliers
for i, cl in enumerate(clade_order, 1):
    sub = df_all[df_all["clade"] == cl]
    out_flags = df_out[df_out["clade"] == cl].set_index("species")["outlier_univariate"]
    ys = sub["total_TE_pct"].values
    xs = np.random.normal(i, 0.08, size=len(ys))
    ax.scatter(xs, ys, color=color_map[cl], s=25, alpha=0.6, zorder=3)
    # label outliers
    for _, r in sub.iterrows():
        if out_flags.get(r["species"], False):
            ax.annotate(r["species"], (i, r["total_TE_pct"]),
                        fontsize=5, ha="center", va="bottom",
                        xytext=(0, 4), textcoords="offset points",
                        color="red", fontweight="bold")
            ax.scatter([i], [r["total_TE_pct"]], color="red", s=60,
                       marker="*", zorder=5)

ax.set_xticks(range(1, len(clade_order) + 1))
ax.set_xticklabels(clade_order, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Total TE content (% genome)", fontsize=11)
ax.set_title(f"Total TE content by clade — outliers (|Z|>{Z_THRESH}) in red", fontsize=12)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plot1 = os.path.join(OUT_DIR, "clade_total_te_boxplot.pdf")
fig.savefig(plot1, dpi=150)
plt.close()
print(f"[OK] Boxplot → {plot1}")

# ── 5. PLOT 2: PCA of TE composition coloured by clade ───────────────────────
comp_mat = df_all[prop_cols].fillna(0).values
mask = ~np.isnan(comp_mat).any(axis=1)
comp_clean = comp_mat[mask]
sp_clean   = df_all["species"].values[mask]
cl_clean   = df_all["clade"].values[mask]

scaler = StandardScaler()
comp_scaled = scaler.fit_transform(comp_clean)

pca = PCA(n_components=2)
coords = pca.fit_transform(comp_scaled)
var_exp = pca.explained_variance_ratio_ * 100

fig, ax = plt.subplots(figsize=(10, 8))

for cl in sorted(set(cl_clean)):
    idx = cl_clean == cl
    ax.scatter(coords[idx, 0], coords[idx, 1],
               color=color_map.get(cl, "grey"), label=cl, s=30, alpha=0.7)

# Mark composition outliers with stars
outlier_sp_set = set(df_out[df_out["outlier_composition"]]["species"])
for j, sp in enumerate(sp_clean):
    if sp in outlier_sp_set:
        ax.scatter(coords[j, 0], coords[j, 1], color="red", s=120,
                   marker="*", zorder=5)
        ax.annotate(sp, (coords[j, 0], coords[j, 1]), fontsize=5,
                    xytext=(4, 2), textcoords="offset points", color="red")

# loadings
loadings = pca.components_.T
feat_labels = [c.replace("prop_", "") for c in prop_cols]
scale = 2.0
for k, lab in enumerate(feat_labels):
    ax.arrow(0, 0, loadings[k, 0]*scale, loadings[k, 1]*scale,
             head_width=0.05, head_length=0.05, fc="0.4", ec="0.4", lw=1)
    ax.text(loadings[k, 0]*scale*1.1, loadings[k, 1]*scale*1.1,
            lab, fontsize=7, color="0.3", ha="center")

ax.axhline(0, color="0.8", lw=0.5)
ax.axvline(0, color="0.8", lw=0.5)
ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)", fontsize=11)
ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)", fontsize=11)
ax.set_title("PCA of TE class composition — composition outliers (red ★)", fontsize=12)

legend_patches = [mpatches.Patch(color=color_map.get(cl, "grey"), label=cl)
                  for cl in sorted(set(cl_clean))]
ax.legend(handles=legend_patches, fontsize=7, ncol=2,
          loc="upper right", framealpha=0.7)
plt.tight_layout()
plot2 = os.path.join(OUT_DIR, "clade_composition_pca.pdf")
fig.savefig(plot2, dpi=150)
plt.close()
print(f"[OK] PCA → {plot2}")

# ── 6. per-clade summary table ────────────────────────────────────────────────
summary = df_out.groupby("clade").agg(
    n_species=("species", "count"),
    mean_TE_pct=("total_TE_pct", "mean"),
    sd_TE_pct=("total_TE_pct", "std"),
    min_TE_pct=("total_TE_pct", "min"),
    max_TE_pct=("total_TE_pct", "max"),
    n_outliers_univ=("outlier_univariate", "sum"),
    n_outliers_comp=("outlier_composition", "sum"),
).round(3)
summary_tsv = os.path.join(OUT_DIR, "clade_summary.tsv")
summary.to_csv(summary_tsv, sep="\t")
print(f"[OK] Clade summary → {summary_tsv}")
print(summary.to_string())
