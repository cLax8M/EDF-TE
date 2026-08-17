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
    # Apply NE_SPECIES override (e.g. CH_Anasp1 → NE)
    import sys as _sys_local, os as _os_local
    _sys_local.path.insert(0, _os_local.path.dirname(_os_local.path.abspath(__file__)))
    from clade_palette import get_clade as _gc_local
    clade = _gc_local(sp)
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
# Drop control yeasts (VL clade) globally — user request 2026-04-21
n_before = len(df_all)
df_all = df_all[~df_all["species"].str.startswith("VL_")].reset_index(drop=True)
if len(df_all) != n_before:
    print(f"[filter] dropped {n_before - len(df_all)} VL species")
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
valid_mask = df_all[prop_cols].notna().any(axis=1)
dom = pd.Series("Unknown", index=df_all.index)
dom[valid_mask] = df_all.loc[valid_mask, prop_cols].idxmax(axis=1).str.replace("prop_", "")
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

# ── 4. PLOT 1: total TE % per clade — half-violin + strip (raincloud) ──────
# Replaces the previous boxplot to avoid visual repetition with the other
# per-clade boxplots in the manuscript (Shannon, Kimura recency, etc.).
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clade_palette import (CLADE_COLORS, CLADE_NAMES,
                            CLADE_ORDER as _CL_PHYLO,
                            get_clade as _get_clade,
                            species_display as _spd)

# Phylogenetic ordering – keep only clades that are actually present
_present = set(df_all["clade"].unique())
clade_order = [c for c in _CL_PHYLO if c in _present] + sorted(_present - set(_CL_PHYLO))
color_map = {c: CLADE_COLORS.get(c, "#cccccc") for c in clade_order}

fig, ax = plt.subplots(figsize=(max(11, len(clade_order) * 0.85), 6.5))

# 1) Half-violin (right side) + strip (left side) for every clade
for i, cl in enumerate(clade_order, 1):
    vals = df_all[df_all["clade"] == cl]["total_TE_pct"].dropna().values
    if len(vals) == 0:
        continue
    # Strip (jittered points) on the left of the tick
    xs = np.random.uniform(i - 0.30, i - 0.05, size=len(vals))
    ax.scatter(xs, vals, color=color_map[cl], s=28, alpha=0.75,
               edgecolor="white", linewidth=0.4, zorder=3)
    # Median crossbar
    med = float(np.median(vals))
    ax.plot([i - 0.32, i - 0.03], [med, med],
            color="black", lw=1.6, zorder=4)
    # Half-violin on the right of the tick
    if len(vals) >= 3:
        try:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(vals, bw_method=0.4)
            ymin = max(0, vals.min() - 2)
            ymax = min(100, vals.max() + 2)
            yy = np.linspace(ymin, ymax, 200)
            dd = kde(yy)
            dd = dd / dd.max() * 0.32  # max half-width 0.32
            ax.fill_betweenx(yy, i, i + dd, facecolor=color_map[cl],
                              alpha=0.45, edgecolor=color_map[cl], linewidth=0.6)
        except Exception:
            pass

# 2) Highlight + label total-TE % outliers on top
for i, cl in enumerate(clade_order, 1):
    out_rows = df_out[(df_out["clade"] == cl) & (df_out["outlier_univariate"])]
    for _, r in out_rows.iterrows():
        y = r["total_TE_pct"]
        ax.scatter([i - 0.15], [y], color="white",
                   edgecolor="black", linewidth=1.6, s=80, zorder=6)
        ax.annotate(_spd(r["species"], italic=True),
                    (i - 0.15, y),
                    fontsize=6.5, fontweight="bold",
                    xytext=(-6, 4), textcoords="offset points",
                    ha="right", color="black",
                    bbox=dict(facecolor="white", edgecolor="none",
                              alpha=0.7, pad=0.5),
                    zorder=7)

ax.set_xticks(range(1, len(clade_order) + 1))
ax.set_xticklabels([CLADE_NAMES.get(c, c) for c in clade_order],
                    rotation=35, ha="right", fontsize=9)
ax.set_ylabel("Total TE content (% genome)", fontsize=11)
ax.set_title(
    f"Total TE content by clade — strip + half-violin\n"
    f"black bar = median  |  outliers within clade (|z| > {Z_THRESH}) ringed and labelled",
    fontsize=11)
ax.grid(axis="y", alpha=0.3)
ax.set_xlim(0.4, len(clade_order) + 0.7)
ax.set_ylim(bottom=0)
plt.tight_layout()
plot1 = os.path.join(OUT_DIR, "clade_total_te_boxplot.pdf")
fig.savefig(plot1, dpi=150)
plt.close()
print(f"[OK] Strip+violin → {plot1}")

# ── 5. PLOT 2: PCA of TE composition coloured by clade ───────────────────────
# Drop the Homo sapiens outlier; its extreme PC coordinates compress the rest
# of the cloud and make the figure unreadable. The species is still kept in
# the per-clade outlier table.
PCA_DROP_SPECIES = {"HO_Homosap"}
_pca_keep_mask = ~df_all["species"].isin(PCA_DROP_SPECIES)
if (~_pca_keep_mask).any():
    print(f"[PCA] dropping for visualisation only: "
          f"{sorted(df_all.loc[~_pca_keep_mask, 'species'].tolist())}")
df_pca = df_all[_pca_keep_mask].reset_index(drop=True)
comp_mat = df_pca[prop_cols].fillna(0).values
mask = ~np.isnan(comp_mat).any(axis=1)
comp_clean = comp_mat[mask]
sp_clean   = df_pca["species"].values[mask]
cl_clean   = df_pca["clade"].values[mask]

scaler = StandardScaler()
comp_scaled = scaler.fit_transform(comp_clean)

pca = PCA(n_components=2)
coords = pca.fit_transform(comp_scaled)
var_exp = pca.explained_variance_ratio_ * 100

# Build per-species lookup of outlier flags for fast access
flag_univ_map = df_out.set_index("species")["outlier_univariate"].to_dict()
flag_comp_map = df_out.set_index("species")["outlier_composition"].to_dict()

# ── Emit two PCA versions ────────────────────────────────────────────────
#   Default (clade_composition_pca.pdf):  ring + label TE-content outliers
#                                          (univariate |z|>2)
#   Wide variant (clade_composition_pca_wide.pdf):  slightly less wide than
#                                          before, and ring + label
#                                          *composition* outliers
#                                          (Mahalanobis from clade centroid)
PCA_VARIANTS = [
    # stem,                          figsize,       s_norm, s_bold, ring_kind
    ("clade_composition_pca",       (12.5, 7.5),   56,     138,    "univariate"),
    ("clade_composition_pca_wide_genomesize",  (13.5, 7.5),   60,     150,    "composition"),
]

import sys as _sys_local
_sys_local.path.insert(0, "/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts")
from clade_palette import species_display as _spd
univ_sp_set = set(df_out[df_out["outlier_univariate"]]["species"])
comp_sp_set = set(df_out[df_out["outlier_composition"]]["species"])

# Force Arial for the PCA so all labels render in the requested font.
import matplotlib as _mpl_pca
_PCA_PREV_FAMILY = _mpl_pca.rcParams["font.family"]
_mpl_pca.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

for _stem, _figsize, _s_norm, _s_bold, _ring_kind in PCA_VARIANTS:
    fig, ax = plt.subplots(figsize=_figsize)

    if _ring_kind == "univariate":
        ring_set = univ_sp_set
        ring_flag_map = flag_univ_map
        ring_legend = "thick black border + label = total-TE % outlier within its clade (|z| > 2)"
    else:
        ring_set = comp_sp_set
        ring_flag_map = flag_comp_map
        ring_legend = "thick black border + label = TE-composition outlier within its clade (Mahalanobis χ²₀.₉₅)"

    # Build genome size lookup for point sizing (from df_pca, aligned with coords)
    sp2gs = dict(zip(df_pca["species"].values[mask], df_pca["genome_size_mb"].values[mask]))
    gs_values = np.array([sp2gs.get(s, 30) for s in sp_clean])

    for cl in sorted(set(cl_clean)):
        idx = cl_clean == cl
        sps_cl = sp_clean[idx]
        coords_cl = coords[idx]
        gs_cl = gs_values[idx]
        bold_mask = np.array([ring_flag_map.get(s, False) for s in sps_cl])
        if (~bold_mask).any():
            # Size proportional to genome_size_mb (base=20, scale=3 per Mb, capped at 300)
            sizes = np.clip(gs_cl[~bold_mask] * 3, 15, 300)
            ax.scatter(coords_cl[~bold_mask, 0], coords_cl[~bold_mask, 1],
                       color=color_map.get(cl, "grey"), label=cl,
                       s=sizes, alpha=0.78,
                       edgecolor="black", linewidth=0.5)
        if bold_mask.any():
            sizes_b = np.clip(gs_cl[bold_mask] * 3, 15, 300)
            ax.scatter(coords_cl[bold_mask, 0], coords_cl[bold_mask, 1],
                       color=color_map.get(cl, "grey"),
                       s=sizes_b, alpha=0.95,
                       edgecolor="black", linewidth=2.4, zorder=4,
                       label=None if (~bold_mask).any() else cl)

    for j, sp in enumerate(sp_clean):
        if sp in ring_set:
            ax.annotate(_spd(sp, italic=True),
                        (coords[j, 0], coords[j, 1]),
                        fontsize=5.5, fontweight="normal",
                        fontstyle="italic",
                        family=["Arial", "DejaVu Sans"],
                        xytext=(4, 2), textcoords="offset points",
                        color="black",
                        alpha=0.7)

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
    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}% of variance)", fontsize=11)
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}% of variance)", fontsize=11)
    ax.set_title(
        "PCA of TE class composition (proportions of LTR / DNA / LINE / SINE / "
        "Rolling-Circle / Penelope / Unclassified)\n"
        "Arrows = variable loadings  |  % = variance explained by each axis  |  "
        "Point size ∝ genome size  |  "
        + ring_legend,
        fontsize=9)

    legend_patches = [mpatches.Patch(color=color_map.get(cl, "grey"),
                                      label=CLADE_NAMES.get(cl, cl))
                      for cl in clade_order if cl in set(cl_clean)]
    ax.legend(handles=legend_patches, fontsize=7, ncol=2,
              loc="upper right", framealpha=0.7)
    plt.tight_layout()
    plot2 = os.path.join(OUT_DIR, f"{_stem}.pdf")
    fig.savefig(plot2, dpi=150)
    fig.savefig(plot2.replace(".pdf", ".png"), dpi=150)
    plt.close(fig)
    print(f"[OK] PCA → {plot2}")

_mpl_pca.rcParams["font.family"] = _PCA_PREV_FAMILY


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
