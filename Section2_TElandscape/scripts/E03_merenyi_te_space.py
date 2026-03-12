#!/usr/bin/env python3
"""
E03_merenyi_te_space.py
=======================
Ordination of all 238 genomes in TE superfamily composition space —
inspired by Merényi et al. 2023, Fig 1b (gene content PCoA).

Here we ask: do Early Diverging Fungi (EDF) occupy a distinct TE space from
Dikarya (Ascomycota + Basidiomycota) and from the Holozoa / Nucleariida
outgroups?

Methods:
  - Parse familyLevelCount.txt to get bp per TE superfamily (NAME#Class/SF)
  - Build species × superfamily proportion matrix (CLR-transformed)
  - PCoA on Aitchison (Euclidean after CLR) distance  →  main ordination
  - UMAP on same matrix  →  non-linear complement
  - Convex hulls per broad group to highlight clustering

Broad group assignment (by species prefix):
  Holozoa:    HO
  Nucleariida: NU
  EDF:        AP, BL, BS, CH, CR, EN, GL, KX, MI, MO, MU, OL, ZO
  Dikarya:    AS, BA

Outputs (in OUT_DIR):
  merenyi_pcoa.pdf   — PCoA PC1/PC2 with convex hulls + per-clade shapes
  merenyi_umap.pdf   — UMAP with same aesthetics
  te_space_matrix.tsv — species × SF proportion matrix (for reproducibility)

Usage:
  python3 E03_merenyi_te_space.py [EARLGREY_DIR] [OUT_DIR]
"""

import sys, os, glob, warnings, re
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ── config ────────────────────────────────────────────────────────────────────
EARLGREY_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/clade_outliers"

os.makedirs(OUT_DIR, exist_ok=True)

# ── broad group mapping ───────────────────────────────────────────────────────
BROAD_GROUP = {}
for code in ["HO"]:
    BROAD_GROUP[code] = "Holozoa"
for code in ["NU"]:
    BROAD_GROUP[code] = "Nucleariida"
for code in ["AS", "BA"]:
    BROAD_GROUP[code] = "Dikarya"
for code in ["AP","BL","BS","CH","CR","EN","GL","KX","MI","MO","MU","OL","ZO"]:
    BROAD_GROUP[code] = "EDF"

GROUP_COLORS = {
    "Holozoa":    "#e63946",
    "Nucleariida":"#f4a261",
    "EDF":        "#2a9d8f",
    "Dikarya":    "#457b9d",
}
GROUP_ALPHA = {
    "Holozoa": 0.85, "Nucleariida": 0.85,
    "EDF": 0.55,     "Dikarya": 0.55,
}

# Per-clade marker shapes (so you can still distinguish clades within groups)
CLADE_MARKERS = {
    "HO":"o","NU":"D",
    "AS":"^","BA":"v",
    "MU":"s","MO":"P","GL":"h","CH":"p","BL":"8","BS":"H",
    "EN":"x","KX":"X","ZO":"*","MI":"+","CR":"d","AP":"<","OL":">",
}

# ── 1. parse familyLevelCount → superfamily matrix ────────────────────────────
def extract_superfamily(name_str):
    """Extract superfamily from 'NAME#Class/Superfamily' or 'NAME#Class'."""
    if "#" not in name_str:
        return "Unknown"
    after_hash = name_str.split("#", 1)[1]
    if "/" in after_hash:
        return after_hash.split("/", 1)[1].strip()
    return after_hash.strip()

records = {}   # sp → {sf: bp}
genome_sizes = {}

# We need genome size — get it from highLevelCount
for hl_fp in sorted(glob.glob(os.path.join(EARLGREY_DIR,
                    "*_EarlGrey","*_summaryFiles","*.highLevelCount.txt"))):
    sp = os.path.basename(hl_fp).replace(".highLevelCount.txt","")
    try:
        df = pd.read_csv(hl_fp, sep="\t")
        df.columns = [c.strip() for c in df.columns]
        size_col = [c for c in df.columns if "Genome Size" in c][0]
        gs = pd.to_numeric(df[size_col], errors="coerce").dropna().iloc[0]
        genome_sizes[sp] = float(gs)
    except Exception:
        pass

# Parse family level counts
fl_files = sorted(glob.glob(os.path.join(EARLGREY_DIR,
                  "*_EarlGrey","*_summaryFiles","*.familyLevelCount.txt")))
print(f"[INFO] Found {len(fl_files)} familyLevelCount files")

for fp in fl_files:
    sp = os.path.basename(fp).replace(".familyLevelCount.txt","")
    sf_bp = {}
    try:
        df = pd.read_csv(fp, sep="\t")
        df.columns = [c.strip() for c in df.columns]
        fam_col = [c for c in df.columns if "Family" in c or "TE" in c][0]
        bp_col  = [c for c in df.columns if "Coverage" in c and "bp" in c][0]
        for _, r in df.iterrows():
            sf = extract_superfamily(str(r[fam_col]).strip())
            bp = float(r[bp_col]) if str(r[bp_col]) not in ("NA","nan","") else 0.0
            sf_bp[sf] = sf_bp.get(sf, 0.0) + bp
        records[sp] = sf_bp
    except Exception as e:
        print(f"[WARN] {sp}: {e}")

print(f"[INFO] Parsed {len(records)} species")

# Build matrix
all_sf = sorted({sf for d in records.values() for sf in d})
sp_list = sorted(records.keys())
mat = pd.DataFrame(
    [[records[sp].get(sf, 0.0) for sf in all_sf] for sp in sp_list],
    index=sp_list, columns=all_sf
)

# Drop Unknown/unclassified superfamilies (too noisy)
keep_sf = [sf for sf in all_sf if sf.lower() not in
           ("unknown","unclassified","simple_repeat","microsatellite","rna","low_complexity","satellite")]
mat = mat[keep_sf]

# Drop superfamilies present in < 5% of species
min_presence = int(len(sp_list) * 0.05)
mat = mat.loc[:, (mat > 0).sum() >= min_presence]
print(f"[INFO] Matrix: {mat.shape[0]} species × {mat.shape[1]} superfamilies (after filtering)")

# Convert to proportions (fraction of TE bp per species)
row_sums = mat.sum(axis=1).replace(0, np.nan)
prop_mat = mat.div(row_sums, axis=0).fillna(0)

# ── 2. CLR transform (Aitchison geometry) ───────────────────────────────────
# Replace zeros with small pseudocount before CLR
pseudo = 1e-6
clr_mat = prop_mat.copy()
clr_mat = clr_mat + pseudo
# CLR: log(x_i) - mean(log(x))
log_mat = np.log(clr_mat.values)
clr_vals = log_mat - log_mat.mean(axis=1, keepdims=True)
clr_df = pd.DataFrame(clr_vals, index=prop_mat.index, columns=prop_mat.columns)

# Save proportions
prop_out = os.path.join(OUT_DIR, "te_space_matrix.tsv")
prop_mat.round(6).to_csv(prop_out, sep="\t")
print(f"[OK] Matrix saved → {prop_out}")

# ── 3. PCoA via PCA on CLR (= Aitchison PCA) ─────────────────────────────────
pca = PCA(n_components=4)
coords_pca = pca.fit_transform(clr_df.values)
var_exp = pca.explained_variance_ratio_ * 100

# ── 4. UMAP ───────────────────────────────────────────────────────────────────
try:
    import umap
    reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.3,
                        metric="euclidean", random_state=42)
    coords_umap = reducer.fit_transform(clr_df.values)
    has_umap = True
except ImportError:
    has_umap = False
    print("[WARN] umap-learn not available, skipping UMAP")

# ── 5. helper: draw one ordination panel ──────────────────────────────────────
def draw_ordination(ax, coords, sp_names, title, xlabel, ylabel):
    groups_present = []
    for sp, xy in zip(sp_names, coords):
        clade = sp[:2]
        group = BROAD_GROUP.get(clade, "EDF")
        if group not in groups_present:
            groups_present.append(group)
        marker = CLADE_MARKERS.get(clade, "o")
        color  = GROUP_COLORS[group]
        alpha  = GROUP_ALPHA[group]
        ax.scatter(xy[0], xy[1], c=color, marker=marker,
                   s=35, alpha=alpha, linewidths=0.3, edgecolors="white", zorder=3)

    # Convex hulls per broad group
    group_coords = {}
    for sp, xy in zip(sp_names, coords):
        g = BROAD_GROUP.get(sp[:2], "EDF")
        group_coords.setdefault(g, []).append(xy)

    for group, pts in group_coords.items():
        pts = np.array(pts)
        if len(pts) < 3:
            continue
        try:
            hull = ConvexHull(pts)
            hull_pts = np.append(hull.vertices, hull.vertices[0])
            ax.fill(pts[hull.vertices, 0], pts[hull.vertices, 1],
                    alpha=0.10, color=GROUP_COLORS[group])
            ax.plot(pts[hull_pts, 0], pts[hull_pts, 1],
                    color=GROUP_COLORS[group], lw=1.5, alpha=0.7)
        except Exception:
            pass

    # Label centroid of each group
    for group, pts in group_coords.items():
        pts = np.array(pts)
        cx, cy = pts.mean(axis=0)
        ax.text(cx, cy, group, ha="center", va="center",
                fontsize=9, fontweight="bold",
                color=GROUP_COLORS[group],
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6))

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.axhline(0, color="0.85", lw=0.5, zorder=0)
    ax.axvline(0, color="0.85", lw=0.5, zorder=0)
    ax.spines[["top","right"]].set_visible(False)

# ── Per-clade secondary legend ────────────────────────────────────────────────
def clade_legend_handles():
    handles = []
    for clade, marker in sorted(CLADE_MARKERS.items()):
        group = BROAD_GROUP.get(clade, "EDF")
        h = Line2D([0],[0], marker=marker, color="w",
                   markerfacecolor=GROUP_COLORS[group],
                   markeredgecolor="0.4", markersize=6, label=clade)
        handles.append(h)
    return handles

# ── 6. PCoA figure ────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

draw_ordination(axes[0], coords_pca[:, :2], sp_list,
    f"Aitchison PCA — PC1 vs PC2",
    f"PC1 ({var_exp[0]:.1f}%)",
    f"PC2 ({var_exp[1]:.1f}%)")

draw_ordination(axes[1], coords_pca[:, 2:4], sp_list,
    f"Aitchison PCA — PC3 vs PC4",
    f"PC3 ({var_exp[2]:.1f}%)",
    f"PC4 ({var_exp[3]:.1f}%)")

# Group colour legend
group_handles = [mpatches.Patch(color=GROUP_COLORS[g], label=g, alpha=0.8)
                 for g in ["Holozoa","Nucleariida","EDF","Dikarya"]]
fig.legend(handles=group_handles + clade_legend_handles(),
           loc="lower center", ncol=10, fontsize=7,
           frameon=False, bbox_to_anchor=(0.5, -0.05))

fig.suptitle("TE superfamily composition space across Holomycota\n"
             "(Aitchison PCA on CLR-transformed proportions)",
             fontsize=13, y=1.01)
plt.tight_layout(rect=[0, 0.06, 1, 1])

out_pcoa = os.path.join(OUT_DIR, "merenyi_pcoa.pdf")
fig.savefig(out_pcoa, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Aitchison PCA → {out_pcoa}")

# ── 7. UMAP figure ────────────────────────────────────────────────────────────
if has_umap:
    fig, ax = plt.subplots(figsize=(9, 8))
    draw_ordination(ax, coords_umap, sp_list,
        "UMAP — TE superfamily composition space",
        "UMAP 1", "UMAP 2")
    group_handles = [mpatches.Patch(color=GROUP_COLORS[g], label=g, alpha=0.8)
                     for g in ["Holozoa","Nucleariida","EDF","Dikarya"]]
    ax.legend(handles=group_handles + clade_legend_handles(),
              loc="upper right", ncol=2, fontsize=7, framealpha=0.8)
    plt.tight_layout()
    out_umap = os.path.join(OUT_DIR, "merenyi_umap.pdf")
    fig.savefig(out_umap, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[OK] UMAP → {out_umap}")

# ── 8. Print variance explained and top loading SFs ──────────────────────────
print("\n[INFO] Variance explained:")
for i, v in enumerate(var_exp[:4], 1):
    print(f"  PC{i}: {v:.1f}%")

loadings = pd.DataFrame(pca.components_.T, index=clr_df.columns,
                        columns=[f"PC{i}" for i in range(1, 5)])
print("\n[INFO] Top 5 superfamilies driving PC1:")
print(loadings["PC1"].abs().nlargest(5).to_string())
print("\n[INFO] Top 5 superfamilies driving PC2:")
print(loadings["PC2"].abs().nlargest(5).to_string())
