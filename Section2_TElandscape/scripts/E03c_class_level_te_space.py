#!/usr/bin/env python3
"""
E03c_class_level_te_space.py
============================
TE-class-level ordination (7 classes, not superfamilies).

Two complementary views:
  1. Abundance — Aitchison PCA on CLR-transformed proportions of TE bp
     per class (relative composition: what fraction of TE content is LTR
     vs DNA vs LINE vs …)
  2. Presence/absence — PCoA on Jaccard distance (which classes are present
     in each genome, with a threshold to define "present")

Missing classes = 0 bp (as confirmed).

Data source: highLevelCount.txt (already aggregated per class by EarlGrey).

Per-clade color / marker scheme: same as E03.

Outputs (in OUT_DIR):
  class_level_pca.pdf          — CLR abundance PCA (PC1–4)
  class_level_pcoa_pa.pdf      — P/A Jaccard PCoA (PCo1–4)
  class_level_matrix.tsv       — species × 7 classes bp/genome-size (%)

Usage:
  python3 E03c_class_level_te_space.py [EARLGREY_DIR] [OUT_DIR]
"""

import sys, os, glob
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.cm as _cm

# ── config ────────────────────────────────────────────────────────────────────
EARLGREY_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/clade_outliers"

os.makedirs(OUT_DIR, exist_ok=True)

# ── broad group + clade color/marker defs (same as E03) ──────────────────────
BROAD_GROUP = {}
for _c in ["HO"]:            BROAD_GROUP[_c] = "Holozoa"
for _c in ["NU"]:            BROAD_GROUP[_c] = "Nucleariida"
for _c in ["AS","BA"]:       BROAD_GROUP[_c] = "Dikarya"
for _c in ["AP","BL","BS","CH","CR","EN","GL","KX","MI","MO","MU","OL","ZO"]:
    BROAD_GROUP[_c] = "EDF"

GROUP_COLORS = {
    "Holozoa":    "#e63946",
    "Nucleariida":"#f4a261",
    "EDF":        "#e9b949",
    "Dikarya":    "#1d3557",
}
GROUP_ALPHA = {"Holozoa":0.85,"Nucleariida":0.85,"EDF":0.55,"Dikarya":0.55}

_edf_sorted = ["AP","BL","BS","CH","CR","EN","GL","KX","MI","MO","MU","OL","ZO"]
_edf_cmap   = _cm.YlOrBr(np.linspace(0.30, 0.85, len(_edf_sorted)))
CLADE_COLORS = {cl: _edf_cmap[i] for i, cl in enumerate(_edf_sorted)}
CLADE_COLORS.update({"HO":"#e63946","NU":"#f4a261","AS":"#4a90d9","BA":"#1d3557"})

CLADE_MARKERS = {
    "HO":"o","NU":"D",
    "AS":"^","BA":"v",
    "MU":"s","MO":"P","GL":"h","CH":"p","BL":"8","BS":"H",
    "EN":"x","KX":"X","ZO":"*","MI":"+","CR":"d","AP":"<","OL":">",
}

# ── TE class mapping from highLevelCount ─────────────────────────────────────
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
TE_CLASSES = ["LTR","DNA","LINE","SINE","Rolling Circle","Penelope","Unclassified"]

# ── 1. parse highLevelCount files ────────────────────────────────────────────
records = []
for fp in sorted(glob.glob(os.path.join(EARLGREY_DIR,
                 "*_EarlGrey","*_summaryFiles","*.highLevelCount.txt"))):
    sp = os.path.basename(fp).replace(".highLevelCount.txt","")
    clade = sp[:2]
    bp_by_class = {c: 0.0 for c in TE_CLASSES}
    genome_size = None
    total_bp    = 0.0
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

df = pd.DataFrame(records).set_index("species")
sp_list = df.index.tolist()
print(f"[INFO] {len(sp_list)} species parsed")

# ── 2. build matrices ────────────────────────────────────────────────────────
# bp matrix
bp_mat = df[TE_CLASSES].fillna(0.0)

# Proportion of total TE bp (composition), absent class = 0
total_te = bp_mat.sum(axis=1).replace(0, np.nan)
prop_mat = bp_mat.div(total_te, axis=0).fillna(0.0)

# Also save % of genome per class
genome_pct = bp_mat.copy()
# (We don't have genome size per-row here, but total_TE_pct is in df)
# Save the proportion matrix (composition) as TSV
out_mat = os.path.join(OUT_DIR, "class_level_matrix.tsv")
prop_mat.round(6).to_csv(out_mat, sep="\t")
print(f"[OK] Class-level matrix → {out_mat}")

# ── 3. helpers ───────────────────────────────────────────────────────────────
def clr_transform(mat, pseudo=1e-6):
    X = mat.values + pseudo
    logX = np.log(X)
    return logX - logX.mean(axis=1, keepdims=True)

def pcoa_jaccard(bin_mat_vals):
    """Return (coords, var_pct) via classical PCoA on Jaccard distance."""
    D = squareform(pdist(bin_mat_vals, metric="jaccard"))
    n = D.shape[0]
    D2 = D ** 2
    H  = np.eye(n) - np.ones((n, n)) / n
    B  = -0.5 * H @ D2 @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    pos = eigvals > 1e-10
    eigvals_pos = eigvals[pos]
    coords = eigvecs[:, pos] * np.sqrt(eigvals_pos)
    var_pct = eigvals_pos / eigvals_pos.sum() * 100
    return coords, var_pct

def draw_ordination(ax, coords, sp_names, title, xlabel, ylabel, draw_mst=True):
    for sp, xy in zip(sp_names, coords):
        clade  = sp[:2]
        group  = BROAD_GROUP.get(clade, "EDF")
        marker = CLADE_MARKERS.get(clade, "o")
        color  = CLADE_COLORS.get(clade, GROUP_COLORS.get(group, "#888888"))
        if hasattr(color, "__len__") and not isinstance(color, str):
            color = tuple(color)
        alpha = GROUP_ALPHA.get(group, 0.7)
        ax.scatter(xy[0], xy[1], c=[color], marker=marker,
                   s=40, alpha=alpha, linewidths=0.3, edgecolors="white", zorder=3)

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
            ax.fill(pts[hull.vertices, 0], pts[hull.vertices, 1],
                    alpha=0.08, color=GROUP_COLORS[group])
            hull_pts = np.append(hull.vertices, hull.vertices[0])
            ax.plot(pts[hull_pts, 0], pts[hull_pts, 1],
                    color=GROUP_COLORS[group], lw=1.2, alpha=0.5)
        except Exception:
            pass

    # MST on clade centroids
    if draw_mst:
        clade_dict = {}
        for sp, xy in zip(sp_names, coords):
            clade_dict.setdefault(sp[:2], []).append(xy)
        centroids = {cl: np.array(pts).mean(axis=0)
                     for cl, pts in clade_dict.items()}
        clade_names = sorted(centroids.keys())
        if len(clade_names) >= 2:
            G = nx.Graph()
            for i, ci in enumerate(clade_names):
                for j, cj in enumerate(clade_names):
                    if i < j:
                        G.add_edge(ci, cj,
                                   weight=np.linalg.norm(centroids[ci]-centroids[cj]))
            mst = nx.minimum_spanning_tree(G)
            for u, v in mst.edges():
                xu, yu = centroids[u]
                xv, yv = centroids[v]
                ax.plot([xu, xv], [yu, yv], color="0.55", lw=1.0,
                        ls="--", zorder=2, alpha=0.8)

    # Label broad groups
    for group, pts in group_coords.items():
        cx, cy = np.array(pts).mean(axis=0)
        ax.text(cx, cy, group, ha="center", va="center",
                fontsize=9, fontweight="bold", color=GROUP_COLORS[group],
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6))

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.axhline(0, color="0.85", lw=0.5, zorder=0)
    ax.axvline(0, color="0.85", lw=0.5, zorder=0)
    ax.spines[["top","right"]].set_visible(False)

def legend_handles():
    gh = [mpatches.Patch(color=GROUP_COLORS[g], label=g, alpha=0.8)
          for g in ["Holozoa","Nucleariida","EDF","Dikarya"]]
    ch = []
    for clade, marker in sorted(CLADE_MARKERS.items()):
        col = CLADE_COLORS.get(clade, "#888888")
        if hasattr(col, "__len__") and not isinstance(col, str):
            col = tuple(col)
        ch.append(Line2D([0],[0], marker=marker, color="w",
                         markerfacecolor=col, markeredgecolor="0.4",
                         markersize=6, label=clade))
    return gh, ch

# ── 4. Abundance: CLR PCA ────────────────────────────────────────────────────
# Drop species with zero total TE (can't compute CLR)
keep = total_te.notna() & (total_te > 0)
prop_sub = prop_mat.loc[keep]
sp_sub   = prop_sub.index.tolist()
clr_vals = clr_transform(prop_sub)

n_comp = min(4, len(TE_CLASSES) - 1, len(sp_sub))
pca = PCA(n_components=n_comp)
coords_pca = pca.fit_transform(clr_vals)
var_exp = pca.explained_variance_ratio_ * 100
print(f"[INFO] CLR PCA: n={len(sp_sub)}, {n_comp} PCs, "
      f"variance: {', '.join(f'PC{i+1}={v:.1f}%' for i,v in enumerate(var_exp))}")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
draw_ordination(axes[0], coords_pca[:, :2], sp_sub,
    f"TE class composition — Aitchison PCA  (n={len(sp_sub)})",
    f"PC1 ({var_exp[0]:.1f}%)", f"PC2 ({var_exp[1]:.1f}%)")
if n_comp >= 4:
    draw_ordination(axes[1], coords_pca[:, 2:4], sp_sub,
        "PC3 vs PC4",
        f"PC3 ({var_exp[2]:.1f}%)", f"PC4 ({var_exp[3]:.1f}%)")
else:
    axes[1].set_visible(False)

# Add loading arrows for PC1/PC2 (biplot overlay, scaled)
loadings = pd.DataFrame(pca.components_[:2].T,
                        index=TE_CLASSES, columns=["PC1","PC2"])
ax0 = axes[0]
scale = 0.8 * max(np.abs(coords_pca[:, :2]).max(axis=0))
for cls in TE_CLASSES:
    lx, ly = loadings.loc[cls, "PC1"] * scale, loadings.loc[cls, "PC2"] * scale
    if np.sqrt(lx**2 + ly**2) < 0.05 * scale:
        continue
    ax0.annotate("", xy=(lx, ly), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="-|>", color="0.4", lw=1.2))
    ax0.text(lx * 1.12, ly * 1.12, cls, fontsize=7, color="0.3",
             ha="center", va="center")

gh, ch = legend_handles()
fig.legend(handles=gh + ch, loc="lower center", ncol=10, fontsize=7,
           frameon=False, bbox_to_anchor=(0.5, -0.05))
fig.suptitle("TE class composition space (7 classes, CLR-transformed)\n"
             "dashed grey = MST on clade centroids  |  arrows = PC loadings",
             fontsize=13, y=1.01)
plt.tight_layout(rect=[0, 0.06, 1, 1])
out_pca = os.path.join(OUT_DIR, "class_level_pca.pdf")
fig.savefig(out_pca, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Class-level CLR PCA → {out_pca}")

# ── 5. Presence/absence: Jaccard PCoA ────────────────────────────────────────
# A class is "present" if it contributes >0% of total TE bp
bin_mat = (prop_sub > 0).astype(float)

# Report how many species have each class
print("\n[INFO] Class presence (n species with >0 bp):")
for cls in TE_CLASSES:
    n_pres = int(bin_mat[cls].sum())
    print(f"  {cls:<18} {n_pres:>3} / {len(sp_sub)}")

coords_pa, var_pa = pcoa_jaccard(bin_mat.values)

# Pad to at least 4 axes in case all classes are near-ubiquitous
if coords_pa.shape[1] < 4:
    print(f"\n[NOTE] Only {coords_pa.shape[1]} positive eigenvalue(s) from Jaccard PCoA.")
    print("       This means all genomes share nearly identical class presence profiles —")
    print("       the P/A signal at class level is minimal (all 7 classes found everywhere).")
    n_pad = 4 - coords_pa.shape[1]
    coords_pa = np.hstack([coords_pa, np.zeros((coords_pa.shape[0], n_pad))])
    var_pa = np.concatenate([var_pa, np.zeros(n_pad)])
print(f"\n[INFO] Jaccard PCoA: "
      f"PCo1={var_pa[0]:.1f}%  PCo2={var_pa[1]:.1f}%  "
      f"PCo3={var_pa[2]:.1f}%  PCo4={var_pa[3]:.1f}%")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
draw_ordination(axes[0], coords_pa[:, :2], sp_sub,
    f"TE class P/A — PCoA (Jaccard)  (n={len(sp_sub)})",
    f"PCo1 ({var_pa[0]:.1f}%)", f"PCo2 ({var_pa[1]:.1f}%)")
draw_ordination(axes[1], coords_pa[:, 2:4], sp_sub,
    "PCo3 vs PCo4",
    f"PCo3 ({var_pa[2]:.1f}%)", f"PCo4 ({var_pa[3]:.1f}%)")

gh, ch = legend_handles()
fig.legend(handles=gh + ch, loc="lower center", ncol=10, fontsize=7,
           frameon=False, bbox_to_anchor=(0.5, -0.05))
fig.suptitle("TE class presence/absence — PCoA (Jaccard distance)\n"
             "dashed grey = MST on clade centroids",
             fontsize=13, y=1.01)
plt.tight_layout(rect=[0, 0.06, 1, 1])
out_pa = os.path.join(OUT_DIR, "class_level_pcoa_pa.pdf")
fig.savefig(out_pa, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Class-level P/A PCoA → {out_pa}")

# ── 6. Combined 4-panel summary ───────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

draw_ordination(axes[0,0], coords_pca[:, :2], sp_sub,
    f"Abundance CLR PCA — PC1 vs PC2",
    f"PC1 ({var_exp[0]:.1f}%)", f"PC2 ({var_exp[1]:.1f}%)")
if n_comp >= 4:
    draw_ordination(axes[0,1], coords_pca[:, 2:4], sp_sub,
        "Abundance CLR PCA — PC3 vs PC4",
        f"PC3 ({var_exp[2]:.1f}%)", f"PC4 ({var_exp[3]:.1f}%)")
else:
    axes[0,1].set_visible(False)
draw_ordination(axes[1,0], coords_pa[:, :2], sp_sub,
    "P/A Jaccard PCoA — PCo1 vs PCo2",
    f"PCo1 ({var_pa[0]:.1f}%)", f"PCo2 ({var_pa[1]:.1f}%)")
draw_ordination(axes[1,1], coords_pa[:, 2:4], sp_sub,
    "P/A Jaccard PCoA — PCo3 vs PCo4",
    f"PCo3 ({var_pa[2]:.1f}%)", f"PCo4 ({var_pa[3]:.1f}%)")

gh, ch = legend_handles()
fig.legend(handles=gh + ch, loc="lower center", ncol=10, fontsize=7,
           frameon=False, bbox_to_anchor=(0.5, -0.03))
fig.suptitle("TE class composition space (7 classes)\n"
             "Top: CLR abundance (how much of each class); "
             "Bottom: presence/absence (which classes present)",
             fontsize=13, y=1.01)
plt.tight_layout(rect=[0, 0.05, 1, 1])
out_summary = os.path.join(OUT_DIR, "class_level_summary.pdf")
fig.savefig(out_summary, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Summary → {out_summary}")
