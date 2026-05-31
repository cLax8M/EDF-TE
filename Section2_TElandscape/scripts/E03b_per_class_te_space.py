#!/usr/bin/env python3
"""
E03b_per_class_te_space.py
==========================
Per-TE-class ordination — one Aitchison PCA per TE class.

For each of the 7 main TE classes (LTR, DNA, LINE, SINE, Rolling Circle,
Penelope, Unclassified) we:
  1. Filter the familyLevelCount matrix to superfamilies of that class only
  2. CLR-transform (with pseudocount)
  3. Run Aitchison PCA (PC1 vs PC2)
  4. Draw scatter + convex hulls + MST on clade centroids

Objective: determine whether EDF occupies an intermediate TE space for each
class, and if the intermediate position is class-specific.

Per-clade color / shape scheme: same as E03.

Usage:
  python3 E03b_per_class_te_space.py [EARLGREY_DIR] [OUT_DIR]
  python3 E03b_per_class_te_space.py  # uses defaults
"""

import sys, os, glob, warnings
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist, squareform
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

# ── TE class → superfamily keyword mapping ───────────────────────────────────
# Keys are how the class appears after '#' in familyLevelCount NAME#Class/SF
# Values are canonical class labels we report
CLASS_KEYWORD_MAP = {
    "LTR":              "LTR",
    "LINE":             "LINE",
    "SINE":             "SINE",
    "DNA":              "DNA",
    "DNA-nested":       "DNA",
    "Rolling Circle":   "Rolling Circle",
    "RC":               "Rolling Circle",
    "Penelope":         "Penelope",
    "PLE":              "Penelope",
    "Unclassified":     "Unclassified",
}

TE_CLASSES_ORDER = ["LTR","DNA","LINE","SINE","Rolling Circle","Penelope","Unclassified"]

# ── 1. parse familyLevelCount files ─────────────────────────────────────────
def extract_class_and_sf(name_str):
    """Return (class_label, superfamily) from 'NAME#Class/SF' or 'NAME#Class'."""
    if "#" not in name_str:
        return None, "Unknown"
    after = name_str.split("#", 1)[1]
    if "/" in after:
        cls_raw, sf = after.split("/", 1)
    else:
        cls_raw, sf = after.strip(), after.strip()
    cls_raw = cls_raw.strip()
    sf      = sf.strip()
    # Normalise nested variants
    base_cls = cls_raw.replace("-nested","").strip()
    mapped = CLASS_KEYWORD_MAP.get(cls_raw) or CLASS_KEYWORD_MAP.get(base_cls)
    return mapped, sf

print("[INFO] Scanning familyLevelCount files …")
fl_files = sorted(glob.glob(os.path.join(EARLGREY_DIR,
                  "*_EarlGrey","*_summaryFiles","*.familyLevelCount.txt")))
print(f"[INFO] Found {len(fl_files)} files")

# sp → {(class_label, sf): bp}
records_raw = {}
for fp in fl_files:
    sp = os.path.basename(fp).replace(".familyLevelCount.txt","")
    try:
        df = pd.read_csv(fp, sep="\t")
        df.columns = [c.strip() for c in df.columns]
        fam_col = [c for c in df.columns if "Family" in c or "TE" in c][0]
        bp_col  = [c for c in df.columns if "Coverage" in c and "bp" in c][0]
        d = {}
        for _, r in df.iterrows():
            cls_lbl, sf = extract_class_and_sf(str(r[fam_col]).strip())
            if cls_lbl is None:
                continue
            bp = float(r[bp_col]) if str(r[bp_col]) not in ("NA","nan","") else 0.0
            key = (cls_lbl, sf)
            d[key] = d.get(key, 0.0) + bp
        records_raw[sp] = d
    except Exception as e:
        print(f"[WARN] {sp}: {e}")

sp_list = sorted(records_raw.keys())
print(f"[INFO] Parsed {len(sp_list)} species")

# ── 2. build one matrix per TE class ────────────────────────────────────────
# For each class, columns = superfamilies from that class only
class_matrices = {}
class_bin_matrices = {}
for cls_label in TE_CLASSES_ORDER:
    sfs = sorted({sf for sp in sp_list
                  for (cl, sf) in records_raw[sp].keys() if cl == cls_label})
    if len(sfs) < 2:
        print(f"[SKIP] {cls_label}: only {len(sfs)} superfamilies")
        continue
    mat = pd.DataFrame(
        [[records_raw[sp].get((cls_label, sf), 0.0) for sf in sfs] for sp in sp_list],
        index=sp_list, columns=sfs
    )
    # Filter SFs present in < 5% of species
    min_pres = max(2, int(len(sp_list) * 0.05))
    mat = mat.loc[:, (mat > 0).sum() >= min_pres]
    if mat.shape[1] < 2:
        print(f"[SKIP] {cls_label}: only {mat.shape[1]} SFs after filtering presence threshold")
        continue
    # Proportion within that class per species
    row_sums = mat.sum(axis=1).replace(0, np.nan)
    prop = mat.div(row_sums, axis=0).fillna(0)
    bin_mat = (mat > 0).astype(float)
    class_matrices[cls_label] = prop
    class_bin_matrices[cls_label] = bin_mat
    print(f"[INFO] {cls_label}: {prop.shape[1]} superfamilies retained")

# ── 3. drawing function ──────────────────────────────────────────────────────
def clr_transform(prop_df, pseudo=1e-6):
    X = prop_df.values + pseudo
    logX = np.log(X)
    return logX - logX.mean(axis=1, keepdims=True)

def pcoa_jaccard(bin_df):
    """Classical PCoA from Jaccard distance matrix on presence/absence data."""
    D = squareform(pdist(bin_df.values, metric="jaccard"))
    n = D.shape[0]
    D2 = D ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ D2 @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pos = eigvals > 1e-10
    eigvals_pos = eigvals[pos]
    if eigvals_pos.size == 0:
        coords = np.zeros((n, 1))
        var = np.array([0.0])
        return coords, var
    coords = eigvecs[:, pos] * np.sqrt(eigvals_pos)
    var = eigvals_pos / eigvals_pos.sum() * 100
    return coords, var

def draw_ordination_class(ax, coords, sp_names, title, xlabel, ylabel):
    for sp, xy in zip(sp_names, coords):
        clade  = sp[:2]
        group  = BROAD_GROUP.get(clade, "EDF")
        marker = CLADE_MARKERS.get(clade, "o")
        color  = CLADE_COLORS.get(clade, GROUP_COLORS.get(group, "#888888"))
        if hasattr(color, "__len__") and not isinstance(color, str):
            color = tuple(color)
        alpha  = GROUP_ALPHA.get(group, 0.7)
        ax.scatter(xy[0], xy[1], c=[color], marker=marker,
                   s=30, alpha=alpha, linewidths=0.3, edgecolors="white", zorder=3)

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
                    alpha=0.07, color=GROUP_COLORS[group])
            hull_pts = np.append(hull.vertices, hull.vertices[0])
            ax.plot(pts[hull_pts, 0], pts[hull_pts, 1],
                    color=GROUP_COLORS[group], lw=1.0, alpha=0.4)
        except Exception:
            pass

    # MST on clade centroids
    clade_coord_map = {}
    for sp, xy in zip(sp_names, coords):
        cl = sp[:2]
        clade_coord_map.setdefault(cl, []).append(xy)
    centroids = {cl: np.array(pts).mean(axis=0)
                 for cl, pts in clade_coord_map.items() if len(pts) >= 1}
    clade_names = sorted(centroids.keys())
    if len(clade_names) >= 2:
        G = nx.Graph()
        for i, ci in enumerate(clade_names):
            for j, cj in enumerate(clade_names):
                if i < j:
                    G.add_edge(ci, cj,
                               weight=np.linalg.norm(centroids[ci] - centroids[cj]))
        mst = nx.minimum_spanning_tree(G)
        for u, v in mst.edges():
            xu, yu = centroids[u]
            xv, yv = centroids[v]
            ax.plot([xu, xv], [yu, yv], color="0.55", lw=0.8,
                    ls="--", zorder=2, alpha=0.8)

    # Label broad groups at centroid
    for group, pts in group_coords.items():
        pts = np.array(pts)
        cx, cy = pts.mean(axis=0)
        ax.text(cx, cy, group, ha="center", va="center",
                fontsize=8, fontweight="bold",
                color=GROUP_COLORS[group],
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.6))

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.axhline(0, color="0.85", lw=0.5, zorder=0)
    ax.axvline(0, color="0.85", lw=0.5, zorder=0)
    ax.spines[["top","right"]].set_visible(False)

# ── 4. legend helpers ────────────────────────────────────────────────────────
def make_legend_handles():
    group_handles = [mpatches.Patch(color=GROUP_COLORS[g], label=g, alpha=0.8)
                     for g in ["Holozoa","Nucleariida","EDF","Dikarya"]]
    clade_handles = []
    for clade, marker in sorted(CLADE_MARKERS.items()):
        col = CLADE_COLORS.get(clade, "#888888")
        if hasattr(col, "__len__") and not isinstance(col, str):
            col = tuple(col)
        clade_handles.append(
            Line2D([0],[0], marker=marker, color="w",
                   markerfacecolor=col, markeredgecolor="0.4",
                   markersize=5, label=clade))
    return group_handles, clade_handles

# ── 5. one figure per TE class ───────────────────────────────────────────────
for cls_label, prop_df in class_matrices.items():
    clr_vals = clr_transform(prop_df)
    # Drop species with all-zero CLR (no TEs of this class)
    row_var = clr_vals.var(axis=1)
    keep_rows = row_var > 0
    clr_sub  = clr_vals[keep_rows]
    sp_sub   = prop_df.index[keep_rows].tolist()
    n_sp = len(sp_sub)
    if n_sp < 5:
        print(f"[SKIP] {cls_label}: only {n_sp} species with non-zero content")
        continue

    n_comp = min(4, clr_sub.shape[1], n_sp)
    pca = PCA(n_components=n_comp)
    coords = pca.fit_transform(clr_sub)
    var_exp = pca.explained_variance_ratio_ * 100

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    draw_ordination_class(
        axes[0], coords[:, :2], sp_sub,
        f"{cls_label} — PC1 vs PC2  (n={n_sp})",
        f"PC1 ({var_exp[0]:.1f}%)", f"PC2 ({var_exp[1]:.1f}%)")

    if n_comp >= 4:
        draw_ordination_class(
            axes[1], coords[:, 2:4], sp_sub,
            f"{cls_label} — PC3 vs PC4",
            f"PC3 ({var_exp[2]:.1f}%)", f"PC4 ({var_exp[3]:.1f}%)")
    else:
        axes[1].set_visible(False)

    group_handles, clade_handles = make_legend_handles()
    fig.legend(handles=group_handles + clade_handles,
               loc="lower center", ncol=10, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.05))
    fig.suptitle(f"TE composition space — {cls_label} elements\n"
                 "(Aitchison PCA on CLR-transformed proportions; "
                 "dashed grey = MST on clade centroids)",
                 fontsize=12, y=1.01)
    plt.tight_layout(rect=[0, 0.06, 1, 1])

    safe_name = cls_label.replace(" ","_").replace("/","_")
    out_fp = os.path.join(OUT_DIR, f"merenyi_class_{safe_name}.pdf")
    fig.savefig(out_fp, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[OK] {cls_label} → {out_fp}")

# ── 5b. one P/A figure per TE class (superfamily presence/absence) ──────────
for cls_label, bin_df in class_bin_matrices.items():
    # Keep species with at least one present superfamily in this class
    keep_rows = bin_df.sum(axis=1) > 0
    bin_sub = bin_df.loc[keep_rows]
    sp_sub = bin_sub.index.tolist()
    n_sp = len(sp_sub)
    if n_sp < 5:
        print(f"[SKIP] {cls_label} P/A: only {n_sp} species with ≥1 SF present")
        continue

    coords_pa, var_pa = pcoa_jaccard(bin_sub)
    if coords_pa.shape[1] < 4:
        pad = 4 - coords_pa.shape[1]
        coords_pa = np.hstack([coords_pa, np.zeros((coords_pa.shape[0], pad))])
        var_pa = np.concatenate([var_pa, np.zeros(pad)])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    draw_ordination_class(
        axes[0], coords_pa[:, :2], sp_sub,
        f"{cls_label} P/A — PCo1 vs PCo2  (n={n_sp})",
        f"PCo1 ({var_pa[0]:.1f}%)", f"PCo2 ({var_pa[1]:.1f}%)")

    draw_ordination_class(
        axes[1], coords_pa[:, 2:4], sp_sub,
        f"{cls_label} P/A — PCo3 vs PCo4",
        f"PCo3 ({var_pa[2]:.1f}%)", f"PCo4 ({var_pa[3]:.1f}%)")

    group_handles, clade_handles = make_legend_handles()
    fig.legend(handles=group_handles + clade_handles,
               loc="lower center", ncol=10, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, -0.05))
    fig.suptitle(f"TE superfamily P/A space — {cls_label} elements\n"
                 "(Jaccard PCoA; dashed grey = MST on clade centroids)",
                 fontsize=12, y=1.01)
    plt.tight_layout(rect=[0, 0.06, 1, 1])

    safe_name = cls_label.replace(" ","_").replace("/","_")
    out_fp = os.path.join(OUT_DIR, f"merenyi_class_{safe_name}_pa.pdf")
    fig.savefig(out_fp, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[OK] {cls_label} P/A → {out_fp}")

# ── 6. multi-panel summary figure (all classes in one PDF) ───────────────────
n_cls = len(class_matrices)
fig_all, axes_all = plt.subplots(2, 4, figsize=(28, 12),
                                  gridspec_kw={"hspace":0.35,"wspace":0.3})
axes_flat = axes_all.flatten()

for idx, (cls_label, prop_df) in enumerate(class_matrices.items()):
    ax = axes_flat[idx]
    clr_vals = clr_transform(prop_df)
    row_var  = clr_vals.var(axis=1)
    keep_rows = row_var > 0
    clr_sub   = clr_vals[keep_rows]
    sp_sub    = prop_df.index[keep_rows].tolist()
    if len(sp_sub) < 5:
        ax.set_visible(False)
        continue
    n_comp = min(2, clr_sub.shape[1], len(sp_sub))
    pca = PCA(n_components=n_comp)
    coords = pca.fit_transform(clr_sub)
    var_exp = pca.explained_variance_ratio_ * 100
    draw_ordination_class(
        ax, coords[:, :2], sp_sub,
        f"{cls_label}  (n={len(sp_sub)})",
        f"PC1 ({var_exp[0]:.1f}%)",
        f"PC2 ({var_exp[1]:.1f}%)" if n_comp > 1 else "PC2")

# Hide unused panels
for idx in range(len(class_matrices), len(axes_flat)):
    axes_flat[idx].set_visible(False)

group_handles, clade_handles = make_legend_handles()
fig_all.legend(handles=group_handles + clade_handles,
               loc="lower center", ncol=12, fontsize=6.5,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
fig_all.suptitle("TE superfamily composition space — by TE class\n"
                  "(Aitchison PCA; dashed grey = MST on clade centroids)",
                  fontsize=14, y=1.01)

out_summary = os.path.join(OUT_DIR, "merenyi_per_class_summary.pdf")
fig_all.savefig(out_summary, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Per-class summary → {out_summary}")

# ── 7. multi-panel summary for per-class P/A (Jaccard) ──────────────────────
fig_pa, axes_pa = plt.subplots(2, 4, figsize=(28, 12),
                               gridspec_kw={"hspace":0.35,"wspace":0.3})
axes_flat = axes_pa.flatten()

for idx, (cls_label, bin_df) in enumerate(class_bin_matrices.items()):
    ax = axes_flat[idx]
    keep_rows = bin_df.sum(axis=1) > 0
    bin_sub = bin_df.loc[keep_rows]
    sp_sub = bin_sub.index.tolist()
    if len(sp_sub) < 5:
        ax.set_visible(False)
        continue
    coords_pa, var_pa = pcoa_jaccard(bin_sub)
    if coords_pa.shape[1] < 2:
        coords_pa = np.hstack([coords_pa, np.zeros((coords_pa.shape[0], 1))])
        var_pa = np.concatenate([var_pa, np.zeros(1)])
    draw_ordination_class(
        ax, coords_pa[:, :2], sp_sub,
        f"{cls_label}  (n={len(sp_sub)})",
        f"PCo1 ({var_pa[0]:.1f}%)",
        f"PCo2 ({var_pa[1]:.1f}%)")

for idx in range(len(class_bin_matrices), len(axes_flat)):
    axes_flat[idx].set_visible(False)

group_handles, clade_handles = make_legend_handles()
fig_pa.legend(handles=group_handles + clade_handles,
              loc="lower center", ncol=12, fontsize=6.5,
              frameon=False, bbox_to_anchor=(0.5, -0.02))
fig_pa.suptitle("TE superfamily presence/absence space — by TE class\n"
                "(Jaccard PCoA; dashed grey = MST on clade centroids)",
                fontsize=14, y=1.01)

out_summary_pa = os.path.join(OUT_DIR, "merenyi_per_class_summary_pa.pdf")
fig_pa.savefig(out_summary_pa, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Per-class P/A summary → {out_summary_pa}")
