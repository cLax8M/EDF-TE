#!/usr/bin/env python3
"""
E03e_te_space_clean_metrics.py
==============================
Clean visualization + robust quantification of TE landscape position.

Data level:
  species x TE superfamily using occupied bp proportions (0 if absent).

Outputs:
  TE_space_global_clean.pdf
  TE_space_centroid_MST.pdf
  TE_space_interpolation_summary.tsv
  TE_space_group_distance_bootstrap.tsv
  TE_space_group_structure_tests.tsv

Methods included:
  1) Bray-Curtis PCoA on occupied-bp genome proportions
  2) Aitchison PCA on CLR-transformed TE superfamily composition
  3) Jaccard PCoA on presence/absence of superfamilies

Usage:
  python3 E03e_te_space_clean_metrics.py [EARLGREY_DIR] [OUT_DIR]
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

EARLGREY_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/clade_outliers"
os.makedirs(OUT_DIR, exist_ok=True)

RNG = np.random.default_rng(42)
N_PERM = 499
N_BOOT = 499

# ---- styles -----------------------------------------------------------------
BROAD_GROUP = {}
for c in ["HO"]:
    BROAD_GROUP[c] = "Holozoa"
for c in ["NU"]:
    BROAD_GROUP[c] = "Nucleariida"
for c in ["AS", "BA"]:
    BROAD_GROUP[c] = "Dikarya"
for c in ["AP", "BL", "BS", "CH", "CR", "EN", "GL", "KX", "MI", "MO", "MU", "OL", "ZO"]:
    BROAD_GROUP[c] = "EDF"

GROUP_COLORS = {
    "Holozoa": "#e63946",
    "Nucleariida": "#f4a261",
    "EDF": "#e9b949",
    "Dikarya": "#1d3557",
}

CLADE_MARKERS = {
    "HO": "o", "NU": "D",
    "AS": "^", "BA": "v",
    "MU": "s", "MO": "P", "GL": "h", "CH": "p", "BL": "8", "BS": "H",
    "EN": "x", "KX": "X", "ZO": "*", "MI": "+", "CR": "d", "AP": "<", "OL": ">",
}

_edf_clades = ["AP", "BL", "BS", "CH", "CR", "EN", "GL", "KX", "MI", "MO", "MU", "OL", "ZO"]
_edf_colors = cm.YlOrBr(np.linspace(0.30, 0.85, len(_edf_clades)))
CLADE_COLORS = {cl: _edf_colors[i] for i, cl in enumerate(_edf_clades)}
CLADE_COLORS.update({"HO": "#e63946", "NU": "#f4a261", "AS": "#4a90d9", "BA": "#1d3557"})

# Override with canonical clade palette (single source of truth)
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clade_palette import CLADE_COLORS as _CANON
CLADE_COLORS.update(_CANON)

# Non-TE repeat classes to exclude from TE-space analysis
NON_TE_CLASSES = frozenset({
    "Simple_repeat", "Low_complexity", "Satellite",
    "rRNA", "tRNA", "snRNA", "scRNA", "sno", "srpRNA",
    "ARTEFACT",
})

# ---- parsing ----------------------------------------------------------------
def extract_class_and_sf(name_str):
    """Return (te_class, superfamily) from NAME#Class/SF or NAME#Class."""
    if "#" not in name_str:
        return "Unknown", "Unknown"
    after = name_str.split("#", 1)[1]
    if "/" in after:
        te_class, sf = after.split("/", 1)
        return te_class.strip(), sf.strip()
    return after.strip(), after.strip()


def read_data(earlgrey_dir):
    genome_size = {}
    hl_files = sorted(glob.glob(os.path.join(earlgrey_dir, "*_EarlGrey", "*_summaryFiles", "*.highLevelCount.txt")))
    for fp in hl_files:
        sp = os.path.basename(fp).replace(".highLevelCount.txt", "")
        try:
            df = pd.read_csv(fp, sep="\t")
            df.columns = [c.strip() for c in df.columns]
            size_col = [c for c in df.columns if "Genome Size" in c][0]
            gs = pd.to_numeric(df[size_col], errors="coerce").dropna().iloc[0]
            genome_size[sp] = float(gs)
        except Exception:
            continue

    rec = {}
    fl_files = sorted(glob.glob(os.path.join(earlgrey_dir, "*_EarlGrey", "*_summaryFiles", "*.familyLevelCount.txt")))
    print(f"[INFO] family files: {len(fl_files)}")
    for fp in fl_files:
        sp = os.path.basename(fp).replace(".familyLevelCount.txt", "")
        if sp not in genome_size:
            continue
        try:
            df = pd.read_csv(fp, sep="\t")
            df.columns = [c.strip() for c in df.columns]
            fam_col = [c for c in df.columns if "Family" in c or "TE" in c][0]
            bp_col = [c for c in df.columns if "Coverage" in c and "bp" in c][0]
            sf_bp = {}
            for _, r in df.iterrows():
                te_class, sf = extract_class_and_sf(str(r[fam_col]).strip())
                if te_class in NON_TE_CLASSES:
                    continue
                bp = float(r[bp_col]) if str(r[bp_col]) not in ("NA", "nan", "") else 0.0
                sf_bp[sf] = sf_bp.get(sf, 0.0) + bp
            rec[sp] = sf_bp
        except Exception as e:
            print(f"[WARN] {sp}: {e}")

    sp_list = sorted(rec.keys())
    sf_all = sorted({sf for d in rec.values() for sf in d})
    mat_bp = pd.DataFrame(
        [[rec[sp].get(sf, 0.0) for sf in sf_all] for sp in sp_list],
        index=sp_list,
        columns=sf_all,
    )
    mat_occ = mat_bp.copy()
    for sp in sp_list:
        gs = genome_size[sp]
        mat_occ.loc[sp] = mat_occ.loc[sp] / gs if gs > 0 else 0.0

    # TE composition (within total TE bp) for CLR space
    row_sums = mat_bp.sum(axis=1).replace(0, np.nan)
    mat_te_comp = mat_bp.div(row_sums, axis=0).fillna(0.0)

    return sp_list, mat_occ, mat_te_comp


# ---- ordination and stats ---------------------------------------------------
def pcoa_from_distance(D):
    n = D.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (D ** 2) @ H
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    pos = evals > 1e-10
    evals_pos = evals[pos]
    if evals_pos.size == 0:
        coords = np.zeros((n, 1))
        var = np.array([0.0])
        return coords, var
    coords = evecs[:, pos] * np.sqrt(evals_pos)
    var = evals_pos / evals_pos.sum() * 100
    return coords, var


def clr_pca(mat_comp):
    X = np.log(mat_comp.values + 1e-6)
    X = X - X.mean(axis=1, keepdims=True)
    pca = PCA(n_components=min(10, X.shape[0], X.shape[1]))
    coords = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100
    return coords, var


def group_labels(species):
    return np.array([BROAD_GROUP.get(s[:2], "EDF") for s in species])


def centroids(coords, labels):
    return {g: coords[labels == g].mean(axis=0) for g in sorted(set(labels))}


def line_projection_metrics(coords, labels, ref_a="Dikarya", ref_b="Holozoa"):
    c = centroids(coords, labels)
    if ref_a not in c or ref_b not in c:
        return pd.DataFrame()
    a = c[ref_a]
    b = c[ref_b]
    ab = b - a
    den = float(ab @ ab)
    rows = []
    for g in sorted(set(labels)):
        p = c[g]
        if den == 0:
            t = np.nan
            perp = np.nan
        else:
            t = float(((p - a) @ ab) / den)
            foot = a + t * ab
            perp = float(np.linalg.norm(p - foot))
        tg = []
        for x in coords[labels == g]:
            if den == 0:
                tg.append(np.nan)
            else:
                tg.append(float(((x - a) @ ab) / den))
        tg = np.array(tg)
        inside = float(np.nanmean((tg >= 0) & (tg <= 1)) * 100)
        rows.append({
            "group": g,
            "t_centroid": t,
            "perp_centroid": perp,
            "inside_segment_species_pct": inside,
            "distance_to_Dikarya_centroid": float(np.linalg.norm(p - a)) if "Dikarya" in c else np.nan,
            "distance_to_Holozoa_centroid": float(np.linalg.norm(p - b)) if "Holozoa" in c else np.nan,
        })
    return pd.DataFrame(rows)


def bootstrap_centroid_distance(coords, labels, g1, g2, n_boot=499):
    idx1 = np.where(labels == g1)[0]
    idx2 = np.where(labels == g2)[0]
    if len(idx1) < 2 or len(idx2) < 2:
        return np.nan, np.nan, np.nan
    dvals = []
    for _ in range(n_boot):
        s1 = RNG.choice(idx1, size=len(idx1), replace=True)
        s2 = RNG.choice(idx2, size=len(idx2), replace=True)
        c1 = coords[s1].mean(axis=0)
        c2 = coords[s2].mean(axis=0)
        dvals.append(float(np.linalg.norm(c1 - c2)))
    dvals = np.array(dvals)
    return float(np.mean(dvals)), float(np.quantile(dvals, 0.025)), float(np.quantile(dvals, 0.975))


def anosim_stat(D, labels, n_perm=499):
    # ANOSIM R statistic from ranked distances
    n = D.shape[0]
    tri = np.triu_indices(n, k=1)
    d = D[tri]
    r = pd.Series(d).rank(method="average").values
    g1 = labels[tri[0]]
    g2 = labels[tri[1]]
    same = g1 == g2
    if same.sum() == 0 or (~same).sum() == 0:
        return np.nan, np.nan
    r_w = r[same].mean()
    r_b = r[~same].mean()
    R = (r_b - r_w) / (len(r) / 2)

    perm = []
    for _ in range(n_perm):
        labp = RNG.permutation(labels)
        g1p = labp[tri[0]]
        g2p = labp[tri[1]]
        samep = g1p == g2p
        if samep.sum() == 0 or (~samep).sum() == 0:
            continue
        r_wp = r[samep].mean()
        r_bp = r[~samep].mean()
        perm.append((r_bp - r_wp) / (len(r) / 2))
    perm = np.array(perm)
    p = float((np.sum(perm >= R) + 1) / (len(perm) + 1))
    return float(R), p


def dispersion_index(coords, labels):
    rows = []
    c = centroids(coords, labels)
    for g in sorted(set(labels)):
        pts = coords[labels == g]
        d = np.linalg.norm(pts - c[g], axis=1)
        rows.append({"group": g, "mean_dist_to_centroid": float(np.mean(d)), "n": len(pts)})
    return pd.DataFrame(rows)


# ---- plotting ----------------------------------------------------------------
def draw_clean_plot(coords_2d, species, labels, var, out_path):
    fig, ax = plt.subplots(figsize=(10, 8))

    # species as muted background
    ax.scatter(coords_2d[:, 0], coords_2d[:, 1], c="0.78", s=16, alpha=0.35, linewidths=0, zorder=1)

    # clade centroids as main layer
    clade_points = {}
    for sp, xy in zip(species, coords_2d):
        clade_points.setdefault(sp[:2], []).append(xy)
    cl_cent = {cl: np.array(v).mean(axis=0) for cl, v in clade_points.items()}

    for cl, xy in sorted(cl_cent.items()):
        group = BROAD_GROUP.get(cl, "EDF")
        marker = CLADE_MARKERS.get(cl, "o")
        color = CLADE_COLORS.get(cl, GROUP_COLORS[group])
        if hasattr(color, "__len__") and not isinstance(color, str):
            color = tuple(color)
        ax.scatter(xy[0], xy[1], c=[color], marker=marker, s=110,
                   edgecolors="black", linewidths=0.5, zorder=3)
        ax.text(xy[0], xy[1], cl, fontsize=8, fontweight="bold", ha="center", va="center", zorder=4)

    # MST on clade centroids
    clades = sorted(cl_cent)
    G = nx.Graph()
    for i, c1 in enumerate(clades):
        for j, c2 in enumerate(clades):
            if i < j:
                G.add_edge(c1, c2, weight=float(np.linalg.norm(cl_cent[c1] - cl_cent[c2])))
    mst = nx.minimum_spanning_tree(G)
    for u, v in mst.edges():
        x1, y1 = cl_cent[u]
        x2, y2 = cl_cent[v]
        ax.plot([x1, x2], [y1, y2], color="0.55", lw=1.0, ls="--", zorder=2)

    # broad group centroids labels
    grp_cent = centroids(coords_2d, labels)
    for g, xy in grp_cent.items():
        ax.text(xy[0], xy[1], g, fontsize=10, fontweight="bold", color=GROUP_COLORS[g],
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8), zorder=5)

    ax.set_title("TE superfamily space (clean view: species background + clade centroids)", fontsize=12)
    ax.set_xlabel(f"Axis 1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"Axis 2 ({var[1]:.1f}%)")
    ax.axhline(0, color="0.9", lw=0.6)
    ax.axvline(0, color="0.9", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close()


def draw_centroid_only(coords_2d, species, labels, var, out_path):
    fig, ax = plt.subplots(figsize=(10, 8))
    clade_points = {}
    for sp, xy in zip(species, coords_2d):
        clade_points.setdefault(sp[:2], []).append(xy)
    cl_cent = {cl: np.array(v).mean(axis=0) for cl, v in clade_points.items()}

    for cl, xy in sorted(cl_cent.items()):
        group = BROAD_GROUP.get(cl, "EDF")
        marker = CLADE_MARKERS.get(cl, "o")
        color = CLADE_COLORS.get(cl, GROUP_COLORS[group])
        if hasattr(color, "__len__") and not isinstance(color, str):
            color = tuple(color)
        ax.scatter(xy[0], xy[1], c=[color], marker=marker, s=120,
                   edgecolors="black", linewidths=0.6, zorder=3)
        ax.text(xy[0], xy[1], f"{cl} (n={len(clade_points[cl])})", fontsize=7,
                ha="left", va="center", zorder=4)

    clades = sorted(cl_cent)
    G = nx.Graph()
    for i, c1 in enumerate(clades):
        for j, c2 in enumerate(clades):
            if i < j:
                G.add_edge(c1, c2, weight=float(np.linalg.norm(cl_cent[c1] - cl_cent[c2])))
    mst = nx.minimum_spanning_tree(G)
    for u, v in mst.edges():
        x1, y1 = cl_cent[u]
        x2, y2 = cl_cent[v]
        ax.plot([x1, x2], [y1, y2], color="0.55", lw=1.0, ls="--", zorder=1)

    ax.set_title("TE superfamily space (clade centroids + MST)", fontsize=12)
    ax.set_xlabel(f"Axis 1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"Axis 2 ({var[1]:.1f}%)")
    ax.axhline(0, color="0.9", lw=0.6)
    ax.axvline(0, color="0.9", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close()


# ---- run --------------------------------------------------------------------
species, mat_occ, mat_te_comp = read_data(EARLGREY_DIR)
labels = group_labels(species)
print(f"[INFO] species={len(species)}  superfamilies={mat_occ.shape[1]}")

# Drop species with zero TE bp (Bray-Curtis undefined for all-zero rows)
_keep = mat_occ.sum(axis=1) > 0
_dropped = mat_occ.index[~_keep].tolist()
if _dropped:
    print(f"[INFO] Dropping {len(_dropped)} zero-TE species: {_dropped}")
    mat_occ = mat_occ.loc[_keep]
    mat_te_comp = mat_te_comp.loc[_keep]
    species = mat_occ.index.tolist()
    labels = group_labels(species)

# spaces
D_bray = squareform(pdist(mat_occ.values, metric="braycurtis"))
coords_bray, var_bray = pcoa_from_distance(D_bray)

coords_ai, var_ai = clr_pca(mat_te_comp)
D_ai = squareform(pdist(coords_ai, metric="euclidean"))

mat_bin = (mat_occ > 0).astype(float)
D_jac = squareform(pdist(mat_bin.values, metric="jaccard"))
coords_jac, var_jac = pcoa_from_distance(D_jac)

# clean plots from the main (Bray) space
out_clean = os.path.join(OUT_DIR, "TE_space_global_clean.pdf")
out_cent = os.path.join(OUT_DIR, "TE_space_centroid_MST.pdf")
draw_clean_plot(coords_bray[:, :2], species, labels, var_bray, out_clean)
draw_centroid_only(coords_bray[:, :2], species, labels, var_bray, out_cent)
print(f"[OK] {out_clean}")
print(f"[OK] {out_cent}")

# interpolation summaries
interp_rows = []
for space_name, coords in [
    ("BrayOcc", coords_bray),
    ("AitchisonCLR", coords_ai),
    ("JaccardPA", coords_jac),
]:
    dfi = line_projection_metrics(coords, labels)
    dfi.insert(0, "space", space_name)
    interp_rows.append(dfi)

interp = pd.concat(interp_rows, ignore_index=True)
out_interp = os.path.join(OUT_DIR, "TE_space_interpolation_summary.tsv")
interp.to_csv(out_interp, sep="\t", index=False)
print(f"[OK] {out_interp}")

# bootstrap centroid distances
pairs = [("EDF", "Dikarya"), ("EDF", "Holozoa"), ("Nucleariida", "Dikarya"), ("Nucleariida", "Holozoa")]
boot_rows = []
for space_name, coords in [
    ("BrayOcc", coords_bray),
    ("AitchisonCLR", coords_ai),
    ("JaccardPA", coords_jac),
]:
    for g1, g2 in pairs:
        mean_d, lo, hi = bootstrap_centroid_distance(coords, labels, g1, g2, n_boot=N_BOOT)
        boot_rows.append({
            "space": space_name,
            "group1": g1,
            "group2": g2,
            "mean_centroid_distance": mean_d,
            "ci95_low": lo,
            "ci95_high": hi,
            "n_boot": N_BOOT,
        })

boot_df = pd.DataFrame(boot_rows)
out_boot = os.path.join(OUT_DIR, "TE_space_group_distance_bootstrap.tsv")
boot_df.to_csv(out_boot, sep="\t", index=False)
print(f"[OK] {out_boot}")

# structure tests: ANOSIM + dispersion
test_rows = []
for space_name, D, coords in [
    ("BrayOcc", D_bray, coords_bray),
    ("AitchisonCLR", D_ai, coords_ai),
    ("JaccardPA", D_jac, coords_jac),
]:
    R, p = anosim_stat(D, labels, n_perm=N_PERM)
    disp = dispersion_index(coords, labels)
    for _, r in disp.iterrows():
        test_rows.append({
            "space": space_name,
            "metric": "dispersion_mean_dist_to_centroid",
            "group": r["group"],
            "value": r["mean_dist_to_centroid"],
            "p_value": np.nan,
        })
    test_rows.append({
        "space": space_name,
        "metric": "ANOSIM_R",
        "group": "ALL",
        "value": R,
        "p_value": p,
    })

# add variance explained rows for first two axes
for space_name, var in [("BrayOcc", var_bray), ("AitchisonCLR", var_ai), ("JaccardPA", var_jac)]:
    test_rows.append({"space": space_name, "metric": "Axis1_var_pct", "group": "ALL", "value": float(var[0]), "p_value": np.nan})
    test_rows.append({"space": space_name, "metric": "Axis2_var_pct", "group": "ALL", "value": float(var[1] if len(var) > 1 else 0.0), "p_value": np.nan})

tests_df = pd.DataFrame(test_rows)
out_tests = os.path.join(OUT_DIR, "TE_space_group_structure_tests.tsv")
tests_df.to_csv(out_tests, sep="\t", index=False)
print(f"[OK] {out_tests}")

print("\n[INFO] Main-space (BrayOcc) variance:")
print(f"  Axis1={var_bray[0]:.2f}%  Axis2={var_bray[1]:.2f}%")
