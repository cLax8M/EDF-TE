#!/usr/bin/env python3
"""
E03d_article_style_bp_superfamilies.py
=====================================
Article-style TE-space plot using base-pair occupancy per superfamily.

Goal:
  Build a species x superfamily matrix using ALL detected superfamilies.
  Value per cell = occupied bp proportion in the genome for that superfamily.
  If a superfamily is absent in a species, value is 0.

Then:
  - Compute PCoA on Bray-Curtis distance (abundance-sensitive, supports zeros)
  - Draw an article-style scatter (PCo1 vs PCo2)
  - Add dashed gray MST (minimum spanning tree) on clade centroids
  - Color by clade (EDF in yellow-green shades, Dikarya blues, HO red, NU orange)

Outputs:
  te_space_bp_superfamily_matrix.tsv
  merenyi_style_bp_superfamily_pcoa.pdf
  merenyi_style_bp_superfamily_pcoa_34.pdf
    merenyi_style_bp_superfamily_intermediate_axis.pdf
    te_space_intermediate_metrics.tsv

Usage:
  python3 E03d_article_style_bp_superfamilies.py [EARLGREY_DIR] [OUT_DIR]
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import matplotlib.cm as cm

# --- config ------------------------------------------------------------------
EARLGREY_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/clade_outliers"

os.makedirs(OUT_DIR, exist_ok=True)

# --- group/clade style -------------------------------------------------------
BROAD_GROUP = {}
for code in ["HO"]:
    BROAD_GROUP[code] = "Holozoa"
for code in ["NU"]:
    BROAD_GROUP[code] = "Nucleariida"
for code in ["AS", "BA"]:
    BROAD_GROUP[code] = "Dikarya"
for code in ["AP", "BL", "BS", "CH", "CR", "EN", "GL", "KX", "MI", "MO", "MU", "OL", "ZO"]:
    BROAD_GROUP[code] = "EDF"

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
CLADE_COLORS.update({
    "HO": "#e63946",
    "NU": "#f4a261",
    "AS": "#4a90d9",
    "BA": "#1d3557",
})

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

# --- parsing -----------------------------------------------------------------
def extract_class_and_sf(name_str):
    """Return (te_class, superfamily) from EarlGrey token NAME#Class/SF or NAME#Class."""
    if "#" not in name_str:
        return "Unknown", "Unknown"
    after = name_str.split("#", 1)[1]
    if "/" in after:
        te_class, sf = after.split("/", 1)
        return te_class.strip(), sf.strip()
    return after.strip(), after.strip()

# genome sizes from highLevelCount
genome_size = {}
for hl_fp in sorted(glob.glob(os.path.join(EARLGREY_DIR, "*_EarlGrey", "*_summaryFiles", "*.highLevelCount.txt"))):
    sp = os.path.basename(hl_fp).replace(".highLevelCount.txt", "")
    try:
        df = pd.read_csv(hl_fp, sep="\t")
        df.columns = [c.strip() for c in df.columns]
        size_col = [c for c in df.columns if "Genome Size" in c][0]
        gs = pd.to_numeric(df[size_col], errors="coerce").dropna().iloc[0]
        genome_size[sp] = float(gs)
    except Exception:
        continue

# family-level bp by superfamily
records = {}  # sp -> {sf: bp}
fl_files = sorted(glob.glob(os.path.join(EARLGREY_DIR, "*_EarlGrey", "*_summaryFiles", "*.familyLevelCount.txt")))
print(f"[INFO] Found {len(fl_files)} familyLevelCount files")

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
        records[sp] = sf_bp
    except Exception as e:
        print(f"[WARN] {sp}: {e}")

sp_list = sorted(records.keys())
all_sf = sorted({sf for d in records.values() for sf in d})
print(f"[INFO] Parsed {len(sp_list)} species with genome size")
print(f"[INFO] Superfamilies (all, no frequency filtering): {len(all_sf)}")

# matrix in bp/genome proportion (0 for absence)
mat_bp = pd.DataFrame(
    [[records[sp].get(sf, 0.0) for sf in all_sf] for sp in sp_list],
    index=sp_list,
    columns=all_sf,
)

mat_prop_genome = mat_bp.copy()
for sp in sp_list:
    gs = genome_size[sp]
    if gs > 0:
        mat_prop_genome.loc[sp] = mat_prop_genome.loc[sp] / gs
    else:
        mat_prop_genome.loc[sp] = 0.0

# save matrix
out_mat = os.path.join(OUT_DIR, "te_space_bp_superfamily_matrix.tsv")
mat_prop_genome.round(10).to_csv(out_mat, sep="\t")
print(f"[OK] Matrix saved -> {out_mat}")

# Drop species with zero TE content (Bray-Curtis is undefined for all-zero rows)
_zero_rows = mat_prop_genome.index[mat_prop_genome.sum(axis=1) == 0].tolist()
if _zero_rows:
    print(f"[INFO] Dropping {len(_zero_rows)} zero-TE species: {_zero_rows}")
    mat_prop_genome = mat_prop_genome.loc[mat_prop_genome.sum(axis=1) > 0]
    sp_list = mat_prop_genome.index.tolist()

# --- PCoA on Bray-Curtis -----------------------------------------------------
D = squareform(pdist(mat_prop_genome.values, metric="braycurtis"))
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
coords = eigvecs[:, pos] * np.sqrt(eigvals_pos)
var_exp = eigvals_pos / eigvals_pos.sum() * 100

if coords.shape[1] < 4:
    pad = 4 - coords.shape[1]
    coords = np.hstack([coords, np.zeros((coords.shape[0], pad))])
    var_exp = np.concatenate([var_exp, np.zeros(pad)])

# --- plotting helper ----------------------------------------------------------
def draw_panel(ax, use_coords, sp_names, title, xlabel, ylabel):
    # points
    for sp, xy in zip(sp_names, use_coords):
        clade = sp[:2]
        group = BROAD_GROUP.get(clade, "EDF")
        marker = CLADE_MARKERS.get(clade, "o")
        color = CLADE_COLORS.get(clade, GROUP_COLORS[group])
        if hasattr(color, "__len__") and not isinstance(color, str):
            color = tuple(color)
        ax.scatter(xy[0], xy[1], c=[color], marker=marker, s=36,
                   alpha=0.75, linewidths=0.3, edgecolors="white", zorder=3)

    # hulls by broad group
    group_coords = {}
    for sp, xy in zip(sp_names, use_coords):
        group_coords.setdefault(BROAD_GROUP.get(sp[:2], "EDF"), []).append(xy)

    for group, pts in group_coords.items():
        pts = np.array(pts)
        if len(pts) < 3:
            continue
        try:
            hull = ConvexHull(pts)
            ax.fill(pts[hull.vertices, 0], pts[hull.vertices, 1],
                    alpha=0.07, color=GROUP_COLORS[group], zorder=1)
        except Exception:
            pass

    # MST on clade centroids
    clade_points = {}
    for sp, xy in zip(sp_names, use_coords):
        clade_points.setdefault(sp[:2], []).append(xy)
    centroids = {cl: np.array(v).mean(axis=0) for cl, v in clade_points.items()}
    clades = sorted(centroids.keys())
    if len(clades) >= 2:
        G = nx.Graph()
        for i, c1 in enumerate(clades):
            for j, c2 in enumerate(clades):
                if i < j:
                    w = np.linalg.norm(centroids[c1] - centroids[c2])
                    G.add_edge(c1, c2, weight=w)
        mst = nx.minimum_spanning_tree(G)
        for u, v in mst.edges():
            x1, y1 = centroids[u]
            x2, y2 = centroids[v]
            ax.plot([x1, x2], [y1, y2], color="0.55", lw=1.0, ls="--", zorder=2)

    # labels for broad-group centroids
    for group, pts in group_coords.items():
        pts = np.array(pts)
        cx, cy = pts.mean(axis=0)
        ax.text(cx, cy, group, ha="center", va="center", fontsize=9,
                fontweight="bold", color=GROUP_COLORS[group],
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.65))

    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.axhline(0, color="0.85", lw=0.5, zorder=0)
    ax.axvline(0, color="0.85", lw=0.5, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

def get_group_centroids(use_coords, sp_names):
    """Return centroid per broad group for the provided coordinate space."""
    group_coords = {}
    for sp, xy in zip(sp_names, use_coords):
        group = BROAD_GROUP.get(sp[:2], "EDF")
        group_coords.setdefault(group, []).append(xy)
    centroids = {g: np.array(v).mean(axis=0) for g, v in group_coords.items() if len(v) > 0}
    return centroids

def projection_metrics(point, a, b):
    """Project point to line a->b and return interpolation t and perpendicular distance."""
    ab = b - a
    ab2 = float(ab @ ab)
    if ab2 == 0:
        return np.nan, np.nan, a
    t = float(((point - a) @ ab) / ab2)
    foot = a + t * ab
    perp = float(np.linalg.norm(point - foot))
    return t, perp, foot

# main figure PCo1/PCo2
fig, ax = plt.subplots(figsize=(10, 8))
draw_panel(
    ax, coords[:, :2], sp_list,
    "TE superfamily space from occupied bp proportions (all SFs)",
    f"PCo1 ({var_exp[0]:.1f}%)",
    f"PCo2 ({var_exp[1]:.1f}%)",
)

group_handles = [mpatches.Patch(color=GROUP_COLORS[g], label=g, alpha=0.8)
                 for g in ["Holozoa", "Nucleariida", "EDF", "Dikarya"]]
clade_handles = []
for cl, mk in sorted(CLADE_MARKERS.items()):
    cc = CLADE_COLORS.get(cl, "#888888")
    if hasattr(cc, "__len__") and not isinstance(cc, str):
        cc = tuple(cc)
    clade_handles.append(Line2D([0], [0], marker=mk, color="w",
                                markerfacecolor=cc, markeredgecolor="0.4",
                                markersize=6, label=cl))

ax.legend(handles=group_handles + clade_handles, loc="upper right", ncol=2,
          fontsize=7, framealpha=0.85)
plt.tight_layout()
out_12 = os.path.join(OUT_DIR, "merenyi_style_bp_superfamily_pcoa.pdf")
fig.savefig(out_12, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Plot PCo1/PCo2 -> {out_12}")

# secondary figure with PCo3/PCo4
fig, ax = plt.subplots(figsize=(10, 8))
draw_panel(
    ax, coords[:, 2:4], sp_list,
    "TE superfamily space from occupied bp proportions (all SFs)",
    f"PCo3 ({var_exp[2]:.1f}%)",
    f"PCo4 ({var_exp[3]:.1f}%)",
)
ax.legend(handles=group_handles + clade_handles, loc="upper right", ncol=2,
          fontsize=7, framealpha=0.85)
plt.tight_layout()
out_34 = os.path.join(OUT_DIR, "merenyi_style_bp_superfamily_pcoa_34.pdf")
fig.savefig(out_34, dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] Plot PCo3/PCo4 -> {out_34}")

# interpretation figure: Dikarya↔Holozoa axis + EDF/NU projections
coords_12 = coords[:, :2]
cent = get_group_centroids(coords_12, sp_list)
has_needed = all(k in cent for k in ["Dikarya", "Holozoa", "EDF", "Nucleariida"])
if has_needed:
    d_c = cent["Dikarya"]
    h_c = cent["Holozoa"]
    e_c = cent["EDF"]
    n_c = cent["Nucleariida"]

    t_edf, perp_edf, foot_edf = projection_metrics(e_c, d_c, h_c)
    t_nu, perp_nu, foot_nu = projection_metrics(n_c, d_c, h_c)

    # species-level interpolation percentages
    ab = h_c - d_c
    ab2 = float(ab @ ab)
    species_rows = []
    for sp, xy in zip(sp_list, coords_12):
        g = BROAD_GROUP.get(sp[:2], "EDF")
        t = float(((xy - d_c) @ ab) / ab2) if ab2 > 0 else np.nan
        species_rows.append((sp, g, t))
    df_t = pd.DataFrame(species_rows, columns=["species", "group", "t_to_HO_from_Dikarya"])

    inside_edf = ((df_t[df_t["group"] == "EDF"]["t_to_HO_from_Dikarya"] >= 0) &
                  (df_t[df_t["group"] == "EDF"]["t_to_HO_from_Dikarya"] <= 1)).mean() * 100
    inside_nu = ((df_t[df_t["group"] == "Nucleariida"]["t_to_HO_from_Dikarya"] >= 0) &
                 (df_t[df_t["group"] == "Nucleariida"]["t_to_HO_from_Dikarya"] <= 1)).mean() * 100

    # save metrics table
    metrics = pd.DataFrame([
        {"group": "EDF", "t_centroid": t_edf, "perp_centroid": perp_edf,
         "inside_segment_species_pct": inside_edf,
         "distance_to_Dikarya_centroid": float(np.linalg.norm(e_c - d_c)),
         "distance_to_Holozoa_centroid": float(np.linalg.norm(e_c - h_c))},
        {"group": "Nucleariida", "t_centroid": t_nu, "perp_centroid": perp_nu,
         "inside_segment_species_pct": inside_nu,
         "distance_to_Dikarya_centroid": float(np.linalg.norm(n_c - d_c)),
         "distance_to_Holozoa_centroid": float(np.linalg.norm(n_c - h_c))},
    ])
    out_metrics = os.path.join(OUT_DIR, "te_space_intermediate_metrics.tsv")
    metrics.to_csv(out_metrics, sep="\t", index=False)
    print(f"[OK] Intermediate metrics -> {out_metrics}")

    fig, ax = plt.subplots(figsize=(10, 8))
    draw_panel(
        ax, coords_12, sp_list,
        "TE space with Dikarya↔Holozoa intermediate-axis interpretation",
        f"PCo1 ({var_exp[0]:.1f}%)",
        f"PCo2 ({var_exp[1]:.1f}%)",
    )

    # axis between Dikarya and Holozoa centroids
    ax.plot([d_c[0], h_c[0]], [d_c[1], h_c[1]], color="black", lw=2.0, zorder=4)
    ax.text(d_c[0], d_c[1], "Dikarya", fontsize=8, fontweight="bold",
            ha="right", va="bottom", color="black")
    ax.text(h_c[0], h_c[1], "Holozoa", fontsize=8, fontweight="bold",
            ha="left", va="bottom", color="black")

    # EDF and NU projection lines
    ax.plot([e_c[0], foot_edf[0]], [e_c[1], foot_edf[1]], color="#2a9d8f", lw=1.4, ls=":", zorder=5)
    ax.plot([n_c[0], foot_nu[0]], [n_c[1], foot_nu[1]], color="#f4a261", lw=1.4, ls=":", zorder=5)
    ax.scatter([foot_edf[0], foot_nu[0]], [foot_edf[1], foot_nu[1]],
               c=["#2a9d8f", "#f4a261"], s=70, marker="X", zorder=6)

    txt = (
        f"EDF: t={t_edf:.2f}, perp={perp_edf:.3f}, inside={inside_edf:.1f}%\n"
        f"NU: t={t_nu:.2f}, perp={perp_nu:.3f}, inside={inside_nu:.1f}%"
    )
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))

    ax.legend(handles=group_handles + clade_handles, loc="upper right", ncol=2,
              fontsize=7, framealpha=0.85)
    plt.tight_layout()
    out_interp = os.path.join(OUT_DIR, "merenyi_style_bp_superfamily_intermediate_axis.pdf")
    fig.savefig(out_interp, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[OK] Intermediate-axis plot -> {out_interp}")

print("\n[INFO] Variance explained (Bray-Curtis PCoA):")
for i in range(min(4, len(var_exp))):
    print(f"  PCo{i+1}: {var_exp[i]:.2f}%")
