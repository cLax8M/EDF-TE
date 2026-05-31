#!/usr/bin/env python3
"""
03_te_landscape_ordination.py  —  TE-space PCoA and diversity analysis.

Computes three complementary ordinations of the TE superfamily composition
matrix (species × superfamily, values = bp/genome):

  1. Bray-Curtis PCoA on Hellinger-transformed relative abundances
     (sensitive to dominant superfamilies)
  2. Aitchison PCoA on CLR-transformed counts
     (compositional, handles zeros via pseudocount)
  3. Jaccard P/A PCoA on binarised presence/absence
     (detects superfamily repertoire differences)

For each ordination, ANOSIM (999 permutations) tests whether broad taxonomic
groups differ more than expected by chance.

Shannon diversity (H) is computed over superfamily-level bp proportions.

Inputs:
  --superfamily_tsv  Long TSV: species, superfamily, bp_total, pct_genome
  --clade_map        Two-column TSV: species_code<TAB>clade_tag
  --tree             Newick tree (for ANOSIM distance weighting, optional)
  --out_dir          Output directory

Outputs:
  pcoa_braycurtis.pdf / .png
  pcoa_aitchison.pdf  / .png
  pcoa_jaccard.pdf    / .png
  shannon_diversity.tsv
  shannon_boxplot.pdf
  pcoa_variance_explained.tsv
  anosim_results.tsv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon
from sklearn.decomposition import PCA

CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

# Broad groups for ANOSIM
BROAD_GROUP = {
    "AS": "Dikarya", "BA": "Dikarya",
    "HO": "Holozoa",
    "NU": "Nucleariida",
}
# all others → "EDF" (Early-Diverging Fungi)


def broad_group(clade: str) -> str:
    return BROAD_GROUP.get(clade, "EDF")


# =============================================================================
# I/O
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--superfamily_tsv", required=True)
    p.add_argument("--clade_map",       required=True)
    p.add_argument("--tree",            default=None)
    p.add_argument("--out_dir",         default="figures")
    p.add_argument("--min_species",     type=int, default=5,
                   help="Min species per superfamily to retain (default: 5)")
    return p.parse_args()


def build_matrix(tsv_path: str, clade_map: dict,
                 min_species: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(tsv_path, sep="\t")
    df["species"] = df["species"].str.replace(r"\.aa$", "", regex=True)
    df["clade"]   = df["species"].map(clade_map)
    df = df.dropna(subset=["clade"])

    # Drop species with zero TE content (Microsporidia outliers)
    total_bp = df.groupby("species")["bp_total"].sum()
    zero_sp  = total_bp[total_bp == 0].index
    df = df[~df["species"].isin(zero_sp)]
    if len(zero_sp):
        print(f"  Dropped {len(zero_sp)} zero-TE species: {list(zero_sp)}")

    # Pivot: species × superfamily (bp/genome)
    mat = df.pivot_table(index="species", columns="superfamily",
                         values="pct_genome", fill_value=0.0)
    # Filter rare superfamilies
    n_sp_per_sf = (mat > 0).sum(axis=0)
    mat = mat.loc[:, n_sp_per_sf >= min_species]
    print(f"  Matrix: {mat.shape[0]} species × {mat.shape[1]} superfamilies")

    meta = df.drop_duplicates("species").set_index("species")[["clade"]]
    return mat, meta


# =============================================================================
# Ordination helpers
# =============================================================================

def cmdscale(D: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Classical metric multidimensional scaling (PCoA)."""
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    # Keep positive eigenvalues only
    pos = eigvals > 0
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    var_explained = eigvals[pos] / eigvals[pos].sum()
    return coords, var_explained


def hellinger(mat: pd.DataFrame) -> np.ndarray:
    row_sums = mat.sum(axis=1).values[:, None]
    row_sums[row_sums == 0] = 1
    return np.sqrt(mat.values / row_sums)


def clr(mat: pd.DataFrame, pseudo: float = 0.001) -> np.ndarray:
    X = mat.values + pseudo
    log_X = np.log(X)
    return log_X - log_X.mean(axis=1, keepdims=True)


def anosim(dist_mat: np.ndarray, groups: list[str], n_perm: int = 999,
           seed: int = 42) -> dict:
    """
    ANOSIM statistic R = (r_B - r_W) / (n*(n-1)/4)
    where r_B, r_W are mean ranks of between- and within-group distances.
    """
    rng = np.random.default_rng(seed)
    n = dist_mat.shape[0]
    # Upper-triangle indices
    idx = np.triu_indices(n, k=1)
    dists = dist_mat[idx]
    # Rank distances (ascending)
    ranks = pd.Series(dists).rank().values
    same  = np.array([groups[i] == groups[j]
                      for i, j in zip(idx[0], idx[1])], dtype=bool)

    def _R(same_mask, ranks):
        rW = ranks[same_mask].mean()  if same_mask.any()  else 0.0
        rB = ranks[~same_mask].mean() if (~same_mask).any() else 0.0
        return (rB - rW) / (n * (n - 1) / 4)

    R_obs = _R(same, ranks)
    R_perm = [_R(same[rng.permutation(len(same))], ranks)
              for _ in range(n_perm)]
    p_val = (np.sum(np.array(R_perm) >= R_obs) + 1) / (n_perm + 1)
    return {"R": R_obs, "p_value": p_val, "n_perm": n_perm}


# =============================================================================
# Shannon diversity
# =============================================================================

def compute_shannon(mat: pd.DataFrame) -> pd.Series:
    rel = mat.div(mat.sum(axis=1), axis=0)
    # Shannon H = -Σ p_i * log(p_i), treating 0 log 0 = 0
    H = -(rel * np.log(rel + 1e-300)).sum(axis=1)
    H.name = "shannon_H"
    return H


# =============================================================================
# Plotting
# =============================================================================

def _convex_hull_patch(x, y, color, alpha=0.12):
    from scipy.spatial import ConvexHull
    pts = np.column_stack([x, y])
    if len(pts) < 3:
        return None
    hull = ConvexHull(pts)
    poly = Polygon(pts[hull.vertices], closed=True,
                   facecolor=color, edgecolor=color,
                   alpha=alpha, linewidth=1.2)
    return poly


def plot_pcoa(coords: np.ndarray, var_exp: np.ndarray,
              meta: pd.DataFrame, title: str,
              out_stem: str, out_dir: Path) -> None:
    species = meta.index.tolist()
    pc1 = coords[:, 0] if coords.shape[1] > 0 else np.zeros(len(species))
    pc2 = coords[:, 1] if coords.shape[1] > 1 else np.zeros(len(species))

    fig, ax = plt.subplots(figsize=(8, 7))
    clades_present = [c for c in CLADE_ORDER if c in meta["clade"].values]

    for clade in clades_present:
        mask = meta["clade"] == clade
        sx, sy = pc1[mask.values], pc2[mask.values]
        color = CLADE_COLORS.get(clade, "#999999")
        ax.scatter(sx, sy, color=color, s=28, alpha=0.85, zorder=3,
                   edgecolors="white", linewidths=0.3)
        poly = _convex_hull_patch(sx, sy, color)
        if poly:
            ax.add_patch(poly)

    v1 = var_exp[0] * 100 if len(var_exp) > 0 else 0
    v2 = var_exp[1] * 100 if len(var_exp) > 1 else 0
    ax.set_xlabel(f"PCo1 ({v1:.1f}%)", fontsize=11)
    ax.set_ylabel(f"PCo2 ({v2:.1f}%)", fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.axhline(0, linewidth=0.5, color="grey", linestyle="--")
    ax.axvline(0, linewidth=0.5, color="grey", linestyle="--")
    ax.spines[["top", "right"]].set_visible(False)

    patches = [mpatches.Patch(color=CLADE_COLORS.get(c, "#999"), label=c)
               for c in clades_present]
    ax.legend(handles=patches, fontsize=7, frameon=False,
              loc="upper right", ncol=2)

    fig.tight_layout()
    fig.savefig(out_dir / f"{out_stem}.pdf")
    fig.savefig(out_dir / f"{out_stem}.png", dpi=200)
    plt.close(fig)
    print(f"Saved: {out_stem}.pdf  (PCo1={v1:.1f}%, PCo2={v2:.1f}%)")


def plot_shannon_boxplot(shannon: pd.Series, meta: pd.DataFrame,
                         out_dir: Path) -> None:
    clades = [c for c in CLADE_ORDER if c in meta["clade"].values]
    data   = [shannon[meta.index[meta["clade"] == c]].dropna().values
              for c in clades]
    colors = [CLADE_COLORS.get(c, "#999") for c in clades]

    fig, ax = plt.subplots(figsize=(13, 4))
    bp = ax.boxplot(data, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.5))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_xticks(range(1, len(clades) + 1))
    ax.set_xticklabels(clades, fontsize=9)
    ax.set_ylabel("Shannon diversity index (H)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "shannon_boxplot.pdf")
    fig.savefig(out_dir / "shannon_boxplot.png", dpi=200)
    plt.close(fig)
    print("Saved: shannon_boxplot.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))

    print(f"Building superfamily matrix from {args.superfamily_tsv}")
    mat, meta = build_matrix(args.superfamily_tsv, clade_map, args.min_species)

    groups_broad = [broad_group(meta.at[s, "clade"]) for s in meta.index]

    anosim_rows = []

    # ── 1. Bray-Curtis PCoA ──────────────────────────────────────────────────
    print("\nBray-Curtis PCoA ...")
    H_mat = hellinger(mat)
    bc = squareform(pdist(H_mat, "braycurtis"))
    coords_bc, var_bc = cmdscale(bc)
    plot_pcoa(coords_bc, var_bc, meta, "TE-space PCoA (Bray-Curtis)",
              "pcoa_braycurtis", out_dir)
    res = anosim(bc, groups_broad)
    anosim_rows.append({"metric": "Bray-Curtis", **res})

    # ── 2. Aitchison PCoA ────────────────────────────────────────────────────
    print("Aitchison (CLR) PCoA ...")
    clr_mat = clr(mat)
    eu = squareform(pdist(clr_mat, "euclidean"))
    coords_ai, var_ai = cmdscale(eu)
    plot_pcoa(coords_ai, var_ai, meta, "TE-space PCoA (Aitchison CLR)",
              "pcoa_aitchison", out_dir)
    res = anosim(eu, groups_broad)
    anosim_rows.append({"metric": "Aitchison", **res})

    # ── 3. Jaccard P/A PCoA ──────────────────────────────────────────────────
    print("Jaccard presence/absence PCoA ...")
    pa = (mat > 0).astype(float)
    jac = squareform(pdist(pa.values, "jaccard"))
    coords_ja, var_ja = cmdscale(jac)
    plot_pcoa(coords_ja, var_ja, meta, "TE-space PCoA (Jaccard P/A)",
              "pcoa_jaccard", out_dir)
    res = anosim(jac, groups_broad)
    anosim_rows.append({"metric": "Jaccard", **res})

    # ── Shannon diversity ─────────────────────────────────────────────────────
    print("Shannon diversity ...")
    shannon = compute_shannon(mat)
    shannon_df = pd.DataFrame({"species": shannon.index, "shannon_H": shannon.values})
    shannon_df["clade"] = shannon_df["species"].map(clade_map)
    shannon_df.to_csv(out_dir / "shannon_diversity.tsv", sep="\t", index=False)
    plot_shannon_boxplot(shannon, meta, out_dir)

    # ── ANOSIM summary ────────────────────────────────────────────────────────
    anosim_df = pd.DataFrame(anosim_rows)
    anosim_df.to_csv(out_dir / "anosim_results.tsv", sep="\t", index=False)
    print("\nANOSIM results (broad groups: Dikarya, EDF, Holozoa, Nucleariida):")
    print(anosim_df.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()
