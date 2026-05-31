#!/usr/bin/env python3
"""
G01_te_landscape_pcoa.py
PCoA of TE-superfamily composition across all annotated species
(Merényi et al. 2023, NEE — Fig 1b style).

Input:
    C01_kimura_superfamily.tsv  (long: species, superfamily, total_bp, pct_of_TE, ...)
    figure_species_data.tsv     (clade per species)

Steps:
    1. Pivot pct_of_TE → species × superfamily matrix.
    2. Hellinger-transform (sqrt of relative abundance) — stabilises Bray-Curtis.
    3. Bray-Curtis distance matrix.
    4. Classical PCoA (via scipy / cmdscale).
    5. Plot dim1 vs dim2, convex hulls per clade.

Output:
    results/Section2_TElandscape/figures/G01_te_landscape_pcoa.{pdf,png}
    results/Section2_TElandscape/G01_pcoa_coords.tsv
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull

BASE = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape")
SF_FILE = BASE / "C01_kimura_superfamily.tsv"
META_FILE = BASE / "figure_species_data.tsv"
OUT_DIR = BASE / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Clade colour palette (matches Section2 cladogram) ───────────────────────
CLADE_COLORS = {
    "Holozoa": "#1f77b4",
    "Cryptomycota": "#aec7e8",
    "Microsporidia": "#9467bd",
    "Aphelida": "#c5b0d5",
    "Olpidiomycota": "#8c564b",
    "Chytridiomycota": "#2ca02c",
    "Neocallimastigomycota": "#98df8a",
    "Blastocladiomycota": "#17becf",
    "Basidiobolales": "#7f7f7f",
    "Entomophthoromycota": "#bcbd22",
    "Zoopagomycota": "#e377c2",
    "Kickxellomycotina": "#ff9896",
    "Mortierellomycota": "#ff7f0e",
    "Mucoromycota": "#d62728",
    "Glomeromycota": "#8c6d31",
    "Basidiomycota": "#e7969c",
    "Ascomycota": "#393b79",
}

def main():
    # Load
    sf = pd.read_csv(SF_FILE, sep="\t")
    meta = pd.read_csv(META_FILE, sep="\t")[["taxon", "clade"]].rename(columns={"taxon": "species"})

    # Pivot → wide matrix species × superfamily, fill 0
    M = sf.pivot_table(index="species", columns="superfamily",
                       values="pct_of_TE", aggfunc="sum", fill_value=0.0)
    print(f"[INFO] Matrix: {M.shape[0]} species × {M.shape[1]} superfamilies", file=sys.stderr)

    # Drop superfamilies present in <2 species (noise)
    keep = (M > 0).sum(axis=0) >= 2
    M = M.loc[:, keep]
    print(f"[INFO] After filter (≥2 species): {M.shape[1]} superfamilies", file=sys.stderr)

    # Re-normalise rows to sum to 1 (relative abundance)
    M = M.div(M.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    # Hellinger transform (sqrt of relative abundance)
    Mh = np.sqrt(M.values)

    # Bray-Curtis distance (on raw relative abundance, not Hellinger)
    D = squareform(pdist(M.values, metric="braycurtis"))

    # Classical PCoA (Gower 1966):
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    # Coords = eigvecs * sqrt(positive eigvals)
    pos = eigvals > 0
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    var_exp = eigvals[pos] / eigvals[pos].sum()
    print(f"[INFO] PCo1 = {var_exp[0]*100:.1f}%   PCo2 = {var_exp[1]*100:.1f}%", file=sys.stderr)

    # DataFrame with coords + clade
    coord_df = pd.DataFrame(coords[:, :3], index=M.index, columns=["PCo1", "PCo2", "PCo3"])
    coord_df = coord_df.merge(meta, left_index=True, right_on="species", how="left")
    coord_df["clade"] = coord_df["clade"].fillna("Unknown")
    out_tsv = BASE / "G01_pcoa_coords.tsv"
    coord_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[INFO] coords → {out_tsv}", file=sys.stderr)

    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11, 9))

    clade_order = sorted(coord_df["clade"].dropna().unique(),
                         key=lambda c: -len(coord_df[coord_df["clade"] == c]))

    for clade in clade_order:
        sub = coord_df[coord_df["clade"] == clade]
        c = CLADE_COLORS.get(clade, "#999999")
        pts = sub[["PCo1", "PCo2"]].values
        # Convex hull (≥3 pts)
        if len(pts) >= 3:
            try:
                hull = ConvexHull(pts)
                hpts = pts[hull.vertices]
                poly = Polygon(hpts, closed=True, facecolor=c, edgecolor=c,
                               alpha=0.18, linewidth=1.2, zorder=1)
                ax.add_patch(poly)
                # Centroid for label
                cx, cy = pts.mean(axis=0)
                ax.text(cx, cy, clade, fontsize=9, fontweight="bold",
                        ha="center", va="center", color=c,
                        bbox=dict(facecolor="white", edgecolor="none",
                                  alpha=0.6, pad=1.2),
                        zorder=4)
            except Exception:
                pass
        # Scatter points
        ax.scatter(pts[:, 0], pts[:, 1], s=22, c=c, edgecolors="white",
                   linewidth=0.5, zorder=3, label=f"{clade} (n={len(sub)})")

    ax.axhline(0, color="#cccccc", lw=0.5, zorder=0)
    ax.axvline(0, color="#cccccc", lw=0.5, zorder=0)
    ax.set_xlabel(f"PCo1 ({var_exp[0]*100:.1f}%)", fontsize=12)
    ax.set_ylabel(f"PCo2 ({var_exp[1]*100:.1f}%)", fontsize=12)
    ax.set_title("TE landscape PCoA (Bray–Curtis on TE superfamily composition)",
                 fontsize=13, pad=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend outside
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=8, frameon=False,
              title=f"Clades (n={len(coord_df)} species)", title_fontsize=9)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"G01_te_landscape_pcoa.{ext}"
        plt.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[OK] {out}", file=sys.stderr)
    plt.close()

if __name__ == "__main__":
    main()
