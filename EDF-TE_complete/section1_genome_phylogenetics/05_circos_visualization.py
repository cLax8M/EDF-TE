#!/usr/bin/env python3
"""
05_circos_visualization.py  —  Circos-style synteny plots for species pairs.

Generates circular chord diagrams showing collinear blocks between two genomes.
Scaffolds are arranged around a circle, coloured by species; chords connect
syntenic regions.

Modes:
  --mode pair    Single pair: requires --sp_a, --sp_b, --clade
  --mode batch   Multiple pairs from a ratio TSV (top-N per clade)
  --mode summary Synteny ratio boxplots and heatmaps (no chord diagrams)

Usage:
  python3 05_circos_visualization.py --mode pair \\
      --sp_a MU_Mucci3 --sp_b MU_Rhimi59 --clade MU \\
      --bed_dir data/bed --collinear_dir data/collinear --out_dir figures/

  python3 05_circos_visualization.py --mode summary \\
      --ratio_dir synteny/ratios --out_dir figures/
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MPath

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", choices=["pair", "batch", "summary"], required=True)
    p.add_argument("--sp_a",          help="Species A code (pair mode)")
    p.add_argument("--sp_b",          help="Species B code (pair mode)")
    p.add_argument("--clade",         help="Clade tag (pair mode)")
    p.add_argument("--bed_dir",       help="Directory with per-species BED files")
    p.add_argument("--collinear_dir", help="Collinear pairs directory")
    p.add_argument("--ratio_dir",     help="Synteny ratio TSV directory (summary/batch)")
    p.add_argument("--out_dir",       default="figures")
    p.add_argument("--max_scaffolds", type=int, default=25,
                   help="Keep top-N scaffolds by size per species (default: 25)")
    p.add_argument("--min_links",     type=int, default=3,
                   help="Drop scaffold pairs with fewer than this many links")
    return p.parse_args()


# =============================================================================
# Circos diagram
# =============================================================================

def _arc(ax, theta1, theta2, r_in, r_out, color, alpha=0.85):
    """Filled arc segment for a chromosome band."""
    n = max(50, int(abs(theta2 - theta1) / np.pi * 200))
    theta = np.linspace(theta1, theta2, n)
    xs = np.concatenate([r_in * np.cos(theta), r_out * np.cos(theta[::-1])])
    ys = np.concatenate([r_in * np.sin(theta), r_out * np.sin(theta[::-1])])
    ax.fill(xs, ys, color=color, alpha=alpha, linewidth=0)


def _bezier_chord(ax, t1_a, t2_a, t1_b, t2_b, r, color, alpha=0.3):
    """Bezier chord connecting two arcs."""
    mid_a = (t1_a + t2_a) / 2
    mid_b = (t1_b + t2_b) / 2
    p0 = np.array([r * np.cos(mid_a), r * np.sin(mid_a)])
    p3 = np.array([r * np.cos(mid_b), r * np.sin(mid_b)])
    p1 = p2 = np.array([0.0, 0.0])  # control point at origin

    verts = [p0, p1, p2, p3]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
    path = MPath(verts, codes)
    patch = mpatches.PathPatch(path, facecolor=color, edgecolor="none", alpha=alpha)
    ax.add_patch(patch)


def _scaffold_angles(sizes: dict, gap_frac: float = 0.02) -> dict:
    """
    Assign angular intervals [theta_start, theta_end] to each scaffold
    proportional to its length.  Returns {scaffold: (t_start, t_end)}.
    """
    total = sum(sizes.values())
    gap = 2 * np.pi * gap_frac
    n_scaf = len(sizes)
    usable = 2 * np.pi - gap * n_scaf
    angles = {}
    theta = 0.0
    for scaf, length in sorted(sizes.items(), key=lambda x: -x[1]):
        span = usable * length / total
        angles[scaf] = (theta, theta + span)
        theta += span + gap
    return angles


def plot_circos_pair(sp_a: str, sp_b: str, bed_dir: Path,
                     collinear_dir: Path, out_dir: Path,
                     max_scaffolds: int, min_links: int,
                     clade: str) -> None:
    """Draw a circos chord diagram for one species pair."""
    clade_a, clade_b = sp_a.split("_")[0], sp_b.split("_")[0]
    color_a = CLADE_COLORS.get(clade_a, "#888888")
    color_b = CLADE_COLORS.get(clade_b, "#555555")

    # Load BED files (chrom, start, end, gene)
    def load_bed(sp):
        bed_file = bed_dir / f"{sp}.bed"
        df = pd.read_csv(bed_file, sep="\t", header=None,
                         names=["chrom", "start", "end", "gene"])
        sizes = df.groupby("chrom")["end"].max().sort_values(ascending=False)
        sizes = sizes.head(max_scaffolds)
        df = df[df["chrom"].isin(sizes.index)]
        return df, sizes.to_dict()

    df_a, sizes_a = load_bed(sp_a)
    df_b, sizes_b = load_bed(sp_b)

    # Load collinear pairs
    pair_key = f"{sp_a}__{sp_b}"
    pair_file = collinear_dir / f"{pair_key}.tsv"
    if not pair_file.exists():
        pair_key = f"{sp_b}__{sp_a}"
        pair_file = collinear_dir / f"{pair_key}.tsv"
    if not pair_file.exists():
        print(f"  No collinear file for {sp_a}×{sp_b}, skipping.")
        return

    links = pd.read_csv(pair_file, sep="\t",
                        names=["gene_a", "gene_b", "chrom_a", "chrom_b",
                               "pos_a", "pos_b"])
    links = links[links["chrom_a"].isin(sizes_a) & links["chrom_b"].isin(sizes_b)]
    # Filter low-link scaffold pairs
    pair_counts = links.groupby(["chrom_a", "chrom_b"]).size()
    valid_pairs = pair_counts[pair_counts >= min_links].index
    links = links[links.set_index(["chrom_a", "chrom_b"]).index.isin(valid_pairs)]

    if len(links) == 0:
        print(f"  No links pass filters for {sp_a}×{sp_b}, skipping.")
        return

    # Assign angles: species A on top half, B on bottom
    angles_a = _scaffold_angles(sizes_a)
    angles_b = {}
    for scaf, (t1, t2) in _scaffold_angles(sizes_b).items():
        angles_b[scaf] = (t1 + np.pi, t2 + np.pi)

    # Gene position → angle
    def gene_angle(df, sizes, angles):
        gene_angles = {}
        for chrom, grp in df.groupby("chrom"):
            if chrom not in angles:
                continue
            t1, t2 = angles[chrom]
            chrom_len = sizes[chrom]
            for _, row in grp.iterrows():
                pos = (row["start"] + row["end"]) / 2
                gene_angles[row["gene"]] = t1 + (t2 - t1) * pos / chrom_len
        return gene_angles

    gangle_a = gene_angle(df_a, sizes_a, angles_a)
    gangle_b = gene_angle(df_b, sizes_b, angles_b)

    # Draw
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect("equal")
    ax.axis("off")

    R_OUT, R_IN = 1.0, 0.88

    for chrom, (t1, t2) in angles_a.items():
        _arc(ax, t1, t2, R_IN, R_OUT, color_a)
    for chrom, (t1, t2) in angles_b.items():
        _arc(ax, t1, t2, R_IN, R_OUT, color_b)

    for _, row in links.iterrows():
        ta = gangle_a.get(row["gene_a"])
        tb = gangle_b.get(row["gene_b"])
        if ta is None or tb is None:
            continue
        _bezier_chord(ax, ta - 0.005, ta + 0.005,
                      tb - 0.005, tb + 0.005, R_IN, "#888888", alpha=0.15)

    # Legend
    patches = [mpatches.Patch(color=color_a, label=sp_a.replace("_", " ")),
               mpatches.Patch(color=color_b, label=sp_b.replace("_", " "))]
    ax.legend(handles=patches, loc="lower center", frameon=False, fontsize=9,
              bbox_to_anchor=(0.5, -0.05))

    n_links = len(links)
    ax.set_title(f"{sp_a} – {sp_b}  ({n_links} collinear pairs)", fontsize=11)

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"circos_{sp_a}_{sp_b}"
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {stem}.pdf")


# =============================================================================
# Summary mode: synteny ratio plots
# =============================================================================

def plot_summary(ratio_dir: Path, out_dir: Path) -> None:
    """Boxplots and heatmap of synteny ratios per clade."""
    all_file = ratio_dir / "all_clades_synteny_ratios.tsv"
    if not all_file.exists():
        print(f"Missing: {all_file}")
        return

    df = pd.read_csv(all_file, sep="\t")
    clades = [c for c in CLADE_ORDER if c in df["clade"].values]

    # Boxplot
    fig, ax = plt.subplots(figsize=(12, 4))
    data_by_clade = [df.loc[df["clade"] == c, "r_ij_norm"].dropna().values
                     for c in clades]
    colors = [CLADE_COLORS.get(c, "#888") for c in clades]
    bp = ax.boxplot(data_by_clade, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.5))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_xticks(range(1, len(clades) + 1))
    ax.set_xticklabels(clades, fontsize=9)
    ax.set_ylabel("Normalized synteny ratio", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "synteny_ratio_boxplot.pdf")
    fig.savefig(out_dir / "synteny_ratio_boxplot.png", dpi=200)
    plt.close(fig)
    print("Saved: synteny_ratio_boxplot.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.mode == "pair":
        if not (args.sp_a and args.sp_b and args.bed_dir and args.collinear_dir):
            print("pair mode requires --sp_a, --sp_b, --bed_dir, --collinear_dir")
            sys.exit(1)
        plot_circos_pair(args.sp_a, args.sp_b,
                         Path(args.bed_dir), Path(args.collinear_dir),
                         out_dir, args.max_scaffolds, args.min_links,
                         args.clade or "")

    elif args.mode == "summary":
        if not args.ratio_dir:
            print("summary mode requires --ratio_dir")
            sys.exit(1)
        plot_summary(Path(args.ratio_dir), out_dir)

    elif args.mode == "batch":
        if not (args.ratio_dir and args.bed_dir and args.collinear_dir):
            print("batch mode requires --ratio_dir, --bed_dir, --collinear_dir")
            sys.exit(1)
        all_file = Path(args.ratio_dir) / "all_clades_synteny_ratios.tsv"
        df = pd.read_csv(all_file, sep="\t")
        for _, row in df.iterrows():
            plot_circos_pair(row["species_A"], row["species_B"],
                             Path(args.bed_dir), Path(args.collinear_dir),
                             out_dir, args.max_scaffolds, args.min_links,
                             row.get("clade", ""))


if __name__ == "__main__":
    main()
