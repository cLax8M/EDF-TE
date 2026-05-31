#!/usr/bin/env python3
"""
10_plot_corrected_synteny.py — Heatmaps from patristic-distance-corrected synteny.

Reads the output of 09_patristic_correction.py and generates per-clade
clustered heatmaps using r_corrected_norm (analogous to 06_plot_synteny.py
but using corrected values). Also produces a side-by-side comparison figure
for each clade showing original vs corrected heatmaps.

Usage:
  python 10_plot_corrected_synteny.py \\
      --corrected /path/to/corrected_synteny_ratios.tsv \\
      --out_dir   /path/to/output/figures \\
      --all
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

CLADES_ALL = ["AP", "AS", "BA", "BL", "BS", "CH", "CR", "EN",
              "GL", "HO", "KX", "MI", "MO", "MU", "NU", "OL", "ZO"]

SPECIES_LABELS = {
    "HO_Acancaste":  "Acanthamoeba castellanii",
    "HO_Capowcza":   "Capsaspora owczarzaki",
    "HO_Ciointe":    "Ciona intestinalis",
    "HO_Dicdisco":   "Dictyostelium discoideum",
    "HO_Drome":      "Drosophila melanogaster",
    "HO_Homosap":    "Homo sapiens",
    "HO_Hydravulg":  "Hydra vulgaris",
    "HO_Magisepte":  "Magicicada septendecula",
    "HO_Mneleidy":   "Mnemiopsis leidyi",
    "HO_Monbrevi":   "Monosiga brevicollis",
    "HO_Naegrub":    "Naegleria gruberi",
    "HO_Pluxyl":     "Plutella xylostella",
    "HO_Salpirose":  "Salpingoeca rosetta",
    "HO_Sphaearti":  "Sphaeroforma arctica",
    "HO_Thecatra":   "Thecamonas trahens",
    "HO_Trichadha":  "Trichoplax adhaerens",
    "NU_Fontalba":   "Fonticula alba",
    "NU_Parvatl":    "Parvularia atlantis",
    "AP_Amoocci":    "Amoeboaphelidium occidentale",
    "AP_Amoprot":    "Amoeboaphelidium protococcarum",
    "AS_Sacce1":     "Saccharomyces cerevisiae",
    "AS_Schpo1":     "Schizosaccharomyces pombe",
    "AS_Ncrass2":    "Neurospora crassa",
    "AS_Bcine1":     "Botrytis cinerea",
    "MU_Rhiarbi":    "Rhizopus arrhizus",
    "MU_Mucci2":     "Mucor circinelloides",
    "MO_Morvert2":   "Mortierella verticillata",
    "CH_Batra1":     "Batrachochytrium dendrobatidis",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--corrected", required=True,
                   help="corrected_synteny_ratios.tsv from step 09")
    p.add_argument("--out_dir", required=True, help="Output directory for figures")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--clade", action="append", choices=CLADES_ALL)
    group.add_argument("--all", action="store_true")
    p.add_argument("--exclude_fragmented", action="store_true",
                   help="Exclude species with r=0 in all pairs but many orthologs")
    p.add_argument("--min_orthologs", type=int, default=2000)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def detect_fragmented(df, min_orthologs=2000):
    fragmented = set()
    all_sp = set(df["species_A"]) | set(df["species_B"])
    for sp in all_sp:
        mask = (df["species_A"] == sp) | (df["species_B"] == sp)
        sub = df[mask]
        if sub.empty:
            continue
        if sub["r_ij_norm"].max() == 0 and sub["n_orthologs"].mean() > min_orthologs:
            fragmented.add(sp)
    return fragmented


def build_matrix(clade_df, value_col="r_corrected_norm", exclude=set()):
    """Build symmetric N×N matrix from pairwise data."""
    df = clade_df.copy()
    if exclude:
        df = df[~df["species_A"].isin(exclude) & ~df["species_B"].isin(exclude)]
    if df.empty:
        return None

    species = sorted(set(df["species_A"]) | set(df["species_B"]))
    n = len(species)
    sp_idx = {sp: i for i, sp in enumerate(species)}

    mat = np.zeros((n, n))
    np.fill_diagonal(mat, 1.0)
    for _, row in df.iterrows():
        i = sp_idx[row["species_A"]]
        j = sp_idx[row["species_B"]]
        val = row[value_col]
        if pd.isna(val):
            val = 0.0
        mat[i, j] = val
        mat[j, i] = val

    labels = []
    for sp in species:
        labels.append(SPECIES_LABELS.get(sp, sp.replace("_", " ")))
    return pd.DataFrame(mat, index=labels, columns=labels)


def plot_single_heatmap(mat_df, title, out_path, dpi=300, cmap="YlOrRd"):
    n = len(mat_df)
    if n < 2:
        fig, ax = plt.subplots(figsize=(4, 4))
        sns.heatmap(mat_df, ax=ax, vmin=0, vmax=1, cmap=cmap,
                    annot=True, fmt=".2f", linewidths=0.5, square=True)
        ax.set_title(title, fontsize=10)
        plt.tight_layout()
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return

    cell_size = 0.55
    fig_size = max(5, n * cell_size + 2)
    g = sns.clustermap(
        mat_df, method="average", metric="euclidean",
        cmap=cmap, vmin=0, vmax=1,
        figsize=(fig_size, fig_size),
        linewidths=0.3, xticklabels=True, yticklabels=True,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        annot=(n <= 20), fmt=".2f" if n <= 20 else "",
    )
    g.ax_heatmap.set_title(title, fontsize=11, pad=12)
    plt.setp(g.ax_heatmap.get_xticklabels(),
             rotation=45, ha="right", fontsize=max(6, 9 - n // 10))
    plt.setp(g.ax_heatmap.get_yticklabels(),
             rotation=0, fontsize=max(6, 9 - n // 10))
    g.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()


def plot_comparison(mat_orig, mat_corr, clade, out_dir, dpi=300):
    """Side-by-side: original vs corrected heatmap."""
    n = len(mat_orig)
    if n < 2:
        return

    cell_size = 0.45
    fig_w = max(10, n * cell_size * 2 + 4)
    fig_h = max(5, n * cell_size + 2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_w, fig_h))

    sns.heatmap(mat_orig, ax=ax1, vmin=0, vmax=1, cmap="YlOrRd",
                linewidths=0.2, square=True, xticklabels=True, yticklabels=True,
                annot=(n <= 15), fmt=".2f" if n <= 15 else "",
                cbar_kws={"shrink": 0.5})
    ax1.set_title(f"{clade} — Original r_ij_norm", fontsize=10)
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", fontsize=max(5, 8 - n // 10))
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=max(5, 8 - n // 10))

    sns.heatmap(mat_corr, ax=ax2, vmin=0, vmax=1, cmap="YlOrRd",
                linewidths=0.2, square=True, xticklabels=True, yticklabels=True,
                annot=(n <= 15), fmt=".2f" if n <= 15 else "",
                cbar_kws={"shrink": 0.5})
    ax2.set_title(f"{clade} — Corrected r (guide tree)", fontsize=10)
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", fontsize=max(5, 8 - n // 10))
    plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=max(5, 8 - n // 10))

    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(Path(out_dir) / f"{clade}_comparison_orig_vs_corrected.{ext}",
                    dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print(f"Loading corrected ratios: {args.corrected}")
    df = pd.read_csv(args.corrected, sep="\t")
    print(f"  {len(df)} pairs across {df['clade'].nunique()} clades")

    target = CLADES_ALL if args.all else args.clade
    exclude = set()
    if args.exclude_fragmented:
        exclude = detect_fragmented(df, args.min_orthologs)
        if exclude:
            print(f"Excluding {len(exclude)} fragmented: {', '.join(sorted(exclude))}")

    for clade in sorted(df["clade"].unique()):
        if clade not in target:
            continue
        clade_df = df[df["clade"] == clade]
        n_sp = len(set(clade_df["species_A"]) | set(clade_df["species_B"]))

        if n_sp < 2:
            print(f"  [{clade}] only {n_sp} species — skipping")
            continue

        # Build both matrices
        mat_orig = build_matrix(clade_df, "r_ij_norm", exclude)
        mat_corr = build_matrix(clade_df, "r_corrected_norm", exclude)

        if mat_corr is None or len(mat_corr) < 2:
            print(f"  [{clade}] no data after filtering — skipping")
            continue

        # Individual corrected heatmap
        for ext in ("pdf", "png"):
            plot_single_heatmap(
                mat_corr,
                f"Gene microsynteny (corrected) — {clade}",
                out / f"{clade}_synteny_corrected.{ext}",
                dpi=args.dpi
            )

        # Comparison figure
        if mat_orig is not None:
            plot_comparison(mat_orig, mat_corr, clade, out, dpi=args.dpi)

        print(f"  [{clade}] {len(mat_corr)}×{len(mat_corr)} heatmap saved")

    print(f"\nFigures saved to {out}")


if __name__ == "__main__":
    main()
