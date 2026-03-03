#!/usr/bin/env python3
"""
06_plot_synteny.py — Plot per-clade pairwise synteny heatmaps.

For each clade, reads the normalized synteny ratios (step 05) and produces:
  - A clustered heatmap (seaborn clustermap) where species are clustered by
    mean synteny ratio, mirroring Figure 3C in Sebé-Pedrós et al. 2013.
  - Output: {FIG_DIR}/{CLADE}_synteny_heatmap.pdf  +  .png

The heatmap diagonal (self-vs-self) is set to 1.0 for visualization.

Assembly-fragmentation filter:
  Species whose r_ij_norm == 0 in ALL pairwise comparisons yet share >2000
  orthologs on average are considered assembly-fragmented (no scaffold-level
  contiguity → synteny cannot be computed).  They are excluded by default when
  --exclude_fragmented is set, and a footnote is added to each affected plot.
  Use --exclude_species to remove specific species manually.

Usage:
  python 06_plot_synteny.py --ratio_dir RATIO_DIR --out_dir FIG_DIR
                            [--clade CLADE ...]  [--all]
                            [--exclude_fragmented] [--min_orthologs 2000]
                            [--exclude_species SP1 SP2 ...]
                            [--flag_species SP1 SP2 ...]
                            [--dpi 300]

Author: EDF-TE project
Last updated: 2026-03-03
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

# Full species names for cleaner axis labels (extend as needed)
SPECIES_LABELS = {
    "HO_Acancaste":  "Acanthamoeba castellanii",
    "HO_Callivi":    "Calliphora vicina",
    "HO_Capowcza":   "Capsaspora owczarzaki",
    "HO_Ciointe":    "Ciona intestinalis",
    "HO_Dicdisco":   "Dictyostelium discoideum",
    "HO_Drome":      "Drosophila melanogaster",
    "HO_Homosap":    "Homo sapiens",
    "HO_Hydravulg":  "Hydra vulgaris",
    "HO_Lucicu":     "Lucilia cuprina",
    "HO_Luciser":    "Lucilia sericata",
    "HO_Magisepte":  "Magicicada septendecula",
    "HO_Mneleidy":   "Mnemiopsis leidyi",
    "HO_Monbrevi":   "Monosiga brevicollis",
    "HO_Musdom":     "Musca domestica",
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
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ratio_dir", required=True,
                   help="Directory with per-clade *_synteny_ratios.tsv files (step 05)")
    p.add_argument("--out_dir",   required=True, help="Output directory for figures")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--clade", action="append", choices=CLADES_ALL)
    group.add_argument("--all", action="store_true")
    p.add_argument("--flag_species", nargs="*", default=[],
                   help="Species to flag with * in labels but keep in the heatmap")
    p.add_argument("--exclude_species", nargs="*", default=[],
                   help="Species to remove from every heatmap (explicit list)")
    p.add_argument("--exclude_fragmented", action="store_true",
                   help="Auto-detect and exclude species with r=0 in all pairs "
                        "but >--min_orthologs orthologs (fragmented assemblies)")
    p.add_argument("--min_orthologs", type=int, default=2000,
                   help="Ortholog threshold for fragmentation detection (default: 2000)")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--min_pairs", type=int, default=1,
                   help="Minimum species in clade to attempt heatmap (default: 2)")
    return p.parse_args()


def detect_fragmented(df, min_orthologs=2000):
    """
    Return set of species where r_ij_norm == 0 in ALL pairwise comparisons
    AND mean n_orthologs > min_orthologs.  These are assembly-fragmented genomes
    that have plenty of shared genes but zero syntenic blocks.
    """
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


def load_ratio_matrix(ratio_file, flag_species, exclude_species=()):
    """
    Read *_synteny_ratios.tsv and build a symmetric NxN DataFrame
    with r_ij_norm values. Diagonal = 1.0.
    Species in exclude_species are dropped before building the matrix.
    Returns (matrix_df, label_list, excluded_list).
    """
    df = pd.read_csv(ratio_file, sep="\t")
    if df.empty:
        return None, [], []

    exclude_set = set(exclude_species)

    # Record which requested exclusions actually appear in this clade
    present = set(df["species_A"]) | set(df["species_B"])
    excluded_here = sorted(present & exclude_set)

    # Drop rows involving excluded species
    if exclude_set:
        df = df[~df["species_A"].isin(exclude_set) & ~df["species_B"].isin(exclude_set)]

    if df.empty:
        return None, [], excluded_here

    species = sorted(set(df["species_A"].tolist() + df["species_B"].tolist()))
    n = len(species)
    sp_idx = {sp: i for i, sp in enumerate(species)}

    mat = np.zeros((n, n))
    np.fill_diagonal(mat, 1.0)

    for _, row in df.iterrows():
        i = sp_idx[row["species_A"]]
        j = sp_idx[row["species_B"]]
        val = float(row["r_ij_norm"])
        mat[i, j] = val
        mat[j, i] = val

    # Build display labels
    labels = []
    for sp in species:
        name = SPECIES_LABELS.get(sp, sp.replace("_", " "))
        if sp in flag_species:
            name += " *"
        labels.append(name)

    mat_df = pd.DataFrame(mat, index=labels, columns=labels)
    return mat_df, labels, excluded_here


def plot_heatmap(mat_df, clade_code, out_dir, dpi=300, excluded=None):
    """Generate a clustered heatmap and save as PDF + PNG.
    If excluded is a non-empty list, a footnote listing removed species is added.
    """
    n = len(mat_df)

    footnote = ""
    if excluded:
        short = [sp.split("_", 1)[-1] if "_" in sp else sp for sp in excluded]
        footnote = "Excluded (fragmented assembly): " + ", ".join(short)

    if n < 2:
        print(f"    [{clade_code}] Only {n} species — cannot cluster, writing flat heatmap")
        fig, ax = plt.subplots(figsize=(max(4, n * 0.6 + 1), max(4, n * 0.6 + 1)))
        sns.heatmap(mat_df, ax=ax, vmin=0, vmax=1,
                    cmap="YlOrRd", annot=True, fmt=".2f",
                    linewidths=0.5, square=True)
        ax.set_title(f"Synteny — {clade_code}", fontsize=11)
        if footnote:
            fig.text(0.5, 0.01, footnote, ha="center", va="bottom",
                     fontsize=7, style="italic", color="#555555")
        plt.tight_layout(rect=[0, 0.04 if footnote else 0, 1, 1])
        for ext in ("pdf", "png"):
            fig.savefig(Path(out_dir) / f"{clade_code}_synteny_heatmap.{ext}",
                        dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return

    cell_size = 0.55
    fig_size  = max(5, n * cell_size + 2)

    g = sns.clustermap(
        mat_df,
        method="average",
        metric="euclidean",
        cmap="YlOrRd",
        vmin=0, vmax=1,
        figsize=(fig_size, fig_size),
        linewidths=0.3,
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        annot=(n <= 20),    # annotate values only for small clades
        fmt=".2f" if n <= 20 else "",
    )

    g.ax_heatmap.set_title(f"Gene microsynteny — {clade_code}",
                            fontsize=11, pad=12)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    plt.setp(g.ax_heatmap.get_xticklabels(),
             rotation=45, ha="right", fontsize=max(6, 9 - n // 10))
    plt.setp(g.ax_heatmap.get_yticklabels(),
             rotation=0, fontsize=max(6, 9 - n // 10))

    if footnote:
        g.fig.text(0.5, 0.005, footnote, ha="center", va="bottom",
                   fontsize=7, style="italic", color="#555555")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        g.savefig(out_dir / f"{clade_code}_synteny_heatmap.{ext}",
                  dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"    [{clade_code}] → heatmap saved ({n}×{n})")


def main():
    args = parse_args()
    target_clades = CLADES_ALL if args.all else args.clade
    flag_set = set(args.flag_species)
    exclude_set = set(args.exclude_species)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Auto-detect fragmented assemblies from the combined ratios file if present
    if args.exclude_fragmented:
        combined = Path(args.ratio_dir) / "all_clades_synteny_ratios.tsv"
        if combined.exists():
            df_all = pd.read_csv(combined, sep="\t")
            auto = detect_fragmented(df_all, min_orthologs=args.min_orthologs)
        else:
            # Fall back to per-clade files
            auto = set()
            for clade in target_clades:
                rf = Path(args.ratio_dir) / clade / f"{clade}_synteny_ratios.tsv"
                if rf.exists():
                    auto |= detect_fragmented(pd.read_csv(rf, sep="\t"),
                                              min_orthologs=args.min_orthologs)
        if auto:
            print(f"Auto-excluded {len(auto)} fragmented-assembly species: "
                  + ", ".join(sorted(auto)))
        exclude_set |= auto

    for clade in target_clades:
        ratio_file = Path(args.ratio_dir) / clade / f"{clade}_synteny_ratios.tsv"
        if not ratio_file.exists():
            print(f"  [{clade}] ratio file not found — skipping")
            continue

        mat_df, labels, excluded = load_ratio_matrix(ratio_file, flag_set,
                                                     exclude_species=exclude_set)
        if mat_df is None or len(mat_df) < args.min_pairs:
            print(f"  [{clade}] insufficient data after filtering — skipping")
            continue

        plot_heatmap(mat_df, clade, out_dir, dpi=args.dpi, excluded=excluded)

    print("Step 06 complete.")


if __name__ == "__main__":
    main()
