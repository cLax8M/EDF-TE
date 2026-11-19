#!/usr/bin/env python3

# Pipeline to test whether MAT loci are enriched in transposable elements (TEs) compared to random genomic regions.

# imports 
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow

# clade order used for consistent plotting across fungal phylogeny
CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

# clade color map for figures
CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

# TE class colors (for snapshot plots)
TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#984ea3",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}


# ARGUMENT PARSING
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mat_coords",  required=True)
    p.add_argument("--te_bed_dir",  required=True)
    p.add_argument("--genome_dir",  required=True)
    p.add_argument("--clade_map",   required=True)
    p.add_argument("--out_dir",     default="sex_locus_out")
    p.add_argument("--window_kb",   type=int, default=10)
    p.add_argument("--n_random",    type=int, default=100)
    p.add_argument("--snapshots",   nargs="*", default=[],
                   help="Species codes for individual MAT snapshots.")
    p.add_argument("--seed",        type=int, default=42)
    return p.parse_args()



# TE content in windows
# compute chromosome lengths from FASTA file (needed for random windows)
def genome_scaffold_lengths(fasta_path: Path) -> dict[str, int]:
    lengths = {}
    cur = None
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                cur = line[1:].split()[0]
                lengths[cur] = 0
            elif cur:
                lengths[cur] += len(line)
    return lengths

# TE COVERAGE IN GENOMIC WINDOWS
# compute how many TE base pairs overlap a given genomic window
# (clips TE intervals to window boundaries)
def te_bp_in_window(te_bed: pd.DataFrame, chrom: str, start: int, end: int) -> int:
    sub = te_bed[(te_bed["chrom"] == chrom) &
                 (te_bed["te_end"] > start) &
                 (te_bed["te_start"] < end)].copy()
    if sub.empty:
        return 0
        
    # compute exact overlap
    sub["ov_start"] = sub["te_start"].clip(lower=start)
    sub["ov_end"]   = sub["te_end"].clip(upper=end)
    return (sub["ov_end"] - sub["ov_start"]).sum()

# generate random genomic windows of fixed size within a chromosome
# used as null model for TE enrichment
def random_windows(chrom_len: int, chrom: str, window: int,
                   n: int, rng: np.random.Generator) -> list[tuple[str, int, int]]:
    out = []
    for _ in range(n * 5):
        s = rng.integers(0, max(1, chrom_len - window))
        out.append((chrom, int(s), int(s + window)))
        if len(out) >= n:
            break
    return out[:n]

# main computation:
# for each species compute TE content in MAT locus vs random windows
def compute_mat_te_content(mat_coords: pd.DataFrame,
                            te_bed_dir: Path, genome_dir: Path,
                            window_bp: int, n_random: int,
                            rng: np.random.Generator,
                            clade_map: dict) -> pd.DataFrame:
    rows = []
    for _, row in mat_coords.iterrows():
        sp    = row["species"]
        chrom = row["chromosome"]
        # MAT window expanded
        mat_s = int(row["mat_start"]) - window_bp
        mat_e = int(row["mat_end"])   + window_bp

        # load TE BED file for species
        te_files = list(Path(te_bed_dir).glob(f"{sp}*.filteredRepeats.bed")) + \
                   list(Path(te_bed_dir).glob(f"{sp}*.bed"))
        if not te_files:
            continue

        te_bed = pd.read_csv(te_files[0], sep="\t", header=None, usecols=[0, 1, 2, 3],
                              names=["chrom", "te_start", "te_end", "family"])
        te_bed = te_bed[te_bed["chrom"] == chrom]

        # estimate chromosome length (fallback if FASTA not available)
        genome_file = list(Path(genome_dir).glob(f"{sp}*.fasta")) + \
                      list(Path(genome_dir).glob(f"{sp}*.fa"))
        if genome_file:
            chrom_lengths = genome_scaffold_lengths(genome_file[0])
            chrom_len = chrom_lengths.get(chrom, mat_e + window_bp * 10)
        else:
            chrom_len = mat_e + window_bp * 10

        # MAT window
        win_len = mat_e - mat_s
        mat_te  = te_bp_in_window(te_bed, chrom, mat_s, mat_e)
        mat_pct = 100.0 * mat_te / win_len if win_len > 0 else 0.0

        # Random comparison windows
        rand_windows = random_windows(chrom_len, chrom, win_len, n_random, rng)
        rand_pcts = [
            100.0 * te_bp_in_window(te_bed, w_chrom, w_s, w_e) / win_len
            for w_chrom, w_s, w_e in rand_windows
        ]

        rows.append({
            "species":      sp,
            "clade":        clade_map.get(sp, sp[:2]),
            "chromosome":   chrom,
            "mat_te_pct":   mat_pct,
            "rand_te_mean": np.mean(rand_pcts),
            "rand_te_sd":   np.std(rand_pcts),
            "n_random":     len(rand_pcts),
        })
    return pd.DataFrame(rows)



# Figures

# plot main comparison:
# MAT TE content vs random genomic expectation
def plot_mat_te_comparison(df: pd.DataFrame, out_dir: Path) -> None:
    """Scatter: MAT TE% vs random mean, grouped by clade."""
    clades = [c for c in CLADE_ORDER if c in df["clade"].values]

    fig, ax = plt.subplots(figsize=(8, 6))
    for clade in clades:
        sub = df[df["clade"] == clade]
        color = CLADE_COLORS.get(clade, "#999")
        ax.scatter(sub["rand_te_mean"], sub["mat_te_pct"],
                   color=color, s=40, alpha=0.8, zorder=3,
                   edgecolors="white", linewidths=0.4, label=clade)
        for _, row in sub.iterrows():
            ax.plot([row["rand_te_mean"], row["rand_te_mean"]],
                    [row["mat_te_pct"], row["mat_te_pct"]],
                    color=color, alpha=0.4, linewidth=0.5)

    lim = max(ax.get_xlim()[1], ax.get_ylim()[1]) + 5
    ax.plot([0, lim], [0, lim], "--", color="grey", linewidth=0.8, zorder=1)
    ax.set_xlabel("Random window TE content (%)", fontsize=10)
    ax.set_ylabel("MAT locus TE content (%)", fontsize=10)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)

    # Wilcoxon test across all species
    if len(df) >= 5:
        stat, pval = wilcoxon(df["mat_te_pct"], df["rand_te_mean"])
        ax.text(0.05, 0.95, f"Wilcoxon p={pval:.3g}", transform=ax.transAxes,
                fontsize=9, va="top")

    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "sex_locus_te_content.pdf")
    fig.savefig(out_dir / "sex_locus_te_content.png", dpi=200)
    plt.close(fig)
    print("Saved: sex_locus_te_content.pdf")


# plot genomic snapshot of MAT locus:
# TE distribution + MAT region context
def plot_mat_snapshot(sp: str, mat_coords_row: pd.Series,
                      te_bed_dir: Path, gene_bed_dir: Path,
                      window_bp: int, out_dir: Path) -> None:
    """Arrow diagram of genes and TEs around a MAT locus."""
    chrom = mat_coords_row["chromosome"]
    mat_s = int(mat_coords_row["mat_start"])
    mat_e = int(mat_coords_row["mat_end"])
    view_s = mat_s - window_bp
    view_e = mat_e + window_bp

    te_files = list(Path(te_bed_dir).glob(f"{sp}*.filteredRepeats.bed"))
    if not te_files:
        return
    te_bed = pd.read_csv(te_files[0], sep="\t", header=None, usecols=[0, 1, 2, 3, 5],
                          names=["chrom", "te_start", "te_end", "family", "strand"])
    te_bed = te_bed[(te_bed["chrom"] == chrom) &
                    (te_bed["te_end"] > view_s) & (te_bed["te_start"] < view_e)]

    fig, ax = plt.subplots(figsize=(14, 3))
    ax.set_xlim(view_s, view_e)
    ax.set_ylim(-1.2, 1.2)
    ax.axis("off")

    # Draw chromosome line
    ax.plot([view_s, view_e], [0, 0], color="#333333", linewidth=2, zorder=1)

    # Highlight MAT region
    ax.axvspan(mat_s, mat_e, alpha=0.12, color="#4dac26", zorder=0)

    # Draw TEs
    def te_class(f):
        f = f.upper()
        if "LTR" in f:    return "LTR"
        if "DNA" in f:    return "DNA"
        if "LINE" in f:   return "LINE"
        if "SINE" in f:   return "SINE"
        return "Unclassified"

    for _, te in te_bed.iterrows():
        color = TE_CLASS_COLORS.get(te_class(te["family"]), "#999")
        s, e = max(te["te_start"], view_s), min(te["te_end"], view_e)
        y = 0.45 if te.get("strand", "+") == "+" else -0.45
        ax.barh(y, e - s, left=s, height=0.3, color=color, alpha=0.8, zorder=2)

    # Scale bar
    sb = (view_e - view_s) / 10
    ax.plot([view_e - sb, view_e], [-1.0, -1.0], color="black", linewidth=2)
    ax.text(view_e - sb / 2, -1.15, f"{sb / 1000:.0f} kb",
            ha="center", fontsize=8)

    ax.set_title(f"{sp.replace('_', ' ')} — MAT locus (±{window_bp // 1000} kb)",
                 fontsize=11)

    # Legend
    patches = [mpatches.Patch(color=c, label=k)
               for k, c in TE_CLASS_COLORS.items()]
    ax.legend(handles=patches, loc="upper right", fontsize=7,
              frameon=False, ncol=3)

    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"mat_snapshot_{sp}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"mat_snapshot_{sp}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: mat_snapshot_{sp}.pdf")



# main


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))

    mat_coords = pd.read_csv(args.mat_coords, sep="\t")
    rng = np.random.default_rng(args.seed)
    window_bp = args.window_kb * 1000

    print(f"Computing TE content at MAT loci ({args.window_kb} kb windows) ...")
    df = compute_mat_te_content(mat_coords, Path(args.te_bed_dir),
                                Path(args.genome_dir), window_bp,
                                args.n_random, rng, clade_map)
    df.to_csv(out_dir / "mat_te_content.tsv", sep="\t", index=False)
    print(f"  {len(df)} species processed")

    plot_mat_te_comparison(df, fig_dir)

    for sp in args.snapshots:
        row = mat_coords[mat_coords["species"] == sp]
        if row.empty:
            print(f"  No MAT coords for {sp}, skipping snapshot.")
            continue
        plot_mat_snapshot(sp, row.iloc[0], Path(args.te_bed_dir),
                          Path(args.gene_bed_dir) if hasattr(args, "gene_bed_dir") else Path("."),
                          window_bp, fig_dir)

    print("Done.")


if __name__ == "__main__":
    main()
