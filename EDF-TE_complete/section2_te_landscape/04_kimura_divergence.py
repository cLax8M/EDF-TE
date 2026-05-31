#!/usr/bin/env python3
"""
04_kimura_divergence.py  —  Kimura divergence analysis and Recency Index.

Processes RepeatMasker .align files to compute per-species, per-superfamily
Kimura substitution divergence distributions and a summary Recency Index (RI).

Recency Index:
  RI = bp with Kimura divergence < 0.1 / total TE bp
  Values near 1 indicate recent/active TE families; near 0 indicate old/eroded.

Analyses:
  A. Per-species divergence landscape (histogram stacked by TE class)
  B. Per-superfamily mean divergence across clades (heatmap)
  C. Recency Index per species per clade (boxplot)
  D. RI per clade × TE superfamily heatmap

Inputs:
  --annot_dir   Root of EarlGrey output (one subfolder per species)
  --clade_map   Two-column TSV: species_code<TAB>clade_tag
  --out_dir     Output directory

Outputs:
  kimura_species.tsv          – per-species RI and mean divergence
  kimura_superfamily.tsv      – per-species × superfamily summary
  figures/kimura_landscape_*.pdf   – per-species landscapes (PDF grid)
  figures/recency_index_boxplot.pdf
  figures/superfamily_divergence_heatmap.pdf
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

KIMURA_RECENT_THRESH = 0.10   # < this → "recent"

TE_CLASS_COLORS = {
    "LTR":             "#e41a1c",
    "DNA":             "#377eb8",
    "LINE":            "#4daf4a",
    "SINE":            "#984ea3",
    "Rolling Circle":  "#ff7f00",
    "Penelope":        "#a65628",
    "Unclassified":    "#999999",
}

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


# =============================================================================
# I/O
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annot_dir",  required=True,
                   help="Root directory containing one EarlGrey folder per species.")
    p.add_argument("--clade_map",  required=True)
    p.add_argument("--out_dir",    default="figures")
    p.add_argument("--bin_width",  type=float, default=1.0,
                   help="Kimura divergence bin width in %% (default: 1.0)")
    return p.parse_args()


def _te_class(family: str) -> str:
    """Extract broad TE class from RepeatMasker family string."""
    family = family.upper()
    if "LTR" in family:    return "LTR"
    if "DNA" in family:    return "DNA"
    if "LINE" in family:   return "LINE"
    if "SINE" in family:   return "SINE"
    if "RC/" in family or "HELITRON" in family: return "Rolling Circle"
    if "PENELOPE" in family: return "Penelope"
    return "Unclassified"


def _kimura_from_align(align_file: Path) -> pd.DataFrame:
    """
    Parse a RepeatMasker .align file and return a DataFrame with columns:
    te_family, te_class, kimura_div (as a fraction 0–1), bp_match.
    """
    rows = []
    opener = gzip.open if str(align_file).endswith(".gz") else open
    pat_header = re.compile(
        r"^\s*\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+"  # score, div, del, ins
        r"\S+\s+\d+\s+\d+\s+\(\d+\)\s+"            # query coords
        r"[+-]\s+"                                   # strand
        r"(\S+)\s+#(\S+)\s+"                        # family # class
        r"\d+\s+\d+\s+\(\d+\)"                      # repeat coords
    )
    pat_kimura = re.compile(r"^Kimura\s+\(with\s+divergence\)=\s*([\d.]+)")

    current = None
    with opener(align_file, "rt", errors="replace") as fh:
        for line in fh:
            m_h = pat_header.match(line)
            if m_h:
                current = {"family": m_h.group(1), "class": m_h.group(2)}
                continue
            if current:
                m_k = pat_kimura.match(line.strip())
                if m_k:
                    div = float(m_k.group(1)) / 100.0
                    rows.append({
                        "te_family": current["family"],
                        "te_class":  _te_class(current["class"]),
                        "kimura_div": div,
                    })
                    current = None

    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["te_family", "te_class", "kimura_div"])


def load_all_species(annot_dir: Path, species_list: list[str]) -> dict[str, pd.DataFrame]:
    data = {}
    for sp in species_list:
        # EarlGrey output structure: {sp}_EarlGrey/{sp}_RepeatMasker/{sp}.fasta.prep.align
        pattern = list((annot_dir / f"{sp}_EarlGrey").glob(
            f"*_RepeatMasker*/{sp}*.align")) + \
                  list((annot_dir / f"{sp}_EarlGrey").glob(
            f"*_RepeatMasker*/{sp}*.align.gz"))
        if not pattern:
            continue
        df = _kimura_from_align(pattern[0])
        if len(df):
            data[sp] = df
    print(f"  Loaded Kimura data for {len(data)} species")
    return data


# =============================================================================
# Summary statistics
# =============================================================================

def compute_species_stats(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for sp, df in data.items():
        total_bp = len(df)
        recent_bp = (df["kimura_div"] < KIMURA_RECENT_THRESH).sum()
        rows.append({
            "species":       sp,
            "total_te_loci": total_bp,
            "recent_loci":   recent_bp,
            "pct_recent":    recent_bp / total_bp if total_bp > 0 else 0.0,
            "mean_kimura":   df["kimura_div"].mean(),
            "median_kimura": df["kimura_div"].median(),
        })
    return pd.DataFrame(rows)


def compute_superfamily_stats(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for sp, df in data.items():
        for sf, grp in df.groupby("te_class"):
            rows.append({
                "species":    sp,
                "te_class":   sf,
                "n_loci":     len(grp),
                "mean_div":   grp["kimura_div"].mean(),
                "pct_recent": (grp["kimura_div"] < KIMURA_RECENT_THRESH).mean(),
            })
    return pd.DataFrame(rows)


# =============================================================================
# Figures
# =============================================================================

def plot_recency_boxplot(stats: pd.DataFrame, clade_map: dict,
                         out_dir: Path) -> None:
    stats = stats.copy()
    stats["clade"] = stats["species"].map(clade_map)
    clades = [c for c in CLADE_ORDER if c in stats["clade"].values]
    data   = [stats.loc[stats["clade"] == c, "pct_recent"].values for c in clades]
    colors = [CLADE_COLORS.get(c, "#999") for c in clades]

    fig, ax = plt.subplots(figsize=(13, 4))
    bp = ax.boxplot(data, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.5))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_xticks(range(1, len(clades) + 1))
    ax.set_xticklabels(clades, fontsize=9)
    ax.set_ylabel("Recency Index (fraction TE bp, Kimura < 0.10)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "recency_index_boxplot.pdf")
    fig.savefig(out_dir / "recency_index_boxplot.png", dpi=200)
    plt.close(fig)
    print("Saved: recency_index_boxplot.pdf")


def plot_landscape_grid(data: dict[str, pd.DataFrame], species_list: list[str],
                        bin_width: float, out_dir: Path,
                        max_per_pdf: int = 50) -> None:
    """One Kimura landscape histogram per species, tiled in a grid PDF."""
    ncols = 5
    batch_size = max_per_pdf
    species_list = [s for s in species_list if s in data]

    for batch_start in range(0, len(species_list), batch_size):
        batch = species_list[batch_start:batch_start + batch_size]
        nrows = (len(batch) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4 * ncols, 2.5 * nrows),
                                 sharex=True)
        axes = np.array(axes).flatten()

        bins = np.arange(0, 0.5 + bin_width / 100, bin_width / 100)

        for ax, sp in zip(axes, batch):
            df = data[sp]
            for te_class, grp in df.groupby("te_class"):
                color = TE_CLASS_COLORS.get(te_class, "#999")
                ax.hist(grp["kimura_div"], bins=bins, color=color,
                        alpha=0.7, label=te_class, density=False)
            ax.set_title(sp, fontsize=6)
            ax.tick_params(labelsize=5)
            ax.spines[["top", "right"]].set_visible(False)

        for ax in axes[len(batch):]:
            ax.set_visible(False)

        fig.tight_layout()
        suffix = f"_batch{batch_start // batch_size + 1}" if len(species_list) > batch_size else ""
        fig.savefig(out_dir / f"kimura_landscape_grid{suffix}.pdf")
        plt.close(fig)
        print(f"Saved: kimura_landscape_grid{suffix}.pdf")


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
    species_list = list(clade_map.keys())

    print(f"Loading .align files from {args.annot_dir}")
    data = load_all_species(Path(args.annot_dir), species_list)

    print("Computing per-species statistics ...")
    stats = compute_species_stats(data)
    stats["clade"] = stats["species"].map(clade_map)
    stats.to_csv(out_dir / "kimura_species.tsv", sep="\t", index=False)

    print("Computing per-class statistics ...")
    sf_stats = compute_superfamily_stats(data)
    sf_stats["clade"] = sf_stats["species"].map(clade_map)
    sf_stats.to_csv(out_dir / "kimura_superfamily.tsv", sep="\t", index=False)

    print("Generating Recency Index boxplot ...")
    plot_recency_boxplot(stats, clade_map, out_dir)

    print("Generating landscape grids ...")
    ordered = sorted(species_list, key=lambda s: CLADE_ORDER.index(
        clade_map.get(s, "AS")) if clade_map.get(s, "AS") in CLADE_ORDER else 99)
    plot_landscape_grid(data, ordered, args.bin_width, out_dir)

    print("Done.")


if __name__ == "__main__":
    main()
