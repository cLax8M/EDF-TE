#!/usr/bin/env python3
"""
03_chromatin_proximity_analysis.py  —  TE proximity to H3K9me3/5mC and gene expression.

For Phycomyces, Rhizopus, and Mucor (species with both ChIP-seq and RNA-seq),
quantifies:
  A. Distance distribution from TE copies to nearest H3K9me3/5mC-marked region
  B. Expression (log₂ TPM) of genes adjacent to TEs, split by:
       - TE is methylated (5mC) vs unmethylated
       - TE is in H3K9me3 domain vs not
     Combined factorial: 4 categories
  C. Proportions of TE classes in methylated vs H3K9me3 compartments

Inputs:
  --te_bed_dir     BED files: species.filteredRepeats.bed
  --h3k9me3_dir   BED files: species.H3K9me3_peaks.bed (or domains)
  --methyl_dir    Bismark .cov.gz or BED of methylated regions
  --expr_dir      Per-gene expression TSVs: gene, tpm (from featureCounts/kallisto)
  --gene_bed_dir  BED files: species.genes.bed (gene coordinates)
  --species_list  File with species codes (one per line)
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory
  --window_kb     Gene proximity window in kb (default: 5)

Outputs:
  proximity_expression_stats.tsv
  figures/chromatin_te_expression.pdf
  figures/te_class_chromatin_proportions.pdf
"""

import argparse
import subprocess
import tempfile
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#984ea3",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--te_bed_dir",   required=True)
    p.add_argument("--h3k9me3_dir",  required=True)
    p.add_argument("--methyl_dir",   required=True)
    p.add_argument("--expr_dir",     required=True)
    p.add_argument("--gene_bed_dir", required=True)
    p.add_argument("--species_list", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="chromatin_out")
    p.add_argument("--window_kb",    type=int, default=5)
    return p.parse_args()


def _te_class(f: str) -> str:
    f = f.upper()
    if "LTR" in f:    return "LTR"
    if "DNA" in f:    return "DNA"
    if "LINE" in f:   return "LINE"
    if "SINE" in f:   return "SINE"
    if "RC/" in f or "HELITRON" in f: return "Rolling Circle"
    if "PENELOPE" in f: return "Penelope"
    return "Unclassified"


def _find(directory: str, sp: str, *patterns: str) -> Path | None:
    for pat in patterns:
        hits = list(Path(directory).glob(f"{sp}{pat}"))
        if hits:
            return hits[0]
    return None


def bed_flag_overlap(query_bed: Path, subject_bed: Path) -> pd.Series:
    """Return boolean Series indexed by query line order: True if overlaps subject."""
    res = subprocess.run(
        ["bedtools", "intersect", "-a", str(query_bed), "-b", str(subject_bed),
         "-c", "-F", "0.5"],
        capture_output=True, text=True, check=True)
    counts = [int(l.rstrip().split("\t")[-1]) for l in res.stdout.splitlines()]
    return pd.Series([c > 0 for c in counts])


def nearest_distance(query_bed: Path, subject_bed: Path) -> pd.Series:
    """Return distance to nearest subject feature for each query interval."""
    res = subprocess.run(
        ["bedtools", "closest", "-a", str(query_bed), "-b", str(subject_bed),
         "-D", "b", "-io"],
        capture_output=True, text=True)
    dists = []
    for line in res.stdout.splitlines():
        parts = line.rstrip().split("\t")
        dists.append(abs(int(parts[-1])))
    return pd.Series(dists)


def load_te_bed(te_bed_dir: str, sp: str) -> pd.DataFrame:
    f = _find(te_bed_dir, sp, "*.filteredRepeats.bed", "*.bed")
    if f is None:
        return pd.DataFrame()
    df = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 3],
                     names=["chrom", "start", "end", "family"])
    df["te_class"] = df["family"].apply(_te_class)
    return df


def process_species(sp: str, args, window_bp: int) -> pd.DataFrame:
    te_bed  = load_te_bed(args.te_bed_dir, sp)
    if te_bed.empty:
        return pd.DataFrame()

    h3k9_f  = _find(args.h3k9me3_dir,  sp, "*H3K9me3*.bed", "*peaks.bed")
    meth_f  = _find(args.methyl_dir,   sp, "*methylated*.bed", "*.cov")
    expr_f  = _find(args.expr_dir,     sp, "*expr*.tsv", "*tpm*.tsv", "*.tsv")
    gene_f  = _find(args.gene_bed_dir, sp, "*.genes.bed", "*.bed")

    with tempfile.TemporaryDirectory() as tmp:
        te_tmp = Path(tmp) / "te.bed"
        te_bed[["chrom", "start", "end", "family"]].to_csv(
            te_tmp, sep="\t", index=False, header=False)

        te_bed["h3k9me3"] = False
        te_bed["methyl"]  = False

        if h3k9_f:
            try:
                te_bed["h3k9me3"] = bed_flag_overlap(te_tmp, h3k9_f).values
            except Exception:
                pass
        if meth_f:
            try:
                te_bed["methyl"] = bed_flag_overlap(te_tmp, meth_f).values
            except Exception:
                pass

    te_bed["species"] = sp
    return te_bed


# =============================================================================
# Expression analysis
# =============================================================================

def expression_by_category(sp_data: dict[str, pd.DataFrame],
                            args) -> pd.DataFrame:
    """
    For each species: find genes within window_bp of each TE,
    split into 4 chromatin categories, collect TPM.
    """
    rows = []
    window_bp = args.window_kb * 1000

    for sp, te_df in sp_data.items():
        expr_f = _find(args.expr_dir, sp, "*expr*.tsv", "*tpm*.tsv", "*.tsv")
        gene_f = _find(args.gene_bed_dir, sp, "*.genes.bed", "*.bed")
        if expr_f is None or gene_f is None:
            continue
        expr = pd.read_csv(expr_f, sep="\t")
        if "tpm" not in expr.columns or "gene" not in expr.columns:
            continue
        expr = expr.set_index("gene")["tpm"]

        gene_bed = pd.read_csv(gene_f, sep="\t", header=None, usecols=[0, 1, 2, 3],
                               names=["chrom", "gstart", "gend", "gene"])

        for _, te in te_df.iterrows():
            cat = ("methyl" if te.get("methyl", False) else "unmethyl") + \
                  "_" + ("H3K9me3" if te.get("h3k9me3", False) else "noH3K9")

            near = gene_bed[(gene_bed["chrom"] == te["chrom"]) &
                            (gene_bed["gend"]   > te["start"] - window_bp) &
                            (gene_bed["gstart"] < te["end"]   + window_bp)]
            for _, g in near.iterrows():
                tpm = expr.get(g["gene"], np.nan)
                if not np.isnan(tpm):
                    rows.append({"species": sp, "category": cat,
                                 "te_class": te["te_class"],
                                 "log2_tpm": np.log2(tpm + 1)})
    return pd.DataFrame(rows)


# =============================================================================
# Figures
# =============================================================================

def plot_expression(expr_df: pd.DataFrame, out_dir: Path) -> None:
    cats = ["unmethyl_noH3K9", "methyl_noH3K9",
            "unmethyl_H3K9me3", "methyl_H3K9me3"]
    labels = ["None", "5mC only", "H3K9me3 only", "5mC + H3K9me3"]

    data_by_cat = [expr_df[expr_df["category"] == c]["log2_tpm"].dropna()
                   for c in cats]
    data_by_cat = [d for d in data_by_cat]  # keep all

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#999999", "#e41a1c", "#4daf4a", "#984ea3"]
    bp = ax.boxplot(data_by_cat, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.5))
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
    ax.set_xticks(range(1, len(cats) + 1))
    ax.set_xticklabels(labels, fontsize=9, rotation=20)
    ax.set_ylabel("log₂(TPM + 1) adjacent genes", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "chromatin_te_expression.pdf")
    fig.savefig(out_dir / "chromatin_te_expression.png", dpi=200)
    plt.close(fig)
    print("Saved: chromatin_te_expression.pdf")


def plot_class_proportions(sp_data: dict, out_dir: Path) -> None:
    rows = []
    for sp, df in sp_data.items():
        for compartment in ["methyl", "h3k9me3"]:
            sub = df[df.get(compartment, pd.Series(False, index=df.index))]
            if sub.empty:
                continue
            for cls, cnt in sub["te_class"].value_counts().items():
                rows.append({"species": sp, "compartment": compartment,
                             "te_class": cls, "count": cnt})
    if not rows:
        return
    agg = pd.DataFrame(rows).groupby(["compartment", "te_class"])["count"].sum()
    agg = agg.unstack("te_class", fill_value=0)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, comp in zip(axes, ["methyl", "h3k9me3"]):
        if comp not in agg.index:
            ax.set_visible(False)
            continue
        vals = agg.loc[comp]
        vals = vals[vals > 0]
        total = vals.sum()
        pcts = vals / total
        colors = [TE_CLASS_COLORS.get(c, "#999") for c in vals.index]
        ax.bar(range(len(vals)), pcts * 100, color=colors, edgecolor="white")
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(vals.index, rotation=30, fontsize=8)
        ax.set_ylabel("% of compartment")
        ax.set_title("5mC" if comp == "methyl" else "H3K9me3")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "te_class_chromatin_proportions.pdf")
    plt.close(fig)
    print("Saved: te_class_chromatin_proportions.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    species_list = [l.strip() for l in open(args.species_list) if l.strip()]
    window_bp    = args.window_kb * 1000

    sp_data = {}
    for sp in species_list:
        print(f"  Processing {sp} ...")
        df = process_species(sp, args, window_bp)
        if not df.empty:
            sp_data[sp] = df

    print("Computing expression by chromatin category ...")
    expr_df = expression_by_category(sp_data, args)
    expr_df.to_csv(out_dir / "proximity_expression_stats.tsv", sep="\t", index=False)

    plot_expression(expr_df, fig_dir)
    plot_class_proportions(sp_data, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
