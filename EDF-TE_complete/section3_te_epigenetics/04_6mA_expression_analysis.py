#!/usr/bin/env python3
"""
04_6mA_expression_analysis.py  —  N6-methyladenine (6mA) and TE expression.

Analyses the relationship between 6mA methylation at TE copies and their
transcriptional activity in species with PacBio SMRT-seq 6mA data.

Strategy:
  - A TE copy is "6mA-marked" if ≥1 motif-context adenine inside the copy
    is called as m6A with IPD ratio ≥ 5× background
  - Expression proxy: normalised coverage from RNA-seq mapped to TE copies
    (ChimeraTE or custom alignment to the TE library)
  - Violin plots: expression split by 6mA-marked vs unmarked, per TE class

Inputs:
  --6ma_dir       6mA calls per species (CSV: chrom, pos, ipd_ratio, strand)
  --te_bed_dir    EarlGrey BED files
  --rnacov_dir    Per-species RNA-seq coverage over TE copies
                  (BED or bedgraph: chrom, start, end, coverage)
  --species_list  File with species codes
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory
  --ipd_thresh    IPD ratio threshold (default: 5.0)

Outputs:
  6mA_te_expression.tsv
  figures/6mA_te_expression_violin.pdf
"""

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IPD_THRESH_DEFAULT = 5.0

TE_CLASSES = ["LTR", "DNA", "LINE", "SINE", "Rolling Circle", "Penelope"]

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
    p.add_argument("--6ma_dir",      required=True, dest="m6a_dir")
    p.add_argument("--te_bed_dir",   required=True)
    p.add_argument("--rnacov_dir",   required=True)
    p.add_argument("--species_list", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="6mA_out")
    p.add_argument("--ipd_thresh",   type=float, default=IPD_THRESH_DEFAULT)
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


def load_6ma(m6a_dir: str, sp: str, ipd_thresh: float) -> pd.DataFrame:
    files = list(Path(m6a_dir).glob(f"{sp}*.csv")) + \
            list(Path(m6a_dir).glob(f"{sp}*.6mA.bed")) + \
            list(Path(m6a_dir).glob(f"{sp}*.bed"))
    if not files:
        return pd.DataFrame()
    df = pd.read_csv(files[0], sep="\t" if str(files[0]).endswith((".bed", ".tsv")) else ",")
    df.columns = [c.lower().replace(" ", "_") for c in df.columns]
    if "ipd_ratio" in df.columns:
        df = df[df["ipd_ratio"] >= ipd_thresh]
    if "chrom" not in df.columns and "chromosome" in df.columns:
        df = df.rename(columns={"chromosome": "chrom"})
    if "pos" not in df.columns and "position" in df.columns:
        df = df.rename(columns={"position": "pos"})
    return df[["chrom", "pos"]].drop_duplicates()


def load_te_bed(te_bed_dir: str, sp: str) -> pd.DataFrame:
    files = list(Path(te_bed_dir).glob(f"{sp}*.filteredRepeats.bed")) + \
            list(Path(te_bed_dir).glob(f"{sp}*.bed"))
    if not files:
        return pd.DataFrame()
    df = pd.read_csv(files[0], sep="\t", header=None, usecols=[0, 1, 2, 3],
                     names=["chrom", "start", "end", "family"])
    df["te_class"] = df["family"].apply(_te_class)
    df["te_id"]    = range(len(df))
    return df


def load_rnacov(rnacov_dir: str, sp: str) -> pd.DataFrame:
    files = list(Path(rnacov_dir).glob(f"{sp}*.bedgraph")) + \
            list(Path(rnacov_dir).glob(f"{sp}*.coverage.bed")) + \
            list(Path(rnacov_dir).glob(f"{sp}*.bed"))
    if not files:
        return pd.DataFrame()
    df = pd.read_csv(files[0], sep="\t", header=None, usecols=[0, 1, 2, 3],
                     names=["chrom", "start", "end", "coverage"])
    return df


def mark_6ma_tes(te_bed: pd.DataFrame, m6a: pd.DataFrame) -> pd.DataFrame:
    """Mark TE copies that contain at least one 6mA site."""
    m6a_by_chrom = {c: g["pos"].values for c, g in m6a.groupby("chrom")}
    flags = []
    for _, te in te_bed.iterrows():
        positions = m6a_by_chrom.get(te["chrom"], np.array([]))
        n_6mA = int(((positions >= te["start"]) & (positions < te["end"])).sum())
        flags.append(n_6mA > 0)
    te_bed = te_bed.copy()
    te_bed["has_6mA"] = flags
    return te_bed


def mean_rnacov_for_tes(te_bed: pd.DataFrame, rna: pd.DataFrame) -> pd.DataFrame:
    """Attach mean RNA-seq coverage to each TE copy."""
    rna_by_chrom = {c: g for c, g in rna.groupby("chrom")}
    covs = []
    for _, te in te_bed.iterrows():
        rna_c = rna_by_chrom.get(te["chrom"], pd.DataFrame())
        if rna_c.empty:
            covs.append(np.nan)
            continue
        ov = rna_c[(rna_c["end"] > te["start"]) & (rna_c["start"] < te["end"])]
        if ov.empty:
            covs.append(0.0)
        else:
            # Weighted mean by bp overlap
            ov = ov.copy()
            ov["ov_len"] = ov["end"].clip(upper=te["end"]) - ov["start"].clip(lower=te["start"])
            covs.append((ov["coverage"] * ov["ov_len"]).sum() / ov["ov_len"].sum())
    te_bed = te_bed.copy()
    te_bed["rna_cov"] = covs
    return te_bed


# =============================================================================
# Figures
# =============================================================================

def plot_violin(data: pd.DataFrame, out_dir: Path) -> None:
    te_classes = [c for c in TE_CLASSES if c in data["te_class"].values]
    if not te_classes:
        print("No TE classes in data — skipping violin.")
        return

    fig, axes = plt.subplots(1, len(te_classes),
                             figsize=(3.5 * len(te_classes), 5), sharey=True)
    if len(te_classes) == 1:
        axes = [axes]

    for ax, cls in zip(axes, te_classes):
        sub = data[(data["te_class"] == cls)].dropna(subset=["rna_cov"])
        marked   = np.log2(sub[sub["has_6mA"]]["rna_cov"] + 1)
        unmarked = np.log2(sub[~sub["has_6mA"]]["rna_cov"] + 1)
        if marked.empty or unmarked.empty:
            ax.set_title(f"{cls}\n(insufficient data)", fontsize=8)
            continue
        parts = ax.violinplot([unmarked, marked], positions=[0, 1],
                              showmedians=True, widths=0.7)
        for pc, col in zip(parts["bodies"], ["#5aade0", "#e05a5a"]):
            pc.set_facecolor(col)
            pc.set_alpha(0.8)
        parts["cmedians"].set_color("black")
        _, pval = mannwhitneyu(unmarked, marked, alternative="two-sided")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["No 6mA", "6mA"], fontsize=8, rotation=15)
        ax.set_title(f"{cls}\np={pval:.2g}", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_ylabel("log₂(RNA-seq coverage + 1)", fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "6mA_te_expression_violin.pdf")
    fig.savefig(out_dir / "6mA_te_expression_violin.png", dpi=200)
    plt.close(fig)
    print("Saved: 6mA_te_expression_violin.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))

    species_list = [l.strip() for l in open(args.species_list) if l.strip()]
    all_data = []

    for sp in species_list:
        print(f"  {sp} ...", end=" ")
        m6a   = load_6ma(args.m6a_dir, sp, args.ipd_thresh)
        te    = load_te_bed(args.te_bed_dir, sp)
        rna   = load_rnacov(args.rnacov_dir, sp)
        if m6a.empty or te.empty:
            print("missing data")
            continue
        te = mark_6ma_tes(te, m6a)
        if not rna.empty:
            te = mean_rnacov_for_tes(te, rna)
        te["species"] = sp
        te["clade"]   = clade_map.get(sp, "?")
        all_data.append(te)
        print(f"{te['has_6mA'].sum()} 6mA-marked TEs")

    if not all_data:
        print("No data loaded.", file=sys.stderr)
        sys.exit(1)

    data = pd.concat(all_data, ignore_index=True)
    data.to_csv(out_dir / "6mA_te_expression.tsv", sep="\t", index=False)
    plot_violin(data, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
