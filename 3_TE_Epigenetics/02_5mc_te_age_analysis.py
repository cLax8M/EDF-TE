#!/usr/bin/env python3

# analyzes relationship between TE age (Kimura divergence)

# and DNA methylation (5mC) from WGBS data


# imports
import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# clade color scheme 
CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

# TE class categories
TE_CLASSES = ["LTR", "DNA", "LINE", "SINE", "Rolling Circle", "Penelope"]

# methylation thresholds
MIN_CPGS    = 3
MIN_COVER   = 5
MIN_METHYL  = 50.0   # %


# parse arguments
def parse_args():
    # input configuration for methylation–TE age analysis
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--te_bed_dir",   required=True)
    p.add_argument("--wgbs_dir",     required=True)
    p.add_argument("--align_dir",    required=True)
    p.add_argument("--species_list", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="5mc_out")
    return p.parse_args()

# classify TE family into broad class
def _te_class(family: str) -> str:
    f = family.upper()
    if "LTR" in f:    return "LTR"
    if "DNA" in f:    return "DNA"
    if "LINE" in f:   return "LINE"
    if "SINE" in f:   return "SINE"
    if "RC/" in f or "HELITRON" in f: return "Rolling Circle"
    if "PENELOPE" in f: return "Penelope"
    return "Unclassified"


# load TE coordinates (BED)
def load_te_bed(te_bed_dir: Path, sp: str) -> pd.DataFrame:
    files = list(Path(te_bed_dir).glob(f"{sp}*.filteredRepeats.bed")) + \
            list(Path(te_bed_dir).glob(f"{sp}*.bed"))
    if not files:
        return pd.DataFrame()
    df = pd.read_csv(files[0], sep="\t", header=None, usecols=[0, 1, 2, 3],
                     names=["chrom", "start", "end", "family"])
    df["te_class"] = df["family"].apply(_te_class)
    df["te_id"]    = range(len(df))
    return df


# load WGBS CpG methylation calls
def load_wgbs(wgbs_dir: Path, sp: str) -> pd.DataFrame:
    """Load bismark .cov.gz file: chrom, start, end, methyl_pct, count_m, count_um."""
    files = list(Path(wgbs_dir).glob(f"{sp}*.bismark.cov.gz")) + \
            list(Path(wgbs_dir).glob(f"{sp}*.cov"))
    if not files:
        return pd.DataFrame()
    with gzip.open(files[0], "rt") if str(files[0]).endswith(".gz") else open(files[0]) as fh:
        df = pd.read_csv(fh, sep="\t", header=None, usecols=[0, 1, 4, 5],
                         names=["chrom", "pos", "count_m", "count_um"])
    df["cover"]      = df["count_m"] + df["count_um"]
    df["methyl_pct"] = 100.0 * df["count_m"] / df["cover"].replace(0, np.nan)
    return df[df["cover"] >= MIN_COVER].dropna(subset=["methyl_pct"])


# parse RepeatMasker alignments (Kimura divergence)
def load_kimura(align_dir: Path, sp: str) -> pd.DataFrame:
    files = list(Path(align_dir).glob(f"{sp}*_EarlGrey/*_RepeatMasker*/*.align")) + \
            list(Path(align_dir).glob(f"{sp}*.align")) + \
            list(Path(align_dir).glob(f"{sp}*.align.gz"))
    if not files:
        return pd.DataFrame()
    rows = []
    pat_h = re.compile(
        r"^\s*\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+(\S+)\s+(\d+)\s+(\d+).*?([+-])\s+(\S+)")
    pat_k = re.compile(r"^Kimura\s+\(with\s+divergence\)=\s*([\d.]+)")
    opener = gzip.open if str(files[0]).endswith(".gz") else open
    current = None
    with opener(files[0], "rt", errors="replace") as fh:
        for line in fh:
            m = pat_h.match(line)
            if m:
                current = {"chrom": m.group(1), "rm_start": int(m.group(2)),
                           "rm_end": int(m.group(3))}
                continue
            if current:
                mk = pat_k.match(line.strip())
                if mk:
                    current["kimura_div"] = float(mk.group(1)) / 100.0
                    rows.append(current)
                    current = None
    return pd.DataFrame(rows) if rows else pd.DataFrame()


# map CpG methylation onto TE copies
def intersect_cpgs_with_tes(te_bed: pd.DataFrame, wgbs: pd.DataFrame) -> pd.DataFrame:
    """
    For each TE, count CpGs covered at ≥5× and compute mean methylation.
    Returns te_bed with added columns: n_cpgs, mean_methyl, is_methylated.
    """
    results = []
    # Efficient: group CpGs by chrom
    wgbs_by_chrom = {c: g for c, g in wgbs.groupby("chrom")}

    for _, te in te_bed.iterrows():
        chrom_cpgs = wgbs_by_chrom.get(te["chrom"], pd.DataFrame())
        if chrom_cpgs.empty:
            results.append({"te_id": te["te_id"], "n_cpgs": 0,
                            "mean_methyl": np.nan, "is_methylated": False})
            continue
        sub = chrom_cpgs[(chrom_cpgs["pos"] >= te["start"]) &
                         (chrom_cpgs["pos"] <  te["end"])]
        n = len(sub)
        if n < MIN_CPGS:
            results.append({"te_id": te["te_id"], "n_cpgs": n,
                            "mean_methyl": np.nan, "is_methylated": False})
        else:
            m = sub["methyl_pct"].mean()
            results.append({"te_id": te["te_id"], "n_cpgs": n,
                            "mean_methyl": m,
                            "is_methylated": m >= MIN_METHYL})
    return te_bed.merge(pd.DataFrame(results), on="te_id")


# map Kimura divergence onto TE copies
def merge_kimura(te_cpg: pd.DataFrame, kimura: pd.DataFrame) -> pd.DataFrame:
    """Join Kimura divergence to TE copies by genomic overlap."""
    rows = []
    k_by_chrom = {c: g for c, g in kimura.groupby("chrom")}
    for _, te in te_cpg.iterrows():
        sub = k_by_chrom.get(te["chrom"], pd.DataFrame())
        if not sub.empty:
            ov = sub[(sub["rm_start"] < te["end"]) & (sub["rm_end"] > te["start"])]
            k_val = ov["kimura_div"].mean() if not ov.empty else np.nan
        else:
            k_val = np.nan
        rows.append(k_val)
    te_cpg = te_cpg.copy()
    te_cpg["kimura_div"] = rows
    return te_cpg



# Figures


# violin plot: methylation vs TE age
def plot_violin(data: pd.DataFrame, out_dir: Path) -> None:
    te_classes = [c for c in TE_CLASSES if c in data["te_class"].values]
    fig, axes = plt.subplots(1, len(te_classes), figsize=(3.5 * len(te_classes), 5),
                             sharey=True)
    if len(te_classes) == 1:
        axes = [axes]
    for ax, cls in zip(axes, te_classes):
        sub = data[data["te_class"] == cls].dropna(subset=["kimura_div"])
        meth   = sub[sub["is_methylated"]]["kimura_div"]
        unmeth = sub[~sub["is_methylated"]]["kimura_div"]
        if meth.empty or unmeth.empty:
            ax.set_title(f"{cls}\n(insufficient data)")
            continue
        parts = ax.violinplot([unmeth, meth], positions=[0, 1],
                              showmedians=True, widths=0.7)
        colors = ["#5aade0", "#e05a5a"]
        for pc, col in zip(parts["bodies"], colors):
            pc.set_facecolor(col)
            pc.set_alpha(0.8)
        parts["cmedians"].set_color("black")
        _, pval = mannwhitneyu(unmeth, meth, alternative="two-sided")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Unmethylated", "Methylated"], fontsize=8, rotation=30)
        ax.set_title(f"{cls}\np={pval:.2g}", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Kimura divergence", fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "5mc_te_age_violin.pdf")
    fig.savefig(out_dir / "5mc_te_age_violin.png", dpi=200)
    plt.close(fig)
    print("Saved: 5mc_te_age_violin.pdf")



# main


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
    print(f"Processing {len(species_list)} species with WGBS data ...")

    all_data = []
    for sp in species_list:
        print(f"  {sp} ...", end=" ")
        te_bed = load_te_bed(Path(args.te_bed_dir), sp)
        if te_bed.empty:
            print("no TE BED")
            continue
        wgbs = load_wgbs(Path(args.wgbs_dir), sp)
        if wgbs.empty:
            print("no WGBS")
            continue
        kimura = load_kimura(Path(args.align_dir), sp)

        te_cpg = intersect_cpgs_with_tes(te_bed, wgbs)
        if not kimura.empty:
            te_cpg = merge_kimura(te_cpg, kimura)
        te_cpg["species"] = sp
        te_cpg["clade"]   = clade_map.get(sp, "?")
        all_data.append(te_cpg)
        n_meth = te_cpg["is_methylated"].sum()
        print(f"{n_meth} methylated TEs")

    if not all_data:
        print("No data loaded — check inputs.", file=sys.stderr)
        sys.exit(1)

    data = pd.concat(all_data, ignore_index=True)
    data.to_csv(out_dir / "5mc_te_kimura.tsv", sep="\t", index=False)

    plot_violin(data, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
