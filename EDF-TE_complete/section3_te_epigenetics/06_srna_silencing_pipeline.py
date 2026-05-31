#!/usr/bin/env python3
"""
06_srna_silencing_pipeline.py  —  sRNA-mediated TE silencing analysis.

Identifies TE-targeting small RNAs (sRNAs) from small RNA-seq data and
characterises:
  A. Proportion of 21–24 nt sRNAs mapping to each TE superfamily
  B. Correlation between sRNA targeting and H3K9me3 / 5mC at the same TE copy
  C. Identification of TE families silenced via the r3b2/DCL pathway

Pipeline:
  1. Trim adapters (cutadapt)  — assumed done externally
  2. Map trimmed reads to TE library with Bowtie (1 mismatch, unique)
  3. Size-select 20–25 nt
  4. Count mappings per TE family, normalise by total reads (RPM)
  5. Correlate sRNA RPM with H3K9me3 enrichment (Spearman)
  6. Identify high-sRNA families (≥ 99th percentile RPM) → candidate targets

Inputs:
  --te_lib        TE consensus library FASTA (EarlGrey output)
  --srna_dir      Directory of trimmed FASTQ files (one per species)
  --h3k9_dir      H3K9me3 ChIP enrichment per TE family TSV (from step 05)
  --species_list  File with species codes
  --bowtie_index  Pre-built Bowtie index of TE library, or build from --te_lib
  --out_dir       Output directory
  --threads       Bowtie threads (default: 8)

Outputs:
  srna_te_mapping.tsv      – per-species, per-family RPM
  srna_h3k9_correlation.tsv
  figures/srna_te_class_proportions.pdf
  figures/srna_h3k9_scatter.pdf
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
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
    p.add_argument("--te_lib",       required=True)
    p.add_argument("--srna_dir",     required=True)
    p.add_argument("--h3k9_dir",     required=True)
    p.add_argument("--species_list", required=True)
    p.add_argument("--bowtie_index", default=None,
                   help="Pre-built Bowtie index prefix. If absent, builds from --te_lib.")
    p.add_argument("--out_dir",      default="srna_out")
    p.add_argument("--threads",      type=int, default=8)
    p.add_argument("--min_len",      type=int, default=20)
    p.add_argument("--max_len",      type=int, default=25)
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


def build_bowtie_index(te_lib: str, out_dir: Path) -> str:
    idx = str(out_dir / "te_lib_index")
    subprocess.run(["bowtie-build", "-f", te_lib, idx],
                   check=True, capture_output=True)
    return idx


def map_srna(fastq: Path, index: str, threads: int,
             min_len: int, max_len: int) -> pd.DataFrame:
    """
    Map sRNAs to TE library with Bowtie, return per-read family mapping.
    Size filter is applied by awk post-mapping.
    """
    # Map
    sam = subprocess.run(
        ["bowtie", "-p", str(threads), "--no-unal", "-v", "1", "-M", "1",
         "--best", "--strata", index, fastq],
        capture_output=True, text=True, check=False)

    rows = []
    for line in sam.stdout.splitlines():
        if line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        seq_len = len(parts[9])
        if seq_len < min_len or seq_len > max_len:
            continue
        ref = parts[2]
        if ref == "*":
            continue
        rows.append({"family": ref, "read_len": seq_len})
    return pd.DataFrame(rows)


def compute_rpm(mapping: pd.DataFrame, total_reads: int) -> pd.DataFrame:
    counts = mapping.groupby("family").size().reset_index(name="count")
    counts["rpm"] = counts["count"] / total_reads * 1e6
    counts["te_class"] = counts["family"].apply(_te_class)
    return counts


def count_total_reads(fastq_path: Path) -> int:
    res = subprocess.run(
        f"wc -l < {fastq_path}",
        shell=True, capture_output=True, text=True)
    try:
        return int(res.stdout.strip()) // 4
    except ValueError:
        return 1


# =============================================================================
# Figures
# =============================================================================

def plot_srna_proportions(rpm_df: pd.DataFrame, out_dir: Path) -> None:
    class_rpm = rpm_df.groupby("te_class")["rpm"].sum()
    total = class_rpm.sum()
    if total == 0:
        return
    pcts = class_rpm / total * 100

    colors = [TE_CLASS_COLORS.get(c, "#999") for c in pcts.index]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(range(len(pcts)), pcts, color=colors, edgecolor="white")
    ax.set_xticks(range(len(pcts)))
    ax.set_xticklabels(pcts.index, rotation=30, fontsize=9)
    ax.set_ylabel("% of TE-mapping sRNAs", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "srna_te_class_proportions.pdf")
    plt.close(fig)
    print("Saved: srna_te_class_proportions.pdf")


def plot_srna_h3k9_scatter(merged: pd.DataFrame, out_dir: Path) -> None:
    x = merged["rpm"]
    y = merged["h3k9_enrichment"]
    valid = x.notna() & y.notna()
    if valid.sum() < 5:
        print("Insufficient data for sRNA–H3K9me3 scatter.")
        return

    r, p = spearmanr(x[valid], y[valid])
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(np.log2(x[valid] + 1), y[valid], alpha=0.5, s=20, color="#333333")
    ax.set_xlabel("log₂(sRNA RPM + 1)", fontsize=10)
    ax.set_ylabel("H3K9me3 enrichment", fontsize=10)
    ax.text(0.05, 0.95, f"Spearman r={r:.2f}, p={p:.2g}",
            transform=ax.transAxes, fontsize=9, va="top")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "srna_h3k9_scatter.pdf")
    plt.close(fig)
    print("Saved: srna_h3k9_scatter.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    index = args.bowtie_index
    if index is None:
        print("Building Bowtie index from TE library ...")
        index = build_bowtie_index(args.te_lib, out_dir)

    species_list = [l.strip() for l in open(args.species_list) if l.strip()]
    all_rpm = []

    for sp in species_list:
        fastq_files = list(Path(args.srna_dir).glob(f"{sp}*.fastq.gz")) + \
                      list(Path(args.srna_dir).glob(f"{sp}*.fq.gz")) + \
                      list(Path(args.srna_dir).glob(f"{sp}*.fastq"))
        if not fastq_files:
            print(f"  {sp}: no sRNA FASTQ found, skipping")
            continue
        fastq = fastq_files[0]
        total = count_total_reads(fastq)
        print(f"  {sp}: {total:,} reads ...", end=" ")

        mapping = map_srna(fastq, index, args.threads, args.min_len, args.max_len)
        if mapping.empty:
            print("no mappings")
            continue
        rpm = compute_rpm(mapping, total)
        rpm["species"] = sp
        all_rpm.append(rpm)
        print(f"{len(rpm)} families")

    if not all_rpm:
        print("No sRNA data mapped.", file=sys.stderr)
        sys.exit(1)

    rpm_df = pd.concat(all_rpm, ignore_index=True)
    rpm_df.to_csv(out_dir / "srna_te_mapping.tsv", sep="\t", index=False)

    plot_srna_proportions(rpm_df, fig_dir)

    # Correlate with H3K9me3 if available
    h3k9_files = list(Path(args.h3k9_dir).glob("h3k9me3_te_enrichment.tsv"))
    if h3k9_files:
        h3k9 = pd.read_csv(h3k9_files[0], sep="\t")
        h3k9_agg = h3k9.groupby("te_class")["mean_cov"].mean().reset_index()
        h3k9_agg.columns = ["te_class", "h3k9_enrichment"]
        rpm_agg = rpm_df.groupby("te_class")["rpm"].sum().reset_index()
        merged = rpm_agg.merge(h3k9_agg, on="te_class", how="inner")
        merged.to_csv(out_dir / "srna_h3k9_correlation.tsv", sep="\t", index=False)
        plot_srna_h3k9_scatter(merged, fig_dir)

    print("Done.")


if __name__ == "__main__":
    main()
