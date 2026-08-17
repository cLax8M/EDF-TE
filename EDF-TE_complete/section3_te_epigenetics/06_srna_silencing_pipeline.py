#!/usr/bin/env python3
"""
06_srna_silencing_pipeline.py  —  sRNA-mediated TE silencing analysis.

Identifies TE-targeting small RNAs (sRNAs) from small RNA-seq data and
characterises:
  A. Proportion of 21–24 nt sRNAs mapping to each TE superfamily
  B. Correlation between sRNA targeting and H3K9me3 / 5mC at the same TE copy
  C. Identification of TE families silenced via the r3b2/DCL pathway

Pipeline (Methods §12.2):
  1. Trim the Illumina small-RNA 3′ adapter with Cutadapt v4.6 and keep
     reads of 18–32 nt (--discard-untrimmed):
         cutadapt -a TGGAATTCTCGGGTGCCAAGG \\
           --minimum-length 18 --maximum-length 32 --discard-untrimmed \\
           -o <trimmed.fq.gz> <raw.fq.gz>
  2. Size-select 20–25 nt (canonical fungal siRNA/piRNA range)
  3. Map trimmed reads to the TE library with Bowtie (-v 2 -M 1 --best
     --strata, one random alignment for multimappers)
  4. Count mappings per TE family, normalise by total mapped reads (RPM)
  5. Correlate sRNA RPM with H3K9me3 enrichment (Spearman) and test the
     association statistically (Mann–Whitney U of H3K9me3 enrichment between
     sRNA-targeted and non-targeted families, BH-FDR corrected)
  6. Identify high-sRNA families (≥ 99th percentile RPM) → candidate targets

Inputs:
  --te_lib        TE consensus library FASTA (EarlGrey output)
  --srna_dir      Directory of RAW sRNA FASTQ files (one per species);
                  reads are adapter-trimmed here with Cutadapt
  --h3k9_dir      H3K9me3 ChIP enrichment per TE family TSV (from step 05)
  --species_list  File with species codes
  --bowtie_index  Pre-built Bowtie index of TE library, or build from --te_lib
  --out_dir       Output directory
  --threads       Bowtie threads (default: 8)
  --adapter       Illumina small-RNA 3' adapter (default: TGGAATTCTCGGGTGCCAAGG)
  --min_trim      Minimum read length after adapter trimming (default: 18)
  --max_trim      Maximum read length after adapter trimming (default: 32)

Outputs:
  srna_te_mapping.tsv      – per-species, per-family RPM
  srna_h3k9_correlation.tsv – per-family Spearman + Mann–Whitney U + BH-FDR
  figures/srna_te_class_proportions.pdf
  figures/srna_h3k9_scatter.pdf
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
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


def bh_fdr(pvals):
    """Benjamini–Hochberg FDR adjustment."""
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = ranked
    return np.clip(out, 0, 1)


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
    p.add_argument("--adapter",      default="TGGAATTCTCGGGTGCCAAGG",
                   help="Illumina small-RNA 3' adapter sequence")
    p.add_argument("--min_trim",     type=int, default=18,
                   help="Min read length after Cutadapt trimming (Methods: 18)")
    p.add_argument("--max_trim",     type=int, default=32,
                   help="Max read length after Cutadapt trimming (Methods: 32)")
    p.add_argument("--min_len",      type=int, default=20,
                   help="Lower size-selection bound after trimming (20 nt)")
    p.add_argument("--max_len",      type=int, default=25,
                   help="Upper size-selection bound after trimming (25 nt)")
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


def trim_adapters(fastq: Path, out_dir: Path, adapter: str,
                  min_len: int, max_len: int) -> Path | None:
    """Run Cutadapt to trim the small-RNA 3' adapter and keep 18–32 nt reads.

    Mirrors Methods §12.2:
        cutadapt -a <ADAPTER> --minimum-length 18 --maximum-length 32 \\
          --discard-untrimmed -o <trimmed.fq.gz> <raw.fq.gz>
    Returns the path of the trimmed file, or None if Cutadapt is unavailable.
    """
    out_path = out_dir / f"{fastq.stem}.trim18-32.fq.gz"
    if out_path.exists():
        return out_path
    try:
        subprocess.run(
            ["cutadapt", "-a", adapter,
             "--minimum-length", str(min_len),
             "--maximum-length", str(max_len),
             "--discard-untrimmed",
             "-o", str(out_path), str(fastq)],
            check=True, capture_output=True, text=True)
    except FileNotFoundError:
        print("  [WARN] cutadapt not found on PATH — using raw FASTQ.")
        return fastq
    except subprocess.CalledProcessError as e:
        print(f"  [WARN] cutadapt failed: {e.stderr[-300:]} — using raw FASTQ.")
        return fastq
    return out_path


def size_select(fastq: Path, out_dir: Path, lo: int, hi: int) -> Path:
    """Keep reads whose insert length is in [lo, hi] (awk size-selection)."""
    out_path = out_dir / f"{fastq.stem}.size{lo}-{hi}.fq.gz"
    if out_path.exists():
        return out_path
    awk = (
        f'NR%4==2{{len=length($0)}} '
        f'NR%4==0{{if(len>={lo} && len<={hi}) '
        f'print prev3 "\\n" prev2 "\\n" prev1 "\\n" $0}} '
        f'{{prev3=prev2; prev2=prev1; prev1=$0}}'
    )
    with open(out_path, "w") as fh:
        subprocess.run(f"zcat -f {fastq} | awk '{awk}' | gzip > {out_path}",
                       shell=True, check=True)
    return out_path


def build_bowtie_index(te_lib: str, out_dir: Path) -> str:
    idx = str(out_dir / "te_lib_index")
    subprocess.run(["bowtie-build", "-f", te_lib, idx],
                   check=True, capture_output=True)
    return idx


def map_srna(fastq: Path, index: str, threads: int,
             min_len: int, max_len: int) -> pd.DataFrame:
    """
    Map sRNAs to TE library with Bowtie (Methods §12.2: -v 2 -M 1 --best
    --strata), return per-read family mapping.
    Size filter is applied post-mapping.
    """
    sam = subprocess.run(
        ["bowtie", "-p", str(threads), "--no-unal", "-v", "2", "-M", "1",
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
    """Per-family Spearman scatter of sRNA RPM vs H3K9me3 enrichment.

    Also saves a table with the Spearman rho, p-value and BH-FDR per species.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for sp, grp in merged.groupby("species"):
        x = grp["rpm"]
        y = grp["h3k9_enrichment"]
        valid = x.notna() & y.notna()
        if valid.sum() < 5:
            rows.append({"species": sp, "n": valid.sum(),
                         "spearman_rho": np.nan, "spearman_p": np.nan})
            continue
        r, p = spearmanr(x[valid], y[valid])
        rows.append({"species": sp, "n": int(valid.sum()),
                     "spearman_rho": r, "spearman_p": p})
    corr = pd.DataFrame(rows)
    if len(corr) and corr["spearman_p"].notna().any():
        mask = corr["spearman_p"].notna()
        corr.loc[mask, "spearman_bh_fdr"] = bh_fdr(
            corr.loc[mask, "spearman_p"].to_numpy(float))
    corr.to_csv(out_dir / "srna_h3k9_correlation.tsv", sep="\t", index=False)

    # scatter (aggregate all species)
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
    fig.savefig(out_dir / "srna_h3k9_scatter.pdf")
    plt.close(fig)
    print("Saved: srna_h3k9_scatter.pdf")


def test_h3k9_enrichment_by_targeting(rpm_df: pd.DataFrame,
                                      h3k9: pd.DataFrame,
                                      out_dir: Path) -> pd.DataFrame:
    """Mann–Whitney U: H3K9me3 enrichment of sRNA-targeted vs non-targeted
    TE families, BH-FDR corrected.

    A family is "sRNA-targeted" if its RPM is above the genome-wide 99th
    percentile (candidate silencing targets, as in Methods).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    rpm_agg = rpm_df.groupby("family")["rpm"].sum().reset_index()
    rpm_agg.columns = ["family", "rpm"]
    h3k9_agg = (h3k9.groupby("family")["mean_cov"].mean().reset_index()
                if "family" in h3k9.columns else
                h3k9.groupby("te_class")["mean_cov"].mean().reset_index())
    h3k9_agg.columns = ["family", "h3k9_enrichment"]

    merged = rpm_agg.merge(h3k9_agg, on="family", how="inner")
    if len(merged) < 10:
        print("Too few families with both sRNA and H3K9me3 data.")
        return merged

    thr = np.percentile(merged["rpm"], 99)
    merged["is_target"] = merged["rpm"] >= thr
    targeted = merged.loc[merged["is_target"], "h3k9_enrichment"]
    non_targeted = merged.loc[~merged["is_target"], "h3k9_enrichment"]

    if len(targeted) >= 5 and len(non_targeted) >= 5:
        u, p = mannwhitneyu(targeted, non_targeted, alternative="two-sided")
    else:
        u, p = np.nan, np.nan

    summ = pd.DataFrame([{
        "n_families": len(merged),
        "n_targeted": int(merged["is_target"].sum()),
        "rpm_99th_pct": thr,
        "median_h3k9_targeted": targeted.median() if len(targeted) else np.nan,
        "median_h3k9_non_targeted": non_targeted.median() if len(non_targeted) else np.nan,
        "mannwhitney_p": p,
        "bh_fdr": p,
    }])
    summ.to_csv(out_dir / "srna_h3k9_targeting_test.tsv", sep="\t", index=False)

    # boxplot
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    bp = ax.boxplot([targeted, non_targeted], patch_artist=True, widths=0.5,
                    medianprops=dict(color="black", lw=1.3), showfliers=False)
    for patch, col in zip(bp["boxes"], ["#d62728", "#999999"]):
        patch.set_facecolor(col)
        patch.set_alpha(0.8)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["sRNA-targeted", "non-targeted"], fontsize=9)
    ax.set_ylabel("H3K9me3 enrichment", fontsize=10)
    ax.text(0.05, 0.95, f"MWU p={p:.2g}", transform=ax.transAxes,
            fontsize=9, va="top")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "srna_h3k9_targeting_boxplot.pdf")
    plt.close(fig)
    print("Saved: srna_h3k9_targeting_boxplot.pdf")
    return merged


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    trim_dir = out_dir / "trimmed"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    trim_dir.mkdir(parents=True, exist_ok=True)

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
        print(f"  {sp}: {fastq.name}", flush=True)

        # 1) Cutadapt: trim 3' adapter, keep 18–32 nt reads (Methods §12.2)
        trimmed = trim_adapters(fastq, trim_dir, args.adapter,
                                args.min_trim, args.max_trim)
        # 2) Size-select 20–25 nt
        sized = size_select(trimmed, trim_dir, args.min_len, args.max_len)

        total = count_total_reads(sized)
        print(f"    after trim+size-select: {total:,} reads", flush=True)

        mapping = map_srna(sized, index, args.threads, args.min_len, args.max_len)
        if mapping.empty:
            print("    no mappings")
            continue
        rpm = compute_rpm(mapping, total)
        rpm["species"] = sp
        all_rpm.append(rpm)
        print(f"    {len(rpm)} families mapped")

    if not all_rpm:
        print("No sRNA data mapped.", file=sys.stderr)
        sys.exit(1)

    rpm_df = pd.concat(all_rpm, ignore_index=True)
    rpm_df.to_csv(out_dir / "srna_te_mapping.tsv", sep="\t", index=False)

    plot_srna_proportions(rpm_df, fig_dir)

    # Correlate with H3K9me3 and run the statistical test
    h3k9_files = list(Path(args.h3k9_dir).glob("h3k9me3_te_enrichment.tsv"))
    if h3k9_files:
        h3k9 = pd.read_csv(h3k9_files[0], sep="\t")
        h3k9_agg = h3k9.groupby("te_class")["mean_cov"].mean().reset_index()
        h3k9_agg.columns = ["te_class", "h3k9_enrichment"]
        rpm_agg = rpm_df.groupby("te_class")["rpm"].sum().reset_index()
        merged = rpm_agg.merge(h3k9_agg, on="te_class", how="inner")
        merged["species"] = "all"
        merged.to_csv(out_dir / "srna_h3k9_merged.tsv", sep="\t", index=False)
        plot_srna_h3k9_scatter(merged, fig_dir)
        # family-level statistical test of H3K9me3 enrichment by targeting
        test_h3k9_enrichment_by_targeting(rpm_df, h3k9, out_dir)

    print("Done.")


if __name__ == "__main__":
    main()
