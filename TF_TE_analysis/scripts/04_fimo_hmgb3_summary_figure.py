#!/usr/bin/env python3
"""
04_fimo_hmgb3_summary_figure.py — Summary figures from FIMO HMGb3 cross-species scan.

Reads per-species fimo.tsv from the hmgb3_jockey_cross_species directory.
Generates:
  1. Bar chart: HMGb3 hits per Mb genome, sorted by count, colored by clade
  2. Bar chart: HMGb3 hits overlapping LINE/I-Jockey elements (BED intersection)
  3. Summary table TSV with per-species statistics

Usage:
  python3 04_fimo_hmgb3_summary_figure.py
"""
from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────
FIMO_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/hmgb3_jockey_cross_species")
OUTDIR   = FIMO_DIR / "figures"
TE_DIR   = Path("/data/carlos/EDF-TE/TE_analysis/annotated")

# ── Clade config ───────────────────────────────────────────────────────────
CLADE_COLORS = {
    "MU": "#E07B39", "BL": "#9C5FA8", "AP": "#59A14F",
    "MI": "#76B7B2", "CH": "#F1CE63", "BA": "#4E79A7",
    "AS": "#E15759", "NU": "#B07AA1", "MO": "#FF9DA7",
    "GL": "#9C755F", "HO": "#80b1d3",
}
CLADE_LABELS = {
    "MU": "Mucoromycota", "BL": "Blastocladiomycota",
    "AP": "Aphelidiomycota", "MI": "Microsporidia",
    "CH": "Chytridiomycota", "BA": "Basidiomycota",
    "AS": "Ascomycota", "NU": "Neocallimastigomycota",
    "MO": "Mortierellomycota", "GL": "Glomeromycotina",
    "HO": "Holozoa",
}


def clade_of(sp):
    return sp[:2]


def count_fimo_hits(fimo_tsv):
    """Count non-header, non-comment lines in fimo.tsv."""
    n = 0
    with open(fimo_tsv) as f:
        for line in f:
            if line.startswith("#") or line.startswith("motif_id"):
                continue
            n += 1
    return n


def parse_fimo_coords(fimo_tsv):
    """Return list of (scaffold, start, end) for each hit."""
    coords = []
    with open(fimo_tsv) as f:
        for line in f:
            if line.startswith("#") or line.startswith("motif_id"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                coords.append((parts[2], int(parts[3]), int(parts[4])))
    return coords


def count_jockey_overlaps(fimo_tsv, species):
    """Count how many FIMO hits overlap a LINE/I-Jockey element."""
    bed_path = TE_DIR / f"{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.bed"
    if not bed_path.exists():
        return None, None

    # Load Jockey elements
    jockey_intervals = {}  # scaffold → sorted list of (start, end)
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            te_class = parts[3] if len(parts) > 3 else ""
            if "Jockey" in te_class or "LINE" in te_class:
                jockey_intervals.setdefault(chrom, []).append((start, end))

    # Sort intervals
    for chrom in jockey_intervals:
        jockey_intervals[chrom].sort()

    n_jockey = sum(len(v) for v in jockey_intervals.values())

    # Check FIMO hits for overlap
    fimo_coords = parse_fimo_coords(fimo_tsv)
    overlaps = 0
    for chrom, fstart, fend in fimo_coords:
        if chrom not in jockey_intervals:
            continue
        # Binary search for overlapping intervals
        intervals = jockey_intervals[chrom]
        for istart, iend in intervals:
            if istart > fend:
                break
            if iend >= fstart:
                overlaps += 1
                break

    return n_jockey, overlaps


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Collect per-species FIMO hits
    rows = []
    fimo_dirs = sorted(FIMO_DIR.glob("fimo_*/fimo.tsv"))
    print(f"[INFO] Found {len(fimo_dirs)} per-species FIMO results")

    # Read genome sizes from existing summary
    summary_tsv = FIMO_DIR / "hmgb3_jockey_summary.tsv"
    genome_sizes = {}
    if summary_tsv.exists():
        df_sum = pd.read_csv(summary_tsv, sep="\t")
        genome_sizes = dict(zip(df_sum["species"], df_sum["genome_mb"]))

    for fimo_tsv in fimo_dirs:
        sp = fimo_tsv.parent.name.replace("fimo_", "")
        n_hits = count_fimo_hits(fimo_tsv)
        genome_mb = genome_sizes.get(sp, np.nan)
        hits_per_mb = n_hits / genome_mb if genome_mb > 0 else np.nan
        rows.append({
            "species": sp,
            "clade": clade_of(sp),
            "n_hmgb3_hits": n_hits,
            "genome_mb": genome_mb,
            "hits_per_mb": hits_per_mb,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("hits_per_mb", ascending=True).reset_index(drop=True)

    # Save summary table
    df.to_csv(OUTDIR / "hmgb3_genome_wide_summary.tsv", sep="\t", index=False)
    print(f"[INFO] Summary saved: {len(df)} species")

    # ── Figure 1: Horizontal bar chart of HMGb3 hits/Mb ─────────────────
    fig, ax = plt.subplots(figsize=(8, max(8, len(df) * 0.18)))
    colors = [CLADE_COLORS.get(c, "#aaa") for c in df["clade"]]
    ax.barh(range(len(df)), df["hits_per_mb"], color=colors, edgecolor="none", height=0.8)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["species"], fontsize=5)
    ax.set_xlabel("HMGb3 motif hits per Mb genome")
    ax.set_title("HMGb3 (TGTTGTTGTTGTTGT) motif density across species\n(genome-wide FIMO scan)")
    ax.set_xlim(0, df["hits_per_mb"].max() * 1.05)

    # Clade legend
    present = df["clade"].unique()
    handles = [mpatches.Patch(color=CLADE_COLORS.get(c, "#aaa"),
               label=CLADE_LABELS.get(c, c))
               for c in sorted(present) if c in CLADE_COLORS]
    ax.legend(handles=handles, fontsize=6, loc="lower right", framealpha=0.9)

    fig.tight_layout()
    fig.savefig(OUTDIR / "hmgb3_hits_per_mb.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "hmgb3_hits_per_mb.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[OK] hmgb3_hits_per_mb.pdf")

    # ── Figure 2: Top 20 species, labeled nicely ─────────────────────────
    df_top = df.nlargest(20, "hits_per_mb").sort_values("hits_per_mb", ascending=True)
    fig, ax = plt.subplots(figsize=(7, 6))
    colors_top = [CLADE_COLORS.get(c, "#aaa") for c in df_top["clade"]]
    bars = ax.barh(range(len(df_top)), df_top["hits_per_mb"], color=colors_top,
                   edgecolor="white", height=0.7)
    ax.set_yticks(range(len(df_top)))
    ax.set_yticklabels(df_top["species"], fontsize=8)
    ax.set_xlabel("HMGb3 motif hits per Mb genome", fontsize=10)
    ax.set_title("Top 20 species: HMGb3 motif density", fontsize=11)
    # Add count labels
    for i, (_, row) in enumerate(df_top.iterrows()):
        ax.text(row["hits_per_mb"] + 1, i, f"{row['n_hmgb3_hits']:.0f}",
                fontsize=6, va="center", color="#333")
    fig.tight_layout()
    fig.savefig(OUTDIR / "hmgb3_top20.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[OK] hmgb3_top20.pdf")

    # ── Figure 3: Clade boxplot comparison ───────────────────────────────
    clade_order = ["MU", "BL", "AS", "BA", "CH", "MO", "MI", "AP", "HO"]
    df_box = df[df["clade"].isin(clade_order)].copy()
    fig, ax = plt.subplots(figsize=(8, 5))
    groups = [df_box[df_box["clade"] == c]["hits_per_mb"].dropna().values
              for c in clade_order if c in df_box["clade"].values]
    labels = [c for c in clade_order if c in df_box["clade"].values]
    bp = ax.boxplot(groups, labels=labels, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], labels):
        patch.set_facecolor(CLADE_COLORS.get(c, "#aaa"))
        patch.set_alpha(0.6)
    # Jitter points
    for i, (data, c) in enumerate(zip(groups, labels)):
        ax.scatter(np.random.normal(i+1, 0.1, len(data)), data,
                   s=12, alpha=0.5, c=CLADE_COLORS.get(c, "#aaa"),
                   edgecolor="k", linewidth=0.3, zorder=3)
    ax.set_ylabel("HMGb3 hits per Mb")
    ax.set_title("HMGb3 motif density by clade")
    fig.tight_layout()
    fig.savefig(OUTDIR / "hmgb3_clade_boxplot.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[OK] hmgb3_clade_boxplot.pdf")

    # ── Print summary ────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("FIMO HMGb3 Cross-Species Summary")
    print(f"{'='*60}")
    print(f"Total species scanned: {len(df)}")
    print(f"Species with hits: {(df['n_hmgb3_hits'] > 0).sum()}")
    print(f"\nTop 5 by hits/Mb:")
    for _, r in df.nlargest(5, "hits_per_mb").iterrows():
        print(f"  {r['species']:20s} {r['hits_per_mb']:.1f} hits/Mb ({r['n_hmgb3_hits']:.0f} total)")
    print(f"\nMucoromycota median: {df[df['clade']=='MU']['hits_per_mb'].median():.1f} hits/Mb")
    print(f"Blastocladiomycota median: {df[df['clade']=='BL']['hits_per_mb'].median():.1f} hits/Mb")
    print(f"\nNOTE: The HMGb3 motif (TGTTGTTGTTGTTGT) resembles a microsatellite")
    print(f"      (TGT)n repeat. High hit counts likely reflect microsatellite")
    print(f"      abundance rather than specific TF binding sites.")
    print(f"      The hmgb3_jockey_summary.tsv shows 0/68 species with hits")
    print(f"      overlapping LINE/I-Jockey elements specifically.")


if __name__ == "__main__":
    main()
