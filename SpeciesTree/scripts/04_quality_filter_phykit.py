#!/usr/bin/env python3
# =============================================================================
# 04_quality_filter_phykit.py — Alignment quality filtering with PhyKIT
#
# PURPOSE : Computes four alignment quality metrics for each trimmed HOG
#           alignment and applies thresholds to retain only high-quality
#           markers for species tree inference.
#
# METRICS AND THRESHOLDS:
#
#   NOTE: This dataset includes 238 diverse eukaryotes (fungi + outgroups).
#   RCV values are therefore higher than in fungi-only studies (Szánthó 2025
#   used 110 fungi with MAX_RCV=0.1; our median RCV ≈ 0.25). Saturation
#   requires a gene tree (computed in step 05), so it is NOT applied here.
#
#   1. Parsimony-Informative Sites (PIS)
#      - Threshold: ≥ MIN_PIS (default 25 sites)
#      - Genes with too few PIS contribute noise, not signal
#
#   2. RCV (Relative Composition Variability) — PhyKIT v2 command: 'rcv'
#      - Threshold: RCV ≤ MAX_RCV (default 0.5, relaxed for diverse eukaryotes)
#      - Lower RCV = more compositionally homogeneous alignment = better signal
#      - Note: PhyKIT v2 renamed 'rcfv' → 'rcv'
#
#   3. Occupancy (re-computed after trimming)
#      - Threshold: ≥ MIN_SPECIES_FRACTION (0.40 = 95/238)
#      - Ensures trimming didn't drastically reduce taxon sampling
#
#   Saturation filter → deferred to step 06 (after gene trees are built)
#
# OUTPUT  :
#   04_filtered/alignments/<HOG>.trim.fa    — passing alignments (symlinks)
#   04_filtered/gene_list.txt               — list of passing HOG IDs
#   04_filtered/quality_metrics.tsv         — all metrics for all HOGs
#   04_filtered/summary.txt                 — filtering summary
#   04_filtered/rejected/                   — rejected HOG details
#
# USAGE   : python 04_quality_filter_phykit.py [--trimmed-dir DIR] [--out-dir DIR]
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# REFS    : Szánthó et al. 2025, Merényi et al. 2023
# =============================================================================

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord

# ─── DEFAULT CONFIGURATION ────────────────────────────────────────────────────
WORK_DIR       = "/data/carlos/EDF-TE/SpeciesTree"
TRIMMED_DIR    = f"{WORK_DIR}/03_trimmed"
OUT_DIR        = f"{WORK_DIR}/04_filtered"
MIN_PIS        = 25       # minimum parsimony-informative sites
MAX_RCV        = 0.5      # max relative composition variability (relaxed for 238 diverse eukaryotes)
MIN_OCCUPANCY  = 0.40     # minimum fraction of 238 species (= 95 species)
N_TOTAL_SPECIES = 238


def parse_args():
    p = argparse.ArgumentParser(
        description="Quality-filter HOG alignments using PhyKIT metrics")
    p.add_argument("--trimmed-dir",       default=TRIMMED_DIR)
    p.add_argument("--out-dir",           default=OUT_DIR)
    p.add_argument("--min-pis",           type=int,   default=MIN_PIS)
    p.add_argument("--max-rcv",           type=float, default=MAX_RCV)
    p.add_argument("--min-occupancy",     type=float, default=MIN_OCCUPANCY)
    p.add_argument("--n-species",         type=int,   default=N_TOTAL_SPECIES)
    p.add_argument("--threads",           type=int,   default=8)
    return p.parse_args()


def run_phykit(cmd: list, capture=True) -> tuple:
    """Run a PhyKIT command and return (stdout, returncode)."""
    result = subprocess.run(
        cmd, capture_output=capture, text=True)
    return result.stdout.strip(), result.returncode


def get_pis(aln_path: str) -> int:
    """
    Calculate parsimony-informative sites using PhyKIT.
    Returns number of PIS or -1 on error.
    """
    stdout, rc = run_phykit(["phykit", "pis", aln_path])
    if rc != 0:
        return -1
    try:
        # PhyKIT pis output: "X\tY" where X = PIS count, Y = fraction
        return int(stdout.split()[0])
    except (ValueError, IndexError):
        return -1


def get_rcv(aln_path: str) -> float:
    """
    Calculate Relative Composition Variability using PhyKIT v2 ('rcv' command).
    In PhyKIT v2, the old 'rcfv' command was renamed to 'rcv'.
    Returns RCV value or -1.0 on error.
    """
    stdout, rc = run_phykit(["phykit", "rcv", aln_path])
    if rc != 0:
        return -1.0
    try:
        # PhyKIT rcv output: single float value
        return float(stdout.strip().split()[0])
    except (ValueError, IndexError):
        return -1.0


def get_occupancy(aln_path: str, n_total: int) -> float:
    """
    Compute occupancy = fraction of total species represented in alignment.
    (Re-computed after trimming, in case trimAl removed sequences.)
    """
    n_seqs = sum(1 for _ in SeqIO.parse(aln_path, "fasta"))
    return n_seqs / n_total


def compute_all_metrics(aln_path: str, n_total: int) -> dict:
    """Compute all four quality metrics for a single alignment."""
    hog_id = Path(aln_path).stem.replace(".trim", "")

    metrics = {
        "og_id"     : hog_id,
        "aln_path"  : aln_path,
        "pis"       : get_pis(aln_path),
        "rcv"       : get_rcv(aln_path),
        "occupancy" : get_occupancy(aln_path, n_total),
    }
    return metrics


def apply_filters(df: pd.DataFrame, args) -> pd.DataFrame:
    """Add pass/fail columns and a final 'pass_all' column."""
    df["pass_pis"]      = df["pis"]       >= args.min_pis
    df["pass_rcv"]      = (df["rcv"]  >= 0) & (df["rcv"]  <= args.max_rcv)
    df["pass_occupancy"]= df["occupancy"] >= args.min_occupancy
    df["pass_all"]      = df["pass_pis"] & df["pass_rcv"] & df["pass_occupancy"]

    # Flag PhyKIT errors (-1 values)
    df["phykit_error"] = (df["pis"] == -1) | (df["rcv"] == -1.0)

    return df


def main():
    args = parse_args()

    print("=" * 70)
    print(" STEP 04: Quality filtering with PhyKIT")
    print("=" * 70)
    print(f" Trimmed dir   : {args.trimmed_dir}")
    print(f" Output dir    : {args.out_dir}")
    print(f" Threads       : {args.threads}")
    print(f" Thresholds:")
    print(f"   Min PIS      : ≥ {args.min_pis}")
    print(f"   Max RCV      : ≤ {args.max_rcv}  (relaxed for 238 diverse eukaryotes)")
    print(f"   Min occupancy: ≥ {args.min_occupancy:.0%} ({int(args.min_occupancy * args.n_species)} of {args.n_species} spp)")
    print(f"   Saturation   : deferred to step 06 (requires gene tree)\n")

    # Find input alignments
    trimmed_files = sorted(Path(args.trimmed_dir).glob("*.trim.fa"))
    print(f"[04] Found {len(trimmed_files)} trimmed alignments")

    if not trimmed_files:
        print("ERROR: No .trim.fa files found in trimmed directory.")
        sys.exit(1)

    # Create output directories
    aln_out_dir  = os.path.join(args.out_dir, "alignments")
    reject_dir   = os.path.join(args.out_dir, "rejected")
    os.makedirs(aln_out_dir,  exist_ok=True)
    os.makedirs(reject_dir,   exist_ok=True)

    # Compute metrics for all alignments
    print(f"\n[04] Computing quality metrics (this may take ~{len(trimmed_files)//60 + 1} min)...")
    all_metrics = []

    for i, aln_path in enumerate(trimmed_files):
        metrics = compute_all_metrics(str(aln_path), args.n_species)
        all_metrics.append(metrics)

        if (i + 1) % 50 == 0:
            print(f"     Processed {i+1}/{len(trimmed_files)} alignments...")

    # Build DataFrame and apply filters
    df = pd.DataFrame(all_metrics)
    df = apply_filters(df, args)

    # Write metrics table
    metrics_path = os.path.join(args.out_dir, "quality_metrics.tsv")
    df.to_csv(metrics_path, sep="\t", index=False, float_format="%.4f")

    # Copy/link passing alignments
    passing = df[df["pass_all"]]
    failing = df[~df["pass_all"]]

    print(f"\n[04] Copying passing alignments...")
    for _, row in passing.iterrows():
        src = row["aln_path"]
        dst = os.path.join(aln_out_dir, os.path.basename(src))
        shutil.copy2(src, dst)

    # Write gene list
    gene_list_path = os.path.join(args.out_dir, "gene_list.txt")
    passing["hog_id"].to_csv(gene_list_path, index=False, header=False)

    # Save rejected gene details
    reject_path = os.path.join(reject_dir, "rejected_genes.tsv")
    failing.to_csv(reject_path, sep="\t", index=False, float_format="%.4f")

    # Failure reason breakdown
    print(f"\n{'=' * 70}")
    print(f" STEP 04 RESULTS — Quality Filter Summary")
    print(f"{'=' * 70}")
    print(f"  Total alignments evaluated : {len(df)}")
    print(f"  Passed all filters         : {len(passing)}")
    print(f"  Rejected                   : {len(failing)}")
    print(f"\n  Rejection breakdown:")
    print(f"    Failed PIS (<{args.min_pis})            : {(~df['pass_pis']).sum()}")
    print(f"    Failed RCV (>{args.max_rcv})            : {(~df['pass_rcv']).sum()}")
    print(f"    Failed occupancy (<{args.min_occupancy:.0%})  : {(~df['pass_occupancy']).sum()}")
    print(f"    PhyKIT errors               : {df['phykit_error'].sum()}")
    print(f"\n  Output:")
    print(f"    Metrics table   : {metrics_path}")
    print(f"    Gene list       : {gene_list_path}")
    print(f"    Alignments dir  : {aln_out_dir}")
    print(f"\n  Next step: sbatch 05_gene_trees_array.sh")
    print(f"{'=' * 70}")

    # Summary file
    summary_path = os.path.join(args.out_dir, "summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("STEP 04 — PhyKIT Quality Filtering — Summary\n")
        fh.write("=" * 45 + "\n")
        fh.write(f"Alignments evaluated  : {len(df)}\n")
        fh.write(f"Passed all filters    : {len(passing)}\n")
        fh.write(f"Rejected              : {len(failing)}\n")
        fh.write(f"\nThresholds used:\n")
        fh.write(f"  Min PIS          >= {args.min_pis}\n")
        fh.write(f"  Max RCV          <= {args.max_rcv}  (relaxed for 238 diverse eukaryotes)\n")
        fh.write(f"  Min occupancy    >= {args.min_occupancy:.0%}\n")
        fh.write(f"  Saturation       → deferred to step 06 (requires gene tree)\n")
        fh.write(f"\nRef: Szánthó et al. 2025, Merényi et al. 2023\n")


if __name__ == "__main__":
    main()
