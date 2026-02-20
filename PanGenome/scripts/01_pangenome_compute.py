#!/usr/bin/env python3
"""
01_pangenome_compute.py
-----------------------
Pangenome analysis from OrthoFinder flat orthogroups (Orthogroups.GeneCount.tsv).

For each taxonomic group (and for the combined dataset), computes:
  - Pangenome accumulation   : gene families present in >= 1 genome
  - Core genome decay        : gene families present in ALL n genomes
  - Shell genome             : gene families present in >= 2 but < n genomes
  - Unique genes             : gene families present in exactly 1 genome

Statistics are estimated via 100 random subsamples at each genome count n,
reporting mean ± standard deviation (matching the approach in Steindorff et al. 2024).

Usage:
    python 01_pangenome_compute.py [--n_iter 100] [--seed 42]

Outputs (to /data/carlos/EDF-TE/PanGenome/results/):
    per_group/{GROUP}_pangenome.tsv
    combined/combined_pangenome.tsv
    summary_final_counts.tsv
"""

import argparse
import os
import sys
import time
import numpy as np
import pandas as pd
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR       = Path("/data/carlos/EDF-TE/PanGenome/data")
RESULTS_DIR    = Path("/data/carlos/EDF-TE/PanGenome/results")
GENECOUNT_TSV  = DATA_DIR / "Orthogroups.GeneCount.tsv"

# ── Target groups (2-letter prefixes) ─────────────────────────────────────────
GROUP_NAMES = {
    "BL": "Blastocladiomycota",
    "CH": "Chytridiomycota",
    "EN": "Entomophthoromycota",
    "GL": "Glomeromycota",
    "KX": "Kickxellomycota",
    "MI": "Microsporidia",
    "MO": "Mortierellomycota",
    "MU": "Mucoromycota",
    "ZO": "Zoopagomycota",
    "CR": "Cryptomycota",
    "BS": "Basidiobolomycota",
}

# ── Argument parsing ───────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Pangenome accumulation analysis")
    p.add_argument("--n_iter", type=int, default=100,
                   help="Number of random subsamples per genome count (default: 100)")
    p.add_argument("--seed", type=int, default=42,
                   help="Base random seed for reproducibility (default: 42)")
    return p.parse_args()


# ── Core computation ───────────────────────────────────────────────────────────
def compute_stats(presence_matrix: np.ndarray, n_iter: int, rng: np.random.Generator) -> pd.DataFrame:
    """
    Given a binary presence/absence matrix (HOGs x species),
    compute pangenome/core/shell/unique accumulation curves via random subsampling.

    Parameters
    ----------
    presence_matrix : np.ndarray, shape (n_hogs, n_species), dtype bool/int
    n_iter          : int, number of random subsamples per n
    rng             : numpy random Generator

    Returns
    -------
    pd.DataFrame with columns:
        n, pangenome_mean, pangenome_sd,
           core_mean,      core_sd,
           shell_mean,     shell_sd,
           unique_mean,    unique_sd
    """
    n_hogs, n_species = presence_matrix.shape
    records = []

    for n in range(1, n_species + 1):
        pan_vals    = np.empty(n_iter, dtype=np.float64)
        core_vals   = np.empty(n_iter, dtype=np.float64)
        shell_vals  = np.empty(n_iter, dtype=np.float64)
        unique_vals = np.empty(n_iter, dtype=np.float64)

        for i in range(n_iter):
            # random subsample of n species (column indices, without replacement)
            idx = rng.choice(n_species, size=n, replace=False)
            sub = presence_matrix[:, idx]  # shape (n_hogs, n)

            row_sums = sub.sum(axis=1)     # number of sampled species having each HOG

            pan_vals[i]    = np.sum(row_sums >= 1)
            core_vals[i]   = np.sum(row_sums == n)
            unique_vals[i] = np.sum(row_sums == 1)
            # shell: present in >=2 AND <n; for n==1, shell is zero by definition
            if n >= 3:
                shell_vals[i] = np.sum((row_sums >= 2) & (row_sums < n))
            elif n == 2:
                # present in exactly 2 but < 2 is empty; shell = present in 2 = core here
                # standard definition: shell = >=2 and <n → for n=2 this is empty
                shell_vals[i] = 0.0
            else:
                shell_vals[i] = 0.0

        records.append({
            "n":              n,
            "pangenome_mean": pan_vals.mean(),
            "pangenome_sd":   pan_vals.std(ddof=0),
            "core_mean":      core_vals.mean(),
            "core_sd":        core_vals.std(ddof=0),
            "shell_mean":     shell_vals.mean(),
            "shell_sd":       shell_vals.std(ddof=0),
            "unique_mean":    unique_vals.mean(),
            "unique_sd":      unique_vals.std(ddof=0),
        })

    return pd.DataFrame(records)


def sanity_check(df: pd.DataFrame, label: str):
    """Basic sanity checks on the accumulation DataFrame."""
    errors = []
    # At n=1: core should equal pangenome
    row1 = df[df["n"] == 1].iloc[0]
    if abs(row1["core_mean"] - row1["pangenome_mean"]) > 1e-3:
        errors.append(f"n=1: core_mean ({row1['core_mean']:.1f}) != pangenome_mean ({row1['pangenome_mean']:.1f})")

    # Pangenome should be non-decreasing
    pan = df["pangenome_mean"].values
    if not np.all(np.diff(pan) >= -1e-6):
        errors.append("Pangenome is not monotonically non-decreasing")

    # Core should be non-increasing
    core = df["core_mean"].values
    if not np.all(np.diff(core) <= 1e-6):
        errors.append("Core is not monotonically non-increasing")

    if errors:
        print(f"  [WARNING] {label}: " + "; ".join(errors))
    else:
        print(f"  [OK] {label}: sanity checks passed")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    rng  = np.random.default_rng(args.seed)

    print(f"[{time.strftime('%H:%M:%S')}] Loading Orthogroups.GeneCount.tsv ...")
    # GeneCount TSV: rows = orthogroups, cols = species (integer gene counts) + 'Total'
    df_hog = pd.read_csv(GENECOUNT_TSV, sep="\t", index_col=0)
    # Drop the 'Total' column if present
    if "Total" in df_hog.columns:
        df_hog = df_hog.drop(columns=["Total"])
    print(f"  Loaded {len(df_hog):,} orthogroups x {df_hog.shape[1]} species")

    species_cols = list(df_hog.columns)

    # Build prefix → species list mapping
    prefix_to_species = {}
    for col in species_cols:
        prefix = col[:2]
        prefix_to_species.setdefault(prefix, []).append(col)

    # Report species counts
    print("\n  Species per group (target groups):")
    for code, name in GROUP_NAMES.items():
        sp_list = prefix_to_species.get(code, [])
        print(f"    {code} ({name}): {len(sp_list)} species")

    # Build binary presence/absence matrix (all species at once)
    print(f"\n[{time.strftime('%H:%M:%S')}] Building binary presence/absence matrix ...")
    # Presence = gene count > 0
    pa_full = (df_hog.values > 0)  # shape (n_orthogroups, n_species), dtype bool

    # Filter to target-group species only
    target_species = []
    for code in GROUP_NAMES:
        target_species.extend(prefix_to_species.get(code, []))

    target_idx = [species_cols.index(sp) for sp in target_species]
    pa_target  = pa_full[:, target_idx]             # (n_hogs, n_target_species)
    # Drop HOGs with zero presence across all target species (irrelevant rows)
    hog_present = pa_target.sum(axis=1) > 0
    pa_target   = pa_target[hog_present]
    print(f"  {pa_target.shape[0]:,} HOGs present in at least one target species; "
          f"{pa_target.shape[1]} target species")

    # ── Per-group analysis ─────────────────────────────────────────────────────
    summary_rows = []
    out_per_group = RESULTS_DIR / "per_group"
    out_per_group.mkdir(parents=True, exist_ok=True)

    for code, name in GROUP_NAMES.items():
        sp_list = prefix_to_species.get(code, [])
        if not sp_list:
            print(f"  [WARNING] No species found for group {code}, skipping")
            continue

        # Column indices within target_species list
        grp_idx  = [target_species.index(sp) for sp in sp_list]
        pa_grp   = pa_target[:, grp_idx]
        # Keep only HOGs present in at least one species of this group
        hog_mask = pa_grp.sum(axis=1) > 0
        pa_grp   = pa_grp[hog_mask]

        n_sp = len(sp_list)
        print(f"\n[{time.strftime('%H:%M:%S')}] Group {code} ({name}): "
              f"{n_sp} species, {pa_grp.shape[0]:,} group-present HOGs")

        df_stats = compute_stats(pa_grp, n_iter=args.n_iter, rng=rng)
        sanity_check(df_stats, f"{code} ({name})")

        # Save per-group TSV
        out_path = out_per_group / f"{code}_pangenome.tsv"
        df_stats.to_csv(out_path, sep="\t", index=False, float_format="%.2f")
        print(f"  Saved → {out_path}")

        # Summary at n = N_species
        last = df_stats.iloc[-1]
        summary_rows.append({
            "group_code":    code,
            "group_name":    name,
            "n_species":     n_sp,
            "n_hogs_total":  pa_grp.shape[0],
            "pangenome":     round(last["pangenome_mean"]),
            "core":          round(last["core_mean"]),
            "shell":         round(last["shell_mean"]),
            "unique":        round(last["unique_mean"]),
        })

    # ── Combined analysis ──────────────────────────────────────────────────────
    print(f"\n[{time.strftime('%H:%M:%S')}] Combined analysis ({pa_target.shape[1]} species) ...")
    out_combined = RESULTS_DIR / "combined"
    out_combined.mkdir(parents=True, exist_ok=True)

    df_combined = compute_stats(pa_target, n_iter=args.n_iter, rng=rng)
    sanity_check(df_combined, "Combined")

    out_path = out_combined / "combined_pangenome.tsv"
    df_combined.to_csv(out_path, sep="\t", index=False, float_format="%.2f")
    print(f"  Saved → {out_path}")

    last_c = df_combined.iloc[-1]
    summary_rows.append({
        "group_code":    "ALL",
        "group_name":    "Combined (all 11 groups)",
        "n_species":     pa_target.shape[1],
        "n_hogs_total":  pa_target.shape[0],
        "pangenome":     round(last_c["pangenome_mean"]),
        "core":          round(last_c["core_mean"]),
        "shell":         round(last_c["shell_mean"]),
        "unique":        round(last_c["unique_mean"]),
    })

    # ── Summary table ──────────────────────────────────────────────────────────
    df_summary = pd.DataFrame(summary_rows)
    summary_path = RESULTS_DIR / "summary_final_counts.tsv"
    df_summary.to_csv(summary_path, sep="\t", index=False)
    print(f"\n[{time.strftime('%H:%M:%S')}] Summary saved → {summary_path}")
    print("\n" + df_summary.to_string(index=False))
    print(f"\n[{time.strftime('%H:%M:%S')}] Done.")


if __name__ == "__main__":
    main()
