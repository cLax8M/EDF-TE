#!/usr/bin/env python3
"""
09_patristic_correction.py
──────────────────────────
Correct within-clade synteny ratios by patristic distance from a guide tree.

Method:
  For each clade, we perform a LOESS/linear regression of raw_score vs
  patristic distance d_ij between species pairs.  The corrected residual
  captures synteny gains or losses beyond what is expected at that divergence.

  Corrected synteny metric per pair:
      1. raw_score  = (c_ij - s_ij) / N_ij        (from step 05)
      2. d_ij       = sum of branch lengths between sp_A and sp_B in tree
      3. r_corrected = raw_score - f(d_ij)          (residual from regression)
         where f() is a LOWESS or polynomial fit of raw_score ~ d_ij
      4. r_norm_corrected = normalized to [0, 1] within each clade

  Per-species summary (analogous to step 08):
      mean_r_corrected, median_r_corrected over all within-clade partners
      (excluding fragmented assemblies)

Inputs:
  --ratios     all_clades_synteny_ratios.tsv (step 05)
  --tree       Guide tree in Newick format (branch lengths required)
  --out_dir    Output directory for corrected results

Outputs:
  {out_dir}/corrected_synteny_ratios.tsv          (all pairs, corrected)
  {out_dir}/perspecies_corrected_synteny.tsv       (per-species summary)
  {out_dir}/figures/                               (diagnostic plots per clade)

Usage:
  python 09_patristic_correction.py \\
      --ratios /data/.../all_clades_synteny_ratios.tsv \\
      --tree   /data/.../guide_rooted_Naegrub.treefile \\
      --out_dir /legserv/Temp/Carlos/EDF-TE/results/Synteny_guide_corrected
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from io import StringIO
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ── Clade definitions ──────────────────────────────────────────────────────
CLADES_ALL = ["AP", "AS", "BA", "BL", "BS", "CH", "CR", "EN",
              "GL", "HO", "KX", "MI", "MO", "MU", "NU", "OL", "ZO"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ratios", required=True,
                   help="all_clades_synteny_ratios.tsv (step 05 output)")
    p.add_argument("--tree", required=True,
                   help="Newick tree with branch lengths (guide tree)")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for corrected results")
    p.add_argument("--degree", type=int, default=2,
                   help="Polynomial degree for regression (default: 2)")
    p.add_argument("--min_pairs", type=int, default=4,
                   help="Minimum pairs per clade for regression (default: 4)")
    p.add_argument("--min_orthologs", type=int, default=2000,
                   help="Ortholog threshold for fragmentation detection")
    return p.parse_args()


def load_tree(tree_path):
    """Load Newick tree and return Bio.Phylo tree object."""
    tree = Phylo.read(tree_path, "newick")
    return tree


def get_tip_names(tree):
    """Return list of all tip (leaf) names in tree."""
    return [t.name for t in tree.get_terminals()]


def patristic_distance(tree, sp_a, sp_b):
    """Compute patristic distance between two tips in a tree."""
    try:
        d = tree.distance(sp_a, sp_b)
        return d
    except Exception:
        return np.nan


def detect_fragmented(df, min_orthologs=2000):
    """Species with r_ij_norm==0 in ALL pairs AND many orthologs."""
    fragmented = set()
    all_sp = set(df["species_A"]) | set(df["species_B"])
    for sp in all_sp:
        mask = (df["species_A"] == sp) | (df["species_B"] == sp)
        sub = df[mask]
        if sub.empty:
            continue
        if sub["r_ij_norm"].max() == 0 and sub["n_orthologs"].mean() > min_orthologs:
            fragmented.add(sp)
    return fragmented


def correct_clade(clade_df, tree, tip_names, degree=2, min_pairs=4):
    """
    Correct synteny ratios for one clade using patristic distances.

    Returns:
        corrected DataFrame with additional columns:
        d_ij, expected_raw, residual, r_corrected_norm
    """
    rows = clade_df.copy()

    # Map tree tip names: try exact match, then with '.aa' suffix
    def resolve_tip(sp):
        if sp in tip_names:
            return sp
        if sp + ".aa" in tip_names:
            return sp + ".aa"
        return None

    # Compute patristic distances
    distances = []
    for _, row in rows.iterrows():
        tip_a = resolve_tip(row["species_A"])
        tip_b = resolve_tip(row["species_B"])
        if tip_a and tip_b:
            d = patristic_distance(tree, tip_a, tip_b)
        else:
            d = np.nan
        distances.append(d)
    rows["d_ij"] = distances

    # Drop rows without valid patristic distances
    valid = rows.dropna(subset=["d_ij"])

    if len(valid) < min_pairs:
        # Not enough data for regression — return uncorrected
        rows["expected_raw"] = np.nan
        rows["residual"] = np.nan
        rows["r_corrected_norm"] = rows["r_ij_norm"]
        return rows

    # Polynomial regression: raw_score ~ poly(d_ij, degree)
    x = valid["d_ij"].values
    y = valid["raw_score"].values

    # Fit polynomial
    coeffs = np.polyfit(x, y, degree)
    poly = np.poly1d(coeffs)

    # Compute expected raw_score and residual for ALL rows
    rows["expected_raw"] = rows["d_ij"].apply(
        lambda d: poly(d) if not np.isnan(d) else np.nan)
    rows["residual"] = rows["raw_score"] - rows["expected_raw"]

    # Normalize residuals to [0, 1] within clade
    valid_residuals = rows["residual"].dropna()
    if len(valid_residuals) > 0:
        r_min = valid_residuals.min()
        r_max = valid_residuals.max()
        r_range = r_max - r_min
        if r_range > 0:
            rows["r_corrected_norm"] = rows["residual"].apply(
                lambda r: (r - r_min) / r_range if not np.isnan(r) else np.nan)
        else:
            rows["r_corrected_norm"] = 0.5
    else:
        rows["r_corrected_norm"] = rows["r_ij_norm"]

    return rows


def plot_clade_correction(clade, clade_df, fig_dir, degree=2):
    """Diagnostic scatter plot: raw_score vs d_ij with regression line."""
    valid = clade_df.dropna(subset=["d_ij", "raw_score"])
    if len(valid) < 3:
        return

    x = valid["d_ij"].values
    y = valid["raw_score"].values

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel 1: Raw score vs patristic distance
    ax = axes[0]
    ax.scatter(x, y, alpha=0.5, s=20, color="#4E79A7", edgecolors="none")
    # Regression line
    x_sort = np.sort(x)
    coeffs = np.polyfit(x, y, degree)
    poly = np.poly1d(coeffs)
    ax.plot(x_sort, poly(x_sort), "r-", linewidth=2, label=f"poly deg={degree}")
    rho, pval = spearmanr(x, y)
    ax.set_xlabel("Patristic distance (subs/site)", fontsize=10)
    ax.set_ylabel("Raw synteny score (c-s)/N", fontsize=10)
    ax.set_title(f"{clade}: Raw score vs distance\n"
                 f"ρ = {rho:.3f}, p = {pval:.2e}, n = {len(x)}",
                 fontsize=10)
    ax.legend(fontsize=8)

    # Panel 2: Residual vs d_ij
    ax = axes[1]
    residuals = valid["residual"].values if "residual" in valid.columns else y - poly(x)
    ax.scatter(x, residuals, alpha=0.5, s=20, color="#E15759", edgecolors="none")
    ax.axhline(0, color="grey", linewidth=0.7, linestyle="--")
    ax.set_xlabel("Patristic distance", fontsize=10)
    ax.set_ylabel("Residual (obs - expected)", fontsize=10)
    ax.set_title(f"{clade}: Distance-corrected residuals", fontsize=10)

    # Panel 3: Corrected normalized vs original normalized
    ax = axes[2]
    if "r_corrected_norm" in valid.columns:
        ax.scatter(valid["r_ij_norm"], valid["r_corrected_norm"],
                   alpha=0.5, s=20, color="#59A14F", edgecolors="none")
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.7)
        ax.set_xlabel("Original r_ij_norm", fontsize=10)
        ax.set_ylabel("Corrected r_corrected_norm", fontsize=10)
        ax.set_title(f"{clade}: Original vs corrected", fontsize=10)
    else:
        ax.text(0.5, 0.5, "No correction applied", ha="center", va="center",
                transform=ax.transAxes)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        path = fig_dir / f"{clade}_patristic_correction.{ext}"
        fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def compute_perspecies(df, min_orthologs=2000):
    """Compute per-species summary of corrected synteny ratios."""
    fragmented = detect_fragmented(df, min_orthologs)

    # Build long form
    rows_a = df[["clade", "species_A", "species_B", "n_orthologs",
                 "r_ij_norm", "r_corrected_norm"]].copy()
    rows_a.columns = ["clade", "species", "partner", "n_orthologs",
                      "r_orig", "r_corrected"]
    rows_b = df[["clade", "species_B", "species_A", "n_orthologs",
                 "r_ij_norm", "r_corrected_norm"]].copy()
    rows_b.columns = ["clade", "species", "partner", "n_orthologs",
                      "r_orig", "r_corrected"]
    long = pd.concat([rows_a, rows_b], ignore_index=True)

    # Exclude fragmented partners
    long = long[~long["partner"].isin(fragmented)]

    records = []
    all_species = sorted(set(df["species_A"]) | set(df["species_B"]))
    for sp in all_species:
        sub = long[long["species"] == sp]
        if sub.empty:
            clade = sp[:2]
            records.append({
                "species": sp, "clade": clade,
                "mean_r_orig": np.nan, "mean_r_corrected": np.nan,
                "median_r_corrected": np.nan, "max_r_corrected": np.nan,
                "n_partners": 0, "fragmented": sp in fragmented,
            })
            continue

        clade = sub["clade"].iloc[0]
        is_frag = sp in fragmented
        records.append({
            "species": sp,
            "clade": clade,
            "mean_r_orig": sub["r_orig"].mean(),
            "mean_r_corrected": sub["r_corrected"].mean(),
            "median_r_corrected": sub["r_corrected"].median(),
            "max_r_corrected": sub["r_corrected"].max(),
            "n_partners": len(sub),
            "fragmented": is_frag,
        })

    return pd.DataFrame(records)


def plot_overview(persp_df, fig_dir):
    """Overview figure: original vs corrected per-species mean synteny."""
    valid = persp_df.dropna(subset=["mean_r_orig", "mean_r_corrected"])
    if valid.empty:
        return

    import sys as _sys, os as _os
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                       "..", "..", "Section2_TElandscape", "scripts"))
    from clade_palette import CLADE_COLORS  # noqa: E402

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: Scatter original vs corrected
    ax = axes[0]
    colors = [CLADE_COLORS.get(c, "#AAAAAA") for c in valid["clade"]]
    ax.scatter(valid["mean_r_orig"], valid["mean_r_corrected"],
               c=colors, s=25, alpha=0.7, edgecolors="none")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("Mean synteny (original)", fontsize=10)
    ax.set_ylabel("Mean synteny (distance-corrected)", fontsize=10)
    ax.set_title("Per-species synteny: original vs patristic-corrected", fontsize=11)

    clades_present = sorted(valid["clade"].unique())
    patches = [mpatches.Patch(color=CLADE_COLORS.get(c, "#AAAAAA"), label=c)
               for c in clades_present]
    ax.legend(handles=patches, fontsize=6, ncol=3, loc="upper left",
              title="Clade", title_fontsize=7)

    # Panel 2: Ranked corrected values
    ax = axes[1]
    ranked = valid.sort_values("mean_r_corrected", ascending=True)
    bar_colors = [CLADE_COLORS.get(c, "#AAAAAA") for c in ranked["clade"]]
    ax.barh(range(len(ranked)), ranked["mean_r_corrected"],
            color=bar_colors, edgecolor="none", height=0.8)
    ax.set_yticks([])
    ax.set_xlabel("Mean corrected synteny ratio", fontsize=10)
    ax.set_title("Species ranked by distance-corrected synteny", fontsize=11)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"overview_corrected_synteny.{ext}",
                    dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    out = Path(args.out_dir)
    fig_dir = out / "figures"
    out.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading synteny ratios: {args.ratios}")
    df = pd.read_csv(args.ratios, sep="\t")
    print(f"  {len(df)} pairs across {df['clade'].nunique()} clades")

    # Load tree
    print(f"Loading guide tree: {args.tree}")
    tree = load_tree(args.tree)
    tip_names = set(get_tip_names(tree))
    print(f"  {len(tip_names)} tips in tree")

    # Correct each clade
    corrected_frames = []
    for clade in sorted(df["clade"].unique()):
        clade_df = df[df["clade"] == clade].copy()
        n_pairs = len(clade_df)
        print(f"  [{clade}] {n_pairs} pairs", end="")

        corrected = correct_clade(clade_df, tree, tip_names,
                                  degree=args.degree, min_pairs=args.min_pairs)
        corrected_frames.append(corrected)

        # Diagnostic figure
        valid_d = corrected.dropna(subset=["d_ij"])
        if len(valid_d) >= 3:
            rho, _ = spearmanr(valid_d["d_ij"], valid_d["raw_score"])
            print(f" → ρ(d, raw) = {rho:.3f} ({len(valid_d)} with distances)")
            plot_clade_correction(clade, corrected, fig_dir, degree=args.degree)
        else:
            print(f" → too few distances ({len(valid_d)})")

    # Combine all clades
    all_corrected = pd.concat(corrected_frames, ignore_index=True)

    # Save corrected ratios
    out_ratios = out / "corrected_synteny_ratios.tsv"
    all_corrected.to_csv(out_ratios, sep="\t", index=False, float_format="%.6f")
    print(f"\nSaved corrected ratios: {out_ratios}")

    # Per-species summary
    persp = compute_perspecies(all_corrected, min_orthologs=args.min_orthologs)
    out_persp = out / "perspecies_corrected_synteny.tsv"
    persp.to_csv(out_persp, sep="\t", index=False, float_format="%.6f")
    print(f"Saved per-species summary: {out_persp}")

    # Overview figure
    plot_overview(persp, fig_dir)
    print(f"Figures saved to: {fig_dir}")

    # Print correlation summary
    print("\n=== Patristic correction summary ===")
    for clade in sorted(all_corrected["clade"].unique()):
        sub = all_corrected[all_corrected["clade"] == clade]
        valid = sub.dropna(subset=["d_ij"])
        if len(valid) >= 3:
            rho_raw, _ = spearmanr(valid["d_ij"], valid["raw_score"])
            rho_corr, _ = spearmanr(valid["d_ij"], valid["residual"])
            print(f"  {clade:3s}: ρ(d, raw)={rho_raw:+.3f} → ρ(d, residual)={rho_corr:+.3f} "
                  f"({len(valid)} pairs)")
        else:
            print(f"  {clade:3s}: insufficient data ({len(valid)} pairs)")


if __name__ == "__main__":
    main()
