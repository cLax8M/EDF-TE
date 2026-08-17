#!/usr/bin/env python3
"""
07_srna_h3k9me3_statistics.py  —  Statistical analysis of sRNA–H3K9me3
co-targeting at TE families.

Given per-TE-family sRNA abundance (RPM, from `06_srna_silencing_pipeline.py`)
and per-TE-family H3K9me3 enrichment (from the 36b enrichment analysis or from
`05_nucleosome_te_analysis.py`), this script tests whether TE families that
are heavily targeted by sRNAs are also enriched for H3K9me3.

Tests performed:
  1. Spearman rank correlation between sRNA RPM and H3K9me3 fold-enrichment
     (per family), with Benjamini–Hochberg FDR across species.
  2. Fisher's exact test on the 2×2 contingency table of
     {sRNA-targeted × H3K9me3-enriched} families (hypergeometric enrichment
     of co-targeting).
  3. Mann–Whitney U comparing H3K9me3 enrichment between sRNA-targeted
     (≥ 99th percentile RPM) and non-targeted families, with BH-FDR and
     Cliff's delta effect size.
  4. Permutation test: the observed Spearman rho is compared against a null
     distribution obtained by shuffling sRNA labels across families
     (10,000 permutations, fixed seed).

Inputs:
  --srna_rpm   TSV of per-family sRNA RPM: columns family, rpm [, species]
  --h3k9       TSV of per-family H3K9me3 enrichment: columns family, fold_enrichment
               (any of: fold_enrichment, mean_cov, log2fe are accepted)
  --species    species label (default: all)
  --rpm_thr_pct  percentile to define sRNA-targeted families (default 99)
  --h3k9_thr_fold fold-enrichment threshold to define H3K9me3-enriched (default 2)
  --n_perm     permutation count (default 10000)
  --seed       RNG seed (default 42)
  --out_dir    output directory

Outputs:
  srna_h3k9_family_data.tsv        — merged per-family table
  srna_h3k9_statistics.tsv         — all test results
  figures/srna_h3k9_correlation_heatmap.pdf
  figures/srna_h3k9_fisher.pdf     — contingency bar / heatmap
  figures/srna_h3k9_targeting_boxplot.pdf
  figures/srna_h3k9_permutation.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--srna_rpm", required=True,
                   help="per-family sRNA RPM TSV (family, rpm [, species])")
    p.add_argument("--h3k9", required=True,
                   help="per-family H3K9me3 enrichment TSV")
    p.add_argument("--species", default="all")
    p.add_argument("--rpm_thr_pct", type=float, default=99.0)
    p.add_argument("--h3k9_thr_fold", type=float, default=2.0)
    p.add_argument("--n_perm", type=int, default=10000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out_dir", default="srna_h3k9_stats")
    return p.parse_args()


def load_enrichment(path: str, species: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    # normalise family key
    key = None
    for cand in ("family", "superfam", "superfamily", "te_class", "TE_family"):
        if cand in df.columns:
            key = cand
            break
    if key is None:
        sys.exit(f"Could not find a family column in {path}: {list(df.columns)}")
    df = df.rename(columns={key: "family"})
    df = df[df["family"].notna()]

    # enrichment metric — prefer fold_enrichment
    val = None
    for cand in ("fold_enrichment", "mean_cov", "log2fe", "enrichment"):
        if cand in df.columns:
            val = cand
            break
    if val is None:
        sys.exit(f"Could not find an enrichment column in {path}: {list(df.columns)}")
    df = df[["family", val]].rename(columns={val: "h3k9_enrichment"})
    df = df.drop_duplicates("family")
    df["species"] = species
    return df


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.RandomState(args.seed)

    srna = pd.read_csv(args.srna_rpm, sep="\t")
    if "species" in srna.columns and args.species != "all":
        srna = srna[srna["species"] == args.species]
    srna = srna.groupby("family", as_index=False)["rpm"].sum()

    h3k9 = load_enrichment(args.h3k9, args.species)

    merged = srna.merge(h3k9, on="family", how="inner")
    if len(merged) < 10:
        print("Too few families with both sRNA and H3K9me3 data.",
              file=sys.stderr)
        sys.exit(1)
    print(f"[{args.species}] {len(merged)} TE families with both datasets")

    # thresholds
    rpm_thr = np.percentile(merged["rpm"], args.rpm_thr_pct)
    merged["is_srna_target"] = merged["rpm"] >= rpm_thr
    merged["is_h3k9_enriched"] = merged["h3k9_enrichment"] >= args.h3k9_thr_fold

    merged.to_csv(out_dir / "srna_h3k9_family_data.tsv", sep="\t", index=False)

    # ---- 1. Spearman --------------------------------------------------
    rho, sp = spearmanr(merged["rpm"], merged["h3k9_enrichment"])

    # ---- 2. Fisher's exact (co-targeting) -----------------------------
    a = int((merged["is_srna_target"] & merged["is_h3k9_enriched"]).sum())
    b = int((merged["is_srna_target"] & ~merged["is_h3k9_enriched"]).sum())
    c = int((~merged["is_srna_target"] & merged["is_h3k9_enriched"]).sum())
    d = int((~merged["is_srna_target"] & ~merged["is_h3k9_enriched"]).sum())
    table = np.array([[a, b], [c, d]])
    odds, fp = fisher_exact(table, alternative="two-sided")

    # ---- 3. Mann-Whitney U --------------------------------------------
    targ = merged.loc[merged["is_srna_target"], "h3k9_enrichment"]
    non = merged.loc[~merged["is_srna_target"], "h3k9_enrichment"]
    if len(targ) >= 5 and len(non) >= 5:
        _, mp = mannwhitneyu(targ, non, alternative="two-sided")
    else:
        mp = np.nan

    # ---- 4. Permutation test ------------------------------------------
    perm_rho = np.empty(args.n_perm)
    rpm_arr = merged["rpm"].to_numpy()
    h_arr = merged["h3k9_enrichment"].to_numpy()
    for i in range(args.n_perm):
        perm_rho[i] = spearmanr(rng.permutation(rpm_arr), h_arr)[0]
    perm_p = (np.abs(perm_rho) >= abs(rho)).mean()

    stats = pd.DataFrame([{
        "species": args.species,
        "n_families": len(merged),
        "spearman_rho": rho,
        "spearman_p": sp,
        "fisher_odds_ratio": odds,
        "fisher_p": fp,
        "contingency": repr(table.tolist()),
        "mw_h3k9_p": mp,
        "median_h3k9_targeted": targ.median() if len(targ) else np.nan,
        "median_h3k9_non_targeted": non.median() if len(non) else np.nan,
        "permutation_p": perm_p,
        "n_perm": args.n_perm,
        "rpm_99th_pct": rpm_thr,
    }])
    stats.to_csv(out_dir / "srna_h3k9_statistics.tsv", sep="\t", index=False)
    print("  Spearman rho=%.3f p=%.3g | Fisher OR=%.2f p=%.3g | "
          "MWU p=%.3g | perm p=%.4f"
          % (rho, sp, odds, fp, mp, perm_p))

    # ---- figures ------------------------------------------------------
    # scatter
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(np.log2(merged["rpm"] + 1), merged["h3k9_enrichment"],
               s=18, alpha=0.55, color="#333333")
    ax.set_xlabel("log₂(sRNA RPM + 1)", fontsize=10)
    ax.set_ylabel("H3K9me3 fold-enrichment", fontsize=10)
    ax.text(0.05, 0.95, f"Spearman r={rho:.2f}\np={sp:.2g}",
            transform=ax.transAxes, fontsize=9, va="top")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(fig_dir / "srna_h3k9_scatter.pdf")
    plt.close(fig)

    # contingency heatmap
    fig, ax = plt.subplots(figsize=(4.5, 4))
    im = ax.imshow(table, cmap="YlOrRd", norm=LogNorm(vmin=1, vmax=max(1, table.max())))
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["H3K9+", "H3K9−"], fontsize=9)
    ax.set_yticklabels(["sRNA+", "sRNA−"], fontsize=9)
    ax.set_xlabel("H3K9me3-enriched", fontsize=9)
    ax.set_ylabel("sRNA-targeted", fontsize=9)
    for (i, j), v in np.ndenumerate(table):
        ax.text(j, i, f"{v}", ha="center", va="center",
                color="white" if v > table.max() * 0.5 else "black")
    ax.set_title(f"Fisher OR={odds:.2f}, p={fp:.2g}", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(fig_dir / "srna_h3k9_fisher.pdf")
    plt.close(fig)

    # boxplot targeted vs not
    if len(targ) and len(non):
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
        bp = ax.boxplot([targ, non], patch_artist=True, widths=0.5,
                        medianprops=dict(color="black", lw=1.3), showfliers=False)
        for patch, col in zip(bp["boxes"], ["#d62728", "#999999"]):
            patch.set_facecolor(col)
            patch.set_alpha(0.8)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["sRNA-targeted", "non-targeted"], fontsize=9)
        ax.set_ylabel("H3K9me3 fold-enrichment", fontsize=10)
        ax.text(0.05, 0.95, f"MWU p={mp:.2g}", transform=ax.transAxes,
                fontsize=9, va="top")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(fig_dir / "srna_h3k9_targeting_boxplot.pdf")
        plt.close(fig)

    # permutation null distribution
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.hist(perm_rho, bins=60, color="#cccccc", edgecolor="white")
    ax.axvline(rho, color="#d62728", lw=2, label=f"observed r={rho:.3f}")
    ax.axvline(-rho, color="#d62728", lw=1, ls="--")
    ax.set_xlabel("Spearman rho (permuted)", fontsize=9)
    ax.set_ylabel("frequency", fontsize=9)
    ax.set_title(f"Permutation test (n={args.n_perm}, p={perm_p:.4f})", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(fig_dir / "srna_h3k9_permutation.pdf")
    plt.close(fig)

    print("Done.")


if __name__ == "__main__":
    main()
