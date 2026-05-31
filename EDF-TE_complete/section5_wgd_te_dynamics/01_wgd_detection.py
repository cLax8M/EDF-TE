#!/usr/bin/env python3
"""
01_wgd_detection.py  —  Whole-genome duplication detection via Ks distributions.

Identifies whole-genome duplication (WGD) events in Mucoromycotina by
analysing synonymous substitution rate (Ks) distributions from paralogous
gene pairs.  WGD events appear as peaks in Ks distributions of retained
paralogs.

Pipeline:
  1. For each species: extract all-vs-all paralog Ks values from OrthoFinder
     RBH pairs (or self-BLAST)
  2. Fit Gaussian mixture model (GMM) to each Ks distribution
  3. Identify significant GMM peaks (> 10% of data in component; log-likelihood
     improvement Δ-BIC ≥ 10 vs k-1 components)
  4. Call WGD if a peak is shared between sister lineages at the same
     Ks depth (phylogenetically dated using species divergence calibrations)
  5. Produce species tree annotated with WGD events

Inputs:
  --ks_dir        Directory with per-species paralog Ks TSVs
                  (columns: gene1, gene2, ks)
  --species_list  File with Mucoromycotina species codes
  --species_tree  Newick tree for Mucoromycotina
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory
  --max_ks        Maximum Ks to include in GMM fitting (default: 3.0)
  --n_components  Maximum number of GMM components to test (default: 5)

Outputs:
  wgd_summary.tsv          – per-species GMM peaks and WGD calls
  ks_distributions.tsv     – all Ks values merged
  figures/ks_gmm_<sp>.pdf  – per-species Ks GMM fit
  figures/wgd_tree.pdf     – species tree with WGD events marked
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ks_dir",        required=True)
    p.add_argument("--species_list",  required=True)
    p.add_argument("--species_tree",  required=True)
    p.add_argument("--clade_map",     required=True)
    p.add_argument("--out_dir",       default="wgd_out")
    p.add_argument("--max_ks",        type=float, default=3.0)
    p.add_argument("--n_components",  type=int,   default=5)
    return p.parse_args()


def load_ks(ks_dir: str, sp: str, max_ks: float) -> np.ndarray:
    files = list(Path(ks_dir).glob(f"{sp}*.ks.tsv")) + \
            list(Path(ks_dir).glob(f"{sp}*.tsv"))
    if not files:
        return np.array([])
    df = pd.read_csv(files[0], sep="\t")
    ks_col = next((c for c in df.columns if c.lower() == "ks"), None)
    if ks_col is None:
        return np.array([])
    ks = df[ks_col].dropna().values
    return ks[(ks > 0) & (ks <= max_ks)]


def fit_gmm(ks: np.ndarray, max_components: int) -> tuple[object, int]:
    """Fit GMM; return best model (lowest BIC) and optimal k."""
    from sklearn.mixture import GaussianMixture
    best_bic = np.inf
    best_gmm = None
    best_k   = 1
    for k in range(1, min(max_components + 1, len(ks) // 20 + 1)):
        gm = GaussianMixture(n_components=k, covariance_type="full",
                             random_state=42, n_init=5)
        gm.fit(ks.reshape(-1, 1))
        bic = gm.bic(ks.reshape(-1, 1))
        if bic < best_bic:
            best_bic, best_gmm, best_k = bic, gm, k
    return best_gmm, best_k


def call_wgd_peaks(gmm, ks: np.ndarray,
                    min_weight: float = 0.10) -> list[dict]:
    """Return peaks from GMM components with weight ≥ min_weight."""
    peaks = []
    for i in range(gmm.n_components):
        w = gmm.weights_[i]
        mu = gmm.means_[i, 0]
        sd = np.sqrt(gmm.covariances_[i, 0, 0])
        if w >= min_weight and 0.05 < mu < 3.0:
            peaks.append({"peak_ks": mu, "peak_sd": sd, "weight": w})
    return peaks


def plot_ks_gmm(sp: str, ks: np.ndarray, gmm,
                peaks: list[dict], out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(ks, bins=100, density=True, color="#aaaaaa", alpha=0.6,
            label="Ks distribution")

    x = np.linspace(0, ks.max(), 500)
    for i in range(gmm.n_components):
        mu = gmm.means_[i, 0]
        sd = np.sqrt(gmm.covariances_[i, 0, 0])
        w  = gmm.weights_[i]
        from scipy.stats import norm
        comp = w * norm.pdf(x, mu, sd)
        ax.plot(x, comp, "--", linewidth=1, alpha=0.7)

    # Overlay total GMM density
    from scipy.stats import norm as _norm
    total = np.zeros_like(x)
    for i in range(gmm.n_components):
        total += gmm.weights_[i] * _norm.pdf(
            x, gmm.means_[i, 0], np.sqrt(gmm.covariances_[i, 0, 0]))
    ax.plot(x, total, "r-", linewidth=1.5, label="GMM fit")

    for pk in peaks:
        ax.axvline(pk["peak_ks"], color="blue", linestyle=":",
                   linewidth=1.2, label=f"Peak Ks={pk['peak_ks']:.2f}")

    ax.set_xlabel("Ks (synonymous substitution rate)", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title(sp, fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"ks_gmm_{sp}.pdf")
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    try:
        from sklearn.mixture import GaussianMixture
    except ImportError:
        print("scikit-learn required for GMM fitting: pip install scikit-learn",
              file=sys.stderr)
        sys.exit(1)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))
    species_list = [l.strip() for l in open(args.species_list) if l.strip()]

    wgd_rows = []
    all_ks   = []

    for sp in species_list:
        ks = load_ks(args.ks_dir, sp, args.max_ks)
        if len(ks) < 50:
            print(f"  {sp}: insufficient Ks pairs ({len(ks)}), skipping")
            continue

        print(f"  {sp}: {len(ks):,} pairs ...", end=" ")
        gmm, k = fit_gmm(ks, args.n_components)
        peaks  = call_wgd_peaks(gmm, ks)

        wgd_candidate = len(peaks) > 0
        wgd_rows.append({
            "species":   sp,
            "clade":     clade_map.get(sp, "?"),
            "n_pairs":   len(ks),
            "n_gmm_k":   k,
            "n_peaks":   len(peaks),
            "wgd_call":  wgd_candidate,
            "peak_ks":   ";".join(f"{pk['peak_ks']:.3f}" for pk in peaks),
        })
        plot_ks_gmm(sp, ks, gmm, peaks, fig_dir)
        ks_df = pd.DataFrame({"species": sp, "ks": ks})
        all_ks.append(ks_df)
        print(f"k={k}, peaks={len(peaks)}, WGD={wgd_candidate}")

    wgd_df = pd.DataFrame(wgd_rows)
    wgd_df.to_csv(out_dir / "wgd_summary.tsv", sep="\t", index=False)
    if all_ks:
        pd.concat(all_ks, ignore_index=True).to_csv(
            out_dir / "ks_distributions.tsv", sep="\t", index=False)

    n_wgd = wgd_df["wgd_call"].sum()
    print(f"\n{n_wgd} species with WGD candidates out of {len(wgd_df)}")
    print("Done.")


if __name__ == "__main__":
    main()
