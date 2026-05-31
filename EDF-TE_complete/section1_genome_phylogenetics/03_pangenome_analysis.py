#!/usr/bin/env python3
"""
03_pangenome_analysis.py  —  Pangenome accumulation curves across fungal clades.

Estimates core, shell, and unique gene-family counts as a function of genome
number for each taxonomic clade and for the combined dataset.  Statistics are
derived from 100 random subsamples at each n (Heap's law approach).

Input:
  --genecount  OrthoFinder Orthogroups.GeneCount.tsv
  --clade_map  two-column TSV: species_code<TAB>clade_tag
  --out_dir    output directory for TSVs and figures

Outputs:
  pangenome/<CLADE>_pangenome.tsv  – mean ± SD accumulation curves per clade
  pangenome/combined_pangenome.tsv – combined across all clades
  pangenome/summary_counts.tsv     – final pan/core/shell/unique at max n
  figures/pangenome_curves.pdf     – accumulation curves (one panel per clade)
  figures/pangenome_summary_bar.pdf– bar chart of final counts
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

GROUP_NAMES = {
    "BL": "Blastocladiomycota", "CH": "Chytridiomycetes",
    "EN": "Entomophthorales",   "GL": "Glomerales",
    "HO": "Holozoa",            "KX": "Kickxellomycotina",
    "MI": "Microsporidia",      "MO": "Mortierellomycota",
    "MU": "Mucoromycotina",     "NE": "Neocallimastigomycota",
    "NU": "Nucleariida",        "ZO": "Zoopagomycota",
    "AS": "Ascomycota",         "BA": "Basidiomycota",
}

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}


# =============================================================================
# I/O
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genecount", required=True,
                   help="OrthoFinder Orthogroups.GeneCount.tsv")
    p.add_argument("--clade_map", required=True,
                   help="Species → clade mapping TSV.")
    p.add_argument("--out_dir",   default="pangenome_out",
                   help="Output directory.")
    p.add_argument("--n_iter",    type=int, default=100,
                   help="Subsampling iterations (default: 100).")
    p.add_argument("--seed",      type=int, default=42)
    return p.parse_args()


def load_genecount(path: str, clade_map: dict) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.drop(columns=["Total"], errors="ignore")
    df.columns = [c.replace(".aa", "") for c in df.columns]

    # Keep only species with a known clade
    known = [c for c in df.columns if c in clade_map]
    df = df[known]
    # Convert to presence/absence (0/1)
    pa = (df > 0).astype(int)
    return pa


# =============================================================================
# Pangenome computation
# =============================================================================

def accumulation_curve(pa: pd.DataFrame, n_iter: int, rng: np.random.Generator):
    """
    Returns a dict with keys 'n', 'pan_mean', 'pan_sd',
    'core_mean', 'core_sd', 'shell_mean', 'shell_sd', 'unique_mean', 'unique_sd'.
    """
    species = pa.columns.tolist()
    n_sp = len(species)
    if n_sp < 2:
        return None

    records = {k: [] for k in
               ("pan", "core", "shell", "unique")}
    ns = list(range(1, n_sp + 1))

    for n in ns:
        pan_vals, core_vals, shell_vals, uniq_vals = [], [], [], []
        for _ in range(n_iter):
            sample = rng.choice(species, size=n, replace=False)
            sub = pa[sample]
            counts = sub.sum(axis=1)
            pan_vals.append((counts >= 1).sum())
            core_vals.append((counts == n).sum())
            shell_vals.append(((counts >= 2) & (counts < n)).sum())
            uniq_vals.append((counts == 1).sum())
        records["pan"].append((np.mean(pan_vals), np.std(pan_vals)))
        records["core"].append((np.mean(core_vals), np.std(core_vals)))
        records["shell"].append((np.mean(shell_vals), np.std(shell_vals)))
        records["unique"].append((np.mean(uniq_vals), np.std(uniq_vals)))

    result = {"n": ns}
    for key in ("pan", "core", "shell", "unique"):
        means, sds = zip(*records[key])
        result[f"{key}_mean"] = list(means)
        result[f"{key}_sd"]   = list(sds)
    return result


def run_pangenome(pa: pd.DataFrame, clade_map: dict, groups: list,
                  n_iter: int, rng: np.random.Generator) -> dict:
    results = {}
    for grp in groups:
        sp = [c for c in pa.columns if clade_map.get(c) == grp]
        if len(sp) < 3:
            continue
        curve = accumulation_curve(pa[sp], n_iter, rng)
        if curve:
            results[grp] = curve

    # combined
    curve_all = accumulation_curve(pa, n_iter, rng)
    if curve_all:
        results["ALL"] = curve_all
    return results


# =============================================================================
# Figures
# =============================================================================

def plot_curves(results: dict, out_dir: Path) -> None:
    n_groups = len(results)
    ncols = 4
    nrows = (n_groups + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axes = axes.flatten()

    for ax, (grp, data) in zip(axes, results.items()):
        ns = data["n"]
        color = CLADE_COLORS.get(grp, "#888888")
        for key, label, ls in [("pan", "Pangenome", "-"),
                                ("core", "Core", "--"),
                                ("shell", "Shell", ":"),
                                ("unique", "Unique", "-.")]:
            mn = data[f"{key}_mean"]
            sd = data[f"{key}_sd"]
            ax.plot(ns, mn, ls=ls, color=color, label=label, linewidth=1.5)
            ax.fill_between(ns,
                            [m - s for m, s in zip(mn, sd)],
                            [m + s for m, s in zip(mn, sd)],
                            alpha=0.15, color=color)
        ax.set_title(GROUP_NAMES.get(grp, grp), fontsize=9)
        ax.set_xlabel("Genomes", fontsize=8)
        ax.set_ylabel("Gene families", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in axes[n_groups:]:
        ax.set_visible(False)

    handles = [
        plt.Line2D([0], [0], ls="-",  color="grey", label="Pangenome"),
        plt.Line2D([0], [0], ls="--", color="grey", label="Core"),
        plt.Line2D([0], [0], ls=":",  color="grey", label="Shell"),
        plt.Line2D([0], [0], ls="-.", color="grey", label="Unique"),
    ]
    fig.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "pangenome_curves.pdf")
    fig.savefig(out_dir / "pangenome_curves.png", dpi=200)
    plt.close(fig)


def save_tsv(results: dict, out_dir: Path) -> None:
    summary_rows = []
    for grp, data in results.items():
        n_final = data["n"][-1]
        rows = {
            "group": grp, "n_genomes": n_final,
            "pan_mean":    data["pan_mean"][-1],    "pan_sd":    data["pan_sd"][-1],
            "core_mean":   data["core_mean"][-1],   "core_sd":   data["core_sd"][-1],
            "shell_mean":  data["shell_mean"][-1],  "shell_sd":  data["shell_sd"][-1],
            "unique_mean": data["unique_mean"][-1], "unique_sd": data["unique_sd"][-1],
        }
        summary_rows.append(rows)
        # full accumulation TSV per group
        curve_df = pd.DataFrame({k: v for k, v in data.items()})
        curve_df.to_csv(out_dir / f"{grp}_pangenome.tsv", sep="\t", index=False)

    pd.DataFrame(summary_rows).to_csv(out_dir / "summary_counts.tsv", sep="\t", index=False)


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clade_df = pd.read_csv(args.clade_map, sep="\t", header=None,
                           names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))

    print(f"Loading gene-count matrix from {args.genecount}")
    pa = load_genecount(args.genecount, clade_map)
    print(f"  {pa.shape[0]} orthogroups × {pa.shape[1]} species")

    groups = [g for g in CLADE_ORDER if any(clade_map.get(c) == g for c in pa.columns)]
    rng = np.random.default_rng(args.seed)

    print(f"Computing accumulation curves ({args.n_iter} iterations) ...")
    results = run_pangenome(pa, clade_map, groups, args.n_iter, rng)

    print("Saving TSVs ...")
    save_tsv(results, out_dir)

    print("Generating figures ...")
    plot_curves(results, fig_dir)

    print("Done.")


if __name__ == "__main__":
    main()
