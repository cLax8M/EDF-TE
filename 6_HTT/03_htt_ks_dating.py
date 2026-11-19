#!/usr/bin/env python3

# HTT Ks dating of horizontal TE transfer candidates
# uses OrthoFinder Ks + species tree calibration

# imports
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# clade colors
CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

# clade order
CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

# args
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--candidates",   required=True)
    p.add_argument("--og_ks_dir",    required=True)
    p.add_argument("--species_tree", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="htt_dating_out")
    return p.parse_args()

# MRCA depth (tree calibration)
def mrca_depth(nwk_file: str, sp_a: str, sp_b: str) -> float:
    """Return the root-to-MRCA distance (branch length sum) for two tips."""
    try:
        from ete3 import Tree
        t   = Tree(nwk_file, format=1)
        anc = t.get_common_ancestor(sp_a, sp_b)
        return t.get_distance(t, anc)
    except Exception:
        return np.nan

# load Ks between clades (OrthoFinder output)
def load_og_ks(og_ks_dir: str, clade_a: str, clade_b: str) -> float:
    """Return median Ks from OrthoFinder orthogroups between two clades."""
    ks_files = list(Path(og_ks_dir).glob(f"{clade_a}_{clade_b}*.tsv")) + \
               list(Path(og_ks_dir).glob(f"{clade_b}_{clade_a}*.tsv"))
    if not ks_files:
        return np.nan
    df = pd.read_csv(ks_files[0], sep="\t")
    ks_col = next((c for c in df.columns if c.lower() == "ks"), None)
    if ks_col is None:
        return np.nan
    vals = df[ks_col].dropna().values
    vals = vals[(vals > 0) & (vals < 10)]
    return float(np.median(vals)) if len(vals) > 0 else np.nan


# Figures

# violin plot (Ks distribution per clade pair)
def plot_timing_violin(dated: pd.DataFrame, out_dir: Path) -> None:
    by_clade_pair = {
        cp: sub["median_ks"].dropna().values
        for cp, sub in dated.groupby("clade_pair")
        if sub["median_ks"].notna().sum() >= 3
    }
    if not by_clade_pair:
        return

    fig, ax = plt.subplots(figsize=(max(6, len(by_clade_pair) * 0.8), 5))
    labels = sorted(by_clade_pair.keys())
    data   = [by_clade_pair[l] for l in labels]
    vp     = ax.violinplot(data, showmedians=True)
    for i, body in enumerate(vp["bodies"]):
        body.set_alpha(0.6)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("MRCA Ks (transfer time proxy)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "htt_timing_violin.pdf")
    plt.close(fig)
    print("Saved: htt_timing_violin.pdf")


# heatmap (HTT counts per clade pair)
def plot_heatmap(heatmap_df: pd.DataFrame, out_dir: Path) -> None:
    """Inter-clade heatmap: transfer counts."""
    clades_present = sorted(set(heatmap_df["q_clade"]) | set(heatmap_df["s_clade"]))
    clades = [c for c in CLADE_ORDER if c in clades_present] + \
             [c for c in clades_present if c not in CLADE_ORDER]
    n = len(clades)
    if n < 2:
        return

    mat = pd.DataFrame(0, index=clades, columns=clades, dtype=float)
    for _, r in heatmap_df.iterrows():
        if r["q_clade"] in mat.index and r["s_clade"] in mat.columns:
            mat.loc[r["q_clade"], r["s_clade"]] = r["n_candidates"]
            mat.loc[r["s_clade"], r["q_clade"]] = r["n_candidates"]

    fig, ax = plt.subplots(figsize=(n * 0.55 + 1.5, n * 0.55 + 1.5))
    im = ax.imshow(mat.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(clades, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(clades, fontsize=8)
    plt.colorbar(im, ax=ax, label="HTT candidates", shrink=0.7)
    ax.set_title("Inter-clade HTT candidate counts", fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "htt_clade_heatmap.pdf")
    plt.close(fig)
    print("Saved: htt_clade_heatmap.pdf")

# main

def main():
    # parse args
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    # clade map
    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))
    # load candidates
    cands = pd.read_csv(args.candidates, sep="\t")
    if cands.empty:
        print("Empty candidates file.", file=sys.stderr)
        sys.exit(0)
    # compute Ks per event
    rows = []
    for _, r in cands.iterrows():
        ks = load_og_ks(args.og_ks_dir, r["q_clade"], r["s_clade"])
        rows.append({
            "q_clade":     r["q_clade"],
            "s_clade":     r["s_clade"],
            "qseqid":      r.get("qseqid", ""),
            "sseqid":      r.get("sseqid", ""),
            "pident":      r.get("pident", np.nan),
            "anomaly":     r.get("anomaly_score", np.nan),
            "median_ks":   ks,
            "clade_pair":  "_".join(sorted([r["q_clade"], r["s_clade"]])),
        })

    dated = pd.DataFrame(rows)
    dated.to_csv(out_dir / "htt_ks_dating.tsv", sep="\t", index=False)

    # aggregate per clade pair
    heatmap_df = dated.groupby(["q_clade", "s_clade"]).agg(
        n_candidates=("pident", "count"),
        mean_ks=("median_ks", "mean"),
    ).reset_index()

    # plots
    plot_timing_violin(dated, fig_dir)
    plot_heatmap(heatmap_df, fig_dir)
    
    # summary
    total = len(dated)
    n_dated = dated["median_ks"].notna().sum()
    print(f"{total} candidate transfers; {n_dated} with Ks dating")
    print("Done.")

# run
if __name__ == "__main__":
    main()
