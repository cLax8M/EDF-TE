#!/usr/bin/env python3

# analyses ChimeraTE output to identify chimeric transcripts and classify TE fusions

# imports
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# define clade structure
CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#984ea3",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}

# parse args
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chimte_dir",   required=True)
    p.add_argument("--species_list", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="chimte_out")
    p.add_argument("--min_reads",    type=int, default=3)
    p.add_argument("--min_samples",  type=int, default=2)
    return p.parse_args()

# Classify TEs
def _te_class(f: str) -> str:
    f = f.upper()
    if "LTR" in f:    return "LTR"
    if "DNA" in f:    return "DNA"
    if "LINE" in f:   return "LINE"
    if "SINE" in f:   return "SINE"
    if "RC/" in f or "HELITRON" in f: return "Rolling Circle"
    if "PENELOPE" in f: return "Penelope"
    return "Unclassified"


# load one species ChimeraTE file
def load_chimte_results(chimte_dir: str, sp: str,
                         min_reads: int) -> pd.DataFrame:
    """
    Load ChimeraTE output for one species.
    Expected file: chimte_dir/<sp>/chimTE.results.tsv
    Columns: gene_id, te_family, chimera_type, n_reads, strand, ...
    """
    candidates = [
        Path(chimte_dir) / sp / "chimTE.results.tsv",
        Path(chimte_dir) / f"{sp}_chimTE.results.tsv",
        *list(Path(chimte_dir).glob(f"{sp}*chimTE*.tsv")),
        *list(Path(chimte_dir).glob(f"{sp}*results*.tsv")),
    ]
    for f in candidates:
        if f.exists():
            df = pd.read_csv(f, sep="\t")
            df.columns = [c.lower().replace(" ", "_") for c in df.columns]
            if "n_reads" in df.columns:
                df = df[df["n_reads"] >= min_reads]
            return df
    return pd.DataFrame()

# aggregate all species into single dataframe
def aggregate_results(chimte_dir: str, species_list: list[str],
                       clade_map: dict, min_reads: int) -> pd.DataFrame:
    all_dfs = []
    for sp in species_list:
        df = load_chimte_results(chimte_dir, sp, min_reads)
        if df.empty:
            continue
        df["species"] = sp
        df["clade"]   = clade_map.get(sp, "?")
        # Normalise TE class column name
        for col in ["te_family", "te_class", "repeat_family"]:
            if col in df.columns:
                df["te_class"] = df[col].apply(_te_class)
                break
        all_dfs.append(df)
    if not all_dfs:
        return pd.DataFrame()
    return pd.concat(all_dfs, ignore_index=True)

# identify conserved chimeras across species
def shared_chimeras(df: pd.DataFrame, min_samples: int) -> pd.DataFrame:
    """Identify chimeric pairs seen in ≥ min_samples species."""
    if "gene_id" not in df.columns or "te_class" not in df.columns:
        return pd.DataFrame()
    pair = df["gene_id"].astype(str) + "||" + df["te_class"].astype(str)
    counts = pair.value_counts()
    conserved = counts[counts >= min_samples]
    return conserved.reset_index().rename(columns={"index": "chimera_pair",
                                                     "count": "n_species"})



# Figures


def plot_per_species_bar(df: pd.DataFrame, out_dir: Path) -> None:
    sp_counts = df.groupby(["species", "clade"]).size().reset_index(name="n_chimeras")
    # sort by clade order
    sp_counts = sp_counts.sort_values(
        "clade", key=lambda s: s.map({c: i for i, c in enumerate(CLADE_ORDER)}))

    fig, ax = plt.subplots(figsize=(max(8, len(sp_counts) * 0.35), 5))
    colors = [CLADE_COLORS.get(c, "#999") for c in sp_counts["clade"]]
    ax.bar(range(len(sp_counts)), sp_counts["n_chimeras"], color=colors, width=0.8)
    ax.set_xticks(range(len(sp_counts)))
    ax.set_xticklabels(sp_counts["species"], rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Number of chimeric transcripts", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    # Clade legend
    seen = {}
    for c in sp_counts["clade"].unique():
        if c not in seen:
            seen[c] = mpatches.Patch(color=CLADE_COLORS.get(c, "#999"), label=c)
            
    # legend by clade
    import matplotlib.patches as mpatches2
    ax.legend(handles=list(seen.values()), fontsize=7, frameon=False,
              ncol=4, loc="upper right")

    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "chimte_per_species_bar.pdf")
    fig.savefig(out_dir / "chimte_per_species_bar.png", dpi=200)
    plt.close(fig)
    print("Saved: chimte_per_species_bar.pdf")


import matplotlib.patches as mpatches


def plot_te_class_proportions(df: pd.DataFrame, out_dir: Path) -> None:
    clades = [c for c in CLADE_ORDER if c in df["clade"].values]
    class_clade = df.groupby(["clade", "te_class"]).size().unstack("te_class", fill_value=0)

    # Normalise rows
    class_clade = class_clade.divide(class_clade.sum(axis=1), axis=0) * 100
    class_clade = class_clade.reindex([c for c in clades if c in class_clade.index])

    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = np.zeros(len(class_clade))
    for cls in class_clade.columns:
        color = TE_CLASS_COLORS.get(cls, "#999")
        ax.bar(range(len(class_clade)), class_clade[cls], bottom=bottom,
               color=color, label=cls, width=0.8)
        bottom += class_clade[cls].values

    ax.set_xticks(range(len(class_clade)))
    ax.set_xticklabels(class_clade.index, fontsize=9, rotation=30)
    ax.set_ylabel("% chimeric transcripts by TE class", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "chimte_te_class_proportions.pdf")
    plt.close(fig)
    print("Saved: chimte_te_class_proportions.pdf")


# main


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    # load clade mapping
    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))
    # load species list
    species_list = [l.strip() for l in open(args.species_list) if l.strip()]

    print(f"Loading ChimeraTE results for {len(species_list)} species ...")
    # aggregate all data
    df = aggregate_results(args.chimte_dir, species_list, clade_map, args.min_reads)
    if df.empty:
        print("No ChimeraTE data loaded — check --chimte_dir.", file=sys.stderr)
        sys.exit(1)

    print(f"  {len(df)} chimeric transcripts across "
          f"{df['species'].nunique()} species")
          
    # save full table
    df.to_csv(out_dir / "chimte_summary.tsv", sep="\t", index=False)
    
    # per clade summary
    clade_sum = df.groupby("clade").size().reset_index(name="n_chimeras")
    clade_sum.to_csv(out_dir / "chimte_per_clade.tsv", sep="\t", index=False)

    # shared events
    shared = shared_chimeras(df, args.min_samples)
    if not shared.empty:
        shared.to_csv(out_dir / "chimte_shared.tsv", sep="\t", index=False)
        print(f"  {len(shared)} chimeras shared ≥{args.min_samples} species")
    # plots
    plot_per_species_bar(df, fig_dir)
    plot_te_class_proportions(df, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
