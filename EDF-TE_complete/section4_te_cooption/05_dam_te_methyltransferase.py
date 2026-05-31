#!/usr/bin/env python3
"""
05_dam_te_methyltransferase.py  —  DAM-like N6-adenine methyltransferase in TEs.

Identifies TE-fused DAM methyltransferase sequences and characterises their
phylogenetic distribution and association with AMT1 (adenine methylation
maintenance) proteins.

Key findings to reproduce:
  - 15 species carry TE-fused DAM methyltransferases (14 DIRS, 1 Gypsy)
  - Fisher exact: DAM-TE × AMT1 co-occurrence (reported OR=5.19, p=0.028)
  - DAM-TE is distributed across multiple fungal clades

Analysis:
  1. Load HMMER hits against DAM methyltransferase HMM profiles
  2. Confirm TE association by overlap with TE BED coordinates
  3. Classify TE class for DAM-carrying TEs
  4. Load AMT1 presence (IPR-based) per species
  5. Fisher test: DAM-TE presence × AMT1 presence across all species
  6. Plot: clade distribution heatmap, TE class composition, AMT1 co-occurrence

Inputs:
  --hmm_hits      HMMER tblout file (all species combined, or per-species dir)
  --te_bed_dir    EarlGrey filtered TE BED files
  --ipr_dir       InterProScan output per species (for AMT1 detection)
  --gene_faa_dir  Per-species protein FASTA (for sequence lengths)
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory
  --hmm_evalue    HMMER e-value cutoff (default: 1e-5)
  --amt1_ipr      IPR ID for AMT1 (default: IPR004550 — AdoMet-MTase type-1)

Outputs:
  dam_te_summary.tsv
  amt1_cooccurrence.tsv
  figures/dam_te_clade_heatmap.pdf
  figures/dam_te_class_bar.pdf
  figures/dam_amt1_fisher.pdf
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AMT1_IPR_DEFAULT = "IPR004550"

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
    "DIRS":           "#984ea3",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hmm_hits",     required=True)
    p.add_argument("--te_bed_dir",   required=True)
    p.add_argument("--ipr_dir",      required=True)
    p.add_argument("--gene_faa_dir", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="dam_te_out")
    p.add_argument("--hmm_evalue",   type=float, default=1e-5)
    p.add_argument("--amt1_ipr",     default=AMT1_IPR_DEFAULT)
    return p.parse_args()


def _te_class(f: str) -> str:
    f = f.upper()
    if "DIRS" in f:   return "DIRS"
    if "LTR" in f:    return "LTR"
    if "GYPSY" in f:  return "LTR/Gypsy"
    if "DNA" in f:    return "DNA"
    if "LINE" in f:   return "LINE"
    if "SINE" in f:   return "SINE"
    if "RC/" in f:    return "Rolling Circle"
    if "PENELOPE" in f: return "Penelope"
    return "Unclassified"


def parse_hmmer_tblout(tblout_path: str, evalue_thresh: float) -> pd.DataFrame:
    """Parse HMMER --tblout format."""
    rows = []
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            target_name = parts[0]
            query_name  = parts[2]
            full_evalue = float(parts[4])
            full_score  = float(parts[5])
            if full_evalue <= evalue_thresh:
                rows.append({"target": target_name, "query": query_name,
                             "evalue": full_evalue, "score": full_score})
    return pd.DataFrame(rows)


def extract_species_from_header(header: str, species_list: list[str]) -> str:
    """Infer species code from a FASTA protein header."""
    for sp in species_list:
        if header.startswith(sp) or f"|{sp}|" in header or f"_{sp}_" in header:
            return sp
    # Fallback: first two segments separated by '_'
    parts = header.split("_")
    if len(parts) >= 2:
        return "_".join(parts[:2])
    return header[:8]


def load_te_bed_by_chrom(te_bed_dir: str, sp: str) -> pd.DataFrame:
    files = list(Path(te_bed_dir).glob(f"{sp}*.filteredRepeats.bed")) + \
            list(Path(te_bed_dir).glob(f"{sp}*.bed"))
    if not files:
        return pd.DataFrame()
    df = pd.read_csv(files[0], sep="\t", header=None, usecols=[0, 1, 2, 3],
                     names=["chrom", "start", "end", "family"])
    df["te_class"] = df["family"].apply(_te_class)
    return df


def load_amt1_presence(ipr_dir: str, species_list: list[str],
                        amt1_ipr: str) -> pd.Series:
    presence = {}
    for sp in species_list:
        files = list(Path(ipr_dir).glob(f"{sp}*.ipr.tab")) + \
                list(Path(ipr_dir).glob(f"{sp}*.tsv"))
        if not files:
            presence[sp] = False
            continue
        df = pd.read_csv(files[0], sep="\t", header=None, usecols=[4],
                         names=["ipr_id"]).dropna()
        presence[sp] = amt1_ipr in df["ipr_id"].values
    return pd.Series(presence)


# =============================================================================
# Figures
# =============================================================================

def plot_clade_heatmap(summary: pd.DataFrame, clade_map: dict,
                        out_dir: Path) -> None:
    summary = summary.copy()
    summary["clade"] = summary["species"].map(clade_map)
    agg = summary.groupby(["clade", "te_class"]).size().reset_index(name="n")
    mat = agg.pivot_table(index="clade", columns="te_class",
                          values="n", fill_value=0, aggfunc="sum")
    mat = mat.reindex([c for c in CLADE_ORDER if c in mat.index])

    import seaborn as sns
    fig, ax = plt.subplots(figsize=(max(5, mat.shape[1] * 1.2),
                                    max(4, mat.shape[0] * 0.6 + 1)))
    sns.heatmap(mat, ax=ax, cmap="YlOrRd", annot=True, fmt="d",
                linewidths=0.3, linecolor="#dddddd",
                cbar_kws={"label": "Number of DAM-TE species"})
    ax.set_xlabel("TE class", fontsize=10)
    ax.set_ylabel("Clade", fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "dam_te_clade_heatmap.pdf")
    plt.close(fig)
    print("Saved: dam_te_clade_heatmap.pdf")


def plot_fisher_result(or_val: float, pval: float, n_dam_amt1: int,
                        n_dam_noamt1: int, n_nodam_amt1: int,
                        n_nodam_noamt1: int, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    cats = ["DAM-TE +\nAMT1 +", "DAM-TE +\nAMT1 −",
            "DAM-TE −\nAMT1 +", "DAM-TE −\nAMT1 −"]
    vals = [n_dam_amt1, n_dam_noamt1, n_nodam_amt1, n_nodam_noamt1]
    colors = ["#e41a1c", "#fb9a99", "#a6cee3", "#1f78b4"]
    ax.bar(range(4), vals, color=colors, edgecolor="white")
    ax.set_xticks(range(4))
    ax.set_xticklabels(cats, fontsize=8)
    ax.set_ylabel("Number of species", fontsize=10)
    ax.set_title(f"DAM-TE × AMT1 co-occurrence\nOR={or_val:.2f}, p={pval:.3g}",
                 fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "dam_amt1_fisher.pdf")
    plt.close(fig)
    print(f"Saved: dam_amt1_fisher.pdf  (OR={or_val:.2f}, p={pval:.3g})")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))
    species_list = list(clade_map.keys())

    print("Parsing HMMER hits ...")
    hmm_hits = parse_hmmer_tblout(args.hmm_hits, args.hmm_evalue)
    if hmm_hits.empty:
        print("No HMMER hits passed e-value threshold.", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(hmm_hits)} hits")

    # Assign species from header
    hmm_hits["species"] = hmm_hits["target"].apply(
        lambda h: extract_species_from_header(h, species_list))

    # Intersect with TE coordinates
    rows = []
    for sp, sp_hits in hmm_hits.groupby("species"):
        if sp not in clade_map:
            continue
        te_bed = load_te_bed_by_chrom(args.te_bed_dir, sp)
        if te_bed.empty:
            continue
        # Load protein→genome mapping — attempt from header naming convention:
        # header format expected: <scaffold>_<start>_<end>|<gene_id>
        for _, hit in sp_hits.iterrows():
            parts = hit["target"].split("|")[0].split("_")
            te_class = "Unknown"
            if len(parts) >= 3:
                try:
                    chrom = "_".join(parts[:-2])
                    s, e  = int(parts[-2]), int(parts[-1])
                    ov = te_bed[(te_bed["chrom"] == chrom) &
                                (te_bed["end"] > s) & (te_bed["start"] < e)]
                    if not ov.empty:
                        te_class = ov.iloc[0]["te_class"]
                except (ValueError, IndexError):
                    pass
            rows.append({"species": sp, "gene": hit["target"],
                          "evalue": hit["evalue"], "score": hit["score"],
                          "te_class": te_class})

    summary = pd.DataFrame(rows)
    summary = summary[summary["te_class"] != "Unknown"]
    print(f"  {len(summary)} DAM hits confirmed in TE coordinates")
    print(f"  {summary['species'].nunique()} species with DAM-TE")

    summary.to_csv(out_dir / "dam_te_summary.tsv", sep="\t", index=False)

    # AMT1 co-occurrence
    print("Testing AMT1 co-occurrence ...")
    amt1 = load_amt1_presence(args.ipr_dir, species_list, args.amt1_ipr)
    dam_species = set(summary["species"].unique())

    n11 = sum(1 for sp in species_list if sp in dam_species and amt1.get(sp, False))
    n10 = sum(1 for sp in species_list if sp in dam_species and not amt1.get(sp, False))
    n01 = sum(1 for sp in species_list if sp not in dam_species and amt1.get(sp, False))
    n00 = sum(1 for sp in species_list if sp not in dam_species and not amt1.get(sp, False))

    od, pval = fisher_exact([[n11, n10], [n01, n00]], alternative="greater")
    print(f"  Fisher exact: OR={od:.2f}, p={pval:.4g}")

    cooc = pd.DataFrame([
        {"category": "DAM+_AMT1+", "n": n11},
        {"category": "DAM+_AMT1-", "n": n10},
        {"category": "DAM-_AMT1+", "n": n01},
        {"category": "DAM-_AMT1-", "n": n00},
        {"category": "odds_ratio",  "n": od},
        {"category": "p_value",     "n": pval},
    ])
    cooc.to_csv(out_dir / "amt1_cooccurrence.tsv", sep="\t", index=False)

    plot_clade_heatmap(summary, clade_map, fig_dir)
    plot_fisher_result(od, pval, n11, n10, n01, n00, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
