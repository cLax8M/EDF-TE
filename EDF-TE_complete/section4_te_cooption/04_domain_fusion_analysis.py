#!/usr/bin/env python3
"""
04_domain_fusion_analysis.py  —  TE domain fusion protein characterisation.

Identifies and characterises proteins that carry both a TE-derived domain
and a non-TE cellular domain in the same polypeptide chain — "domain fusion"
proteins where a TE domain has been co-opted into a new multi-domain context.

Focal fusion classes (based on observed data):
  - FLYWCH (DNA-binding) fused to MULE transposase
  - CotH (kinase) fused to DNA transposon TIR domains
  - NACHT (NTPase) fused to reverse transcriptase
  - PARP (poly-ADP ribose) fused to chromo domain

Analysis:
  1. Screen InterProScan output for proteins with both TE and non-TE domains
  2. Classify fusion type by domain combination
  3. Count fusions per species; compute clade distribution
  4. For each fusion type: summarise phylogenetic distribution (heatmap)

Inputs:
  --ipr_dir       Directory with per-species InterPro TSV annotation
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory
  --te_domains    Comma-separated IPR IDs to flag as TE-related
  --fusion_defs   TSV defining fusion classes: class, te_ipr, partner_ipr, description

Outputs:
  domain_fusions.tsv
  fusion_clade_distribution.tsv
  figures/fusion_domain_heatmap.pdf
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Default fusion class definitions (if --fusion_defs not provided)
FUSION_DEFAULTS = [
    {"class": "FLYWCH-MULE",    "te_ipr": "IPR021967", "partner_ipr": "IPR003893",
     "description": "FLYWCH DNA-binding + MULE transposase"},
    {"class": "CotH-DDE",       "te_ipr": "IPR007721", "partner_ipr": "IPR009102",
     "description": "CotH kinase + DDE transposase"},
    {"class": "NACHT-RT",       "te_ipr": "IPR005093", "partner_ipr": "IPR002104",
     "description": "NACHT ATPase + chromo domain"},
    {"class": "PARP-Chromo",    "te_ipr": "IPR002104", "partner_ipr": "IPR012317",
     "description": "Chromo domain + PARP catalytic"},
]

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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ipr_dir",     required=True)
    p.add_argument("--clade_map",   required=True)
    p.add_argument("--out_dir",     default="domain_fusion_out")
    p.add_argument("--te_domains",  default=None)
    p.add_argument("--fusion_defs", default=None,
                   help="TSV defining fusion classes (optional; uses built-in defaults).")
    return p.parse_args()


def load_gene_domains(ipr_dir: Path, species_list: list[str]) -> pd.DataFrame:
    """Load gene-level domain lists across all species."""
    dfs = []
    for sp in species_list:
        files = list(Path(ipr_dir).glob(f"{sp}*.ipr.tab")) + \
                list(Path(ipr_dir).glob(f"{sp}*.tsv"))
        if not files:
            continue
        df = pd.read_csv(files[0], sep="\t", header=None, usecols=[0, 4, 5],
                         names=["gene", "ipr_id", "ipr_desc"]).dropna(subset=["ipr_id"])
        df["species"] = sp
        dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def detect_fusions(gene_domains: pd.DataFrame,
                   fusion_defs: list[dict]) -> pd.DataFrame:
    """
    For each fusion class, find genes carrying both te_ipr and partner_ipr.
    Returns a table of: species, gene, fusion_class, description.
    """
    # Build per-gene domain sets
    gene_ipr = gene_domains.groupby(["species", "gene"])["ipr_id"].apply(set).reset_index()
    gene_ipr.columns = ["species", "gene", "ipr_set"]

    rows = []
    for fdef in fusion_defs:
        te_ipr      = fdef["te_ipr"]
        partner_ipr = fdef["partner_ipr"]
        cls         = fdef["class"]
        desc        = fdef["description"]

        mask = gene_ipr["ipr_set"].apply(
            lambda s: te_ipr in s and partner_ipr in s)
        hits = gene_ipr[mask].copy()
        if hits.empty:
            continue
        hits["fusion_class"]  = cls
        hits["description"]   = desc
        hits["te_ipr"]        = te_ipr
        hits["partner_ipr"]   = partner_ipr
        rows.append(hits[["species", "gene", "fusion_class", "description",
                           "te_ipr", "partner_ipr"]])
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


# =============================================================================
# Figures
# =============================================================================

def plot_fusion_heatmap(fusions: pd.DataFrame, clade_map: dict,
                         out_dir: Path) -> None:
    fusions = fusions.copy()
    fusions["clade"] = fusions["species"].map(clade_map)

    # Count fusions per class × clade
    agg = fusions.groupby(["fusion_class", "clade"]).size().reset_index(name="n")
    mat = agg.pivot_table(index="fusion_class", columns="clade",
                          values="n", fill_value=0, aggfunc="sum")
    # Reorder columns
    mat = mat.reindex([c for c in CLADE_ORDER if c in mat.columns], axis=1)

    if mat.empty:
        print("No fusion data for heatmap.")
        return

    fig, ax = plt.subplots(figsize=(max(7, mat.shape[1] * 0.7),
                                    max(3, mat.shape[0] * 0.7 + 1)))
    sns.heatmap(mat, ax=ax, cmap="YlOrRd", linewidths=0.3, linecolor="#dddddd",
                annot=True, fmt="d", annot_kws={"fontsize": 7},
                cbar_kws={"label": "Number of fusion proteins"})
    ax.set_xlabel("Clade", fontsize=10)
    ax.set_ylabel("Fusion class", fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "fusion_domain_heatmap.pdf")
    fig.savefig(out_dir / "fusion_domain_heatmap.png", dpi=200)
    plt.close(fig)
    print("Saved: fusion_domain_heatmap.pdf")


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

    fusion_defs = FUSION_DEFAULTS
    if args.fusion_defs and Path(args.fusion_defs).exists():
        fd = pd.read_csv(args.fusion_defs, sep="\t")
        fusion_defs = fd.to_dict("records")

    print(f"Loading domain annotations for {len(species_list)} species ...")
    gene_domains = load_gene_domains(Path(args.ipr_dir), species_list)
    print(f"  {len(gene_domains):,} gene-domain entries")

    print("Detecting domain fusions ...")
    fusions = detect_fusions(gene_domains, fusion_defs)
    if fusions.empty:
        print("No domain fusions detected.", file=sys.stderr)
        sys.exit(0)

    print(f"  {len(fusions)} fusion protein instances across "
          f"{fusions['species'].nunique()} species")
    fusions["clade"] = fusions["species"].map(clade_map)
    fusions.to_csv(out_dir / "domain_fusions.tsv", sep="\t", index=False)

    clade_dist = fusions.groupby(["fusion_class", "clade"]).size().reset_index(name="n")
    clade_dist.to_csv(out_dir / "fusion_clade_distribution.tsv", sep="\t", index=False)

    plot_fusion_heatmap(fusions, clade_map, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
