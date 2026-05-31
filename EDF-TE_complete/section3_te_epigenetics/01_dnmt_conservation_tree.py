#!/usr/bin/env python3
"""
01_dnmt_conservation_tree.py  —  DNMT1/DIM2 conservation across Mucorales.

Builds a phylogenetic tree decorated with presence/absence of DNMT1 and DIM2
(RID2/DMAP-domain proteins) inferred from OrthoFinder orthogroups and
InterProScan domain calls.

Analysis:
  1. Load OrthoFinder gene-count matrix; identify DNMT1 orthogroup(s) by
     IPR005105 (Cytosine-5 DNA methyltransferase catalytic domain) and
     IPR018117 presence
  2. Load DIM2/RID2 orthogroups identified by IPR022686 (AdoMet-dependent MTase)
     combined with DMAP-interaction domain (IPR024773)
  3. Score each species: DNMT1_present (bool), DIM2_present (bool),
     copy_number (int)
  4. Map scores onto a pruned species tree (Mucoromycotina + outgroups)
  5. Plot cladogram with tip bar marks and dot matrix for DNMT1 / DIM2

Inputs:
  --orthogroups   OrthoFinder Gene_Counts.tsv (rows = orthogroups)
  --ipr_dir       Directory with per-species IPR annotation TSVs
  --dnmt_og       Comma-separated orthogroup IDs for DNMT1 (e.g. OG0000123)
  --dim2_og       Comma-separated orthogroup IDs for DIM2/RID2
  --species_tree  Newick species tree (full)
  --mucor_list    File listing Mucoromycotina species codes to prune to
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory

Output:
  dnmt_presence_matrix.tsv
  figures/dnmt_conservation_tree.pdf
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

DNMT1_IPR  = "IPR005105"
DIM2_IPR   = "IPR022686"
DIM2_DMAP  = "IPR024773"

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
    p.add_argument("--orthogroups", required=True)
    p.add_argument("--ipr_dir",     required=True)
    p.add_argument("--dnmt_og",     required=True,
                   help="Comma-separated DNMT1 orthogroup IDs")
    p.add_argument("--dim2_og",     required=True,
                   help="Comma-separated DIM2 orthogroup IDs")
    p.add_argument("--species_tree", required=True)
    p.add_argument("--mucor_list",  required=True)
    p.add_argument("--clade_map",   required=True)
    p.add_argument("--out_dir",     default="epigenetics_out")
    return p.parse_args()


def load_ipr_presence(ipr_dir: Path, species_list: list[str],
                       target_ipr: set[str]) -> pd.DataFrame:
    rows = []
    for sp in species_list:
        files = list(Path(ipr_dir).glob(f"{sp}*.ipr.tab")) + \
                list(Path(ipr_dir).glob(f"{sp}*.tsv"))
        if not files:
            rows.append({"species": sp, **{k: False for k in target_ipr}})
            continue
        df = pd.read_csv(files[0], sep="\t", header=None, usecols=[0, 4],
                         names=["gene", "ipr_id"]).dropna()
        row = {"species": sp}
        for ipr in target_ipr:
            row[ipr] = ipr in df["ipr_id"].values
        rows.append(row)
    return pd.DataFrame(rows)


def load_og_counts(og_file: Path, og_ids: list[str],
                   species_list: list[str]) -> pd.Series:
    df = pd.read_csv(og_file, sep="\t", index_col=0)
    df.index = df.index.str.strip()
    sub = df.loc[df.index.isin(og_ids), :]
    # Keep only columns (species) in our list
    cols = [c for c in df.columns if c in species_list]
    if sub.empty or not cols:
        return pd.Series(0, index=species_list)
    return sub[cols].sum(axis=0).reindex(species_list, fill_value=0)


def plot_dnmt_tree(tree_nwk: str, mucor_species: list[str],
                   presence: pd.DataFrame, out_dir: Path,
                   clade_map: dict) -> None:
    try:
        from ete3 import Tree, TreeStyle, NodeStyle, TextFace, RectFace
    except ImportError:
        print("ete3 not available — skipping tree plot; writing TSV only.")
        return

    t = Tree(tree_nwk, format=1)
    # Prune to Mucoromycotina + a few outgroup leaves
    keep = set(mucor_species)
    t.prune([n for n in t.get_leaves() if n.name in keep], preserve_branch_length=True)

    # Map presence to leaf nodes
    pres_dict = presence.set_index("species").to_dict("index")

    for leaf in t.get_leaves():
        sp = leaf.name
        clade = clade_map.get(sp, "MU")
        nst = NodeStyle()
        nst["fgcolor"] = CLADE_COLORS.get(clade, "#333333")
        nst["size"] = 6
        leaf.set_style(nst)

        info = pres_dict.get(sp, {})
        for marker, col in [("DNMT1", "#e41a1c"), ("DIM2", "#377eb8")]:
            key = "dnmt_present" if marker == "DNMT1" else "dim2_present"
            if info.get(key, False):
                leaf.add_face(RectFace(8, 8, col, col), column=1 if marker == "DNMT1" else 2)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "r"
    ts.scale = 120
    t.render(str(out_dir / "dnmt_conservation_tree.pdf"), tree_style=ts)
    print("Saved: dnmt_conservation_tree.pdf")


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))

    mucor_species = [l.strip() for l in open(args.mucor_list) if l.strip()]

    dnmt_ogs = [o.strip() for o in args.dnmt_og.split(",")]
    dim2_ogs = [o.strip() for o in args.dim2_og.split(",")]

    print("Loading OrthoFinder counts ...")
    dnmt_counts = load_og_counts(Path(args.orthogroups), dnmt_ogs, mucor_species)
    dim2_counts  = load_og_counts(Path(args.orthogroups), dim2_ogs, mucor_species)

    print("Loading IPR domain presence ...")
    ipr_pres = load_ipr_presence(
        Path(args.ipr_dir), mucor_species,
        target_ipr={DNMT1_IPR, DIM2_IPR, DIM2_DMAP})

    # Combine OG count and IPR evidence
    pres = ipr_pres.copy()
    pres["dnmt_present"] = (dnmt_counts > 0).values | pres[DNMT1_IPR]
    pres["dim2_present"]  = (dim2_counts  > 0).values | (pres[DIM2_IPR] & pres[DIM2_DMAP])
    pres["dnmt_copies"]   = dnmt_counts.values
    pres["dim2_copies"]   = dim2_counts.values

    pres.to_csv(out_dir / "dnmt_presence_matrix.tsv", sep="\t", index=False)

    tree_nwk = open(args.species_tree).read()
    plot_dnmt_tree(tree_nwk, mucor_species, pres, fig_dir, clade_map)

    print("Done.")


if __name__ == "__main__":
    main()
