#!/usr/bin/env python3
"""
00_prep_data.py
---------------
Prepare input matrices for treeWAS analysis.

Outputs (written to --outdir, default: TreeWAS/results/):
  phenotypes.tsv   — 183 × 8 phenotype table (Y)
  X_OG.tsv         — 183 × N_OG binary presence/absence matrix
  X_IPR.tsv        — 183 × N_IPR binary presence/absence matrix
  tree_pruned.nwk  — Newick tree pruned to 183 species

Usage:
  python 00_prep_data.py [options]

  --genomic_features  Path to ALL_SPECIES_ALL_SUMMARY.txt
  --duplication       Path to orthogroup_stats.tsv
  --unknown           Path to unknown_function_summary.tsv
  --synteny           Path to perspecies_synteny.tsv
  --orthogroups       Path to Orthogroups.GeneCount.tsv
  --ipr_dir           Directory containing *_IPR.tab files
  --tree              Path to SpeciesTree_rooted.txt (Newick)
  --outdir            Output directory (default: ../results)
  --min_freq          Minimum frequency for OG/IPR filter (default: 0.05)
"""

import argparse
import os
import sys
import math
import glob

import pandas as pd

# ── attempt to import ete3; graceful error if missing ──────────────────────
try:
    from ete3 import Tree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False

# ── species to exclude ──────────────────────────────────────────────────────
EXCLUDE_PREFIXES = ("AS_", "BA_", "HO_", "OL_")
EXCLUDE_FRAGMENTED = {
    "BA_Rhisol1", "CH_Batsa1", "GL_Rhicl1",
    "MU_Bifad1", "MU_Mucam1", "MU_Muccirc1",
}

# ── TE-related Pfam IDs (for annotation in script 02) ──────────────────────
TE_PFAM = {"PF00078", "PF03184", "PF13456", "PF14214", "PF13683", "PF01497", "PF13952"}


def parse_args():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_base = "/legserv/Temp/Carlos/EDF-TE"
    p = argparse.ArgumentParser(description="Prepare treeWAS input matrices")
    p.add_argument("--genomic_features",
                   default=f"{data_base}/results/GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY.txt")
    p.add_argument("--duplication",
                   default=f"{data_base}/results/DuplicationCG/orthogroup_stats.tsv")
    p.add_argument("--unknown",
                   default=f"{data_base}/results/unknown_function/unknown_function_summary.tsv")
    p.add_argument("--synteny",
                   default=f"{data_base}/results/Synteny/perspecies_synteny.tsv")
    p.add_argument("--orthogroups",
                   default=f"{data_base}/results/Orthofinder_run/Orthogroups/Orthogroups.GeneCount.tsv")
    p.add_argument("--ipr_dir",
                   default=f"{data_base}/data/functional_interpro")
    p.add_argument("--tree",
                   default=f"{data_base}/results/Orthofinder_run/Species_Tree/SpeciesTree_rooted.txt")
    p.add_argument("--outdir",
                   default=os.path.join(base, "results"))
    p.add_argument("--min_freq", type=float, default=0.05,
                   help="Minimum species frequency to retain an OG/IPR (default: 0.05)")
    return p.parse_args()


def should_exclude(species: str) -> bool:
    """Return True if this species should be excluded from the analysis."""
    if any(species.startswith(pref) for pref in EXCLUDE_PREFIXES):
        return True
    if species in EXCLUDE_FRAGMENTED:
        return True
    return False


# ── 1. PHENOTYPE MATRIX ─────────────────────────────────────────────────────

def load_phenotypes(args) -> pd.DataFrame:
    print("[1/5] Loading phenotype data …")

    # --- genomic features ---
    gf = pd.read_csv(args.genomic_features, sep="\t")
    gf = gf[["Species", "Genome_Size_bp", "Avg_Exons_Per_Gene",
              "Avg_Introns_Multiexon", "Pct_Intergenic"]].copy()
    gf.rename(columns={"Species": "species"}, inplace=True)

    # --- duplication ---
    dup = pd.read_csv(args.duplication, sep="\t")
    dup = dup[["Species", "Pct_paralog_genes", "Pct_OGs_multicopy"]].copy()
    dup.rename(columns={"Species": "species"}, inplace=True)

    # --- unknown function ---
    unk = pd.read_csv(args.unknown, sep="\t")
    unk = unk[["Species", "Percent_unknown"]].copy()
    unk.rename(columns={"Species": "species"}, inplace=True)

    # --- synteny ---
    syn = pd.read_csv(args.synteny, sep="\t")
    # for fragmented species mean_r_intra should be NaN, not 0
    syn.loc[syn["fragmented"].astype(str).str.lower() == "true", "mean_r_intra"] = float("nan")
    syn = syn[["species", "mean_r_intra"]].copy()

    # --- merge ---
    df = gf.merge(dup, on="species", how="outer")
    df = df.merge(unk, on="species", how="outer")
    df = df.merge(syn, on="species", how="outer")

    # --- apply exclusions ---
    before = len(df)
    df = df[~df["species"].apply(should_exclude)].copy()
    print(f"    Excluded {before - len(df)} species → {len(df)} remaining")

    # --- log10 transform genome size ---
    df["log10_Genome_Size"] = df["Genome_Size_bp"].apply(
        lambda x: math.log10(x) if pd.notna(x) and x > 0 else float("nan")
    )
    df.drop(columns=["Genome_Size_bp"], inplace=True)

    # --- reorder columns ---
    cols = ["species", "log10_Genome_Size", "Pct_Intergenic",
            "Avg_Introns_Multiexon", "Avg_Exons_Per_Gene",
            "Pct_paralog_genes", "Pct_OGs_multicopy",
            "Percent_unknown", "mean_r_intra"]
    df = df[cols]
    df.set_index("species", inplace=True)

    print(f"    Phenotype table: {df.shape[0]} species × {df.shape[1]} traits")
    missing = df.isnull().sum()
    if missing.any():
        print("    Missing values per trait:")
        for col, n in missing[missing > 0].items():
            print(f"      {col}: {n}")
    return df


# ── 2. OG BINARY MATRIX ─────────────────────────────────────────────────────

def load_og_matrix(args, species_set: set, min_count: int) -> pd.DataFrame:
    print("[2/5] Loading OrthoGroup matrix …")
    # OG file is large: read in chunks to manage memory
    chunks = []
    reader = pd.read_csv(args.orthogroups, sep="\t", index_col=0, chunksize=5000)
    for chunk in reader:
        # strip .aa/.faa suffix from column names before filtering
        chunk.columns = [c.replace(".aa", "").replace(".faa", "") for c in chunk.columns]
        cols_keep = [c for c in chunk.columns if c in species_set]
        if not cols_keep:
            continue
        sub = chunk[cols_keep]
        # binarise: ≥1 → 1
        sub = (sub >= 1).astype("int8")
        # filter rows: present in ≥ min_count species in this chunk
        sub = sub[sub.sum(axis=1) >= min_count]
        if not sub.empty:
            chunks.append(sub)

    if not chunks:
        raise RuntimeError("No OrthoGroups passed the frequency filter — check paths and species names.")

    og = pd.concat(chunks)
    # recompute filter after concatenation (chunks filtered independently)
    og = og[[c for c in og.columns if c in species_set]]
    og = og[og.sum(axis=1) >= min_count]
    # transpose: species × OGs
    og = og.T
    print(f"    OG matrix: {og.shape[0]} species × {og.shape[1]} orthogroups (freq ≥ {min_count}/{len(species_set)})")
    return og


# ── 3. IPR BINARY MATRIX ────────────────────────────────────────────────────

def load_ipr_matrix(args, species_set: set, min_count: int) -> pd.DataFrame:
    print("[3/5] Loading InterPro domain matrix …")
    ipr_files = glob.glob(os.path.join(args.ipr_dir, "*_IPR.tab"))
    if not ipr_files:
        raise FileNotFoundError(f"No *_IPR.tab files found in {args.ipr_dir}")

    records = {}
    for fpath in ipr_files:
        # derive species name from filename: AP_Amoocci_IPR.tab → AP_Amoocci
        basename = os.path.basename(fpath)
        species = basename.replace("_IPR.tab", "")
        if species not in species_set:
            continue
        try:
            df = pd.read_csv(fpath, sep="\t", header=None,
                             usecols=[1],  # col 1 = iprId
                             names=["iprId"])
            unique_iprs = set(df["iprId"].dropna().unique())
            records[species] = unique_iprs
        except Exception as e:
            print(f"    WARNING: could not read {fpath}: {e}")

    if not records:
        raise RuntimeError("No IPR files matched included species.")

    # get universe of IPR terms
    all_iprs = sorted(set().union(*records.values()))
    # build binary matrix
    data = {sp: pd.Series({ipr: 1 for ipr in iprs}) for sp, iprs in records.items()}
    ipr_df = pd.DataFrame(data, index=all_iprs).T.fillna(0).astype("int8")
    # filter columns (IPR terms) by minimum frequency
    ipr_df = ipr_df.loc[:, ipr_df.sum(axis=0) >= min_count]
    print(f"    IPR matrix: {ipr_df.shape[0]} species × {ipr_df.shape[1]} IPR terms (freq ≥ {min_count}/{len(species_set)})")
    return ipr_df


# ── 4. PRUNE TREE ───────────────────────────────────────────────────────────

def prune_tree(args, species_set: set, outdir: str):
    print("[4/5] Pruning phylogenetic tree …")
    if not HAS_ETE3:
        print("    WARNING: ete3 not installed. Skipping tree pruning.")
        print("    Install with: pip install ete3")
        print("    Or: conda install -c etetoolkit ete3")
        return

    t = Tree(args.tree, format=1)

    # strip .aa suffix from all leaf names; handle CH_Batde5 (no suffix)
    for leaf in t.get_leaves():
        leaf.name = leaf.name.replace(".aa", "").replace(".faa", "")

    # prune to included species
    leaves_in_tree = {leaf.name for leaf in t.get_leaves()}
    to_keep = list(species_set & leaves_in_tree)
    missing_from_tree = species_set - leaves_in_tree
    if missing_from_tree:
        print(f"    WARNING: {len(missing_from_tree)} species not in tree: {sorted(missing_from_tree)}")

    t.prune(to_keep, preserve_branch_length=True)
    out_path = os.path.join(outdir, "tree_pruned.nwk")
    t.write(outfile=out_path, format=1)
    print(f"    Pruned tree with {len(to_keep)} leaves → {out_path}")


# ── 5. WRITE OUTPUTS ────────────────────────────────────────────────────────

def write_outputs(pheno: pd.DataFrame, og: pd.DataFrame, ipr: pd.DataFrame, outdir: str):
    print("[5/5] Writing output files …")
    os.makedirs(outdir, exist_ok=True)

    pheno_path = os.path.join(outdir, "phenotypes.tsv")
    pheno.to_csv(pheno_path, sep="\t")
    print(f"    {pheno_path}")

    # align OG and IPR rows to phenotype species
    species = list(pheno.index)

    og_aligned = og.reindex(species)
    og_path = os.path.join(outdir, "X_OG.tsv")
    og_aligned.to_csv(og_path, sep="\t")
    print(f"    {og_path}  ({og_aligned.shape[1]} OGs)")

    ipr_aligned = ipr.reindex(species)
    ipr_path = os.path.join(outdir, "X_IPR.tsv")
    ipr_aligned.to_csv(ipr_path, sep="\t")
    print(f"    {ipr_path}  ({ipr_aligned.shape[1]} IPR terms)")


# ── MAIN ────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    pheno = load_phenotypes(args)
    species_set = set(pheno.index)
    min_count = max(1, round(args.min_freq * len(species_set)))
    print(f"    Frequency filter: retain features present in ≥{min_count} species "
          f"({args.min_freq*100:.0f}% of {len(species_set)})")

    og = load_og_matrix(args, species_set, min_count)
    ipr = load_ipr_matrix(args, species_set, min_count)

    os.makedirs(args.outdir, exist_ok=True)   # ensure outdir exists before tree write
    prune_tree(args, species_set, args.outdir)
    write_outputs(pheno, og, ipr, args.outdir)

    print("\nDone. Verify:")
    print(f"  Phenotypes: {len(pheno)} species")
    og_path = os.path.join(args.outdir, "X_OG.tsv")
    ipr_path = os.path.join(args.outdir, "X_IPR.tsv")
    print(f"  OG matrix:  see {og_path}")
    print(f"  IPR matrix: see {ipr_path}")
    print(f"  Tree:       see {os.path.join(args.outdir, 'tree_pruned.nwk')}")
    print("\nNext step: sbatch TreeWAS/scripts/01_submit_treewas.sh")


if __name__ == "__main__":
    main()
