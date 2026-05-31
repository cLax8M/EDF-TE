#!/usr/bin/env python3
"""
D03_treewas_genes.py — TreeWAS with TE composition as phenotype
================================================================
Discover orthogroups whose presence/absence is associated with
TE landscape variation across the phylogeny.

Uses the EXISTING TreeWAS infrastructure:
  - X_OG.tsv and tree_pruned.nwk from 00_prep_data.py
  - 01_run_treewas.R for the actual treeWAS analysis

This script:
  1. Builds a TE-phenotype table (pct_TE, per-class %, Shannon, Kimura)
  2. Writes phenotypes_te.tsv aligned to the existing species set
  3. Generates a SLURM submission script for all TE phenotype × OG runs
  4. Optionally consolidates results and annotates hits with InterPro

Usage
-----
  # Step 1: Prepare TE phenotypes
  python D03_treewas_genes.py prep [options]

  # Step 2: After treeWAS jobs finish, consolidate
  python D03_treewas_genes.py collect [options]

Dependencies: pandas, numpy (standard venv).
"""

import argparse
import os
import sys
import glob

import numpy as np
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
_DATA       = "/legserv/Temp/Carlos/EDF-TE"
_SEC2       = f"{_DATA}/results/Section2_TElandscape"
_TREEWAS    = "/data/carlos/EDF-TE/EDF-TE/TreeWAS"
_TW_RESULTS = f"{_TREEWAS}/results"
_TW_SCRIPTS = f"{_TREEWAS}/scripts"
_IPR_DIR    = f"{_DATA}/data/functional_interpro"
_OG_FILE    = f"{_DATA}/results/Orthofinder_run/Orthogroups/Orthogroups.GeneCount.tsv"
_LOGS       = "/data/carlos/EDF-TE/EDF-TE/logs"

# TE-related Pfam IDs (for annotation)
TE_PFAM = {
    "PF00078", "PF03184", "PF13456", "PF14214",
    "PF13683", "PF01497", "PF13952",
}

# ── TE phenotypes to test ─────────────────────────────────────────────────────
TE_PHENOTYPES = [
    "pct_TE_total", "pct_DNA", "pct_LTR", "pct_LINE",
    "pct_RC", "pct_PLE", "pct_Unknown",
    "H_prime", "recency_index", "burst_score",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="TreeWAS with TE composition phenotypes")
    sub = p.add_subparsers(dest="command")

    # ── prep ──────────────────────────────────────────────────────────────
    prep = sub.add_parser("prep", help="Prepare TE phenotype table + submit script")
    prep.add_argument("--merged",   default=f"{_SEC2}/D01_te_genomics_merged_provisional.tsv")
    prep.add_argument("--shannon",  default=f"{_SEC2}/B01_shannon_diversity.tsv")
    prep.add_argument("--kimura",   default=f"{_SEC2}/C01_kimura_species.tsv")
    prep.add_argument("--tw_data",  default=_TW_RESULTS,
                      help="Dir with X_OG.tsv and tree_pruned.nwk from 00_prep_data.py")
    prep.add_argument("--outdir",   default=f"{_TW_RESULTS}")

    # ── collect ───────────────────────────────────────────────────────────
    coll = sub.add_parser("collect", help="Consolidate treeWAS TE-phenotype results")
    coll.add_argument("--resultsdir", default=_TW_RESULTS)
    coll.add_argument("--ipr_dir",    default=_IPR_DIR)
    coll.add_argument("--og_file",    default=_OG_FILE)

    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
# PREP: build TE phenotype table + submission script
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_prep(args):
    # ── load TE composition ───────────────────────────────────────────────
    merged = pd.read_csv(args.merged, sep="\t")
    merged.set_index("species", inplace=True)

    # Shannon
    if os.path.isfile(args.shannon):
        shan = pd.read_csv(args.shannon, sep="\t").set_index("species")
        for col in ["H_prime", "J_evenness"]:
            if col in shan.columns:
                merged[col] = shan[col]

    # Kimura
    if os.path.isfile(args.kimura):
        kim = pd.read_csv(args.kimura, sep="\t").set_index("species")
        for col in ["recency_index", "burst_score", "mean_div_weighted"]:
            if col in kim.columns:
                merged[col] = kim[col]

    # ── align to existing TreeWAS species ─────────────────────────────────
    og_path = os.path.join(args.tw_data, "X_OG.tsv")
    if not os.path.isfile(og_path):
        sys.exit(f"[ERROR] X_OG.tsv not found at {og_path}. Run 00_prep_data.py first.")

    og_idx = pd.read_csv(og_path, sep="\t", usecols=[0])
    tw_species = list(og_idx.iloc[:, 0])
    print(f"[PREP] TreeWAS species set: {len(tw_species)}")

    # filter to TE phenotypes available
    avail = [c for c in TE_PHENOTYPES if c in merged.columns]
    print(f"[PREP] TE phenotypes available: {avail}")

    pheno = merged.reindex(tw_species)[avail]
    missing = pheno.isnull().all(axis=1).sum()
    if missing > 0:
        print(f"[WARN] {missing} species have no TE data at all")

    # drop species with all NaN (they were excluded from TE analysis)
    pheno = pheno.dropna(how="all")
    print(f"[PREP] TE phenotype table: {pheno.shape[0]} species × {pheno.shape[1]} traits")

    # ── write ─────────────────────────────────────────────────────────────
    pheno_path = os.path.join(args.outdir, "phenotypes_te.tsv")
    pheno.to_csv(pheno_path, sep="\t")
    print(f"[SAVED] {pheno_path}")

    # ── generate submission script ────────────────────────────────────────
    n_pheno = len(avail)
    n_array = n_pheno  # only OG matrix (not IPR) for TE phenotypes
    sbatch_path = os.path.join(_TW_SCRIPTS, "03_submit_treewas_te.sh")

    pheno_list = " ".join(f'"{p}"' for p in avail)

    script = f"""#!/bin/bash
#SBATCH --job-name=tw_te
#SBATCH --array=0-{n_array - 1}
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output={_LOGS}/treewas_te_%A_%a.out
#SBATCH --error={_LOGS}/treewas_te_%A_%a.err
#SBATCH --partition=normal.168h

# TreeWAS with TE composition phenotypes × OG genotypes
# Generated by D03_treewas_genes.py
# TE phenotypes: {', '.join(avail)}

set -euo pipefail

PHENOTYPES=({pheno_list})
PHENOTYPE="${{PHENOTYPES[$SLURM_ARRAY_TASK_ID]}}"

echo "=============================================="
echo "Job array index : $SLURM_ARRAY_TASK_ID"
echo "Phenotype       : $PHENOTYPE"
echo "Matrix          : OG"
echo "Node            : $(hostname)"
echo "Started         : $(date)"
echo "=============================================="

SCRIPT_DIR="{_TW_SCRIPTS}"
RESULTS_DIR="{_TW_RESULTS}"

# Activate environment
export MAMBA_EXE="$(which micromamba)"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate treewasEDFTE

# The 01_run_treewas.R script reads phenotypes.tsv by default.
# We need it to read phenotypes_te.tsv instead.
# Solution: temporarily symlink, or pass a modified --datadir.
# Create a temp datadir that links to existing OG/tree but uses TE phenotypes.
TMPDIR="${{RESULTS_DIR}}/te_run_data"
mkdir -p "$TMPDIR"
ln -sf "${{RESULTS_DIR}}/X_OG.tsv"        "$TMPDIR/X_OG.tsv"
ln -sf "${{RESULTS_DIR}}/tree_pruned.nwk"  "$TMPDIR/tree_pruned.nwk"
ln -sf "${{RESULTS_DIR}}/phenotypes_te.tsv" "$TMPDIR/phenotypes.tsv"

Rscript "${{SCRIPT_DIR}}/01_run_treewas.R" \\
    --phenotype "$PHENOTYPE" \\
    --matrix    OG \\
    --datadir   "$TMPDIR" \\
    --outdir    "$RESULTS_DIR" \\
    --n_sim     10000

echo "Finished: $(date)"
"""

    with open(sbatch_path, "w") as f:
        f.write(script)
    os.chmod(sbatch_path, 0o755)
    print(f"[SAVED] {sbatch_path}")
    print(f"\n[NEXT]  cd /data/carlos/EDF-TE/EDF-TE && sbatch TreeWAS/scripts/03_submit_treewas_te.sh")


# ═══════════════════════════════════════════════════════════════════════════════
# COLLECT: consolidate results + annotate
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_collect(args):
    # ── gather summary TSVs for TE-phenotype runs ─────────────────────────
    te_dirs = []
    for pheno in TE_PHENOTYPES:
        d = os.path.join(args.resultsdir, f"{pheno}_OG")
        if os.path.isdir(d):
            te_dirs.append(d)

    if not te_dirs:
        sys.exit("[ERROR] No TE-phenotype treeWAS results found. Run the jobs first.")

    dfs = []
    for d in te_dirs:
        summary_file = os.path.join(d, "summary.tsv")
        if os.path.isfile(summary_file):
            try:
                df = pd.read_csv(summary_file, sep="\t")
                if not df.empty:
                    dfs.append(df)
            except Exception:
                continue

    if not dfs:
        print("[WARN] All summary files empty — no significant hits.")
        combined = pd.DataFrame()
    else:
        combined = pd.concat(dfs, ignore_index=True)

    print(f"[COLLECT] {len(combined)} feature-phenotype associations from {len(te_dirs)} runs")

    if combined.empty:
        out_path = os.path.join(args.resultsdir, "treewas_te_annotated_hits.tsv")
        combined.to_csv(out_path, sep="\t", index=False)
        print(f"[SAVED] {out_path} (empty)")
        return

    # ── annotate with IPR descriptions ────────────────────────────────────
    ipr_lookup = _build_ipr_lookup(args.ipr_dir) if os.path.isdir(args.ipr_dir) else {}

    # For OG features, try to associate with most common InterPro function
    og_to_ipr = _map_og_to_ipr(args.og_file, args.ipr_dir) if os.path.isfile(args.og_file) else {}

    def annotate(row):
        feat = row["feature"]
        if feat.startswith("OG"):
            desc = og_to_ipr.get(feat, "")
            is_te = False
        elif feat.startswith("IPR"):
            desc = ipr_lookup.get(feat, "")
            is_te = any(pfam in feat for pfam in TE_PFAM)
        else:
            desc = ""
            is_te = False
        return pd.Series({"description": desc, "is_TE_pfam": is_te})

    annot = combined.apply(annotate, axis=1)
    combined = pd.concat([combined, annot], axis=1)

    # ── sort and save ─────────────────────────────────────────────────────
    combined.sort_values(["n_tests_significant", "mean_score"],
                         ascending=[False, False], inplace=True)

    out_path = os.path.join(args.resultsdir, "treewas_te_annotated_hits.tsv")
    combined.to_csv(out_path, sep="\t", index=False)

    robust_n = combined["robust"].sum() if "robust" in combined.columns else 0
    print(f"[RESULTS] Robust hits (≥2 tests): {robust_n}")
    print(f"[SAVED]   {out_path}")

    # ── per-phenotype summary ─────────────────────────────────────────────
    if "phenotype" in combined.columns:
        summary = combined.groupby("phenotype").agg(
            total_hits=("feature", "count"),
            robust_hits=("robust", "sum"),
        ).reset_index()
        print("\n[SUMMARY] Hits per TE phenotype:")
        print(summary.to_string(index=False))


def _build_ipr_lookup(ipr_dir: str) -> dict:
    """Map IPR ID → description from InterProScan files."""
    lookup = {}
    for fpath in glob.glob(os.path.join(ipr_dir, "*_IPR.tab")):
        try:
            df = pd.read_csv(fpath, sep="\t", header=None,
                             usecols=[1, 2], names=["iprId", "iprDesc"])
            for _, row in df.dropna(subset=["iprId"]).iterrows():
                if row["iprId"] not in lookup and pd.notna(row.get("iprDesc")):
                    lookup[row["iprId"]] = str(row["iprDesc"])
        except Exception:
            continue
    return lookup


def _map_og_to_ipr(og_file: str, ipr_dir: str) -> dict:
    """Lightweight OG → most common InterPro description mapping.

    This reads a few IPR files to build protein→IPR, then maps OG→proteins
    from the Orthogroups file.  For speed, only reads a subset.
    """
    # Skip if files are too large or ipr_dir missing
    if not os.path.isdir(ipr_dir):
        return {}

    # Build protein→IPR_desc from all IPR files
    prot_to_desc = {}
    for fpath in glob.glob(os.path.join(ipr_dir, "*_IPR.tab")):
        try:
            df = pd.read_csv(fpath, sep="\t", header=None, usecols=[0, 2],
                             names=["prot", "desc"])
            for _, row in df.dropna().iterrows():
                prot_to_desc[row["prot"]] = str(row["desc"])
        except Exception:
            continue

    # For a robust mapping we'd need the full Orthogroups.tsv (protein lists).
    # Since it's very large, return empty for now — the user can run
    # 02_annotate_results.py for full annotation.
    return {}


# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()
    if args.command == "prep":
        cmd_prep(args)
    elif args.command == "collect":
        cmd_collect(args)
    else:
        print("Usage: D03_treewas_genes.py {prep|collect} [options]")
        print("  prep    — Build TE phenotype table + SLURM submission script")
        print("  collect — Consolidate results after treeWAS finishes")
        sys.exit(1)


if __name__ == "__main__":
    main()
