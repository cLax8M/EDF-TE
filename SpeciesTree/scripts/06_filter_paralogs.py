#!/usr/bin/env python3
# =============================================================================
# 06_filter_paralogs.py — Detect and remove ancient paralogs from gene trees
#
# PURPOSE : Inspects gene trees from Step 05 to identify HOGs that contain
#           ancient paralogs (where gene duplication predates the speciation
#           events of interest). Two criteria are used:
#
#   1. BRANCH LENGTH Z-SCORE (primary filter)
#      Branch lengths that deviate > PARALOG_ZSCORE (default 3.5) standard
#      deviations from the tree's mean branch length are flagged as potentially
#      paralogous. This detects ancient, divergent gene copies that inflate
#      distances.
#      Ref: Merényi et al. 2023 (Nat. Ecol. Evol.)
#
#   2. MULTI-COPY RESOLUTION (secondary)
#      For HOGs with exactly 2 copies per species (allowed in Step 01):
#      retain only the sequence with the branch length closest to the gene
#      tree's median, discarding the divergent copy. This produces a pruned
#      FASTA with exactly 1 sequence per species.
#
# OUTCOME:
#   - HOGs with >PARALOG_ZSCORE_THRESHOLD anomalous branches WHERE position
#     suggests true ancient paralogy → entire HOG discarded
#   - HOGs with localized long branches (terminal duplications) → pruned to
#     one sequence per species → retained
#
# OUTPUT  :
#   06_final_genes/alignments/<HOG>.final.fa  — one sequence per species
#   06_final_genes/gene_list.txt              — final HOG list
#   06_final_genes/paralog_report.tsv         — per-HOG paralog analysis
#   06_final_genes/summary.txt                — final gene count summary
#
# USAGE   : python 06_filter_paralogs.py [--gene-trees-dir DIR] [--filtered-dir DIR]
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# REFS    : Merényi et al. 2023, Liu et al. 2024
# =============================================================================

import os
import sys
import argparse
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from Bio import Phylo

# ─── DEFAULT CONFIGURATION ────────────────────────────────────────────────────
WORK_DIR        = "/data/carlos/EDF-TE/SpeciesTree"
GENE_TREES_DIR  = f"{WORK_DIR}/05_gene_trees"
FILTERED_ALN_DIR = f"{WORK_DIR}/04_filtered/alignments"
OUT_DIR         = f"{WORK_DIR}/06_final_genes"
PARALOG_ZSCORE  = 3.5    # branch length z-score threshold (Merényi et al. 2023)


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter ancient paralogs from gene trees")
    p.add_argument("--gene-trees-dir",  default=GENE_TREES_DIR)
    p.add_argument("--filtered-aln-dir", default=FILTERED_ALN_DIR)
    p.add_argument("--out-dir",         default=OUT_DIR)
    p.add_argument("--zscore",          type=float, default=PARALOG_ZSCORE,
                   help="Branch length z-score threshold for paralog detection")
    return p.parse_args()


def get_branch_lengths(tree) -> np.ndarray:
    """Return array of all branch lengths in the tree (excluding root).
    Bio.Phylo: branch_length is None for the root clade.
    """
    bls = [c.branch_length for c in tree.find_clades()
           if c.branch_length is not None and c.branch_length > 0]
    return np.array(bls)


def detect_anomalous_branches(tree, zscore_threshold: float) -> list:
    """
    Detect leaves (sequences) connected by anomalously long branches.
    Returns list of (leaf_name, branch_length, z_score) for anomalous leaves.
    Bio.Phylo: tree.get_terminals() for leaves.
    """
    bls = get_branch_lengths(tree)
    if len(bls) < 3:
        return []

    mean_bl = np.mean(bls)
    std_bl  = np.std(bls)

    if std_bl == 0:
        return []

    anomalous = []
    for leaf in tree.get_terminals():
        bl = leaf.branch_length or 0.0
        z = (bl - mean_bl) / std_bl
        if z > zscore_threshold:
            anomalous.append((leaf.name, bl, z))

    return anomalous


def get_species_from_seq_id(seq_id: str) -> str:
    """
    Extract species name from sequence header.
    Format: SPECIES_NAME|gene_id or similar.
    """
    return seq_id.split("|")[0] if "|" in seq_id else seq_id


def prune_to_one_per_species(tree, aln_records: dict,
                             aln_records_sanitized: dict = None) -> tuple:
    """
    For multi-copy HOGs: keep one sequence per species (closest to median bl).
    Returns (pruned_sequences dict, action_log list).
    Bio.Phylo: tree.get_terminals() and leaf.branch_length.
    aln_records_sanitized: fallback lookup dict where keys have '%' replaced by '_',
    matching the IQTree newick sanitization (%252E -> _252E etc.).
    """
    by_species = defaultdict(list)
    for leaf in tree.get_terminals():
        sp = get_species_from_seq_id(leaf.name)
        by_species[sp].append((leaf.name, leaf.branch_length or 0.0))

    pruned_seqs = {}
    action_log  = []

    bls = get_branch_lengths(tree)
    median_bl = float(np.median(bls)) if len(bls) > 0 else 0.0

    for sp, copies in by_species.items():
        if len(copies) == 1:
            seq_id = copies[0][0]
        else:
            # Choose copy with branch length closest to tree median
            seq_id = min(copies, key=lambda x: abs(x[1] - median_bl))[0]
            discarded = [c[0] for c in copies if c[0] != seq_id]
            action_log.append({
                "species": sp, "kept": seq_id,
                "discarded": "; ".join(discarded),
                "reason": "multi_copy_pruned_to_median"
            })

        rec = aln_records.get(seq_id)
        if rec is None and aln_records_sanitized is not None:
            # IQTree replaces '%' with '_' in newick leaf names;
            # seq_id from the tree has '_252E' while the FASTA key has '%252E'.
            rec = aln_records_sanitized.get(seq_id)
        if rec is not None:
            pruned_seqs[sp] = rec

    return pruned_seqs, action_log


def analyze_tree_for_paralogs(tree_path: str, aln_path: str,
                               zscore_threshold: float) -> dict:
    """
    Main analysis function for a single HOG.
    Returns a result dictionary with status and (if kept) pruned sequences.
    """
    hog_id = Path(tree_path).stem

    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        return {"hog_id": hog_id, "status": "TREE_PARSE_ERROR",
                "reason": str(e), "n_seq_final": 0, "records": None}

    # Load alignment sequences
    try:
        aln_records = {rec.id: rec for rec in SeqIO.parse(aln_path, "fasta")}
        # IQTree sanitizes sequence names by replacing '%' with '_' in newick output.
        # Build a fallback lookup keyed by the sanitized (IQTree) name so sequences
        # with URL-encoded characters (e.g. %252E) are still found.
        aln_records_sanitized = {rid.replace('%', '_'): rec
                                  for rid, rec in aln_records.items()}
    except Exception as e:
        return {"hog_id": hog_id, "status": "ALN_READ_ERROR",
                "reason": str(e), "n_seq_final": 0, "records": None}

    # Detect anomalous branches
    anomalous = detect_anomalous_branches(tree, zscore_threshold)
    bls = get_branch_lengths(tree)
    mean_bl = float(np.mean(bls)) if len(bls) > 0 else 0.0
    max_bl  = float(np.max(bls))  if len(bls) > 0 else 0.0

    # Decision logic:
    # If >20% of leaves are anomalous → ancient paralogy → REJECT whole HOG
    n_leaves = len(tree.get_terminals())
    n_anomalous = len(anomalous)
    anomaly_fraction = n_anomalous / n_leaves if n_leaves > 0 else 0

    if anomaly_fraction > 0.20:
        return {
            "hog_id": hog_id,
            "status": "REJECTED_ANCIENT_PARALOG",
            "reason": f"{n_anomalous}/{n_leaves} leaves ({anomaly_fraction:.0%}) "
                      f"have z-score > {zscore_threshold}",
            "n_anomalous": n_anomalous, "n_leaves": n_leaves,
            "mean_bl": mean_bl, "max_bl": max_bl,
            "n_seq_final": 0, "records": None
        }

    # Prune to one sequence per species
    pruned_seqs, action_log = prune_to_one_per_species(tree, aln_records,
                                                       aln_records_sanitized)

    return {
        "hog_id": hog_id,
        "status": "PASS",
        "reason": "OK",
        "n_anomalous": n_anomalous, "n_leaves": n_leaves,
        "mean_bl": mean_bl, "max_bl": max_bl,
        "n_seq_final": len(pruned_seqs),
        "records": pruned_seqs,
        "pruning_log": action_log
    }


def main():
    args = parse_args()

    print("=" * 70)
    print(" STEP 06: Ancient paralog detection and removal")
    print("=" * 70)
    print(f" Gene trees dir : {args.gene_trees_dir}")
    print(f" Alignments dir : {args.filtered_aln_dir}")
    print(f" Output dir     : {args.out_dir}")
    print(f" Z-score cutoff : {args.zscore}  (Merényi et al. 2023)\n")

    # Find tree files
    tree_files = sorted(Path(args.gene_trees_dir).glob("*.treefile"))
    print(f"[06] Found {len(tree_files)} gene trees")

    if not tree_files:
        print("ERROR: No .treefile files found.")
        sys.exit(1)

    # Create output directories
    aln_out_dir = os.path.join(args.out_dir, "alignments")
    os.makedirs(aln_out_dir, exist_ok=True)

    # Process each HOG
    results = []
    passed  = 0
    rejected = 0

    for i, tree_path in enumerate(tree_files):
        hog_id = tree_path.stem

        # Find corresponding alignment
        aln_candidates = list(Path(args.filtered_aln_dir).glob(f"{hog_id}*.fa"))
        if not aln_candidates:
            results.append({"hog_id": hog_id, "status": "ALN_NOT_FOUND",
                             "reason": "alignment not found in filtered dir",
                             "n_seq_final": 0})
            rejected += 1
            continue

        aln_path = str(aln_candidates[0])
        result   = analyze_tree_for_paralogs(str(tree_path), aln_path, args.zscore)
        results.append(result)

        if result["status"] == "PASS" and result["records"]:
            # Write pruned FASTA (one sequence per species)
            out_fa = os.path.join(aln_out_dir, f"{hog_id}.final.fa")
            records_to_write = list(result["records"].values())
            SeqIO.write(records_to_write, out_fa, "fasta")
            passed += 1
        else:
            rejected += 1

        if (i + 1) % 50 == 0:
            print(f"     Processed {i+1}/{len(tree_files)} "
                  f"({passed} passed, {rejected} rejected)...")

    # Write results table
    report_df = pd.DataFrame([
        {k: v for k, v in r.items() if k not in ("records", "pruning_log")}
        for r in results
    ])
    report_path = os.path.join(args.out_dir, "paralog_report.tsv")
    report_df.to_csv(report_path, sep="\t", index=False, float_format="%.4f")

    # Write final gene list
    gene_list = [r["hog_id"] for r in results if r["status"] == "PASS"]
    gene_list_path = os.path.join(args.out_dir, "gene_list.txt")
    with open(gene_list_path, "w") as fh:
        fh.write("\n".join(gene_list) + "\n")

    # Summary
    print(f"\n{'=' * 70}")
    print(f" STEP 06 RESULTS — Paralog Filtering Summary")
    print(f"{'=' * 70}")
    print(f"  Gene trees analyzed    : {len(tree_files)}")
    print(f"  Passed (no paralogs)   : {passed}")
    print(f"  Rejected (ancient par) : {rejected}")
    status_counts = report_df["status"].value_counts()
    for status, count in status_counts.items():
        print(f"    {status:35s}: {count}")
    print(f"\n  Output:")
    print(f"    Final alignments : {aln_out_dir}")
    print(f"    Gene list        : {gene_list_path}")
    print(f"    Report           : {report_path}")
    print(f"\n  Next step: python 07_concatenate.py")
    print(f"{'=' * 70}")

    # Summary file
    summary_path = os.path.join(args.out_dir, "summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("STEP 06 — Paralog Filtering — Summary\n")
        fh.write("=" * 40 + "\n")
        fh.write(f"Gene trees analyzed    : {len(tree_files)}\n")
        fh.write(f"Passed                 : {passed}\n")
        fh.write(f"Rejected               : {rejected}\n")
        fh.write(f"Final gene list        : {gene_list_path}\n")
        fh.write(f"\nMethod: Branch length z-score threshold = {args.zscore}\n")
        fh.write(f"Ref: Merényi et al. 2023 (Nat. Ecol. Evol.)\n")


if __name__ == "__main__":
    main()
