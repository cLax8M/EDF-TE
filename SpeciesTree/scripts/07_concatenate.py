#!/usr/bin/env python3
# =============================================================================
# 07_concatenate.py — Build supermatrix for species tree inference
#
# PURPOSE : Concatenates all single-copy (paralog-filtered) HOG alignments
#           from Step 06 into a single supermatrix FASTA file. Also generates
#           a partition file in IQTree format specifying the start and end
#           position of each gene in the concatenated alignment.
#
#           Missing sequences (species absent from a HOG) are represented as
#           gap-only rows (---...) to maintain a complete matrix with defined
#           missing data positions.
#
# OUTPUT  :
#   07_concatenated/supermatrix.fa          — concatenated alignment (FASTA)
#   07_concatenated/supermatrix.phy         — same in PHYLIP format
#   07_concatenated/partitions.txt          — IQTree partition file
#   07_concatenated/supermatrix_stats.tsv   — per-species completeness
#   07_concatenated/summary.txt             — overall supermatrix statistics
#
# USAGE   : python 07_concatenate.py [--final-genes-dir DIR] [--out-dir DIR]
#
# NOTES   :
#   - Species names are extracted from the first field before '|' in seq headers
#   - All 238 species are included even if absent from some genes (gaps added)
#   - Completeness (non-gap fraction) is reported per species and overall
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# =============================================================================

import os
import sys
import argparse
from pathlib import Path
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ─── DEFAULT CONFIGURATION ────────────────────────────────────────────────────
WORK_DIR         = "/data/carlos/EDF-TE/SpeciesTree"
FINAL_GENES_DIR  = f"{WORK_DIR}/06_final_genes/alignments"
OUT_DIR          = f"{WORK_DIR}/07_concatenated"
# Master species list: GeneCount TSV header (all 238 species)
SPECIES_REF      = "/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups/Orthogroups.GeneCount.tsv"


def parse_args():
    p = argparse.ArgumentParser(
        description="Concatenate paralog-filtered HOG alignments into supermatrix")
    p.add_argument("--final-genes-dir", default=FINAL_GENES_DIR,
                   help="Directory with .final.fa alignments from Step 06")
    p.add_argument("--out-dir",         default=OUT_DIR)
    p.add_argument("--min-genes-per-species", type=int, default=10,
                   help="Species with fewer genes than this are flagged (not removed)")
    p.add_argument("--species-ref",     default=SPECIES_REF,
                   help="GeneCount TSV — header used as authoritative species list. "
                        "Species absent from all alignments get all-gap rows.")
    return p.parse_args()


def load_reference_species(ref_tsv: str) -> set:
    """
    Read the header of Orthogroups.GeneCount.tsv and return the set of all
    species names (column headers, stripped of '.aa' suffix for consistency
    with sequence header format used in the alignments).
    """
    with open(ref_tsv) as fh:
        header = fh.readline().rstrip('\n').split('\t')
    # First column is 'Orthogroup', drop it and 'Total' if present
    species = set()
    for col in header[1:]:
        if col == 'Total':
            continue
        # Strip '.aa' suffix to match seq header format (e.g. MU_Jimfl.aa → MU_Jimfl.aa)
        # Actually keep as-is; seq headers use the full column name (e.g. MU_Jimfl.aa)
        species.add(col)
    return species


def collect_all_species(aln_files: list, ref_species: set = None) -> set:
    """Collect the union of all species across all alignments.
    If ref_species is provided, it is used as the master list (any species
    in ref_species but absent from all alignments will get all-gap rows).
    """
    all_species = set()
    for aln_path in aln_files:
        for rec in SeqIO.parse(aln_path, "fasta"):
            sp = rec.id.split("|")[0]
            all_species.add(sp)
    if ref_species:
        missing = ref_species - all_species
        if missing:
            print(f"     WARNING: {len(missing)} species absent from all alignments "
                  f"— added as all-gap rows: {', '.join(sorted(missing))}")
        all_species = all_species | ref_species
    return all_species


def concatenate_alignments(aln_files: list, all_species: set) -> tuple:
    """
    Build supermatrix as dict: {species: concatenated_sequence_string}.
    Also builds partition data: list of (hog_id, start, end).
    """
    # Initialize sequences for all species
    supermatrix = {sp: [] for sp in all_species}
    partitions  = []   # list of (hog_id, start_1based, end_1based)
    gene_lengths = {}  # hog_id → alignment length in AA columns

    current_pos = 1    # 1-based

    for aln_path in aln_files:
        hog_id = Path(aln_path).stem.replace(".final", "")
        records = {rec.id.split("|")[0]: str(rec.seq)
                   for rec in SeqIO.parse(aln_path, "fasta")}

        # Get alignment length from first sequence
        if not records:
            continue
        aln_len = len(next(iter(records.values())))

        # Add sequence (or gaps) for each species
        for sp in all_species:
            if sp in records:
                supermatrix[sp].append(records[sp])
            else:
                supermatrix[sp].append("-" * aln_len)

        gene_lengths[hog_id] = aln_len
        partitions.append((hog_id, current_pos, current_pos + aln_len - 1))
        current_pos += aln_len

    # Join all fragments
    supermatrix_joined = {sp: "".join(frags) for sp, frags in supermatrix.items()}
    return supermatrix_joined, partitions, gene_lengths


def compute_completeness(supermatrix: dict, partitions: list) -> pd.DataFrame:
    """Compute per-species completeness (fraction of non-gap positions)."""
    total_len = len(next(iter(supermatrix.values())))
    rows = []
    for sp, seq in supermatrix.items():
        n_total = len(seq)
        n_gaps  = seq.count("-")
        n_data  = n_total - n_gaps
        frac    = n_data / n_total if n_total > 0 else 0.0
        rows.append({
            "species": sp,
            "total_positions": n_total,
            "data_positions": n_data,
            "gap_positions": n_gaps,
            "completeness": frac
        })
    df = pd.DataFrame(rows).sort_values("completeness", ascending=False)
    return df


def write_fasta(supermatrix: dict, out_path: str):
    """Write supermatrix as FASTA."""
    records = [
        SeqRecord(Seq(seq), id=sp, description="")
        for sp, seq in sorted(supermatrix.items())
    ]
    SeqIO.write(records, out_path, "fasta")


def write_phylip(supermatrix: dict, out_path: str):
    """Write supermatrix in PHYLIP format."""
    n_taxa = len(supermatrix)
    aln_len = len(next(iter(supermatrix.values())))
    with open(out_path, "w") as fh:
        fh.write(f" {n_taxa} {aln_len}\n")
        for sp in sorted(supermatrix.keys()):
            # PHYLIP requires name ≤ 10 chars padded — use full name with space
            fh.write(f"{sp:<30}{supermatrix[sp]}\n")


def write_partitions(partitions: list, out_path: str):
    """
    Write IQTree partition file (-spp / nexus format).
    Format: LG, GENE_NAME = start-end
    """
    with open(out_path, "w") as fh:
        fh.write("#nexus\n")
        fh.write("begin sets;\n")
        for hog_id, start, end in partitions:
            fh.write(f"  charset {hog_id} = {start}-{end};\n")
        fh.write("end;\n")


def main():
    args = parse_args()

    print("=" * 70)
    print(" STEP 07: Build supermatrix for species tree inference")
    print("=" * 70)
    print(f" Input dir   : {args.final_genes_dir}")
    print(f" Output dir  : {args.out_dir}\n")

    # Find alignment files
    aln_files = sorted(Path(args.final_genes_dir).glob("*.final.fa"))
    print(f"[07] Found {len(aln_files)} final HOG alignments")

    if len(aln_files) == 0:
        print("ERROR: No .final.fa files found.")
        sys.exit(1)

    if len(aln_files) < 100:
        print(f"WARNING: Only {len(aln_files)} genes. Consider relaxing earlier filters.")
    elif len(aln_files) < 225:
        print(f"WARNING: {len(aln_files)} genes (below typical target of 225-300).")
    else:
        print(f"INFO: {len(aln_files)} genes — within target range (225-300+).")

    # Collect all species (use GeneCount header as reference to catch any missing)
    print(f"\n[07] Collecting species list from all alignments...")
    ref_species = None
    if args.species_ref and os.path.exists(args.species_ref):
        ref_species = load_reference_species(args.species_ref)
        print(f"     Reference species list: {len(ref_species)} taxa from GeneCount TSV")
    all_species = collect_all_species([str(f) for f in aln_files], ref_species)
    print(f"     Total species in supermatrix: {len(all_species)}")

    # Build supermatrix
    print(f"\n[07] Concatenating {len(aln_files)} alignments...")
    supermatrix, partitions, gene_lengths = concatenate_alignments(
        [str(f) for f in aln_files], all_species)

    total_len = len(next(iter(supermatrix.values())))
    print(f"     Supermatrix dimensions: {len(all_species)} taxa × {total_len:,} AA columns")

    # Compute completeness
    completeness_df = compute_completeness(supermatrix, partitions)
    mean_completeness = completeness_df["completeness"].mean()
    print(f"     Mean completeness     : {mean_completeness:.1%}")
    print(f"     Min completeness      : {completeness_df['completeness'].min():.1%} "
          f"({completeness_df.iloc[-1]['species']})")

    # Flag flagged species
    flagged = completeness_df[
        completeness_df["data_positions"] < args.min_genes_per_species * 50]
    if len(flagged) > 0:
        print(f"\n  WARNING: {len(flagged)} species with very low data:")
        for _, row in flagged.head(10).iterrows():
            print(f"    {row['species']}: {row['completeness']:.1%}")

    # Write outputs
    os.makedirs(args.out_dir, exist_ok=True)

    fasta_path      = os.path.join(args.out_dir, "supermatrix.fa")
    phylip_path     = os.path.join(args.out_dir, "supermatrix.phy")
    partition_path  = os.path.join(args.out_dir, "partitions.txt")
    stats_path      = os.path.join(args.out_dir, "supermatrix_stats.tsv")
    summary_path    = os.path.join(args.out_dir, "summary.txt")

    print(f"\n[07] Writing outputs...")
    write_fasta(supermatrix, fasta_path)
    write_phylip(supermatrix, phylip_path)
    write_partitions(partitions, partition_path)
    completeness_df.to_csv(stats_path, sep="\t", index=False, float_format="%.4f")

    # Summary
    with open(summary_path, "w") as fh:
        fh.write("STEP 07 — Supermatrix — Summary\n")
        fh.write("=" * 40 + "\n")
        fh.write(f"Number of genes     : {len(aln_files)}\n")
        fh.write(f"Number of taxa      : {len(all_species)}\n")
        fh.write(f"Alignment length    : {total_len:,} AA columns\n")
        fh.write(f"Mean completeness   : {mean_completeness:.1%}\n")
        fh.write(f"\nFiles:\n")
        fh.write(f"  FASTA       : {fasta_path}\n")
        fh.write(f"  PHYLIP      : {phylip_path}\n")
        fh.write(f"  Partitions  : {partition_path}\n")
        fh.write(f"  Stats       : {stats_path}\n")

    print(f"\n{'=' * 70}")
    print(f" STEP 07 COMPLETE")
    print(f"   Genes            : {len(aln_files)}")
    print(f"   Taxa             : {len(all_species)}")
    print(f"   Alignment length : {total_len:,} AA columns")
    print(f"   Mean completeness: {mean_completeness:.1%}")
    print(f"   FASTA            : {fasta_path}")
    print(f"   Partitions       : {partition_path}")
    print(f"\n   Next step: sbatch 08_species_tree_guide.sh")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
