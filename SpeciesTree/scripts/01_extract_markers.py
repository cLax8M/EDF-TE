#!/usr/bin/env python3
# =============================================================================
# 01_extract_markers.py — Extract phylogenetic markers from OrthoFinder orthogroups
#
# PURPOSE : Reads Orthogroups.GeneCount.tsv to filter orthogroups by occupancy
#           and copy-number, then extracts sequences from Orthogroup_Sequences/
#           with species-labelled headers (using Orthogroups.tsv for the mapping).
#
# DATA SOURCES:
#   - Orthogroups.GeneCount.tsv : gene counts per species per OG → fast filtering
#   - Orthogroups.tsv           : gene IDs per species per OG → sequence mapping
#   - Orthogroup_Sequences/*.fa : unaligned per-OG FASTA sequences
#
# NOTE: N1.tsv and other HOG files contain HOGs rooted at specific tree nodes
#       and have very low per-HOG species coverage. Orthogroups.tsv is the
#       correct file for broad-coverage marker selection.
#
# CRITERIA:
#   - Occupancy  : OG present in ≥ MIN_SPECIES (default 95 = 40% of 238)
#   - Copy number: ≤ MAX_COPIES genes per species in that OG
#
# OUTPUT  :
#   01_markers/<OG_ID>.fa        — one FASTA per passing OG (headers: SPECIES|gene_id)
#   01_markers/marker_stats.tsv  — per-OG statistics
#   01_markers/summary.txt       — overall extraction summary
#
# USAGE   : python 01_extract_markers.py
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# REFS    : Szánthó et al. 2025, Liu et al. 2024, Merényi et al. 2023
# =============================================================================

import os
import sys
import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ─── DEFAULT PATHS ─────────────────────────────────────────────────────────────
ORTHO_DIR      = "/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run"
OG_COUNTS_FILE = f"{ORTHO_DIR}/Orthogroups/Orthogroups.GeneCount.tsv"
OG_TABLE_FILE  = f"{ORTHO_DIR}/Orthogroups/Orthogroups.tsv"
OG_SEQ_DIR     = f"{ORTHO_DIR}/Orthogroup_Sequences"
WORK_DIR       = "/data/carlos/EDF-TE/SpeciesTree"
OUT_DIR        = f"{WORK_DIR}/01_markers"
MIN_SPECIES    = 95   # ≥40% of 238 species
MAX_COPIES     = 2    # max gene copies per species per OG


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract phylogenetic markers from OrthoFinder orthogroups")
    p.add_argument("--counts-file",  default=OG_COUNTS_FILE,
                   help="Path to Orthogroups.GeneCount.tsv")
    p.add_argument("--table-file",   default=OG_TABLE_FILE,
                   help="Path to Orthogroups.tsv (gene IDs per species)")
    p.add_argument("--og-seq-dir",   default=OG_SEQ_DIR,
                   help="Path to Orthogroup_Sequences/")
    p.add_argument("--out-dir",      default=OUT_DIR,
                   help="Output directory for FASTAs")
    p.add_argument("--min-species",  type=int, default=MIN_SPECIES,
                   help="Minimum species occupancy")
    p.add_argument("--max-copies",   type=int, default=MAX_COPIES,
                   help="Max gene copies per species per OG")
    p.add_argument("--dry-run",      action="store_true",
                   help="Report OGs to extract without writing files")
    return p.parse_args()


def filter_orthogroups(counts_file: str, min_species: int,
                       max_copies: int) -> pd.DataFrame:
    """
    Load Orthogroups.GeneCount.tsv and filter by occupancy + copy number.
    Returns DataFrame with passing OG IDs and per-OG stats.
    """
    print(f"[01] Loading gene count table: {counts_file}")
    counts = pd.read_csv(counts_file, sep="\t", index_col=0)

    # Drop 'Total' column if present (OrthoFinder adds it)
    if "Total" in counts.columns:
        counts = counts.drop(columns=["Total"])

    n_total   = len(counts)
    n_species = counts.shape[1]
    print(f"     Orthogroups   : {n_total:,}")
    print(f"     Species       : {n_species}")

    # Count species with 1..max_copies genes (the "good" copy-number species).
    # An OG is selected if at least min_species taxa have 1-max_copies copies.
    # Species with 0 copies are absent; species with >max_copies are treated as
    # absent during FASTA extraction (extract_og_fasta enforces max_copies per sp).
    # This avoids silently dropping species that are single-copy in most OGs but
    # happen to share those OGs with a few high-copy species in other lineages.
    good_count = ((counts >= 1) & (counts <= max_copies)).sum(axis=1)
    max_cp     = counts.max(axis=1)

    mask_good  = good_count >= min_species
    passing    = counts[mask_good].copy()
    passing["n_species_present"] = good_count[mask_good]
    passing["max_copies"]        = max_cp[mask_good]

    n_old_strict = (((counts > 0).sum(axis=1) >= min_species) & (max_cp <= max_copies)).sum()
    print(f"\n     Filter: ≥ {min_species} species with 1-{max_copies} copies "
          f"({min_species/n_species:.0%} occupancy)")
    print(f"     Pass (relaxed) : {len(passing):,}  "
          f"[old strict filter would give: {n_old_strict:,}]")

    if len(passing) == 0:
        print("\n  ERROR: 0 OGs passed. Check input file and thresholds.")
        sys.exit(1)
    elif len(passing) < 200:
        print(f"\n  WARNING: Only {len(passing)} OGs. Consider relaxing MIN_SPECIES.")
    else:
        print(f"\n  INFO: {len(passing)} OGs — within expected range (400-700).")

    return passing


def build_gene_to_species_map(table_file: str, og_ids: list) -> dict:
    """
    Read Orthogroups.tsv for the passing OG_IDs only.
    Returns nested dict: {og_id: {gene_id: species_name}}
    Species column names are keys; gene IDs are comma-separated within each cell.
    """
    print(f"\n[01] Reading gene-to-species mapping for {len(og_ids)} OGs...")
    print(f"     Source: {table_file}")

    og_id_set = set(og_ids)
    gene_map  = {}   # {og_id: {gene_id: species_name}}

    # Read in chunks for memory efficiency (Orthogroups.tsv can be >50MB)
    chunk_size = 5000
    rows_read  = 0
    found      = 0

    for chunk in pd.read_csv(table_file, sep="\t", index_col=0,
                              chunksize=chunk_size, low_memory=False):
        chunk_hits = chunk[chunk.index.isin(og_id_set)]
        for og_id, row in chunk_hits.iterrows():
            species_map = {}
            for sp_col, cell in row.items():
                if pd.isna(cell) or str(cell).strip() == "":
                    continue
                for gene_id in str(cell).split(", "):
                    gene_id = gene_id.strip()
                    if gene_id:
                        species_map[gene_id] = sp_col
            gene_map[og_id] = species_map
            found += 1
        rows_read += len(chunk)
        if found >= len(og_ids):
            break  # all found, stop reading

    print(f"     OGs mapped: {found} / {len(og_ids)}")
    return gene_map


def extract_og_fasta(og_id: str, gene_to_species: dict,
                     og_seq_dir: str, out_dir: str,
                     max_copies: int) -> dict:
    """
    Extract and rename sequences for a single OG.
    Headers become: SPECIES_NAME|gene_id
    Only keeps up to max_copies sequences per species.
    """
    fa_path = os.path.join(og_seq_dir, f"{og_id}.fa")
    if not os.path.exists(fa_path):
        return {"og_id": og_id, "status": "FASTA_NOT_FOUND",
                "n_sequences": 0, "n_species": 0}

    records_to_write = []
    species_count    = {}   # species → number of seqs already added

    try:
        for rec in SeqIO.parse(fa_path, "fasta"):
            gene_id = rec.id
            if gene_id not in gene_to_species:
                continue
            sp = gene_to_species[gene_id]
            # Enforce max_copies per species
            if species_count.get(sp, 0) >= max_copies:
                continue
            new_rec = SeqRecord(
                rec.seq,
                id=f"{sp}|{gene_id}",
                description=""
            )
            records_to_write.append(new_rec)
            species_count[sp] = species_count.get(sp, 0) + 1
    except Exception as e:
        return {"og_id": og_id, "status": f"FASTA_READ_ERROR: {e}",
                "n_sequences": 0, "n_species": 0}

    if not records_to_write:
        return {"og_id": og_id, "status": "NO_SEQUENCES_EXTRACTED",
                "n_sequences": 0, "n_species": 0}

    out_path = os.path.join(out_dir, f"{og_id}.fa")
    SeqIO.write(records_to_write, out_path, "fasta")

    return {"og_id": og_id, "status": "OK",
            "n_sequences": len(records_to_write),
            "n_species": len(species_count)}


def main():
    args = parse_args()

    print("=" * 70)
    print(" STEP 01: Extract phylogenetic markers from OrthoFinder")
    print("=" * 70)
    print(f" Counts file : {args.counts_file}")
    print(f" Table file  : {args.table_file}")
    print(f" OG seq dir  : {args.og_seq_dir}")
    print(f" Output dir  : {args.out_dir}")
    print(f" Min species : {args.min_species} (≥40% of 238)")
    print(f" Max copies  : {args.max_copies}")
    if args.dry_run:
        print(" DRY RUN     : will not write files")
    print()

    # Step 1: Filter orthogroups by occupancy and copy number
    passing_df = filter_orthogroups(
        args.counts_file, args.min_species, args.max_copies)
    og_ids = list(passing_df.index)

    if args.dry_run:
        passing_df[["n_species_present", "max_copies"]].to_csv(
            "dry_run_markers.tsv", sep="\t")
        print(f"\n[DRY RUN] {len(og_ids)} OGs would be extracted.")
        print(f"  List written to: dry_run_markers.tsv")
        sys.exit(0)

    # Step 2: Build gene-to-species mapping for passing OGs
    gene_map = build_gene_to_species_map(args.table_file, og_ids)

    # Step 3: Extract FASTAs
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"\n[01] Extracting FASTAs for {len(og_ids)} OGs...")
    stats_rows = []
    success = 0
    failed  = 0

    for i, og_id in enumerate(og_ids):
        gts = gene_map.get(og_id, {})
        stat = extract_og_fasta(
            og_id, gts, args.og_seq_dir, args.out_dir, args.max_copies)
        # Attach counts stats
        stat["n_species_filter"] = int(passing_df.loc[og_id, "n_species_present"])
        stat["max_copies_filter"] = int(passing_df.loc[og_id, "max_copies"])
        stats_rows.append(stat)

        if stat["status"] == "OK":
            success += 1
        else:
            failed += 1

        if (i + 1) % 100 == 0:
            print(f"     Processed {i+1}/{len(og_ids)} OGs "
                  f"({success} OK, {failed} failed)...")

    # Write stats
    stats_df   = pd.DataFrame(stats_rows)
    stats_path = os.path.join(args.out_dir, "marker_stats.tsv")
    stats_df.to_csv(stats_path, sep="\t", index=False)

    # Summary file
    summary_path = os.path.join(args.out_dir, "summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("STEP 01 — Extract Markers — Summary\n")
        fh.write("=" * 40 + "\n")
        fh.write(f"Orthogroups in GeneCount.tsv  : {len(passing_df):,} (after filter)\n")
        fh.write(f"FASTAs written successfully   : {success}\n")
        fh.write(f"FASTAs failed                 : {failed}\n")
        fh.write(f"Filters: min_species={args.min_species}, "
                 f"max_copies={args.max_copies}\n")
        fh.write(f"Output directory              : {args.out_dir}\n")
        fh.write(f"\nDETAILS: {stats_path}\n")

    print(f"\n{'=' * 70}")
    print(f" STEP 01 COMPLETE")
    print(f"   OGs extracted : {success}")
    print(f"   Failed        : {failed}")
    print(f"   Stats file    : {stats_path}")
    print(f"\n Next: sbatch 02_align_mafft_array.sh")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
