#!/usr/bin/env python3
"""
01_gff_to_bed.py — Convert GFF/GFF3 annotation files to 4-column BED and
                   genome-size files required for the synteny analysis.

Output per species:
  {BED_DIR}/{species}.bed        — chrom  start  end  gene_id  (0-based, sorted)
  {BED_DIR}/{species}.genome     — chrom  size   (for bedtools shuffle)

Usage:
  python 01_gff_to_bed.py --gff_dir GFF_DIR --out_dir BED_DIR [--species SP1 SP2 ...]
  
  If --species is omitted, all *.gff and *.gff3 files in GFF_DIR are processed.

Author: EDF-TE project
Last updated: 2026-02-20
"""

import argparse
import os
import sys
import re
from pathlib import Path
from collections import defaultdict


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff_dir", required=True,
                   help="Directory containing {CLADE}_{Species}.gff[3] files")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for .bed and .genome files")
    p.add_argument("--species", nargs="*",
                   help="Species prefixes to process (e.g. HO_Naegrub). "
                        "Default: all files in gff_dir")
    p.add_argument("--min_scaffold_genes", type=int, default=3,
                   help="Minimum genes on a scaffold to include (default: 3)")
    p.add_argument("--dry_run", action="store_true",
                   help="Print species list without processing")
    return p.parse_args()


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def parse_attributes(attr_str):
    """Parse GFF3 attribute string into a dict."""
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if "=" in field:
            k, v = field.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def extract_gene_id(attrs, source_line):
    """
    Extract a clean gene identifier from GFF attributes.
    Priority: locus_tag → Name → ID (stripped of 'gene-' prefix).
    """
    if "locus_tag" in attrs:
        return attrs["locus_tag"]
    if "Name" in attrs:
        return attrs["Name"]
    if "ID" in attrs:
        gid = attrs["ID"]
        # Strip common prefixes from NCBI GFF3
        gid = re.sub(r"^gene-", "", gid)
        gid = re.sub(r"^gene:", "", gid)
        return gid
    return None


def gff_to_bed(gff_path, min_scaffold_genes=3):
    """
    Parse a GFF/GFF3 file and return:
      - records: list of (chrom, start, end, gene_id)  [0-based BED coords]
      - genome_sizes: dict  chrom -> max_end seen (proxy for scaffold length)
      - skipped: number of 'gene' lines without a parseable ID

    Only 'gene' feature lines are extracted. CDS/exon lines are ignored.
    Scaffolds with fewer than min_scaffold_genes genes are filtered out.
    """
    genes_by_chrom = defaultdict(list)  # chrom -> [(start0, end, gene_id)]
    genome_sizes = {}
    skipped = 0

    # Try UTF-8 first, fall back to latin-1 for legacy GFF files with non-ASCII bytes
    for encoding in ("utf-8", "latin-1"):
        try:
            with open(gff_path, encoding=encoding) as fh:
                fh.read(1024)  # probe
            break
        except UnicodeDecodeError:
            continue
    else:
        encoding = "latin-1"  # last resort

    with open(gff_path, encoding=encoding, errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                # Capture ##sequence-region directives for genome sizes
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            genome_sizes[parts[1]] = int(parts[3])
                        except ValueError:
                            pass
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            feature = cols[2]
            if feature != "gene":
                continue

            chrom = cols[0]
            try:
                start1 = int(cols[3])  # GFF is 1-based inclusive
                end1   = int(cols[4])
            except ValueError:
                continue

            start0 = start1 - 1  # convert to 0-based
            attrs = parse_attributes(cols[8])
            gene_id = extract_gene_id(attrs, line)

            if gene_id is None:
                skipped += 1
                continue

            genes_by_chrom[chrom].append((start0, end1, gene_id))

            # Update genome size fallback from max coordinate seen
            if chrom not in genome_sizes or end1 > genome_sizes[chrom]:
                genome_sizes[chrom] = end1

    # Filter scaffolds with too few genes, then sort by position
    records = []
    filtered_chroms = set()
    for chrom, genes in genes_by_chrom.items():
        if len(genes) < min_scaffold_genes:
            filtered_chroms.add(chrom)
            continue
        genes_sorted = sorted(genes, key=lambda x: x[0])
        for start0, end, gene_id in genes_sorted:
            records.append((chrom, start0, end, gene_id))

    return records, genome_sizes, skipped, filtered_chroms


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def find_gff_files(gff_dir, species_list=None):
    """Return list of (species_prefix, path) tuples."""
    gff_dir = Path(gff_dir)
    results = []
    for ext in ("*.gff", "*.gff3"):
        for fp in sorted(gff_dir.glob(ext)):
            # Skip macOS AppleDouble hidden metadata files (._*)
            if fp.name.startswith("._"):
                continue
            # Species prefix = filename without extension
            prefix = fp.stem
            # Strip trailing _c (cleaned) suffix used by OrthoFinder inputs
            prefix_clean = re.sub(r"_c$", "", prefix)
            if species_list and prefix_clean not in species_list and prefix not in species_list:
                continue
            results.append((prefix_clean, fp))
    return results


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gff_files = find_gff_files(args.gff_dir, args.species)

    if not gff_files:
        print(f"ERROR: no GFF files found in {args.gff_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(gff_files)} GFF file(s) to process.")

    if args.dry_run:
        for prefix, fp in gff_files:
            print(f"  {prefix}  <-  {fp}")
        return

    summary_rows = []

    for prefix, gff_path in gff_files:
        print(f"  Processing {prefix} ...", end=" ", flush=True)

        records, genome_sizes, skipped, filtered_chroms = gff_to_bed(
            gff_path, min_scaffold_genes=args.min_scaffold_genes
        )

        if not records:
            print(f"WARNING: no usable genes — skipped.")
            continue

        # Write BED file
        bed_out = out_dir / f"{prefix}.bed"
        with open(bed_out, "w") as fh:
            for chrom, start0, end, gene_id in records:
                fh.write(f"{chrom}\t{start0}\t{end}\t{gene_id}\n")

        # Write genome-size file (for bedtools shuffle)
        # Include ALL scaffolds (even small ones) so shuffle can place genes anywhere
        genome_out = out_dir / f"{prefix}.genome"
        with open(genome_out, "w") as fh:
            for chrom, size in sorted(genome_sizes.items()):
                fh.write(f"{chrom}\t{size}\n")

        n_genes = len(records)
        n_chroms = len(set(r[0] for r in records))
        print(f"{n_genes} genes on {n_chroms} scaffolds "
              f"(filtered {len(filtered_chroms)} short scaffolds, {skipped} IDs skipped)")

        summary_rows.append({
            "species": prefix,
            "n_genes": n_genes,
            "n_scaffolds": n_chroms,
            "n_filtered_scaffolds": len(filtered_chroms),
            "n_skipped_ids": skipped
        })

    # Write summary TSV
    summary_path = out_dir / "01_gff_to_bed_summary.tsv"
    with open(summary_path, "w") as fh:
        fh.write("species\tn_genes\tn_scaffolds\tn_filtered_scaffolds\tn_skipped_ids\n")
        for row in summary_rows:
            fh.write(
                f"{row['species']}\t{row['n_genes']}\t{row['n_scaffolds']}\t"
                f"{row['n_filtered_scaffolds']}\t{row['n_skipped_ids']}\n"
            )
    print(f"\nSummary written to {summary_path}")
    print("Step 01 complete.")


if __name__ == "__main__":
    main()
