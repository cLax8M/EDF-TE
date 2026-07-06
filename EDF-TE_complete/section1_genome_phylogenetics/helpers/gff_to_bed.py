#!/usr/bin/env python3
"""Minimal GFF/GFF3 -> BED converter for the EDF-TE synteny pipeline."""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gff", required=True, help="Input GFF/GFF3 file")
    p.add_argument("--out", required=True, help="Output BED file")
    p.add_argument("--min_scaffold_genes", type=int, default=3)
    p.add_argument("--dry_run", action="store_true")
    return p.parse_args()


def open_text(path):
    for enc in ("utf-8", "latin-1"):
        try:
            return open(path, encoding=enc, errors="replace")
        except UnicodeDecodeError:
            continue
    return open(path, encoding="latin-1", errors="replace")


def parse_gff(gff_path, min_scaffold_genes=3):
    raw_entries = defaultdict(list)
    genome_sizes = {}
    has_gene_features = False

    with open_text(gff_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 3 and cols[2] == "gene":
                    has_gene_features = True
                    break

    with open_text(gff_path) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        genome_sizes[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            chrom = cols[0]
            try:
                start0 = int(cols[3]) - 1
                end = int(cols[4])
            except ValueError:
                continue
            attrs = {}
            attr_str = cols[8]
            if "=" in attr_str:
                for field in attr_str.split(";"):
                    field = field.strip()
                    if "=" in field:
                        k, v = field.split("=", 1)
                        attrs[k.strip()] = v.strip()

            if has_gene_features and feature == "gene":
                keys = []
                for key in ("locus_tag", "ID", "Name", "gene", "proteinId"):
                    if key in attrs:
                        keys.append(attrs[key])
                if "Parent" in attrs:
                    keys.append(attrs["Parent"])
                if not keys:
                    keys = [f"gene_{len(raw_entries[chrom]) + 1}"]
                raw_entries[chrom].append((start0, end, keys))
            elif not has_gene_features and feature in {"CDS", "exon"}:
                key_candidates = []
                if "Name" in attrs:
                    key_candidates.append(attrs["Name"])
                if "gene" in attrs:
                    key_candidates.append(attrs["gene"])
                if "locus_tag" in attrs:
                    key_candidates.append(attrs["locus_tag"])
                if "proteinId" in attrs:
                    key_candidates.append(attrs["proteinId"])
                if key_candidates:
                    raw_entries[chrom].append((start0, end, key_candidates))
            elif has_gene_features and feature in {"CDS", "exon"}:
                # Use the parent gene if available
                parent = attrs.get("Parent") or attrs.get("ID")
                if parent:
                    raw_entries[chrom].append((start0, end, [parent]))

    gene_entries = defaultdict(list)
    for chrom, entries in raw_entries.items():
        if len(entries) < min_scaffold_genes:
            continue
        entries_sorted = sorted(entries, key=lambda x: x[0])
        for rank, (start0, end, keys) in enumerate(entries_sorted):
            gene_id = None
            for candidate in keys:
                if candidate:
                    gene_id = candidate
                    break
            if not gene_id:
                gene_id = f"{chrom}_g{rank+1}"
            gene_entries[chrom].append((start0, end, gene_id))

    return gene_entries, genome_sizes


def write_output(gene_entries, genome_sizes, out_path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    bed_path = out_path
    genome_path = out_path.with_suffix(".genome")
    summary_path = out_path.with_suffix(".id_map.tsv")

    bed_lines = []
    for chrom in sorted(gene_entries):
        entries = sorted(gene_entries[chrom], key=lambda x: x[0])
        for start0, end, gene_id in entries:
            bed_lines.append(f"{chrom}\t{start0}\t{end}\t{gene_id}\n")

    with open(bed_path, "w") as fh:
        fh.writelines(bed_lines)

    with open(genome_path, "w") as fh:
        for chrom in sorted(genome_sizes):
            fh.write(f"{chrom}\t{genome_sizes[chrom]}\n")

    with open(summary_path, "w") as fh:
        fh.write("gene_id\tchrom\tstart0\tend\n")
        for chrom in sorted(gene_entries):
            for start0, end, gene_id in sorted(gene_entries[chrom], key=lambda x: x[0]):
                fh.write(f"{gene_id}\t{chrom}\t{start0}\t{end}\n")

    print(f"Wrote BED: {bed_path}")
    print(f"Wrote genome: {genome_path}")


def main():
    args = parse_args()
    if args.dry_run:
        print(f"Would convert {args.gff} -> {args.out}")
        return
    gene_entries, genome_sizes = parse_gff(Path(args.gff), min_scaffold_genes=args.min_scaffold_genes)
    write_output(gene_entries, genome_sizes, Path(args.out))


if __name__ == "__main__":
    main()
