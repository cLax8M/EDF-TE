#!/usr/bin/env python3
"""Identify collinear ortholog pairs from BED files and 1:1 ortholog tables."""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bed_dir", required=True, help="Directory with *.bed files")
    p.add_argument("--ortho_dir", required=True, help="Directory with *.pairs.tsv files")
    p.add_argument("--out_dir", required=True, help="Output directory for collinear results")
    return p.parse_args()


def load_bed(bed_path):
    chrom_genes = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            chrom, start, end, gene_id = cols[0], int(cols[1]), int(cols[2]), cols[3]
            chrom_genes[chrom].append((start, gene_id))
    gene_info = {}
    for chrom in chrom_genes:
        entries = sorted(chrom_genes[chrom], key=lambda x: x[0])
        for rank, (_, gene_id) in enumerate(entries):
            gene_info[gene_id] = (chrom, rank)
    return gene_info


def load_pairs(pairs_path):
    pairs = []
    with open(pairs_path) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                pairs.append((cols[1], cols[2]))
    return pairs


def find_collinear(gene_info_a, gene_info_b, pairs):
    a_to_b = {}
    b_to_a = {}
    for ga, gb in pairs:
        if ga in gene_info_a and gb in gene_info_b:
            a_to_b[ga] = gb
            b_to_a[gb] = ga
    chrom_orthologs_a = defaultdict(list)
    for ga in a_to_b:
        chrom, rank = gene_info_a[ga]
        chrom_orthologs_a[chrom].append((rank, ga))
    for chrom in chrom_orthologs_a:
        chrom_orthologs_a[chrom].sort(key=lambda x: x[0])

    chrom_orthologs_b = defaultdict(list)
    for gb in b_to_a:
        chrom, rank = gene_info_b[gb]
        chrom_orthologs_b[chrom].append((rank, gb))
    for chrom in chrom_orthologs_b:
        chrom_orthologs_b[chrom].sort(key=lambda x: x[0])

    gene_b_ortholog_rank = {}
    for chrom, entries in chrom_orthologs_b.items():
        for ortho_rank, gb in enumerate([g for _, g in entries]):
            gene_b_ortholog_rank[gb] = (chrom, ortho_rank)

    collinear = []
    n_adj = 0
    for chrom_a, entries in chrom_orthologs_a.items():
        if len(entries) < 2:
            continue
        n_adj += len(entries) - 1
        for i in range(len(entries) - 1):
            _, ga1 = entries[i]
            _, ga2 = entries[i + 1]
            gb1 = a_to_b[ga1]
            gb2 = a_to_b[ga2]
            cb1, rb1 = gene_b_ortholog_rank[gb1]
            cb2, rb2 = gene_b_ortholog_rank[gb2]
            if cb1 == cb2 and abs(rb1 - rb2) == 1:
                collinear.append((chrom_a, ga1, ga2, cb1, gb1, gb2))
    return collinear, n_adj


def main():
    args = parse_args()
    bed_dir = Path(args.bed_dir)
    ortho_dir = Path(args.ortho_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pair_files = sorted(ortho_dir.glob("*.pairs.tsv")) + sorted(ortho_dir.rglob("*.pairs.tsv"))
    pair_files = sorted(set(pair_files))
    if not pair_files:
        print("No pair files found", file=sys.stderr)
        return

    for pairs_path in pair_files:
        if pairs_path.name == "*.pairs.tsv":
            continue
        name = pairs_path.stem.replace(".pairs", "")
        if "__" not in name:
            continue
        sp_a, sp_b = name.split("__", 1)
        bed_a = bed_dir / f"{sp_a}.bed"
        bed_b = bed_dir / f"{sp_b}.bed"
        if not bed_a.exists() or not bed_b.exists():
            continue
        gene_info_a = load_bed(bed_a)
        gene_info_b = load_bed(bed_b)
        pairs = load_pairs(pairs_path)
        collinear, n_adj = find_collinear(gene_info_a, gene_info_b, pairs)
        clade = sp_a.split("_")[0]
        clade_out = out_dir / clade
        clade_out.mkdir(parents=True, exist_ok=True)
        out_path = clade_out / f"{sp_a}__{sp_b}.collinear.tsv"
        with open(out_path, "w") as fh:
            fh.write("scaffold_A\tgene_A1\tgene_A2\tscaffold_B\tgene_B1\tgene_B2\n")
            for row in collinear:
                fh.write("\t".join(map(str, row)) + "\n")
        summary_path = clade_out / f"{clade}_collinear_summary.tsv"
        with open(summary_path, "a") as fh:
            if fh.tell() == 0:
                fh.write("species_A\tspecies_B\tn_orthologs\tc_ij\tN_ij\n")
            fh.write(f"{sp_a}\t{sp_b}\t{len(pairs)}\t{len(collinear)}\t{n_adj}\n")
    print(f"Wrote collinear outputs under {out_dir}")


if __name__ == "__main__":
    main()
