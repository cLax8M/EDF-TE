#!/usr/bin/env python3
"""Extract 1:1 ortholog pairs from OrthoFinder pairwise tables."""

import argparse
import sys
from pathlib import Path
from itertools import combinations


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ortho_dir", required=True, help="Path to OrthoFinder Orthologues/ directory")
    p.add_argument("--out_dir", required=True, help="Directory for .pairs.tsv files")
    p.add_argument("--dry_run", action="store_true")
    return p.parse_args()


def find_orthologues_files(orth_dir):
    orth_dir = Path(orth_dir)
    files = []
    for d in orth_dir.glob("Orthologues_*.aa"):
        files.extend(d.glob("*.tsv"))
    for d in orth_dir.glob("Orthologues_*"):
        if d.is_dir() and d.name != "Orthologues_*.aa":
            files.extend(d.glob("*.tsv"))
    return sorted(set(files))


def parse_species_from_path(path):
    name = path.name
    if "__v__" in name:
        left, right = name.split("__v__", 1)
        sp_a = left.replace(".aa", "")
        sp_b = right.replace(".aa", "").replace(".tsv", "")
        return sp_a, sp_b
    if ".tsv" in name:
        stem = name[:-4]
        if "__v__" in stem:
            sp_a, sp_b = stem.split("__v__", 1)
            return sp_a, sp_b
    return None, None


def extract_one_to_one(tsv_path):
    pairs = []
    with open(tsv_path) as fh:
        header = fh.readline()
        if not header:
            return pairs
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            og = cols[0]
            a = cols[1].strip()
            b = cols[2].strip()
            if not a or not b:
                continue
            if "," in a or "," in b:
                continue
            pairs.append((og, a, b))
    return pairs


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = find_orthologues_files(args.ortho_dir)
    if args.dry_run:
        print(f"Would process {len(files)} Orthofinder files")
        return

    processed = 0
    for path in files:
        sp_a, sp_b = parse_species_from_path(path)
        if not sp_a or not sp_b:
            continue
        pairs = extract_one_to_one(path)
        out_path = out_dir / f"{sp_a}__{sp_b}.pairs.tsv"
        with open(out_path, "w") as fh:
            fh.write("orthogroup\tgene_A\tgene_B\n")
            for og, ga, gb in pairs:
                fh.write(f"{og}\t{ga}\t{gb}\n")
        processed += 1
    print(f"Wrote {processed} pair files to {out_dir}")


if __name__ == "__main__":
    main()
