#!/usr/bin/env python3
"""
02_collect_consensus.py — Collect FULL-LENGTH-present consensus sequences from
                          all species into a single FASTA, with headers tagged
                          with the species code, ready for cross-species clustering.

Header format: >{SPECIES}|{family_name}

Inputs:
  - {fl_dir} with per-species *_fl_families.tsv (from step 01)
  - {annotated_dir} with EarlGrey outputs (combined_library.fasta sources)

Output:
  - {out_fasta}  multi-FASTA with one entry per (species, family-with-FL-copy)

Usage:
  python 02_collect_consensus.py \
      --fl_dir /path/to/fl_families \
      --annotated_dir /data/.../annotated \
      --out_fasta /path/to/all_fl_consensus.fa

Author: EDF-TE project
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fl_dir", required=True)
    p.add_argument("--annotated_dir", required=True)
    p.add_argument("--out_fasta", required=True)
    return p.parse_args()


def load_fl_families(tsv):
    """Return set of family_names with fl_copies >= 1."""
    fams = set()
    with open(tsv) as fh:
        next(fh)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            try:
                if int(cols[3]) >= 1:
                    fams.add(cols[0])  # already uppercase from step 01
            except ValueError:
                pass
    return fams


def stream_library(lib_path, wanted, species, out_fh):
    """Write entries from lib_path with name in `wanted` to out_fh, retagged."""
    keep = False
    n = 0
    with open(lib_path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].rstrip().split()[0].split("#")[0].upper()
                if name in wanted:
                    out_fh.write(f">{species}|{name}\n")
                    keep = True
                    n += 1
                else:
                    keep = False
            else:
                if keep:
                    out_fh.write(line)
    return n


def main():
    args = parse_args()
    fl_dir = Path(args.fl_dir)
    annot_dir = Path(args.annotated_dir)
    out_fa = Path(args.out_fasta)
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    tsvs = sorted(fl_dir.glob("*_fl_families.tsv"))
    print(f"Collecting from {len(tsvs)} species → {out_fa}")

    total = 0
    with open(out_fa, "w") as out_fh:
        for tsv in tsvs:
            species = tsv.stem.replace("_fl_families", "")
            wanted = load_fl_families(tsv)
            if not wanted:
                continue
            lib = (annot_dir / f"{species}_EarlGrey" /
                   f"{species}_summaryFiles" /
                   f"{species}_combined_library.fasta")
            if not lib.exists():
                print(f"  SKIP {species}: library missing", file=sys.stderr)
                continue
            n = stream_library(lib, wanted, species, out_fh)
            total += n
            print(f"  {species}: {n} FL-present consensus written")

    print(f"Total: {total} sequences written to {out_fa}")


if __name__ == "__main__":
    main()
