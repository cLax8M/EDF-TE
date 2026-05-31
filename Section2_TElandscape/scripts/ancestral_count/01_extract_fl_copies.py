#!/usr/bin/env python3
"""
01_extract_fl_copies.py — Per-species, parse EarlGrey filteredRepeats.gff and
                          combined_library.fasta to produce a per-family table
                          of full-length copies (copy_length >= FL_FRAC * consensus_length).

Output: one TSV per species:
  {out_dir}/{species}_fl_families.tsv
  Columns: family_name  consensus_length  total_copies  fl_copies  max_copy_len

A family is considered "full-length present" if fl_copies >= 1.

Usage:
  python 01_extract_fl_copies.py --annotated_dir /data/.../TE_analysis/annotated \
                                  --out_dir /data/.../ancestral_count/data/fl_families \
                                  --fl_frac 0.8

Author: EDF-TE project
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotated_dir", required=True,
                   help="Directory with {SPECIES}_EarlGrey/ subdirs")
    p.add_argument("--out_dir", required=True)
    p.add_argument("--fl_frac", type=float, default=0.8,
                   help="Minimum copy_len/consensus_len to count as full-length (default 0.8)")
    p.add_argument("--species", default=None,
                   help="Optional: restrict to this single species code")
    return p.parse_args()


def load_consensus_lengths(fasta_path):
    """Return dict consensus_name -> length (bp). Strips '#class' suffix."""
    lens = {}
    name = None
    L = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = L
                name = line[1:].rstrip().split()[0]
                # strip any '#class/subclass' suffix to match GFF Name=
                name = name.split("#")[0].upper()  # uppercase to match GFF
                L = 0
            else:
                L += len(line.rstrip())
        if name is not None:
            lens[name] = L
    return lens


_NAME_RE = re.compile(r'Name=([^;]+)')


def parse_gff_copies(gff_path):
    """Yield (family_name, copy_length) from EarlGrey filteredRepeats.gff."""
    with open(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            m = _NAME_RE.search(cols[8])
            if not m:
                continue
            fam = m.group(1).split("#")[0].upper()  # uppercase to match consensus
            yield fam, end - start + 1


def process_species(eg_dir, out_dir, fl_frac):
    species = eg_dir.name.replace("_EarlGrey", "")
    summary_dir = eg_dir / f"{species}_summaryFiles"
    gff = summary_dir / f"{species}.filteredRepeats.gff"
    lib = summary_dir / f"{species}_combined_library.fasta"

    if not gff.exists() or not lib.exists():
        print(f"  SKIP {species}: missing gff or library", file=sys.stderr)
        return None

    cons_len = load_consensus_lengths(lib)

    n_total = defaultdict(int)
    n_fl = defaultdict(int)
    max_len = defaultdict(int)

    for fam, copy_len in parse_gff_copies(gff):
        n_total[fam] += 1
        if copy_len > max_len[fam]:
            max_len[fam] = copy_len
        cl = cons_len.get(fam)
        if cl is None or cl == 0:
            continue
        if copy_len >= fl_frac * cl:
            n_fl[fam] += 1

    out_file = out_dir / f"{species}_fl_families.tsv"
    with open(out_file, "w") as fh:
        fh.write("family_name\tconsensus_length\ttotal_copies\tfl_copies\tmax_copy_len\n")
        for fam in sorted(n_total):
            cl = cons_len.get(fam, 0)
            fh.write(f"{fam}\t{cl}\t{n_total[fam]}\t{n_fl[fam]}\t{max_len[fam]}\n")

    n_fams = len(n_total)
    n_present = sum(1 for f in n_total if n_fl[f] >= 1)
    print(f"  {species}: {n_fams} families | {n_present} full-length present")
    return out_file


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    annot_dir = Path(args.annotated_dir)

    eg_dirs = sorted([d for d in annot_dir.iterdir()
                      if d.is_dir() and d.name.endswith("_EarlGrey")])
    if args.species:
        eg_dirs = [d for d in eg_dirs if d.name == f"{args.species}_EarlGrey"]

    print(f"Processing {len(eg_dirs)} species (fl_frac={args.fl_frac})")
    n_done = 0
    for eg in eg_dirs:
        if process_species(eg, out_dir, args.fl_frac) is not None:
            n_done += 1
    print(f"Done: {n_done}/{len(eg_dirs)} species")


if __name__ == "__main__":
    main()
