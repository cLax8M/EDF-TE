#!/usr/bin/env python3
"""
build_fl_beds.py — For each EarlGrey species, emit a FL-filtered BED
(copy_len >= fl_frac * consensus_len, default 0.8) parallel to
{SP}.filteredRepeats.bed.

Usage:
    python3 build_fl_beds.py <annotated_dir> <out_dir> [--fl_frac 0.8]
"""
import os, sys, re, argparse
from collections import defaultdict

_NAME_RE = re.compile(r'Name=([^;]+)')


def load_consensus(lib_fa):
    out = {}
    name, L = None, 0
    with open(lib_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = L
                hdr = line[1:].rstrip().split()[0]
                fam = hdr.split("#", 1)[0]
                name = fam.upper()
                L = 0
            else:
                L += len(line.rstrip())
        if name is not None:
            out[name] = L
    return out


def process(eg_dir, out_dir, fl_frac):
    sp = os.path.basename(eg_dir).replace("_EarlGrey", "")
    summary = os.path.join(eg_dir, f"{sp}_summaryFiles")
    gff = os.path.join(summary, f"{sp}.filteredRepeats.gff")
    bed = os.path.join(summary, f"{sp}.filteredRepeats.bed")
    lib = os.path.join(summary, f"{sp}_combined_library.fasta")
    if not (os.path.isfile(gff) and os.path.isfile(bed) and os.path.isfile(lib)):
        return None

    cons = load_consensus(lib)

    # Collect (chrom, start, end) of FL copies from GFF
    fl_set = set()
    with open(gff) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            try:
                start = int(cols[3]); end = int(cols[4])
            except ValueError:
                continue
            m = _NAME_RE.search(cols[8])
            if not m:
                continue
            fam = m.group(1).split("#")[0].upper()
            cl = cons.get(fam, 0)
            if cl == 0:
                continue
            if (end - start + 1) >= fl_frac * cl:
                # Store as (chrom, 0-based start, end) key
                fl_set.add((cols[0], start - 1, end))

    # Walk BED and keep only rows matching FL set
    out_bed = os.path.join(out_dir, f"{sp}.filteredRepeats.fl.bed")
    n_kept = n_total = 0
    with open(bed) as fh, open(out_bed, "w") as out:
        for ln in fh:
            if not ln or ln.startswith("track") or ln.startswith("#"):
                out.write(ln); continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            n_total += 1
            try:
                chrom = p[0]; bs = int(p[1]); be = int(p[2])
            except ValueError:
                continue
            if (chrom, bs, be) in fl_set:
                out.write(ln); n_kept += 1
    return sp, n_kept, n_total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("annotated_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--fl_frac", type=float, default=0.8)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    for entry in sorted(os.listdir(args.annotated_dir)):
        if not entry.endswith("_EarlGrey"):
            continue
        eg = os.path.join(args.annotated_dir, entry)
        if not os.path.isdir(eg):
            continue
        r = process(eg, args.out_dir, args.fl_frac)
        if r is None:
            print(f"[SKIP] {entry}", file=sys.stderr); continue
        sp, kept, total = r
        frac = (kept / total * 100) if total else 0
        print(f"{sp}\t{kept}\t{total}\t{frac:.1f}%")


if __name__ == "__main__":
    main()
