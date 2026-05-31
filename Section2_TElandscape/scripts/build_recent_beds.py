#!/usr/bin/env python3
"""
build_recent_beds.py — For each EarlGrey species, emit two BEDs in parallel
to {SP}.filteredRepeats.bed:

    {SP}.filteredRepeats.recent.bed     -- copies with KIMURA80 < THRESH
    {SP}.filteredRepeats.fl_recent.bed  -- copies that are FL AND recent
                                          (FL = copy_len >= fl_frac * consensus_len)

Source of divergence: KIMURA80 attribute in
    {EG}/{SP}_RepeatLandscape/{SP}.filteredRepeats.withDivergence.gff
Source of consensus length: {SP}_combined_library.fasta in summaryFiles/.

Usage:
    python3 build_recent_beds.py <annotated_dir> <out_dir> \\
        [--fl_frac 0.8] [--kimura_thresh 0.10]
"""
import os, sys, re, argparse

_NAME_RE = re.compile(r'NAME=([^;]+)')
_KIM_RE  = re.compile(r'KIMURA80=([0-9.eE+-]+)')


def load_consensus(lib_fa):
    out, name, L = {}, None, 0
    with open(lib_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = L
                hdr = line[1:].rstrip().split()[0]
                name = hdr.split("#", 1)[0].upper()
                L = 0
            else:
                L += len(line.rstrip())
        if name is not None:
            out[name] = L
    return out


def process(eg_dir, out_dir, fl_frac, kim_thresh):
    sp = os.path.basename(eg_dir).replace("_EarlGrey", "")
    summary = os.path.join(eg_dir, f"{sp}_summaryFiles")
    landscape = os.path.join(eg_dir, f"{sp}_RepeatLandscape")
    bed = os.path.join(summary, f"{sp}.filteredRepeats.bed")
    lib = os.path.join(summary, f"{sp}_combined_library.fasta")
    div_gff = os.path.join(landscape, f"{sp}.filteredRepeats.withDivergence.gff")

    if not (os.path.isfile(bed) and os.path.isfile(lib) and os.path.isfile(div_gff)):
        return None

    cons = load_consensus(lib)

    # Map (chrom, 0-based start, end) -> (kim, is_FL)
    info = {}
    with open(div_gff) as fh:
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
            mk = _KIM_RE.search(cols[8])
            if not mk:
                continue
            kim = float(mk.group(1))
            mn = _NAME_RE.search(cols[8])
            fam = mn.group(1).split("#")[0].upper() if mn else ""
            cl = cons.get(fam, 0)
            is_fl = (cl > 0) and ((end - start + 1) >= fl_frac * cl)
            info[(cols[0], start - 1, end)] = (kim, is_fl)

    out_recent    = os.path.join(out_dir, f"{sp}.filteredRepeats.recent.bed")
    out_fl_recent = os.path.join(out_dir, f"{sp}.filteredRepeats.fl_recent.bed")
    n_total = n_recent = n_fl_recent = 0
    with open(bed) as fh, open(out_recent, "w") as orec, open(out_fl_recent, "w") as ofl:
        for ln in fh:
            if not ln or ln.startswith("track") or ln.startswith("#"):
                orec.write(ln); ofl.write(ln); continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            n_total += 1
            try:
                chrom = p[0]; bs = int(p[1]); be = int(p[2])
            except ValueError:
                continue
            v = info.get((chrom, bs, be))
            if v is None:
                continue
            kim, is_fl = v
            if kim < kim_thresh:
                orec.write(ln); n_recent += 1
                if is_fl:
                    ofl.write(ln); n_fl_recent += 1
    return sp, n_recent, n_fl_recent, n_total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("annotated_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--fl_frac",       type=float, default=0.8)
    ap.add_argument("--kimura_thresh", type=float, default=0.10)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print("species\trecent\tfl_recent\ttotal\trecent%\tfl_recent%")
    for entry in sorted(os.listdir(args.annotated_dir)):
        if not entry.endswith("_EarlGrey"):
            continue
        eg = os.path.join(args.annotated_dir, entry)
        if not os.path.isdir(eg):
            continue
        r = process(eg, args.out_dir, args.fl_frac, args.kimura_thresh)
        if r is None:
            print(f"[SKIP] {entry}", file=sys.stderr); continue
        sp, n_rec, n_flr, n_tot = r
        f1 = (n_rec / n_tot * 100) if n_tot else 0
        f2 = (n_flr / n_tot * 100) if n_tot else 0
        print(f"{sp}\t{n_rec}\t{n_flr}\t{n_tot}\t{f1:.1f}%\t{f2:.1f}%")


if __name__ == "__main__":
    main()
