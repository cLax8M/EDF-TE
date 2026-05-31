#!/usr/bin/env python3
"""
F9_scaled_body_meta.py
Per-species meta-TE profile with scale-regions semantics like deepTools
computeMatrix scale-regions:

  [-FLANK ......... TSS] [BODY scaled to BODY_BINS] [TES ......... +FLANK]

Each TE >= MIN_LEN is scaled to BODY_BINS bins; flanks are sampled at
1 bp/bin and binned into FLANK_BINS bins each.

Output columns (TSV):
  species, te_class, region, bin_idx, mean_occupancy, n_insertions

`region` ∈ {upstream, body, downstream}; bin_idx is local to the region.
"""
import argparse
import math
import os
import statistics
import sys
from collections import defaultdict

# -- WIG parser (same as F6) -------------------------------------------------
def parse_wig(path):
    scores = defaultdict(dict)
    chrom = None
    step = None
    pos = None
    span = 1
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            if line.startswith("variableStep"):
                step = "var"
                for tok in line.split()[1:]:
                    k, v = tok.split("=")
                    if k == "chrom": chrom = v
                    elif k == "span": span = int(v)
                continue
            if line.startswith("fixedStep"):
                step = "fix"
                start = 1; sstep = 1; span = 1
                for tok in line.split()[1:]:
                    k, v = tok.split("=")
                    if k == "chrom": chrom = v
                    elif k == "start": start = int(v)
                    elif k == "step": sstep = int(v)
                    elif k == "span": span = int(v)
                pos = start
                continue
            try:
                if step == "var":
                    p, v = line.split()
                    p = int(p); v = float(v)
                    scores[chrom][p] = v
                elif step == "fix":
                    v = float(line)
                    scores[chrom][pos] = v
                    pos += sstep
            except Exception:
                continue
    return scores


def parse_class(name: str) -> str:
    if "#" in name:
        name = name.split("#", 1)[1]
    head = name.split("/", 1)[0].rstrip("?")
    return {
        "DNA": "DNA", "LTR": "LTR", "LINE": "LINE", "SINE": "SINE",
        "RC": "Rolling Circle", "Helitron": "Rolling Circle",
        "PLE": "PLE", "Penelope": "PLE",
        "Retroposon": "LINE",
    }.get(head, "Other")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("wig")
    p.add_argument("bed")
    p.add_argument("species_id")
    p.add_argument("output_tsv")
    p.add_argument("--flank", type=int, default=500)
    p.add_argument("--flank_bins", type=int, default=50)
    p.add_argument("--body_bins", type=int, default=100)
    p.add_argument("--min_len",  type=int, default=200)
    p.add_argument("--max_len",  type=int, default=20000)
    a = p.parse_args()

    print(f"[{a.species_id}] loading WIG {a.wig}")
    raw = parse_wig(a.wig)
    prefix = a.species_id + "."
    scores = {(c[len(prefix):] if c.startswith(prefix) else c): v
              for c, v in raw.items()}

    # acc[(class, region, bin)] -> list of values
    acc = defaultdict(list)
    n_te = defaultdict(int)
    flank_step = a.flank / a.flank_bins  # bp per bin

    with open(a.bed) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"): continue
            f = ln.split("\t")
            if len(f) < 4: continue
            chrom = f[0]
            try:
                s = int(f[1]); e = int(f[2])
            except ValueError: continue
            tlen = e - s
            if tlen < a.min_len or tlen > a.max_len: continue
            cls = parse_class(f[3])
            if cls == "Other": continue
            cd = scores.get(chrom)
            if not cd: continue
            strand = f[5] if len(f) >= 6 else "+"
            n_te[cls] += 1

            # ── Upstream flank
            for b in range(a.flank_bins):
                offset = b * flank_step + flank_step / 2
                if strand == "-":
                    pos = e + int(round(offset)) + 1   # 1-based wig
                else:
                    pos = s - a.flank + int(round(offset)) + 1
                v = cd.get(pos)
                if v is not None:
                    acc[(cls, "upstream", b if strand != "-" else (a.flank_bins - 1 - b))].append(v)

            # ── Body scaled
            body_step = tlen / a.body_bins
            for b in range(a.body_bins):
                offset = b * body_step + body_step / 2
                if strand == "-":
                    pos = e - int(round(offset)) + 1
                else:
                    pos = s + int(round(offset)) + 1
                v = cd.get(pos)
                if v is not None:
                    acc[(cls, "body", b)].append(v)

            # ── Downstream flank
            for b in range(a.flank_bins):
                offset = b * flank_step + flank_step / 2
                if strand == "-":
                    pos = s - int(round(offset))   # subtract since reversed
                else:
                    pos = e + int(round(offset)) + 1
                v = cd.get(pos)
                if v is not None:
                    acc[(cls, "downstream", b if strand != "-" else (a.flank_bins - 1 - b))].append(v)

    os.makedirs(os.path.dirname(a.output_tsv), exist_ok=True)
    with open(a.output_tsv, "w") as oh:
        oh.write("species\tte_class\tregion\tbin_idx\tmean_occupancy\tn_insertions\tn_obs\n")
        for (cls, reg, b), vals in sorted(acc.items()):
            if not vals: continue
            m = statistics.mean(vals)
            oh.write(f"{a.species_id}\t{cls}\t{reg}\t{b}\t{m:.4f}\t{n_te[cls]}\t{len(vals)}\n")
    print(f"[{a.species_id}] DONE  classes: {dict(n_te)}  -> {a.output_tsv}")


if __name__ == "__main__":
    main()
