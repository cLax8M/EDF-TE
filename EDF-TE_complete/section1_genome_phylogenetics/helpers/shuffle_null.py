#!/usr/bin/env python3
"""Compute one shuffle replicate for the synteny null model."""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--collinear_dir", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--bed_dir", required=True)
    p.add_argument("--ortho_dir", required=True)
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args()


def read_summary(summary_path):
    rows = []
    with open(summary_path) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 5:
                rows.append((cols[0], cols[1], int(cols[2]), int(cols[3]), int(cols[4])))
    return rows


def compute_c_shuffled(sp_a, sp_b, bed_dir, ortho_dir, out_dir, seed):
    tmp_dir = Path(tempfile.mkdtemp(prefix="shuffle_", dir="/tmp"))
    try:
        bed_a = Path(bed_dir) / f"{sp_a}.bed"
        bed_b = Path(bed_dir) / f"{sp_b}.bed"
        geno_a = Path(bed_dir) / f"{sp_a}.genome"
        geno_b = Path(bed_dir) / f"{sp_b}.genome"
        if not bed_a.exists() or not bed_b.exists() or not geno_a.exists() or not geno_b.exists():
            return 0

        bedtools = shutil.which("bedtools")
        if not bedtools:
            return 0

        tmp_beds = tmp_dir / "beds"
        tmp_beds.mkdir(exist_ok=True)
        tmp_pairs = tmp_dir / "pairs"
        tmp_pairs.mkdir(exist_ok=True)
        tmp_out = tmp_dir / "collinear"
        tmp_out.mkdir(exist_ok=True)

        shuffled_a = tmp_beds / f"{sp_a}.bed"
        shuffled_b = tmp_beds / f"{sp_b}.bed"
        subprocess.run([bedtools, "shuffle", "-i", str(bed_a), "-g", str(geno_a), "-chrom", "-noOverlapping", "-seed", str(seed)], check=True, stdout=open(shuffled_a, "w"))
        subprocess.run([bedtools, "shuffle", "-i", str(bed_b), "-g", str(geno_b), "-chrom", "-noOverlapping", "-seed", str(seed + 1)], check=True, stdout=open(shuffled_b, "w"))

        pair_file = None
        for candidate in Path(ortho_dir).glob("*.pairs.tsv"):
            if candidate.stem.startswith(sp_a + "__") and candidate.stem.endswith(sp_b):
                pair_file = candidate
                break
        if pair_file is None:
            for candidate in Path(ortho_dir).rglob("*.pairs.tsv"):
                if candidate.stem.startswith(sp_a + "__") and candidate.stem.endswith(sp_b):
                    pair_file = candidate
                    break
        if pair_file is None:
            return 0
        shutil.copy2(pair_file, tmp_pairs / pair_file.name)

        script = Path(__file__).resolve()
        # call the sibling collinear_pairs script
        subprocess.run([
            sys.executable, str(script.with_name("collinear_pairs.py")),
            "--bed_dir", str(tmp_beds),
            "--ortho_dir", str(tmp_pairs),
            "--out_dir", str(tmp_out),
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        summary_path = tmp_out / sp_a.split("_")[0] / f"{sp_a.split('_')[0]}_collinear_summary.tsv"
        if not summary_path.exists():
            return 0
        rows = read_summary(summary_path)
        for r_sp_a, r_sp_b, n_orthologs, c_ij, n_adj in rows:
            if r_sp_a == sp_a and r_sp_b == sp_b:
                return c_ij
        return 0
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def append_result(out_dir, clade, sp_a, sp_b, c_ij, seed):
    clade_out = Path(out_dir) / clade
    clade_out.mkdir(parents=True, exist_ok=True)
    summary_path = clade_out / f"{sp_a}__{sp_b}.shuffle_summary.tsv"
    with open(summary_path, "a") as fh:
        if fh.tell() == 0:
            fh.write("replicate\tc_shuffled\n")
        fh.write(f"{seed}\t{c_ij}\n")

    # Recompute means and SD for this pair
    rows = []
    with open(summary_path) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                rows.append(float(cols[1]))
    if not rows:
        return
    mean = sum(rows) / len(rows)
    var = sum((x - mean) ** 2 for x in rows) / len(rows) if len(rows) > 1 else 0.0
    sd = var ** 0.5
    null_means = clade_out / f"{clade}_null_means.tsv"
    # rewrite the whole file for simplicity
    pairs = []
    if null_means.exists():
        with open(null_means) as fh:
            first = fh.readline()
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 4:
                    pairs.append((cols[0], cols[1], cols[2], cols[3]))
    new_pairs = []
    found = False
    for a, b, m, s in pairs:
        if a == sp_a and b == sp_b:
            new_pairs.append((a, b, f"{mean:.4f}", f"{sd:.4f}"))
            found = True
        else:
            new_pairs.append((a, b, m, s))
    if not found:
        new_pairs.append((sp_a, sp_b, f"{mean:.4f}", f"{sd:.4f}"))
    with open(null_means, "w") as fh:
        fh.write("species_A\tspecies_B\ts_ij_mean\ts_ij_sd\n")
        for a, b, m, s in new_pairs:
            fh.write(f"{a}\t{b}\t{m}\t{s}\n")


def main():
    args = parse_args()
    collinear_dir = Path(args.collinear_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for summary_path in sorted(collinear_dir.glob("*/ *_collinear_summary.tsv")):
        pass
    for clade_dir in sorted(collinear_dir.iterdir()):
        if not clade_dir.is_dir():
            continue
        summary_path = clade_dir / f"{clade_dir.name}_collinear_summary.tsv"
        if not summary_path.exists():
            continue
        for sp_a, sp_b, n_orthologs, c_ij, n_adj in read_summary(summary_path):
            c_shuffled = compute_c_shuffled(sp_a, sp_b, args.bed_dir, args.ortho_dir, out_dir, args.seed)
            append_result(out_dir, clade_dir.name, sp_a, sp_b, c_shuffled, args.seed)
    print(f"Wrote shuffle results under {out_dir}")


if __name__ == "__main__":
    main()
