#!/usr/bin/env python3
"""
00_genome_qc_filter.py — QC filter for fragmented genomes used in synteny analysis.

Excludes species from downstream synteny computations if their genome assembly
is too fragmented to yield meaningful collinearity. Criteria (PASS if EITHER):
  - N50 >= 50,000 bp (computed from genome FASTA), OR
  - BUSCO completeness >= 80% (read from optional BUSCO summary table).

Inputs:
  --genome_dir DIR     Directory with per-species genome FASTA files
                       (filenames like {SPECIES_PREFIX}.{fasta,fna})
  --busco_table FILE   (optional) TSV with columns: species_prefix  busco_complete_pct
  --out_dir DIR        Output directory

Outputs:
  {out_dir}/species_whitelist.txt   one species_prefix per line (PASS)
  {out_dir}/species_qc_summary.tsv  full QC table (species, n50, busco_pct, status, reason)

Usage:
  python 00_genome_qc_filter.py --genome_dir /legserv/Temp/Carlos/EDF-TE/data/genomes \
      --busco_table /legserv/Temp/Carlos/EDF-TE/results/QC/busco_summary.tsv \
      --out_dir /data/carlos/EDF-TE/EDF-TE/Synteny/data/qc

Author: EDF-TE project
"""

import argparse
import sys
from pathlib import Path


N50_THRESHOLD = 50_000
BUSCO_THRESHOLD = 80.0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genome_dir", required=True)
    p.add_argument("--busco_table", default=None)
    p.add_argument("--out_dir", required=True)
    return p.parse_args()


def compute_n50(fasta_path):
    """Compute N50 from a FASTA file (no biopython dependency)."""
    lengths = []
    cur = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur > 0:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.rstrip())
        if cur > 0:
            lengths.append(cur)

    if not lengths:
        return 0

    lengths.sort(reverse=True)
    total = sum(lengths)
    cum = 0
    for L in lengths:
        cum += L
        if cum >= total / 2:
            return L
    return 0


def load_busco_table(path):
    """Load optional BUSCO table: species_prefix -> busco_complete_pct (float)."""
    data = {}
    if not path or not Path(path).exists():
        return data
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_sp = header.index("species_prefix")
            i_b = header.index("busco_complete_pct")
        except ValueError:
            print(f"  WARNING: BUSCO table missing required columns: {path}",
                  file=sys.stderr)
            return data
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(i_sp, i_b):
                continue
            try:
                data[cols[i_sp]] = float(cols[i_b])
            except ValueError:
                pass
    return data


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    genome_dir = Path(args.genome_dir)
    if not genome_dir.is_dir():
        sys.exit(f"ERROR: genome_dir not found: {genome_dir}")

    busco = load_busco_table(args.busco_table)

    whitelist = []
    rows = []

    fastas = sorted(list(f for f in genome_dir.glob("*.fasta") if not f.name.startswith("._")) + list(f for f in genome_dir.glob("*.fna") if not f.name.startswith("._")))
    print(f"Scanning {len(fastas)} genome files...")

    for fa in fastas:
        if fa.is_symlink():
            continue
        species = fa.stem
        print(f"Processing {fa}")
        n50 = compute_n50(fa)
        busco_pct = busco.get(species)

        pass_n50 = n50 >= N50_THRESHOLD
        pass_busco = (busco_pct is not None) and (busco_pct >= BUSCO_THRESHOLD)
        status = "PASS" if (pass_n50 or pass_busco) else "EXCLUDE"
        reasons = []
        if pass_n50:
            reasons.append(f"N50={n50}>={N50_THRESHOLD}")
        if pass_busco:
            reasons.append(f"BUSCO={busco_pct:.1f}>={BUSCO_THRESHOLD}")
        if status == "EXCLUDE":
            reasons.append(f"N50={n50}<{N50_THRESHOLD}")
            if busco_pct is not None:
                reasons.append(f"BUSCO={busco_pct:.1f}<{BUSCO_THRESHOLD}")
            else:
                reasons.append("BUSCO=NA")

        if status == "PASS":
            whitelist.append(species)

        rows.append((species, n50, busco_pct, status, ";".join(reasons)))

    # Write whitelist
    whitelist_file = out_dir / "species_whitelist.txt"
    with open(whitelist_file, "w") as fh:
        for sp in whitelist:
            fh.write(sp + "\n")

    # Write summary
    summary_file = out_dir / "species_qc_summary.tsv"
    with open(summary_file, "w") as fh:
        fh.write("species_prefix\tn50\tbusco_complete_pct\tstatus\treason\n")
        for sp, n50, busco_pct, status, reason in rows:
            b = "" if busco_pct is None else f"{busco_pct:.2f}"
            fh.write(f"{sp}\t{n50}\t{b}\t{status}\t{reason}\n")

    n_pass = len(whitelist)
    n_excl = len(rows) - n_pass
    print(f"  PASS:    {n_pass}")
    print(f"  EXCLUDE: {n_excl}")
    print(f"  Whitelist: {whitelist_file}")
    print(f"  Summary:   {summary_file}")


if __name__ == "__main__":
    main()
