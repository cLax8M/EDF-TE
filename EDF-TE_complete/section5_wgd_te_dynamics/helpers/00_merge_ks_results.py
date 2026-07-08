#!/usr/bin/env python3
import argparse
import glob
import os
import sys


def main():
    parser = argparse.ArgumentParser(description="Merge per-species Ks tables into a global dataset.")
    parser.add_argument("--root", default="results/ks/paml_yn00", help="Root directory with <species>/<species>.ks.tsv")
    parser.add_argument("--out", default="results/ks/paml_yn00/all_species_ks.tsv", help="Output merged TSV")
    args = parser.parse_args()

    root = args.root.rstrip("/")
    pattern = os.path.join(root, "*", "*.ks.tsv")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[merge_ks] ERROR: no Ks tables found under {pattern}", file=sys.stderr)
        sys.exit(1)

    out_path = args.out
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w") as out_fh:
        out_fh.write("species\tgene1\tgene2\tKs\n")
        for fpath in files:
            try:
                with open(fpath, "r") as fh:
                    header = fh.readline()
                    for line in fh:
                        line = line.rstrip("\n")
                        if not line:
                            continue
                        out_fh.write(line + "\n")
            except OSError as e:
                print(f"[merge_ks] WARNING: cannot read {fpath}: {e}", file=sys.stderr)
                continue

    print(f"[merge_ks] Wrote merged Ks table to {out_path} (from {len(files)} species tables)", file=sys.stderr)


if __name__ == "__main__":
    main()
