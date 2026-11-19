#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
from typing import Optional, Tuple


def parse_ds_from_yn00(path: Path) -> Optional[float]:
    """Parse dS (Ks) from a PAML yn00 output file.

    We take the Yang & Nielsen (2000) section (B), and from the first
    data line under that header we read the third value from the end,
    which corresponds to dS (before '+-' and its SE).
    """

    try:
        with path.open() as fh:
            in_section_b = False
            for line in fh:
                if "(B) Yang & Nielsen" in line:
                    in_section_b = True
                    continue
                if not in_section_b:
                    continue
                stripped = line.strip()
                if not stripped:
                    continue
                # Expect a line starting with two integers: seq. seq.
                parts = stripped.split()
                if len(parts) < 5:
                    continue
                if not (parts[0].isdigit() and parts[1].isdigit()):
                    continue
                # Take the third-to-last token as dS
                try:
                    ds = float(parts[-3])
                except ValueError:
                    return None
                return ds
    except OSError:
        return None
    return None


def process_species(sp_dir: Path, species: str) -> int:
    pairs_dir = sp_dir / "pairs"
    mapping = pairs_dir / "pairs_index.tsv"
    if not mapping.is_file():
        return 0

    out_path = sp_dir / f"{species}.ks.tsv"
    n_written = 0
    with mapping.open() as mf, out_path.open("w") as out:
        out.write("species\tgene1\tgene2\tKs\n")
        header = next(mf, None)
        for line in mf:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            pair_id, g1, g2 = parts[0], parts[1], parts[2]
            yn_path = pairs_dir / f"{pair_id}.yn00.out"
            if not yn_path.is_file():
                continue
            ks = parse_ds_from_yn00(yn_path)
            if ks is None or ks <= 0:
                continue
            out.write(f"{species}\t{g1}\t{g2}\t{ks}\n")
            n_written += 1
    return n_written


def main() -> None:
    ap = argparse.ArgumentParser(description="Rebuild per-species Ks tables by parsing existing yn00 outputs.")
    ap.add_argument("--root", default="results/ks/paml_yn00", help="Root directory with <species>/pairs/*.yn00.out")
    ap.add_argument("--species-list", help="Optional file with species names (one per line). If omitted, infer from subdirectories of root.")
    args = ap.parse_args()

    root = Path(args.root)
    if not root.is_dir():
        raise SystemExit(f"Root directory not found: {root}")

    if args.species_list:
        species = [l.strip() for l in Path(args.species_list).read_text().splitlines() if l.strip()]
    else:
        species = [p.name for p in root.iterdir() if p.is_dir()]

    total_pairs = 0
    sp_with_any = 0
    for sp in species:
        sp_dir = root / sp
        if not sp_dir.is_dir():
            continue
        n = process_species(sp_dir, sp)
        if n > 0:
            sp_with_any += 1
            total_pairs += n
            print(f"[rebuild_ks] {sp}: wrote {n} pairs")
        else:
            print(f"[rebuild_ks] {sp}: no Ks pairs written (no yn00 outputs or all invalid)")

    print(f"[rebuild_ks] Species with Ks: {sp_with_any}; total pairs: {total_pairs}")


if __name__ == "__main__":
    main()
