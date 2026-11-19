#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path
from typing import Dict, Tuple


def parse_diag(path: Path) -> Tuple[int, int, int, int, int]:
    if not path.is_file():
        return 0, 0, 0, 0, 0
    total = ok = miss_prot = miss_cds = miss_both = 0
    with path.open() as fh:
        header = next(fh, None)
        for line in fh:
            line = line.strip()
            if not line:
                continue
            total += 1
            parts = line.split("\t")
            status = parts[-1] if parts else ""
            if status == "OK":
                ok += 1
            elif status == "MISSING_PROTEIN":
                miss_prot += 1
            elif status == "MISSING_CDS":
                miss_cds += 1
            elif status == "MISSING_BOTH":
                miss_both += 1
    return total, ok, miss_prot, miss_cds, miss_both


def count_ks_rows(path: Path) -> int:
    if not path.is_file():
        return 0
    with path.open() as fh:
        # subtract header if any
        n = sum(1 for _ in fh)
    return max(0, n - 1)


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize which species have Ks results and diagnostics.")
    ap.add_argument("--species-list", default="results/ks/species_for_paml_ks.txt")
    ap.add_argument("--root", default="results/ks/paml_yn00")
    ap.add_argument("--out", default="results/ks/paml_yn00/ks_species_status.tsv")
    args = ap.parse_args()

    species_file = Path(args.species_list)
    root = Path(args.root)
    out_path = Path(args.out)

    if not species_file.is_file():
        print(f"[status] ERROR: species list not found at {species_file}", file=sys.stderr)
        sys.exit(1)

    species = [l.strip() for l in species_file.read_text().splitlines() if l.strip()]
    if not species:
        print(f"[status] ERROR: species list {species_file} is empty", file=sys.stderr)
        sys.exit(1)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out:
        out.write(
            "species\thas_ks\tn_ks_pairs\tn_diag_genes\tn_diag_ok\t"
            "n_missing_protein\tn_missing_cds\tn_missing_both\treason\n"
        )
        for sp in species:
            sp_dir = root / sp
            ks_path = sp_dir / f"{sp}.ks.tsv"
            diag_path = sp_dir / f"{sp}.id_diagnostics.tsv"

            n_ks = count_ks_rows(ks_path)
            has_ks = 1 if n_ks > 0 else 0

            total, ok, miss_prot, miss_cds, miss_both = parse_diag(diag_path)

            if has_ks:
                reason = "OK"
            else:
                if total == 0:
                    reason = "NO_DIAGNOSTICS_OR_INPUT"
                elif ok == 0 and total > 0:
                    if miss_cds > 0 and miss_prot == 0 and miss_both == 0:
                        reason = "NO_IDS_WITH_CDS"  # like AP_*: only CDS IDs mismatched
                    elif miss_prot > 0 and miss_cds == 0 and miss_both == 0:
                        reason = "NO_IDS_WITH_PROTEIN"
                    else:
                        reason = "NO_IDS_WITH_BOTH_PROT_AND_CDS"
                else:
                    reason = "UNKNOWN_NO_KS"

            out.write(
                f"{sp}\t{has_ks}\t{n_ks}\t{total}\t{ok}\t{miss_prot}\t"
                f"{miss_cds}\t{miss_both}\t{reason}\n"
            )

    print(f"[status] Wrote species Ks status to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
