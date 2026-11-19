#!/usr/bin/env python3
import argparse
import os
import sys
from typing import Set


def read_fasta_ids(path: str) -> Set[str]:
    ids: Set[str] = set()
    try:
        with open(path, "r") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                if not header:
                    continue
                seq_id = header.split()[0]
                ids.add(seq_id)
    except OSError as e:
        print(f"[diag_ids] ERROR: cannot read FASTA {path}: {e}", file=sys.stderr)
        sys.exit(1)
    return ids


def read_pairs_ids(path: str) -> Set[str]:
    ids: Set[str] = set()
    try:
        with open(path, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                ids.add(parts[0])
                ids.add(parts[1])
    except OSError as e:
        print(f"[diag_ids] ERROR: cannot read pairs file {path}: {e}", file=sys.stderr)
        sys.exit(1)
    return ids


def main() -> None:
    p = argparse.ArgumentParser(description="Diagnose ID coverage between pairs, protein FASTA and CDS FASTA for one species.")
    p.add_argument("--species", required=True, help="Species name (for reporting only)")
    p.add_argument("--proteins", required=True, help="Protein FASTA path for the species")
    p.add_argument("--cds", required=True, help="CDS FASTA path for the species")
    p.add_argument("--pairs", required=True, help="TSV with gene1, gene2 (per species)")
    p.add_argument("--out", required=True, help="Output TSV with per-gene coverage diagnostics")
    args = p.parse_args()

    sp = args.species

    prot_ids = read_fasta_ids(args.proteins)
    cds_ids = read_fasta_ids(args.cds)
    pair_ids = read_pairs_ids(args.pairs)

    if not pair_ids:
        print(f"[diag_ids] WARNING: no IDs in pairs file for {sp}; nothing to diagnose.", file=sys.stderr)
        # still write empty diagnostics file
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w") as out_fh:
            out_fh.write("gene\tin_pairs\thas_protein\thas_cds\tstatus\n")
        return

    ok_ids = set()
    missing_prot = set()
    missing_cds = set()
    missing_both = set()

    for gid in pair_ids:
        has_p = gid in prot_ids
        has_c = gid in cds_ids
        if has_p and has_c:
            ok_ids.add(gid)
        elif has_p and not has_c:
            missing_cds.add(gid)
        elif not has_p and has_c:
            missing_prot.add(gid)
        else:
            missing_both.add(gid)

    total_pairs_ids = len(pair_ids)
    print(
        f"[diag_ids] {sp}: unique IDs in pairs={total_pairs_ids}, ok={len(ok_ids)}, "
        f"missing_protein_only={len(missing_prot)}, missing_cds_only={len(missing_cds)}, "
        f"missing_both={len(missing_both)}",
        file=sys.stderr,
    )

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as out_fh:
        out_fh.write("gene\tin_pairs\thas_protein\thas_cds\tstatus\n")
        for gid in sorted(pair_ids):
            has_p = gid in prot_ids
            has_c = gid in cds_ids
            if has_p and has_c:
                status = "OK"
            elif has_p and not has_c:
                status = "MISSING_CDS"
            elif not has_p and has_c:
                status = "MISSING_PROTEIN"
            else:
                status = "MISSING_BOTH"
            out_fh.write(
                f"{gid}\t1\t{int(has_p)}\t{int(has_c)}\t{status}\n"
            )

    # Exit code 0 if at least one ID has both protein and CDS; 2 otherwise
    if ok_ids:
        sys.exit(0)
    else:
        sys.exit(2)


if __name__ == "__main__":
    main()
