#!/usr/bin/env python3

#Extract CDS sequences from an NCBI *cds_from_genomic.fna* file.

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, Set


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build CDS FASTA keyed by protein_id from NCBI cds_from_genomic file",
    )
    p.add_argument(
        "--cds-from-genomic",
        required=True,
        help="Path to GCF_xxx_cds_from_genomic.fna (decompressed)",
    )
    p.add_argument(
        "--proteins",
        required=True,
        help="Proteome FASTA used in OrthoFinder (e.g. MU_Rhior3.aa.fasta)",
    )
    p.add_argument(
        "--out-cds",
        required=True,
        help=(
            "Output CDS FASTA; headers will be protein IDs that also occur "
            "in the proteome FASTA"
        ),
    )
    return p.parse_args()


def load_protein_ids(prot_fasta: Path) -> Set[str]:
    ids: Set[str] = set()
    with prot_fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    continue
                prot_id = header.split()[0]
                ids.add(prot_id)
    return ids


def iter_fasta_records(path: Path) -> Iterable[tuple[str, str]]:
    """Yield (header, seq) for each FASTA record (header without '>')."""

    with path.open() as fh:
        header: str | None = None
        seq_chunks: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def extract_protein_id(cds_header: str) -> str | None:
    """Parse protein_id from an NCBI cds_from_genomic header.

    Example header:

        lcl|NW_027137854.1_cds_1 [gene=RO3G_17256] ... [protein_id=XP_067527781.1] ...
    """

    tag = "protein_id="
    idx = cds_header.find(tag)
    if idx == -1:
        return None

    sub = cds_header[idx + len(tag) :]
    # Typically ends at the closing bracket or whitespace
    # e.g. "XP_067527781.1] [location=..." -> take up to ']' or space
    for sep in "] ":
        sub = sub.split(sep, 1)[0]
    return sub.strip()


def build_protein_to_cds(
    cds_from_genomic: Path, wanted_proteins: Set[str]
) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for header, seq in iter_fasta_records(cds_from_genomic):
        prot_id = extract_protein_id(header)
        if prot_id is None:
            continue
        if prot_id in wanted_proteins and prot_id not in mapping:
            mapping[prot_id] = seq
    return mapping


def write_cds_fasta(out_path: Path, prot_to_cds: Dict[str, str]) -> None:
    with out_path.open("w") as out:
        for prot_id, seq in sorted(prot_to_cds.items()):
            out.write(f">{prot_id}\n")
            # Wrap sequences at 60–80 chars per line for readability
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")


def main() -> None:
    args = parse_args()

    cds_from_genomic = Path(args.cds_from_genomic).resolve()
    proteins = Path(args.proteins).resolve()
    out_cds = Path(args.out_cds).resolve()

    if not cds_from_genomic.is_file():
        raise SystemExit(f"cds_from_genomic file not found: {cds_from_genomic}")
    if not proteins.is_file():
        raise SystemExit(f"Proteome FASTA not found: {proteins}")

    wanted_ids = load_protein_ids(proteins)
    prot_to_cds = build_protein_to_cds(cds_from_genomic, wanted_ids)

    if not prot_to_cds:
        raise SystemExit(
            "No overlapping protein IDs found between cds_from_genomic and proteome."
        )

    out_cds.parent.mkdir(parents=True, exist_ok=True)
    write_cds_fasta(out_cds, prot_to_cds)

    import sys

    sys.stderr.write(
        f"[WGD][Ks] Wrote {len(prot_to_cds)} CDS sequences to {out_cds} "
        f"(from {cds_from_genomic.name}, matched to {proteins.name})\n"
    )


if __name__ == "__main__":
    main()
