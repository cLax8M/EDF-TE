#!/usr/bin/env python3
import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Set


def parse_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    current = None
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs.setdefault(current, [])
            else:
                if current is not None:
                    seqs[current].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def load_protein_ids(protein_fasta: Path) -> Set[str]:
    ids: Set[str] = set()
    with protein_fasta.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            hdr = line[1:].strip().split()[0]
            if hdr:
                ids.add(hdr)
    return ids


def parse_ncbi_gff(path: Path, wanted_proteins: Set[str]) -> Dict[str, List[Tuple[str, int, int, str]]]:
    """Parse an NCBI-style GFF and collect CDS exons per protein_id.

    Returns: protein_id -> list of (seqid, start, end, strand).
    Coordinates are 1-based inclusive as in GFF.
    """

    cds_by_prot: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attrs = parts
            if feature != "CDS":
                continue
            prot_id = None
            for field in attrs.split(";"):
                field = field.strip()
                if not field:
                    continue
                if field.startswith("protein_id="):
                    # e.g. protein_id=KAI3655558.1
                    prot_id = field.split("=", 1)[1].strip()
                    break
            if prot_id is None:
                continue
            if prot_id not in wanted_proteins:
                # Only build CDS for proteins present in the OrthoFinder proteome
                continue
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            cds_by_prot[prot_id].append((seqid, s, e, strand))
    return cds_by_prot


def write_cds_fasta_ncbi(
    genome_fasta: Path,
    protein_fasta: Path,
    gff_file: Path,
    out_fasta: Path,
) -> int:
    genome = parse_fasta(genome_fasta)
    prot_ids = load_protein_ids(protein_fasta)
    cds_by_prot = parse_ncbi_gff(gff_file, prot_ids)

    written = 0
    with out_fasta.open("w") as out:
        for prot_id, exons in cds_by_prot.items():
            # order exons along the genome by start coordinate
            exons_sorted = sorted(exons, key=lambda x: x[1])
            seqid, _s0, _e0, strand = exons_sorted[0]
            if seqid not in genome:
                continue
            seq = genome[seqid]
            nuc_frags: List[str] = []
            for seqid2, s, e, strand2 in exons_sorted:
                if seqid2 != seqid:
                    continue
                start0 = max(s - 1, 0)
                end0 = min(e, len(seq))
                if start0 >= end0:
                    continue
                nuc_frags.append(seq[start0:end0])
            if not nuc_frags:
                continue
            nuc_seq = "".join(nuc_frags)
            if strand == "-":
                nuc_seq = revcomp(nuc_seq)
            out.write(f">{prot_id}\n")
            for i in range(0, len(nuc_seq), 80):
                out.write(nuc_seq[i : i + 80] + "\n")
            written += 1
    return written


def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Extract CDS sequences from an NCBI-style GFF (using protein_id) "
            "and a genome FASTA, writing a CDS FASTA keyed by protein IDs that "
            "match the OrthoFinder proteome."
        )
    )
    p.add_argument("--genome", required=True, help="Genome FASTA file")
    p.add_argument("--gff", required=True, help="NCBI-style GFF/GFF3 file")
    p.add_argument(
        "--proteins",
        required=True,
        help="Protein FASTA used by OrthoFinder (for protein_id set)",
    )
    p.add_argument("--out-cds", required=True, help="Output CDS FASTA path")
    args = p.parse_args()

    genome = Path(args.genome)
    gff = Path(args.gff)
    proteins = Path(args.proteins)
    out_cds = Path(args.out_cds)

    if not genome.is_file():
        raise SystemExit(f"Genome FASTA not found: {genome}")
    if not gff.is_file():
        raise SystemExit(f"GFF file not found: {gff}")
    if not proteins.is_file():
        raise SystemExit(f"Protein FASTA not found: {proteins}")

    out_cds.parent.mkdir(parents=True, exist_ok=True)
    written = write_cds_fasta_ncbi(genome, proteins, gff, out_cds)

    import sys

    sys.stderr.write(
        f"[WGD][Ks] NCBI-GFF extractor wrote {written} CDS sequences to {out_cds}\n"
    )


if __name__ == "__main__":
    main()
