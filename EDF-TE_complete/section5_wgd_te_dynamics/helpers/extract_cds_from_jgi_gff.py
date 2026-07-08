#!/usr/bin/env python

import argparse
from collections import defaultdict

from pathlib import Path

from typing import Dict, List, Tuple


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


def load_protein_map(protein_fasta: Path) -> Dict[str, str]:

	"""Map numeric proteinId -> full protein header, for JGI-style FASTA.

	E.g. header:
	  >jgi|ChetHAW225_4|179880|e_gw1.1.255.1
	we map "179880" -> "jgi|ChetHAW225_4|179880|e_gw1.1.255.1".
	"""

	map_id: Dict[str, str] = {}
	with protein_fasta.open() as fh:
		for line in fh:
			if not line.startswith(">"):
				continue
			h = line[1:].strip().split()[0]
			parts = h.split("|")
			if len(parts) >= 3 and parts[2].isdigit():
				map_id[parts[2]] = h
	return map_id


def parse_jgi_gff(path: Path) -> Dict[str, List[Tuple[str, str, int, int, str]]]:

	"""Parse a JGI-style GFF and collect CDS exons per proteinId.

	Returns: proteinId -> list of (seqid, name, start, end, strand).
	Coordinates are 1-based inclusive as in GFF.
	"""

	cds_by_prot: Dict[str, List[Tuple[str, str, int, int, str]]] = defaultdict(list)
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
			# parse attributes like:
			# name "e_gw1.1.255.1"; proteinId 179880; exonNumber 1; ...
			name = None
			prot_id = None
			for field in attrs.split(";"):
				field = field.strip()
				if not field:
					continue
				if field.startswith("name "):
					val = field.split(" ", 1)[1].strip().strip("\"")
					name = val
				elif field.startswith("proteinId "):
					val = field.split(" ", 1)[1].strip().strip("\"")
					prot_id = val
			if prot_id is None:
				continue
			try:
				s = int(start)
				e = int(end)
			except ValueError:
				continue
			cds_by_prot[prot_id].append((seqid, name or prot_id, s, e, strand))
	return cds_by_prot


def write_cds_fasta(
	genome_fasta: Path,
	protein_fasta: Path,
	gff_file: Path,
	out_fasta: Path,
) -> int:

	genome = parse_fasta(genome_fasta)
	prot_map = load_protein_map(protein_fasta)
	cds_by_prot = parse_jgi_gff(gff_file)

	written = 0
	with out_fasta.open("w") as out:
		for prot_id, exons in cds_by_prot.items():
			if prot_id not in prot_map:
				continue
			# order exons along the genome
			exons_sorted = sorted(exons, key=lambda x: x[2])
			seqid, _name, _s, _e, strand = exons_sorted[0]
			if seqid not in genome:
				continue
			seq = genome[seqid]
			nuc_fragments: List[str] = []
			for _seqid2, _name2, s, e, _strand2 in exons_sorted:
				if _seqid2 != seqid:
					# skip mismatched contig just in case
					continue
				# convert to 0-based slice
				start0 = max(s - 1, 0)
				end0 = min(e, len(seq))
				if start0 >= end0:
					continue
				nuc_fragments.append(seq[start0:end0])
			if not nuc_fragments:
				continue
			nuc_seq = "".join(nuc_fragments)
			if strand == "-":
				nuc_seq = revcomp(nuc_seq)
			hdr = prot_map[prot_id]
			out.write(f">{hdr}\n")
			for i in range(0, len(nuc_seq), 80):
				out.write(nuc_seq[i : i + 80] + "\n")
			written += 1
	return written


def main() -> None:
	parser = argparse.ArgumentParser(
		description=(
			"Extract CDS sequences from a JGI-style GFF (using proteinId/name) "
			"and a genome FASTA, writing a CDS FASTA keyed by protein FASTA IDs."
		)
	)
	parser.add_argument("--genome", required=True, help="Genome FASTA file")
	parser.add_argument("--gff", required=True, help="JGI-style GFF file")
	parser.add_argument(
		"--proteins",
		required=True,
		help="Protein FASTA used by OrthoFinder (for mapping proteinId)",
	)
	parser.add_argument("--out-cds", required=True, help="Output CDS FASTA path")
	args = parser.parse_args()

	written = write_cds_fasta(
		Path(args.genome), Path(args.proteins), Path(args.gff), Path(args.out_cds)
	)
	print(
		f"[WGD][Ks] JGI-GFF extractor wrote {written} CDS sequences to {args.out_cds}",
	)


if __name__ == "__main__":
	main()
