#!/usr/bin/env python3
import argparse
import os
import sys


def read_fasta(path):
    seqs = {}
    try:
        with open(path, "r") as fh:
            header = None
            chunks = []
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        seqs[header] = "".join(chunks)
                    header = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if header is not None:
                seqs[header] = "".join(chunks)
    except OSError as e:
        print(f"[prepare_paml_pair_fastas] ERROR: cannot read FASTA {path}: {e}", file=sys.stderr)
        sys.exit(1)
    return seqs


def main():
    parser = argparse.ArgumentParser(description="Prepare per-pair protein/CDS FASTA files for PAML/yn00.")
    parser.add_argument("--proteins", required=True, help="Protein FASTA for the species")
    parser.add_argument("--cds", required=True, help="CDS FASTA for the species")
    parser.add_argument("--pairs", required=True, help="TSV with two columns: gene1<TAB>gene2")
    parser.add_argument("--out-dir", required=True, help="Output directory for per-pair FASTAs")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    prot_seqs = read_fasta(args.proteins)
    cds_seqs = read_fasta(args.cds)

    mapping_path = os.path.join(args.out_dir, "pairs_index.tsv")
    skipped = 0
    total = 0

    try:
        with open(args.pairs, "r") as pf, open(mapping_path, "w") as mf:
            mf.write("pair_id\tgene1\tgene2\n")
            idx = 0
            for line in pf:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                g1, g2 = parts[0], parts[1]
                total += 1
                pep1 = prot_seqs.get(g1)
                pep2 = prot_seqs.get(g2)
                cds1 = cds_seqs.get(g1)
                cds2 = cds_seqs.get(g2)
                if pep1 is None or pep2 is None or cds1 is None or cds2 is None:
                    skipped += 1
                    print(f"[prepare_paml_pair_fastas] WARNING: missing sequence(s) for pair {g1},{g2}; skipping.", file=sys.stderr)
                    continue
                idx += 1
                pair_id = f"pair_{idx:06d}"
                pep_path = os.path.join(args.out_dir, f"{pair_id}.pep.fa")
                cds_path = os.path.join(args.out_dir, f"{pair_id}.cds.fa")
                with open(pep_path, "w") as fh:
                    fh.write(f">{g1}\n{pep1}\n>{g2}\n{pep2}\n")
                with open(cds_path, "w") as fh:
                    fh.write(f">{g1}\n{cds1}\n>{g2}\n{cds2}\n")
                mf.write(f"{pair_id}\t{g1}\t{g2}\n")
    except OSError as e:
        print(f"[prepare_paml_pair_fastas] ERROR: cannot process pairs file {args.pairs}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"[prepare_paml_pair_fastas] Processed {total} pairs, wrote {total - skipped} pairs, skipped {skipped} due to missing sequences.", file=sys.stderr)


if __name__ == "__main__":
    main()
