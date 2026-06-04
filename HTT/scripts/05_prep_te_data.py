#!/usr/bin/env python3
"""
05_prep_te_data.py
Build step-05 outputs for the HTT pipeline from EarlGrey BED + genome FASTAs.

Produces:
  <OUT_DIR>/all_tes_combined.fasta        – all TE sequences, header = {genome}|{te_id}
  <OUT_DIR>/all_tes_info.tsv              – te_id, genome, chr, start, end, strand, family, class
  <OUT_DIR>/filtered_tes/<genome>_filtered_tes.tsv  – stub domain files (real TEs ≥ MIN_BP)

TE ids inside sequences and domain files: {chr}:{start}-{end}:{strand}
Full ids (as in FASTA header and info TSV): {genome}|{chr}:{start}-{end}:{strand}

Usage:
  python3 05_prep_te_data.py [--out_dir DIR] [--min_bp INT] [--threads INT]
"""

import os
import sys
import gzip
import argparse
import subprocess
import concurrent.futures
from pathlib import Path
from collections import defaultdict

# ── Paths ────────────────────────────────────────────────────────────────────
TE_ANNOT_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/HTT/04_te_annotation")
GENOMES_DIR  = Path("/legserv/Temp/Carlos/EDF-TE/data/genomes")

# TE classes considered real autonomous/non-repetitive TEs (allowed in domain stubs)
REAL_TE_CLASSES = {
    "LTR", "LINE", "SINE", "DNA", "RC", "Penelope",
    "Unknown", "Unclassified",
}

# Classes to EXCLUDE entirely (even from all_tes_combined.fasta)
SKIP_CLASSES = {
    "Simple_repeat", "Satellite", "Low_complexity", "rRNA",
    "snRNA", "srpRNA", "tRNA", "scRNA", "RNA", "ARTEFACT",
    "Tandem_repeat",
}


def parse_bed(bed_path: Path):
    """Yield (chr, start0, end, family, score, strand, te_class) tuples."""
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chr_    = parts[0]
            start0  = int(parts[1])   # 0-based
            end     = int(parts[2])   # 1-based end
            family  = parts[3]        # e.g. LTR/Gypsy or Simple_repeat
            score   = parts[4] if len(parts) > 4 else "."
            strand  = parts[5] if len(parts) > 5 else "+"
            # Derive class from family
            te_class = family.split("/")[0] if "/" in family else family
            yield chr_, start0, end, family, score, strand, te_class


def read_genome_index(fai_path: Path):
    """Return dict chr -> (length, byte_offset, bases_per_line, bytes_per_line)."""
    idx = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chr_ = parts[0]
            length = int(parts[1])
            offset = int(parts[2])
            bpl    = int(parts[3])   # bases per line
            bypl   = int(parts[4])   # bytes per line
            idx[chr_] = (length, offset, bpl, bypl)
    return idx


def extract_seq(genome_fh, idx, chr_, start0, end, strand):
    """Extract subsequence [start0, end) from indexed genome file, revcomp if strand==-."""
    if chr_ not in idx:
        return None
    length, offset, bpl, bypl = idx[chr_]
    end = min(end, length)
    if end <= start0:
        return None

    n_before = start0
    full_lines_before = n_before // bpl
    remainder_before  = n_before % bpl
    seek_pos = offset + full_lines_before * bypl + remainder_before

    genome_fh.seek(seek_pos)
    seq_len = end - start0
    chars_read = 0
    seq_parts = []
    while chars_read < seq_len:
        chunk = genome_fh.read(512)
        if not chunk:
            break
        for c in chunk:
            # In Python 3, iterating bytes gives ints; convert to char
            ch = chr(c) if isinstance(c, int) else c
            if ch in "ACGTNacgtnWSRYKMBDHVwsrykmbdhv":
                seq_parts.append(ch)
                chars_read += 1
                if chars_read >= seq_len:
                    break

    seq = "".join(seq_parts).upper()

    if strand == "-":
        comp = str.maketrans("ACGTN", "TGCAN")
        seq = seq.translate(comp)[::-1]

    return seq


def process_species(species: str, out_dir: Path, min_bp: int):
    """
    Extract TE sequences for one species, return list of:
      (full_te_id, chr_, start0, end, strand, family, te_class, seq)
    """
    bed_path    = TE_ANNOT_DIR / f"{species}.filteredRepeats.bed"
    genome_path = GENOMES_DIR / f"{species}.fasta"
    fai_path    = GENOMES_DIR / f"{species}.fasta.fai"

    if not bed_path.exists():
        print(f"  SKIP {species}: no BED file", flush=True)
        return []
    if not genome_path.exists():
        print(f"  SKIP {species}: no genome FASTA", flush=True)
        return []
    if not fai_path.exists():
        # Try to create index
        try:
            subprocess.run(["samtools", "faidx", str(genome_path)], check=True,
                           capture_output=True)
        except Exception:
            print(f"  SKIP {species}: cannot index genome", flush=True)
            return []

    idx = read_genome_index(fai_path)
    records = []
    n_skip = 0

    with open(genome_path, "rb") as gfh:
        for chr_, start0, end, family, score, strand, te_class in parse_bed(bed_path):
            # Skip non-TE classes
            if str(te_class) in SKIP_CLASSES:
                n_skip += 1
                continue
            length = end - start0
            if length < 1:
                n_skip += 1
                continue

            seq = extract_seq(gfh, idx, chr_, start0, end, strand)
            if seq is None or len(seq) < 1:
                n_skip += 1
                continue

            te_id      = f"{chr_}:{start0}-{end}:{strand}"
            full_te_id = f"{species}|{te_id}"
            records.append((full_te_id, te_id, chr_, start0, end, strand, family, te_class, seq))

    print(f"  {species}: {len(records)} TEs extracted ({n_skip} skipped)", flush=True)
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", default="/legserv/Temp/Carlos/EDF-TE/results/HTT/05_te_classification")
    parser.add_argument("--min_bp",  type=int, default=300,
                        help="Min TE length (bp) to include in domain stub files")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    domains_dir = out_dir / "filtered_tes"
    out_dir.mkdir(parents=True, exist_ok=True)
    domains_dir.mkdir(exist_ok=True)

    # Collect species list
    species_list = sorted({
        p.stem.replace(".filteredRepeats", "")
        for p in TE_ANNOT_DIR.glob("*.filteredRepeats.bed")
    })
    print(f"Found {len(species_list)} species in {TE_ANNOT_DIR}", flush=True)

    fasta_out  = out_dir / "all_tes_combined.fasta"
    info_out   = out_dir / "all_tes_info.tsv"

    total_seqs  = 0
    total_domains = 0

    with open(fasta_out, "w") as fa_fh, open(info_out, "w") as info_fh:
        info_fh.write("te_id\tgenome\tchr\tstart\tend\tstrand\tfamily\tte_class\n")

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = {
                pool.submit(process_species, sp, out_dir, args.min_bp): sp
                for sp in species_list
            }
            for fut in concurrent.futures.as_completed(futures):
                sp = futures[fut]
                try:
                    records = fut.result()
                except Exception as e:
                    print(f"  ERROR {sp}: {e}", flush=True)
                    records = []

                if not records:
                    continue

                # Domain stub file for this species
                dom_rows = []

                for (full_te_id, te_id, chr_, start0, end, strand, family, te_class, seq) in records:
                    # Write FASTA
                    fa_fh.write(f">{full_te_id}\n")
                    for i in range(0, len(seq), 80):
                        fa_fh.write(seq[i:i+80] + "\n")
                    total_seqs += 1

                    # Write info
                    info_fh.write(f"{full_te_id}\t{sp}\t{chr_}\t{start0}\t{end}\t{strand}\t{family}\t{te_class}\n")

                    # Domain stub: only for real TEs ≥ min_bp
                    te_len = end - start0
                    if te_len >= args.min_bp and te_class in REAL_TE_CLASSES:
                        dom_rows.append(
                            f"{te_id}\tRepeatMasker\t{family}\t1\t{te_len}\t1\t{te_len}\t{strand}\n"
                        )
                        total_domains += 1

                # Write domain file for this species
                if dom_rows:
                    dom_file = domains_dir / f"{sp}_filtered_tes.tsv"
                    with open(dom_file, "w") as dom_fh:
                        dom_fh.write("te_id\tsource\tdomain\tte_nt_start\tte_nt_end\tprofile_start\tprofile_end\tstrand\n")
                        dom_fh.writelines(dom_rows)

    print(f"\nDone.", flush=True)
    print(f"  Total TE sequences: {total_seqs:,}", flush=True)
    print(f"  Domain stub entries: {total_domains:,}", flush=True)
    print(f"  FASTA: {fasta_out}", flush=True)
    print(f"  Info:  {info_out}", flush=True)
    print(f"  Domains: {domains_dir}/", flush=True)


if __name__ == "__main__":
    main()
