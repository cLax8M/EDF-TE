#!/usr/bin/env python3
"""
03b_fimo_cross_species_per_sp.py
─────────────────────────────────
Per-species FIMO scan for key TF motif × TE class associations
identified from the Rhizopus microsporus enrichment analysis.

For each (TE class, motif) pair, this script:
  1. Filters the TE BED for that class (≥ MIN_TE_LEN)
  2. Extracts FASTA with bedtools getfasta
  3. Runs FIMO to find motif occurrences
  4. Counts the fraction of elements with ≥1 motif hit

Output: one TSV row per (species, te_class, motif)
Columns: species, te_class, motif, n_elements, n_with_motif, fraction_with_motif

Usage
  python 03b_fimo_cross_species_per_sp.py \\
      --species MU_Rhimi59 \\
      --bed /path/to/MU_Rhimi59.filteredRepeats.bed \\
      --genome /path/to/MU_Rhimi59.fasta \\
      --motif_dir /path/to/motifs \\
      --outdir /path/to/output \\
      --out_tsv /path/to/cross_species_results.tsv

Requirements:
  module load BEDTools/2.31.1-GCC-13.3.0
  module load MEME/5.5.7-gompi-2023b
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

# ── Key TE class × motif pairs to scan ─────────────────────────────────────
# Based on significant enrichments in MU_Rhimi59 (fimo_enrichment.tsv)
TE_MOTIF_PAIRS = [
    # LINE / HMGb associations — the main cross-species hypothesis
    ("LINE/I-Jockey", "HMGb3"),
    ("LINE/I-Jockey", "HMGb2"),
    ("LINE/R1",        "HMGb3"),
    ("LINE/R1",        "HMGb2"),
    ("LINE",           "HMGb2"),
    ("LINE",           "HMGb3"),
    ("LINE/I",         "HMGb2"),
    ("LINE/I",         "HMGb3"),
    ("LINE/RTE-X",     "HMGb2"),
    ("LINE/RTE-X",     "HMGb3"),
    ("LINE/L1",        "HMGb3"),
    ("LINE/L1",        "HMGb2"),
    # LTR / GATA_Znf associations
    ("LTR/Gypsy",      "GATA_Znf9"),
    ("LTR/Gypsy",      "GATA_Znf10"),
    ("LTR/DIRS",       "GATA_Znf10"),
    ("LTR/Ngaro",      "HMGb3"),
    # DNA / HLH_Myc associations
    ("DNA/CMC-EnSpm",  "HLH_Myc6"),
    ("DNA/CMC-EnSpm",  "HMGb2"),
    ("DNA/Merlin",     "bZip3"),
    ("DNA/Merlin",     "HLH_Myc6"),
    ("DNA/Crypton-F",  "HLH_Myc15"),
    # SINE
    ("SINE/Alu",       "HMGb3"),
]

MIN_TE_LEN  = 200    # bp minimum element length
FIMO_THRESH = 1e-4   # p-value threshold for FIMO hits
MODULE_LOAD = (
    "source /etc/profile.d/lmod.sh 2>/dev/null; "
    "module load BEDTools/2.31.1-GCC-13.3.0; "
    "module load MEME/5.5.7-gompi-2023b"
)

SKIP_TE_CLASSES = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA",
                   "snRNA", "Unknown", "Unspecified", "Other"}


def run_cmd(cmd: str, desc: str = "") -> str:
    full = f"bash -c '{MODULE_LOAD}; {cmd}'"
    r = subprocess.run(full, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Command failed [{desc}]:\n{r.stderr[:600]}")
    return r.stdout


def safe_name(s: str) -> str:
    return s.replace("/", "_").replace(" ", "_")


def extract_fasta_for_class(te_bed: Path, te_class: str, genome_fa: Path,
                             work_dir: Path, min_len: int) -> tuple[Path | None, int]:
    """
    Extract FASTA sequences for a given TE class.
    Returns (fasta_path, n_elements). fasta_path is None if no elements found.
    """
    sname = safe_name(te_class)
    bed_path = work_dir / f"te_{sname}.bed"
    fa_path  = work_dir / f"te_{sname}.fa"

    # Write per-class BED
    entries = []
    with open(te_bed) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            cls = parts[3].strip()
            if cls == te_class and (end - start) >= min_len:
                entries.append(f"{chrom}\t{start}\t{end}\t{cls}\n")

    if not entries:
        return None, 0

    with open(bed_path, "w") as fout:
        fout.writelines(entries)

    # Extract FASTA
    if not fa_path.exists() or fa_path.stat().st_size == 0:
        run_cmd(
            f"bedtools getfasta -fi {genome_fa} -bed {bed_path} -fo {fa_path} -name",
            f"getfasta {te_class}"
        )

    n_seqs = sum(1 for line in open(fa_path) if line.startswith(">"))
    return fa_path if n_seqs > 0 else None, n_seqs


def run_fimo_and_count(motif_file: Path, fasta_path: Path,
                       fimo_out_dir: Path, n_total: int,
                       threshold: float = FIMO_THRESH) -> tuple[int, int]:
    """
    Run FIMO and return (n_seqs_with_hit, n_total_seqs).
    """
    fimo_out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = fimo_out_dir / "fimo.tsv"

    if not (tsv_path.exists() and tsv_path.stat().st_size > 0):
        run_cmd(
            f"fimo --thresh {threshold} --oc {fimo_out_dir} "
            f"{motif_file} {fasta_path}",
            f"fimo {motif_file.stem}"
        )

    if not tsv_path.exists():
        return 0, n_total

    seqs_hit: set[str] = set()
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("motif_id"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                seqs_hit.add(parts[2])
    return len(seqs_hit), n_total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species",   required=True, help="Species identifier")
    ap.add_argument("--bed",       required=True, help="TE BED file")
    ap.add_argument("--genome",    required=True, help="Genome FASTA")
    ap.add_argument("--motif_dir", required=True, help="Directory with .meme motif files")
    ap.add_argument("--outdir",    required=True, help="Working directory for this species")
    ap.add_argument("--out_tsv",   required=True, help="Output summary TSV (append mode)")
    args = ap.parse_args()

    species    = args.species
    te_bed     = Path(args.bed)
    genome_fa  = Path(args.genome)
    motif_dir  = Path(args.motif_dir)
    work_dir   = Path(args.outdir) / species
    out_tsv    = Path(args.out_tsv)
    work_dir.mkdir(parents=True, exist_ok=True)

    # Verify inputs
    if not te_bed.exists():
        print(f"[ERROR] BED not found: {te_bed}", file=sys.stderr)
        sys.exit(1)
    if not genome_fa.exists():
        print(f"[ERROR] Genome not found: {genome_fa}", file=sys.stderr)
        sys.exit(1)

    # Index genome if needed
    fai = Path(str(genome_fa) + ".fai")
    if not fai.exists():
        run_cmd(f"samtools faidx {genome_fa}", "samtools faidx")

    # Collect the TE classes present in this species
    present_classes: set[str] = set()
    with open(te_bed) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                present_classes.add(p[3].strip())

    # FASTA cache per class (avoid re-extracting for same class × multiple motifs)
    fasta_cache: dict[str, tuple[Path | None, int]] = {}

    rows = []
    for (te_class, motif_name) in TE_MOTIF_PAIRS:
        if te_class not in present_classes:
            continue

        # Get motif file
        motif_file = motif_dir / f"{motif_name}.meme"
        if not motif_file.exists():
            # Try alternate capitalisation
            alts = list(motif_dir.glob(f"{motif_name}*.meme"))
            if not alts:
                print(f"[WARN] Motif not found: {motif_name}", file=sys.stderr)
                continue
            motif_file = alts[0]

        # Extract FASTA for this TE class (cached)
        if te_class not in fasta_cache:
            fasta_cache[te_class] = extract_fasta_for_class(
                te_bed, te_class, genome_fa, work_dir, MIN_TE_LEN
            )
        fasta_path, n_elements = fasta_cache[te_class]

        if fasta_path is None or n_elements == 0:
            continue

        # FIMO output directory per (species, te_class, motif)
        fimo_out = work_dir / "fimo" / safe_name(te_class) / motif_name
        n_with_motif, n_total = run_fimo_and_count(
            motif_file, fasta_path, fimo_out, n_elements
        )
        fraction = n_with_motif / n_total if n_total > 0 else 0.0

        rows.append((species, te_class, motif_name, n_total, n_with_motif,
                     f"{fraction:.4f}"))
        print(f"[{species}] {te_class} × {motif_name}: "
              f"{n_with_motif}/{n_total} = {fraction:.3f}")

    # Append rows to shared TSV (use lock file to avoid race conditions)
    if rows:
        lock = Path(str(out_tsv) + ".lock")
        import time
        # Simple spin-lock
        for _ in range(60):
            try:
                fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.close(fd)
                break
            except FileExistsError:
                time.sleep(1)

        try:
            write_header = not out_tsv.exists() or out_tsv.stat().st_size == 0
            with open(out_tsv, "a") as fout:
                if write_header:
                    fout.write("species\tte_class\tmotif\tn_elements\t"
                               "n_with_motif\tfraction_with_motif\n")
                for row in rows:
                    fout.write("\t".join(str(x) for x in row) + "\n")
        finally:
            try:
                lock.unlink()
            except FileNotFoundError:
                pass

    print(f"[DONE] {species} — {len(rows)} rows written")


if __name__ == "__main__":
    main()
