#!/usr/bin/env python3
"""
03a_extract_cds_for_ogs.py
--------------------------
Extract CDS nucleotide sequences for each OrthoFinder orthogroup × species pair,
creating a fake busco_runs/ directory structure compatible with HTT step 03.

Output structure (mirrors BUSCO single_copy output expected by 03_ks):
  BUSCO_RUNS_DIR / {SPECIES}_dothideomycetes_odb10 /
      run_dothideomycetes_odb10 / busco_sequences / single_copy_busco_sequences /
          {OG_ID}.fna          ← one CDS sequence for that species

Sources:
  - OrthoFinder MSA:  results/Orthofinder_run/MultipleSequenceAlignments/{OG}.fa
  - GFF files:        data/genes_annot_gff/{SPECIES}.gff[3]
  - Genome FASTAs:    data/genomes/{SPECIES}.fasta (.fai must exist)

Usage:
    micromamba run -n syntenyEDFTE python3 03a_extract_cds_for_ogs.py [--threads N]

Requirements: Python ≥ 3.8 (stdlib only); genome .fasta.fai files must exist.
"""

import os, re, sys, argparse, logging
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# ── Paths ─────────────────────────────────────────────────────────────────────
_ROOT        = Path("/legserv/Temp/Carlos/EDF-TE")
GENOMES_DIR  = _ROOT / "data/genomes"
GFF_DIR      = _ROOT / "data/genes_annot_gff"
MSA_DIR      = _ROOT / "results/Orthofinder_run/MultipleSequenceAlignments"
OG_LIST_FILE = _ROOT / "results/HTT/01_phylo_tree/selected_ogs.txt"
BUSCO_RUNS   = _ROOT / "results/HTT/01_phylo_tree/busco_runs"
LOG_FILE     = Path("/data/carlos/EDF-TE/logs/HTT/03a_cds_extraction.log")

DB_SUFFIX    = "dothideomycetes_odb10"   # matches step 03 glob pattern

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger(__name__)


# ── Sequence utilities ────────────────────────────────────────────────────────
_COMP = str.maketrans("ACGTacgtNnRrYySsWwKkMmBbDdHhVv",
                       "TGCAtgcaNnYyRrSsWwMmKkVvHhDdBb")

def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ── FAI-indexed random access ─────────────────────────────────────────────────
def load_fai(fasta_path: Path):
    """Return dict: chrom -> (length, byte_offset, bpl, byl)."""
    for fai_path in [fasta_path.with_suffix(".fasta.fai"),
                     Path(str(fasta_path) + ".fai"),
                     Path(str(fasta_path) + ".prep.fai"),
                     fasta_path.parent / (fasta_path.stem + ".fasta.fai")]:
        if fai_path.exists():
            idx = {}
            with open(fai_path) as f:
                for line in f:
                    p = line.strip().split("\t")
                    if len(p) >= 5:
                        idx[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
            return idx
    return {}


def get_seq(fasta_path: Path, fai: dict, chrom: str, start: int, end: int) -> str:
    """Extract sequence [start:end) (0-based) using FAI random access."""
    if chrom not in fai:
        return ""
    length, byte_off, bpl, byl = fai[chrom]
    start = max(0, start)
    end   = min(length, end)
    if start >= end:
        return ""
    # byte positions
    b_start = byte_off + (start // bpl) * byl + (start % bpl)
    b_end   = byte_off + ((end - 1) // bpl) * byl + ((end - 1) % bpl) + 1
    n_newlines = (end - 1) // bpl - start // bpl   # newlines between start and end
    with open(fasta_path, "rb") as fh:
        fh.seek(b_start)
        raw = fh.read(b_end - b_start + n_newlines + 4)
    seq = raw.decode("ascii", errors="replace").replace("\n", "").replace("\r", "")
    return seq[: end - start].upper()


def extract_cds(fasta_path: Path, fai: dict, intervals: list, strand: str) -> str:
    """Extract + concatenate CDS exons; reverse-complement if minus strand."""
    if not intervals:
        return ""
    # intervals: list of (chrom, start, end) with 0-based coords
    sorted_ivs = sorted(intervals, key=lambda x: x[1])
    pieces = []
    for chrom, s, e in sorted_ivs:
        pieces.append(get_seq(fasta_path, fai, chrom, s, e))
    cds = "".join(pieces)
    return revcomp(cds) if strand == "-" else cds


# ── GFF parsing ───────────────────────────────────────────────────────────────
def build_gff_index(gff_path: Path) -> dict:
    """
    Return dict: protein_id -> [(chrom, start, end)] with 0-based coords.
    strand stored as dict: protein_id -> '+'/'-'
    Returns (cds_index, strand_index).
    """
    cds_idx    = defaultdict(list)
    strand_idx = {}
    mrna_protid = {}  # mRNA feature ID -> proteinId  (JGI GFF3)

    try:
        fh = open(gff_path)
    except OSError:
        return {}, {}

    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        chrom, _, ftype, start_s, end_s, _, strand, phase_s, attrs = parts
        try:
            start = int(start_s) - 1   # 0-based
            end   = int(end_s)          # exclusive
        except ValueError:
            continue

        if ftype == "mRNA":
            # JGI GFF3: mRNA has proteinId=N
            m = re.search(r"[Pp]rotein[Ii]d[= ](\d+)", attrs)
            if m:
                mid = re.search(r"ID=([^;]+)", attrs)
                if mid:
                    mrna_protid[mid.group(1)] = m.group(1)

        elif ftype == "CDS":
            pid = None

            # NCBI GFF3: protein_id=ACC
            m = re.search(r"protein_id=([^;,\s]+)", attrs)
            if m:
                pid = m.group(1)

            # JGI GFF2: proteinId N (space-separated, old format)
            if pid is None:
                m = re.search(r"proteinId (\d+)", attrs)
                if m:
                    pid = m.group(1)

            # JGI GFF3: Parent=mRNA_ID → look up in mrna_protid
            if pid is None:
                m = re.search(r"Parent=([^;]+)", attrs)
                if m:
                    pid = mrna_protid.get(m.group(1))

            if pid:
                cds_idx[pid].append((chrom, start, end))
                if pid not in strand_idx:
                    strand_idx[pid] = strand

    fh.close()
    return dict(cds_idx), strand_idx


# ── MSA header parsing ────────────────────────────────────────────────────────
def parse_msa_header(header: str):
    """
    Parse OrthoFinder MSA header → (species_code, protein_id).

    Formats encountered:
      >SPECIES_aa_jgi|ShortName|proteinId|geneName   → JGI, with _aa_
      >SPECIES_aa_KAI3657220.1                        → NCBI, with _aa_
      >SPECIES_jgi|ShortName|proteinId|geneName       → JGI, without _aa_
      >SPECIES                                        → species only (step01 processed)
    """
    h = header.lstrip(">").strip()

    # Extract species code
    m = re.match(r"^([A-Z]{2,3}_[A-Za-z0-9]+)", h)
    if not m:
        return None, None
    species = m.group(1)
    rest = h[len(species):]

    # Format: _aa_jgi|ShortName|proteinId|...
    m2 = re.match(r"^_aa_jgi\|[^|]+\|(\d+)\|", rest)
    if m2:
        return species, m2.group(1)

    # Format: _jgi|ShortName|proteinId|... (without _aa_)
    m3 = re.match(r"^_jgi\|[^|]+\|(\d+)\|", rest)
    if m3:
        return species, m3.group(1)

    # Format: _aa_NCBI_ACC (e.g. KAI3657220.1)
    m4 = re.match(r"^_aa_([A-Z0-9_]+\.\d+)", rest)
    if m4:
        return species, m4.group(1)

    # No protein ID parseable
    return species, None


# ── Find GFF file for a species ───────────────────────────────────────────────
def find_gff(species: str) -> Path | None:
    for ext in (".gff3", ".gff", ".gff3.gz", ".gff.gz"):
        p = GFF_DIR / (species + ext)
        if p.exists():
            return p
    return None


def find_genome(species: str) -> Path | None:
    for ext in (".fasta", ".fa", ".fna"):
        p = GENOMES_DIR / (species + ext)
        if p.exists():
            return p
    return None


# ── Main processing ───────────────────────────────────────────────────────────
def process_species(species: str, og_pid_map: dict):
    """
    og_pid_map: {og_id: protein_id}   for this species
    Returns: {og_id: cds_sequence_or_None}
    """
    gff_path    = find_gff(species)
    genome_path = find_genome(species)
    if gff_path is None or genome_path is None:
        return {}

    cds_idx, strand_idx = build_gff_index(gff_path)
    fai = load_fai(genome_path)
    if not fai:
        log.warning(f"  [{species}] FAI not found for {genome_path}")
        return {}

    results = {}
    for og_id, prot_id in og_pid_map.items():
        if prot_id is None or prot_id not in cds_idx:
            continue
        intervals = cds_idx[prot_id]
        strand    = strand_idx.get(prot_id, "+")
        cds = extract_cds(genome_path, fai, intervals, strand)
        if cds and len(cds) >= 3:
            results[og_id] = cds
    return results


def write_busco_fna(species: str, og_id: str, cds: str):
    """Write CDS to busco_runs-compatible path."""
    out_dir = (BUSCO_RUNS / f"{species}_{DB_SUFFIX}"
               / f"run_{DB_SUFFIX}" / "busco_sequences" / "single_copy_busco_sequences")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{og_id}.fna"
    with open(out_file, "w") as f:
        f.write(f">{species}\n{cds}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of parallel threads (default: 4)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing .fna files")
    args = parser.parse_args()

    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(LOG_FILE)
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(file_handler)

    # Load OG list
    ogs = [l.strip() for l in OG_LIST_FILE.read_text().splitlines() if l.strip()]
    log.info(f"Loaded {len(ogs)} orthogroups")

    # Build species → {og_id → protein_id} mapping from original OrthoFinder MSAs
    # structure: sp_og_map[species][og_id] = protein_id
    sp_og_map = defaultdict(dict)
    n_headers_parsed = 0
    n_headers_missing_pid = 0

    for og_id in ogs:
        msa_file = MSA_DIR / f"{og_id}.fa"
        if not msa_file.exists():
            continue
        with open(msa_file) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                species, prot_id = parse_msa_header(line.rstrip())
                if species is None:
                    continue
                n_headers_parsed += 1
                if prot_id is None:
                    n_headers_missing_pid += 1
                # Take the FIRST protein per species per OG (should be 1 in SC OGs)
                if og_id not in sp_og_map[species]:
                    sp_og_map[species][og_id] = prot_id

    log.info(f"Parsed {n_headers_parsed} headers "
             f"({n_headers_missing_pid} missing protein ID)")
    log.info(f"Species to process: {len(sp_og_map)}")

    # Count how many CDS are already written
    already_done = 0
    if not args.overwrite:
        for sp in sp_og_map:
            sc_dir = (BUSCO_RUNS / f"{sp}_{DB_SUFFIX}" / f"run_{DB_SUFFIX}"
                      / "busco_sequences" / "single_copy_busco_sequences")
            if sc_dir.exists():
                already_done += len(list(sc_dir.glob("*.fna")))
    if already_done:
        log.info(f"  {already_done} CDS files already written (use --overwrite to redo)")

    # Process species in parallel
    total_written = 0
    total_failed  = 0
    species_list  = sorted(sp_og_map.keys())

    def _process(sp):
        og_pid = sp_og_map[sp]
        if not args.overwrite:
            sc_dir = (BUSCO_RUNS / f"{sp}_{DB_SUFFIX}" / f"run_{DB_SUFFIX}"
                      / "busco_sequences" / "single_copy_busco_sequences")
            # Skip OGs already extracted
            og_pid = {og: pid for og, pid in og_pid.items()
                      if not (sc_dir / f"{og}.fna").exists()}
        if not og_pid:
            return 0, 0
        cds_results = process_species(sp, og_pid)
        written = 0
        failed  = 0
        for og_id, cds in cds_results.items():
            write_busco_fna(sp, og_id, cds)
            written += 1
        failed = len(og_pid) - written
        return written, failed

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {pool.submit(_process, sp): sp for sp in species_list}
        done = 0
        for fut in as_completed(futures):
            sp = futures[fut]
            done += 1
            try:
                w, f = fut.result()
                total_written += w
                total_failed  += f
                if done % 20 == 0 or done == len(species_list):
                    log.info(f"  [{done}/{len(species_list)}] "
                             f"written={total_written} failed={total_failed}")
            except Exception as exc:
                log.error(f"  [{sp}] ERROR: {exc}")

    log.info(f"Done! Total CDS written: {total_written} | Failed: {total_failed}")
    log.info(f"busco_runs dir: {BUSCO_RUNS}")


if __name__ == "__main__":
    main()
