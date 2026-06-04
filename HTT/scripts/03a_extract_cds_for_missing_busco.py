#!/usr/bin/env python3
"""
03a_extract_cds_for_missing_busco.py

For 6 species missing BUSCO nucleotide sequence files, extract CDS from
GFF+genome using OrthoFinder protein assignments.
Writes fake BUSCO nucleotide files at the paths expected by step 03 Ks calculation.

Missing species:
  AS_AspfumUD1, HO_Homosap, HO_Magisepte_c, MI_MicG_I_3, MU_Rhior3, NU_Parvatl_c

Output: {BUSCO_RUNS_DIR}/{SP}_dothideomycetes_odb10/run_dothideomycetes_odb10/
           busco_sequences/single_copy_busco_sequences/{OG_ID}.fna
"""

import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# ─── Configuration ────────────────────────────────────────────────────────────
RESULTS_ROOT   = Path('/legserv/Temp/Carlos/EDF-TE/results/HTT')
DATA_ROOT      = Path('/legserv/Temp/Carlos/EDF-TE/data')
SEQUENCES_DIR  = RESULTS_ROOT / '01_phylo_tree/sequences'
BUSCO_RUNS_DIR = RESULTS_ROOT / '01_phylo_tree/busco_runs'
ORTHOGROUPS_TSV = Path(
    '/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups/Orthogroups.tsv'
)
OG_LIST_FILE   = RESULTS_ROOT / 'orthofinder_near_universal_80pct_exactly1.txt'
GFF_DIR        = DATA_ROOT / 'genes_annot_gff'
GENOME_DIR     = DATA_ROOT / 'genomes'
BUSCO_LINEAGE  = 'dothideomycetes_odb10'

MISSING_SPECIES = [
    'AS_AspfumUD1',
    'HO_Homosap',
    'HO_Magisepte_c',
    'MI_MicG_I_3',
    'MU_Rhior3',
    'NU_Parvatl_c',
]

# GFF format per species: 'NCBI', 'JGI', or 'AUGUSTUS'
# and the file extension (some species use .gff instead of .gff3)
SPECIES_CONFIG = {
    'AS_AspfumUD1':   {'gff_fmt': 'JGI',      'gff_ext': 'gff3'},
    'HO_Homosap':     {'gff_fmt': 'NCBI',     'gff_ext': 'gff3'},
    'HO_Magisepte_c': {'gff_fmt': 'AUGUSTUS', 'gff_ext': 'gff3'},
    'MI_MicG_I_3':    {'gff_fmt': 'JGI',      'gff_ext': 'gff3'},
    'MU_Rhior3':      {'gff_fmt': 'NCBI',     'gff_ext': 'gff'},
    'NU_Parvatl_c':   {'gff_fmt': 'AUGUSTUS', 'gff_ext': 'gff3'},
}

# ─── Helper functions ─────────────────────────────────────────────────────────
_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def parse_gff_attrs(attr_str: str) -> dict:
    attrs = {}
    for field in attr_str.rstrip(';\n').split(';'):
        field = field.strip()
        if '=' in field:
            k, v = field.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


def extract_protein_id(ortho_id: str, gff_fmt: str) -> str:
    """Return the ID that can be found in GFF attributes."""
    ortho_id = ortho_id.strip()
    if gff_fmt == 'JGI':
        # e.g. 'jgi|AspfumUD1_1|662927|estExt_...' → '662927'
        parts = ortho_id.split('|')
        if len(parts) >= 3:
            return parts[2]
    # NCBI: 'NP_006018.4' or 'XP_067511036.1' → return as-is
    # AUGUSTUS: 'g72696.t1' → return as-is
    return ortho_id.split()[0]


def build_gff_cds_index(gff_path: Path, gff_fmt: str) -> dict:
    """
    Parse GFF and build:
      protein_id → [(chrom, start0, end1, strand), ...]
    Coordinates are 0-based, end is exclusive (BED convention).
    """
    index = defaultdict(list)
    mrna_to_protein = {}  # mRNA_ID → protein_ID  (used for JGI and AUGUSTUS)

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, _src, feat, start, end, _sc, strand, _fr, attrs_str = parts
            attrs  = parse_gff_attrs(attrs_str)
            start0 = int(start) - 1   # GFF3 1-based → 0-based
            end1   = int(end)          # already exclusive in 0-based coords

            if gff_fmt == 'NCBI':
                if feat == 'CDS' and 'protein_id' in attrs:
                    index[attrs['protein_id']].append((chrom, start0, end1, strand))

            elif gff_fmt == 'JGI':
                if feat in ('mRNA', 'transcript'):
                    mrna_id  = attrs.get('ID', '')
                    prot_id  = attrs.get('proteinId', '')
                    if mrna_id and prot_id:
                        mrna_to_protein[mrna_id] = prot_id
                elif feat == 'CDS':
                    parent = attrs.get('Parent', '')
                    if parent in mrna_to_protein:
                        index[mrna_to_protein[parent]].append(
                            (chrom, start0, end1, strand)
                        )

            elif gff_fmt == 'AUGUSTUS':
                if feat in ('mRNA', 'transcript'):
                    tid = attrs.get('ID', '')
                    if tid:
                        mrna_to_protein[tid] = tid  # transcript ID == protein ID
                elif feat == 'CDS':
                    parent = attrs.get('Parent', '')
                    if parent in mrna_to_protein:
                        index[mrna_to_protein[parent]].append(
                            (chrom, start0, end1, strand)
                        )

    return index


SAMTOOLS = '/opt/easybuild/software/SAMtools/1.22.1-GCC-14.2.0/bin/samtools'


def ensure_fai(genome_path: Path) -> None:
    fai = Path(str(genome_path) + '.fai')
    if not fai.exists():
        print(f'  Indexing {genome_path.name} …', flush=True)
        subprocess.run([SAMTOOLS, 'faidx', str(genome_path)], check=True)


def build_chrom_name_map(genome_path: Path, gff_path: Path) -> dict:
    """
    Build mapping: gff_chrom_name → fasta_seq_name.

    Handles three cases:
     1. Exact match (most species).
     2. GFF chrom has underscore-joined long NCBI name, FASTA uses short accession
        (e.g. HO_Magisepte_c, NU_Parvatl_c).
     3. GFF uses NCBI RefSeq accessions (NC_/NW_), FASTA uses different names
        (e.g. HO_Homosap: NC_000001.11 → 1, MU_Rhior3: NW_xxx → supercontig).
        Resolved by matching sequence lengths from ##sequence-region comments vs .fai.
    """
    fai_path = Path(str(genome_path) + '.fai')
    # fasta_name → length
    fai_lengths: dict[str, int] = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.split('\t')
            if len(parts) >= 2:
                fai_lengths[parts[0]] = int(parts[1])

    # length → fasta_name (for length-based matching; assumes unique lengths)
    length_to_fasta: dict[int, str] = {}
    for name, length in fai_lengths.items():
        length_to_fasta[length] = name

    # Single-pass: collect ##sequence-region lengths and all unique GFF chroms
    # NOTE: NCBI INSDC format has ##sequence-region scattered (one per chromosome).
    gff_region_lengths: dict[str, int] = {}
    gff_chroms: set[str] = set()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('##sequence-region'):
                parts = line.split()
                if len(parts) >= 4:
                    gff_region_lengths[parts[1]] = int(parts[3]) - int(parts[2]) + 1
            elif not line.startswith('#'):
                gff_chroms.add(line.split('\t', 1)[0])

    chrom_map: dict[str, str] = {}

    # Determine all unique GFF chroms used in the file
    for gc in gff_chroms:
        # Case 1: exact match
        if gc in fai_lengths:
            chrom_map[gc] = gc
            continue
        # Case 2: strip trailing long NCBI description
        #  e.g. "CM068367.1_Magicicada_..." → first token is "CM068367.1"
        short = gc.split('_')[0]
        if short in fai_lengths:
            chrom_map[gc] = short
            continue
        # Case 3: length-based match
        gc_len = gff_region_lengths.get(gc)
        if gc_len is not None and gc_len in length_to_fasta:
            chrom_map[gc] = length_to_fasta[gc_len]
            continue
        # Fallback: use as-is (may fail later, will produce a warning)
        chrom_map[gc] = gc

    # Report any unresolved
    unresolved = [gc for gc, fc in chrom_map.items() if fc not in fai_lengths]
    if unresolved:
        print(
            f'  WARNING: {len(unresolved)} GFF chroms could not be mapped '
            f'to FASTA seqs (first 3: {unresolved[:3]})',
            flush=True,
        )
    return chrom_map


def extract_cds_sequence(genome_path: Path, cds_intervals: list,
                          chrom_map: dict) -> str:
    """
    Extract and concatenate CDS exons using samtools faidx.
    Handles multi-exon genes, strand, and chromosome name remapping.
    Returns empty string on failure.
    """
    if not cds_intervals:
        return ''

    strand = cds_intervals[0][3]
    # Sort exons by position: ascending for +, descending for -
    exons = sorted(cds_intervals, key=lambda x: x[1], reverse=(strand == '-'))

    pieces = []
    for gff_chrom, start0, end1, _strand in exons:
        fasta_chrom = chrom_map.get(gff_chrom, gff_chrom)
        region = f'{fasta_chrom}:{start0 + 1}-{end1}'   # samtools: 1-based inclusive
        try:
            res = subprocess.run(
                [SAMTOOLS, 'faidx', str(genome_path), region],
                capture_output=True, text=True, check=True
            )
            seq = ''.join(res.stdout.split('\n')[1:]).replace(' ', '').upper()
            if not seq:
                print(
                    f'    WARNING: empty sequence for {region} '
                    f'(gff_chrom={gff_chrom!r})',
                    flush=True,
                )
                return ''
            pieces.append(seq)
        except subprocess.CalledProcessError as e:
            print(
                f'    WARNING samtools failed for {region}: '
                f'{e.stderr.strip()[:120]}',
                flush=True,
            )
            return ''

    cds_seq = ''.join(pieces)
    if strand == '-':
        cds_seq = reverse_complement(cds_seq)
    return cds_seq


def load_og_list(og_list_file: Path) -> list:
    with open(og_list_file) as fh:
        return [line.strip() for line in fh if line.strip()]


def load_orthogroups(orthogroups_tsv: Path, og_set: set, species_list: list) -> dict:
    """
    Returns: {og_id: {species: [protein_id, ...]}}
    Only for OGs in og_set and species in species_list.
    """
    result = defaultdict(lambda: defaultdict(list))
    col_map = {}   # column_index → species_name

    with open(orthogroups_tsv) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        for i, col in enumerate(header):
            sp = col.replace('.aa', '')
            if sp in species_list:
                col_map[i] = sp

        for row in reader:
            if not row:
                continue
            og_id = row[0]
            if og_id not in og_set:
                continue
            for col_idx, sp in col_map.items():
                if col_idx < len(row) and row[col_idx].strip():
                    proteins = [p.strip() for p in row[col_idx].split(',')
                                if p.strip()]
                    result[og_id][sp] = proteins

    return result


# ─── Main ─────────────────────────────────────────────────────────────────────
def main():
    print('=== CDS Extraction for Missing BUSCO Species ===', flush=True)

    og_list = load_og_list(OG_LIST_FILE)
    og_set  = set(og_list)
    print(f'OGs to process: {len(og_set)}', flush=True)

    print('Loading OrthoFinder Orthogroups (this may take a moment) …', flush=True)
    orthogroups = load_orthogroups(ORTHOGROUPS_TSV, og_set, MISSING_SPECIES)
    print(f'Loaded OrthoGroup data for {len(orthogroups)} OGs', flush=True)

    total_written = 0
    total_failed  = 0

    for sp in MISSING_SPECIES:
        cfg      = SPECIES_CONFIG[sp]
        gff_fmt  = cfg['gff_fmt']
        gff_ext  = cfg['gff_ext']

        gff_path    = GFF_DIR    / f'{sp}.{gff_ext}'
        genome_path = GENOME_DIR / f'{sp}.fasta'

        if not gff_path.exists():
            print(f'\n[{sp}] ERROR: GFF not found: {gff_path}', flush=True)
            continue
        if not genome_path.exists():
            print(f'\n[{sp}] ERROR: Genome not found: {genome_path}', flush=True)
            continue

        # Output directory (fake BUSCO structure)
        busco_out_dir = (
            BUSCO_RUNS_DIR / f'{sp}_{BUSCO_LINEAGE}'
            / f'run_{BUSCO_LINEAGE}'
            / 'busco_sequences'
            / 'single_copy_busco_sequences'
        )
        busco_out_dir.mkdir(parents=True, exist_ok=True)

        print(f'\n[{sp}] GFF format: {gff_fmt}', flush=True)
        print(f'  Building GFF CDS index …', flush=True)
        cds_index = build_gff_cds_index(gff_path, gff_fmt)
        print(f'  CDS index: {len(cds_index)} proteins indexed', flush=True)

        ensure_fai(genome_path)
        chrom_map = build_chrom_name_map(genome_path, gff_path)
        print(f'  Chrom map: {len(chrom_map)} entries', flush=True)

        sp_written  = 0
        sp_failed   = 0
        sp_no_og    = 0
        sp_no_cds   = 0

        for og_id in og_list:
            og_proteins = orthogroups.get(og_id, {}).get(sp, [])

            # Require exactly one copy (single-copy OG for this species)
            if len(og_proteins) != 1:
                sp_no_og += 1
                continue

            protein_id = extract_protein_id(og_proteins[0], gff_fmt)

            cds_intervals = cds_index.get(protein_id, [])
            if not cds_intervals:
                sp_no_cds += 1
                if sp_no_cds <= 5:   # only print the first few warnings
                    print(
                        f'  WARNING: no CDS for {sp} protein '
                        f'"{protein_id}" (OG {og_id})',
                        flush=True,
                    )
                continue

            cds_seq = extract_cds_sequence(genome_path, cds_intervals, chrom_map)
            if not cds_seq:
                sp_failed += 1
                continue

            out_file = busco_out_dir / f'{og_id}.fna'
            with open(out_file, 'w') as fout:
                fout.write(f'>{sp}\n{cds_seq}\n')
            sp_written += 1

        print(
            f'  Written: {sp_written}  |  No-CDS-in-GFF: {sp_no_cds}  |  '
            f'Samtools-fail: {sp_failed}  |  Multi/absent-copy: {sp_no_og}',
            flush=True,
        )
        total_written += sp_written
        total_failed  += sp_failed

    print(f'\n=== Total CDS files written: {total_written}', flush=True)
    print(f'=== Total samtools failures: {total_failed}', flush=True)


if __name__ == '__main__':
    main()
