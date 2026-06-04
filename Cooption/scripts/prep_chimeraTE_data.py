#!/usr/bin/env python3
"""
prep_chimeraTE_data.py
Prepare data for ChimeraTE Mode 1 for MUC species:
  - De-interleave paired-end reads (interleaved -> R1 + R2)
  - Convert JGI GFF2 -> GTF (gene_id / transcript_id attributes)
  - Convert EarlGrey filteredRepeats.bed -> TE GTF (filter Simple_repeat/Satellite/etc.)
  - Write ChimeraTE input table (tab-delimited: R1, R2, replicate_name)

Usage:
    python3 prep_chimeraTE_data.py --species MU_Mucci3 --reads_dir /path/to/reads
"""

import os
import sys
import gzip
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict

# ─── Paths ──────────────────────────────────────────────────────────────────
BASE_DATA   = "/legserv/Temp/Carlos/EDF-TE/data"
BASE_RES    = "/legserv/Temp/Carlos/EDF-TE/results"
CHIMERA_OUT = f"{BASE_RES}/Co-option/chimerateTE"
READS_BASE_MUC  = "/legserv/Temp/Carlos/EDF-TE/data/Extra/RNAseqRawMuc"
READS_BASE_EDF  = "/legserv/Temp/Carlos/EDF-TE/data/rnaseq"
TE_BED_BASE_EDF = "/legserv/Temp/Carlos/EDF-TE/results/HTT/04_te_annotation"
TE_BED_BASE_MUC = "/data/carlos/EDF-TE/TE_analysis/annotated"

# Map species code -> (reads_dir, reads_type)
# reads_type: 'interleaved' = one file per rep to de-interleave
#             'paired'      = *_1.fastq.gz + *_2.fastq.gz already split
SPECIES_CONFIG = {
    # ── Original 3 Mucoromycotina ─────────────────────────────────────────
    "MU_Mucci3":  {"reads_dir": f"{READS_BASE_MUC}/Mucor",      "reads_type": "interleaved"},
    "MU_Rhimi59": {"reads_dir": f"{READS_BASE_MUC}/Rhizopus",   "reads_type": "interleaved"},
    "MU_Phybl4":  {"reads_dir": f"{READS_BASE_MUC}/Phycomyces", "reads_type": "interleaved"},
    # ── EDF additional species (paired-end, split reads from SRA) ─────────
    "AP_Amoprot": {"reads_dir": f"{READS_BASE_EDF}/AP_Amoprot", "reads_type": "paired"},
    "BL_Blabri1": {"reads_dir": f"{READS_BASE_EDF}/BL_Blabri1", "reads_type": "paired"},
    "BS_Basme92": {"reads_dir": f"{READS_BASE_EDF}/BS_Basme92", "reads_type": "paired"},
    "CH_Batde5":  {"reads_dir": f"{READS_BASE_EDF}/CH_Batde5",  "reads_type": "paired"},
    "CH_Batsa1":  {"reads_dir": f"{READS_BASE_EDF}/CH_Batsa1",  "reads_type": "paired"},
    "CH_Pirfi3":  {"reads_dir": f"{READS_BASE_EDF}/CH_Pirfi3",  "reads_type": "paired"},
    "EN_Entmai1": {"reads_dir": f"{READS_BASE_EDF}/EN_Entmai1", "reads_type": "paired"},
    "EN_Zoorad1": {"reads_dir": f"{READS_BASE_EDF}/EN_Zoorad1", "reads_type": "paired"},
    "GL_Geopyr1": {"reads_dir": f"{READS_BASE_EDF}/GL_Geopyr1", "reads_type": "paired"},
    "GL_Rhicl1":  {"reads_dir": f"{READS_BASE_EDF}/GL_Rhicl1",  "reads_type": "paired"},
    "GL_RhiirA4": {"reads_dir": f"{READS_BASE_EDF}/GL_RhiirA4", "reads_type": "paired"},
    "KX_Coemoj1": {"reads_dir": f"{READS_BASE_EDF}/KX_Coemoj1", "reads_type": "paired"},
    "KX_Kicala1": {"reads_dir": f"{READS_BASE_EDF}/KX_Kicala1", "reads_type": "paired"},
    "MI_Enccu1":  {"reads_dir": f"{READS_BASE_EDF}/MI_Enccu1",  "reads_type": "paired"},
    "MI_Sprlo1":  {"reads_dir": f"{READS_BASE_EDF}/MI_Sprlo1",  "reads_type": "paired"},
    "MI_Vavcu1":  {"reads_dir": f"{READS_BASE_EDF}/MI_Vavcu1",  "reads_type": "paired"},
    "MO_Disdec1": {"reads_dir": f"{READS_BASE_EDF}/MO_Disdec1", "reads_type": "paired"},
    "MO_Morel2":  {"reads_dir": f"{READS_BASE_EDF}/MO_Morel2",  "reads_type": "paired"},
    "NU_Fontalba": {"reads_dir": f"{READS_BASE_EDF}/NU_Fontalba","reads_type": "paired"},
    "ZO_Pipcy3":  {"reads_dir": f"{READS_BASE_EDF}/ZO_Pipcy3",  "reads_type": "paired"},
}
# Legacy alias for backward compatibility
SPECIES_READS_MAP = {k: v["reads_dir"] for k, v in SPECIES_CONFIG.items()}

# ─── TE classes to EXCLUDE from ChimeraTE (same as RepeatMasker convention) ──
TE_EXCLUDE_CLASSES = {
    "Simple_repeat", "Satellite", "Low_complexity", "rRNA",
    "snRNA", "srpRNA", "tRNA", "scRNA", "RNA", "ARTEFACT",
}

# ─── Genome SA index sizes (min(14, floor(log2(genomeSize)/2 - 1))) ──────────
GENOME_SA_SIZES = {
    # Original 3
    "MU_Mucci3":  11,
    "MU_Rhimi59": 11,
    "MU_Phybl4":  12,
    # EDF additional
    "AP_Amoprot": 11,
    "BL_Blabri1": 11,
    "BS_Basme92": 12,
    "CH_Batde5":  11,
    "CH_Batsa1":  11,
    "CH_Pirfi3":  11,
    "EN_Entmai1": 14,
    "EN_Zoorad1": 13,
    "GL_Geopyr1": 12,
    "GL_Rhicl1":  12,
    "GL_RhiirA4": 12,
    "KX_Coemoj1": 10,
    "KX_Kicala1": 11,
    "MI_Enccu1":   9,
    "MI_Sprlo1":  10,
    "MI_Vavcu1":  10,
    "MO_Disdec1": 11,
    "MO_Morel2":  11,
    "NU_Fontalba": 11,
    "ZO_Pipcy3":  10,
}


def deinterleave_fastq(input_gz: Path, out_r1: Path, out_r2: Path):
    """Split interleaved paired-end FASTQ (both mates in one file) into R1 and R2."""
    print(f"  De-interleaving {input_gz.name} ...", flush=True)
    n_pairs = 0
    with gzip.open(input_gz, "rt") as fh, \
         gzip.open(out_r1, "wt") as r1, \
         gzip.open(out_r2, "wt") as r2:
        while True:
            # Read 4 lines = 1 FASTQ record for R1
            h1 = fh.readline()
            if not h1:
                break
            s1 = fh.readline()
            p1 = fh.readline()
            q1 = fh.readline()
            # Read 4 lines = 1 FASTQ record for R2
            h2 = fh.readline()
            s2 = fh.readline()
            p2 = fh.readline()
            q2 = fh.readline()
            if not h2:
                print(f"  WARNING: odd number of records in {input_gz.name}", flush=True)
                r1.write(h1); r1.write(s1); r1.write(p1); r1.write(q1)
                break
            r1.write(h1); r1.write(s1); r1.write(p1); r1.write(q1)
            r2.write(h2); r2.write(s2); r2.write(p2); r2.write(q2)
            n_pairs += 1
    print(f"    Wrote {n_pairs:,} pairs -> {out_r1.name} / {out_r2.name}", flush=True)
    return n_pairs


def gff2_to_gtf(gff_path: Path, gtf_path: Path):
    """Convert JGI GFF2 format to GTF with gene_id and transcript_id attributes."""
    print(f"  Converting GFF2 -> GTF: {gff_path.name} ...", flush=True)

    # First pass: collect exon coordinates per transcript to build gene/transcript bounds
    transcripts = defaultdict(lambda: {"chr": None, "source": None, "strand": None,
                                        "start": float("inf"), "end": 0, "exons": []})

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chr_, source, feat, start, end, score, strand, frame, attrs = parts
            if feat not in ("exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR",
                            "stop_codon", "start_codon"):
                continue
            # Extract 'name' attribute as transcript ID
            tx_id = None
            for token in attrs.replace('"', '').split(";"):
                token = token.strip()
                if token.lower().startswith("name "):
                    tx_id = token.split(" ", 1)[1].strip()
                    break
            if tx_id is None:
                continue
            start_i = int(start)
            end_i   = int(end)
            t = transcripts[tx_id]
            t["chr"]    = chr_
            t["source"] = source
            t["strand"] = strand
            t["start"]  = min(t["start"], start_i)
            t["end"]    = max(t["end"],   end_i)
            if feat == "exon":
                t["exons"].append((start_i, end_i, score, strand, frame))

    # Second pass: write GTF
    n_genes = 0
    with open(gff_path) as fh_in, open(gtf_path, "w") as fh_out:
        written_tx = set()

        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                chr_, source, feat, start, end, score, strand, frame, attrs = parts
                if feat not in ("exon", "CDS", "UTR", "five_prime_UTR",
                                "three_prime_UTR", "stop_codon", "start_codon"):
                    continue
                # Extract transcript id
                tx_id = None
                for token in attrs.replace('"', '').split(";"):
                    token = token.strip()
                    if token.lower().startswith("name "):
                        tx_id = token.split(" ", 1)[1].strip()
                        break
                if tx_id is None:
                    continue

                # Emit gene + transcript features on first encounter
                if tx_id not in written_tx:
                    written_tx.add(tx_id)
                    t = transcripts[tx_id]
                    attr_str = f'gene_id "{tx_id}"; transcript_id "{tx_id}";'
                    fh_out.write(f'{t["chr"]}\t{t["source"]}\tgene\t{t["start"]}\t{t["end"]}\t.\t{t["strand"]}\t.\t{attr_str}\n')
                    fh_out.write(f'{t["chr"]}\t{t["source"]}\ttranscript\t{t["start"]}\t{t["end"]}\t.\t{t["strand"]}\t.\t{attr_str}\n')
                    n_genes += 1

                attr_str = f'gene_id "{tx_id}"; transcript_id "{tx_id}";'
                fh_out.write(f'{chr_}\t{source}\t{feat}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attr_str}\n')

    print(f"    Wrote {n_genes:,} transcripts -> {gtf_path.name}", flush=True)


def gff3_to_gtf(gff_path: Path, gtf_path: Path):
    """Convert JGI GFF3 format (ID=mRNA_N;Parent=gene_N; exon Parent=mRNA_N) to GTF."""
    print(f"  Converting GFF3 -> GTF: {gff_path.name} ...", flush=True)

    def parse_attrs(attrs_str):
        d = {}
        for token in attrs_str.split(";"):
            token = token.strip()
            if "=" in token:
                k, v = token.split("=", 1)
                d[k.strip()] = v.strip()
        return d

    # First pass: build mRNA_id -> gene bounds & strand from mRNA records
    mrna_info = {}  # mRNA_id -> {chr, source, strand, start, end, gene_id}
    exon_records = []  # (chr, source, start, end, score, strand, frame, mRNA_id)

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chr_, source, feat, start, end, score, strand, frame, attrs = parts
            d = parse_attrs(attrs)
            if feat == "mRNA":
                mrna_id = d.get("ID")
                gene_id = d.get("Parent", mrna_id)
                if mrna_id:
                    mrna_info[mrna_id] = {
                        "chr": chr_, "source": source, "strand": strand,
                        "start": int(start), "end": int(end), "gene_id": gene_id,
                    }
            elif feat == "exon":
                parent = d.get("Parent")
                if parent:
                    exon_records.append((chr_, source, int(start), int(end),
                                         score, strand, frame, parent))

    # Second pass: write GTF
    n_genes = 0
    written_mrna = set()
    with open(gtf_path, "w") as fh_out:
        for (chr_, source, start, end, score, strand, frame, mrna_id) in exon_records:
            info = mrna_info.get(mrna_id)
            if info is None:
                continue
            gene_id = info["gene_id"]
            attr_str = f'gene_id "{gene_id}"; transcript_id "{mrna_id}";'
            if mrna_id not in written_mrna:
                written_mrna.add(mrna_id)
                fh_out.write(f'{info["chr"]}\t{info["source"]}\tgene\t{info["start"]}\t{info["end"]}'
                             f'\t.\t{info["strand"]}\t.\t{attr_str}\n')
                fh_out.write(f'{info["chr"]}\t{info["source"]}\ttranscript\t{info["start"]}\t{info["end"]}'
                             f'\t.\t{info["strand"]}\t.\t{attr_str}\n')
                n_genes += 1
            fh_out.write(f'{chr_}\t{source}\texon\t{start}\t{end}'
                         f'\t{score}\t{strand}\t{frame}\t{attr_str}\n')

    print(f"    Wrote {n_genes:,} transcripts -> {gtf_path.name}", flush=True)


def bed_to_te_gtf(bed_path: Path, gtf_path: Path):
    """Convert EarlGrey filteredRepeats BED (0-based) -> TE GTF (1-based)
    filtering out non-TE classes (Simple_repeat, Satellite, etc.)"""
    print(f"  Converting TE BED -> GTF: {bed_path.name} ...", flush=True)
    n_written = 0
    n_skipped = 0
    with open(bed_path) as fh_in, open(gtf_path, "w") as fh_out:
        for line in fh_in:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chr_   = parts[0]
            start  = int(parts[1]) + 1   # BED 0-based -> GTF 1-based
            end    = int(parts[2])
            te_name = parts[3]
            score  = parts[4] if len(parts) > 4 else "."
            strand = parts[5] if len(parts) > 5 else "."

            # Filter unwanted classes
            te_class = te_name.split("/")[0] if "/" in te_name else te_name
            if te_class in TE_EXCLUDE_CLASSES:
                n_skipped += 1
                continue

            fh_out.write(f'{chr_}\tRepeatMasker\tsimilarity\t{start}\t{end}\t{score}\t{strand}\t.\t{te_name}\n')
            n_written += 1

    print(f"    Wrote {n_written:,} TE entries (skipped {n_skipped:,} non-TE)", flush=True)


def prepare_species(species: str, force: bool = False):
    cfg = SPECIES_CONFIG.get(species)
    if cfg is None:
        print(f"ERROR: unknown species '{species}'. Add to SPECIES_CONFIG.", flush=True)
        sys.exit(1)

    reads_dir    = Path(cfg["reads_dir"])
    reads_type   = cfg["reads_type"]   # 'interleaved' or 'paired'
    genome_path  = Path(BASE_DATA) / "genomes" / f"{species}.fasta"
    gff_path     = Path(BASE_DATA) / "genes_annot_gff" / f"{species}.gff"
    if not gff_path.exists():
        gff_path = Path(BASE_DATA) / "genes_annot_gff" / f"{species}.gff3"
    # TE BED: try EDF path first, then EarlGrey summaryFiles
    te_bed_path = Path(TE_BED_BASE_EDF) / f"{species}.filteredRepeats.bed"
    if not te_bed_path.exists():
        te_bed_path = Path(TE_BED_BASE_MUC) / f"{species}_EarlGrey" / f"{species}_summaryFiles" / f"{species}.filteredRepeats.bed"

    out_dir  = Path(CHIMERA_OUT) / species
    prep_dir = out_dir / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}", flush=True)
    print(f"Preparing {species}", flush=True)
    print(f"{'='*60}", flush=True)

    # Sanity checks
    for p, label in [(genome_path, "Genome"), (gff_path, "GFF"), (te_bed_path, "TE BED"), (reads_dir, "Reads dir")]:
        if not p.exists():
            print(f"ERROR: {label} not found: {p}", flush=True)
            sys.exit(1)
        print(f"  {label}: {p} OK", flush=True)

    # ── 1. Convert GFF -> GTF ─────────────────────────────────────────────
    gtf_genes = prep_dir / f"{species}_genes.gtf"
    if force or not gtf_genes.exists() or gtf_genes.stat().st_size == 0:
        if gff_path.suffix == ".gff3":
            gff3_to_gtf(gff_path, gtf_genes)
        else:
            gff2_to_gtf(gff_path, gtf_genes)
    else:
        print(f"  GTF (genes) already exists, skipping: {gtf_genes.name}", flush=True)

    # ── 2. Convert TE BED -> GTF ──────────────────────────────────────────
    gtf_te = prep_dir / f"{species}_TEs.gtf"
    if force or not gtf_te.exists():
        bed_to_te_gtf(te_bed_path, gtf_te)
    else:
        print(f"  GTF (TEs) already exists, skipping: {gtf_te.name}", flush=True)

    # ── 3. De-interleave reads ────────────────────────────────────────────
    reads_out_dir = prep_dir / "reads"
    reads_out_dir.mkdir(exist_ok=True)

    if reads_type == "interleaved":
        # Original Mucoromycotina: one interleaved file per rep
        input_reads = sorted(reads_dir.glob("*.filter-RNA.fastq.gz"))
        if not input_reads:
            input_reads = sorted(reads_dir.glob("*.fastq.gz"))
            raw_names = {r.name.replace(".filter-RNA.fastq.gz", ".fastq.gz") for r in input_reads}
            input_reads = [r for r in input_reads if r.name not in raw_names]

        print(f"  Found {len(input_reads)} interleaved read file(s)", flush=True)
        input_table_rows = []
        for i, gz_file in enumerate(input_reads, 1):
            rep_name = f"rep{i:02d}"
            out_r1   = reads_out_dir / f"{rep_name}_R1.fastq.gz"
            out_r2   = reads_out_dir / f"{rep_name}_R2.fastq.gz"
            if force or not out_r1.exists():
                deinterleave_fastq(gz_file, out_r1, out_r2)
            else:
                print(f"  {rep_name} already de-interleaved, skipping", flush=True)
            input_table_rows.append(f"{out_r1.resolve()}\t{out_r2.resolve()}\t{rep_name}")
    else:
        # EDF species: *_1.fastq.gz + *_2.fastq.gz already split (single replicate)
        r1_files = sorted(reads_dir.glob("*_1.fastq.gz"))
        r2_files = sorted(reads_dir.glob("*_2.fastq.gz"))
        # Filter out wrong SRR6057017 reads if present
        r1_files = [f for f in r1_files if "SRR6057017" not in f.name]
        r2_files = [f for f in r2_files if "SRR6057017" not in f.name]
        if not r1_files:
            print(f"  [WARN] No paired reads found in {reads_dir} (skipping prep)", flush=True)
            return None
        print(f"  Found {len(r1_files)} paired-end read set(s)", flush=True)
        input_table_rows = []
        for i, (r1, r2) in enumerate(zip(sorted(r1_files), sorted(r2_files)), 1):
            rep_name = f"rep{i:02d}"
            # Symlink or copy into reads subdir with standard names
            out_r1 = reads_out_dir / f"{rep_name}_R1.fastq.gz"
            out_r2 = reads_out_dir / f"{rep_name}_R2.fastq.gz"
            if force or not out_r1.exists():
                if not out_r1.exists():
                    out_r1.symlink_to(r1.resolve())
                if not out_r2.exists():
                    out_r2.symlink_to(r2.resolve())
            print(f"  {rep_name}: {r1.name} + {r2.name}", flush=True)
            input_table_rows.append(f"{out_r1.resolve()}\t{out_r2.resolve()}\t{rep_name}")

    # ── 4. Write ChimeraTE input table ────────────────────────────────────
    input_table = prep_dir / f"{species}_input_table.tsv"
    with open(input_table, "w") as fh:
        fh.write("\n".join(input_table_rows) + "\n")
    print(f"  Input table: {input_table} ({len(input_table_rows)} replicates)", flush=True)

    print(f"\nPreparation complete for {species}:", flush=True)
    print(f"  Gene GTF : {gtf_genes}", flush=True)
    print(f"  TE GTF   : {gtf_te}", flush=True)
    print(f"  Input tbl: {input_table}", flush=True)

    return {
        "species":      species,
        "genome":       str(genome_path),
        "gtf_genes":    str(gtf_genes),
        "gtf_te":       str(gtf_te),
        "input_table":  str(input_table),
        "out_dir":      str(out_dir),
        "sa_index":     GENOME_SA_SIZES.get(species, 11),
    }


def main():
    parser = argparse.ArgumentParser(description="Prepare data for ChimeraTE Mode 1")
    parser.add_argument("--species", nargs="+",
                        default=list(SPECIES_CONFIG.keys()),
                        help="Species code(s) to process (default: all)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if output files already exist")
    args = parser.parse_args()

    for sp in args.species:
        result = prepare_species(sp, force=args.force)
        if result is None:
            print(f"  Skipped {sp} (no reads yet?)", flush=True)

    print("\nAll species prepared.", flush=True)


if __name__ == "__main__":
    main()
