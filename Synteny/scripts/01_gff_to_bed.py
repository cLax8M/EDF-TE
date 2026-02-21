#!/usr/bin/env python3
"""
01_gff_to_bed.py — Convert GFF/GFF3 annotation files to 4-column BED and
                   genome-size files required for the synteny analysis.

Gene IDs in the BED file match the identifiers used by OrthoFinder (i.e., the
first word of each protein FASTA header). This requires --faa_dir so that the
protein accession -> genomic position mapping can be built from GFF + FASTA.

Protein FASTA header formats handled:
  NCBI  : >KAI3655145.1 hypothetical protein LOCUS_TAG [species]
  JGI   : >jgi|Taxon|num|GENE_NAME
  Plain : >GENE_ID description

Output per species:
  {out_dir}/{species}.bed        — chrom  start  end  protein_id  (0-based, sorted)
  {out_dir}/{species}.genome     — chrom  size
  {out_dir}/{species}.id_map.tsv — protein_id  gene_lookup_key  matched  (diagnostics)

Usage:
  python 01_gff_to_bed.py --gff_dir GFF_DIR --faa_dir FAA_DIR --out_dir BED_DIR
                          [--species SP1 SP2 ...] [--min_scaffold_genes 3]

Author: EDF-TE project
Last updated: 2026-02-20
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff_dir", required=True,
                   help="Directory containing {CLADE}_{Species}.gff[3] files")
    p.add_argument("--faa_dir", required=True,
                   help="Directory containing {CLADE}_{Species}.aa.fasta protein files")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for .bed, .genome and .id_map.tsv files")
    p.add_argument("--species", nargs="*",
                   help="Species prefixes to process (e.g. HO_Naegrub). "
                        "Default: all GFF files found in gff_dir")
    p.add_argument("--min_scaffold_genes", type=int, default=3,
                   help="Min genes on a scaffold to be included (default: 3)")
    p.add_argument("--dry_run", action="store_true",
                   help="Print species list without processing")
    return p.parse_args()


# ---------------------------------------------------------------------------
# GFF parsing — build gene_info dict keyed by ALL possible IDs
# ---------------------------------------------------------------------------

def open_gff(gff_path):
    """Open GFF with UTF-8, falling back to latin-1."""
    for enc in ("utf-8", "latin-1"):
        try:
            with open(gff_path, encoding=enc) as fh:
                fh.read(512)
            return open(gff_path, encoding=enc, errors="replace")
        except UnicodeDecodeError:
            continue
    return open(gff_path, encoding="latin-1", errors="replace")


def parse_gff(gff_path, min_scaffold_genes=3):
    """
    Parse GFF/GFF2/GFF3 -> returns:
      gene_info      : dict  key -> (chrom, start0, end, rank)
                       Multiple keys per gene: locus_tag, Name, stripped ID,
                       proteinId, name (GFF2), and jgi.p|TAXON|PROTEINID form.
      chrom_gene_list: dict  chrom -> [(start0, end, keys_set)] sorted
      genome_sizes   : dict  chrom -> size

    Handles three annotation styles:
      1. NCBI GFF3  : 'gene' features with key=value; locus_tag attribute
      2. JGI GFF3   : 'gene' features with key=value; Name=jgi.p|...|NUM
      3. JGI GFF2   : NO 'gene' features; 'CDS'/'exon' with space-separated
                      attrs: name "GENENAME"; proteinId NUM
    """
    raw_genes = defaultdict(list)   # chrom -> [(start0, end, keys_set)]
    genome_sizes = {}

    # --- Detect style: does the file have any 'gene' type features? ---
    has_gene_features = False
    with open_gff(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 3 and cols[2] == "gene":
                has_gene_features = True
                break

    # --- GFF2 style: aggregate CDS intervals per (chrom, gene_name) ---
    gene_per_chrom_name = {}  # (chrom, name) -> [min_start0, max_end, keys]

    with open_gff(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            genome_sizes[parts[1]] = int(parts[3])
                        except ValueError:
                            pass
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            feature = cols[2]

            # ---- GFF3 path: proper 'gene' features ----
            if has_gene_features and feature == "gene":
                chrom = cols[0]
                try:
                    start0 = int(cols[3]) - 1
                    end    = int(cols[4])
                except ValueError:
                    continue

                attr_str = cols[8]
                keys = set()

                if "=" in attr_str:
                    attrs = {}
                    for field in attr_str.strip().split(";"):
                        field = field.strip()
                        if "=" in field:
                            k, v = field.split("=", 1)
                            attrs[k.strip()] = v.strip()

                    if "locus_tag" in attrs:
                        keys.add(attrs["locus_tag"])
                    if "Name" in attrs:
                        keys.add(attrs["Name"])
                    if "ID" in attrs:
                        raw_id = attrs["ID"]
                        keys.add(raw_id)
                        keys.add(re.sub(r"^gene[-:]", "", raw_id))
                    # JGI GFF3: store proteinId and jgi.p|TAXON|PROTEINID form
                    if "proteinId" in attrs:
                        pid = attrs["proteinId"]
                        keys.add(pid)
                        if "Name" in attrs:
                            nm_parts = attrs["Name"].split("|")
                            if len(nm_parts) >= 2:
                                taxon = nm_parts[1]
                                keys.add(f"jgi.p|{taxon}|{pid}")
                                keys.add(f"jgi|{taxon}|{pid}")

                if keys:
                    raw_genes[chrom].append((start0, end, keys))
                    if chrom not in genome_sizes or end > genome_sizes[chrom]:
                        genome_sizes[chrom] = end

            # ---- GFF2 path: CDS or exon features (no gene features) ----
            elif not has_gene_features and feature in ("CDS", "exon"):
                chrom = cols[0]
                try:
                    start0 = int(cols[3]) - 1
                    end    = int(cols[4])
                except ValueError:
                    continue

                attr_str = cols[8]
                gene_name = None
                protein_id_num = None

                if "=" in attr_str:
                    # GFF3-style attrs even without gene features
                    attrs = {}
                    for field in attr_str.strip().split(";"):
                        if "=" in field:
                            k, v = field.split("=", 1)
                            attrs[k.strip()] = v.strip()
                    gene_name = attrs.get("Name") or attrs.get("name")
                    protein_id_num = attrs.get("proteinId")
                else:
                    # Old GFF2 space-separated: name "YAL069W"; proteinId 1
                    for tok in re.split(r";\s*", attr_str):
                        tok = tok.strip()
                        m = re.match(r'^name\s+"([^"]+)"', tok)
                        if m:
                            gene_name = m.group(1)
                        m2 = re.match(r'^proteinId\s+(\S+)', tok)
                        if m2:
                            protein_id_num = m2.group(1)

                if gene_name is None:
                    continue

                cg_key = (chrom, gene_name)
                if cg_key not in gene_per_chrom_name:
                    gene_per_chrom_name[cg_key] = [start0, end, set()]
                entry = gene_per_chrom_name[cg_key]
                entry[0] = min(entry[0], start0)
                entry[1] = max(entry[1], end)
                entry[2].add(gene_name)
                if protein_id_num:
                    entry[2].add(protein_id_num)
                if chrom not in genome_sizes or end > genome_sizes[chrom]:
                    genome_sizes[chrom] = end

    # Promote GFF2 genes into raw_genes
    if not has_gene_features:
        for (chrom, _), (s, e, keys) in gene_per_chrom_name.items():
            raw_genes[chrom].append((s, e, keys))

    # Filter short scaffolds, build flat gene_info dict
    gene_info = {}
    chrom_gene_list = {}
    for chrom, entries in raw_genes.items():
        if len(entries) < min_scaffold_genes:
            continue
        entries_sorted = sorted(entries, key=lambda x: x[0])
        chrom_gene_list[chrom] = entries_sorted
        for rank, (start0, end, keys) in enumerate(entries_sorted):
            for key in keys:
                gene_info[key] = (chrom, start0, end, rank)

    # ---- Extra pass for NCBI GFF3: map CDS protein_id → gene locus ----
    # Many NCBI proteins (NP_/XP_) have no locus_tag in the FAA description
    # but the GFF CDS feature carries protein_id=XP_... and gene=LOC...
    # which lets us directly add the accession as a lookup key.
    if has_gene_features:
        with open_gff(gff_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[2] != "CDS":
                    continue
                attr_str = cols[8]
                if "=" not in attr_str:
                    continue
                attrs = {}
                for field in attr_str.strip().split(";"):
                    if "=" in field:
                        k, v = field.split("=", 1)
                        attrs[k.strip()] = v.strip()
                cds_pid = attrs.get("protein_id")  # e.g. XP_048489190.1
                if not cds_pid:
                    continue
                if cds_pid in gene_info:
                    continue  # already known
                # Find the gene position via any gene-level key
                gene_pos = None
                for key_attr in ("gene", "locus_tag", "Name"):
                    candidate = attrs.get(key_attr)
                    if candidate and (gene_info.get(candidate)
                                      or gene_info.get(candidate.upper())
                                      or gene_info.get(candidate.lower())):
                        gene_pos = (gene_info.get(candidate)
                                    or gene_info.get(candidate.upper())
                                    or gene_info.get(candidate.lower()))
                        break
                if gene_pos:
                    gene_info[cds_pid] = gene_pos

    return gene_info, chrom_gene_list, genome_sizes


# ---------------------------------------------------------------------------
# Protein FASTA parsing — build protein_id -> gene_lookup_key mapping
# ---------------------------------------------------------------------------

def parse_faa_header(header):
    """
    Return (protein_id, [lookup_keys]) from a FASTA header line.
    Returns a LIST of lookup keys to try (most specific first).

    Formats:
      NCBI : >KAI3655145.1 hypothetical protein LOCUS_TAG [species]
             protein_id="KAI3655145.1", lookup=["LOCUS_TAG"]

      JGI GFF3: >jgi|Taxon|PROTEINID|MODEL_NAME
             protein_id=full_tag
             lookup=["PROTEINID", "MODEL_NAME", "jgi|Taxon|PROTEINID",
                     "jgi.p|Taxon|PROTEINID"]

      JGI GFF2: >jgi|Taxon|NUM|GENE_NAME  (NUM is sequential, matches proteinId)
             lookup=["GENE_NAME", "NUM", "jgi.p|Taxon|NUM"]

      Plain: >GENE_ID ...
             protein_id=gene_lookup_key=["GENE_ID"]
    """
    line = header.lstrip(">").rstrip()
    parts = line.split(None, 1)
    protein_id  = parts[0]
    description = parts[1] if len(parts) > 1 else ""

    if protein_id.startswith("jgi|"):
        fields = protein_id.split("|")
        taxon      = fields[1] if len(fields) >= 2 else ""
        num        = fields[2] if len(fields) >= 3 else ""
        model_name = fields[3] if len(fields) >= 4 else ""
        # Strip .tN transcript suffix from model name
        model_name_clean = re.sub(r"\.t\d+.*$", "", model_name)
        lookup_keys = []
        # Add all candidate keys (deduplicated, preserving order)
        for k in [model_name_clean, model_name, num,
                  f"jgi.p|{taxon}|{num}", f"jgi|{taxon}|{num}"]:
            if k and k not in lookup_keys:
                lookup_keys.append(k)
        return protein_id, lookup_keys

    # NCBI or plain: try to extract locus_tag from description
    if "[" in description:
        before_bracket = description[: description.index("[")].strip()
        before_bracket = re.sub(r",\s*\S.*$", "", before_bracket).strip()
        words = before_bracket.split()
        if words:
            candidate = words[-1]
            if re.match(r"^[A-Za-z0-9_:.\-]+$", candidate) and len(candidate) > 2:
                # Locus tag first; protein_id itself as fallback (matched by CDS extra pass)
                return protein_id, [candidate, protein_id]

    # Gene.transcript format: >g211176.t1 (no description) — strip .tN suffix
    m = re.match(r'^(.+)\.t\d+$', protein_id)
    if m:
        gene_id = m.group(1)
        return protein_id, [gene_id, protein_id]

    return protein_id, [protein_id]  # fallback


def parse_faa(faa_path):
    """Read a protein FASTA -> list of (protein_id, [lookup_keys])."""
    proteins = []
    for enc in ("utf-8", "latin-1"):
        try:
            with open(faa_path, encoding=enc) as fh:
                for line in fh:
                    if line.startswith(">"):
                        proteins.append(parse_faa_header(line))
            break
        except UnicodeDecodeError:
            continue
    return proteins


# ---------------------------------------------------------------------------
# Build BED for one species
# ---------------------------------------------------------------------------

def build_bed(prefix, gff_path, faa_path, min_scaffold_genes):
    """
    Returns:
      bed_records  : list of (chrom, start0, end, protein_id)  sorted
      genome_sizes : dict chrom -> size
      id_map       : list of (protein_id, gene_lookup_key, matched:bool)
      stats        : summary dict
    """
    gene_info, chrom_gene_list, genome_sizes = parse_gff(gff_path, min_scaffold_genes)
    n_gff_genes = sum(len(v) for v in chrom_gene_list.values())

    if faa_path is None:
        # No FAA available: use GFF IDs directly (first key of each gene)
        records = []
        for chrom, entries in sorted(chrom_gene_list.items()):
            for start0, end, keys in entries:
                key = sorted(keys)[0]
                records.append((chrom, start0, end, key))
        return records, genome_sizes, [], {
            "n_gff_genes": n_gff_genes, "n_faa_proteins": 0,
            "n_matched": len(records), "n_unmatched": 0,
        }

    proteins = parse_faa(faa_path)
    n_faa = len(proteins)
    id_map = []
    bed_positions = {}  # protein_id -> (chrom, start0, end)
    unmatched = 0

    for protein_id, lookup_keys in proteins:
        pos = None
        matched_key = None
        for key in lookup_keys:
            pos = (gene_info.get(key)
                   or gene_info.get(key.upper())
                   or gene_info.get(key.lower()))
            if pos is not None:
                matched_key = key
                break
        if pos is not None:
            bed_positions[protein_id] = (pos[0], pos[1], pos[2])
            id_map.append((protein_id, matched_key, True))
        else:
            id_map.append((protein_id, lookup_keys[0] if lookup_keys else "", False))
            unmatched += 1

    # Group by chrom, apply min_scaffold_genes filter on matched proteins
    bed_by_chrom = defaultdict(list)
    for pid, (chrom, start0, end) in bed_positions.items():
        bed_by_chrom[chrom].append((start0, end, pid))

    records = []
    for chrom in sorted(bed_by_chrom):
        entries = sorted(bed_by_chrom[chrom])
        if len(entries) < min_scaffold_genes:
            continue
        for start0, end, pid in entries:
            records.append((chrom, start0, end, pid))

    return records, genome_sizes, id_map, {
        "n_gff_genes": n_gff_genes, "n_faa_proteins": n_faa,
        "n_matched": len(bed_positions), "n_unmatched": unmatched,
    }


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def find_files(gff_dir, faa_dir, species_list=None):
    """Return list of (species_prefix, gff_path, faa_path) for all species."""
    gff_dir = Path(gff_dir)
    faa_dir = Path(faa_dir)
    results = []
    seen = set()

    for ext in ("*.gff", "*.gff3"):
        for gff_fp in sorted(gff_dir.glob(ext)):
            if gff_fp.name.startswith("._"):  # skip macOS hidden files
                continue
            prefix = re.sub(r"_c$", "", gff_fp.stem)
            if prefix in seen:
                continue
            if species_list and prefix not in species_list:
                continue
            seen.add(prefix)

            # Locate FAA (try both plain and _c variants)
            faa_fp = None
            for candidate in (f"{prefix}.aa.fasta", f"{prefix}_c.aa.fasta",
                               f"{prefix}.aa.fa",    f"{prefix}_c.aa.fa"):
                fp = faa_dir / candidate
                if fp.exists():
                    faa_fp = fp
                    break

            results.append((prefix, gff_fp, faa_fp))

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    species_files = find_files(args.gff_dir, args.faa_dir, args.species)

    if not species_files:
        print(f"ERROR: no GFF files found in {args.gff_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(species_files)} species to process.")

    if args.dry_run:
        for prefix, gff_fp, faa_fp in species_files:
            faa_status = faa_fp.name if faa_fp else "NO FAA FOUND"
            print(f"  {prefix}  GFF={gff_fp.name}  FAA={faa_status}")
        return

    summary_rows = []

    for prefix, gff_path, faa_path in species_files:
        print(f"  {prefix} ...", end=" ", flush=True)

        records, genome_sizes, id_map, stats = build_bed(
            prefix, gff_path, faa_path, args.min_scaffold_genes
        )

        n_faa = stats["n_faa_proteins"]
        n_matched = stats["n_matched"]
        pct = 100 * n_matched / n_faa if n_faa > 0 else 0

        if not records:
            print(f"WARNING: 0 BED entries "
                  f"(GFF={stats['n_gff_genes']}, FAA={n_faa}, "
                  f"matched={n_matched} [{pct:.0f}%]) — skipped")
            continue

        # Write BED
        bed_out = out_dir / f"{prefix}.bed"
        with open(bed_out, "w") as fh:
            for chrom, start0, end, pid in records:
                fh.write(f"{chrom}\t{start0}\t{end}\t{pid}\n")

        # Write genome sizes (all scaffolds, for bedtools shuffle)
        genome_out = out_dir / f"{prefix}.genome"
        with open(genome_out, "w") as fh:
            for chrom, size in sorted(genome_sizes.items()):
                fh.write(f"{chrom}\t{size}\n")

        # Write ID map (diagnostics)
        if id_map:
            idmap_out = out_dir / f"{prefix}.id_map.tsv"
            with open(idmap_out, "w") as fh:
                fh.write("protein_id\tgene_lookup_key\tmatched\n")
                for pid, key, found in id_map:
                    fh.write(f"{pid}\t{key}\t{'yes' if found else 'no'}\n")

        n_scaffolds = len(set(r[0] for r in records))
        print(f"{len(records)} genes | {n_scaffolds} scaffolds | "
              f"{n_matched}/{n_faa} FAA matched ({pct:.0f}%)")

        summary_rows.append({
            "species": prefix, "n_bed_genes": len(records),
            "n_scaffolds": n_scaffolds, "n_gff_genes": stats["n_gff_genes"],
            "n_faa_proteins": n_faa, "n_matched": n_matched,
            "pct_matched": f"{pct:.1f}",
        })

    summary_path = out_dir / "01_gff_to_bed_summary.tsv"
    with open(summary_path, "w") as fh:
        fh.write("species\tn_bed_genes\tn_scaffolds\tn_gff_genes\t"
                 "n_faa_proteins\tn_matched\tpct_matched\n")
        for row in summary_rows:
            fh.write(
                f"{row['species']}\t{row['n_bed_genes']}\t{row['n_scaffolds']}\t"
                f"{row['n_gff_genes']}\t{row['n_faa_proteins']}\t"
                f"{row['n_matched']}\t{row['pct_matched']}\n"
            )
    print(f"\nSummary written to {summary_path}")
    print("Step 01 complete.")


if __name__ == "__main__":
    main()
