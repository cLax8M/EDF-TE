#!/usr/bin/env python3
"""
07_circos_prep.py — Prepare Circos v0.67 input files for selected species pairs.

For a chosen species pair (SpA, SpB), generates:
  karyotype.{SpA}.txt     — scaffold lengths for SpA (Circos karyotype format)
  karyotype.{SpB}.txt     — scaffold lengths for SpB
  links.{SpA}__{SpB}.txt  — syntenic gene-pair links (collinear pairs from step 03)
  circos.conf             — minimal Circos configuration file

Only scaffolds that carry at least one collinear pair are included in the
karyotype, to keep the plot readable. If a species has many scaffolds, only
the top N by collinear-pair count are retained (--max_scaffolds, default 30).

The final Circos plot can be rendered with:
  circos -conf circos.conf

Usage:
  python 07_circos_prep.py --sp_a SpA --sp_b SpB  --clade CLADE
                           --bed_dir BED_DIR
                           --collinear_dir COLLINEAR_DIR
                           --out_dir CIRCOS_DIR
                           [--max_scaffolds 30]

Author: EDF-TE project
Last updated: 2026-02-20
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


# Circos color palette (one color per SpA scaffold)
COLORS = [
    "red", "orange", "yellow", "green", "blue", "purple", "grey",
    "vdred", "vdorange", "vdyellow", "vdgreen", "vdblue",
    "lblue", "lgreen", "lpurple", "lorange",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sp_a",          required=True, help="Species A prefix (e.g. HO_Lucicu)")
    p.add_argument("--sp_b",          required=True, help="Species B prefix (e.g. HO_Luciser)")
    p.add_argument("--clade",         required=True, help="Clade code (e.g. HO)")
    p.add_argument("--bed_dir",       required=True, help="BED directory (step 01 output)")
    p.add_argument("--collinear_dir", required=True, help="Collinear directory (step 03 output)")
    p.add_argument("--out_dir",       required=True, help="Output directory for Circos files")
    p.add_argument("--max_scaffolds", type=int, default=30,
                   help="Max scaffolds per species to include (default: 30)")
    return p.parse_args()


def load_genome_sizes(bed_dir, sp):
    """Read {sp}.genome -> dict chrom -> size."""
    genome_file = Path(bed_dir) / f"{sp}.genome"
    sizes = {}
    with open(genome_file) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                try:
                    sizes[cols[0]] = int(cols[1])
                except ValueError:
                    pass
    return sizes


def load_bed_positions(bed_dir, sp):
    """Read {sp}.bed -> dict gene_id -> (chrom, start, end)."""
    gene_pos = {}
    bed_file = Path(bed_dir) / f"{sp}.bed"
    with open(bed_file) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 4:
                gene_pos[cols[3]] = (cols[0], int(cols[1]), int(cols[2]))
    return gene_pos


def load_collinear(collinear_dir, clade, sp_a, sp_b):
    """
    Load collinear pairs. Returns list of (chrom_a,ga1,ga2,chrom_b,gb1,gb2).
    Tries both orderings.
    """
    col_dir = Path(collinear_dir) / clade
    f1 = col_dir / f"{sp_a}__{sp_b}.collinear.tsv"
    f2 = col_dir / f"{sp_b}__{sp_a}.collinear.tsv"

    pairs = []
    target = f1 if f1.exists() else (f2 if f2.exists() else None)
    if target is None:
        print(f"ERROR: collinear file not found for {sp_a} × {sp_b}", file=sys.stderr)
        sys.exit(1)

    reversed_order = (target == f2)
    with open(target) as fh:
        next(fh)  # header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            if reversed_order:
                pairs.append((cols[3], cols[4], cols[5], cols[0], cols[1], cols[2]))
            else:
                pairs.append(tuple(cols[:6]))
    return pairs


def select_top_scaffolds(pairs, pos_a, pos_b, max_scaffolds):
    """
    Count collinear pairs per scaffold for each species and select top N,
    ensuring both scaffolds of a pair are retained together.
    """
    scaf_count_a = defaultdict(int)
    scaf_count_b = defaultdict(int)

    for ca, ga1, ga2, cb, gb1, gb2 in pairs:
        scaf_count_a[ca] += 1
        scaf_count_b[cb] += 1

    top_a = set(sorted(scaf_count_a, key=lambda x: -scaf_count_a[x])[:max_scaffolds])
    top_b = set(sorted(scaf_count_b, key=lambda x: -scaf_count_b[x])[:max_scaffolds])

    # Filter pairs to only those where both scaffolds are in top sets
    filtered = [(ca, ga1, ga2, cb, gb1, gb2) for ca, ga1, ga2, cb, gb1, gb2 in pairs
                if ca in top_a and cb in top_b]

    # Final scaffold sets based on filtered pairs
    retained_a = set(p[0] for p in filtered)
    retained_b = set(p[3] for p in filtered)

    return filtered, retained_a, retained_b


def write_karyotype(sp, scaffolds, genome_sizes, out_dir, prefix, color_map):
    """Write Circos karyotype file."""
    kary_file = Path(out_dir) / f"karyotype.{sp}.txt"
    with open(kary_file, "w") as fh:
        for scaf in sorted(scaffolds, key=lambda x: -genome_sizes.get(x, 0)):
            size = genome_sizes.get(scaf, 1)
            color = color_map.get(scaf, "grey")
            # Circos karyotype: chr - ID LABEL start end color
            chr_id = f"{prefix}_{scaf.replace('|', '_').replace('.', '_')}"
            fh.write(f"chr - {chr_id} {scaf} 0 {size} {color}\n")
    return kary_file


def write_links(pairs, pos_a, pos_b, sp_a, sp_b, out_dir, scaf_id_a, scaf_id_b):
    """Write Circos links file (one link per collinear gene pair)."""
    links_file = Path(out_dir) / f"links.{sp_a}__{sp_b}.txt"
    with open(links_file, "w") as fh:
        for ca, ga1, ga2, cb, gb1, gb2 in pairs:
            # Write midpoint of each gene pair as the link anchor
            def midpoint(gene, pos_dict, scaf_id):
                if gene in pos_dict:
                    chrom, s, e = pos_dict[gene]
                    cid = scaf_id.get(chrom, chrom.replace("|","_").replace(".","_"))
                    mid = (s + e) // 2
                    return f"{cid} {mid} {mid}"
                return None

            cid_a = f"A_{ca.replace('|','_').replace('.','_')}"
            cid_b = f"B_{cb.replace('|','_').replace('.','_')}"

            pa1 = pos_a.get(ga1)
            pa2 = pos_a.get(ga2)
            pb1 = pos_b.get(gb1)
            pb2 = pos_b.get(gb2)

            if not all([pa1, pa2, pb1, pb2]):
                continue

            # Link A-region to B-region (use start of first gene to end of second)
            start_a = min(pa1[1], pa2[2])
            end_a   = max(pa1[2], pa2[2])
            start_b = min(pb1[1], pb2[2])
            end_b   = max(pb1[2], pb2[2])

            fh.write(f"{cid_a} {start_a} {end_a} {cid_b} {start_b} {end_b}\n")

    return links_file


def write_circos_conf(sp_a, sp_b, kary_a, kary_b, links_file, out_dir):
    """Write a minimal Circos configuration file."""
    conf_file = Path(out_dir) / "circos.conf"
    with open(conf_file, "w") as fh:
        fh.write(f"""\
# Circos configuration for {sp_a} vs {sp_b}
# Generated by 07_circos_prep.py
# Run with: circos -conf circos.conf

karyotype = {kary_a},{kary_b}

<ideogram>
  <spacing>
    default = 0.005r
  </spacing>
  radius           = 0.85r
  thickness        = 20p
  fill             = yes
  stroke_color     = dgrey
  stroke_thickness = 2p
  show_label       = yes
  label_font       = default
  label_radius     = dims(ideogram,radius) + 0.05r
  label_size       = 18
  label_parallel   = yes
</ideogram>

<ticks>
  radius           = 1r
  color            = black
  thickness        = 2p
  multiplier       = 1e-6
  format           = %d Mb
  <tick>
    spacing        = 10Mb
    size           = 10p
    show_label     = yes
    label_size     = 12p
  </tick>
</ticks>

<links>
  <link>
    file          = {links_file}
    radius        = 0.95r
    bezier_radius = 0.1r
    color         = grey_a3
    thickness     = 1
  </link>
</links>

<image>
  <<include etc/image.conf>>
  file* = synteny_{sp_a}__{sp_b}.png
  radius* = 1500p
  background* = white
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
""")
    return conf_file


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading gene positions for {args.sp_a} and {args.sp_b} ...")
    pos_a = load_bed_positions(args.bed_dir, args.sp_a)
    pos_b = load_bed_positions(args.bed_dir, args.sp_b)
    sizes_a = load_genome_sizes(args.bed_dir, args.sp_a)
    sizes_b = load_genome_sizes(args.bed_dir, args.sp_b)

    print(f"Loading collinear pairs ...")
    pairs = load_collinear(args.collinear_dir, args.clade, args.sp_a, args.sp_b)
    print(f"  {len(pairs)} collinear pairs loaded")

    if not pairs:
        print("No pairs — cannot generate Circos input.")
        sys.exit(0)

    print(f"Selecting top {args.max_scaffolds} scaffolds per species ...")
    pairs_filt, scaffolds_a, scaffolds_b = select_top_scaffolds(
        pairs, pos_a, pos_b, args.max_scaffolds
    )
    print(f"  {len(pairs_filt)} pairs retained on "
          f"{len(scaffolds_a)} SpA scaffolds, {len(scaffolds_b)} SpB scaffolds")

    # Assign colors to SpA scaffolds (for visual link grouping by anchor)
    sorted_scafs_a = sorted(scaffolds_a, key=lambda x: -sizes_a.get(x, 0))
    color_map_a = {scaf: COLORS[i % len(COLORS)] for i, scaf in enumerate(sorted_scafs_a)}

    # Build scaffold ID maps for Circos (dots/pipes not allowed in IDs)
    scaf_id_a = {s: f"A_{s.replace('|','_').replace('.','_')}" for s in scaffolds_a}
    scaf_id_b = {s: f"B_{s.replace('|','_').replace('.','_')}" for s in scaffolds_b}

    print("Writing karyotype files ...")
    kary_a = write_karyotype(args.sp_a, scaffolds_a, sizes_a, out_dir, "A", color_map_a)
    # SpB scaffolds get grey (links colored by SpA scaffold)
    kary_b = write_karyotype(args.sp_b, scaffolds_b, sizes_b, out_dir, "B",
                             {s: "lgrey" for s in scaffolds_b})

    print("Writing links file ...")
    links_f = write_links(pairs_filt, pos_a, pos_b, args.sp_a, args.sp_b, out_dir,
                          scaf_id_a, scaf_id_b)

    print("Writing Circos configuration ...")
    conf_f = write_circos_conf(args.sp_a, args.sp_b, kary_a, kary_b, links_f, out_dir)

    print(f"\nCircos files written to {out_dir}")
    print(f"  Karyotype A : {kary_a}")
    print(f"  Karyotype B : {kary_b}")
    print(f"  Links       : {links_f}")
    print(f"  Config      : {conf_f}")
    print(f"\nRender with: circos -conf {conf_f}")
    print("Step 07 complete.")


if __name__ == "__main__":
    main()
