#!/usr/bin/env python3
"""
graubove_validate.py
Grau-Bové-style contiguous expression validation for TE co-opted genes.

For each co-opted gene: checks RNA-seq coverage (>= min_cov reads, default 2)
across all CDS exon regions, using genome-mapped BAMs (STAR accepted_hits.bam).
Approximates the "contiguous expression" criterion of Grau-Bové et al. 2022
(Nature Ecol. Evol.) without requiring domain-level coordinate conversion.

GFF format auto-detected:
  - JGI GFF2 (Rhimi59, Mucci3): CDS line has 'proteinId NNNN' in attr col
  - GFF3 (Phybl4): mRNA line has 'proteinId=NNNN'; CDS lines linked via Parent

Usage:
    python3 graubove_validate.py \\
        --species MU_Rhimi59 \\
        --gff  /path/to/MU_Rhimi59.gff \\
        --domains_tsv /path/to/MU_Rhimi59.domains.tsv \\
        --coopted_tsv /path/to/MU_Rhimi59.true_TE_coopted_genes.tsv \\
        --bam_dir /path/to/chimerateTE/MU_Rhimi59/projects/chimeraTE_out \\
        --out /path/to/MU_Rhimi59.graubove_validation.tsv

Output columns:
    protein_id  gene_name  n_exons  n_exon_bp  bp_covered  frac_covered
    min_depth   median_depth  mean_depth  validated(YES/PARTIAL/NO)
    has_TE_domain  has_nonTE_domain
"""

import sys
import os
import re
import argparse
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

# TE keywords (same as co-option pipeline script 02)
TE_KEYWORDS = [
    "retrotrans", "integrase", "reverse transcriptase",
    "gag", "pol", "endonuclease", "ribonuclease h",
    "helicase", "protease", "transposase",
]


def is_te_domain(name):
    name = str(name).lower()
    return any(k in name for k in TE_KEYWORDS)


# --------------------------------------------------------------------------- #
# GFF parsing
# --------------------------------------------------------------------------- #

def parse_gff2(gff_path):
    """
    Parse JGI GFF2 format (Rhimi59, Mucci3).
    Extracts CDS features; maps proteinId -> [(chr, start0, end, strand)].
    start0 = 0-based (GFF is 1-based, we subtract 1).
    """
    gene_exons = defaultdict(list)   # proteinId -> list of (chr, start0, end, strand)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chr_, source, feat, start, end, score, strand, frame, attrs = parts
            if feat != 'CDS':
                continue
            m = re.search(r'proteinId\s+(\d+)', attrs)
            if m:
                pid = int(m.group(1))
                gene_exons[pid].append((chr_, int(start) - 1, int(end), strand))
    return gene_exons


def parse_gff3(gff_path):
    """
    Parse JGI GFF3 format (Phybl4).
    mRNA features carry proteinId=NNNN; CDS features are linked via Parent.
    Returns proteinId -> [(chr, start0, end, strand)].
    """
    # Step 1: mRNA_id -> proteinId
    mrna_to_prot = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feat = parts[2]
            if feat not in ('mRNA', 'transcript'):
                continue
            attrs = parts[8]
            id_m = re.search(r'ID=([^;]+)', attrs)
            prot_m = re.search(r'proteinId=(\d+)', attrs)
            if id_m and prot_m:
                mrna_to_prot[id_m.group(1)] = int(prot_m.group(1))

    # Step 2: CDS features -> map to proteinId via Parent
    gene_exons = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chr_, source, feat, start, end, score, strand, frame, attrs = parts
            if feat != 'CDS':
                continue
            parent_m = re.search(r'Parent=([^;]+)', attrs)
            if parent_m:
                parent = parent_m.group(1)
                # Parent may be comma-separated (multiple parents)
                for p in parent.split(','):
                    if p in mrna_to_prot:
                        pid = mrna_to_prot[p]
                        gene_exons[pid].append((chr_, int(start) - 1, int(end), strand))
    return gene_exons


def parse_gff(gff_path):
    """Auto-detect GFF2 vs GFF3 and parse accordingly."""
    # Peek at first few non-comment lines
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if 'gff-version 3' in line or 'gff-version\t3' in line:
                break
            # Check for GFF3 attributes style (key=value)
            parts = line.split('\t')
            if len(parts) >= 9 and '=' in parts[8]:
                return parse_gff3(gff_path)
            else:
                return parse_gff2(gff_path)
    return parse_gff3(gff_path)


# --------------------------------------------------------------------------- #
# Coverage analysis
# --------------------------------------------------------------------------- #

def build_merged_bam(bam_dir, tmp_dir, threads=8, samtools_bin="samtools",
                     use_star_bam=False):
    """
    Merge RNA-seq BAM files found in repXX/alignment/.
    By default uses accepted_hits.bam (ChimeraTE strand-filtered).
    With use_star_bam=True, uses the full STAR alignment
    (*_Aligned.sortedByCoord.out.bam) which has higher genome coverage.
    Returns path to merged BAM (sorted, indexed).
    """
    if use_star_bam:
        bam_files = sorted(Path(bam_dir).glob(
            "rep*/alignment/*_Aligned.sortedByCoord.out.bam"))
        bam_label = "Aligned.sortedByCoord.out.bam"
    else:
        bam_files = sorted(Path(bam_dir).glob("rep*/alignment/accepted_hits.bam"))
        bam_label = "accepted_hits.bam"
    bam_files = [str(b) for b in bam_files if Path(b).exists()]
    if not bam_files:
        raise FileNotFoundError(
            f"No {bam_label} found under {bam_dir}")
    print(f"  Merging {len(bam_files)} BAM file(s): "
          f"{[Path(b).parent.parent.name for b in bam_files]}")

    merged = os.path.join(tmp_dir, "merged_all_reps.bam")
    if len(bam_files) == 1:
        # Just symlink / copy
        subprocess.run(f"cp '{bam_files[0]}' '{merged}'", shell=True, check=True)
    else:
        bam_args = " ".join(f"'{b}'" for b in bam_files)
        subprocess.run(
            f"'{samtools_bin}' merge -@ {threads} -f '{merged}' {bam_args}",
            shell=True, check=True
        )
    subprocess.run(f"'{samtools_bin}' index '{merged}'", shell=True, check=True)
    return merged


def compute_genomecov(merged_bam, tmp_dir, bedtools_bin="bedtools"):
    """
    Run bedtools genomecov -bga to get per-base coverage across all covered bases.
    Returns path to bedgraph file.
    """
    bg = os.path.join(tmp_dir, "genomecov.bg")
    print(f"  Computing genome-wide coverage (bedtools genomecov)...")
    subprocess.run(
        f"'{bedtools_bin}' genomecov -bga -ibam '{merged_bam}' > '{bg}'",
        shell=True, check=True
    )
    return bg


def load_coverage_bedgraph(bg_path):
    """
    Load bedgraph into dict: chr -> sorted list of (start, end, depth).
    """
    print(f"  Loading coverage bedgraph...")
    cov = defaultdict(list)
    with open(bg_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            chr_, start, end, depth = parts[0], int(parts[1]), int(parts[2]), int(float(parts[3]))
            cov[chr_].append((start, end, depth))
    return cov


def query_coverage(cov, chrom, start, end):
    """
    Given the coverage dict and a genomic interval [start, end),
    return (n_bases, bp_covered_>=2, min_depth, sum_depth).
    Uses linear scan over overlapping bedgraph intervals.
    """
    intervals = cov.get(chrom, [])
    n_bases = end - start
    bp_covered = 0
    min_depth = None
    sum_depth = 0

    for (iv_start, iv_end, depth) in intervals:
        # Overlap with [start, end)
        ov_start = max(start, iv_start)
        ov_end = min(end, iv_end)
        if ov_start >= ov_end:
            continue
        bp = ov_end - ov_start
        if depth >= 2:
            bp_covered += bp
        if min_depth is None or depth < min_depth:
            min_depth = depth
        sum_depth += depth * bp

    if min_depth is None:
        min_depth = 0
    mean_depth = sum_depth / n_bases if n_bases > 0 else 0
    frac = bp_covered / n_bases if n_bases > 0 else 0
    return n_bases, bp_covered, min_depth, mean_depth, frac


def gene_coverage(gene_exons_list, cov):
    """
    Compute aggregate coverage stats for a gene (list of exon tuples).
    Merges overlapping exons first.
    Returns: (n_exons_raw, n_exon_bp, bp_covered, min_depth, median_depth, mean_depth, frac).
    """
    if not gene_exons_list:
        return 0, 0, 0, 0, 0, 0, 0.0

    # Deduplicate and merge overlapping exons
    chrom = gene_exons_list[0][0]
    strand = gene_exons_list[0][3]
    intervals = sorted(set((s, e) for (_, s, e, _) in gene_exons_list))
    merged = []
    for s, e in intervals:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append([s, e])

    n_exon_bp = sum(e - s for s, e in merged)
    total_bp_cov = 0
    all_depths = []
    global_min = None
    sum_depth = 0

    for (s, e) in merged:
        n_bases, bp_cov, min_d, mean_d, frac = query_coverage(cov, chrom, s, e)
        total_bp_cov += bp_cov
        sum_depth += mean_d * n_bases
        if global_min is None or min_d < global_min:
            global_min = min_d

        # Collect per-base depths for median (approximate from bedgraph)
        for iv in cov.get(chrom, []):
            ov_s = max(s, iv[0])
            ov_e = min(e, iv[1])
            if ov_s < ov_e:
                all_depths.extend([iv[2]] * (ov_e - ov_s))

    if global_min is None:
        global_min = 0
    mean_depth = sum_depth / n_exon_bp if n_exon_bp > 0 else 0
    median_depth = sorted(all_depths)[len(all_depths) // 2] if all_depths else 0
    frac_covered = total_bp_cov / n_exon_bp if n_exon_bp > 0 else 0.0

    return len(gene_exons_list), n_exon_bp, total_bp_cov, global_min, median_depth, mean_depth, frac_covered


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--species', required=True)
    parser.add_argument('--gff', required=True)
    parser.add_argument('--domains_tsv', required=True)
    parser.add_argument('--coopted_tsv', required=True)
    parser.add_argument('--bam_dir', required=True,
                        help='chimTE out dir with repXX/alignment/accepted_hits.bam')
    parser.add_argument('--out', required=True)
    parser.add_argument('--min_cov', type=int, default=2,
                        help='Minimum read depth for coverage (default: 2)')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--bedtools', default='bedtools',
                        help='Path to bedtools binary (default: bedtools in PATH)')
    parser.add_argument('--samtools', default='samtools',
                        help='Path to samtools binary (default: samtools in PATH)')
    parser.add_argument('--use-star-bam', action='store_true',
                        help='Use full STAR Aligned.sortedByCoord.out.bam instead of '
                             'accepted_hits.bam (useful for species with sparse coverage)')
    args = parser.parse_args()

    print(f"\n=== Grau-Bové contiguous expression validation ===")
    print(f"  Species:   {args.species}")
    print(f"  GFF:       {args.gff}")
    print(f"  min_cov:   {args.min_cov}")

    # ---- Parse GFF ----
    print("\n[1] Parsing GFF...")
    gene_exons = parse_gff(args.gff)
    print(f"    {len(gene_exons)} genes with CDS exons parsed.")

    # ---- Load domains ----
    print("\n[2] Loading domain annotations...")
    prot_te_domains = defaultdict(list)    # proteinId -> [domain_names TE]
    prot_nonte_domains = defaultdict(list)  # proteinId -> [domain_names non-TE]
    with open(args.domains_tsv) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            try:
                pid = int(parts[0])
            except ValueError:
                continue
            dname = parts[2]
            if is_te_domain(dname):
                prot_te_domains[pid].append(dname)
            else:
                prot_nonte_domains[pid].append(dname)
    print(f"    Domains loaded for {len(prot_te_domains) + len(prot_nonte_domains)} proteins.")

    # ---- Load co-opted genes ----
    print("\n[3] Loading co-opted genes...")
    coopted = []
    with open(args.coopted_tsv) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if parts and parts[0].isdigit():
                coopted.append(int(parts[0]))
    print(f"    {len(coopted)} co-opted genes to validate.")

    # ---- Build merged BAM and compute genome-wide coverage ----
    print("\n[4] Building coverage profile from RNA-seq BAMs...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        use_star = getattr(args, 'use_star_bam', False)
        merged_bam = build_merged_bam(args.bam_dir, tmp_dir, args.threads,
                                      samtools_bin=args.samtools,
                                      use_star_bam=use_star)
        bg_path = compute_genomecov(merged_bam, tmp_dir, bedtools_bin=args.bedtools)
        cov = load_coverage_bedgraph(bg_path)

        # ---- Validate each co-opted gene ----
        print("\n[5] Validating co-opted genes...")
        results = []
        for i, prot_id in enumerate(sorted(coopted)):
            if i > 0 and i % 20 == 0:
                print(f"    {i}/{len(coopted)}...", flush=True)

            exons = gene_exons.get(prot_id, [])
            has_te = len(prot_te_domains[prot_id]) > 0
            has_nonte = len(prot_nonte_domains[prot_id]) > 0

            if not exons:
                results.append({
                    'protein_id': prot_id, 'n_exons': 0, 'n_exon_bp': 0,
                    'bp_covered': 0, 'frac_covered': 0.0,
                    'min_depth': 0, 'median_depth': 0, 'mean_depth': 0.0,
                    'validated': 'NO_GFF',
                    'has_TE_domain': has_te, 'has_nonTE_domain': has_nonte,
                })
                continue

            n_exons_raw, n_exon_bp, bp_cov, min_d, med_d, mean_d, frac = \
                gene_coverage(exons, cov)

            # Validation criterion: Grau-Bové requires >= 2 reads across the
            # region between domains. We approximate as: >= 2 reads over
            # >= 95% of exon bases.
            if frac >= 0.95 and min_d >= args.min_cov:
                validated = 'YES'
            elif frac >= 0.5:
                validated = 'PARTIAL'
            else:
                validated = 'NO'

            results.append({
                'protein_id': prot_id, 'n_exons': n_exons_raw,
                'n_exon_bp': n_exon_bp, 'bp_covered': bp_cov,
                'frac_covered': round(frac, 4),
                'min_depth': min_d,
                'median_depth': med_d,
                'mean_depth': round(mean_d, 2),
                'validated': validated,
                'has_TE_domain': has_te, 'has_nonTE_domain': has_nonte,
            })

    # ---- Write output ----
    cols = ['protein_id', 'n_exons', 'n_exon_bp', 'bp_covered', 'frac_covered',
            'min_depth', 'median_depth', 'mean_depth', 'validated',
            'has_TE_domain', 'has_nonTE_domain']
    with open(args.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for row in results:
            fh.write('\t'.join(str(row[c]) for c in cols) + '\n')

    # ---- Summary ----
    n_yes = sum(1 for r in results if r['validated'] == 'YES')
    n_partial = sum(1 for r in results if r['validated'] == 'PARTIAL')
    n_no = sum(1 for r in results if r['validated'] in ('NO', 'NO_GFF'))
    total = len(results)

    print(f"\n=== Validation summary for {args.species} ===")
    print(f"  Co-opted genes:  {total}")
    print(f"  YES (>=95% exon bases covered, min_depth>={args.min_cov}): "
          f"{n_yes} ({100*n_yes/max(1,total):.1f}%)")
    print(f"  PARTIAL (50-95% exon bases covered): "
          f"{n_partial} ({100*n_partial/max(1,total):.1f}%)")
    print(f"  NO / NO_GFF: {n_no} ({100*n_no/max(1,total):.1f}%)")
    print(f"\n  Output: {args.out}")


if __name__ == '__main__':
    main()
