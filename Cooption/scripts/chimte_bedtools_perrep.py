#!/usr/bin/env python3
"""
chimte_bedtools_perrep.py
Generate per-rep chimTE-format TSVs using bedtools (no pyranges).
Uses existing intermediate BED files from chimTE's alignment step.

Usage:
    python3 chimte_bedtools_perrep.py <aln_dir> <tmp_dir> <group>

Example:
    python3 chimte_bedtools_perrep.py \\
        /path/to/chimerateTE/MU_Mucci3/projects/chimeraTE_out/rep02/alignment \\
        /path/to/chimerateTE/MU_Mucci3/projects/chimeraTE_out/tmp \\
        rep02

Output: writes TSVs to <tmp_dir>/TE-initiated-<group>.tsv,
        <tmp_dir>/TE-terminated-<group>.tsv, etc.
"""

import sys
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def run_bedtools_intersect(a, b, extra=""):
    """Run bedtools intersect -a A -b B -wa -wb and return stdout."""
    cmd = f"bedtools intersect -a '{a}' -b '{b}' -wa -wb {extra}"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [WARN] bedtools error: {r.stderr[:200]}", file=sys.stderr)
    return r.stdout


def strip_mate(name):
    """Strip /1 or /2 suffix from read name."""
    if name.endswith("/1") or name.endswith("/2"):
        return name[:-2]
    return name


def write_temp_bed(records, tmpdir):
    """Write list of (chr, start, end, name, dot, strand) to a temp BED file."""
    fh = tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", dir=tmpdir, delete=False
    )
    for rec in records:
        fh.write("\t".join(str(x) for x in rec) + "\n")
    fh.close()
    return fh.name


def make_window_bed(expressed_bed, window, direction, outpath):
    """
    Create upstream (direction='up') or downstream (direction='down') window BED
    from expressed genes BED (chr, start, end, gene_id, dot, strand).
    """
    if os.path.exists(outpath) and os.path.getsize(outpath) > 0:
        print(f"  {outpath} exists, skipping.")
        return
    lines = []
    with open(expressed_bed) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            chrom, start, end, gene_id, dot, strand = parts[:6]
            start, end = int(start), int(end)
            if direction == "up":
                if strand == "+":
                    w_start = max(0, start - window)
                    w_end = start
                else:
                    w_start = end
                    w_end = end + window
            else:  # downstream
                if strand == "+":
                    w_start = end
                    w_end = end + window
                else:
                    w_start = max(0, start - window)
                    w_end = start
            lines.append(f"{chrom}\t{w_start}\t{w_end}\t{gene_id}\t.\t{strand}\n")
    with open(outpath, "w") as fh:
        fh.writelines(lines)
    print(f"  Written {len(lines)} windows to {outpath}")


# --------------------------------------------------------------------------- #
# Chimeric TE-gene analysis functions
# --------------------------------------------------------------------------- #

def parse_gene_te_pairs(intersect_out):
    """
    Parse output of: bedtools intersect -a genes_window.bed -b chimeric_TEs.bed -wa -wb
    Cols: chr start end gene_id . strand TE_chr TE_start TE_end TE_class . TE_strand
    Returns:
      gene_te_info: dict (gene_id, te_locus) → (gene_strand, gene_window, te_class, te_strand, te_chr, te_start, te_end)
      te_regions: list of (te_chr, te_start, te_end, te_locus)  for making temp BED
    """
    gene_te_info = {}
    te_regions = []
    seen_te = set()
    for line in intersect_out.splitlines():
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        # gene cols: 0=chr, 1=start, 2=end, 3=gene_id, 4=dot, 5=strand
        # TE cols:   6=chr, 7=start, 8=end, 9=class, 10=dot, 11=strand
        gene_chr, gene_start, gene_end = parts[0], int(parts[1]), int(parts[2])
        gene_id, gene_dot, gene_strand = parts[3], parts[4], parts[5]
        te_chr, te_start, te_end = parts[6], int(parts[7]), int(parts[8])
        te_class, te_dot, te_strand = parts[9], parts[10], parts[11]
        te_locus = f"{te_chr}:{te_start}-{te_end}"
        gene_window = f"{gene_chr}:{gene_start}-{gene_end}"
        key = (gene_id, te_locus)
        if key not in gene_te_info:
            gene_te_info[key] = (gene_strand, gene_window, te_class, te_strand,
                                 te_chr, te_start, te_end)
        if te_locus not in seen_te:
            seen_te.add(te_locus)
            te_regions.append((te_chr, te_start, te_end, te_locus, ".", te_strand))
    return gene_te_info, te_regions


def build_reads_by_feature(intersect_out, key_col):
    """
    Parse bedtools intersect -wa -wb output.
    key_col: column index in the intersect output that identifies the feature
             (e.g. col 10 = TE locus name, or col 10 = gene_id in chim_exons)
    Returns: dict {feature_name → set of stripped read names}
    """
    d = defaultdict(set)
    for line in intersect_out.splitlines():
        parts = line.split("\t")
        if len(parts) <= key_col:
            continue
        read_name = strip_mate(parts[3])  # col 3 = read name
        feature_name = parts[key_col]
        d[feature_name].add(read_name)
    return d


def count_chimeric_reads(te_reads_dict, exon_reads_dict, gene_id, te_locus, gene_strand):
    """
    Count reads present in BOTH te_reads_dict[te_locus] AND exon_reads_dict[gene_id].
    te_reads_dict keys = te_locus (from bulk intersect with TE upstream BED)
    exon_reads_dict keys = gene_id (from bulk intersect with chim_exons.bed)
    """
    te_reads = te_reads_dict.get(te_locus, set())
    exon_reads = exon_reads_dict.get(gene_id, set())
    return len(te_reads & exon_reads)


def generate_te_initiated(aln_dir, tmp_dir, group, window=3000):
    """Generate TE-initiated TSV."""
    output_tsv = os.path.join(tmp_dir, f"TE-initiated-{group}.tsv")
    if os.path.exists(output_tsv) and os.path.getsize(output_tsv) > 0:
        print(f"  {output_tsv} already exists, skipping.")
        return

    genes_up = os.path.join(aln_dir, "genes_UP_window.bed")
    chimeric_tes = os.path.join(aln_dir, "chimeric_TEs.bed")
    te_reads_fwd = os.path.join(aln_dir, "TE_reads_fwd.bed")
    te_reads_rev = os.path.join(aln_dir, "TE_reads_rev.bed")
    gene_reads_rev = os.path.join(aln_dir, "gene_reads_rev.bed")
    chim_exons = os.path.join(aln_dir, "chim_exons.bed")

    for f in [genes_up, chimeric_tes, te_reads_fwd, te_reads_rev, gene_reads_rev, chim_exons]:
        if not os.path.exists(f):
            print(f"  [SKIP TE-initiated] Missing file: {f}")
            return

    print(f"  TE-initiated: intersecting upstream windows with chimeric TEs...")
    pair_out = run_bedtools_intersect(genes_up, chimeric_tes)
    if not pair_out.strip():
        print(f"  No gene-TE upstream pairs found.")
        return

    gene_te_info, te_regions = parse_gene_te_pairs(pair_out)
    print(f"  Found {len(gene_te_info)} gene-TE pairs (upstream)")

    if not te_regions:
        return

    # Write temp BED of TE upstream regions (named by te_locus for grouping)
    te_tmp_bed = write_temp_bed(te_regions, tmp_dir)

    print(f"  Intersecting TE reads with upstream TE regions (fwd)...")
    fwd_out = run_bedtools_intersect(te_reads_fwd, te_tmp_bed)
    te_reads_fwd_dict = build_reads_by_feature(fwd_out, key_col=9)  # col 10 = te_locus name

    print(f"  Intersecting TE reads with upstream TE regions (rev)...")
    rev_out = run_bedtools_intersect(te_reads_rev, te_tmp_bed)
    te_reads_rev_dict = build_reads_by_feature(rev_out, key_col=9)

    os.unlink(te_tmp_bed)

    print(f"  Intersecting gene reads with exons...")
    exon_out = run_bedtools_intersect(gene_reads_rev, chim_exons)
    exon_reads_dict = build_reads_by_feature(exon_out, key_col=9)  # col 10 = gene_id

    print(f"  Counting chimeric reads per gene-TE pair...")
    results = []
    for (gene_id, te_locus), (gene_strand, gene_window, te_class, te_strand,
                               te_chr, te_start, te_end) in gene_te_info.items():
        if gene_strand == "+":
            te_reads_dict = te_reads_fwd_dict
        else:
            te_reads_dict = te_reads_rev_dict
        cov = count_chimeric_reads(te_reads_dict, exon_reads_dict, gene_id, te_locus, gene_strand)
        if cov > 0:
            results.append(
                f"{gene_id}\t{gene_strand}\t{gene_window}\t{te_class}\t{te_strand}\t{te_locus}\t{cov}"
            )

    if results:
        with open(output_tsv, "w") as fh:
            fh.write("\n".join(results) + "\n")
        print(f"  Written {len(results)} TE-initiated pairs → {output_tsv}")
    else:
        print(f"  No TE-initiated chimeras found (all zero coverage).")


def generate_te_terminated(aln_dir, tmp_dir, group, window=3000):
    """Generate TE-terminated TSV (TE downstream of gene TES)."""
    output_tsv = os.path.join(tmp_dir, f"TE-terminated-{group}.tsv")
    if os.path.exists(output_tsv) and os.path.getsize(output_tsv) > 0:
        print(f"  {output_tsv} already exists, skipping.")
        return

    genes_down = os.path.join(aln_dir, "genes_DOWN_window.bed")
    chimeric_tes = os.path.join(aln_dir, "chimeric_TEs.bed")
    te_reads_fwd = os.path.join(aln_dir, "TE_reads_fwd.bed")
    te_reads_rev = os.path.join(aln_dir, "TE_reads_rev.bed")
    gene_reads_rev = os.path.join(aln_dir, "gene_reads_rev.bed")
    chim_exons = os.path.join(aln_dir, "chim_exons.bed")

    for f in [genes_down, chimeric_tes, te_reads_fwd, te_reads_rev, gene_reads_rev, chim_exons]:
        if not os.path.exists(f):
            print(f"  [SKIP TE-terminated] Missing file: {f}")
            return

    print(f"  TE-terminated: intersecting downstream windows with chimeric TEs...")
    pair_out = run_bedtools_intersect(genes_down, chimeric_tes)
    if not pair_out.strip():
        print(f"  No gene-TE downstream pairs found.")
        return

    gene_te_info, te_regions = parse_gene_te_pairs(pair_out)
    print(f"  Found {len(gene_te_info)} gene-TE pairs (downstream)")

    if not te_regions:
        return

    te_tmp_bed = write_temp_bed(te_regions, tmp_dir)

    fwd_out = run_bedtools_intersect(te_reads_fwd, te_tmp_bed)
    te_reads_fwd_dict = build_reads_by_feature(fwd_out, key_col=9)

    rev_out = run_bedtools_intersect(te_reads_rev, te_tmp_bed)
    te_reads_rev_dict = build_reads_by_feature(rev_out, key_col=9)

    os.unlink(te_tmp_bed)

    exon_out = run_bedtools_intersect(gene_reads_rev, chim_exons)
    exon_reads_dict = build_reads_by_feature(exon_out, key_col=9)

    results = []
    for (gene_id, te_locus), (gene_strand, gene_window, te_class, te_strand,
                               te_chr, te_start, te_end) in gene_te_info.items():
        if gene_strand == "+":
            te_reads_dict = te_reads_rev_dict  # reads for TE downstream of + gene are in rev
        else:
            te_reads_dict = te_reads_fwd_dict
        cov = count_chimeric_reads(te_reads_dict, exon_reads_dict, gene_id, te_locus, gene_strand)
        if cov > 0:
            results.append(
                f"{gene_id}\t{gene_strand}\t{gene_window}\t{te_class}\t{te_strand}\t{te_locus}\t{cov}"
            )

    if results:
        with open(output_tsv, "w") as fh:
            fh.write("\n".join(results) + "\n")
        print(f"  Written {len(results)} TE-terminated pairs → {output_tsv}")
    else:
        print(f"  No TE-terminated chimeras found.")


def generate_te_exonized(aln_dir, tmp_dir, group):
    """
    Generate TE-exonized TSVs (embedded, overlapped, intronic).
    Uses gene_TE_INSIDE-ALL_exons.bed if available, else chim_exons + chimeric_TEs intersection.
    """
    # TE-exonized_embedded: TE completely inside a gene exon with chimeric reads
    output_emb = os.path.join(tmp_dir, f"TE-exonized_embedded-{group}.tsv")
    output_ovl = os.path.join(tmp_dir, f"TE-exonized_overlapped-{group}.tsv")
    output_int = os.path.join(tmp_dir, f"TE-exonized_intronic-{group}.tsv")

    chimeric_tes = os.path.join(aln_dir, "chimeric_TEs.bed")
    chim_exons = os.path.join(aln_dir, "chim_exons.bed")
    chim_reads = os.path.join(aln_dir, "chim_reads.bed")

    for f in [chimeric_tes, chim_exons, chim_reads]:
        if not os.path.exists(f):
            print(f"  [SKIP TE-exonized] Missing file: {f}")
            return

    # Embedded: TE fully inside exon
    if not (os.path.exists(output_emb) and os.path.getsize(output_emb) > 0):
        print(f"  TE-exonized embedded: TEs fully inside exons...")
        # TEs completely within exons: TE is subset of exon
        emb_out = run_bedtools_intersect(chimeric_tes, chim_exons, "-F 1.0")
        results_emb = []
        for line in emb_out.splitlines():
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            te_chr, te_start, te_end = parts[0], int(parts[1]), int(parts[2])
            te_class, te_strand = parts[3], parts[5]
            gene_id, gene_strand = parts[9], parts[11]
            te_locus = f"{te_chr}:{te_start}-{te_end}"
            gene_coord = f"{parts[6]}:{parts[7]}-{parts[8]}"

            # Count chimeric reads supporting this TE in this exon
            te_tmp = write_temp_bed([(te_chr, te_start, te_end, te_locus, ".", te_strand)], tmp_dir)
            exon_tmp = write_temp_bed([(parts[6], int(parts[7]), int(parts[8]), gene_id, ".", gene_strand)], tmp_dir)
            reads_in_te = run_bedtools_intersect(chim_reads, te_tmp)
            te_names = {strip_mate(l.split("\t")[3]) for l in reads_in_te.splitlines() if "\t" in l}
            reads_in_exon = run_bedtools_intersect(chim_reads, exon_tmp)
            exon_names = {strip_mate(l.split("\t")[3]) for l in reads_in_exon.splitlines() if "\t" in l}
            os.unlink(te_tmp)
            os.unlink(exon_tmp)
            cov = len(te_names & exon_names)
            if cov > 0:
                results_emb.append(
                    f"{gene_id}\t{gene_strand}\t{gene_coord}\t{te_class}\t{te_strand}\t{te_locus}\t{cov}\tEmbedded"
                )
        if results_emb:
            with open(output_emb, "w") as fh:
                fh.write("\n".join(results_emb) + "\n")
            print(f"  Written {len(results_emb)} embedded pairs → {output_emb}")
        else:
            open(output_emb, "w").close()
            print(f"  No embedded TE-exonized chimeras.")

    # Overlapped: TE partially overlaps exon boundary
    if not (os.path.exists(output_ovl) and os.path.getsize(output_ovl) > 0):
        print(f"  TE-exonized overlapped: TEs partially overlapping exons...")
        # TEs overlapping but not fully inside exons: overlap > 0 but F < 1
        # Use simple overlap (any) minus embedded
        ovl_out = run_bedtools_intersect(chimeric_tes, chim_exons, "-f 0.5 -F 0.5")
        results_ovl = []
        for line in ovl_out.splitlines():
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            te_chr, te_start, te_end = parts[0], int(parts[1]), int(parts[2])
            te_class, te_strand = parts[3], parts[5]
            gene_id, gene_strand = parts[9], parts[11]
            te_locus = f"{te_chr}:{te_start}-{te_end}"
            gene_coord = f"{parts[6]}:{parts[7]}-{parts[8]}"
            te_tmp = write_temp_bed([(te_chr, te_start, te_end, te_locus, ".", te_strand)], tmp_dir)
            exon_tmp = write_temp_bed([(parts[6], int(parts[7]), int(parts[8]), gene_id, ".", gene_strand)], tmp_dir)
            reads_in_te = run_bedtools_intersect(chim_reads, te_tmp)
            te_names = {strip_mate(l.split("\t")[3]) for l in reads_in_te.splitlines() if "\t" in l}
            reads_in_exon = run_bedtools_intersect(chim_reads, exon_tmp)
            exon_names = {strip_mate(l.split("\t")[3]) for l in reads_in_exon.splitlines() if "\t" in l}
            os.unlink(te_tmp)
            os.unlink(exon_tmp)
            cov = len(te_names & exon_names)
            if cov > 0:
                results_ovl.append(
                    f"{gene_id}\t{gene_strand}\t{gene_coord}\t{te_class}\t{te_strand}\t{te_locus}\t{cov}\tOverlapped"
                )
        if results_ovl:
            with open(output_ovl, "w") as fh:
                fh.write("\n".join(results_ovl) + "\n")
            print(f"  Written {len(results_ovl)} overlapped pairs → {output_ovl}")
        else:
            open(output_ovl, "w").close()
            print(f"  No overlapped TE-exonized chimeras.")

    # Intronic: placeholder (complex, skip for now)
    if not (os.path.exists(output_int) and os.path.getsize(output_int) > 0):
        open(output_int, "w").close()
        print(f"  TE-exonized intronic: skipped (placeholder empty file).")


# --------------------------------------------------------------------------- #
# Missing intermediate file generation (for Phybl4 rep08)
# --------------------------------------------------------------------------- #

def generate_missing_intermediates(aln_dir, tmp_dir, threads=8):
    """
    Generate missing intermediate BED files for reps where prep is incomplete.
    Required inputs: accepted_hits.bam, fwd.bam (or fwd1_f.bam + fwd2_f.bam),
                     genes_total_expressed.bed, expressed_TEs.bed, tmp/TE_file.bed
    """
    # 1. Generate rev.bam if missing
    rev_bam = os.path.join(aln_dir, "rev.bam")
    if not os.path.exists(rev_bam) or os.path.getsize(rev_bam) == 0:
        print(f"  Generating rev.bam from accepted_hits.bam...")
        accepted = os.path.join(aln_dir, "accepted_hits.bam")
        rev1 = os.path.join(aln_dir, "rev1_f.bam")
        rev2 = os.path.join(aln_dir, "rev2_f.bam")
        subprocess.run(
            f"samtools view -@ {threads} -b -f 144 '{accepted}' -o '{rev1}'",
            shell=True, check=True
        )
        subprocess.run(
            f"samtools view -@ {threads} -b -f 64 -F 16 '{accepted}' -o '{rev2}'",
            shell=True, check=True
        )
        subprocess.run(
            f"samtools merge -@ {threads} -f '{rev_bam}' '{rev1}' '{rev2}'",
            shell=True, check=True
        )
        print(f"    rev.bam generated.")
    else:
        print(f"  rev.bam exists.")

    # 2. Generate rev.bed if missing
    rev_bed = os.path.join(aln_dir, "rev.bed")
    if not os.path.exists(rev_bed) or os.path.getsize(rev_bed) == 0:
        print(f"  Generating rev.bed from rev.bam...")
        subprocess.run(
            f"bedtools bamtobed -split -i '{rev_bam}' > '{rev_bed}'",
            shell=True, check=True
        )
        print(f"    rev.bed generated.")
    else:
        print(f"  rev.bed exists.")

    te_file = os.path.join(tmp_dir, "TE_file.bed")
    fwd_bed = os.path.join(aln_dir, "fwd.bed")
    genes_exp = os.path.join(aln_dir, "genes_total_expressed.bed")

    # 2.5. Generate genes_total_expressed.bed if missing (from FPKM tracking)
    if not os.path.exists(genes_exp) or os.path.getsize(genes_exp) == 0:
        print(f"  Generating genes_total_expressed.bed from FPKM data...")
        fpkm_file = os.path.join(aln_dir, "fpkm", "genes.fpkm_tracking")
        gene_coord_bed = os.path.join(tmp_dir, "gene_coord.bed")
        expressed_ids = set()
        with open(fpkm_file) as fh:
            next(fh)  # skip header
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 10 and parts[9]:
                    try:
                        if float(parts[9]) > 0:
                            expressed_ids.add(parts[0])
                    except ValueError:
                        pass
        n_written = 0
        with open(gene_coord_bed) as fin, open(genes_exp, 'w') as fout:
            for line in fin:
                parts = line.split('\t')
                if len(parts) >= 4 and parts[3].strip() in expressed_ids:
                    fout.write(line)
                    n_written += 1
        print(f"    Written {n_written} expressed genes.")

    # 2.6. Generate expressed_TEs.bed if missing
    expressed_tes = os.path.join(aln_dir, "expressed_TEs.bed")
    if not os.path.exists(expressed_tes) or os.path.getsize(expressed_tes) == 0:
        print(f"  Generating expressed_TEs.bed...")
        subprocess.run(
            f"cat '{fwd_bed}' '{rev_bed}' | bedtools intersect -a '{te_file}' -b stdin -wa -u > '{expressed_tes}'",
            shell=True, check=True
        )

    # 3. Get chimeric read names
    chim_names_file = os.path.join(aln_dir, "chim_names.txt")
    if not os.path.exists(chim_names_file):
        print(f"  Identifying chimeric reads (TE ∩ gene reads)...")
        # Get TE read names (from fwd + rev)
        fwd_te_names = os.path.join(aln_dir, "_fwd_te_names.txt")
        rev_te_names = os.path.join(aln_dir, "_rev_te_names.txt")
        all_te_names = os.path.join(aln_dir, "_all_te_names.txt")
        subprocess.run(
            f"bedtools intersect -a '{fwd_bed}' -b '{te_file}' -wa "
            f"| awk '{{n=$4; gsub(/\\/[12]$/,\"\",n); print n}}' | sort -u > '{fwd_te_names}'",
            shell=True, check=True
        )
        subprocess.run(
            f"bedtools intersect -a '{rev_bed}' -b '{te_file}' -wa "
            f"| awk '{{n=$4; gsub(/\\/[12]$/,\"\",n); print n}}' | sort -u > '{rev_te_names}'",
            shell=True, check=True
        )
        subprocess.run(
            f"sort -mu '{fwd_te_names}' '{rev_te_names}' > '{all_te_names}'",
            shell=True, check=True
        )
        # Get gene read names
        fwd_gene_names = os.path.join(aln_dir, "_fwd_gene_names.txt")
        rev_gene_names = os.path.join(aln_dir, "_rev_gene_names.txt")
        all_gene_names = os.path.join(aln_dir, "_all_gene_names.txt")
        subprocess.run(
            f"bedtools intersect -a '{fwd_bed}' -b '{genes_exp}' -wa "
            f"| awk '{{n=$4; gsub(/\\/[12]$/,\"\",n); print n}}' | sort -u > '{fwd_gene_names}'",
            shell=True, check=True
        )
        subprocess.run(
            f"bedtools intersect -a '{rev_bed}' -b '{genes_exp}' -wa "
            f"| awk '{{n=$4; gsub(/\\/[12]$/,\"\",n); print n}}' | sort -u > '{rev_gene_names}'",
            shell=True, check=True
        )
        subprocess.run(
            f"sort -mu '{fwd_gene_names}' '{rev_gene_names}' > '{all_gene_names}'",
            shell=True, check=True
        )
        # Chimeric = in both
        subprocess.run(
            f"comm -12 '{all_te_names}' '{all_gene_names}' > '{chim_names_file}'",
            shell=True, check=True
        )
        n_chim = int(subprocess.run(f"wc -l < '{chim_names_file}'", shell=True, capture_output=True, text=True).stdout.strip())
        print(f"    {n_chim} chimeric read pairs identified.")
        # Clean temp files
        for f in [fwd_te_names, rev_te_names, all_te_names, fwd_gene_names, rev_gene_names, all_gene_names]:
            try:
                os.unlink(f)
            except:
                pass
    else:
        n_chim = int(subprocess.run(f"wc -l < '{chim_names_file}'", shell=True, capture_output=True, text=True).stdout.strip())
        print(f"  chim_names.txt exists ({n_chim} chimeric reads).")

    # Early exit if no chimeric reads — create empty downstream files
    if os.path.getsize(chim_names_file) == 0:
        print(f"  No chimeric reads found; creating empty downstream files...")
        for fname in ["TE_reads_fwd.bed", "TE_reads_rev.bed", "gene_reads_fwd.bed",
                      "gene_reads_rev.bed", "chim_reads.bed", "chim_exons.bed", "chimeric_TEs.bed"]:
            path = os.path.join(aln_dir, fname)
            if not os.path.exists(path):
                Path(path).touch()
        make_window_bed(genes_exp, 3000, "up", os.path.join(aln_dir, "genes_UP_window.bed"))
        make_window_bed(genes_exp, 3000, "down", os.path.join(aln_dir, "genes_DOWN_window.bed"))
        print(f"  All intermediate files ready (no chimeric reads).")
        return

    # 4. Generate TE_reads_fwd.bed and TE_reads_rev.bed if missing
    te_reads_fwd = os.path.join(aln_dir, "TE_reads_fwd.bed")
    te_reads_rev = os.path.join(aln_dir, "TE_reads_rev.bed")

    if not os.path.exists(te_reads_fwd) or os.path.getsize(te_reads_fwd) == 0:
        print(f"  Generating TE_reads_fwd.bed...")
        subprocess.run(
            f"bedtools intersect -a '{fwd_bed}' -b '{te_file}' -wa "
            f"| awk 'NR==FNR{{a[$1]=1; next}} "
            f"{{n=$4; gsub(/\\/[12]$/,\"\",n); if(n in a) print}}' "
            f"'{chim_names_file}' - > '{te_reads_fwd}'",
            shell=True, check=True
        )

    if not os.path.exists(te_reads_rev) or os.path.getsize(te_reads_rev) == 0:
        print(f"  Generating TE_reads_rev.bed...")
        subprocess.run(
            f"bedtools intersect -a '{rev_bed}' -b '{te_file}' -wa "
            f"| awk 'NR==FNR{{a[$1]=1; next}} "
            f"{{n=$4; gsub(/\\/[12]$/,\"\",n); if(n in a) print}}' "
            f"'{chim_names_file}' - > '{te_reads_rev}'",
            shell=True, check=True
        )

    # 5. Generate gene_reads_fwd.bed and gene_reads_rev.bed
    gene_reads_fwd = os.path.join(aln_dir, "gene_reads_fwd.bed")
    gene_reads_rev = os.path.join(aln_dir, "gene_reads_rev.bed")

    if not os.path.exists(gene_reads_fwd) or os.path.getsize(gene_reads_fwd) == 0:
        print(f"  Generating gene_reads_fwd.bed...")
        subprocess.run(
            f"bedtools intersect -a '{fwd_bed}' -b '{genes_exp}' -wa "
            f"| awk 'NR==FNR{{a[$1]=1; next}} "
            f"{{n=$4; gsub(/\\/[12]$/,\"\",n); if(n in a) print}}' "
            f"'{chim_names_file}' - > '{gene_reads_fwd}'",
            shell=True, check=True
        )

    if not os.path.exists(gene_reads_rev) or os.path.getsize(gene_reads_rev) == 0:
        print(f"  Generating gene_reads_rev.bed...")
        subprocess.run(
            f"bedtools intersect -a '{rev_bed}' -b '{genes_exp}' -wa "
            f"| awk 'NR==FNR{{a[$1]=1; next}} "
            f"{{n=$4; gsub(/\\/[12]$/,\"\",n); if(n in a) print}}' "
            f"'{chim_names_file}' - > '{gene_reads_rev}'",
            shell=True, check=True
        )

    # 6. Generate chim_reads.bed
    accepted_bed = os.path.join(aln_dir, "accepted_hits.bed")
    chim_reads = os.path.join(aln_dir, "chim_reads.bed")
    if not os.path.exists(chim_reads) or os.path.getsize(chim_reads) == 0:
        print(f"  Generating chim_reads.bed...")
        subprocess.run(
            f"awk 'NR==FNR{{a[$1]=1; next}} "
            f"{{n=$4; gsub(/\\/[12]$/,\"\",n); if(n in a) print}}' "
            f"'{chim_names_file}' '{accepted_bed}' > '{chim_reads}'",
            shell=True, check=True
        )

    # 7. Generate chim_exons.bed
    exon_file = os.path.join(tmp_dir, "exon_file.bed")
    chim_exons = os.path.join(aln_dir, "chim_exons.bed")
    if not os.path.exists(chim_exons) or os.path.getsize(chim_exons) == 0:
        print(f"  Generating chim_exons.bed...")
        subprocess.run(
            f"bedtools intersect -a '{exon_file}' -b '{chim_reads}' -wa -u > '{chim_exons}'",
            shell=True, check=True
        )

    # 8. Generate chimeric_TEs.bed
    expressed_tes = os.path.join(aln_dir, "expressed_TEs.bed")
    chimeric_tes = os.path.join(aln_dir, "chimeric_TEs.bed")
    if not os.path.exists(chimeric_tes) or os.path.getsize(chimeric_tes) == 0:
        print(f"  Generating chimeric_TEs.bed...")
        subprocess.run(
            f"bedtools intersect -a '{expressed_tes}' -b '{chim_reads}' -wa -u > '{chimeric_tes}'",
            shell=True, check=True
        )

    # 9. Generate genes_UP_window.bed
    genes_up = os.path.join(aln_dir, "genes_UP_window.bed")
    if not os.path.exists(genes_up) or os.path.getsize(genes_up) == 0:
        print(f"  Generating genes_UP_window.bed...")
        make_window_bed(genes_exp, 3000, "up", genes_up)

    # 10. Generate genes_DOWN_window.bed
    genes_down = os.path.join(aln_dir, "genes_DOWN_window.bed")
    if not os.path.exists(genes_down) or os.path.getsize(genes_down) == 0:
        print(f"  Generating genes_DOWN_window.bed...")
        make_window_bed(genes_exp, 3000, "down", genes_down)

    print(f"  All intermediate files ready.")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    aln_dir = sys.argv[1].rstrip("/")
    tmp_dir = sys.argv[2].rstrip("/")
    group = sys.argv[3]

    print(f"\n=== chimte_bedtools_perrep for {group} ===")
    print(f"  aln_dir: {aln_dir}")
    print(f"  tmp_dir: {tmp_dir}")

    # Check if we need to generate missing intermediate files
    te_reads_fwd = os.path.join(aln_dir, "TE_reads_fwd.bed")
    if not os.path.exists(te_reads_fwd) or os.path.getsize(te_reads_fwd) == 0:
        print(f"\n[PREP] Generating missing intermediate files...")
        generate_missing_intermediates(aln_dir, tmp_dir)

    # Generate genes_DOWN_window.bed if missing (needed for TE-terminated)
    genes_down = os.path.join(aln_dir, "genes_DOWN_window.bed")
    if not os.path.exists(genes_down) or os.path.getsize(genes_down) == 0:
        genes_exp = os.path.join(aln_dir, "genes_total_expressed.bed")
        if os.path.exists(genes_exp):
            print(f"\n[PREP] Generating genes_DOWN_window.bed...")
            make_window_bed(genes_exp, 3000, "down", genes_down)

    # Remove partial TE-initiated_tmp.tsv to avoid duplicate entries
    partial_tmp = os.path.join(aln_dir, "TE-initiated_tmp.tsv")
    if os.path.exists(partial_tmp):
        print(f"\n[CLEANUP] Removing partial {partial_tmp}")
        os.remove(partial_tmp)

    # Generate per-rep TSVs
    print(f"\n[ANALYSIS] Generating chimera TSVs for {group}...")
    generate_te_initiated(aln_dir, tmp_dir, group)
    generate_te_terminated(aln_dir, tmp_dir, group)
    generate_te_exonized(aln_dir, tmp_dir, group)

    print(f"\n=== Done: {group} ===")
    print(f"TSVs written to {tmp_dir}/")


if __name__ == "__main__":
    main()
