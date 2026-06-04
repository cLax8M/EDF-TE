#!/usr/bin/env python3
"""
16_multicopy_te_proximity.py  — v1.0
=====================================
TEs near multicopy gene families across EDF fungi.

Motivation: Argonaute (Ago) gene expansions in Rhizophagus irregularis are
linked to nearby MULE-class TEs.  This script tests whether high-copy
orthogroups systematically harbour more flanking TE content than single-copy
orthologues, and whether MULE/DNA TEs are enriched in their vicinity.

Approach:
  1. Load OrthoFinder gene-count matrix (Orthogroups.GeneCount.tsv).
  2. For each species, split orthogroups into:
       - single_copy  : count == 1
       - multicopy    : count >= MULTI_MIN  (default 4)
  3. Extract gene coordinates from GFF and add ±FLANK_BP flanking windows.
  4. Intersect with filteredRepeats.bed (bedtools intersect).
  5. Compute TE bp and TE superfamily composition per orthogroup/copy-class.
  6. Statistical comparison: Mann-Whitney U (multi vs single) per species.
  7. Highlight Rhizophagus GL_RhiirA4 MULE enrichment in high-copy OGs.

Outputs:
  OUT_DIR/16_multicopy_te_summary.tsv           (per-species stats)
  OUT_DIR/16_multicopy_te_orthogroup.tsv        (per-orthogroup TE bp)
  OUT_DIR/16_multicopy_te_enrichment.pdf        (figure: multi vs single)
  OUT_DIR/16_mule_proximity_rhizophagus.pdf     (spotlight: Rhizophagus MULE)

Usage:
  python3 16_multicopy_te_proximity.py [--species SP1 SP2 ...]

Dependencies: python3, bedtools (in PATH or syntenyEDFTE env)
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ─── PATHS ────────────────────────────────────────────────────────────────────
OG_COUNT      = Path("/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/"
                     "Orthogroups/Orthogroups.GeneCount.tsv")
GENE_GFF_DIR  = Path("/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff")
TE_ANNOT_BASE = Path("/data/carlos/EDF-TE/TE_analysis/annotated")
OUT_DIR       = Path("/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures")

# ─── PARAMETERS ──────────────────────────────────────────────────────────────
FLANK_BP   = 5_000          # flanking window around each gene
MULTI_MIN  = 4              # minimum copies to call "multicopy"
PVAL_THRES = 0.05
HIGHLIGHT_SP = "GL_RhiirA4"  # Rhizophagus irregularis A4
MULE_PATTERN = re.compile(r"MULE|Mutator|MuDR", re.IGNORECASE)
DNA_PATTERN  = re.compile(r"^DNA", re.IGNORECASE)

# ─── HELPERS ─────────────────────────────────────────────────────────────────
def find_te_bed(sp: str) -> Path | None:
    """Locate filteredRepeats.bed for a species (handles .aa suffix)."""
    sp_clean = sp.rstrip(".aa") if sp.endswith(".aa") else sp
    for suffix in ("", "_EarlGrey"):
        d = TE_ANNOT_BASE / f"{sp_clean}{suffix}" / f"{sp_clean}_summaryFiles"
        bed = d / f"{sp_clean}.filteredRepeats.bed"
        if bed.exists():
            return bed
    return None


def find_gene_gff(sp: str) -> Path | None:
    sp_clean = sp.rstrip(".aa") if sp.endswith(".aa") else sp
    for ext in (".gff3", ".gff"):
        p = GENE_GFF_DIR / f"{sp_clean}{ext}"
        if p.exists():
            return p
    return None


def gff_to_bed(gff: Path, tmp_dir: str, flank: int, sp: str = "") -> Path | None:
    """Extract gene loci from GFF → BED with ±flank window.

    Handles three formats:
      - GFF3 with gene/mRNA features (standard)
      - GFF3 CDS-only with Parent= attribute
      - JGI GFF2 CDS-only with 'proteinId N' (Rhizophagus/Glomeromycota)
        → gene IDs become jgi|{sp_short}|{proteinId} to match OrthoFinder IDs
    """
    bed_path = Path(tmp_dir) / f"{gff.stem}.flanked.bed"
    sp_short = sp.split("_", 1)[1] if "_" in sp else sp   # "GL_RhiirA4" → "RhiirA4"

    def _parse_id(attrs: str) -> str:
        m = re.search(r'ID=([^;]+)', attrs)
        return m.group(1) if m else ""

    def _parse_parent(attrs: str) -> str:
        m = re.search(r'[Pp]arent=([^;]+)', attrs)
        return m.group(1) if m else ""

    def _parse_proteinid(attrs: str):
        m = re.search(r'proteinId\s+(\d+)', attrs)
        return m.group(1) if m else None

    # First pass: determine format
    feature_types: set = set()
    raw_lines: list = []
    has_proteinid = False
    with (gzip.open(gff, "rt") if str(gff).endswith(".gz") else open(gff)) as fh:
        for i, line in enumerate(fh):
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature_types.add(parts[2])
            raw_lines.append(parts)
            if i < 500 and "proteinId" in parts[8]:
                has_proteinid = True

    # Choose strategy
    if "gene" in feature_types:
        primary = "gene"
    elif "mRNA" in feature_types:
        primary = "mRNA"
    else:
        primary = None

    lines = []

    if primary:
        for parts in raw_lines:
            if parts[2] != primary:
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end   = int(parts[4])
            gid   = _parse_id(parts[8]) or primary
            lines.append(f"{chrom}\t{max(0, start - flank)}\t{end + flank}\t{gid}")

    elif has_proteinid:
        # JGI GFF2: group CDS by proteinId → gene span
        # gene ID = "jgi|{sp_short}|{proteinId}" matches OrthoFinder prefix
        gene_spans: dict = {}
        for parts in raw_lines:
            if parts[2] != "CDS":
                continue
            pid = _parse_proteinid(parts[8])
            if pid is None:
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end   = int(parts[4])
            gid   = f"jgi|{sp_short}|{pid}"
            if gid not in gene_spans:
                gene_spans[gid] = [chrom, start, end]
            else:
                gene_spans[gid][1] = min(gene_spans[gid][1], start)
                gene_spans[gid][2] = max(gene_spans[gid][2], end)
        for gid, (chrom, g_start, g_end) in gene_spans.items():
            lines.append(
                f"{chrom}\t{max(0, g_start - flank)}\t{g_end + flank}\t{gid}"
            )

    else:
        # GFF3 CDS-only with Parent attribute
        gene_spans2: dict = {}
        for parts in raw_lines:
            if parts[2] != "CDS":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end   = int(parts[4])
            parent = _parse_parent(parts[8])
            if not parent:
                parent = _parse_id(parts[8]) or f"{chrom}:{start}"
            if parent not in gene_spans2:
                gene_spans2[parent] = [chrom, start, end]
            else:
                gene_spans2[parent][1] = min(gene_spans2[parent][1], start)
                gene_spans2[parent][2] = max(gene_spans2[parent][2], end)
        for gid, (chrom, g_start, g_end) in gene_spans2.items():
            lines.append(
                f"{chrom}\t{max(0, g_start - flank)}\t{g_end + flank}\t{gid}"
            )

    if not lines:
        return None
    with open(bed_path, "w") as out:
        out.write("\n".join(lines) + "\n")
    return bed_path


def bedtools_intersect(gene_bed: Path, te_bed: Path, tmp_dir: str) -> dict[str, dict]:
    """Return {gene_id: {'te_bp': int, 'superfams': list[str]}} via bedtools intersect."""
    out_path = Path(tmp_dir) / "intersect.bed"
    cmd = ["bedtools", "intersect", "-a", str(gene_bed), "-b", str(te_bed),
           "-wo", "-sorted"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError:
        # Fall back without -sorted if chromosomes differ
        cmd = ["bedtools", "intersect", "-a", str(gene_bed), "-b", str(te_bed), "-wo"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    hits: dict[str, dict] = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        # cols: chrom, start, end, gene_id, te_chrom, te_start, te_end, te_class, ..., overlap_bp
        if len(parts) < 5:
            continue
        gene_id  = parts[3]
        te_class = parts[7] if len(parts) > 7 else "Unknown"
        overlap  = int(parts[-1]) if parts[-1].isdigit() else 0
        if gene_id not in hits:
            hits[gene_id] = {"te_bp": 0, "superfams": []}
        hits[gene_id]["te_bp"]       += overlap
        hits[gene_id]["superfams"].append(te_class)
    return hits


def mann_whitney(a, b):
    """Two-sided Mann-Whitney U; returns (U, p)."""
    from scipy.stats import mannwhitneyu
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    return stat, p


# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="*",
                    help="Limit to these species (use full name e.g. GL_RhiirA4)")
    args = ap.parse_args()

    # Load gene count matrix
    print("[INFO] Loading OrthoFinder gene count matrix …")
    with open(OG_COUNT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    # header[0] = "Orthogroup", rest = species labels (with .aa)
    sp_cols = {col: idx for idx, col in enumerate(header[1:], 1)}

    # Read full matrix into dict: og → array of counts
    og_counts = {}    # og_id → {sp_col: count}
    with open(OG_COUNT) as fh:
        fh.readline()    # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            og_id = parts[0]
            og_counts[og_id] = [int(x) for x in parts[1:]]

    all_species_raw = list(sp_cols.keys())
    if args.species:
        # accept with or without .aa suffix
        requested = {s.rstrip(".aa") for s in args.species}
        all_species_raw = [s for s in all_species_raw
                           if s.rstrip(".aa") in requested or s in requested]

    sp_labels = [s.replace(".aa", "") for s in all_species_raw]
    print(f"[INFO] Processing {len(sp_labels)} species")

    summary_rows  = []
    og_rows       = []
    # For the spotlight plot
    rhi_data = {"single": [], "multi": [], "mule_multi": [], "mule_single": []}

    with tempfile.TemporaryDirectory() as tmp_dir:
        for sp_raw, sp in zip(all_species_raw, sp_labels):
            col_idx = sp_cols[sp_raw] - 1   # 0-based index into count arrays

            te_bed  = find_te_bed(sp)
            gff     = find_gene_gff(sp)
            if te_bed is None or gff is None:
                print(f"  [SKIP] {sp}: missing TE bed or GFF")
                continue

            # Partition orthogroups by copy number for this species
            single_ogs = [og for og, counts in og_counts.items()
                          if counts[col_idx] == 1]
            multi_ogs  = [og for og, counts in og_counts.items()
                          if counts[col_idx] >= MULTI_MIN]

            if not multi_ogs:
                continue

            gene_bed_path = gff_to_bed(gff, tmp_dir, FLANK_BP, sp)
            if gene_bed_path is None:
                print(f"  [SKIP] {sp}: empty GFF")
                continue

            try:
                hits = bedtools_intersect(gene_bed_path, te_bed, tmp_dir)
            except Exception as e:
                print(f"  [WARN] {sp}: bedtools error: {e}")
                continue

            # We need gene_id → orthogroup mapping
            # Build from GFF: gene_id → OG
            # (OrthoFinder uses "species.aa_geneID" format in Orthogroups.tsv)
            # For simplicity: intersect all gene windows, aggregate by OG later
            # Actually: gene IDs in GFF match what OrthoFinder used in its input
            # Orthogroups.GeneCount only gives per-species counts; we need the actual
            # gene lists from Orthogroups.tsv to map gene→OG
            # Quick approach: load Orthogroups.tsv for this species
            og_tsv = OG_COUNT.parent / "Orthogroups.tsv"
            gene_to_og: dict[str, str] = {}
            if og_tsv.exists():
                try:
                    with open(og_tsv) as fh:
                        og_header = fh.readline().rstrip("\n").split("\t")
                        try:
                            sp_tsv_col = og_header.index(sp_raw)
                        except ValueError:
                            try:
                                sp_tsv_col = og_header.index(sp + ".aa")
                            except ValueError:
                                sp_tsv_col = None
                        if sp_tsv_col is not None:
                            for line in fh:
                                parts = line.rstrip("\n").split("\t")
                                og_id = parts[0]
                                if sp_tsv_col < len(parts) and parts[sp_tsv_col].strip():
                                    for full_gid in parts[sp_tsv_col].split(", "):
                                        full_gid = full_gid.strip()
                                        gene_to_og[full_gid] = og_id
                                        # Also index by jgi|sp|proteinId prefix so
                                        # GFF2-derived IDs (without gene name suffix) match
                                        if full_gid.startswith("jgi|"):
                                            jgi_parts = full_gid.split("|")
                                            if len(jgi_parts) >= 3:
                                                prefix = "|".join(jgi_parts[:3])
                                                gene_to_og[prefix] = og_id
                except Exception:
                    pass

            # Aggregate TE bp per OG class
            single_te_bps = []
            multi_te_bps  = []
            mule_in_multi = []
            mule_in_single= []

            for gene_id, hit_info in hits.items():
                og_id = gene_to_og.get(gene_id)
                if og_id is None:
                    continue
                count = og_counts.get(og_id, [0] * len(sp_labels))[col_idx]
                te_bp = hit_info["te_bp"]
                has_mule = any(MULE_PATTERN.search(sf) for sf in hit_info["superfams"])

                if count == 1:
                    single_te_bps.append(te_bp)
                    mule_in_single.append(int(has_mule))
                elif count >= MULTI_MIN:
                    multi_te_bps.append(te_bp)
                    mule_in_multi.append(int(has_mule))

            if not single_te_bps or not multi_te_bps:
                continue

            U, p = mann_whitney(multi_te_bps, single_te_bps)
            summary_rows.append({
                "species":        sp,
                "n_single":       len(single_te_bps),
                "n_multi":        len(multi_te_bps),
                "median_single":  float(np.median(single_te_bps)),
                "median_multi":   float(np.median(multi_te_bps)),
                "mw_U":           U,
                "mw_p":           p,
                "signif":         p < PVAL_THRES if p == p else False,
                "pct_mule_single": np.mean(mule_in_single) * 100 if mule_in_single else 0,
                "pct_mule_multi":  np.mean(mule_in_multi) * 100 if mule_in_multi else 0,
            })

            if sp == HIGHLIGHT_SP:
                rhi_data["single"]      = single_te_bps
                rhi_data["multi"]       = multi_te_bps
                rhi_data["mule_multi"]  = mule_in_multi
                rhi_data["mule_single"] = mule_in_single

            print(f"  {sp}: single={len(single_te_bps)} median={np.median(single_te_bps):.0f} bp | "
                  f"multi={len(multi_te_bps)} median={np.median(multi_te_bps):.0f} bp  p={p:.3g}")

    if not summary_rows:
        print("[WARN] No data to plot")
        return

    # ── Write TSV ────────────────────────────────────────────────────────────
    tsv_path = OUT_DIR / "16_multicopy_te_summary.tsv"
    keys = list(summary_rows[0].keys())
    with open(tsv_path, "w") as out:
        out.write("\t".join(keys) + "\n")
        for row in summary_rows:
            out.write("\t".join(str(row[k]) for k in keys) + "\n")
    print(f"[DONE] {tsv_path}")

    # ── Figure 1: enrichment summary across species ───────────────────────────
    signif = [r for r in summary_rows if r["signif"]]
    n_sig  = len(signif)

    fig, ax = plt.subplots(figsize=(10, 6))
    species  = [r["species"] for r in summary_rows]
    med_s    = np.array([r["median_single"] for r in summary_rows])
    med_m    = np.array([r["median_multi"]  for r in summary_rows])
    ratio    = np.where(med_s > 0, med_m / med_s, np.nan)
    colors   = ["#C0392B" if r["signif"] else "#7f8c8d" for r in summary_rows]

    x = np.arange(len(species))
    ax.bar(x, ratio, color=colors, edgecolor="white", linewidth=0.4)
    ax.axhline(1, color="black", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(species, rotation=90, fontsize=5)
    ax.set_ylabel(f"Median TE bp (multi / single), ±{FLANK_BP//1000} kb flank")
    ax.set_title(f"TE enrichment near multicopy (≥{MULTI_MIN} copies) vs single-copy genes\n"
                 f"n={len(species)} species, {n_sig} significant (p<{PVAL_THRES}, MW-U)")
    sig_patch = mpatches.Patch(color="#C0392B", label="p < 0.05")
    ns_patch  = mpatches.Patch(color="#7f8c8d", label="n.s.")
    ax.legend(handles=[sig_patch, ns_patch], loc="upper right", fontsize=8)

    fig_path = OUT_DIR / "16_multicopy_te_enrichment.pdf"
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    fig.savefig(str(fig_path).replace(".pdf", ".png"), dpi=150)
    plt.close(fig)
    print(f"[DONE] {fig_path}")

    # ── Figure 2: Rhizophagus MULE spotlight ─────────────────────────────────
    if rhi_data["single"] or rhi_data["multi"]:
        fig, axes = plt.subplots(1, 2, figsize=(9, 5))

        # Panel A: TE bp distribution
        ax = axes[0]
        data_s = np.array(rhi_data["single"])
        data_m = np.array(rhi_data["multi"])
        vp = ax.violinplot([data_s, data_m], positions=[0, 1], showmedians=True)
        for pc in vp["bodies"]:
            pc.set_alpha(0.7)
        vp["bodies"][0].set_facecolor("#3498db")
        vp["bodies"][1].set_facecolor("#e74c3c")
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Single-copy\n(n={len(data_s)})",
                             f"Multicopy (≥{MULTI_MIN})\n(n={len(data_m)})"])
        ax.set_ylabel(f"TE bp within ±{FLANK_BP//1000} kb")
        ax.set_title(f"GL_RhiirA4 TE proximity\nby gene copy number")
        if rhi_data["single"] and rhi_data["multi"]:
            _, p_rhi = mann_whitney(rhi_data["multi"], rhi_data["single"])
            ax.text(0.5, 0.97, f"p={p_rhi:.3g}", ha="center", va="top",
                    transform=ax.transAxes, fontsize=9)

        # Panel B: MULE prevalence
        ax = axes[1]
        mule_pct = [np.mean(rhi_data["mule_single"]) * 100,
                    np.mean(rhi_data["mule_multi"])  * 100]
        ax.bar([0, 1], mule_pct, color=["#3498db", "#e74c3c"],
               edgecolor="white", linewidth=0.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Single-copy\n(n={len(rhi_data['mule_single'])})",
                             f"Multicopy\n(n={len(rhi_data['mule_multi'])})"])
        ax.set_ylabel("% genes with nearby MULE TE")
        ax.set_title("MULE prevalence near gene copies\n(GL_RhiirA4 Rhizophagus irregularis A4)")

        fig.suptitle("TE-gene co-occurrence: multicopy expansion and MULE association\n"
                     "Rhizophagus irregularis A4", fontsize=10, fontweight="bold")
        fig.tight_layout()
        rhi_path = OUT_DIR / "16_mule_proximity_rhizophagus.pdf"
        fig.savefig(rhi_path, dpi=150)
        fig.savefig(str(rhi_path).replace(".pdf", ".png"), dpi=150)
        plt.close(fig)
        print(f"[DONE] {rhi_path}")
    else:
        print(f"[NOTE] Rhizophagus ({HIGHLIGHT_SP}) data not available in this run")


if __name__ == "__main__":
    main()
