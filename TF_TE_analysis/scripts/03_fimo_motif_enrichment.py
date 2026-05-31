#!/usr/bin/env python3
"""
03_fimo_motif_enrichment.py
───────────────────────────
Test whether TF binding motifs (MEME format) are enriched in TE sequences
relative to the genomic background.

Strategy
  1. Extract TE sequences per class using bedtools getfasta
  2. Shuffle those sequences 3× using bedtools shuffle → background
  3. Run FIMO (from MEME suite) for each motif file against
     - TE sequences
     - shuffled background sequences
  4. Fisher's exact test: (hits in TE | no-hits in TE) vs
                          (hits in background | no-hits in background)
  5. Output enrichment table + bubble plot

Outputs
  <outdir>/fimo_results/         ← raw FIMO TSV files
  <outdir>/fimo_enrichment.tsv   ← enrichment log2OR, p-value per motif × TE class
  <outdir>/fimo_enrichment_bubble.pdf
  <outdir>/fimo_enrichment_heatmap.pdf

Usage
  python 03_fimo_motif_enrichment.py [--outdir /path/to/out]

Requirements
  module load BEDTools/2.31.1-GCC-13.3.0
  module load MEME/5.5.7-gompi-2023b
  pip install pandas scipy matplotlib biopython

Note: This script is CPU-intensive and is best submitted to the cluster.
      A companion SLURM script (run_fimo_enrichment.sh) is provided.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from scipy.stats import fisher_exact

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── CONFIG ──────────────────────────────────────────────────────────────────
MOTIF_DIR  = Path("/legserv/Temp/Carlos/EDF-TE/data/TFsDAP/motifs")
TE_BED     = Path("/data/carlos/EDF-TE/TE_analysis/annotated/MU_Rhimi59_EarlGrey"
                  "/MU_Rhimi59_mergedRepeats/looseMerge/MU_Rhimi59.filteredRepeats.bed")
GENOME_FA  = Path("/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Rhimi59.fasta")
GENOME_FAI = Path("/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Rhimi59.fasta.fai")
GFF_PATH   = Path("/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff/MU_Rhimi59.gff")
OUTDIR     = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/03_fimo_enrichment")

SKIP_TE_CLASSES = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA",
                   "snRNA", "Unknown", "Unspecified", "Other"}

MIN_TE_LEN = 200   # bp
N_SHUFFLE  = 3     # multiply background size
FIMO_THRESH = 1e-4 # FIMO p-value threshold

MODULE_LOAD = ("source /etc/profile.d/lmod.sh 2>/dev/null; "
               "module load BEDTools/2.31.1-GCC-13.3.0; "
               "module load MEME/5.5.7-gompi-2023b")
# ────────────────────────────────────────────────────────────────────────────


def run_cmd(cmd: str, desc: str = "") -> str:
    full = f"bash -c '{MODULE_LOAD}; {cmd}'"
    r = subprocess.run(full, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Command failed [{desc}]:\n{r.stderr[:600]}")
    return r.stdout


def extract_te_fastas_by_class(te_bed: Path, genome_fa: Path,
                                genome_fai: Path, skip: set, min_len: int,
                                out_root: Path) -> dict[str, Path]:
    """
    Write one FASTA per TE class (filtered), return {te_class: fasta_path}.
    Also writes a per-class BED for bedtools shuffle.
    """
    # Group TE entries by class
    class_beds = {}
    all_entries = []
    with open(te_bed) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            chrom, start, end = p[0], int(p[1]), int(p[2])
            te_class = p[3].strip()
            if te_class in skip or (end - start) < min_len:
                continue
            all_entries.append((te_class, chrom, start, end))

    # Write BED files per class
    beds_by_class = {}
    for te_class, chrom, start, end in all_entries:
        if te_class not in beds_by_class:
            bed_path = out_root / f"te_{_safe_name(te_class)}.bed"
            beds_by_class[te_class] = open(bed_path, "w")
        beds_by_class[te_class].write(f"{chrom}\t{start}\t{end}\t{te_class}\n")

    for cls, fh in beds_by_class.items():
        fh.close()

    # Extract FASTA per class
    fa_by_class = {}
    for te_class, bed_fh in beds_by_class.items():
        bed_path = out_root / f"te_{_safe_name(te_class)}.bed"
        fa_path  = out_root / f"te_{_safe_name(te_class)}.fa"
        if not fa_path.exists() or fa_path.stat().st_size == 0:
            run_cmd(
                f"bedtools getfasta -fi {genome_fa} -bed {bed_path} -fo {fa_path} -name",
                f"getfasta {te_class}"
            )
        fa_by_class[te_class] = fa_path

    return fa_by_class


def shuffle_fasta(bed_path: Path, genome_fa: Path, genome_fai: Path,
                  n_shuffle: int, out_path: Path) -> Path:
    """
    Generate shuffled background sequences by shuffling BED n times and
    extracting FASTA.
    """
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    # Write combined shuffled BED
    shuffled_bed = out_path.with_suffix(".shuffled.bed")
    with open(shuffled_bed, "w") as fout:
        for _ in range(n_shuffle):
            out = run_cmd(
                f"bedtools shuffle -i {bed_path} -g {genome_fai} -seed 42",
                "shuffle"
            )
            fout.write(out)
    # Extract FASTA
    run_cmd(
        f"bedtools getfasta -fi {genome_fa} -bed {shuffled_bed} -fo {out_path} -name",
        "getfasta shuffle"
    )
    return out_path


def count_fasta_seqs(fasta_path: Path) -> int:
    """Count sequences in a FASTA file."""
    n = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def run_fimo(motif_file: Path, fasta_path: Path, out_dir: Path,
             threshold: float = FIMO_THRESH) -> Path:
    """Run FIMO and return path to fimo.tsv. Returns None if output exists."""
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / "fimo.tsv"
    if tsv_path.exists() and tsv_path.stat().st_size > 0:
        return tsv_path
    run_cmd(
        f"fimo --thresh {threshold} --oc {out_dir} {motif_file} {fasta_path}",
        f"fimo {motif_file.stem}"
    )
    return tsv_path if tsv_path.exists() else None


def count_fimo_hits(tsv_path: Path) -> tuple[int, int]:
    """Return (n_seqs_with_hit, n_seqs_total) from a FIMO output TSV."""
    if tsv_path is None or not tsv_path.exists():
        return 0, 0
    seqs_with_hit = set()
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("motif"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                # sequence_name is col 2 (0-indexed)
                # Strip position from name like "seq::scaffold_1:100-200(+)"
                seqs_with_hit.add(parts[2])
    return len(seqs_with_hit), len(seqs_with_hit)


def fimo_hit_sequences(tsv_path: Path, n_total_seqs: int) -> tuple[int, int]:
    """
    Return (n_seqs_with any hit, n_total_seqs).
    """
    if tsv_path is None or not tsv_path.exists():
        return 0, n_total_seqs
    seqs_hit = set()
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("motif_id"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                seqs_hit.add(parts[2])
    return len(seqs_hit), n_total_seqs


def motif_family_from_name(meme_name: str) -> str:
    """Extract short family label from motif filename."""
    for prefix in ["HLH_Myc", "GATA_Znf", "HOX_Myb", "HOX_KN", "bZip", "bzip",
                   "HMGb", "Copf", "Cp2", "MADS", "HSF", "HOX", "HLH", "FK", "Znf"]:
        if meme_name.lower().startswith(prefix.lower()):
            return prefix
    return meme_name


def list_te_classes_and_exit(outdir: Path):
    """Print filtered TE classes, save to te_class_list.txt, and exit."""
    outdir.mkdir(parents=True, exist_ok=True)
    classes = set()
    with open(TE_BED) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            te_class = p[3].strip()
            length = int(p[2]) - int(p[1])
            if te_class not in SKIP_TE_CLASSES and length >= MIN_TE_LEN:
                classes.add(te_class)
    sorted_classes = sorted(classes)
    out_file = outdir / "te_class_list.txt"
    with open(out_file, "w") as fh:
        for cls in sorted_classes:
            fh.write(cls + "\n")
    print(f"Found {len(sorted_classes)} TE classes:")
    for i, cls in enumerate(sorted_classes):
        print(f"  [{i}] {cls}")
    print(f"\nSaved to: {out_file}")
    print(f"Use --array=0-{len(sorted_classes)-1} in SLURM script")
    sys.exit(0)


def merge_results_and_plot(outdir: Path):
    """Merge per-class enrichment TSVs and generate final figures."""
    fimo_root = outdir / "fimo_results"
    per_class_tsvs = sorted(outdir.glob("fimo_enrichment_class_*.tsv"))
    if not per_class_tsvs:
        print("No per-class TSV files found. Run individual tasks first.")
        sys.exit(1)
    dfs = [pd.read_csv(tsv, sep="\t") for tsv in per_class_tsvs]
    merged = pd.concat(dfs, ignore_index=True)
    out_tsv = outdir / "fimo_enrichment.tsv"
    merged.to_csv(out_tsv, sep="\t", index=False)
    print(f"Merged {len(per_class_tsvs)} files → {out_tsv} ({len(merged)} rows)")
    _plot_bubble(merged, outdir)
    _plot_heatmap(merged, outdir)
    sys.exit(0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=str(OUTDIR))
    parser.add_argument("--te-class", default=None,
                        help="Process only this TE class (for SLURM array use)")
    parser.add_argument("--list-classes", action="store_true",
                        help="List available TE classes and exit")
    parser.add_argument("--merge-only", action="store_true",
                        help="Merge per-class results and generate figures, then exit")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.list_classes:
        list_te_classes_and_exit(outdir)

    if args.merge_only:
        merge_results_and_plot(outdir)

    fimo_root = outdir / "fimo_results"
    fimo_root.mkdir(parents=True, exist_ok=True)
    seq_dir = outdir / "sequences"
    seq_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Extract TE sequences per class ────────────────────────────────────
    print("Extracting TE sequences per class …")
    fa_by_class = extract_te_fastas_by_class(
        TE_BED, GENOME_FA, GENOME_FAI, SKIP_TE_CLASSES, MIN_TE_LEN, seq_dir
    )
    if args.te_class:
        # Filter to single class for SLURM array
        if args.te_class not in fa_by_class:
            print(f"TE class '{args.te_class}' not found. Available: {list(fa_by_class)}")
            sys.exit(1)
        fa_by_class = {args.te_class: fa_by_class[args.te_class]}

    print(f"  TE classes: {sorted(fa_by_class)}")

    # ── 2. Build shuffled background per class ───────────────────────────────
    print("Building shuffled backgrounds …")
    bg_by_class = {}
    for te_class, fa_path in fa_by_class.items():
        bed_path = seq_dir / f"te_{_safe_name(te_class)}.bed"
        bg_fa    = seq_dir / f"bg_{_safe_name(te_class)}.fa"
        bg_by_class[te_class] = shuffle_fasta(bed_path, GENOME_FA, GENOME_FAI,
                                               N_SHUFFLE, bg_fa)

    # ── 3. Run FIMO for each motif × TE class ───────────────────────────────
    motif_files = sorted(MOTIF_DIR.glob("*.meme")) + sorted(MOTIF_DIR.glob("*.MEME"))
    print(f"Running FIMO: {len(motif_files)} motif files × {len(fa_by_class)} TE classes …")

    results = []
    for te_class, fa_path in sorted(fa_by_class.items()):
        n_te_seqs = count_fasta_seqs(fa_path)
        bg_fa = bg_by_class[te_class]
        n_bg_seqs = count_fasta_seqs(bg_fa)
        print(f"  TE class: {te_class} ({n_te_seqs} seqs, bg={n_bg_seqs})")

        for mf in motif_files:
            motif_id = mf.stem
            motif_family = motif_family_from_name(motif_id)
            safe_cls = _safe_name(te_class)

            # FIMO on TE sequences
            te_fimo_dir = fimo_root / f"{safe_cls}__{motif_id}__te"
            te_tsv = run_fimo(mf, fa_path, te_fimo_dir)
            te_hits, _ = fimo_hit_sequences(te_tsv, n_te_seqs)

            # FIMO on background
            bg_fimo_dir = fimo_root / f"{safe_cls}__{motif_id}__bg"
            bg_tsv = run_fimo(mf, bg_fa, bg_fimo_dir)
            bg_hits, _ = fimo_hit_sequences(bg_tsv, n_bg_seqs)

            # Fisher's exact test
            # [hits_TE, no_hit_TE]
            # [hits_bg, no_hit_bg]
            table = [[te_hits, n_te_seqs - te_hits],
                     [bg_hits, n_bg_seqs - bg_hits]]
            try:
                odds_ratio, pval = fisher_exact(table, alternative="greater")
                log2or = np.log2(odds_ratio + 1e-9)
            except Exception:
                odds_ratio, pval, log2or = float("nan"), float("nan"), float("nan")

            results.append({
                "te_class": te_class,
                "motif_id": motif_id,
                "motif_family": motif_family,
                "n_te_seqs": n_te_seqs,
                "te_hits": te_hits,
                "te_hit_frac": round(te_hits / max(n_te_seqs, 1), 4),
                "n_bg_seqs": n_bg_seqs,
                "bg_hits": bg_hits,
                "bg_hit_frac": round(bg_hits / max(n_bg_seqs, 1), 4),
                "odds_ratio": round(odds_ratio, 4),
                "log2_odds_ratio": round(log2or, 4),
                "fisher_pval": round(pval, 8),
            })

    # ── 4. Save results ───────────────────────────────────────────────────────
    res_df = pd.DataFrame(results)
    if args.te_class:
        # Save per-class TSV for later merging
        safe = _safe_name(args.te_class)
        out_tsv = outdir / f"fimo_enrichment_class_{safe}.tsv"
        res_df.to_csv(out_tsv, sep="\t", index=False)
        print(f"Saved: {out_tsv}")
    else:
        # Single-run mode: save combined and plot immediately
        out_tsv = outdir / "fimo_enrichment.tsv"
        res_df.to_csv(out_tsv, sep="\t", index=False)
        print(f"Saved: {out_tsv}")
        _plot_bubble(res_df, outdir)
        _plot_heatmap(res_df, outdir)

    print("\nDone.")


def _safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", s)


def _plot_bubble(res_df: pd.DataFrame, outdir: Path):
    """Bubble plot: log2OR (x) vs –log10p (y), color=motif family, size=te_hit_frac."""
    if res_df.empty:
        return
    # Keep only significant (p < 0.1) or OR > 1
    plot_df = res_df[res_df["odds_ratio"] > 1].copy()
    if plot_df.empty:
        plot_df = res_df.copy()

    families = sorted(plot_df["motif_family"].unique())
    cmap = plt.get_cmap("tab20")
    color_map = {f: cmap(i / max(len(families), 1)) for i, f in enumerate(families)}

    te_classes = sorted(plot_df["te_class"].unique())
    n_cls = len(te_classes)
    cols = min(3, n_cls)
    rows = (n_cls + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4))
    axes = np.array(axes).flatten()

    for i, te_class in enumerate(te_classes):
        ax = axes[i]
        sub = plot_df[plot_df["te_class"] == te_class]
        added = set()
        for _, row in sub.iterrows():
            pval = max(row["fisher_pval"], 1e-300)
            y = -np.log10(pval)
            x = row["log2_odds_ratio"]
            sz = max(15, row["te_hit_frac"] * 500)
            fam = row["motif_family"]
            label = fam if fam not in added else "_nolegend_"
            added.add(fam)
            ax.scatter(x, y, s=sz, color=color_map[fam], alpha=0.75,
                       edgecolors="k", linewidths=0.3, label=label)
            if y > 1.3 and x > 0.5:
                ax.annotate(row["motif_id"], (x, y), fontsize=5,
                            xytext=(2, 2), textcoords="offset points")
        ax.axhline(-np.log10(0.05), ls="--", color="grey", lw=0.7)
        ax.axvline(0, ls=":", color="grey", lw=0.7)
        ax.set_title(te_class, fontsize=9, fontweight="bold")
        ax.set_xlabel("log₂ OR", fontsize=7)
        ax.set_ylabel("–log₁₀ p", fontsize=7)
        if len(added) > 0:
            ax.legend(fontsize=5, ncol=2, loc="upper left", framealpha=0.5)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("Motif enrichment in TE classes (FIMO vs shuffled background)", fontsize=11)
    plt.tight_layout()
    out = outdir / "fimo_enrichment_bubble.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def _plot_heatmap(res_df: pd.DataFrame, outdir: Path):
    """Heatmap: motif family × TE class, color=mean log2OR (filtered p<0.05)."""
    if res_df.empty:
        return

    # Aggregate: median log2OR per (motif_family, te_class)
    sig_df = res_df.copy()
    sig_df.loc[sig_df["fisher_pval"] > 0.05, "log2_odds_ratio"] = 0
    pivot = sig_df.pivot_table(
        index="motif_family", columns="te_class",
        values="log2_odds_ratio", aggfunc="median"
    ).fillna(0)

    pivot = pivot.loc[(pivot.abs() > 0).any(axis=1)]

    if pivot.empty:
        print("  No significant enrichments for heatmap")
        return

    fig_h = max(4, len(pivot) * 0.4)
    fig_w = max(6, len(pivot.columns) * 0.7)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmax = max(1.0, np.nanpercentile(np.abs(pivot.values), 95))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=9)
    ax.set_xlabel("TE class")
    ax.set_ylabel("TF motif family")
    ax.set_title("Median log₂OR of motif enrichment in TEs\n(cells with p > 0.05 set to 0)")
    plt.colorbar(im, ax=ax, label="log₂ Odds Ratio")
    plt.tight_layout()
    out = outdir / "fimo_enrichment_heatmap.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
