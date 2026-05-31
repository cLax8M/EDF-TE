#!/usr/bin/env python3
"""
04_htt_locus_conservation.py  —  Genomic locus conservation at HTT sites.

Assesses whether horizontally transferred TE copies are inserted into
syntenic loci across donor and recipient genomes.  Syntenic conservation
of insertion sites provides additional support for horizontal transfer by
ruling out convergent insertions at mutational hotspots.

Analysis:
  1. For each HTT candidate pair (donor: sp_A, recipient: sp_B):
     - Retrieve flanking 10 kb on each side of TE copy in both genomes
     - BLAST flanking sequences cross-species
     - Score syntenic conservation: fraction of flanking length with
       significant BLAST hits (≥ 50 bp, ≥ 60% identity)
  2. Stratify candidates by conservation score: conserved (≥ 30%),
     partial (10–30%), non-conserved (< 10%)
  3. Produce locus snapshots for top candidates (matplotlib — syntenic
     block representation)
  4. Summary: Monosiga/Linnaemania LINE as exemplar of conserved locus

Inputs:
  --candidates    htt_candidates.tsv from 02_htt_candidate_detection.py
  --genome_dir    Directory with per-species genome FASTAs
  --te_bed_dir    Per-species TE BED files (EarlGrey)
  --out_dir       Output directory
  --flank_kb      Flanking size in kb (default: 10)
  --n_threads     BLAST threads (default: 4)
  --top_n         Number of candidates for locus snapshots (default: 20)

Outputs:
  htt_locus_conservation.tsv
  figures/htt_locus_conservation_bar.pdf
  figures/htt_locus_snapshot_<qseqid>.pdf  (top_n candidates)
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--candidates",  required=True)
    p.add_argument("--genome_dir",  required=True)
    p.add_argument("--te_bed_dir",  required=True)
    p.add_argument("--out_dir",     default="htt_locus_out")
    p.add_argument("--flank_kb",    type=int, default=10)
    p.add_argument("--n_threads",   type=int, default=4)
    p.add_argument("--top_n",       type=int, default=20)
    return p.parse_args()


def sp_from_header(header: str) -> str:
    return header.split("__")[0] if "__" in header else header.split("_")[0]


def genome_path(genome_dir: str, sp: str) -> Path | None:
    for ext in [".fa", ".fasta", ".fna"]:
        p = Path(genome_dir) / f"{sp}{ext}"
        if p.exists():
            return p
    return None


def te_copy_coords(te_bed_dir: str, sp: str, family: str) -> pd.DataFrame:
    """Return BED rows for a given TE family in a species."""
    beds = list(Path(te_bed_dir).glob(f"{sp}*.bed"))
    if not beds:
        return pd.DataFrame()
    df = pd.read_csv(beds[0], sep="\t", header=None,
                     names=["chrom", "start", "end", "family", "score", "strand"],
                     usecols=[0, 1, 2, 3])
    return df[df["family"].str.contains(family, case=False, na=False)]


def extract_flanks_bedtools(genome_fa: Path, chrom: str, start: int, end: int,
                             flank_bp: int) -> str | None:
    """Return a FASTA string of both flanking regions concatenated."""
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as bed_f:
        # left flank
        l_start = max(0, start - flank_bp)
        bed_f.write(f"{chrom}\t{l_start}\t{start}\tleft\t0\t+\n")
        # right flank
        bed_f.write(f"{chrom}\t{end}\t{end + flank_bp}\tright\t0\t+\n")
        bed_path = bed_f.name

    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=True) as fa_f:
        ret = subprocess.run(
            ["bedtools", "getfasta", "-fi", str(genome_fa), "-bed", bed_path,
             "-fo", fa_f.name], capture_output=True)
        if ret.returncode != 0:
            return None
        return Path(fa_f.name).read_text()


def blast_flank_vs_genome(flank_fa: str, subject_fa: Path,
                           n_threads: int) -> float:
    """Return fraction of flank query length with significant BLAST hits."""
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as fq:
        fq.write(flank_fa)
        query_path = fq.name

    with tempfile.NamedTemporaryFile(suffix=".tsv", mode="w", delete=True) as out:
        cmd = ["blastn", "-query", query_path, "-subject", str(subject_fa),
               "-outfmt", "6 qseqid sseqid pident length qlen",
               "-word_size", "11", "-num_threads", str(n_threads), "-out", out.name]
        subprocess.run(cmd, capture_output=True)
        try:
            df = pd.read_csv(out.name, sep="\t", header=None,
                             names=["q", "s", "pident", "length", "qlen"])
        except Exception:
            return 0.0

    df = df[(df["pident"] >= 60) & (df["length"] >= 50)]
    if df.empty:
        return 0.0
    total_qlen = df["qlen"].max()
    covered    = df["length"].sum()
    return min(1.0, covered / total_qlen) if total_qlen > 0 else 0.0


def conservation_category(frac: float) -> str:
    if frac >= 0.30:
        return "conserved"
    if frac >= 0.10:
        return "partial"
    return "non-conserved"


# =============================================================================
# Figures
# =============================================================================

def plot_conservation_bar(result: pd.DataFrame, out_dir: Path) -> None:
    counts = result["conservation_category"].value_counts()
    cats   = ["conserved", "partial", "non-conserved"]
    vals   = [counts.get(c, 0) for c in cats]
    colors = ["#1f78b4", "#a6cee3", "#d9d9d9"]

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(cats, vals, color=colors)
    ax.set_ylabel("HTT candidates", fontsize=10)
    ax.set_title("Locus conservation at HTT sites", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "htt_locus_conservation_bar.pdf")
    plt.close(fig)
    print("Saved: htt_locus_conservation_bar.pdf")


def plot_locus_snapshot(q_id: str, q_sp: str, s_sp: str,
                         cons_frac: float, out_dir: Path) -> None:
    """Simple schematic locus diagram (schematic only — no genome coordinates)."""
    fig, axes = plt.subplots(2, 1, figsize=(8, 3), sharex=False)
    for ax, sp, label in zip(axes, [q_sp, s_sp], ["Donor", "Recipient"]):
        ax.barh(0, 10, height=0.4, color="#aaaaaa", left=0, label="Flanking region")
        ax.barh(0,  2, height=0.6, color="#e31a1c", left=4, label="TE copy")
        ax.set_xlim(0, 10)
        ax.set_yticks([])
        ax.set_title(f"{label}: {sp}", fontsize=9)
        ax.spines[["top", "right", "left"]].set_visible(False)
    axes[-1].set_xlabel("Relative genomic position (schematic)", fontsize=8)
    fig.suptitle(f"{q_id}  |  Conservation: {cons_frac:.0%}", fontsize=9)
    fig.tight_layout()
    safe_id = q_id.replace("/", "_").replace(":", "_")[:50]
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"htt_locus_snapshot_{safe_id}.pdf")
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    snap_dir = out_dir / "snapshots"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    cands = pd.read_csv(args.candidates, sep="\t")
    flank_bp = args.flank_kb * 1000

    rows = []
    for _, r in cands.iterrows():
        q_sp = r.get("q_sp", sp_from_header(r.get("qseqid", "")))
        s_sp = r.get("s_sp", sp_from_header(r.get("sseqid", "")))
        family = r.get("qseqid", "").split("__")[-1] if "__" in str(r.get("qseqid", "")) else ""

        q_genome = genome_path(args.genome_dir, q_sp)
        s_genome = genome_path(args.genome_dir, s_sp)
        if q_genome is None or s_genome is None:
            cons_frac = np.nan
        else:
            q_copies = te_copy_coords(args.te_bed_dir, q_sp, family)
            if q_copies.empty:
                cons_frac = np.nan
            else:
                # Use first matching copy
                copy = q_copies.iloc[0]
                flank_fa = extract_flanks_bedtools(
                    q_genome, copy["chrom"], copy["start"], copy["end"], flank_bp)
                if flank_fa is None:
                    cons_frac = np.nan
                else:
                    cons_frac = blast_flank_vs_genome(flank_fa, s_genome, args.n_threads)

        rows.append({
            "qseqid":                r.get("qseqid", ""),
            "sseqid":                r.get("sseqid", ""),
            "q_sp":                  q_sp,
            "s_sp":                  s_sp,
            "q_clade":               r.get("q_clade", ""),
            "s_clade":               r.get("s_clade", ""),
            "pident":                r.get("pident", np.nan),
            "conservation_fraction": cons_frac,
            "conservation_category": conservation_category(cons_frac)
                                     if not np.isnan(cons_frac) else "unknown",
        })

    result = pd.DataFrame(rows)
    result.to_csv(out_dir / "htt_locus_conservation.tsv", sep="\t", index=False)
    plot_conservation_bar(result, fig_dir)

    # Locus snapshots for top candidates by pident
    top = result.dropna(subset=["pident"]).nlargest(args.top_n, "pident")
    for _, r in top.iterrows():
        plot_locus_snapshot(str(r["qseqid"]), r["q_sp"], r["s_sp"],
                            r["conservation_fraction"] if not np.isnan(r["conservation_fraction"]) else 0,
                            snap_dir)

    print(f"Conservation scores computed for {result['conservation_fraction'].notna().sum()} "
          f"of {len(result)} candidates")
    print("Done.")


if __name__ == "__main__":
    main()
