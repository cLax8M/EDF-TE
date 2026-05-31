#!/usr/bin/env python3
"""
01_peak_TE_overlap.py
─────────────────────
Overlap DAP-seq peaks (58 TFs, Rhizopus microsporus ATCC11559) with the
EarlGrey TE annotation.  For each TF:

  - Count how many peaks overlap any TE
  - Break them down by TE class/superfamily
  - Test enrichment vs genomic TE coverage (Fisher's exact test)

Outputs
  <outdir>/peak_TE_overlap_stats.tsv      ← per-TF summary
  <outdir>/peak_TE_overlap_matrix.tsv     ← TF × TE-class n_overlapping_peaks
  <outdir>/peak_TE_fraction_heatmap.pdf   ← heatmap (fraction of peaks on TEs)
  <outdir>/peak_TE_enrichment_bubble.pdf  ← bubble plot (log2OR vs TE class)

Usage
  python 01_peak_TE_overlap.py   [adjust CONFIG below if needed]

Requirements
  module load BEDTools/2.31.1-GCC-13.3.0   (or 13.2)
  pip install pandas scipy matplotlib
"""

import os
import re
import subprocess
import tempfile
from pathlib import Path
from scipy.stats import fisher_exact

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# ── CONFIG ──────────────────────────────────────────────────────────────────
PEAKS_DIR  = Path("/legserv/Temp/Carlos/EDF-TE/data/TFsDAP/peaks")
TE_BED     = Path("/data/carlos/EDF-TE/TE_analysis/annotated/MU_Rhimi59_EarlGrey"
                  "/MU_Rhimi59_mergedRepeats/looseMerge/MU_Rhimi59.filteredRepeats.bed")
IPR_TAB    = Path("/legserv/Temp/Carlos/EDF-TE/data/functional_interpro/MU_Rhimi59_IPR.tab")
GENOME_FAI = Path("/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Rhimi59.fasta.fai")
OUTDIR     = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/01_peak_overlap")

# TE classes to keep (skip simple/low-complexity repeats)
SKIP_TE_CLASSES = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA",
                   "snRNA", "Unknown", "Unspecified", "Other"}

# Minimum q-value (–log10) for a peak to be considered
MIN_QVAL = 2.0   # q < 0.01

# Bed module load command (empty string = assume already in PATH)
MODULE_LOAD = "module load BEDTools/2.31.1-GCC-13.3.0"
# ────────────────────────────────────────────────────────────────────────────


def run_bedtools(cmd: str) -> str:
    """Run a bedtools command via bash (with module load) and return stdout."""
    full_cmd = f"bash -c 'source /etc/profile.d/lmod.sh 2>/dev/null; {MODULE_LOAD}; {cmd}'"
    result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bedtools error:\n{result.stderr[:500]}")
    return result.stdout


def load_ipr(ipr_tab: Path) -> dict[str, str]:
    """Return {protein_id: short_family_label} from an IPR .tab file."""
    ipr_map = {}
    domain_to_family = {
        "IPR000232": "HSF",
        "IPR001766": "ForkHead",
        "IPR007087": "Znf_C2H2",
        "IPR002100": "MADS",
        "IPR011598": "bHLH_Myc",
        "IPR001356": "HOX",
        "IPR007604": "CP2",
        "IPR001005": "SANT_Myb",
        "IPR008422": "HOX_KN",
        "IPR000679": "GATA_Znf",
        "IPR004827": "bZIP",
        "IPR009071": "HMGb",
        "IPR001083": "CopperFist",
    }
    with open(ipr_tab) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            pid = parts[0].strip()
            ipr = parts[1].strip()
            if pid not in ipr_map and ipr in domain_to_family:
                ipr_map[pid] = domain_to_family[ipr]
    return ipr_map


def parse_peaks_to_bed(xls_path: Path, min_qval: float) -> list[tuple]:
    """Parse a MACS3 XLS peaks file → list of (chr, start, end) tuples."""
    records = []
    with open(xls_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chr\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            try:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                qval = float(parts[8])
                if qval >= min_qval:
                    records.append((chrom, start, end))
            except (ValueError, IndexError):
                continue
    return records


def genome_size(fai: Path) -> int:
    """Total genome size from .fai index."""
    total = 0
    with open(fai) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 2:
                total += int(parts[1])
    return total


def te_coverage_by_class(te_bed: Path, skip: set) -> dict[str, int]:
    """Return {te_class: total_bp} from TE BED, skipping unwanted classes."""
    cov = {}
    with open(te_bed) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            te_class = parts[3].strip()
            if te_class in skip:
                continue
            length = int(parts[2]) - int(parts[1])
            cov[te_class] = cov.get(te_class, 0) + length
    return cov


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── 1. Load TF family annotations ──────────────────────────────────────
    print("Loading IPR annotations …")
    ipr_map = load_ipr(IPR_TAB)

    # ── 2. Genome stats ─────────────────────────────────────────────────────
    print("Computing genome stats …")
    genome_bp = genome_size(GENOME_FAI)
    te_class_bp = te_coverage_by_class(TE_BED, SKIP_TE_CLASSES)
    total_te_bp = sum(te_class_bp.values())
    print(f"  Genome: {genome_bp:,} bp | Annotated TE (filtered): {total_te_bp:,} bp "
          f"({100*total_te_bp/genome_bp:.1f}%)")

    # ── 3. Process each TF ──────────────────────────────────────────────────
    peak_files = sorted(PEAKS_DIR.glob("*_peaks.xls"))
    print(f"Processing {len(peak_files)} peak files …")

    # Extract protein ID from filename: Rhimi59_2_<PROTID>_R_microsporus_...
    def protein_id(p: Path) -> str:
        m = re.search(r"Rhimi59_\d+_(\d+)_R_microsporus", p.stem)
        return m.group(1) if m else p.stem

    # Extract plate-well code from filename for short labels
    def plate_well(p: Path) -> str:
        m = re.search(r"DAPi\d+([A-H]\d+)-", p.stem)
        return m.group(1) if m else p.stem

    summary_rows = []
    # TF × TE class matrix: counts of overlapping peaks
    matrix_data = {}   # {tf_id: {te_class: count}}

    # Prepare filtered TE BED (no simple/unknown)
    te_filtered_tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
    with open(TE_BED) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and parts[3].strip() not in SKIP_TE_CLASSES:
                te_filtered_tmp.write(line)
    te_filtered_tmp.close()

    for pf in peak_files:
        pid = protein_id(pf)
        pw = plate_well(pf)
        tf_family = ipr_map.get(pid, "Unknown")
        tf_label = f"{tf_family}_{pid[-4:]}"  # e.g. "HSF_2580"

        # Parse peaks
        peaks = parse_peaks_to_bed(pf, MIN_QVAL)
        n_peaks = len(peaks)
        if n_peaks == 0:
            print(f"  WARNING: no peaks passing qval threshold in {pf.name}")
            continue

        # Write peaks to temp BED
        peaks_tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
        for c, s, e in peaks:
            peaks_tmp.write(f"{c}\t{s}\t{e}\n")
        peaks_tmp.close()

        # bedtools intersect -a peaks -b TE_BED -wa -wb (overlapping peaks with TE class info)
        intersect_out = run_bedtools(
            f"bedtools intersect -a {peaks_tmp.name} -b {te_filtered_tmp.name} -wa -wb"
        )
        os.unlink(peaks_tmp.name)

        # Count overlapping peaks (unique peak coords) and by TE class
        overlapping_peaks = set()
        te_class_counts = {}
        for line in intersect_out.splitlines():
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            peak_coord = (parts[0], parts[1], parts[2])
            overlapping_peaks.add(peak_coord)
            te_class = parts[6].strip()  # col 6 = TE class from filtered BED (0-based col 3 + 3)
            te_class_counts[te_class] = te_class_counts.get(te_class, 0) + 1

        n_overlap = len(overlapping_peaks)
        frac_overlap = n_overlap / n_peaks if n_peaks > 0 else 0

        # Fisher's exact test: peaks on TEs vs expected from TE density
        # 2×2 table:
        #   peaks_on_TE   | peaks_not_on_TE
        #   expected_on_TE| expected_not_on_TE
        # Use genomic TE fraction to estimate expected
        te_frac_genome = total_te_bp / genome_bp
        expected_on_te = n_peaks * te_frac_genome
        table = [[n_overlap, n_peaks - n_overlap],
                 [round(expected_on_te), round(n_peaks - expected_on_te)]]
        try:
            odds_ratio, pval = fisher_exact(table, alternative="greater")
        except Exception:
            odds_ratio, pval = float("nan"), float("nan")

        summary_rows.append({
            "protein_id": pid,
            "plate_well": pw,
            "tf_family": tf_family,
            "tf_label": tf_label,
            "n_peaks": n_peaks,
            "n_peaks_on_TE": n_overlap,
            "frac_peaks_on_TE": round(frac_overlap, 4),
            "genomic_TE_frac": round(te_frac_genome, 4),
            "log2_odds_ratio": round(np.log2(odds_ratio + 1e-9), 3),
            "fisher_pval": round(pval, 6),
        })
        matrix_data[tf_label] = te_class_counts
        print(f"  {tf_label}: {n_peaks} peaks, {n_overlap} on TEs "
              f"({100*frac_overlap:.1f}%), OR={odds_ratio:.2f} p={pval:.2e}")

    os.unlink(te_filtered_tmp.name)

    # ── 4. Save tables ───────────────────────────────────────────────────────
    summary_df = pd.DataFrame(summary_rows).sort_values("frac_peaks_on_TE", ascending=False)
    summary_df.to_csv(OUTDIR / "peak_TE_overlap_stats.tsv", sep="\t", index=False)
    print(f"\nSaved: {OUTDIR}/peak_TE_overlap_stats.tsv")

    # Matrix: TF × TE class
    matrix_df = pd.DataFrame(matrix_data).fillna(0).astype(int).T
    matrix_df.to_csv(OUTDIR / "peak_TE_overlap_matrix.tsv", sep="\t")
    print(f"Saved: {OUTDIR}/peak_TE_overlap_matrix.tsv")

    # ── 5. Heatmap: fraction of peaks on TE per TF × TE class ───────────────
    print("Plotting heatmap …")
    _plot_heatmap(summary_df, matrix_df, OUTDIR)

    # ── 6. Bubble plot: enrichment ───────────────────────────────────────────
    _plot_enrichment_bubble(summary_df, OUTDIR)

    print("\nDone.")


def _plot_heatmap(summary_df: pd.DataFrame, matrix_df: pd.DataFrame, outdir: Path):
    """Heatmap: rows = TFs (ordered by family), cols = TE classes."""
    if matrix_df.empty:
        return

    # Normalize by total peaks per TF to get fraction
    totals = summary_df.set_index("tf_label")["n_peaks"]
    frac_df = matrix_df.div(totals, axis=0).fillna(0)

    # Only keep TE classes with at least some signal
    frac_df = frac_df.loc[:, frac_df.max() > 0.005]

    # Sort columns by total signal
    col_order = frac_df.sum(axis=0).sort_values(ascending=False).index
    frac_df = frac_df[col_order]

    # Sort rows by TF family
    row_order = summary_df.set_index("tf_label")["tf_family"].sort_values().index
    frac_df = frac_df.reindex(row_order.intersection(frac_df.index))

    fig_h = max(6, len(frac_df) * 0.35)
    fig_w = max(8, len(frac_df.columns) * 0.65)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(frac_df.values, aspect="auto", cmap="YlOrRd",
                   vmin=0, vmax=frac_df.values.max())
    ax.set_xticks(range(len(frac_df.columns)))
    ax.set_xticklabels(frac_df.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(frac_df.index)))
    ax.set_yticklabels(frac_df.index, fontsize=8)
    ax.set_xlabel("TE class / superfamily")
    ax.set_ylabel("TF (family_protID)")
    ax.set_title("Fraction of DAP-seq peaks overlapping TE class")
    plt.colorbar(im, ax=ax, label="Fraction of peaks on TE")
    plt.tight_layout()
    out_path = outdir / "peak_TE_fraction_heatmap.pdf"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def _plot_enrichment_bubble(summary_df: pd.DataFrame, outdir: Path):
    """Bubble plot: log2 OR (x) vs –log10 p (y), size=n_peaks, color=TF family."""
    fig, ax = plt.subplots(figsize=(8, 6))
    families = summary_df["tf_family"].unique()
    cmap = plt.get_cmap("tab20")
    color_map = {f: cmap(i / max(len(families), 1)) for i, f in enumerate(sorted(families))}

    added_legends = set()
    for _, row in summary_df.iterrows():
        pval = row["fisher_pval"]
        if pval <= 0:
            pval = 1e-300
        y = -np.log10(pval)
        x = row["log2_odds_ratio"]
        sz = max(20, row["n_peaks"] / 5)
        fam = row["tf_family"]
        color = color_map[fam]
        label = fam if fam not in added_legends else "_nolegend_"
        added_legends.add(fam)
        ax.scatter(x, y, s=sz, color=color, alpha=0.75, edgecolors="k",
                   linewidths=0.4, label=label)

    ax.axhline(-np.log10(0.05), ls="--", color="grey", lw=0.8, label="p=0.05")
    ax.axvline(0, ls=":", color="grey", lw=0.8)
    ax.set_xlabel("log₂ Odds Ratio (peaks on TE vs genomic TE density)")
    ax.set_ylabel("–log₁₀ Fisher p-value")
    ax.set_title("Enrichment of DAP-seq peaks on TEs")
    ax.legend(fontsize=7, ncol=2, loc="upper left", framealpha=0.6)
    plt.tight_layout()
    out_path = outdir / "peak_TE_enrichment_bubble.pdf"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
