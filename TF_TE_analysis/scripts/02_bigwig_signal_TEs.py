#!/usr/bin/env python3
"""
02_bigwig_signal_TEs.py
───────────────────────
For each of the 58 DAP-seq BigWigs, compute the mean normalised signal over
every TE element annotated by EarlGrey.  Aggregate signals by TE class and
produce a TF × TE-class heatmap (z-scored per TF).

Outputs
  <outdir>/bw_TE_signal_raw.tsv          ← mean signal per TF per TE element
  <outdir>/bw_TE_signal_class_mean.tsv   ← mean per TF × TE class (log2)
  <outdir>/bw_TE_signal_heatmap.pdf      ← clustermap
  <outdir>/bw_TE_top_TFs_dotplot.pdf     ← per-class top-signal TFs

Usage
  python 02_bigwig_signal_TEs.py

Requirements
  pip install pyBigWig pandas scipy matplotlib
"""

import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import zscore

# ── CONFIG ──────────────────────────────────────────────────────────────────
BW_DIR     = Path("/legserv/Temp/Carlos/EDF-TE/data/TFsDAP/bigwigs")
TE_BED     = Path("/data/carlos/EDF-TE/TE_analysis/annotated/MU_Rhimi59_EarlGrey"
                  "/MU_Rhimi59_mergedRepeats/looseMerge/MU_Rhimi59.filteredRepeats.bed")
IPR_TAB    = Path("/legserv/Temp/Carlos/EDF-TE/data/functional_interpro/MU_Rhimi59_IPR.tab")
OUTDIR     = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/02_bigwig_signal")

SKIP_TE_CLASSES = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA",
                   "snRNA", "Unknown", "Unspecified", "Other"}

# Minimum TE length to consider (bp)
MIN_TE_LEN = 200
# ────────────────────────────────────────────────────────────────────────────


def load_ipr(ipr_tab: Path) -> dict[str, str]:
    domain_to_family = {
        "IPR000232": "HSF",       "IPR001766": "ForkHead",
        "IPR007087": "Znf_C2H2",  "IPR002100": "MADS",
        "IPR011598": "bHLH_Myc",  "IPR001356": "HOX",
        "IPR007604": "CP2",        "IPR001005": "SANT_Myb",
        "IPR008422": "HOX_KN",    "IPR000679": "GATA_Znf",
        "IPR004827": "bZIP",       "IPR009071": "HMGb",
        "IPR001083": "CopperFist",
    }
    ipr_map = {}
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


def load_te_bed(te_bed: Path, skip: set, min_len: int) -> pd.DataFrame:
    """Load TE BED, filter, return DataFrame."""
    rows = []
    with open(te_bed) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            te_class = parts[3].strip()
            if te_class in skip:
                continue
            if (end - start) < min_len:
                continue
            rows.append({"chrom": chrom, "start": start, "end": end, "te_class": te_class})
    return pd.DataFrame(rows)


def protein_id_from_bw(bw_path: Path) -> str:
    m = re.search(r"Rhimi59_\d+_(\d+)_R_microsporus", bw_path.stem)
    return m.group(1) if m else bw_path.stem


def mean_signal_over_intervals(bw_path: Path, te_df: pd.DataFrame) -> np.ndarray:
    """Return array of mean BigWig values for each TE interval."""
    signals = np.zeros(len(te_df), dtype=np.float32)
    try:
        bw = pyBigWig.open(str(bw_path))
        chroms_in_bw = set(bw.chroms().keys())
        for i, row in te_df.iterrows():
            if row["chrom"] not in chroms_in_bw:
                continue
            chrom_len = bw.chroms()[row["chrom"]]
            end = min(row["end"], chrom_len)
            if row["start"] >= end:
                continue
            val = bw.stats(row["chrom"], row["start"], end, type="mean")[0]
            if val is not None:
                signals[i] = float(val)
        bw.close()
    except Exception as e:
        print(f"  WARNING: {bw_path.name}: {e}")
    return signals


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("Loading TE annotations …")
    te_df = load_te_bed(TE_BED, SKIP_TE_CLASSES, MIN_TE_LEN)
    print(f"  {len(te_df)} TE elements retained (len ≥ {MIN_TE_LEN} bp)")

    print("Loading IPR annotations …")
    ipr_map = load_ipr(IPR_TAB)

    bw_files = sorted(BW_DIR.glob("*_dedup_norm.bw"))
    print(f"Processing {len(bw_files)} BigWig files …")

    # {tf_label: mean_signal_per_TE[]}
    signal_matrix = {}

    for bwf in bw_files:
        pid = protein_id_from_bw(bwf)
        tf_family = ipr_map.get(pid, "Unknown")
        tf_label = f"{tf_family}_{pid[-4:]}"
        print(f"  {tf_label} …", end=" ", flush=True)
        sig = mean_signal_over_intervals(bwf, te_df)
        signal_matrix[tf_label] = sig
        print(f"mean={sig.mean():.3f}")

    # Build element-level DataFrame and save
    sig_df = pd.DataFrame(signal_matrix, index=te_df.index)
    sig_df["te_class"] = te_df["te_class"].values
    sig_df.to_csv(OUTDIR / "bw_TE_signal_raw.tsv", sep="\t")
    print(f"Saved: {OUTDIR}/bw_TE_signal_raw.tsv")

    # ── Aggregate by TE class (median) ──────────────────────────────────────
    tf_cols = [c for c in sig_df.columns if c != "te_class"]
    class_agg = sig_df.groupby("te_class")[tf_cols].median()

    # log2 transform (add small pseudocount)
    log2_agg = np.log2(class_agg + 1e-4)

    log2_agg.to_csv(OUTDIR / "bw_TE_signal_class_mean.tsv", sep="\t")
    print(f"Saved: {OUTDIR}/bw_TE_signal_class_mean.tsv")

    # ── Heatmap: TE class × TF (z-scored per TF) ─────────────────────────────
    _plot_heatmap(log2_agg, OUTDIR, ipr_map)
    _plot_top_TFs(log2_agg, OUTDIR)

    print("\nDone.")


def _plot_heatmap(log2_df: pd.DataFrame, outdir: Path, ipr_map: dict):
    """Clustered heatmap of z-scored log2 mean BigWig signal per TF × TE class."""
    # Z-score per TF (column)
    z_df = log2_df.apply(zscore, axis=0, ddof=1)

    # Remove TE classes with all-zero variance
    z_df = z_df.loc[z_df.std(axis=1) > 0]

    if z_df.empty:
        print("  WARNING: z_df empty, skipping heatmap")
        return

    fig_h = max(5, len(z_df) * 0.4)
    fig_w = max(10, len(z_df.columns) * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmax = np.nanpercentile(np.abs(z_df.values), 95)
    im = ax.imshow(z_df.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)

    ax.set_xticks(range(len(z_df.columns)))
    ax.set_xticklabels(z_df.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(z_df.index)))
    ax.set_yticklabels(z_df.index, fontsize=8)
    ax.set_xlabel("TF (family_protID)")
    ax.set_ylabel("TE class / superfamily")
    ax.set_title("Mean DAP-seq signal over TE classes (z-scored log₂, normalised BigWig)")
    plt.colorbar(im, ax=ax, label="z-score (log₂ signal)")
    plt.tight_layout()

    out = outdir / "bw_TE_signal_heatmap.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def _plot_top_TFs(log2_df: pd.DataFrame, outdir: Path, top_n: int = 5):
    """For each TE class, bar plot of top-N TFs by mean signal."""
    n_classes = len(log2_df.index)
    if n_classes == 0:
        return

    cols = 4
    rows = (n_classes + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.5, rows * 3))
    axes = np.array(axes).flatten()

    for i, te_class in enumerate(log2_df.index):
        ax = axes[i]
        row = log2_df.loc[te_class]
        top = row.nlargest(top_n)
        ax.barh(range(len(top)), top.values[::-1], color="steelblue")
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top.index[::-1], fontsize=7)
        ax.set_xlabel("log₂ signal", fontsize=7)
        ax.set_title(te_class, fontsize=8, fontweight="bold")
        ax.invert_xaxis()

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f"Top-{top_n} TFs by mean DAP-seq signal per TE class", fontsize=10, y=1.01)
    plt.tight_layout()
    out = outdir / "bw_TE_top_TFs_dotplot.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
