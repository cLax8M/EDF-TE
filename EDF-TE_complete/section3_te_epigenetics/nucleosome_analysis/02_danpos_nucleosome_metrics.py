#!/usr/bin/env python3
"""
02_danpos_nucleosome_metrics.py  —  Per-nucleosome DANPOS metrics.

For every nucleosome called by DANPOS ``dpos`` (``positions.xls``) in the
in-vivo MNase / MAINE datasets, classifies the dyad into one of three
genomic contexts — TE locus (dyad inside an EarlGrey-annotated TE),
protein-coding gene body (dyad inside a GFF ``gene`` feature), or intergenic
(neither) — and compares three per-nucleosome metrics:

  (i)   fuzziness_score  — positional SD of the read population (bp);
                          lower = more precisely positioned
  (ii)  smt_value        — smoothed peak height; higher = more nucleosome-like
  (iii) nrl_bp           — distance to the next dyad on the same chromosome
                          (nucleosome repeat length proxy; 100–400 bp filter)

Statistics: two-sided Mann–Whitney U tests for TE vs Gene and TE vs Inter,
with Benjamini–Hochberg FDR correction across all comparisons, plus
Cliff's delta effect sizes.

Inputs:
  --positions   DANPOS positions.xls
  --te_bed      EarlGrey filteredRepeats.bed
  --gff         gene annotation GFF/GFF3
  --label       sample label (e.g. AS_Ncrass2_MNase)
  --out_dir     output directory

Outputs:
  <label>_nucleosomes.tsv.gz    — per-nucleosome metrics + context
  <label>_summary.tsv           — context medians + Mann–Whitney U + BH-FDR
  figures/<label>_boxplots.pdf
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent / "helpers"))
from nucleosome_common import (  # noqa: E402
    load_danpos_positions,
    add_nrl,
    load_te_bed,
    te_intervals_by_chrom,
    load_gene_intervals,
    bh_fdr,
    cliffs_delta,
    style_ax,
)

CONTEXT_ORDER = ["TE", "Gene", "Inter"]
CONTEXT_COLORS = {"TE": "#d62728", "Gene": "#1f77b4", "Inter": "#999999"}
METRICS = {
    "fuzziness_score": "Fuzziness score (bp)",
    "smt_value": "Peak height (smt_value)",
    "nrl_bp": "NRL (bp)",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--positions", required=True, help="DANPOS positions.xls")
    p.add_argument("--te_bed",    required=True)
    p.add_argument("--gff",       required=True)
    p.add_argument("--label",     required=True)
    p.add_argument("--out_dir",   default="danpos_out")
    return p.parse_args()


def classify_context(df, te_dict, gene_dict):
    """Vectorised context classification via searchsorted per chromosome."""
    ctx = np.array(["Inter"] * len(df), dtype=object)
    te_hit = np.zeros(len(df), dtype=bool)
    gene_hit = np.zeros(len(df), dtype=bool)
    for chrom, grp in df.groupby("chr"):
        pos = grp["smt_pos"].to_numpy()
        idx = grp.index

        for label, ivs, arr in (("te", te_dict, te_hit),
                                ("gene", gene_dict, gene_hit)):
            iv = ivs.get(chrom)
            if iv is None or len(iv) == 0:
                continue
            j = np.searchsorted(iv[:, 0], pos, side="right") - 1
            in_iv = np.where(j >= 0, pos < iv[np.clip(j, 0, len(iv) - 1), 1], False)
            arr[idx] = in_iv
    ctx[te_hit] = "TE"
    ctx[~te_hit & gene_hit] = "Gene"
    return ctx


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{args.label}] Loading DANPOS positions ...", flush=True)
    df = load_danpos_positions(args.positions)
    if df.empty:
        sys.exit("No nucleosomes parsed from positions.xls")
    df = add_nrl(df)
    print(f"[{args.label}] {len(df):,} nucleosomes", flush=True)

    te_dict = te_intervals_by_chrom(load_te_bed(args.te_bed))
    gene_dict = load_gene_intervals(args.gff)
    df["context"] = classify_context(df, te_dict, gene_dict)

    with gzip.open(out_dir / f"{args.label}_nucleosomes.tsv.gz", "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)

    # ---- statistics --------------------------------------------------
    summary_rows = []
    for metric in METRICS:
        groups = {c: df.loc[df["context"] == c, metric].dropna()
                  for c in CONTEXT_ORDER}
        for c in CONTEXT_ORDER:
            vals = groups[c]
            summary_rows.append({
                "metric": metric,
                "context": c,
                "n": len(vals),
                "median": vals.median() if len(vals) else np.nan,
                "mean": vals.mean() if len(vals) else np.nan,
                "sd": vals.std() if len(vals) else np.nan,
            })
        # TE vs Gene / TE vs Inter Mann-Whitney U
        for comp, other in (("Gene", "Gene"), ("Inter", "Inter")):
            x = groups["TE"].to_numpy(float)
            y = groups[other].to_numpy(float)
            if len(x) >= 5 and len(y) >= 5:
                _, p = mannwhitneyu(x, y, alternative="two-sided")
            else:
                p = np.nan
            summary_rows.append({
                "metric": metric,
                "context": f"TE_vs_{other}",
                "n": len(x) + len(y),
                "median": np.nan,
                "mean": np.nan,
                "sd": np.nan,
                "mannwhitney_p": p,
                "cliffs_delta": cliffs_delta(x, y),
            })

    summ = pd.DataFrame(summary_rows)
    # BH-FDR over all Mann-Whitney p-values
    mw = summ["mannwhitney_p"].notna()
    if mw.any():
        summ.loc[mw, "bh_fdr"] = bh_fdr(summ.loc[mw, "mannwhitney_p"].to_numpy(float))
    summ.to_csv(out_dir / f"{args.label}_summary.tsv", sep="\t", index=False)

    # ---- figures -----------------------------------------------------
    for metric, ylab in METRICS.items():
        data = [df.loc[df["context"] == c, metric].dropna()
                for c in CONTEXT_ORDER]
        fig, ax = plt.subplots(figsize=(5.5, 4))
        bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                        medianprops=dict(color="black", lw=1.2),
                        showfliers=False)
        for patch, c in zip(bp["boxes"], CONTEXT_ORDER):
            patch.set_facecolor(CONTEXT_COLORS[c])
            patch.set_alpha(0.8)
        ax.set_xticks(range(1, 4))
        ax.set_xticklabels(CONTEXT_ORDER, fontsize=9)
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_title(f"{args.label} — {metric}", fontsize=10)
        style_ax(ax)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{args.label}_{metric}_boxplot.pdf")
        plt.close(fig)

    print(f"[{args.label}] Done. Summary saved to {out_dir / (args.label + '_summary.tsv')}")


if __name__ == "__main__":
    main()
