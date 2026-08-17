#!/usr/bin/env python3
"""
06_te_fuzziness.py  —  Nucleosome fuzziness at TE loci.

Analyses the fuzziness (positional dispersion) of DANPOS-called
nucleosomes at TE loci compared to genes and intergenic regions, and — when
an H3K9me3 broad-peak BED is provided — tests whether TE-associated
nucleosomes within H3K9me3 domains are more fuzzy than TE nucleosomes
outside H3K9me3 domains (consistent with the heterochromatin-proxy state).

Metrics:
  * fuzziness_score — positional SD (bp); lower = better positioned
  * NRL among TE nucleosomes (100–400 bp filter)

Statistics:
  * Kruskal–Wallis across TE / Gene / Inter contexts, followed by pairwise
    Mann–Whitney U (TE vs Gene, TE vs Inter) with BH-FDR correction and
    Cliff's delta.
  * Mann–Whitney U: TE∩H3K9me3 vs TE∩no-H3K9me3 fuzziness.

Inputs:
  --positions   DANPOS positions.xls
  --te_bed      EarlGrey filteredRepeats.bed
  --gff         gene annotation GFF/GFF3
  --h3k9_bed    (optional) H3K9me3 broad-peaks BED
  --label       sample label
  --out_dir     output directory

Outputs:
  <label>_fuzziness_nucleosomes.tsv.gz
  <label>_fuzziness_summary.tsv      — group medians + tests
  figures/<label>_fuzziness_contexts.pdf
  figures/<label>_fuzziness_h3k9.pdf
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--positions", required=True)
    p.add_argument("--te_bed", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--h3k9_bed", default=None)
    p.add_argument("--label", required=True)
    p.add_argument("--out_dir", default="fuzziness_out")
    return p.parse_args()


def load_bed_intervals(path):
    """Load BED intervals -> {chrom: Nx2 int32 array}."""
    d = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(("#", "track", "browser")):
                continue
            f = ln.split("\t")
            if len(f) < 3:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            d.setdefault(f[0], []).append((s, e))
    return {c: np.array(sorted(v), dtype=np.int32) for c, v in d.items()}


def inside(ivs, chrom, pos):
    if ivs is None or chrom not in ivs or len(ivs[chrom]) == 0:
        return False
    arr = ivs[chrom]
    j = np.searchsorted(arr[:, 0], pos, side="right") - 1
    return j >= 0 and pos < arr[j, 1]


def classify_context(df, te_dict, gene_dict):
    ctx = np.array(["Inter"] * len(df), dtype=object)
    for chrom, grp in df.groupby("chr"):
        pos = grp["smt_pos"].to_numpy()
        idx = grp.index
        tv = te_dict.get(chrom)
        gv = gene_dict.get(chrom)
        for label, iv, arr in (("te", tv, None), ("gene", gv, None)):
            if iv is None or len(iv) == 0:
                continue
            j = np.searchsorted(iv[:, 0], pos, side="right") - 1
            hit = np.where(j >= 0, pos < iv[np.clip(j, 0, len(iv) - 1), 1], False)
            if label == "te":
                ctx[idx[hit]] = "TE"
            else:
                ctx[idx[hit & (ctx[idx] != "TE")]] = "Gene"
    return ctx


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = load_danpos_positions(args.positions)
    if df.empty:
        sys.exit("No nucleosomes parsed from positions.xls")
    df = add_nrl(df)

    te_dict = te_intervals_by_chrom(load_te_bed(args.te_bed))
    gene_dict = load_gene_intervals(args.gff)
    df["context"] = classify_context(df, te_dict, gene_dict)

    # H3K9me3 overlap for TE nucleosomes
    df["in_h3k9"] = False
    if args.h3k9_bed:
        h3k9 = load_bed_intervals(args.h3k9_bed)
        te_mask = df["context"] == "TE"
        df.loc[te_mask, "in_h3k9"] = [
            inside(h3k9, chrom, pos)
            for chrom, pos in df.loc[te_mask, ["chr", "smt_pos"]].itertuples(index=False)]

    with gzip.open(out_dir / f"{args.label}_fuzziness_nucleosomes.tsv.gz", "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)

    # ---- Kruskal-Wallis + pairwise MWU across contexts ----------------
    groups = {c: df.loc[df["context"] == c, "fuzziness_score"].dropna()
              for c in ("TE", "Gene", "Inter")}
    groups = {c: v for c, v in groups.items() if len(v) >= 5}

    rows = []
    if len(groups) >= 3:
        _, kp = kruskal(*[groups[c].to_numpy(float) for c in groups])
        rows.append({"test": "kruskal_wallis_contexts", "p": kp})

    for other in ("Gene", "Inter"):
        if other not in groups or "TE" not in groups:
            continue
        x = groups["TE"].to_numpy(float)
        y = groups[other].to_numpy(float)
        _, p = mannwhitneyu(x, y, alternative="two-sided")
        rows.append({
            "test": f"MWU_TE_vs_{other}",
            "p": p,
            "median_TE": groups["TE"].median(),
            f"median_{other}": groups[other].median(),
            "cliffs_delta": cliffs_delta(x, y),
        })

    # ---- H3K9me3 within-TE comparison ---------------------------------
    if args.h3k9_bed and (df["in_h3k9"]).any():
        x = df.loc[(df["context"] == "TE") & df["in_h3k9"], "fuzziness_score"].dropna()
        y = df.loc[(df["context"] == "TE") & ~df["in_h3k9"], "fuzziness_score"].dropna()
        if len(x) >= 5 and len(y) >= 5:
            _, p = mannwhitneyu(x, y, alternative="two-sided")
            rows.append({
                "test": "MWU_TE_h3k9_vs_TE_no_h3k9",
                "p": p,
                "median_TE_in_h3k9": x.median(),
                "median_TE_out_h3k9": y.median(),
                "cliffs_delta": cliffs_delta(x, y),
            })

    summ = pd.DataFrame(rows)
    mw = summ["test"].str.startswith("MWU")
    if mw.any():
        summ.loc[mw, "bh_fdr"] = bh_fdr(summ.loc[mw, "p"].to_numpy(float))
    summ.to_csv(out_dir / f"{args.label}_fuzziness_summary.tsv", sep="\t", index=False)

    # ---- figures -------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.5, 4))
    data = [groups[c] for c in groups]
    labels = list(groups)
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False,
                    medianprops=dict(color="black", lw=1.2))
    colors = {"TE": "#d62728", "Gene": "#1f77b4", "Inter": "#999999"}
    for patch, lab in zip(bp["boxes"], labels):
        patch.set_facecolor(colors.get(lab, "#999"))
        patch.set_alpha(0.8)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Fuzziness score (bp)", fontsize=9)
    ax.set_title(f"{args.label} — nucleosome fuzziness", fontsize=10)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.label}_fuzziness_contexts.pdf")
    plt.close(fig)

    if args.h3k9_bed and (df["in_h3k9"]).any():
        fig, ax = plt.subplots(figsize=(5, 4))
        data = [df.loc[(df["context"] == "TE") & df["in_h3k9"], "fuzziness_score"].dropna(),
                df.loc[(df["context"] == "TE") & ~df["in_h3k9"], "fuzziness_score"].dropna()]
        bp = ax.boxplot(data, patch_artist=True, widths=0.5, showfliers=False,
                        medianprops=dict(color="black", lw=1.2))
        for patch, col in zip(bp["boxes"], ["#8B0000", "#e41a1c"]):
            patch.set_facecolor(col)
            patch.set_alpha(0.8)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["TE ∩ H3K9me3", "TE ∩ no-H3K9me3"], fontsize=9)
        ax.set_ylabel("Fuzziness score (bp)", fontsize=9)
        ax.set_title(f"{args.label} — TE fuzziness vs H3K9me3", fontsize=10)
        style_ax(ax)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{args.label}_fuzziness_h3k9.pdf")
        plt.close(fig)

    print(f"[{args.label}] Done.")


if __name__ == "__main__":
    main()
