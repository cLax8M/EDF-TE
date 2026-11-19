#!/usr/bin/env python3
"""
01_nuclpred_analysis.py  —  Genome-wide NuclPred occupancy analysis.

For each species with a NuclPred per-base predicted-occupancy WIG
(``Section2_TElandscape/nucleosomes/nuclpred/{SP}_nuclpred.wig``):

  A. Tiles the genome into 500 bp windows and computes, per window,
     mean occupancy and autocorrelation peak (140–240 bp lag) — the two
     metrics used to define the four chromatin states (Open / Phased /
     Fuzzy_mid / Fuzzy_high; Fuzzy_high ≈ heterochromatin proxy).
  B. Classifies every window into one of the four chromatin states using
     the genome-wide q33/q67 occupancy thresholds and the median phasing.
  C. Reports TE-vs-background occupancy: for each TE superfamily the mean
     predicted occupancy over TE-overlapping windows vs non-TE, non-genic
     (background) windows, together with the Mann–Whitney U test and
     Cliff's delta effect size.

Inputs (CLI):
  --wig          NuclPred WIG (fixedStep/variableStep, gzip allowed)
  --fai          genome .fasta.fai (chromosome lengths)
  --te_bed       EarlGrey filteredRepeats.bed
  --gff          (optional) gene annotation GFF/GFF3 to mask genic windows
  --species      species code (output labelling)
  --win          window size (default 500)
  --out_dir      output directory

Outputs:
  <species>_window_metrics.tsv.gz   — per-window metrics + state
  <species>_chromatin_states.tsv    — per-state summary
  <species>_te_bias.tsv             — per-superfamily TE vs background
  figures/<species>_states_summary.pdf
  figures/<species>_te_bias.pdf
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
    load_wig,
    load_chrom_sizes,
    load_te_bed,
    window_metrics,
    cliffs_delta,
    style_ax,
)

STATE_COLORS = {
    "Open":       "#2ca02c",
    "Phased":     "#1f77b4",
    "Fuzzy_mid":  "#ffbb78",
    "Fuzzy_high": "#d62728",
}
STATE_ORDER = ["Open", "Phased", "Fuzzy_mid", "Fuzzy_high"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--wig",      required=True, help="NuclPred WIG")
    p.add_argument("--fai",      required=True, help="genome .fasta.fai")
    p.add_argument("--te_bed",   required=True, help="EarlGrey filteredRepeats.bed")
    p.add_argument("--gff",      default=None, help="gene GFF/GFF3 (optional)")
    p.add_argument("--species",  required=True)
    p.add_argument("--win",      type=int, default=500)
    p.add_argument("--out_dir",  default="nuclpred_out")
    return p.parse_args()


def classify_states(occ: np.ndarray, ac_val: np.ndarray):
    """Assign the 4 chromatin states from occupancy + phasing thresholds."""
    valid = np.isfinite(occ)
    occ_q33, occ_q67 = np.nanpercentile(occ[valid], [33, 67])
    ac_med = np.nanmedian(ac_val[valid])
    state = np.full(occ.shape, "", dtype=object)
    for i in range(len(occ)):
        if not valid[i]:
            state[i] = "NA"
            continue
        if occ[i] < occ_q33:
            state[i] = "Open"
        elif ac_val[i] >= ac_med:
            state[i] = "Phased"
        elif occ[i] < occ_q67:
            state[i] = "Fuzzy_mid"
        else:
            state[i] = "Fuzzy_high"
    return state, occ_q33, occ_q67, ac_med


def window_bed(metrics: dict, win: int, n_chrom_windows: int,
               chrom_names, chrom_sizes):
    """Build per-window coordinate table for overlap computations."""
    rows = []
    for c in chrom_names:
        n = n_chrom_windows[c]
        for i in range(n):
            rows.append((c, i * win, min((i + 1) * win, chrom_sizes[c])))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    df["mean_occ"] = metrics["mean_occ"]
    df["cv"] = metrics["cv"]
    df["ac_lag"] = metrics["ac_lag"]
    df["ac_val"] = metrics["ac_val"]
    df["state"] = metrics["state"]
    return df


def overlap_frac(starts, ends, ivs):
    """Fraction of each interval covered by sorted intervals `ivs` (Nx2)."""
    out = np.zeros(len(starts), dtype=np.float32)
    if ivs is None or len(ivs) == 0:
        return out
    s0 = ivs[:, 0]
    e0 = ivs[:, 1]
    for i in range(len(starts)):
        s, e = starts[i], ends[i]
        if e <= s:
            continue
        lo = np.searchsorted(s0, s, side="right") - 1
        hi = np.searchsorted(e0, e, side="left")
        lo = max(lo, 0)
        hi = min(hi, len(s0) - 1)
        if lo > hi:
            # still check adjacency
            idx = np.searchsorted(s0, s, side="right")
            cand = []
            if idx > 0:
                cand.append(idx - 1)
            if idx < len(s0):
                cand.append(idx)
            for j in cand:
                if s0[j] < e and e0[j] > s:
                    ov = min(e, e0[j]) - max(s, s0[j])
                    out[i] += ov
        else:
            seg = np.clip(np.minimum(e0[lo:hi + 1], e) -
                          np.maximum(s0[lo:hi + 1], s), 0, None)
            out[i] = seg.sum()
        out[i] = min(out[i], e - s) / (e - s)
    return out


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    chrom_sizes = load_chrom_sizes(args.fai)
    print(f"[{args.species}] Loading WIG ...", flush=True)
    occ = load_wig(args.wig, chrom_sizes)
    te = load_te_bed(args.te_bed)

    # --- tiled windows per chromosome --------------------------------
    win = args.win
    w_rows = []
    m_rows = []
    for chrom, size in chrom_sizes.items():
        arr = occ.get(chrom)
        if arr is None or size < win:
            continue
        m = window_metrics(arr, win)
        n = len(m["mean_occ"])
        if n == 0:
            continue
        for i in range(n):
            w_rows.append((chrom, i * win, min((i + 1) * win, size)))
            m_rows.append((m["mean_occ"][i], m["cv"][i],
                           m["ac_lag"][i], m["ac_val"][i], bool(m["valid"][i])))

    wdf = pd.DataFrame(w_rows, columns=["chrom", "start", "end"])
    mdf = pd.DataFrame(m_rows, columns=["mean_occ", "cv", "ac_lag", "ac_val", "valid"])
    wdf = pd.concat([wdf, mdf], axis=1)
    wdf = wdf[wdf["valid"]].copy()
    print(f"[{args.species}] {len(wdf):,} valid {win}-bp windows", flush=True)

    # --- chromatin-state classification ------------------------------
    state, q33, q67, ac_med = classify_states(
        wdf["mean_occ"].to_numpy(float), wdf["ac_val"].to_numpy(float))
    wdf["state"] = state

    state_summary = (wdf.groupby("state").agg(
        n_windows=("mean_occ", "size"),
        mean_occ=("mean_occ", "mean"),
        median_ac=("ac_val", "median")).reset_index())
    state_summary["frac"] = state_summary["n_windows"] / len(wdf)
    state_summary.to_csv(out_dir / f"{args.species}_chromatin_states.tsv",
                         sep="\t", index=False)

    with gzip.open(out_dir / f"{args.species}_window_metrics.tsv.gz", "wt") as fh:
        wdf.to_csv(fh, sep="\t", index=False)
    print(f"[{args.species}] states: " +
          ", ".join(f"{s}:{state_summary.loc[state_summary.state==s,'n_windows'].iloc[0] if (state_summary.state==s).any() else 0}"
                    for s in STATE_ORDER), flush=True)

    # --- TE vs background occupancy per superfamily ------------------
    te_by_chrom = {c: g[["start", "end"]].to_numpy(np.int32)
                   for c, g in te.groupby("chrom")}
    gene_by_chrom = None
    if args.gff:
        from nucleosome_common import load_gene_intervals
        gene_by_chrom = load_gene_intervals(args.gff)

    te_frac = np.zeros(len(wdf), dtype=np.float32)
    gene_frac = np.zeros(len(wdf), dtype=np.float32)
    for chrom, grp in wdf.groupby("chrom"):
        idx = grp.index
        te_frac[idx] = overlap_frac(grp["start"].to_numpy(),
                                    grp["end"].to_numpy(),
                                    te_by_chrom.get(chrom))
        if gene_by_chrom:
            gene_frac[idx] = overlap_frac(grp["start"].to_numpy(),
                                          grp["end"].to_numpy(),
                                          gene_by_chrom.get(chrom))
    wdf["te_frac"] = te_frac
    wdf["gene_frac"] = gene_frac

    # per-TE-copy occupancy (associate each insertion with its window)
    te["win_mean_occ"] = np.nan
    for chrom, grp in te.groupby("chrom"):
        w = wdf[wdf["chrom"] == chrom]
        if w.empty:
            continue
        for _, row in grp.iterrows():
            mid = (row["start"] + row["end"]) // 2
            hit = w[(w["start"] <= mid) & (w["end"] > mid)]
            if not hit.empty:
                te.loc[row.name, "win_mean_occ"] = hit.iloc[0]["mean_occ"]
    te_occ = te[te["win_mean_occ"].notna()].copy()
    te_occ = te_occ[te_occ["length"] >= 50]

    bg = wdf[(wdf["te_frac"] < 0.2) & (wdf["gene_frac"] < 0.3)]
    bg_occ = bg["mean_occ"].to_numpy(float)

    rows = []
    for cls, grp in te_occ.groupby("te_class"):
        vals = grp["win_mean_occ"].to_numpy(float)
        if len(vals) < 5 or len(bg_occ) < 5:
            continue
        try:
            u, p = mannwhitneyu(vals, bg_occ, alternative="two-sided")
        except ValueError:
            u, p = np.nan, np.nan
        rows.append({
            "species": args.species,
            "te_class": cls,
            "n_insertions": len(vals),
            "mean_te_occupancy": vals.mean(),
            "median_te_occupancy": np.median(vals),
            "background_occupancy": bg_occ.mean(),
            "delta_occupancy": vals.mean() - bg_occ.mean(),
            "mannwhitney_p": p,
            "cliffs_delta": cliffs_delta(vals, bg_occ),
        })
    bias = pd.DataFrame(rows)
    if len(bias):
        bias.to_csv(out_dir / f"{args.species}_te_bias.tsv", sep="\t", index=False)
        print(f"[{args.species}] TE bias computed for {len(bias)} classes", flush=True)

    # --- figures -----------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 3.5))
    counts = [state_summary.loc[state_summary.state == s, "n_windows"].sum()
              for s in STATE_ORDER]
    ax.bar(STATE_ORDER, counts, color=[STATE_COLORS[s] for s in STATE_ORDER])
    ax.set_ylabel("windows", fontsize=9)
    ax.set_title(f"{args.species} — chromatin states", fontsize=10)
    style_ax(ax)
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.species}_states_summary.pdf")
    plt.close(fig)

    if len(bias):
        bias_sorted = bias.sort_values("delta_occupancy")
        fig, ax = plt.subplots(figsize=(6, 0.45 * len(bias_sorted) + 1))
        cols = [("#d62728" if d > 0 else "#1f77b4")
                for d in bias_sorted["delta_occupancy"]]
        ax.barh(bias_sorted["te_class"], bias_sorted["delta_occupancy"], color=cols)
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_xlabel("Δ occupancy (TE − background)", fontsize=9)
        ax.set_title(f"{args.species} — TE nucleosome bias", fontsize=10)
        style_ax(ax)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{args.species}_te_bias.pdf")
        plt.close(fig)

    print(f"[{args.species}] Done.")


if __name__ == "__main__":
    main()
