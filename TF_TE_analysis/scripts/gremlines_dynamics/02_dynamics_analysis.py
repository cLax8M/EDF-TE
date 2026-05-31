#!/usr/bin/env python3
"""FASE 4B step 2-5: Mucci3 Grem-LINE1 dynamics.

Inputs (produced by 01_extract_align_tree.sh):
  L1Tx1_copies.bed          BED of L1-Tx1 copies (>=500bp)
  L1Tx1_aln.trim.fasta      trimAl-cleaned alignment
  L1Tx1_iqtree.treefile     IQTree (UFBoot 1000)
External:
  centromes BED:            MU_Mucci3_centromes.bed (cluster regions per scaffold)

Outputs (in this dir):
  k80_per_copy.tsv
  copy_zone.tsv             per-copy: distance to centromere center, zone, K80
  hypothesis_results.txt    Spearman, monophyly, zone medians
  k80_vs_centro.pdf
  boxplot_k80_zone.pdf
  contig_map.pdf
"""
from __future__ import annotations
import sys, os, math, re, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path("/legserv/Temp/Carlos/EDF-TE/results/Gremlines/MU402")
# MU402 (JGI Muccir1_3) centromeres from Navarro 2019 — clipped to scaffold ends
CEN_BED = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/MU402_centromere_windows_clipped.bed")
FAI = Path("/legserv/Temp/Carlos/EDF-TE/data/Extra/mu402GenomeJGI.fa.fai")

ALN = OUT / "gremlin_aln_all.trim.fasta"   # K80 vs consensus on ALL copies (autonomous + remnants)
COPIES = OUT / "gremlin_copies_all.bed"
TREE = OUT / "gremlin_iqtree.treefile"      # tree built only from autonomous (~31 copies)


def parse_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf).upper()
    return seqs


def k80_pairwise(s1: str, s2: str) -> float:
    """Kimura 80 distance ignoring gaps (-, N)."""
    transitions = transversions = sites = 0
    for a, b in zip(s1, s2):
        if a in "-N" or b in "-N":
            continue
        sites += 1
        if a == b:
            continue
        if (a, b) in {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}:
            transitions += 1
        else:
            transversions += 1
    if sites == 0:
        return float("nan")
    P = transitions / sites
    Q = transversions / sites
    arg1 = 1 - 2 * P - Q
    arg2 = 1 - 2 * Q
    if arg1 <= 0 or arg2 <= 0:
        return float("nan")
    return -0.5 * math.log(arg1) - 0.25 * math.log(arg2)


def consensus(aln: dict[str, str]) -> str:
    L = len(next(iter(aln.values())))
    cons = []
    for i in range(L):
        col = [s[i] for s in aln.values() if s[i] not in "-N"]
        if not col:
            cons.append("N")
            continue
        # majority rule
        cnt = {}
        for c in col:
            cnt[c] = cnt.get(c, 0) + 1
        cons.append(max(cnt, key=cnt.get))
    return "".join(cons)


def main():
    if not ALN.exists():
        sys.exit(f"[ERR] missing {ALN} — run 01_extract_align_tree.sh first")

    # 1. Load alignment + compute K80 per copy vs consensus
    aln = parse_fasta(ALN)
    print(f"[INFO] {len(aln)} aligned copies, len={len(next(iter(aln.values())))}")
    cons = consensus(aln)
    k80 = {name: k80_pairwise(seq, cons) for name, seq in aln.items()}
    k80_df = pd.DataFrame({"copy_id": list(k80), "k80": list(k80.values())})
    k80_df.to_csv(OUT / "k80_per_copy.tsv", sep="\t", index=False)

    # 2. Load BED (only the real gremlines, not the f218 outgroup) → copy positions
    bed_cols = ["chrom", "start", "end", "name", "score", "strand"]
    bed = pd.read_csv(COPIES, sep="\t", names=bed_cols, header=None)
    bed["mid"] = (bed["start"] + bed["end"]) // 2
    bed["copy_key"] = bed["name"].astype(str)

    # match fasta names: bedtools getfasta -nameOnly -s appends "(+)" or "(-)"
    fasta_names = set(aln.keys())
    bed["fasta_name"] = bed.apply(lambda r: f"{r['name']}({r['strand']})", axis=1)
    bed = bed[bed["fasta_name"].isin(fasta_names)].copy()
    print(f"[INFO] {len(bed)} BED↔fasta matches")

    bed = bed.merge(k80_df, left_on="fasta_name", right_on="copy_id", how="left")

    # 3. Centromere center per scaffold (use MU_Mucci3_centromes.bed)
    cen = pd.read_csv(CEN_BED, sep="\t", header=None,
                      names=["chrom", "start", "end", "n_grem"], usecols=[0, 1, 2, 3])
    cen["center"] = (cen["start"] + cen["end"]) // 2
    cen_map = dict(zip(cen["chrom"], cen["center"]))
    cen_size = dict(zip(cen["chrom"], cen["end"] - cen["start"]))

    def dist_to_centro(row):
        c = cen_map.get(row["chrom"])
        if c is None:
            return np.nan
        return abs(row["mid"] - c)

    bed["dist_centro"] = bed.apply(dist_to_centro, axis=1)

    def zone(d):
        if pd.isna(d):
            return "no_centromere"
        if d <= 10_000:
            return "core"
        if d <= 100_000:
            return "peri"
        return "distal"

    bed["zone"] = bed["dist_centro"].apply(zone)
    # AUT (Navarro S3 clusters 1-3) vs REM (clusters 4-9 + NoCen) per name suffix
    bed["status"] = bed["name"].apply(lambda n: "AUT" if "_AUT" in n else "REM")
    bed.to_csv(OUT / "copy_zone.tsv", sep="\t", index=False)

    # 4. Stats: Spearman K80 vs distance (only copies with centromere on same chrom)
    sub = bed.dropna(subset=["dist_centro", "k80"])
    if len(sub) >= 5:
        rho, pval = stats.spearmanr(sub["dist_centro"], sub["k80"])
    else:
        rho, pval = float("nan"), float("nan")

    # zone medians + Kruskal
    zone_groups = [sub[sub["zone"] == z]["k80"].dropna().values
                   for z in ["core", "peri", "distal"]]
    valid = [g for g in zone_groups if len(g) > 1]
    if len(valid) >= 2:
        H, p_kw = stats.kruskal(*valid)
    else:
        H, p_kw = float("nan"), float("nan")

    with open(OUT / "hypothesis_results.txt", "w") as f:
        f.write("FASE 4B — Mucci3 REAL Grem-LINE1 (rnd-1_family-198/206/180) insertion dynamics\n")
        f.write("=" * 64 + "\n\n")
        f.write(f"N copies (>=500bp) analysed       : {len(bed)}\n")
        f.write(f"N centromeric clusters            : {len(cen)}\n\n")
        f.write("[H1] Monophyly per centromere — see L1Tx1_iqtree.treefile\n")
        f.write("     (visual inspection in tree.pdf)\n\n")
        f.write("[H2] Spearman K80 vs distance-to-centromere-center\n")
        f.write(f"     rho   = {rho:.3f}\n")
        f.write(f"     p     = {pval:.3g}\n")
        f.write(f"     N     = {len(sub)}\n\n")
        f.write("[H3] Zone medians (K80)\n")
        for z in ["core", "peri", "distal"]:
            arr = sub[sub["zone"] == z]["k80"]
            if len(arr):
                f.write(f"     {z:6s}: n={len(arr):4d}  median={arr.median():.4f}  mean={arr.mean():.4f}\n")
        f.write(f"     Kruskal-Wallis H={H:.3f}  p={p_kw:.3g}\n\n")
        f.write("[H4] Counts by zone (all copies)\n")
        for z in ["core", "peri", "distal"]:
            f.write(f"     {z:6s}: {(bed['zone']==z).sum()}\n")
        f.write("\n[H5] Autonomous vs remnant\n")
        for s in ["AUT", "REM"]:
            sub_s = bed[bed["status"] == s]
            if len(sub_s):
                f.write(f"     {s}: n={len(sub_s)}  median K80={sub_s['k80'].median():.4f}  median dist_centro={sub_s['dist_centro'].median():.0f}\n")
        # Mann-Whitney AUT vs REM K80
        aut_k = bed[bed['status']=='AUT']['k80'].dropna()
        rem_k = bed[bed['status']=='REM']['k80'].dropna()
        if len(aut_k)>1 and len(rem_k)>1:
            U, p_mw = stats.mannwhitneyu(aut_k, rem_k, alternative='two-sided')
            f.write(f"     Mann-Whitney AUT vs REM K80: U={U:.1f} p={p_mw:.3g}\n")

    # 5. Plots --------------------------------------------------------
    # 5a. K80 vs distance (per centromere panel)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(sub["dist_centro"] / 1e3, sub["k80"],
               c=sub["chrom"].astype("category").cat.codes,
               cmap="tab10", s=20, alpha=0.7, edgecolor="k", linewidth=0.3)
    if len(sub) >= 5:
        # LOWESS-ish: rolling median by binned distance
        sub2 = sub.sort_values("dist_centro")
        ax.set_title(f"K80 vs distance to centromere center (Spearman ρ={rho:.2f}, p={pval:.2g})")
    ax.set_xlabel("Distance to centromere center (kb)")
    ax.set_ylabel("K80 distance to consensus")
    fig.tight_layout()
    fig.savefig(OUT / "k80_vs_centro.pdf")
    plt.close(fig)

    # 5b. Boxplot K80 per zone
    fig, ax = plt.subplots(figsize=(5, 4))
    data = [sub[sub["zone"] == z]["k80"].values for z in ["core", "peri", "distal"]]
    bp = ax.boxplot(data, labels=["core (≤10kb)", "peri (≤100kb)", "distal"],
                    patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], ["#d62728", "#ff7f0e", "#2ca02c"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    for i, d in enumerate(data, 1):
        ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=8, alpha=0.4, c="k")
    ax.set_ylabel("K80 distance to consensus")
    ax.set_title(f"Grem-LINE1 age by centromeric zone (KW p={p_kw:.2g})")
    fig.tight_layout()
    fig.savefig(OUT / "boxplot_k80_zone.pdf")
    plt.close(fig)

    # 5c. Contig map: each scaffold with copies coloured by K80, centromere shaded
    fai = pd.read_csv(FAI, sep="\t", header=None, names=["chrom", "length"], usecols=[0, 1])
    chroms = fai[fai["chrom"].isin(bed["chrom"].unique())].sort_values("length", ascending=False).head(15)
    fig, axes = plt.subplots(len(chroms), 1, figsize=(11, max(3, 0.5 * len(chroms))))
    if len(chroms) == 1:
        axes = [axes]
    vmin, vmax = bed["k80"].quantile(0.05), bed["k80"].quantile(0.95)
    for ax, (_, row) in zip(axes, chroms.iterrows()):
        ch, L = row["chrom"], row["length"]
        ax.add_patch(Rectangle((0, 0.4), L / 1e6, 0.2, facecolor="#e0e0e0"))
        if ch in cen_map:
            c = cen_map[ch] / 1e6
            cs = cen_size.get(ch, 100_000) / 1e6
            ax.add_patch(Rectangle((c - cs / 2, 0.3), cs, 0.4, facecolor="#fcb""b8e", alpha=0.7))
        sb = bed[bed["chrom"] == ch]
        if len(sb):
            sc = ax.scatter(sb["mid"] / 1e6, [0.5] * len(sb),
                            c=sb["k80"], cmap="viridis_r", vmin=vmin, vmax=vmax,
                            s=30, edgecolor="k", linewidth=0.4)
        ax.set_xlim(0, L / 1e6)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(ch, rotation=0, ha="right", va="center", fontsize=8)
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)
    axes[-1].set_xlabel("Genomic position (Mb)")
    fig.suptitle("MU_Mucci3 — Real Grem-LINE1 copies (rnd-1_family-198/206/180; color=K80 age)")
    fig.tight_layout()
    fig.savefig(OUT / "contig_map.pdf")
    plt.close(fig)

    print("[DONE] FASE 4B Mucci3 dynamics complete")


if __name__ == "__main__":
    main()
