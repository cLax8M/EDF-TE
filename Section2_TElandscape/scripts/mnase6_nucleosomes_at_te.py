#!/usr/bin/env python3
"""
mnase6_nucleosomes_at_te.py
===========================
Same descriptive analysis as rhimi59_nucleosomes_at_h3k9_te.py but for the
6 species with REAL DANPOS MNase-seq data, comparing TE+ vs TE- windows
(no H3K9me3 because only Rhizopus has it).

Per species, 500-bp windows are split into:
  TE+   window has >=50% overlap with EarlGrey TE annotation
  TE-   background (no TE)

Features (raw, not rank-normalised, no model):
  mean_occ   nucleosome occupancy
  cv         sd / mean (rugosity)
  ac_val     autocorr peak in 140-240 bp (phasing strength)
  ac_lag     position of that peak (NRL)

Outputs in
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/mnase6_te/
    per_species/<SP>_te_summary.tsv
    per_species/<SP>_te_features.pdf
    mnase6_te_summary_all.tsv         (all 6 + Rhimi59 reference)
    mnase6_te_features_grid.pdf       (4x6 boxplots, grid)
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics

WIN = 500
MIN_OVERLAP = 0.5

# 6 species with real MNase data
SPECIES = {
    "MU_Rhimi59":  "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig",
    "AS_Usnflo1":  "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/AS_Usnflo1/AS_Usnflo1_danpos/pooled/AS_Usnflo1_dedup.smooth.wig",
    "BA_Mereb1":   "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Mereb1/BA_Mereb1_danpos/pooled/BA_Mereb1_dedup.smooth.wig",
    "BA_Mrafri1":  "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Mrafri1/BA_Mrafri1_danpos/pooled/BA_Mrafri1_dedup.smooth.wig",
    "BA_Rhoto":    "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Rhoto/BA_Rhoto_danpos/pooled/BA_Rhoto_dedup.smooth.wig",
    "BL_Catan2":   "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BL_Catan2/BL_Catan2_danpos/pooled/BL_Catan2_dedup.smooth.wig",
}
FAI_FMT = "/legserv/Temp/Carlos/EDF-TE/data/genomes/{sp}.fasta.fai"
TE_FMT  = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"

OUT  = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/mnase6_te"
os.makedirs(os.path.join(OUT, "per_species"), exist_ok=True)


def chrom_sizes_from_fai(p):
    fai = pd.read_csv(p, sep="\t", header=None,
                      names=["chrom","L","off","lb","lby"])
    return dict(zip(fai["chrom"], fai["L"].astype(int)))


def overlap_array(bed_path, sizes, win):
    bp = {c: np.zeros(L, dtype=np.uint8) for c, L in sizes.items()}
    if not os.path.exists(bed_path):
        return {c: np.zeros(L // win, dtype=bool) for c, L in sizes.items()}
    with open(bed_path) as fh:
        for ln in fh:
            f = ln.split("\t")
            if len(f) < 3: continue
            c = f[0]
            if c not in bp: continue
            try:
                s = int(f[1]); e = int(f[2])
            except ValueError: continue
            e = min(e, len(bp[c])); s = max(0, s)
            if s < e: bp[c][s:e] = 1
    out = {}
    for c, L in sizes.items():
        n = L // win
        if n == 0: out[c] = np.zeros(0, dtype=bool); continue
        a = bp[c][:n*win].reshape(n, win)
        out[c] = (a.mean(axis=1) >= MIN_OVERLAP)
    return out


def build_table(sp, wig):
    fai = FAI_FMT.format(sp=sp)
    bed = TE_FMT.format(sp=sp)
    print(f"[{sp}] sizes/wig/te")
    sizes = chrom_sizes_from_fai(fai)
    arrs  = load_nuclpred_wig(wig, sizes)
    rows  = []
    for c, occ in arrs.items():
        n = len(occ) // WIN
        if n == 0: continue
        means, cvs, aclags, acvals, _ = window_metrics(occ, WIN)
        rows.append(pd.DataFrame({
            "chrom": c,
            "start": np.arange(n) * WIN,
            "end":   np.arange(n) * WIN + WIN,
            "mean_occ": means, "cv": cvs,
            "ac_val": acvals, "ac_lag": aclags,
        }))
    if not rows: return None
    df = pd.concat(rows, ignore_index=True)
    te_ov = overlap_array(bed, sizes, WIN)
    te_flag = np.zeros(len(df), dtype=bool)
    i0 = 0
    for c, sub in df.groupby("chrom", sort=False):
        n = len(sub)
        te_flag[i0:i0+n] = te_ov.get(c, np.zeros(n, dtype=bool))[:n]
        i0 += n
    df["te"] = te_flag
    df = df[df["mean_occ"] > 1e-3].copy()
    df["category"] = np.where(df["te"], "TE+", "TE-")
    return df


def per_species_summary(sp, df):
    g = df.groupby("category")[["mean_occ","cv","ac_val","ac_lag"]]
    s = g.agg(["mean","median"])
    s["n"] = df["category"].value_counts().reindex(s.index)
    s.to_csv(os.path.join(OUT, "per_species", f"{sp}_te_summary.tsv"), sep="\t")
    return s


def plot_species(sp, df, out_pdf):
    ORDER  = ["TE-", "TE+"]
    COLORS = {"TE-": "#7F7F7F", "TE+": "#1F77B4"}
    feats = [
        ("mean_occ", "Mean occupancy", True),
        ("cv",       "CV (rugosity)",  False),
        ("ac_val",   "AC peak (phasing)", False),
        ("ac_lag",   "AC peak position (bp, NRL)", False),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
    for ax, (col, ttl, log) in zip(axes, feats):
        data = [df.loc[df["category"]==c, col].values for c in ORDER]
        bp = ax.boxplot(data, positions=np.arange(len(ORDER)), widths=0.65,
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=1.5))
        for patch, c in zip(bp["boxes"], ORDER):
            patch.set_facecolor(COLORS[c]); patch.set_edgecolor("black")
        ax.set_xticks(np.arange(len(ORDER)))
        ax.set_xticklabels([f"{c}\n(n={(df['category']==c).sum()})" for c in ORDER],
                            fontsize=9)
        ax.set_title(ttl, fontsize=10)
        if log: ax.set_yscale("log")
        if col == "ac_lag": ax.set_ylim(140, 240)
    fig.suptitle(f"{sp} — nucleosome features at TE vs background (real MNase, 500 bp windows)",
                 fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.93])
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf",".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    all_sum = {}
    species_df = {}
    for sp, wig in SPECIES.items():
        if not os.path.exists(wig):
            print(f"[SKIP] {sp}: missing {wig}"); continue
        df = build_table(sp, wig)
        if df is None or df.empty:
            print(f"[SKIP] {sp}: no informative windows"); continue
        species_df[sp] = df
        s = per_species_summary(sp, df)
        plot_species(sp, df, os.path.join(OUT, "per_species", f"{sp}_te_features.pdf"))
        all_sum[sp] = s

    # combined summary table: per species, mean of mean_occ and cv per category, plus TE+ vs TE- ratio/diff
    rows = []
    for sp, s in all_sum.items():
        try:
            rows.append({
                "species": sp,
                "n_TEpos": int(s.loc["TE+","n"].iloc[0]) if hasattr(s.loc["TE+","n"], "iloc") else int(s.loc["TE+","n"]),
                "n_TEneg": int(s.loc["TE-","n"].iloc[0]) if hasattr(s.loc["TE-","n"], "iloc") else int(s.loc["TE-","n"]),
                "mean_occ_TEpos":  s.loc["TE+",("mean_occ","mean")],
                "mean_occ_TEneg":  s.loc["TE-",("mean_occ","mean")],
                "delta_mean_occ":  s.loc["TE+",("mean_occ","mean")] - s.loc["TE-",("mean_occ","mean")],
                "cv_TEpos":        s.loc["TE+",("cv","mean")],
                "cv_TEneg":        s.loc["TE-",("cv","mean")],
                "delta_cv":        s.loc["TE+",("cv","mean")] - s.loc["TE-",("cv","mean")],
                "ac_val_TEpos":    s.loc["TE+",("ac_val","mean")],
                "ac_val_TEneg":    s.loc["TE-",("ac_val","mean")],
                "delta_ac_val":    s.loc["TE+",("ac_val","mean")] - s.loc["TE-",("ac_val","mean")],
                "ac_lag_TEpos":    s.loc["TE+",("ac_lag","mean")],
                "ac_lag_TEneg":    s.loc["TE-",("ac_lag","mean")],
            })
        except Exception as ex:
            print(f"[WARN] summary {sp}: {ex}")
    summ_df = pd.DataFrame(rows)
    summ_df.to_csv(os.path.join(OUT, "mnase6_te_summary_all.tsv"), sep="\t", index=False)
    print("\n=== Combined summary (TE+ vs TE-) ===")
    print(summ_df.to_string(index=False))

    # 4x6 grid
    feats = [("mean_occ","Mean occupancy", True),
             ("cv","CV (rugosity)", False),
             ("ac_val","AC peak (phasing)", False),
             ("ac_lag","AC lag (NRL, bp)", False)]
    sps = list(species_df.keys())
    fig, axes = plt.subplots(len(feats), len(sps),
                             figsize=(2.4*len(sps), 2.8*len(feats)),
                             squeeze=False)
    ORDER  = ["TE-", "TE+"]
    COLORS = {"TE-": "#7F7F7F", "TE+": "#1F77B4"}
    for j, sp in enumerate(sps):
        df = species_df[sp]
        for i, (col, ttl, log) in enumerate(feats):
            ax = axes[i][j]
            data = [df.loc[df["category"]==c, col].values for c in ORDER]
            bp = ax.boxplot(data, positions=[0,1], widths=0.65,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.2))
            for patch, c in zip(bp["boxes"], ORDER):
                patch.set_facecolor(COLORS[c]); patch.set_edgecolor("black")
            ax.set_xticks([0,1]); ax.set_xticklabels(ORDER, fontsize=8)
            if log: ax.set_yscale("log")
            if col == "ac_lag": ax.set_ylim(140, 240)
            if i == 0: ax.set_title(sp, fontsize=10, fontweight="bold")
            if j == 0: ax.set_ylabel(ttl, fontsize=9)
    fig.suptitle("Real MNase: nucleosome features at TE vs background (6 species)",
                 fontsize=12)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(os.path.join(OUT, "mnase6_te_features_grid.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "mnase6_te_features_grid.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[DONE] outputs in {OUT}")


if __name__ == "__main__":
    main()
