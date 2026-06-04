#!/usr/bin/env python3
"""
te_rich_vs_poor_chromatin.py
============================
Compare nucleosome / chromatin features in TE-RICH vs TE-POOR regions across
the 9 species with real chromatin data (4 in vivo + 5 in vitro Nuc-seq).

Why a separate analysis from `mnase_invivo_nucleosomes_at_te.py`?
  - That script splits 500-bp windows by direct TE overlap (TE+ vs TE-).
  - This script aggregates 10-kb windows and ranks them by TE density,
    contrasting TE-rich (top quartile TE%) vs TE-poor (bottom quartile)
    *neighbourhoods*. This captures large heterochromatic domains where TEs
    cluster, more comparable to the "TE compartment" in published studies.

Per species:
  1. Tile genome in 10 kb windows.
  2. Compute TE% per window (bp covered by EarlGrey TE / 10 000).
  3. Aggregate the WIG occupancy (or accessibility) per window: mean signal,
     CV, autocorrelation peak (140-240 bp) value/lag.
  4. Quartiles by TE%: Q1 (TE-poor) vs Q4 (TE-rich); also Q2/Q3 saved.
  5. Mann-Whitney U Q4 vs Q1 for each feature; effect size = Cliff's delta.

Outputs in
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/te_rich_vs_poor/
    per_species/<LABEL>_windows.tsv.gz   (one row per 10 kb window)
    per_species/<LABEL>_quartiles.tsv    (mean/median per quartile)
    per_species/<LABEL>_features.{pdf,png}
    te_rich_vs_poor_summary.tsv          (per species, Q4 vs Q1 stats)
    te_rich_vs_poor_grid.{pdf,png}       (4 features x 9 species, Q1..Q4)
    invivo_only_grid.{pdf,png}
"""
import os, sys, gzip
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics

WIN = 10_000

MV = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes"
RHIMI_WIG = "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"

PANEL = [
    dict(label="MU_Rhimi59",     sp="MU_Rhimi59",
         wig=RHIMI_WIG, assay="MNase",  in_vivo=True),
    dict(label="AS_Ncrass2",     sp="AS_Ncrass2",
         wig=f"{MV}/nuc-seq/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.wig",
         assay="MNase",  in_vivo=True),
    dict(label="AS_Bcine1",      sp="AS_Bcine1",
         wig=f"{MV}/nuc-seq/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.wig",
         assay="MAINE",  in_vivo=True),
    dict(label="BA_Mereb1_atac", sp="BA_Mereb1",
         wig=f"{MV}/atac-seq/BA_Mereb1/BA_Mereb1_atac.wig",
         assay="ATAC",   in_vivo=True),
    dict(label="AS_Usnflo1",     sp="AS_Usnflo1",
         wig=f"{MV}/nuc-seq/AS_Usnflo1/AS_Usnflo1_danpos/pooled/AS_Usnflo1_dedup.smooth.wig",
         assay="NucSeq", in_vivo=False),
    dict(label="BA_Mereb1",      sp="BA_Mereb1",
         wig=f"{MV}/nuc-seq/BA_Mereb1/BA_Mereb1_danpos/pooled/BA_Mereb1_dedup.smooth.wig",
         assay="NucSeq", in_vivo=False),
    dict(label="BA_Mrafri1",     sp="BA_Mrafri1",
         wig=f"{MV}/nuc-seq/BA_Mrafri1/BA_Mrafri1_danpos/pooled/BA_Mrafri1_dedup.smooth.wig",
         assay="NucSeq", in_vivo=False),
    dict(label="BA_Rhoto",       sp="BA_Rhoto",
         wig=f"{MV}/nuc-seq/BA_Rhoto/BA_Rhoto_danpos/pooled/BA_Rhoto_dedup.smooth.wig",
         assay="NucSeq", in_vivo=False),
    dict(label="BL_Catan2",      sp="BL_Catan2",
         wig=f"{MV}/nuc-seq/BL_Catan2/BL_Catan2_danpos/pooled/BL_Catan2_dedup.smooth.wig",
         assay="NucSeq", in_vivo=False),
]

FAI_FMT = "/legserv/Temp/Carlos/EDF-TE/data/genomes/{sp}.fasta.fai"
TE_DATA = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
TE_MV   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"

OUT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/te_rich_vs_poor"
os.makedirs(os.path.join(OUT, "per_species"), exist_ok=True)


def te_bed(sp):
    a = TE_DATA.format(sp=sp); b = TE_MV.format(sp=sp)
    return a if os.path.exists(a) else b


def chrom_sizes_from_fai(p):
    fai = pd.read_csv(p, sep="\t", header=None,
                      names=["chrom","L","off","lb","lby"])
    return dict(zip(fai["chrom"], fai["L"].astype(int)))


def te_pct_per_window(bed_path, sizes, win):
    """Return dict chrom -> array of TE bp coverage per window (0..win)."""
    bp = {c: np.zeros(L, dtype=np.uint8) for c, L in sizes.items()}
    if os.path.exists(bed_path):
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
        if n == 0:
            out[c] = np.zeros(0, dtype=np.float32); continue
        a = bp[c][:n*win].reshape(n, win)
        out[c] = a.sum(axis=1).astype(np.float32) / win  # 0..1
    return out


def cliffs_delta(x, y):
    """Non-parametric effect size, robust to non-normality."""
    x = np.asarray(x); y = np.asarray(y)
    if len(x) == 0 or len(y) == 0: return np.nan
    # subsample for tractability
    nmax = 5000
    if len(x) > nmax:
        x = np.random.RandomState(0).choice(x, nmax, replace=False)
    if len(y) > nmax:
        y = np.random.RandomState(1).choice(y, nmax, replace=False)
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return float(gt - lt) / (len(x) * len(y))


def build_table(label, sp, wig):
    fai = FAI_FMT.format(sp=sp); bed = te_bed(sp)
    print(f"[{label}] wig={os.path.basename(wig)}  bed={os.path.basename(bed)}")
    sizes = chrom_sizes_from_fai(fai)
    arrs  = load_nuclpred_wig(wig, sizes)
    te_pct = te_pct_per_window(bed, sizes, WIN)
    rows = []
    for c, occ in arrs.items():
        n = len(occ) // WIN
        if n == 0: continue
        means, cvs, aclags, acvals, _ = window_metrics(occ, WIN)
        tp = te_pct.get(c, np.zeros(n, dtype=np.float32))[:n]
        rows.append(pd.DataFrame({
            "chrom": c,
            "start": np.arange(n) * WIN,
            "end":   np.arange(n) * WIN + WIN,
            "te_pct":  tp,
            "mean_sig": means, "cv": cvs,
            "ac_val": acvals, "ac_lag": aclags,
        }))
    if not rows: return None
    df = pd.concat(rows, ignore_index=True)
    # drop empty/no-signal windows
    df = df[df["mean_sig"] > 1e-3].copy()
    if df.empty: return df
    # quartiles by TE%
    q = pd.qcut(df["te_pct"].rank(method="first"), 4, labels=["Q1","Q2","Q3","Q4"])
    df["quartile"] = q.astype(str)
    return df


def quartile_summary(label, df):
    g = df.groupby("quartile")[["te_pct","mean_sig","cv","ac_val","ac_lag"]]
    s = g.agg(["mean","median","std","count"])
    s.to_csv(os.path.join(OUT, "per_species", f"{label}_quartiles.tsv"), sep="\t")
    return s


def stats_q4_vs_q1(df):
    out = {}
    for col in ["mean_sig","cv","ac_val","ac_lag"]:
        x = df.loc[df["quartile"]=="Q4", col].dropna().values
        y = df.loc[df["quartile"]=="Q1", col].dropna().values
        if len(x) < 3 or len(y) < 3:
            out[f"{col}_p"] = np.nan; out[f"{col}_d"] = np.nan; continue
        try:
            _, p = mannwhitneyu(x, y, alternative="two-sided")
        except Exception:
            p = np.nan
        out[f"{col}_p"] = float(p)
        out[f"{col}_d"] = cliffs_delta(x, y)
        out[f"{col}_Q4mean"] = float(np.mean(x))
        out[f"{col}_Q1mean"] = float(np.mean(y))
        out[f"{col}_delta_Q4_Q1"] = float(np.mean(x) - np.mean(y))
    return out


def plot_species(label, df, assay, in_vivo, out_pdf):
    ORDER = ["Q1","Q2","Q3","Q4"]
    cmap = plt.cm.viridis
    COLORS = {q: cmap(i/3) for i, q in enumerate(ORDER)}
    feats = [
        ("mean_sig", "Mean signal" + (" (acc.)" if assay=="ATAC" else " (occ.)"), True),
        ("cv",       "CV (rugosity)",  False),
        ("ac_val",   "AC peak (phasing)", False),
        ("ac_lag",   "AC lag (NRL, bp)", False),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
    for ax, (col, ttl, log) in zip(axes, feats):
        data = [df.loc[df["quartile"]==q, col].values for q in ORDER]
        bp = ax.boxplot(data, positions=np.arange(len(ORDER)), widths=0.65,
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=1.5))
        for patch, q in zip(bp["boxes"], ORDER):
            patch.set_facecolor(COLORS[q]); patch.set_edgecolor("black")
        ax.set_xticks(np.arange(len(ORDER)))
        ax.set_xticklabels([f"{q}\n(n={(df['quartile']==q).sum()})" for q in ORDER],
                            fontsize=9)
        ax.set_title(ttl, fontsize=10)
        if log: ax.set_yscale("log")
        if col == "ac_lag": ax.set_ylim(140, 240)
    suffix = f" — {assay} ({'in vivo' if in_vivo else 'in vitro'})"
    fig.suptitle(f"{label}{suffix} — chromatin features by TE-density quartile (10 kb windows)",
                 fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.93])
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf",".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def grid_plot(species_meta, species_df, fname_base, suptitle):
    feats = [("mean_sig","Mean signal", True),
             ("cv","CV (rugosity)", False),
             ("ac_val","AC peak (phasing)", False),
             ("ac_lag","AC lag (NRL, bp)", False)]
    sps = list(species_df.keys())
    ORDER = ["Q1","Q2","Q3","Q4"]
    cmap = plt.cm.viridis
    COLORS = {q: cmap(i/3) for i, q in enumerate(ORDER)}
    fig, axes = plt.subplots(len(feats), len(sps),
                             figsize=(2.6*len(sps), 2.8*len(feats)),
                             squeeze=False)
    for j, lab in enumerate(sps):
        df = species_df[lab]; meta = species_meta[lab]
        for i, (col, ttl, log) in enumerate(feats):
            ax = axes[i][j]
            data = [df.loc[df["quartile"]==q, col].values for q in ORDER]
            bp = ax.boxplot(data, positions=np.arange(4), widths=0.65,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.2))
            for patch, q in zip(bp["boxes"], ORDER):
                patch.set_facecolor(COLORS[q]); patch.set_edgecolor("black")
            ax.set_xticks(np.arange(4)); ax.set_xticklabels(ORDER, fontsize=8)
            if log: ax.set_yscale("log")
            if col == "ac_lag": ax.set_ylim(140, 240)
            if i == 0:
                tag = f"{lab}\n{meta['assay']} ({'in vivo' if meta['in_vivo'] else 'in vitro'})"
                ax.set_title(tag, fontsize=9, fontweight="bold")
            if j == 0: ax.set_ylabel(ttl, fontsize=9)
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(os.path.join(OUT, f"{fname_base}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{fname_base}.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    species_df = {}; species_meta = {}; rows = []
    for entry in PANEL:
        label = entry["label"]; sp = entry["sp"]; wig = entry["wig"]
        if not os.path.exists(wig):
            print(f"[SKIP] {label}: missing {wig}"); continue
        df = build_table(label, sp, wig)
        if df is None or df.empty:
            print(f"[SKIP] {label}: no usable windows"); continue
        # save windows table
        df.to_csv(os.path.join(OUT, "per_species", f"{label}_windows.tsv.gz"),
                  sep="\t", index=False, compression="gzip")
        quartile_summary(label, df)
        plot_species(label, df, entry["assay"], entry["in_vivo"],
                     os.path.join(OUT, "per_species", f"{label}_features.pdf"))
        st = stats_q4_vs_q1(df)
        st["species"]   = label
        st["assay"]     = entry["assay"]
        st["in_vivo"]   = entry["in_vivo"]
        st["te_pct_Q1"] = float(df.loc[df["quartile"]=="Q1","te_pct"].mean())
        st["te_pct_Q4"] = float(df.loc[df["quartile"]=="Q4","te_pct"].mean())
        st["n_windows"] = int(len(df))
        rows.append(st)
        species_df[label] = df; species_meta[label] = entry

    summ = pd.DataFrame(rows)
    cols_first = ["species","assay","in_vivo","n_windows","te_pct_Q1","te_pct_Q4"]
    rest = [c for c in summ.columns if c not in cols_first]
    summ = summ[cols_first + rest]
    summ = summ.sort_values(["in_vivo","species"], ascending=[False, True])
    summ.to_csv(os.path.join(OUT, "te_rich_vs_poor_summary.tsv"), sep="\t", index=False)
    print("\n=== Q4 vs Q1 summary (TE-rich vs TE-poor 10 kb windows) ===")
    show = summ[["species","assay","in_vivo","n_windows","te_pct_Q1","te_pct_Q4",
                 "mean_sig_delta_Q4_Q1","mean_sig_d","mean_sig_p",
                 "cv_delta_Q4_Q1","cv_d",
                 "ac_val_delta_Q4_Q1","ac_val_d",
                 "ac_lag_delta_Q4_Q1","ac_lag_d"]].copy()
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # ordered grids: in vivo first
    ordered = list(summ["species"])
    sd = {k: species_df[k] for k in ordered if k in species_df}
    sm = {k: species_meta[k] for k in ordered if k in species_meta}
    grid_plot(sm, sd, "te_rich_vs_poor_grid",
              "Chromatin features by TE-density quartile (10 kb) — 4 in vivo + 5 in vitro Nuc-seq")

    sd2 = {k: v for k, v in sd.items() if sm[k]["in_vivo"]}
    sm2 = {k: v for k, v in sm.items() if sm[k]["in_vivo"]}
    if sd2:
        grid_plot(sm2, sd2, "invivo_only_grid",
                  "In vivo native chromatin only — features by TE-density quartile (10 kb)")
    print(f"\n[DONE] outputs in {OUT}")


if __name__ == "__main__":
    main()
