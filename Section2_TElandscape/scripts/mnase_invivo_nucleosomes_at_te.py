#!/usr/bin/env python3
"""
mnase_invivo_nucleosomes_at_te.py
=================================
Extension of mnase6_nucleosomes_at_te.py: adds 3 truly in vivo datasets
(plus the original Rhimi59 in vivo MNase as positive control), keeping the
5 in vitro Nuc-seq species (Mondo et al. JGI) for comparison.

Datasets:
  IN VIVO  MU_Rhimi59       MNase-seq (Rhimi)         WIG = pre-existing
  IN VIVO  AS_Ncrass2       MNase-seq (N. crassa)     WIG = new pipeline
  IN VIVO  AS_Bcine1        MAINE-seq (B. cinerea)    WIG = new pipeline
  IN VIVO  BA_Mereb1_atac   ATAC-seq (M. eburnea)     WIG = new pipeline (sign inverted)
  IN VITRO AS_Usnflo1       Nuc-seq                    WIG = Mondo et al.
  IN VITRO BA_Mereb1        Nuc-seq                    WIG = Mondo et al.
  IN VITRO BA_Mrafri1       Nuc-seq                    WIG = Mondo et al.
  IN VITRO BA_Rhoto         Nuc-seq                    WIG = Mondo et al.
  IN VITRO BL_Catan2        Nuc-seq                    WIG = Mondo et al.

For ATAC the signal is open-chromatin (high = accessible). To keep the same
"compaction" semantic (high = compacted) we invert the per-window mean
(occ_compact = max - mean_open). CV/AC features are kept as-is because they
measure rugosity/phasing regardless of sign.

Outputs:
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/mnase_invivo_te/
    per_species/<SP>_te_summary.tsv
    per_species/<SP>_te_features.pdf
    mnase_invivo_te_summary_all.tsv
    mnase_invivo_te_features_grid.pdf
    invivo_only_grid.pdf       (only the 4 truly in vivo, side-by-side)
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

# (species_key, wig_path, assay_label, in_vivo_bool, invert_sign)
# species_key is the EarlGrey key (used for FAI/BED) -- for the ATAC of Mereb1
# we re-use the same TE annotation but tag with a different display key.
SPECIES = [
    # truly in vivo
    ("MU_Rhimi59",     "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig",
        "MNase",  True,  False, "MU_Rhimi59"),
    ("AS_Ncrass2",     "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.wig",
        "MNase",  True,  False, "AS_Ncrass2"),
    ("AS_Bcine1",      "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.wig",
        "MAINE",  True,  False, "AS_Bcine1"),
    ("BA_Mereb1",      "/data/carlos/EDF-TE/TE_analysis/nucleosomes/atac-seq/BA_Mereb1/BA_Mereb1_atac.wig",
        "ATAC",   True,  True,  "BA_Mereb1_atac"),
    # in vitro Nuc-seq (Mondo)
    ("AS_Usnflo1",     "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/AS_Usnflo1/AS_Usnflo1_danpos/pooled/AS_Usnflo1_dedup.smooth.wig",
        "NucSeq", False, False, "AS_Usnflo1"),
    ("BA_Mereb1",      "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Mereb1/BA_Mereb1_danpos/pooled/BA_Mereb1_dedup.smooth.wig",
        "NucSeq", False, False, "BA_Mereb1"),
    ("BA_Mrafri1",     "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Mrafri1/BA_Mrafri1_danpos/pooled/BA_Mrafri1_dedup.smooth.wig",
        "NucSeq", False, False, "BA_Mrafri1"),
    ("BA_Rhoto",       "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Rhoto/BA_Rhoto_danpos/pooled/BA_Rhoto_dedup.smooth.wig",
        "NucSeq", False, False, "BA_Rhoto"),
    ("BL_Catan2",      "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BL_Catan2/BL_Catan2_danpos/pooled/BL_Catan2_dedup.smooth.wig",
        "NucSeq", False, False, "BL_Catan2"),
]
FAI_FMT = "/legserv/Temp/Carlos/EDF-TE/data/genomes/{sp}.fasta.fai"
TE_FMT  = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"

OUT  = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/mnase_invivo_te"
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


def build_table(sp_key, wig, invert):
    fai = FAI_FMT.format(sp=sp_key)
    bed = TE_FMT.format(sp=sp_key)
    print(f"[{sp_key}] sizes/wig/te (invert={invert})")
    sizes = chrom_sizes_from_fai(fai)
    arrs  = load_nuclpred_wig(wig, sizes)
    if invert:
        # global max over the genome to flip sign while keeping non-negative
        gmax = max((float(np.nanmax(v)) if v.size else 0.0) for v in arrs.values())
        arrs = {c: (gmax - v) for c, v in arrs.items()}
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


def per_species_summary(disp, df):
    g = df.groupby("category")[["mean_occ","cv","ac_val","ac_lag"]]
    s = g.agg(["mean","median"])
    s["n"] = df["category"].value_counts().reindex(s.index)
    s.to_csv(os.path.join(OUT, "per_species", f"{disp}_te_summary.tsv"), sep="\t")
    return s


def plot_species(disp, df, assay, in_vivo, out_pdf):
    ORDER  = ["TE-", "TE+"]
    COLORS = {"TE-": "#7F7F7F", "TE+": "#1F77B4" if in_vivo else "#FF7F0E"}
    feats = [
        ("mean_occ", "Mean occupancy" + (" (1-open)" if assay=="ATAC" else ""), True),
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
    tag = "in vivo" if in_vivo else "in vitro"
    fig.suptitle(f"{disp} — {assay} ({tag}) — TE vs background (500 bp windows)",
                 fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.93])
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf",".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    all_sum = {}
    species_df = {}
    meta = {}    # disp -> (assay, in_vivo)
    for sp_key, wig, assay, in_vivo, invert, disp in SPECIES:
        if not os.path.exists(wig):
            print(f"[SKIP] {disp}: missing {wig}"); continue
        df = build_table(sp_key, wig, invert)
        if df is None or df.empty:
            print(f"[SKIP] {disp}: no informative windows"); continue
        species_df[disp] = df
        meta[disp] = (assay, in_vivo)
        s = per_species_summary(disp, df)
        plot_species(disp, df, assay, in_vivo,
                     os.path.join(OUT, "per_species", f"{disp}_te_features.pdf"))
        all_sum[disp] = s

    # combined summary
    rows = []
    for disp, s in all_sum.items():
        assay, in_vivo = meta[disp]
        try:
            n_tp = s.loc["TE+","n"]
            n_tn = s.loc["TE-","n"]
            rows.append({
                "species": disp,
                "assay": assay,
                "in_vivo": in_vivo,
                "n_TEpos": int(n_tp.iloc[0]) if hasattr(n_tp,"iloc") else int(n_tp),
                "n_TEneg": int(n_tn.iloc[0]) if hasattr(n_tn,"iloc") else int(n_tn),
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
            print(f"[WARN] summary {disp}: {ex}")
    summ_df = pd.DataFrame(rows)
    summ_df.to_csv(os.path.join(OUT, "mnase_invivo_te_summary_all.tsv"),
                   sep="\t", index=False)
    print("\n=== Combined summary (TE+ vs TE-) ===")
    print(summ_df.to_string(index=False))

    # Order: in vivo first, then in vitro
    order = sorted(species_df.keys(),
                   key=lambda d: (0 if meta[d][1] else 1, d))

    feats = [("mean_occ","Mean occupancy", True),
             ("cv","CV (rugosity)", False),
             ("ac_val","AC peak (phasing)", False),
             ("ac_lag","AC lag (NRL, bp)", False)]
    ORDER  = ["TE-", "TE+"]

    def grid(sps, fname, title):
        fig, axes = plt.subplots(len(feats), len(sps),
                                 figsize=(2.4*len(sps), 2.8*len(feats)),
                                 squeeze=False)
        for j, disp in enumerate(sps):
            df = species_df[disp]
            assay, in_vivo = meta[disp]
            COLORS = {"TE-": "#7F7F7F", "TE+": "#1F77B4" if in_vivo else "#FF7F0E"}
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
                if i == 0:
                    tag = "in vivo" if in_vivo else "in vitro"
                    ax.set_title(f"{disp}\n{assay} ({tag})", fontsize=9,
                                 fontweight="bold" if in_vivo else "normal")
                if j == 0: ax.set_ylabel(ttl, fontsize=9)
        fig.suptitle(title, fontsize=12)
        fig.tight_layout(rect=[0,0,1,0.96])
        fig.savefig(os.path.join(OUT, fname), bbox_inches="tight")
        fig.savefig(os.path.join(OUT, fname.replace(".pdf",".png")),
                    dpi=180, bbox_inches="tight")
        plt.close(fig)

    grid(order, "mnase_invivo_te_features_grid.pdf",
         "Nucleosome features at TE vs background — in vivo (blue) vs in vitro (orange)")

    invivo_only = [d for d in order if meta[d][1]]
    if invivo_only:
        grid(invivo_only, "invivo_only_grid.pdf",
             f"In vivo chromatin assays — TE vs background ({len(invivo_only)} species)")

    print(f"\n[DONE] outputs in {OUT}")


if __name__ == "__main__":
    main()
