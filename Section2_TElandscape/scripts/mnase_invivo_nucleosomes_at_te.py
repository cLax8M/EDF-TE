#!/usr/bin/env python3
"""
mnase_invivo_nucleosomes_at_te.py
=================================
Same descriptive analysis as mnase6_nucleosomes_at_te.py but extended with
3 NEW in-vivo datasets to test the "TE = compact chromatin" hypothesis with
proper native chromatin (the previous 5 Mondo et al. JGI Nuc-seq tracks are
in vitro reconstituted nucleosomes and consistently show TE+ < TE-, opposite
to the in vivo Rhimi59 MNase signature).

Datasets (9):
  In vivo (target panel)
    MU_Rhimi59      MNase    Rhizopus microsporus  (legacy, in vivo)
    AS_Ncrass2      MNase    Neurospora crassa OR74A   SRR13067301-303
    AS_Bcine1       MAINE    Botrytis cinerea B05.10   SRR13224644/47/48
    BA_Mereb1_atac  ATAC     Microbotryum eburnea      SRR28069093/94/106
  In vitro reference (Mondo et al., JGI Nuc-seq)
    AS_Usnflo1      NucSeq
    BA_Mereb1       NucSeq
    BA_Mrafri1      NucSeq
    BA_Rhoto        NucSeq
    BL_Catan2       NucSeq

Per species, 500-bp windows split TE+ (>=50% EarlGrey overlap) vs TE-.
Features: mean_occ, cv, ac_val (140-240 bp peak), ac_lag (NRL).

ATAC note: signal is accessibility (high = open). For comparability with the
nucleosome occupancy panels, the raw AC/occupancy metrics are still reported
"as is" (interpretation is opposite). For the supervised compaction Δ scores
downstream, ATAC species are flagged so the sign can be inverted.

Outputs in
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/mnase_invivo_te/
    per_species/<SP>_te_summary.tsv
    per_species/<SP>_te_features.{pdf,png}
    mnase_invivo_te_summary_all.tsv     (all 9, with assay + in_vivo flags)
    mnase_invivo_te_features_grid.{pdf,png}    (4 features x 9 species)
    invivo_only_grid.{pdf,png}                 (only the 4 in-vivo tracks)
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

# Per-species (label, wig path, assay, in_vivo flag, source species ID for TE/FAI)
MV = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes"
RHIMI_WIG = "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"

PANEL = [
    # in vivo
    dict(label="MU_Rhimi59",     sp="MU_Rhimi59",
         wig=RHIMI_WIG,                                             assay="MNase",  in_vivo=True),
    dict(label="AS_Ncrass2",     sp="AS_Ncrass2",
         wig=f"{MV}/nuc-seq/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.wig",
         assay="MNase",  in_vivo=True),
    dict(label="AS_Bcine1",      sp="AS_Bcine1",
         wig=f"{MV}/nuc-seq/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.wig",
         assay="MAINE",  in_vivo=True),
    dict(label="BA_Mereb1_atac", sp="BA_Mereb1",
         wig=f"{MV}/atac-seq/BA_Mereb1/BA_Mereb1_atac.wig",
         assay="ATAC",   in_vivo=True),
    # in vitro Nuc-seq (Mondo et al., JGI)
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

FAI_FMT  = "/legserv/Temp/Carlos/EDF-TE/data/genomes/{sp}.fasta.fai"
TE_DATA  = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
TE_MV    = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"

OUT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/mnase_invivo_te"
os.makedirs(os.path.join(OUT, "per_species"), exist_ok=True)


def te_bed(sp):
    a = TE_DATA.format(sp=sp); b = TE_MV.format(sp=sp)
    return a if os.path.exists(a) else b


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
    bed = te_bed(sp)
    print(f"[{sp}] sizes/wig/te ({os.path.basename(wig)})")
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


def per_species_summary(label, df):
    g = df.groupby("category")[["mean_occ","cv","ac_val","ac_lag"]]
    s = g.agg(["mean","median"])
    s["n"] = df["category"].value_counts().reindex(s.index)
    s.to_csv(os.path.join(OUT, "per_species", f"{label}_te_summary.tsv"), sep="\t")
    return s


def plot_species(label, df, assay, in_vivo, out_pdf):
    ORDER  = ["TE-", "TE+"]
    COLORS = {"TE-": "#7F7F7F", "TE+": "#1F77B4"}
    feats = [
        ("mean_occ", "Mean signal" + (" (accessibility)" if assay=="ATAC" else " (occupancy)"), True),
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
    suffix = f" — {assay} ({'in vivo' if in_vivo else 'in vitro'})"
    fig.suptitle(f"{label}{suffix} — features at TE vs background (500 bp windows)",
                 fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.93])
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf",".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def grid_plot(species_meta, species_df, fname_base, suptitle):
    feats = [("mean_occ","Mean signal", True),
             ("cv","CV (rugosity)", False),
             ("ac_val","AC peak (phasing)", False),
             ("ac_lag","AC lag (NRL, bp)", False)]
    sps = list(species_df.keys())
    fig, axes = plt.subplots(len(feats), len(sps),
                             figsize=(2.4*len(sps), 2.8*len(feats)),
                             squeeze=False)
    ORDER  = ["TE-", "TE+"]
    COLORS = {"TE-": "#7F7F7F", "TE+": "#1F77B4"}
    for j, lab in enumerate(sps):
        df = species_df[lab]; meta = species_meta[lab]
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
                tag = f"{lab}\n{meta['assay']} ({'in vivo' if meta['in_vivo'] else 'in vitro'})"
                ax.set_title(tag, fontsize=9, fontweight="bold")
            if j == 0: ax.set_ylabel(ttl, fontsize=9)
    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(os.path.join(OUT, f"{fname_base}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{fname_base}.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    all_sum = {}
    species_df = {}
    species_meta = {}
    for entry in PANEL:
        label = entry["label"]; wig = entry["wig"]; sp = entry["sp"]
        if not os.path.exists(wig):
            print(f"[SKIP] {label}: missing {wig}"); continue
        df = build_table(sp, wig)
        if df is None or df.empty:
            print(f"[SKIP] {label}: no informative windows"); continue
        species_df[label] = df
        species_meta[label] = entry
        s = per_species_summary(label, df)
        plot_species(label, df, entry["assay"], entry["in_vivo"],
                     os.path.join(OUT, "per_species", f"{label}_te_features.pdf"))
        all_sum[label] = s

    # combined summary
    rows = []
    for label, s in all_sum.items():
        meta = species_meta[label]
        try:
            def gv(idx, col, agg):
                v = s.loc[idx, (col, agg)]
                return float(v.iloc[0]) if hasattr(v, "iloc") else float(v)
            def gn(idx):
                v = s.loc[idx, "n"]
                return int(v.iloc[0]) if hasattr(v, "iloc") else int(v)
            rows.append({
                "species": label,
                "assay":   meta["assay"],
                "in_vivo": meta["in_vivo"],
                "n_TEpos": gn("TE+"),
                "n_TEneg": gn("TE-"),
                "mean_occ_TEpos": gv("TE+","mean_occ","mean"),
                "mean_occ_TEneg": gv("TE-","mean_occ","mean"),
                "delta_mean_occ": gv("TE+","mean_occ","mean") - gv("TE-","mean_occ","mean"),
                "cv_TEpos":       gv("TE+","cv","mean"),
                "cv_TEneg":       gv("TE-","cv","mean"),
                "delta_cv":       gv("TE+","cv","mean") - gv("TE-","cv","mean"),
                "ac_val_TEpos":   gv("TE+","ac_val","mean"),
                "ac_val_TEneg":   gv("TE-","ac_val","mean"),
                "delta_ac_val":   gv("TE+","ac_val","mean") - gv("TE-","ac_val","mean"),
                "ac_lag_TEpos":   gv("TE+","ac_lag","mean"),
                "ac_lag_TEneg":   gv("TE-","ac_lag","mean"),
            })
        except Exception as ex:
            print(f"[WARN] summary {label}: {ex}")
    summ_df = pd.DataFrame(rows)
    # in vivo first
    summ_df = summ_df.sort_values(["in_vivo","species"], ascending=[False, True])
    summ_df.to_csv(os.path.join(OUT, "mnase_invivo_te_summary_all.tsv"), sep="\t", index=False)
    print("\n=== Combined summary (TE+ vs TE-) ===")
    print(summ_df.to_string(index=False))

    # Full grid (9 species)
    ordered = list(summ_df["species"])
    sd = {k: species_df[k] for k in ordered if k in species_df}
    sm = {k: species_meta[k] for k in ordered if k in species_meta}
    grid_plot(sm, sd, "mnase_invivo_te_features_grid",
              "Nucleosome / chromatin features at TE vs background — in vivo (4) + in vitro Nuc-seq (5)")

    # in-vivo only grid
    sd2 = {k: v for k, v in sd.items() if sm[k]["in_vivo"]}
    sm2 = {k: v for k, v in sm.items() if sm[k]["in_vivo"]}
    if sd2:
        grid_plot(sm2, sd2, "invivo_only_grid",
                  "In vivo native chromatin only — TE vs background")
    print(f"\n[DONE] outputs in {OUT}")


if __name__ == "__main__":
    main()
