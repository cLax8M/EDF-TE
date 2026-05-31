#!/usr/bin/env python3
"""
delta_score_simple2f.py
=======================
Panel A: per phylum, mean Δ score (TE − genome background) ± SE across species.
Panel B: per TE class, mean Δ score across species.

Uses cross_species_2state/<SP>_states.tsv.gz from simple_2feature_compaction.py.
score = occ_pct - cv_pct. Δ > 0 => TE in compacted regions.
"""
import os, glob, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT   = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/simple2f"
CS_DIR = os.path.join(ROOT, "cross_species_2state")
TE_DIR = "/data/carlos/EDF-TE/TE_analysis/annotated"
WIN = 500

PHYLUM = {
    "AS":"Ascomycota","BA":"Basidiomycota","BL":"Blastocladiomycota",
    "BS":"Basidiobolomycota","CH":"Chytridiomycota","CR":"Cryptomycota",
    "EN":"Entomophthoromycota","GL":"Glomeromycota","KX":"Kickxellomycota",
    "MI":"Microsporidia","MO":"Mortierellomycota","MU":"Mucoromycota",
    "NU":"Nuclear-other","OL":"Olpidiomycota","AP":"Aphelidiomycota",
    "ZO":"Zoopagomycota",
}

def te_class(s):
    s = s or ""
    if s.startswith("LTR"): return "LTR"
    if s.startswith("LINE"): return "LINE"
    if s.startswith("SINE"): return "SINE"
    if s.startswith("DNA"): return "DNA"
    if s.startswith("RC"): return "RC"
    if s.startswith("Satellite"): return "Satellite"
    if s.startswith("Simple"): return "Simple"
    if s.startswith("Low"): return "Low"
    if s == "rRNA": return "rRNA"
    if s == "tRNA": return "tRNA"
    if s == "Unknown": return "Unknown"
    if s == "Other": return "Other"
    return "Unknown"

CLASS_ORDER = ["LTR","LINE","SINE","DNA","RC","Satellite","Simple","Low","rRNA","tRNA","Unknown"]


def process(args):
    sp, ts, bed = args
    try:
        df = pd.read_csv(ts, sep="\t",
                         usecols=["chrom","start","end","mean_occ","score","state"])
        df = df[df["state"] != "NoSig"].copy()
        if df.empty: return sp, None, None
        bg = float(df["score"].mean())
        # build chrom -> array of score (NaN for uncovered windows)
        ca = {}
        for c, sub in df.groupby("chrom", sort=False):
            idx = (sub["start"].values // WIN).astype(np.int64)
            n = int(idx.max()) + 1
            arr = np.full(n, np.nan, dtype=np.float32)
            arr[idx] = sub["score"].values.astype(np.float32)
            ca[c] = arr
        if not os.path.exists(bed):
            return sp, bg, None
        te = pd.read_csv(bed, sep="\t", header=None,
                         names=["chrom","start","end","superfam","score","strand"],
                         usecols=["chrom","start","end","superfam"])
        if te.empty: return sp, bg, None
        te["te_class"] = te["superfam"].map(te_class)
        sums={}; counts={}
        for r in te.itertuples(index=False):
            arr = ca.get(r.chrom)
            if arr is None: continue
            i0 = r.start // WIN
            i1 = min((r.end - 1)//WIN, len(arr)-1)
            if i0 > i1: continue
            seg = arr[i0:i1+1]; mask = ~np.isnan(seg)
            if not mask.any(): continue
            sums[r.te_class]   = sums.get(r.te_class, 0.0) + float(seg[mask].sum())
            counts[r.te_class] = counts.get(r.te_class, 0)   + int(mask.sum())
        rows=[]
        for cls,v in sums.items():
            n=counts[cls]
            if n<5: continue
            rows.append(dict(species=sp, te_class=cls, mean_score_TE=v/n,
                             n_windows=n, bg_score=bg, delta=v/n - bg))
        return sp, bg, pd.DataFrame(rows)
    except Exception as ex:
        print(f"[WARN] {sp}: {ex}", file=sys.stderr)
        return sp, None, None


def main():
    files = sorted(glob.glob(os.path.join(CS_DIR, "*_states.tsv.gz")))
    jobs=[]
    for f in files:
        sp = os.path.basename(f).replace("_states.tsv.gz","")
        bed = os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                           f"{sp}.filteredRepeats.bed")
        jobs.append((sp, f, bed))
    print(f"[RUN] {len(jobs)} species")
    rows=[]; n=0
    with ProcessPoolExecutor(max_workers=16) as ex:
        futs=[ex.submit(process, j) for j in jobs]
        for fut in as_completed(futs):
            sp, bg, d = fut.result(); n+=1
            if d is not None and len(d): rows.append(d)
            if n%25==0: print(f"  [{n}/{len(jobs)}]")
    df = pd.concat(rows, ignore_index=True)
    df["phylum_code"] = df["species"].str.split("_", n=1).str[0]
    df["phylum"] = df["phylum_code"].map(PHYLUM).fillna(df["phylum_code"])
    df.to_csv(os.path.join(ROOT, "delta_score_TE_vs_background.tsv"), sep="\t", index=False)

    sp_pool = (df.groupby(["species","phylum"])
                 .apply(lambda d: pd.Series({
                     "delta_TE_all": np.average(d["delta"], weights=d["n_windows"])
                 })).reset_index())
    sp_pool.to_csv(os.path.join(ROOT, "delta_score_per_species.tsv"), sep="\t", index=False)

    phy = (sp_pool.groupby("phylum")
                  .agg(mean=("delta_TE_all","mean"),
                       sd  =("delta_TE_all","std"),
                       n   =("delta_TE_all","count"))
                  .reset_index())
    phy["se"] = phy["sd"] / np.sqrt(phy["n"].clip(lower=1))
    phy = phy[phy["n"]>=2].sort_values("mean")
    phy.to_csv(os.path.join(ROOT, "delta_score_by_phylum.tsv"), sep="\t", index=False)

    cls_per_sp = (df.groupby(["te_class","species"])
                    .apply(lambda d: np.average(d["delta"], weights=d["n_windows"]))
                    .reset_index().rename(columns={0:"delta"}))
    cls = (cls_per_sp.groupby("te_class")
                     .agg(mean=("delta","mean"),
                          sd  =("delta","std"),
                          n   =("delta","count"))
                     .reset_index())
    cls["se"] = cls["sd"]/np.sqrt(cls["n"].clip(lower=1))
    cls = cls[cls["te_class"].isin(CLASS_ORDER)].sort_values("mean")
    cls.to_csv(os.path.join(ROOT, "delta_score_by_TEclass.tsv"), sep="\t", index=False)

    # -------- plot --------
    fig, axes = plt.subplots(1, 2, figsize=(15, 0.32*max(len(phy),len(cls))+2),
                             gridspec_kw={"width_ratios":[1.05,1.0]})
    fig.patch.set_facecolor("#F5F8FB")

    ax = axes[0]; ax.set_facecolor("#F5F8FB")
    palette = plt.cm.Set2.colors + plt.cm.Set3.colors
    colors = [palette[i % len(palette)] for i in range(len(phy))]
    yy = np.arange(len(phy))
    ax.barh(yy, phy["mean"].values, xerr=phy["se"].values,
            color=colors, edgecolor="black", linewidth=0.4,
            error_kw=dict(ecolor="black", capsize=3, lw=0.8))
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels(phy["phylum"], fontsize=10)
    ax.set_xlabel("Δ score (TE − genome background)", fontsize=10)
    ax.set_title("A   Compaction signal at TE insertions vs Genome Background",
                 loc="left", fontsize=11, fontweight="bold")
    ax.text(0.98,0.04,
            "score = occ_pct − cv_pct\nΔ > 0 → TEs sit in more-compacted regions\nΔ < 0 → TEs sit in more-open regions",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#9CB4CC"))

    ax = axes[1]; ax.set_facecolor("#F5F8FB")
    color_map = {"LTR":"#D6604D","LINE":"#F4A582","SINE":"#B2182B",
                 "DNA":"#92C5DE","RC":"#1B7837"}
    bar_colors = [color_map.get(c,"#888888") for c in cls["te_class"]]
    yy = np.arange(len(cls))
    ax.barh(yy, cls["mean"].values, xerr=cls["se"].values,
            color=bar_colors, edgecolor="black", linewidth=0.4,
            error_kw=dict(ecolor="black", capsize=3, lw=0.8))
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels(cls["te_class"], fontsize=10)
    ax.set_xlabel("Δ score (TE − genome background)", fontsize=10)
    ax.set_title("B   Compaction Bias by TE Class (averaged across all species)",
                 loc="left", fontsize=11, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "delta_score_TE_vs_background.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(ROOT, "delta_score_TE_vs_background.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[DONE]")


if __name__ == "__main__":
    main()
