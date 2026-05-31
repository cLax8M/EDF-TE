#!/usr/bin/env python3
"""
delta_compaction_te_vs_background.py
====================================
For every species, compute per-window p_het (compaction probability from the
supervised model). Then for each TE annotation, compute the mean p_het over
windows overlapping the TE vs the mean p_het over the whole genome (informative
windows only, i.e. mean_occ > 1e-3).

Δ p_het (TE − background) per species per TE class.
Aggregate to:
  Panel A: per phylum (mean ± SE across species), TE-pooled (all classes).
  Panel B: per TE class, mean ± SE across species (all phyla pooled).
"""
import os, glob, sys, re, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT     = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin"
CS_DIR   = os.path.join(ROOT, "cross_species")
TE_DIR   = "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_ROOT = ROOT
WIN = 500

# Phylum mapping (prefix -> full name)
PHYLUM = {
    "AS": "Ascomycota",
    "BA": "Basidiomycota",
    "BL": "Blastocladiomycota",
    "BS": "Basidiobolomycota",
    "CH": "Chytridiomycota",
    "CR": "Cryptomycota",
    "EN": "Entomophthoromycota",
    "GL": "Glomeromycota",
    "KX": "Kickxellomycota",
    "MI": "Microsporidia",
    "MO": "Mortierellomycota",
    "MU": "Mucoromycota",
    "NU": "Nuclear-other",
    "OL": "Olpidiomycota",
    "AP": "Aphelidiomycota",
    "ZO": "Zoopagomycota",
}

# TE class mapping (collapse superfam -> high-level class for panel B)
def te_class(superfam: str) -> str:
    s = superfam or ""
    if s.startswith("LTR"):           return "LTR"
    if s.startswith("LINE"):          return "LINE"
    if s.startswith("SINE"):          return "SINE"
    if s.startswith("DNA"):           return "DNA"
    if s.startswith("RC"):            return "RC"
    if s.startswith("Satellite"):     return "Satellite"
    if s.startswith("Simple"):        return "Simple"
    if s.startswith("Low"):           return "Low"
    if s in ("rRNA",):                return "rRNA"
    if s in ("tRNA",):                return "tRNA"
    if s == "Unknown":                return "Unknown"
    if s == "Other":                  return "Other"
    return "Unknown"

CLASS_ORDER = ["LTR","LINE","SINE","DNA","RC","Satellite","Simple","Low",
               "rRNA","tRNA","Unknown"]


def process_species(args):
    sp, ts, bed = args
    try:
        df = pd.read_csv(ts, sep="\t",
                         usecols=["chrom","start","end","mean_occ","p_het"])
        df = df[df["mean_occ"].values > 1e-3].copy()
        if df.empty:
            return sp, None, None
        bg_mean = float(df["p_het"].mean())

        # window index per chrom
        # build dict: chrom -> (start_array_sorted by index, p_het_array)
        # because windows are uniform WIN, we can index directly:
        chrom_arrays = {}
        for c, sub in df.groupby("chrom", sort=False):
            # build a full array of p_het indexed by window number; gaps -> NaN
            idx = (sub["start"].values // WIN).astype(np.int64)
            n = int(idx.max()) + 1
            arr = np.full(n, np.nan, dtype=np.float32)
            arr[idx] = sub["p_het"].values.astype(np.float32)
            chrom_arrays[c] = arr

        if not os.path.exists(bed):
            return sp, bg_mean, None
        te = pd.read_csv(bed, sep="\t", header=None,
                         names=["chrom","start","end","superfam","score","strand"],
                         usecols=["chrom","start","end","superfam"])
        if te.empty:
            return sp, bg_mean, None
        te["te_class"] = te["superfam"].map(te_class)

        # accumulate sum/n by class using window-overlap weighting
        sums = {}; counts = {}
        for r in te.itertuples(index=False):
            arr = chrom_arrays.get(r.chrom)
            if arr is None: continue
            i0 = r.start // WIN
            i1 = (r.end - 1) // WIN
            i1 = min(i1, len(arr) - 1)
            if i0 > i1: continue
            seg = arr[i0:i1+1]
            mask = ~np.isnan(seg)
            if not mask.any(): continue
            v = float(seg[mask].sum()); n = int(mask.sum())
            sums[r.te_class] = sums.get(r.te_class, 0.0) + v
            counts[r.te_class] = counts.get(r.te_class, 0) + n

        rows = []
        for cls, v in sums.items():
            n = counts[cls]
            if n < 5: continue
            rows.append(dict(species=sp, te_class=cls,
                             mean_pHet_TE=v/n, n_windows=n,
                             bg_pHet=bg_mean,
                             delta=v/n - bg_mean))
        return sp, bg_mean, pd.DataFrame(rows)
    except Exception as ex:
        print(f"[WARN] {sp}: {ex}", file=sys.stderr)
        return sp, None, None


def main():
    ts_files = sorted(glob.glob(os.path.join(CS_DIR, "*_supervised.tsv.gz")))
    jobs = []
    for f in ts_files:
        sp = os.path.basename(f).replace("_supervised.tsv.gz","")
        bed = os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                           f"{sp}.filteredRepeats.bed")
        jobs.append((sp, f, bed))
    print(f"[RUN] {len(jobs)} species")

    rows = []; bg_rows = []
    n_done = 0
    with ProcessPoolExecutor(max_workers=16) as ex:
        futs = [ex.submit(process_species, j) for j in jobs]
        for fut in as_completed(futs):
            sp, bg, d = fut.result()
            n_done += 1
            if bg is not None:
                bg_rows.append((sp, bg))
            if d is not None and len(d):
                rows.append(d)
            if n_done % 25 == 0:
                print(f"  [{n_done}/{len(jobs)}]")
    if not rows:
        print("nothing"); return
    df = pd.concat(rows, ignore_index=True)
    df["phylum_code"] = df["species"].str.split("_", n=1).str[0]
    df["phylum"] = df["phylum_code"].map(PHYLUM).fillna(df["phylum_code"])
    df.to_csv(os.path.join(OUT_ROOT, "delta_pHet_TE_vs_background.tsv"),
              sep="\t", index=False)
    print(f"[OUT] rows={len(df)}  species={df['species'].nunique()}")

    # ---------------- per-species TE-pooled (weighted mean over classes) ----
    sp_pool = (df.groupby(["species","phylum"])
                 .apply(lambda d: pd.Series({
                     "delta_TE_all": np.average(d["delta"],
                                                weights=d["n_windows"]),
                     "n": d["n_windows"].sum(),
                 }))
                 .reset_index())
    sp_pool.to_csv(os.path.join(OUT_ROOT, "delta_pHet_per_species.tsv"),
                   sep="\t", index=False)

    # phylum aggregate
    phy = (sp_pool.groupby("phylum")
                  .agg(mean=("delta_TE_all","mean"),
                       sd  =("delta_TE_all","std"),
                       n   =("delta_TE_all","count"))
                  .reset_index())
    phy["se"] = phy["sd"] / np.sqrt(phy["n"].clip(lower=1))
    # keep phyla with >=2 species
    phy = phy[phy["n"] >= 2].sort_values("mean")
    phy.to_csv(os.path.join(OUT_ROOT, "delta_pHet_by_phylum.tsv"),
               sep="\t", index=False)

    # class aggregate (mean over species per class)
    cls_per_sp = (df.groupby(["te_class","species"])
                    .apply(lambda d: np.average(d["delta"],
                                                weights=d["n_windows"]))
                    .reset_index().rename(columns={0:"delta"}))
    cls = (cls_per_sp.groupby("te_class")
                     .agg(mean=("delta","mean"),
                          sd  =("delta","std"),
                          n   =("delta","count"))
                     .reset_index())
    cls["se"] = cls["sd"] / np.sqrt(cls["n"].clip(lower=1))
    cls = cls[cls["te_class"].isin(CLASS_ORDER)]
    cls = cls.sort_values("mean")
    cls.to_csv(os.path.join(OUT_ROOT, "delta_pHet_by_TEclass.tsv"),
               sep="\t", index=False)

    # ---------------- plot ----------------
    fig, axes = plt.subplots(1, 2, figsize=(15, 0.32 * max(len(phy), len(cls)) + 2),
                             gridspec_kw={"width_ratios":[1.05, 1.0]})
    fig.patch.set_facecolor("#F5F8FB")

    # -------- Panel A: per phylum --------
    ax = axes[0]
    ax.set_facecolor("#F5F8FB")
    palette = plt.cm.Set2.colors + plt.cm.Set3.colors
    colors = [palette[i % len(palette)] for i in range(len(phy))]
    yy = np.arange(len(phy))
    ax.barh(yy, phy["mean"].values, xerr=phy["se"].values,
            color=colors, edgecolor="black", linewidth=0.4,
            error_kw=dict(ecolor="black", capsize=3, lw=0.8))
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels(phy["phylum"], fontsize=10)
    ax.set_xlabel("Δ p_het  (TE − genome background)", fontsize=10)
    ax.set_title("A   Compaction signal at TE insertions vs Genome Background",
                 loc="left", fontsize=11, fontweight="bold")
    # legend box
    txt = ("Δ < 0  →  TEs sit in less-compacted (open) regions\n"
           "Δ > 0  →  TEs sit in more-compacted (heterochromatic) regions")
    ax.text(0.98, 0.04, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, bbox=dict(boxstyle="round,pad=0.4",
                                  fc="white", ec="#9CB4CC"))

    # -------- Panel B: per TE class --------
    ax = axes[1]
    ax.set_facecolor("#F5F8FB")
    cls = cls.set_index("te_class").reindex(
        [c for c in CLASS_ORDER if c in cls["te_class"].values][::-1]
        if False else cls.sort_values("mean")["te_class"].values
    ).reset_index()
    yy = np.arange(len(cls))
    # color map: colored bars for major TE classes, grey for repeats/unknown
    color_map = {"LTR":"#D6604D", "LINE":"#F4A582", "SINE":"#B2182B",
                 "DNA":"#92C5DE", "RC":"#1B7837"}
    bar_colors = [color_map.get(c, "#888888") for c in cls["te_class"]]
    ax.barh(yy, cls["mean"].values, xerr=cls["se"].values,
            color=bar_colors, edgecolor="black", linewidth=0.4,
            error_kw=dict(ecolor="black", capsize=3, lw=0.8))
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels(cls["te_class"], fontsize=10)
    ax.set_xlabel("Δ p_het  (TE − genome background)", fontsize=10)
    ax.set_title("B   Compaction Bias by TE Class (averaged across all species)",
                 loc="left", fontsize=11, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "delta_pHet_TE_vs_background.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(OUT_ROOT, "delta_pHet_TE_vs_background.png"),
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[DONE]")


if __name__ == "__main__":
    main()
