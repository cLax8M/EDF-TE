#!/usr/bin/env python3
"""
delta_score_subgroups.py
========================
Re-uses cross_species_2state/<SP>_states.tsv.gz (simple 2-feature score) to
compute Δ score (TE − background) under several SUBGROUP filters:

  ALL                : every TE annotation in the GFF (baseline; previous figure)
  FL                 : full-length copies only
                       cov = (tend - tstart + 1) / consensus_len >= 0.80
  FL_recent          : full-length AND kimura80 < 5
  MNASE6             : restrict species to the 6 with real DANPOS MNase
  RHIMI59            : restrict to Rhizopus microsporus only

For each subgroup we generate a 2-panel figure (Panel A: per phylum; Panel B:
per TE class) like the original.

Outputs in
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/simple2f/subgroups/
"""
import os, sys, glob, re, gzip, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT     = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/simple2f"
CS_DIR   = os.path.join(ROOT, "cross_species_2state")
TE_DIR   = "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_DIR  = os.path.join(ROOT, "subgroups")
os.makedirs(OUT_DIR, exist_ok=True)

WIN = 500
FL_FRAC = 0.80
RECENT_KIMURA = 5.0

PHYLUM = {
    "AS":"Ascomycota","BA":"Basidiomycota","BL":"Blastocladiomycota",
    "BS":"Basidiobolomycota","CH":"Chytridiomycota","CR":"Cryptomycota",
    "EN":"Entomophthoromycota","GL":"Glomeromycota","KX":"Kickxellomycota",
    "MI":"Microsporidia","MO":"Mortierellomycota","MU":"Mucoromycota",
    "NU":"Nuclear-other","OL":"Olpidiomycota","AP":"Aphelidiomycota",
    "ZO":"Zoopagomycota",
}

MNASE6 = {"MU_Rhimi59","AS_Usnflo1","BA_Mereb1","BA_Mrafri1","BA_Rhoto","BL_Catan2"}

CLASS_ORDER = ["LTR","LINE","SINE","DNA","RC","Satellite","Simple","Low",
               "rRNA","tRNA","Unknown"]
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

# ---------------- helpers ----------------
NAME_RE  = re.compile(r"Name=([^;]+)")
KIM_RE   = re.compile(r"kimura80=(-?[\d.]+)")
TSTART   = re.compile(r"tstart=(\d+)")
TEND     = re.compile(r"tend=(\d+)")

def consensus_lengths(fa_path):
    """Return dict UPPER(name without #class) -> length."""
    if not os.path.exists(fa_path):
        return {}
    out = {}; name = None; L = 0
    with open(fa_path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None: out[name] = L
                name = ln[1:].rstrip().split()[0].split("#")[0].upper()
                L = 0
            else:
                L += len(ln.rstrip())
        if name is not None: out[name] = L
    return out


def parse_gff_with_filter(gff_path, fa_path, mode):
    """Yield (chrom,start,end,superfam,te_class) under filter mode."""
    if not os.path.exists(gff_path):
        return []
    cons = consensus_lengths(fa_path) if mode in ("FL","FL_recent") else None
    rows = []
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip(): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9: continue
            try:
                s = int(f[3]) - 1; e = int(f[4])
            except ValueError: continue
            sf = f[2]
            attr = f[8]
            # filters
            if mode in ("FL","FL_recent"):
                m = NAME_RE.search(attr); ts = TSTART.search(attr); te = TEND.search(attr)
                if not (m and ts and te): continue
                fam = m.group(1).upper()
                clen = cons.get(fam, 0)
                if clen <= 0: continue
                cov = (int(te.group(1)) - int(ts.group(1)) + 1) / float(clen)
                if cov < FL_FRAC: continue
            if mode == "FL_recent":
                k = KIM_RE.search(attr)
                if not k: continue
                try:
                    if float(k.group(1)) >= RECENT_KIMURA: continue
                except ValueError: continue
            rows.append((f[0], s, e, sf, te_class(sf)))
    return rows


def load_score_array(ts_path):
    """Return dict chrom -> np.float32 array (NaN for non-informative)."""
    df = pd.read_csv(ts_path, sep="\t",
                     usecols=["chrom","start","end","mean_occ","score","state"])
    df = df[df["state"] != "NoSig"]
    if df.empty: return None, None
    bg = float(df["score"].mean())
    out = {}
    for c, sub in df.groupby("chrom", sort=False):
        idx = (sub["start"].values // WIN).astype(np.int64)
        n = int(idx.max()) + 1
        arr = np.full(n, np.nan, dtype=np.float32)
        arr[idx] = sub["score"].values.astype(np.float32)
        out[c] = arr
    return bg, out


def process_species(args):
    sp, ts, gff, fa, mode = args
    try:
        bg, ca = load_score_array(ts)
        if ca is None: return sp, None, None
        rows = parse_gff_with_filter(gff, fa, mode)
        if not rows: return sp, bg, None
        sums={}; counts={}
        for c,s,e,_sf,cls in rows:
            arr = ca.get(c)
            if arr is None: continue
            i0 = s // WIN; i1 = min((e-1)//WIN, len(arr)-1)
            if i0 > i1: continue
            seg = arr[i0:i1+1]; mask = ~np.isnan(seg)
            if not mask.any(): continue
            sums[cls] = sums.get(cls, 0.0) + float(seg[mask].sum())
            counts[cls] = counts.get(cls, 0) + int(mask.sum())
        recs=[]
        for cls,v in sums.items():
            n = counts[cls]
            if n < 5: continue
            recs.append(dict(species=sp, te_class=cls, mean_score_TE=v/n,
                             n_windows=n, bg_score=bg, delta=v/n - bg))
        return sp, bg, pd.DataFrame(recs) if recs else None
    except Exception as ex:
        print(f"[WARN] {sp} {mode}: {ex}", file=sys.stderr)
        return sp, None, None


def aggregate_and_plot(df, label, out_pdf, out_png):
    if df is None or df.empty:
        print(f"[SKIP] {label}: no data"); return None, None
    df = df.copy()
    df["phylum_code"] = df["species"].str.split("_", n=1).str[0]
    df["phylum"] = df["phylum_code"].map(PHYLUM).fillna(df["phylum_code"])

    # per-species TE-pooled (weighted)
    sp_pool = (df.groupby(["species","phylum"])
                 .apply(lambda d: pd.Series({
                     "delta_TE_all": np.average(d["delta"], weights=d["n_windows"])
                 })).reset_index())
    phy = (sp_pool.groupby("phylum")
                  .agg(mean=("delta_TE_all","mean"),
                       sd  =("delta_TE_all","std"),
                       n   =("delta_TE_all","count"))
                  .reset_index())
    phy["se"] = phy["sd"] / np.sqrt(phy["n"].clip(lower=1))
    phy = phy[phy["n"] >= 1].sort_values("mean")

    cls_per_sp = (df.groupby(["te_class","species"])
                    .apply(lambda d: np.average(d["delta"], weights=d["n_windows"]))
                    .reset_index().rename(columns={0:"delta"}))
    cls = (cls_per_sp.groupby("te_class")
                     .agg(mean=("delta","mean"),
                          sd  =("delta","std"),
                          n   =("delta","count"))
                     .reset_index())
    cls["se"] = cls["sd"] / np.sqrt(cls["n"].clip(lower=1))
    cls = cls[cls["te_class"].isin(CLASS_ORDER)].sort_values("mean")

    # plot
    fig, axes = plt.subplots(1, 2, figsize=(15, 0.32*max(len(phy),len(cls))+2.2),
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
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{p}  (n={n})" for p,n in zip(phy["phylum"], phy["n"])],
                       fontsize=10)
    ax.set_xlabel("Δ score (TE − background)", fontsize=10)
    ax.set_title(f"A   Compaction signal at TE insertions  —  {label}",
                 loc="left", fontsize=11, fontweight="bold")
    ax.text(0.98,0.04,
            "score = occ_pct − cv_pct\nΔ > 0 → more compacted\nΔ < 0 → more open",
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
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{c}  (n={n})" for c,n in zip(cls["te_class"], cls["n"])],
                       fontsize=10)
    ax.set_xlabel("Δ score (TE − background)", fontsize=10)
    ax.set_title(f"B   By TE class (mean ± SE across species)  —  {label}",
                 loc="left", fontsize=11, fontweight="bold")

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return phy, cls


def run_subgroup(name, mode, sp_filter):
    """sp_filter: None | set of species names."""
    files = sorted(glob.glob(os.path.join(CS_DIR, "*_states.tsv.gz")))
    jobs=[]
    for f in files:
        sp = os.path.basename(f).replace("_states.tsv.gz","")
        if sp_filter is not None and sp not in sp_filter:
            continue
        gff = os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                           f"{sp}.filteredRepeats.gff")
        fa  = os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                           f"{sp}_combined_library.fasta")
        jobs.append((sp, f, gff, fa, mode))
    print(f"[{name}] mode={mode}  species={len(jobs)}", flush=True)
    rows=[]; n=0
    nproc = min(16, max(1, len(jobs)))
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        futs=[ex.submit(process_species, j) for j in jobs]
        for fut in as_completed(futs):
            sp, bg, d = fut.result(); n+=1
            if d is not None and len(d): rows.append(d)
            if n % max(1, len(jobs)//10) == 0:
                print(f"  [{name} {n}/{len(jobs)}]", flush=True)
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    df.to_csv(os.path.join(OUT_DIR, f"delta_{name}.tsv"), sep="\t", index=False)
    phy, cls = aggregate_and_plot(df, name,
                                  os.path.join(OUT_DIR, f"delta_{name}.pdf"),
                                  os.path.join(OUT_DIR, f"delta_{name}.png"))
    if phy is not None:
        phy.to_csv(os.path.join(OUT_DIR, f"delta_{name}_byPhylum.tsv"), sep="\t", index=False)
        cls.to_csv(os.path.join(OUT_DIR, f"delta_{name}_byClass.tsv"), sep="\t", index=False)
        print(f"[{name}] phylum=\n{phy.to_string(index=False)}", flush=True)
        print(f"[{name}] class=\n{cls.to_string(index=False)}", flush=True)


def main():
    # 1) Full-length only, ALL species
    run_subgroup("FL_allSP",       mode="FL",        sp_filter=None)
    # 2) Full-length AND recent (kimura80<5), ALL species
    run_subgroup("FLrecent_allSP", mode="FL_recent", sp_filter=None)
    # 3) ALL TEs, only the 6 MNase species
    run_subgroup("ALL_MNASE6",     mode="ALL",       sp_filter=MNASE6)
    # 4) ALL TEs, only Rhizopus
    run_subgroup("ALL_RHIMI59",    mode="ALL",       sp_filter={"MU_Rhimi59"})


if __name__ == "__main__":
    main()
