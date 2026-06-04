#!/usr/bin/env python3
"""
delta_nrl_by_clade_and_superfamily.py
======================================
Two complementary figures — same layout as the existing NuclPred bar charts:

  Panel A (left) : Δ metric by clade (mean ± SE across species in clade)
  Panel B (right): Δ metric by TE superfamily + PCG (mean ± SE across species)

Figure 1 — ALL FUNGI (216 species)
  Metric : Δ score  = mean_score(class windows) − mean_score(background windows)
  Source : _states.tsv.gz  (NuclPred-based 2-feature score, 500 bp windows)
  Classes: LTR, LINE, DNA, Helitron, PLE, SINE, Unknown, PCG

Figure 2 — IN VIVO ONLY (4 species with real chromatin data)
  Metric : Δ NRL (bp)  = median_NRL(class windows) − median_NRL(background windows)
  Source : real MNase / MAINE / ATAC WIG files  (5 kb windows, ac_lag)
  Classes: LTR, LINE, DNA, Helitron, PLE, SINE, Unknown, PCG

Panel A in Figure 2 shows per-species bars (no clade pooling — only 4 species).

Background = windows with no TE coverage ≥ MIN_TE_FRAC AND no gene coverage ≥ MIN_GENE_FRAC.

Outputs:
  invivo_landscape/figE1_all_fungi_score.{pdf,png}
  invivo_landscape/figE2_invivo_nrl.{pdf,png}
  invivo_landscape/figE_summary_allfungi.tsv
  invivo_landscape/figE_summary_invivo.tsv
"""
import os, sys, gzip, re
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics

# ----------------------------------------------------------------- paths ---
CS_DIR   = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/simple2f/cross_species_2state"
GFF_DIR  = "/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff"
FAI_DIR  = "/legserv/Temp/Carlos/EDF-TE/data/genomes"
TE_DATA  = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
TE_MV    = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
MV_NUC   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes/nuc-seq"
MV_ATC   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes/atac-seq"
OUT      = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/invivo_landscape"
os.makedirs(OUT, exist_ok=True)

WIN_ALLFUNGI = 500      # bp (matches states TSV)
WIN_INVIVO   = 5_000    # bp (for reliable ac_lag)
MIN_TE_FRAC  = 0.20     # fraction of window covered to classify as TE class
MIN_GENE_FRAC= 0.30     # fraction covered to classify as PCG
MIN_WINDOWS  = 10       # per class per species minimum
N_WORKERS    = 8

# ---------------------------------------------------------------- taxonomy ---
PHYLUM = {
    "AS":"Ascomycota","BA":"Basidiomycota","BL":"Blastocladiomycota",
    "BS":"Basidiobolomycota","CH":"Chytridiomycota","CR":"Cryptomycota",
    "EN":"Entomophthoromycota","GL":"Glomeromycota","KX":"Kickxellomycota",
    "MI":"Microsporidia","MO":"Mortierellomycota","MU":"Mucoromycota",
    "NU":"Nuclear-other","OL":"Olpidiomycota","AP":"Aphelidiomycota",
    "ZO":"Zoopagomycota",
}

IN_VIVO_PANEL = [
    dict(label="MU_Rhimi59",     sp="MU_Rhimi59",  assay="MNase", atac=False,
         wig="/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"),
    dict(label="AS_Ncrass2",     sp="AS_Ncrass2",  assay="MNase", atac=False,
         wig=f"{MV_NUC}/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.wig"),
    dict(label="AS_Bcine1",      sp="AS_Bcine1",   assay="MAINE", atac=False,
         wig=f"{MV_NUC}/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.wig"),
    dict(label="BA_Mereb1_atac", sp="BA_Mereb1",   assay="ATAC",  atac=True,
         wig=f"{MV_ATC}/BA_Mereb1/BA_Mereb1_atac.wig"),
]

# TE superfamily classes (mobile elements only)
SF_CLASSES = ["LTR","LINE","DNA","Helitron","PLE","SINE","Unknown"]
SF_ORDER   = SF_CLASSES + ["PCG"]

SF_COLORS  = {
    "LTR":      "#D6604D",
    "LINE":     "#F4A582",
    "DNA":      "#4393C3",
    "Helitron": "#9B59B6",
    "PLE":      "#1ABC9C",
    "SINE":     "#E67E22",
    "Unknown":  "#7F8C8D",
    "PCG":      "#2ECC71",
}

def classify_te_sf(col4):
    """Return simplified SF name or None to exclude non-TE repeats."""
    raw = col4.strip(); c = raw.split("/")[0].rstrip("?").strip()
    if c.startswith("LTR"):                                        return "LTR"
    if c.startswith("LINE"):                                       return "LINE"
    if c.startswith("DNA"):                                        return "DNA"
    if c.startswith("RC") or c.lower().startswith("helitron"):     return "Helitron"
    if c.startswith("PLE") or c.startswith("Penelope"):            return "PLE"
    if c.startswith("SINE"):                                       return "SINE"
    if c in ("Unknown","Retroposon","Unspecified","Other"):         return "Unknown"
    return None  # exclude Simple_repeat, Low_complexity, Satellite, tRNA, rRNA, snRNA

def te_bed_path(sp):
    a = TE_DATA.format(sp=sp); b = TE_MV.format(sp=sp)
    return a if os.path.exists(a) else (b if os.path.exists(b) else None)

def gff_path(sp):
    for ext in (".gff",".gff3",".gff2"):
        p = os.path.join(GFF_DIR, sp + ext)
        if os.path.exists(p): return p
    return None

# ============================================================== FIGURE 1 ==
# All fungi — score from NuclPred states TSV

def load_dominant_class_arrays(bed_path, gff_p, chroms_sizes, win):
    """
    Returns dict: chrom -> int8 array of length n_windows.
    Values encode dominant class index:
      0 = background
      1..len(SF_CLASSES) = TE class index (0-based in SF_CLASSES)
      len(SF_CLASSES)+1  = PCG
    """
    nsf = len(SF_CLASSES)
    PCG_IDX = nsf + 1
    BG_IDX  = 0

    # Accumulate bp coverage per class per window
    cov = {}
    for c, L in chroms_sizes.items():
        n = L // win
        if n == 0: continue
        cov[c] = np.zeros((n, nsf + 2), dtype=np.float32)  # [bg, sf0..sfN-1, pcg]

    # TE coverage
    if bed_path and os.path.exists(bed_path):
        with open(bed_path) as fh:
            for ln in fh:
                f = ln.split("\t")
                if len(f) < 4: continue
                c = f[0]
                if c not in cov: continue
                try: s, e = int(f[1]), int(f[2])
                except ValueError: continue
                sf = classify_te_sf(f[3])
                if sf is None: continue
                sf_idx = SF_CLASSES.index(sf) + 1  # 1-based
                n = len(cov[c])
                e = min(e, n * win)
                i0 = max(0, s // win); i1 = min((e - 1) // win, n - 1)
                if i0 > i1: continue
                # fractional coverage per window
                for i in range(i0, i1 + 1):
                    ws = i * win; we = ws + win
                    ov = min(e, we) - max(s, ws)
                    if ov > 0:
                        cov[c][i, sf_idx] += ov / win

    # Gene coverage (PCG)
    if gff_p and os.path.exists(gff_p):
        gene_feat = "gene"
        # auto-detect if gene feature exists
        with open(gff_p) as fh:
            for ln in fh:
                if ln.startswith("#"): continue
                f = ln.split("\t")
                if len(f) < 3: continue
                if f[2] == "gene": break
                if f[2] in ("exon","CDS"):
                    gene_feat = f[2]; break
        with open(gff_p) as fh:
            for ln in fh:
                if ln.startswith("#"): continue
                f = ln.split("\t")
                if len(f) < 5 or f[2] != gene_feat: continue
                c = f[0]
                if c not in cov: continue
                try: s, e = int(f[3]) - 1, int(f[4])
                except ValueError: continue
                n = len(cov[c])
                e = min(e, n * win)
                i0 = max(0, s // win); i1 = min((e - 1) // win, n - 1)
                if i0 > i1: continue
                for i in range(i0, i1 + 1):
                    ws = i * win; we = ws + win
                    ov = min(e, we) - max(s, ws)
                    if ov > 0:
                        cov[c][i, PCG_IDX] += ov / win

    # Assign dominant class
    dom = {}
    for c, arr in cov.items():
        n = arr.shape[0]
        cls = np.full(n, BG_IDX, dtype=np.int8)
        # TE dominant: any TE class > MIN_TE_FRAC  → pick max-coverage TE
        te_cols = arr[:, 1:nsf+1]
        te_max  = te_cols.max(axis=1)
        te_which = te_cols.argmax(axis=1) + 1   # 1-based SF index
        te_mask  = te_max >= MIN_TE_FRAC
        cls[te_mask] = te_which[te_mask]
        # PCG: gene coverage > MIN_GENE_FRAC AND NOT TE-dominant
        pcg_mask = (arr[:, PCG_IDX] >= MIN_GENE_FRAC) & (~te_mask)
        cls[pcg_mask] = PCG_IDX
        dom[c] = cls
    return dom


def process_species_score(args):
    sp, ts_gz, bed_p, gff_p = args
    try:
        df = pd.read_csv(ts_gz, sep="\t",
                         usecols=["chrom","start","score","state"],
                         compression="gzip")
        df = df[df["state"].notna() & (df["state"] != "NoSig") & df["score"].notna()]
        if df.empty: return sp, None

        chroms_sizes = {}
        for c, sub in df.groupby("chrom", sort=False):
            max_end = (sub["start"].max() + WIN_ALLFUNGI)
            chroms_sizes[c] = int(max_end)

        dom = load_dominant_class_arrays(bed_p, gff_p, chroms_sizes, WIN_ALLFUNGI)

        # Build score arrays
        sc = {}
        for c, sub in df.groupby("chrom", sort=False):
            n = chroms_sizes[c] // WIN_ALLFUNGI
            arr = np.full(n, np.nan, dtype=np.float32)
            idx = (sub["start"].values // WIN_ALLFUNGI).astype(int)
            mask = idx < n
            arr[idx[mask]] = sub["score"].values[mask].astype(np.float32)
            sc[c] = arr

        nsf = len(SF_CLASSES)
        PCG_IDX = nsf + 1

        # Per-class scores
        class_vals = {sf: [] for sf in SF_CLASSES + ["PCG", "BG"]}
        for c, darr in dom.items():
            sarr = sc.get(c)
            if sarr is None: continue
            n = min(len(darr), len(sarr))
            for i in range(n):
                v = sarr[i]
                if np.isnan(v): continue
                d = darr[i]
                if d == 0:
                    class_vals["BG"].append(v)
                elif 1 <= d <= nsf:
                    class_vals[SF_CLASSES[d-1]].append(v)
                elif d == PCG_IDX:
                    class_vals["PCG"].append(v)

        bg_mean = np.mean(class_vals["BG"]) if len(class_vals["BG"]) >= MIN_WINDOWS else np.nan
        if np.isnan(bg_mean): return sp, None

        rows = []
        for cls in SF_CLASSES + ["PCG"]:
            vals = np.array(class_vals[cls], dtype=np.float32)
            if len(vals) < MIN_WINDOWS: continue
            rows.append({
                "species": sp, "sf_class": cls,
                "mean_score": float(np.mean(vals)),
                "bg_score":   float(bg_mean),
                "delta":      float(np.mean(vals)) - float(bg_mean),
                "n_windows":  len(vals),
                "n_bg":       len(class_vals["BG"]),
            })
        return sp, pd.DataFrame(rows) if rows else None
    except Exception as ex:
        import traceback
        print(f"[WARN] {sp}: {ex}", file=sys.stderr)
        return sp, None


def run_all_fungi():
    print("\n" + "="*60)
    print("FIGURE 1: All fungi — Δ score by clade and superfamily")
    print("="*60)

    # Gather species
    tsv_files = sorted(
        f for f in os.listdir(CS_DIR) if f.endswith("_states.tsv.gz"))
    args_list = []
    for fn in tsv_files:
        sp  = fn.replace("_states.tsv.gz","")
        ts  = os.path.join(CS_DIR, fn)
        bed = te_bed_path(sp)
        gff = gff_path(sp)
        args_list.append((sp, ts, bed, gff))

    print(f"  Processing {len(args_list)} species …")

    all_rows = []
    done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(process_species_score, a): a[0] for a in args_list}
        for fut in as_completed(futs):
            sp, df = fut.result()
            done += 1
            if df is not None:
                all_rows.append(df)
            if done % 20 == 0:
                print(f"  … {done}/{len(args_list)}")

    if not all_rows:
        print("No data!"); return None

    df_all = pd.concat(all_rows, ignore_index=True)
    df_all["phylum_code"] = df_all["species"].str.split("_", n=1).str[0]
    df_all["phylum"]      = df_all["phylum_code"].map(PHYLUM).fillna(df_all["phylum_code"])
    print(f"  {df_all['species'].nunique()} species with data")
    df_all.to_csv(os.path.join(OUT, "figE_summary_allfungi.tsv"), sep="\t", index=False)
    return df_all


def make_figure1(df_all):
    if df_all is None or df_all.empty: return

    # --- Panel A: by clade ---
    sp_pool = (df_all[df_all["sf_class"].isin(SF_CLASSES)]
               .groupby(["species","phylum"])
               .apply(lambda d: np.average(d["delta"], weights=d["n_windows"]))
               .reset_index().rename(columns={0:"delta_TE"}))
    phy = (sp_pool.groupby("phylum")
           .agg(mean=("delta_TE","mean"), sd=("delta_TE","std"), n=("delta_TE","count"))
           .reset_index())
    phy["se"] = phy["sd"] / np.sqrt(phy["n"].clip(lower=1))
    phy = phy[phy["n"] >= 1].sort_values("mean")

    # --- Panel B: by TE superfamily + PCG ---
    cls_per_sp = (df_all.groupby(["sf_class","species"])
                  .apply(lambda d: np.average(d["delta"], weights=d["n_windows"]))
                  .reset_index().rename(columns={0:"delta"}))
    cls = (cls_per_sp.groupby("sf_class")
           .agg(mean=("delta","mean"), sd=("delta","std"), n=("delta","count"))
           .reset_index())
    cls["se"] = cls["sd"] / np.sqrt(cls["n"].clip(lower=1))
    order_present = [c for c in SF_ORDER if c in cls["sf_class"].values]
    cls = cls[cls["sf_class"].isin(order_present)]
    cls = cls.set_index("sf_class").loc[order_present].reset_index()
    cls = cls.sort_values("mean")

    fig, axes = plt.subplots(1, 2, figsize=(16, max(5.5, 0.42*max(len(phy), len(cls))+2.5)),
                             gridspec_kw={"width_ratios":[1.1, 0.9]})
    fig.patch.set_facecolor("#F5F8FB")

    # Panel A
    ax = axes[0]; ax.set_facecolor("#F5F8FB")
    palette = list(plt.cm.Set2.colors) + list(plt.cm.Set3.colors)
    colors = [palette[i % len(palette)] for i in range(len(phy))]
    yy = np.arange(len(phy))
    ax.barh(yy, phy["mean"].values, xerr=phy["se"].values,
            color=colors, edgecolor="black", linewidth=0.5,
            error_kw=dict(ecolor="black", capsize=3, lw=0.9))
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{p}  (n={n})" for p,n in zip(phy["phylum"], phy["n"])], fontsize=10)
    ax.set_xlabel("Δ NuclPred score (TE regions − genome background)", fontsize=10)
    ax.set_title("A   All fungi — by clade", loc="left", fontsize=11, fontweight="bold")
    ax.text(0.97, 0.03,
            "Δ < 0 → TEs in nucleosome-depleted regions\nΔ > 0 → TEs in nucleosome-enriched regions",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#9CB4CC", alpha=0.9))

    # Panel B
    ax = axes[1]; ax.set_facecolor("#F5F8FB")
    bar_colors = [SF_COLORS.get(c, "#888888") for c in cls["sf_class"]]
    yy2 = np.arange(len(cls))
    ax.barh(yy2, cls["mean"].values, xerr=cls["se"].values,
            color=bar_colors, edgecolor="black", linewidth=0.5,
            error_kw=dict(ecolor="black", capsize=3, lw=0.9))
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy2)
    ax.set_yticklabels([f"{c}  (n={n})" for c, n in zip(cls["sf_class"], cls["n"])], fontsize=10)
    ax.set_xlabel("Δ NuclPred score (class − background)", fontsize=10)
    ax.set_title("B   By TE superfamily + PCG  (mean ± SE across species)",
                 loc="left", fontsize=11, fontweight="bold")

    fig.suptitle(f"Nucleosome positioning at TE insertions vs genome background\n"
                 f"NuclPred-based score, {WIN_ALLFUNGI} bp windows, "
                 f"{df_all['species'].nunique()} species",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    for ext in ("pdf","png"):
        kw = {"bbox_inches":"tight"} if ext=="pdf" else {"dpi":200,"bbox_inches":"tight"}
        fig.savefig(os.path.join(OUT, f"figE1_all_fungi_score.{ext}"), **kw)
    plt.close(fig)
    print("  → figE1_all_fungi_score saved")


# ============================================================== FIGURE 2 ==
# In vivo — NRL (ac_lag) from real MNase/ATAC WIG

def chrom_sizes_from_wig_or_fai(sp, wig_path):
    fai = os.path.join(FAI_DIR, f"{sp}.fasta.fai")
    if os.path.exists(fai):
        df = pd.read_csv(fai, sep="\t", header=None, names=["c","L","o","lb","lby"])
        return dict(zip(df["c"], df["L"].astype(int)))
    # fallback: infer from WIG
    from chromatin_states_from_nuclpred import load_nuclpred_wig
    arrs = load_nuclpred_wig(wig_path, None)
    return {c: len(a) for c, a in arrs.items()}


def run_invivo():
    print("\n" + "="*60)
    print("FIGURE 2: In vivo — Δ NRL by species and superfamily (5 kb)")
    print("="*60)

    all_rows = []
    species_data = {}

    for entry in IN_VIVO_PANEL:
        lab = entry["label"]; sp = entry["sp"]; wig = entry["wig"]
        if not os.path.exists(wig):
            print(f"[SKIP] {lab}: {wig}"); continue
        print(f"\n[{lab}]")

        sz   = chrom_sizes_from_wig_or_fai(sp, wig)
        arrs = load_nuclpred_wig(wig, sz)
        bed  = te_bed_path(sp)
        gff  = gff_path(sp)

        dom  = load_dominant_class_arrays(bed, gff, sz, WIN_INVIVO)

        rows_win = []
        for c, occ in arrs.items():
            n = len(occ) // WIN_INVIVO
            if n == 0: continue
            means, cvs, aclags, acvals, _ = window_metrics(occ, WIN_INVIVO)
            darr = dom.get(c, np.zeros(n, dtype=np.int8))[:n]
            rows_win.append(pd.DataFrame({
                "chrom": c, "win_idx": np.arange(n),
                "ac_lag": aclags, "mean_sig": means,
                "dom_class": darr,
            }))
        if not rows_win: continue

        df = pd.concat(rows_win, ignore_index=True)
        df = df[df["mean_sig"] > 1e-3].copy()
        species_data[lab] = (df, entry)

        nsf = len(SF_CLASSES); PCG_IDX = nsf + 1
        bg_vals = df.loc[df["dom_class"]==0, "ac_lag"].dropna().values
        if len(bg_vals) < MIN_WINDOWS: continue
        bg_med = float(np.median(bg_vals))

        print(f"  BG NRL median = {bg_med:.1f} bp  (n={len(bg_vals)})")
        for ci, cls in enumerate(SF_CLASSES + ["PCG"], start=1):
            idx = ci if cls != "PCG" else PCG_IDX
            vals = df.loc[df["dom_class"]==idx, "ac_lag"].dropna().values
            if len(vals) < MIN_WINDOWS: continue
            med = float(np.median(vals))
            print(f"    {cls:12s}: NRL={med:.1f} bp  Δ={med-bg_med:+.1f}  n={len(vals)}")
            all_rows.append({
                "species": lab, "assay": entry["assay"],
                "phylum_code": sp.split("_")[0],
                "phylum": PHYLUM.get(sp.split("_")[0], sp.split("_")[0]),
                "sf_class": cls,
                "median_nrl": med,
                "bg_nrl": bg_med,
                "delta_nrl": med - bg_med,
                "n_windows": len(vals),
                "n_bg": len(bg_vals),
            })

    if not all_rows:
        print("No in-vivo data!"); return None, None

    df_iv = pd.DataFrame(all_rows)
    df_iv.to_csv(os.path.join(OUT, "figE_summary_invivo.tsv"), sep="\t", index=False)
    return df_iv, species_data


def make_figure2(df_iv, species_data):
    if df_iv is None or df_iv.empty: return

    # --- Panel A: per species (4 bars, showing TE-pooled Δ NRL) ---
    sp_pool = (df_iv[df_iv["sf_class"].isin(SF_CLASSES)]
               .groupby(["species","phylum","phylum_code","assay"])
               .apply(lambda d: np.average(d["delta_nrl"], weights=d["n_windows"]))
               .reset_index().rename(columns={0:"delta_nrl_pooled"}))
    sp_pool = sp_pool.sort_values("delta_nrl_pooled")

    # --- Panel B: by superfamily + PCG ---
    cls = (df_iv.groupby("sf_class")
           .agg(mean=("delta_nrl","mean"), sd=("delta_nrl","std"), n=("delta_nrl","count"))
           .reset_index())
    cls["se"] = cls["sd"] / np.sqrt(cls["n"].clip(lower=1))
    present = [c for c in SF_ORDER if c in cls["sf_class"].values]
    cls = cls[cls["sf_class"].isin(present)].set_index("sf_class").loc[present].reset_index()
    cls = cls.sort_values("mean")

    fig, axes = plt.subplots(1, 2, figsize=(15, max(5.5, 0.55*max(len(sp_pool), len(cls))+3.0)),
                             gridspec_kw={"width_ratios":[0.8, 1.1]})
    fig.patch.set_facecolor("#F5F8FB")

    # Panel A — per species
    ax = axes[0]; ax.set_facecolor("#F5F8FB")
    sp_colors = {"MU":"#27AE60","AS":"#E74C3C","BA":"#3498DB","BL":"#9B59B6"}
    colors_a = [sp_colors.get(r["phylum_code"], "#888888")
                for _, r in sp_pool.iterrows()]
    yy = np.arange(len(sp_pool))
    ax.barh(yy, sp_pool["delta_nrl_pooled"].values,
            color=colors_a, edgecolor="black", linewidth=0.5)
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{r['species']}\n({r['assay']})"
                        for _, r in sp_pool.iterrows()], fontsize=10)
    ax.set_xlabel("Δ NRL (bp)  (TE windows − background)", fontsize=10)
    ax.set_title("A   In vivo — by species\n(TE-pooled Δ NRL)",
                 loc="left", fontsize=11, fontweight="bold")
    ax.text(0.97, 0.03,
            "Δ < 0 → shorter NRL over TEs\n(nucleosomes more closely packed)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#9CB4CC", alpha=0.9))

    # Panel B — by superfamily + PCG
    ax = axes[1]; ax.set_facecolor("#F5F8FB")
    bar_colors = [SF_COLORS.get(c, "#888888") for c in cls["sf_class"]]
    yy2 = np.arange(len(cls))
    ax.barh(yy2, cls["mean"].values, xerr=cls["se"].values,
            color=bar_colors, edgecolor="black", linewidth=0.5,
            error_kw=dict(ecolor="black", capsize=3, lw=0.9))
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(yy2)
    ax.set_yticklabels([f"{c}  (n={n})" for c, n in zip(cls["sf_class"], cls["n"])], fontsize=10)
    ax.set_xlabel("Δ NRL (bp)  (class − background)", fontsize=10)
    ax.set_title("B   By TE superfamily + PCG  (mean ± SE, 4 in vivo species)",
                 loc="left", fontsize=11, fontweight="bold")

    fig.suptitle(f"Nucleosome repeat length (NRL) at TE insertions vs genome background\n"
                 f"Real chromatin data (MNase/MAINE/ATAC), {WIN_INVIVO//1000} kb windows",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    for ext in ("pdf","png"):
        kw = {"bbox_inches":"tight"} if ext=="pdf" else {"dpi":200,"bbox_inches":"tight"}
        fig.savefig(os.path.join(OUT, f"figE2_invivo_nrl.{ext}"), **kw)
    plt.close(fig)
    print("  → figE2_invivo_nrl saved")


# ============================================================== FIGURE 3 ==
# Side-by-side comparison: NuclPred (all fungi) vs NRL (in vivo)
# Panel B only (by superfamily + PCG), two bars per class

def make_figure3(df_all, df_iv):
    if df_all is None or df_iv is None: return

    # NuclPred (all fungi, Panel B data)
    cls_all = (df_all.groupby(["sf_class","species"])
               .apply(lambda d: np.average(d["delta"], weights=d["n_windows"]))
               .reset_index().rename(columns={0:"delta_score"}))
    agg_all = (cls_all.groupby("sf_class")
               .agg(mean=("delta_score","mean"), se=("delta_score", lambda x: x.std()/np.sqrt(len(x))))
               .reset_index())
    agg_all["source"] = "All fungi\n(NuclPred score)"

    # NRL in vivo
    agg_iv = (df_iv.groupby("sf_class")
              .agg(mean=("delta_nrl","mean"), se=("delta_nrl", lambda x: x.std()/np.sqrt(len(x))))
              .reset_index())
    agg_iv["source"] = "In vivo\n(Δ NRL bp)"

    present = [c for c in SF_ORDER
               if c in agg_all["sf_class"].values and c in agg_iv["sf_class"].values]
    if not present: return

    fig, ax = plt.subplots(figsize=(8, 5.5))
    fig.patch.set_facecolor("#F5F8FB"); ax.set_facecolor("#F5F8FB")

    n = len(present); y = np.arange(n)
    bh = 0.38

    for agg, offset, hatch, alpha in [
            (agg_all, -bh/2, "", 0.85),
            (agg_iv,  +bh/2, "///", 0.65),
    ]:
        vals = []; errs = []
        for cls in present:
            row = agg[agg["sf_class"]==cls]
            if row.empty:
                vals.append(0.0); errs.append(0.0)
            else:
                vals.append(float(row["mean"].iloc[0]))
                errs.append(float(row["se"].iloc[0]))

        # Normalise to [-1,1] for visual comparison (z-score by range)
        arr = np.array(vals)
        rng = max(abs(arr.max()), abs(arr.min()), 1e-9)
        arr_norm = arr / rng

        colors = [SF_COLORS.get(c, "#888888") for c in present]
        ax.barh(y + offset, arr_norm, height=bh,
                color=colors, edgecolor="black", linewidth=0.5,
                hatch=hatch, alpha=alpha,
                label=agg["source"].iloc[0])

    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(present, fontsize=11)
    ax.set_xlabel("Δ metric (normalised to max per dataset)", fontsize=10)
    ax.set_title("NuclPred (all fungi, solid) vs NRL in vivo (hatched)\n"
                 "Δ < 0 → depleted / shorter NRL over TE class vs background",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="x", alpha=0.3, lw=0.6)

    fig.tight_layout()
    for ext in ("pdf","png"):
        kw = {"bbox_inches":"tight"} if ext=="pdf" else {"dpi":200,"bbox_inches":"tight"}
        fig.savefig(os.path.join(OUT, f"figE3_comparison.{ext}"), **kw)
    plt.close(fig)
    print("  → figE3_comparison saved")


# ======================================================================= main
def main():
    df_all = run_all_fungi()
    make_figure1(df_all)

    df_iv, species_data = run_invivo()
    make_figure2(df_iv, species_data)

    make_figure3(df_all, df_iv)

    print(f"\n[DONE] outputs in {OUT}")

if __name__ == "__main__":
    main()
