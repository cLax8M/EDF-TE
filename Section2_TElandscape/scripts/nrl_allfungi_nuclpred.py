#!/usr/bin/env python3
"""
nrl_allfungi_nuclpred.py
=========================
Compute Δ NRL (nucleosome repeat length, from autocorrelation of NuclPred
occupancy in 10 kb windows) for all ~217 fungi and plot as horizontal bar
charts by clade and by TE superfamily + PCG.

This is the in silico companion to the in vivo MNase/MAINE analysis:
the signal is nucleosome-occupancy predicted from DNA sequence alone
(NuclPred), so any Δ NRL reflects the DNA-sequence-encoded chromatin
landscape — what the sequence "wants" independent of epigenetic marks.

  Δ NRL = median(ac_lag in class windows) − median(ac_lag in BG windows)
  where BG = windows with <20 % TE coverage AND <30 % gene coverage.

Outputs (invivo_landscape/):
  figG1_allfungi_nrl.{pdf,png}
  figG1_summary_allfungi_nrl.tsv
"""
import os, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics
from clade_palette import CLADE_COLORS as _CP, CLADE_NAMES, CLADE_ORDER as _CO, TE_CLASS_COLORS, get_clade

# ----------------------------------------------------------------- paths ---
NUC_DIR  = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"
GFF_DIR  = "/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff"
FAI_DIR  = "/legserv/Temp/Carlos/EDF-TE/data/genomes"
TE_DATA  = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
TE_MV    = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
OUT      = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/invivo_landscape"
os.makedirs(OUT, exist_ok=True)

WIN      = 10_000      # window size (bp)
MIN_TE_FRAC  = 0.20   # window classified as TE-dominant if this fraction covered
MIN_GENE_FRAC= 0.30   # window classified as PCG if this fraction covered
MIN_WINDOWS  = 10     # min windows per class to include
N_WORKERS    = 8

# ---------------------------------------------------------------- taxonomy ---
PHYLUM = {
    "AS":"Ascomycota","BA":"Basidiomycota","BL":"Blastocladiomycota",
    "BS":"Basidiobolomycota","CH":"Chytridiomycota","CR":"Cryptomycota",
    "EN":"Entomophthoromycota","GL":"Glomeromycota","KX":"Kickxellomycota",
    "MI":"Microsporidia","MO":"Mortierellomycota","MU":"Mucoromycota",
    "NU":"Nuclear-other","OL":"Olpidiomycota","AP":"Aphelidiomycota",
    "ZO":"Zoopagomycota","HO":"Holozoa",
}

PHYLUM_ORDER = ["Ascomycota","Basidiomycota","Basidiobolomycota","Blastocladiomycota",
                "Chytridiomycota","Cryptomycota","Entomophthoromycota","Glomeromycota",
                "Kickxellomycota","Microsporidia","Mortierellomycota","Mucoromycota",
                "Zoopagomycota","Aphelidiomycota","Olpidiomycota","Holozoa"]

CLADE_COLORS = {
    "Ascomycota":        "#E74C3C","Basidiomycota":      "#3498DB",
    "Blastocladiomycota":"#9B59B6","Basidiobolomycota":  "#1ABC9C",
    "Chytridiomycota":   "#F39C12","Cryptomycota":       "#7F8C8D",
    "Entomophthoromycota":"#E91E63","Glomeromycota":     "#27AE60",
    "Kickxellomycota":   "#FF5722","Microsporidia":      "#795548",
    "Mortierellomycota": "#607D8B","Mucoromycota":       "#4CAF50",
    "Zoopagomycota":     "#FF9800","Aphelidiomycota":    "#00BCD4",
    "Olpidiomycota":     "#8BC34A","Holozoa":            "#BDC3C7",
}

SF_CLASSES = ["LTR","LINE","DNA","Helitron","PLE","SINE","Unknown"]
SF_ORDER   = SF_CLASSES + ["PCG"]
SF_COLORS  = {
    "LTR":      TE_CLASS_COLORS["LTR"],
    "LINE":     TE_CLASS_COLORS["LINE"],
    "DNA":      TE_CLASS_COLORS["DNA"],
    "Helitron": TE_CLASS_COLORS.get("Rolling Circle","#ff7f00"),
    "PLE":      TE_CLASS_COLORS["PLE"],
    "SINE":     TE_CLASS_COLORS["SINE"],
    "Unknown":  TE_CLASS_COLORS["Unknown"],
    "PCG":      "#333333",
}

# ---------------------------------------------------------------- helpers ---
def classify_te_sf(col4):
    raw = col4.strip(); c = raw.split("/")[0].rstrip("?").strip()
    if c.startswith("LTR"):                                      return "LTR"
    if c.startswith("LINE"):                                     return "LINE"
    if c.startswith("DNA"):                                      return "DNA"
    if c.startswith("RC") or c.lower().startswith("helitron"):  return "Helitron"
    if c.startswith("PLE") or c.startswith("Penelope"):         return "PLE"
    if c.startswith("SINE"):                                     return "SINE"
    if c in ("Unknown","Retroposon","Unspecified","Other"):      return "Unknown"
    return None

def te_bed_path(sp):
    a = TE_DATA.format(sp=sp); b = TE_MV.format(sp=sp)
    return a if os.path.exists(a) else (b if os.path.exists(b) else None)

def gff_path(sp):
    for ext in (".gff",".gff3",".gff2"):
        p = os.path.join(GFF_DIR, sp + ext)
        if os.path.exists(p): return p
    return None

def chrom_sizes_fai(sp):
    fai = os.path.join(FAI_DIR, f"{sp}.fasta.fai")
    if not os.path.exists(fai): return None
    df = pd.read_csv(fai, sep="\t", header=None, names=["c","L","o","lb","lby"])
    return dict(zip(df["c"], df["L"].astype(int)))


def build_dominant_class_vec(bed_path, gff_p, chrom_windows, win):
    """
    Returns dict chrom -> int8 array (length = n_windows).
    0=BG, 1..7=SF_CLASSES index+1, 8=PCG
    """
    nsf = len(SF_CLASSES); PCG_IDX = nsf + 1
    cov = {c: np.zeros((n, nsf + 2), dtype=np.float32)
           for c, n in chrom_windows.items() if n > 0}

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
                sf_idx = SF_CLASSES.index(sf) + 1
                n_win = cov[c].shape[0]
                e = min(e, n_win * win)
                i0 = max(0, s // win); i1 = min((e - 1) // win, n_win - 1)
                if i0 > i1: continue
                for i in range(i0, i1 + 1):
                    ws = i * win; we = ws + win
                    ov = min(e, we) - max(s, ws)
                    if ov > 0: cov[c][i, sf_idx] += ov / win

    if gff_p and os.path.exists(gff_p):
        gene_feat = "gene"
        with open(gff_p) as fh:
            for ln in fh:
                if ln.startswith("#"): continue
                f = ln.split("\t")
                if len(f) < 3: continue
                if f[2] == "gene": break
                if f[2] in ("exon","CDS"): gene_feat = f[2]; break
        with open(gff_p) as fh:
            for ln in fh:
                if ln.startswith("#"): continue
                f = ln.split("\t")
                if len(f) < 5 or f[2] != gene_feat: continue
                c = f[0]
                if c not in cov: continue
                try: s, e = int(f[3]) - 1, int(f[4])
                except ValueError: continue
                n_win = cov[c].shape[0]
                e = min(e, n_win * win)
                i0 = max(0, s // win); i1 = min((e - 1) // win, n_win - 1)
                if i0 > i1: continue
                for i in range(i0, i1 + 1):
                    ws = i * win; we = ws + win
                    ov = min(e, we) - max(s, ws)
                    if ov > 0: cov[c][i, PCG_IDX] += ov / win

    dom = {}
    for c, arr in cov.items():
        n = arr.shape[0]; cls = np.zeros(n, dtype=np.int8)
        te_cols = arr[:, 1:nsf+1]; te_max = te_cols.max(axis=1)
        te_which = te_cols.argmax(axis=1) + 1
        te_mask = te_max >= MIN_TE_FRAC
        cls[te_mask] = te_which[te_mask]
        pcg_mask = (arr[:, PCG_IDX] >= MIN_GENE_FRAC) & (~te_mask)
        cls[pcg_mask] = PCG_IDX
        dom[c] = cls
    return dom


# ======================================================= per-species work ===

def process_one(args):
    sp, wig_path, bed_p, gff_p = args
    nsf = len(SF_CLASSES); PCG_IDX = nsf + 1

    try:
        sz = chrom_sizes_fai(sp)
        arrs = load_nuclpred_wig(wig_path, sz)  # sz may be None → inferred inside
        if sz is None:
            sz = {c: len(a) for c, a in arrs.items()}

        chrom_windows = {c: L // WIN for c, L in sz.items() if L >= WIN}
        if not chrom_windows: return sp, None

        dom = build_dominant_class_vec(bed_p, gff_p, chrom_windows, WIN)

        # collect 4 metrics per class
        MKEYS = ("ac_lag", "cv", "ac_val", "mean_sig")
        class_vals = {cls: {m: [] for m in MKEYS}
                      for cls in SF_CLASSES + ["PCG","BG"]}

        for c, occ in arrs.items():
            n_win = len(occ) // WIN
            if n_win == 0: continue
            darr = dom.get(c)
            if darr is None: continue
            n_win = min(n_win, len(darr))

            means, cvs, aclags, acvals, valids = window_metrics(occ, WIN)
            means   = means[:n_win];  cvs    = cvs[:n_win]
            aclags  = aclags[:n_win]; acvals = acvals[:n_win]; valids = valids[:n_win]

            for i in range(n_win):
                if not valids[i]: continue
                vals_i = {"ac_lag": aclags[i], "cv": cvs[i],
                          "ac_val": acvals[i], "mean_sig": means[i]}
                if any(np.isnan(v) for v in vals_i.values()): continue
                d = darr[i]
                if d == 0:             bucket = "BG"
                elif 1 <= d <= nsf:    bucket = SF_CLASSES[d-1]
                elif d == PCG_IDX:     bucket = "PCG"
                else:                  continue
                for m, v in vals_i.items():
                    class_vals[bucket][m].append(v)

        bg_n = len(class_vals["BG"]["ac_lag"])
        if bg_n < MIN_WINDOWS: return sp, None
        bg = {m: float(np.median(class_vals["BG"][m])) for m in MKEYS}

        rows = []
        for cls in SF_CLASSES + ["PCG"]:
            n = len(class_vals[cls]["ac_lag"])
            if n < MIN_WINDOWS: continue
            row = {"species": sp, "sf_class": cls, "n_windows": n}
            for m in MKEYS:
                v = float(np.median(class_vals[cls][m]))
                row[f"median_{m}"] = v
                row[f"bg_{m}"]     = bg[m]
                row[f"delta_{m}"]  = v - bg[m]
            rows.append(row)

        return sp, pd.DataFrame(rows) if rows else None

    except Exception as ex:
        print(f"[WARN] {sp}: {ex}", file=sys.stderr)
        return sp, None


# ============================================================= main loop ===

def run():
    print("\n" + "="*60)
    print("NRL for all fungi — NuclPred WIG → 10 kb windows")
    print("="*60)

    wig_files = sorted(f for f in os.listdir(NUC_DIR) if f.endswith("_nuclpred.wig"))
    args_list = []
    for wf in wig_files:
        sp = wf.replace("_nuclpred.wig","")
        args_list.append((sp,
                          os.path.join(NUC_DIR, wf),
                          te_bed_path(sp),
                          gff_path(sp)))

    print(f"  {len(args_list)} species to process …")
    all_rows = []; done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(process_one, a): a[0] for a in args_list}
        for fut in as_completed(futs):
            sp, df = fut.result(); done += 1
            if df is not None: all_rows.append(df)
            if done % 30 == 0: print(f"    … {done}/{len(args_list)}")

    if not all_rows: print("No data!"); return None
    df_all = pd.concat(all_rows, ignore_index=True)
    df_all["phylum_code"] = df_all["species"].apply(get_clade)
    df_all["phylum"]      = df_all["phylum_code"].map(PHYLUM).fillna(df_all["phylum_code"])
    print(f"  {df_all['species'].nunique()} species with data")

    tsv = os.path.join(OUT, "figG1_summary_allfungi_multimetric.tsv")
    df_all.to_csv(tsv, sep="\t", index=False)
    print(f"  TSV saved: {tsv}")
    return df_all


def make_figure(df_all):
    if df_all is None or df_all.empty: return

    METRICS = [
        ("delta_ac_lag",   "Δ NRL (bp)",            "Nucleosome spacing\n(Δ < 0 = closer packing)"),
        ("delta_cv",       "Δ CV (rugosity)",        "Chromatin irregularity\n(Δ > 0 = more heterogeneous)"),
        ("delta_ac_val",   "Δ Phasing (AC peak)",    "Array regularity\n(Δ > 0 = stronger phasing)"),
        ("delta_mean_sig", "Δ Mean occupancy",        "NuclPred occupancy\n(Δ < 0 = less nucleosome-favorable)"),
    ]

    # ── Clade aggregation ─────────────────────────────────────────────────
    # Per-species weighted mean Δ across TE SF classes
    phy_agg = {}
    for m, _, _ in METRICS:
        sp_pool = (df_all[df_all["sf_class"].isin(SF_CLASSES)]
                   .groupby(["species","phylum"])
                   .apply(lambda d: np.average(d[m], weights=d["n_windows"].clip(1)))
                   .reset_index().rename(columns={0: "val"}))
        phy = (sp_pool.groupby("phylum")
               .agg(mean=("val","mean"), sd=("val","std"), n=("val","count"))
               .reset_index())
        phy["se"] = phy["sd"] / np.sqrt(phy["n"].clip(1))
        phy_agg[m] = phy

    # common phylum order: sort by delta_ac_lag
    base_phy = phy_agg["delta_ac_lag"]
    all_phyla = [p for p in PHYLUM_ORDER if p in base_phy["phylum"].values] + \
                [p for p in base_phy["phylum"].values if p not in PHYLUM_ORDER]
    order_map  = {p: float(base_phy.loc[base_phy["phylum"]==p,"mean"].values[0])
                  for p in all_phyla if p in base_phy["phylum"].values}
    present_phyla = sorted(order_map.keys(), key=lambda p: order_map[p])
    # Exclude clades with very few species or uncertain phylogenetic placement
    EXCLUDE_CLADES = {"Cryptomycota", "Aphelidiomycota"}
    present_phyla = [p for p in present_phyla if p not in EXCLUDE_CLADES]

    # ── Superfamily aggregation ───────────────────────────────────────────
    sf_agg = {}
    for m, _, _ in METRICS:
        cls_per_sp = (df_all.groupby(["sf_class","species"])
                      .apply(lambda d: float(np.average(d[m], weights=d["n_windows"].clip(1))))
                      .reset_index().rename(columns={0:"val"}))
        sf = (cls_per_sp.groupby("sf_class")
              .agg(mean=("val","mean"), sd=("val","std"), n=("val","count"))
              .reset_index())
        sf["se"] = sf["sd"] / np.sqrt(sf["n"].clip(1))
        sf_nwin = df_all.groupby("sf_class")["n_windows"].sum().reset_index()\
                        .rename(columns={"n_windows":"n_win_total"})
        sf = sf.merge(sf_nwin, on="sf_class", how="left")
        sf_agg[m] = sf

    base_sf = sf_agg["delta_ac_lag"]
    sf_ord_map = {s: float(base_sf.loc[base_sf["sf_class"]==s,"mean"].values[0])
                  for s in SF_ORDER if s in base_sf["sf_class"].values}
    present_sf = sorted(sf_ord_map.keys(), key=lambda s: sf_ord_map[s])
    present_sf = [s for s in present_sf if s != "PLE"]   # remove Penelope
    SF_DISPLAY = {s: ("Rolling Circle" if s == "Helitron" else s) for s in SF_ORDER}

    # ── Figure: 4 rows × 2 cols ───────────────────────────────────────────
    n_metrics = len(METRICS)
    n_phyla   = len(present_phyla)
    n_sf      = len(present_sf)
    row_h     = max(0.42 * max(n_phyla, n_sf), 2.5)

    fig = plt.figure(figsize=(18, row_h * n_metrics + 1.5))
    fig.patch.set_facecolor("#F6F8FB")
    gs  = gridspec.GridSpec(n_metrics, 2, figure=fig,
                            hspace=0.55, wspace=0.38,
                            left=0.17, right=0.97, top=0.93, bottom=0.03)

    for row, (m, xlabel, subtitle) in enumerate(METRICS):
        # ── Panel A: by clade ──
        ax = fig.add_subplot(gs[row, 0])
        ax.set_facecolor("#F6F8FB")
        phy = phy_agg[m]
        vals_a = []; errs_a = []; cols_a = []; labs_a = []
        for p in present_phyla:
            sub = phy[phy["phylum"]==p]
            if sub.empty: continue
            vals_a.append(float(sub["mean"].iloc[0]))
            errs_a.append(float(sub["se"].iloc[0]))
            cols_a.append(CLADE_COLORS.get(p, "#888888"))
            labs_a.append(f"{p}  (n={int(sub['n'].iloc[0])} sp)")
        yy = np.arange(len(vals_a))
        ax.barh(yy, vals_a, xerr=errs_a, color=cols_a,
                edgecolor="black", linewidth=0.4,
                error_kw=dict(ecolor="#333333", capsize=3, lw=0.8))
        ax.axvline(0, color="black", lw=1, ls="--")
        ax.set_yticks(yy); ax.set_yticklabels(labs_a, fontsize=8)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.grid(axis="x", alpha=0.25, lw=0.6)
        if row == 0:
            ax.set_title("A   By fungal clade  (mean ± SE across species)",
                         loc="left", fontsize=10, fontweight="bold")
        ax.text(0.98, 0.02, subtitle, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#9CB4CC", alpha=0.9))

        # ── Panel B: by TE class + PCG ──
        ax2 = fig.add_subplot(gs[row, 1])
        ax2.set_facecolor("#F6F8FB")
        sf = sf_agg[m]
        vals_b = []; errs_b = []; cols_b = []; labs_b = []
        for s in present_sf:
            sub = sf[sf["sf_class"]==s]
            if sub.empty: continue
            vals_b.append(float(sub["mean"].iloc[0]))
            errs_b.append(float(sub["se"].iloc[0]))
            cols_b.append(SF_COLORS.get(s, "#888888"))
            n_win = int(sub["n_win_total"].iloc[0]); n_sp = int(sub["n"].iloc[0])
            labs_b.append(f"{SF_DISPLAY.get(s, s)}  ({n_win:,} win, {n_sp} sp)")
        yy2 = np.arange(len(vals_b))
        ax2.barh(yy2, vals_b, xerr=errs_b, color=cols_b,
                 edgecolor="black", linewidth=0.4,
                 error_kw=dict(ecolor="#333333", capsize=3, lw=0.8))
        ax2.axvline(0, color="black", lw=1, ls="--")
        ax2.set_yticks(yy2); ax2.set_yticklabels(labs_b, fontsize=8.5)
        ax2.set_xlabel(xlabel + "  (mean ± SE across species)", fontsize=9)
        ax2.grid(axis="x", alpha=0.25, lw=0.6)
        if row == 0:
            ax2.set_title("B   By TE superfamily + PCG  (mean ± SE across species)",
                          loc="left", fontsize=10, fontweight="bold")

    fig.suptitle(
        f"NuclPred-predicted chromatin organisation at TE insertions vs background\n"
        f"4 chromatin metrics — 10 kb windows — {df_all['species'].nunique()} fungal species",
        fontsize=13, fontweight="bold", y=0.975)

    for ext in ("pdf","png"):
        kw = {"bbox_inches":"tight"} if ext=="pdf" else {"dpi":200,"bbox_inches":"tight"}
        fig.savefig(os.path.join(OUT, f"figG1_allfungi_multimetric.{ext}"), **kw)
    plt.close(fig)
    print("  → figG1_allfungi_multimetric saved")


if __name__ == "__main__":
    df = run()
    make_figure(df)
    print(f"\n[DONE] outputs in {OUT}")
