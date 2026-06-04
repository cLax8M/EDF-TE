#!/usr/bin/env python3
"""
chromatin_by_class_multimetric.py
==================================
Two publication-quality figures — horizontal bar layout (mean ± SE):

  FIGURE 1 — All fungi (205+ species, NuclPred-based)
    Rows: Δ NuclPred score | Δ mean occupancy | Δ CV (rugosity)
    Left  panel: by fungal clade
    Right panel: by dominant TE superfamily + PCG

  FIGURE 2 — In vivo validation (4 species, real chromatin)
    Rows: Δ NRL (bp) | Δ CV | Δ phasing (ac_val) | Δ mean signal
    Left  panel: by species
    Right panel: by TE superfamily + PCG

In both cases windows are classified by dominant occupancy class:
  BG (background) = no TE ≥ MIN_TE_FRAC and no PCG ≥ MIN_GENE_FRAC
  LTR / LINE / DNA / Helitron / PLE / SINE / Unknown = dominant TE class
  PCG = gene coverage ≥ MIN_GENE_FRAC and NOT TE-dominant

Δ metric = mean/median(class windows) − mean/median(BG windows) per species,
then aggregated (mean ± SE) across species within clade / superfamily.

Source files:
  All fungi  : _states.tsv.gz  (500 bp windows, NuclPred-based)
  In vivo    : WIG → window_metrics at 10 kb

Outputs (all in invivo_landscape/ dir):
  figF1_allfungi_multimetric.{pdf,png}
  figF2_invivo_multimetric.{pdf,png}
  figF_summary_allfungi.tsv
  figF_summary_invivo.tsv
"""
import os, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics
from clade_palette import CLADE_COLORS as _CP, CLADE_NAMES, CLADE_ORDER as _CO, TE_CLASS_COLORS, get_clade

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

WIN_ALLFUNGI = 500
WIN_INVIVO   = 10_000
MIN_TE_FRAC  = 0.20
MIN_GENE_FRAC= 0.30
MIN_WINDOWS  = 10
N_WORKERS    = 8

# ---------------------------------------------------------------- taxonomy ---
# Use canonical clade codes (2-letter) → names and colors from clade_palette
PHYLUM       = CLADE_NAMES
PHYLUM_ORDER = [CLADE_NAMES[c] for c in _CO if c in CLADE_NAMES]
CLADE_COLORS = {CLADE_NAMES[c]: _CP[c] for c in _CO if c in CLADE_NAMES}

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
    if c.startswith("LTR"):                                    return "LTR"
    if c.startswith("LINE"):                                   return "LINE"
    if c.startswith("DNA"):                                    return "DNA"
    if c.startswith("RC") or c.lower().startswith("helitron"): return "Helitron"
    if c.startswith("PLE") or c.startswith("Penelope"):        return "PLE"
    if c.startswith("SINE"):                                   return "SINE"
    if c in ("Unknown","Retroposon","Unspecified","Other"):     return "Unknown"
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
                sf_idx = SF_CLASSES.index(sf) + 1
                n_win = cov[c].shape[0]
                e = min(e, n_win * win)
                i0 = max(0, s // win); i1 = min((e - 1) // win, n_win - 1)
                if i0 > i1: continue
                for i in range(i0, i1 + 1):
                    ws = i * win; we = ws + win
                    ov = min(e, we) - max(s, ws)
                    if ov > 0: cov[c][i, sf_idx] += ov / win

    # Gene (PCG) coverage
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


# ========================================================== ALL FUNGI ====

def process_species_allfungi(args):
    sp, ts_gz, bed_p, gff_p = args
    nsf = len(SF_CLASSES); PCG_IDX = nsf + 1
    try:
        df = pd.read_csv(ts_gz, sep="\t",
                         usecols=["chrom","start","mean_occ","cv","score","state"],
                         compression="gzip")
        df = df[df["state"].notna() & (df["state"] != "NoSig")]
        df = df[df["score"].notna() & df["mean_occ"].notna() & df["cv"].notna()]
        if df.empty: return sp, None

        chrom_windows = {}
        for c, sub in df.groupby("chrom", sort=False):
            chrom_windows[c] = int(sub["start"].max() // WIN_ALLFUNGI) + 1

        dom = build_dominant_class_vec(bed_p, gff_p, chrom_windows, WIN_ALLFUNGI)

        # Build per-window arrays for each metric
        metric_arrs = {m: {} for m in ("score","mean_occ","cv")}
        for c, sub in df.groupby("chrom", sort=False):
            n = chrom_windows[c]
            for m in ("score","mean_occ","cv"):
                arr = np.full(n, np.nan, dtype=np.float32)
                idx = (sub["start"].values // WIN_ALLFUNGI).astype(int)
                mask = idx < n
                arr[idx[mask]] = sub[m].values[mask].astype(np.float32)
                metric_arrs[m][c] = arr

        # Accumulate per class
        class_vals = {m: {cls: [] for cls in SF_CLASSES + ["PCG","BG"]}
                      for m in ("score","mean_occ","cv")}
        for c, darr in dom.items():
            n = len(darr)
            for m in ("score","mean_occ","cv"):
                sarr = metric_arrs[m].get(c)
                if sarr is None: continue
                nn = min(n, len(sarr))
                for i in range(nn):
                    v = sarr[i]
                    if np.isnan(v): continue
                    d = darr[i]
                    if d == 0:
                        class_vals[m]["BG"].append(v)
                    elif 1 <= d <= nsf:
                        class_vals[m][SF_CLASSES[d-1]].append(v)
                    elif d == PCG_IDX:
                        class_vals[m]["PCG"].append(v)

        # Compute backgrounds
        bg = {}
        for m in ("score","mean_occ","cv"):
            bg_vals = class_vals[m]["BG"]
            bg[m] = float(np.mean(bg_vals)) if len(bg_vals) >= MIN_WINDOWS else np.nan
        if any(np.isnan(v) for v in bg.values()): return sp, None

        rows = []
        for cls in SF_CLASSES + ["PCG"]:
            row = {"species": sp, "sf_class": cls}
            enough = True
            for m in ("score","mean_occ","cv"):
                vals = np.array(class_vals[m][cls], dtype=np.float32)
                if len(vals) < MIN_WINDOWS: enough = False; break
                row[f"mean_{m}"] = float(np.mean(vals))
                row[f"bg_{m}"]   = float(bg[m])
                row[f"delta_{m}"] = float(np.mean(vals)) - float(bg[m])
                row[f"n_{m}"]    = len(vals)
            if not enough: continue
            rows.append(row)

        return sp, pd.DataFrame(rows) if rows else None
    except Exception as ex:
        print(f"[WARN] {sp}: {ex}", file=sys.stderr)
        return sp, None


def run_all_fungi():
    print("\n" + "="*60)
    print("FIGURE 1: All fungi — 3 metrics × (clade + superfamily)")
    print("="*60)

    tsv_files = sorted(f for f in os.listdir(CS_DIR) if f.endswith("_states.tsv.gz"))
    args_list = [(fn.replace("_states.tsv.gz",""),
                  os.path.join(CS_DIR, fn),
                  te_bed_path(fn.replace("_states.tsv.gz","")),
                  gff_path(fn.replace("_states.tsv.gz","")))
                 for fn in tsv_files]

    print(f"  Processing {len(args_list)} species …")
    all_rows = []; done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(process_species_allfungi, a): a[0] for a in args_list}
        for fut in as_completed(futs):
            sp, df = fut.result(); done += 1
            if df is not None: all_rows.append(df)
            if done % 30 == 0: print(f"    … {done}/{len(args_list)}")

    if not all_rows: print("No data!"); return None
    df_all = pd.concat(all_rows, ignore_index=True)
    df_all["phylum_code"] = df_all["species"].str.split("_",n=1).str[0]
    df_all["phylum"]      = df_all["phylum_code"].map(PHYLUM).fillna(df_all["phylum_code"])
    print(f"  {df_all['species'].nunique()} species with data")
    df_all.to_csv(os.path.join(OUT,"figF_summary_allfungi.tsv"), sep="\t", index=False)
    return df_all


def make_figure1(df_all):
    if df_all is None or df_all.empty: return

    METRICS = [
        ("delta_score",    "Δ NuclPred score",      "Compacted vs open\n(higher = more nucleosome-rich)"),
        ("delta_mean_occ", "Δ Mean occupancy",       "Nucleosome occupancy\n(higher = more packed)"),
        ("delta_cv",       "Δ CV (rugosity)",        "Chromatin irregularity\n(higher = more heterogeneous)"),
    ]

    # Pre-aggregate per-species
    sp_records = {}
    for m, _, _ in METRICS:
        sp_pool = (df_all[df_all["sf_class"].isin(SF_CLASSES)]
                   .groupby(["species","phylum"])
                   .apply(lambda d: np.average(d[m], weights=d["n_score"].clip(1)))
                   .reset_index().rename(columns={0: m}))
        sp_records[m] = sp_pool

    # Clade aggregation
    phy_agg = {}
    for m, _, _ in METRICS:
        phy = (sp_records[m].groupby("phylum")
               .agg(mean=(m,"mean"), sd=(m,"std"), n=(m,"count"))
               .reset_index())
        phy["se"] = phy["sd"] / np.sqrt(phy["n"].clip(lower=1))
        phy_agg[m] = phy

    # Common clade order: sort by mean score delta
    base_phy = phy_agg["delta_score"].copy()
    common_phyla = [p for p in PHYLUM_ORDER if p in base_phy["phylum"].values] + \
                   [p for p in base_phy["phylum"].values if p not in PHYLUM_ORDER]
    # Filter to present
    present_phyla = [p for p in common_phyla if p in base_phy["phylum"].values]
    # Sort by score delta for final order
    order_map = {p: float(base_phy.loc[base_phy["phylum"]==p,"mean"].values[0])
                 if p in base_phy["phylum"].values else 0 for p in present_phyla}
    present_phyla = sorted(present_phyla, key=lambda p: order_map[p])

    # Superfamily aggregation
    sf_agg = {}
    for m, _, _ in METRICS:
        cls_per_sp = (df_all.groupby(["sf_class","species"])
                     .apply(lambda d: np.average(d[m], weights=d["n_score"].clip(1)))
                     .reset_index().rename(columns={0: m}))
        cls = (cls_per_sp.groupby("sf_class")
               .agg(mean=(m,"mean"), sd=(m,"std"), n=(m,"count"))
               .reset_index())
        cls["se"] = cls["sd"] / np.sqrt(cls["n"].clip(lower=1))
        sf_agg[m] = cls

    # SF order: sort by score delta
    base_sf = sf_agg["delta_score"]
    present_sf = [s for s in SF_ORDER if s in base_sf["sf_class"].values]
    sf_order_map = {s: float(base_sf.loc[base_sf["sf_class"]==s,"mean"].values[0])
                    if s in base_sf["sf_class"].values else 0 for s in present_sf}
    present_sf = sorted(present_sf, key=lambda s: sf_order_map[s])
    present_sf = [s for s in present_sf if s != "PLE"]   # remove Penelope

    # ---- FIGURE ----
    n_metrics  = len(METRICS)
    n_phyla    = len(present_phyla)
    n_sf       = len(present_sf)
    row_h      = max(0.45 * max(n_phyla, n_sf), 2.5)

    fig = plt.figure(figsize=(18, row_h * n_metrics + 1.5))
    fig.patch.set_facecolor("#F6F8FB")
    gs  = gridspec.GridSpec(n_metrics, 2, figure=fig,
                            hspace=0.45, wspace=0.35,
                            left=0.13, right=0.97, top=0.93, bottom=0.04)

    palette = list(plt.cm.Set2.colors) + list(plt.cm.Set3.colors)

    for row, (m, xlabel, subtitle) in enumerate(METRICS):
        # Panel A — by clade
        ax = fig.add_subplot(gs[row, 0])
        ax.set_facecolor("#F6F8FB")
        phy = phy_agg[m]
        vals = []; errs = []; colors_a = []; labels_a = []
        for p in present_phyla:
            sub = phy[phy["phylum"]==p]
            if sub.empty: continue
            vals.append(float(sub["mean"].iloc[0])); errs.append(float(sub["se"].iloc[0]))
            colors_a.append(CLADE_COLORS.get(p, "#888888"))
            n = int(sub["n"].iloc[0])
            labels_a.append(f"{p}  (n={n})")
        yy = np.arange(len(vals))
        ax.barh(yy, vals, xerr=errs, color=colors_a,
                edgecolor="black", linewidth=0.4,
                error_kw=dict(ecolor="#333333", capsize=3, lw=0.8))
        ax.axvline(0, color="black", lw=1, ls="--")
        ax.set_yticks(yy); ax.set_yticklabels(labels_a, fontsize=8.5)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.grid(axis="x", alpha=0.25, lw=0.6)
        if row == 0:
            ax.set_title("A   By fungal clade", loc="left", fontsize=10, fontweight="bold")
        ax.text(0.98, 0.02, subtitle, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#9CB4CC", alpha=0.9))

        # Panel B — by superfamily + PCG
        ax2 = fig.add_subplot(gs[row, 1])
        ax2.set_facecolor("#F6F8FB")
        sf_data = sf_agg[m]
        vals2 = []; errs2 = []; colors2 = []; labels2 = []
        for s in present_sf:
            sub = sf_data[sf_data["sf_class"]==s]
            if sub.empty: continue
            vals2.append(float(sub["mean"].iloc[0])); errs2.append(float(sub["se"].iloc[0]))
            colors2.append(SF_COLORS.get(s, "#888888"))
            labels2.append(f"{'Rolling Circle' if s == 'Helitron' else s}  (n={int(sub['n'].iloc[0])})")
        yy2 = np.arange(len(vals2))
        ax2.barh(yy2, vals2, xerr=errs2, color=colors2,
                 edgecolor="black", linewidth=0.4,
                 error_kw=dict(ecolor="#333333", capsize=3, lw=0.8))
        ax2.axvline(0, color="black", lw=1, ls="--")
        ax2.set_yticks(yy2); ax2.set_yticklabels(labels2, fontsize=9)
        ax2.set_xlabel(xlabel + "  (averaged across species)", fontsize=9)
        ax2.grid(axis="x", alpha=0.25, lw=0.6)
        if row == 0:
            ax2.set_title("B   By TE superfamily + PCG  (mean ± SE across species)",
                          loc="left", fontsize=10, fontweight="bold")

    fig.suptitle(f"Chromatin organisation at TE insertions vs genome background\n"
                 f"NuclPred-based scores (500 bp windows) — {df_all['species'].nunique()} fungal species",
                 fontsize=13, fontweight="bold", y=0.975)

    for ext in ("pdf","png"):
        kw = {"bbox_inches":"tight"} if ext=="pdf" else {"dpi":200,"bbox_inches":"tight"}
        fig.savefig(os.path.join(OUT, f"figF1_allfungi_multimetric.{ext}"), **kw)
    plt.close(fig)
    print("  → figF1_allfungi_multimetric saved")


# ========================================================== IN VIVO ======

INVIVO_METRICS = [
    ("ac_lag",   "Δ NRL (bp)",         "Nucleosome repeat length\n(Δ < 0 = nucleosomes more closely packed)"),
    ("cv",       "Δ CV (rugosity)",     "Chromatin irregularity\n(Δ > 0 = more heterogeneous)"),
    ("ac_val",   "Δ Phasing (AC peak)","Nucleosome array regularity\n(Δ > 0 = stronger periodic phasing)"),
    ("mean_sig", "Δ Mean signal",       "Chromatin accessibility/occupancy\n(MNase: Δ < 0 = less accessible)"),
]

SP_COLORS_IV = {"MU":"#4CAF50","AS":"#E74C3C","BA":"#3498DB","BL":"#9B59B6"}

def run_invivo():
    print("\n" + "="*60)
    print("FIGURE 2: In vivo — 4 metrics × (species + superfamily) @ 10 kb")
    print("="*60)

    all_rows = []

    for entry in IN_VIVO_PANEL:
        lab = entry["label"]; sp = entry["sp"]; wig = entry["wig"]
        if not os.path.exists(wig): print(f"[SKIP] {lab}"); continue
        print(f"\n[{lab}]")

        sz  = chrom_sizes_fai(sp)
        if sz is None:
            from chromatin_states_from_nuclpred import load_nuclpred_wig as _lw
            _tmp = _lw(wig, None)
            sz = {c: len(a) for c, a in _tmp.items()}
        arrs = load_nuclpred_wig(wig, sz)
        bed  = te_bed_path(sp)
        gff  = gff_path(sp)

        chrom_windows = {c: L // WIN_INVIVO for c, L in sz.items() if L >= WIN_INVIVO}
        dom = build_dominant_class_vec(bed, gff, chrom_windows, WIN_INVIVO)

        rows_win = []
        for c, occ in arrs.items():
            n = len(occ) // WIN_INVIVO
            if n == 0: continue
            means, cvs, aclags, acvals, _ = window_metrics(occ, WIN_INVIVO)
            darr = dom.get(c, np.zeros(min(n, chrom_windows.get(c,0)), dtype=np.int8))
            n = min(n, len(darr))
            rows_win.append(pd.DataFrame({
                "chrom": c, "win_idx": np.arange(n),
                "ac_lag": aclags[:n], "mean_sig": means[:n],
                "cv": cvs[:n], "ac_val": acvals[:n],
                "dom_class": darr[:n],
            }))
        if not rows_win: continue

        df = pd.concat(rows_win, ignore_index=True)
        df = df[df["mean_sig"] > 1e-3].copy()

        nsf = len(SF_CLASSES); PCG_IDX = nsf + 1

        # Background medians
        bg = {}
        for m in ("ac_lag","cv","ac_val","mean_sig"):
            vals = df.loc[df["dom_class"]==0, m].dropna().values
            bg[m] = float(np.median(vals)) if len(vals) >= MIN_WINDOWS else np.nan

        if any(np.isnan(v) for v in bg.values()):
            print(f"  [WARN] insufficient BG windows"); continue

        print(f"  BG: NRL={bg['ac_lag']:.1f}  CV={bg['cv']:.3f}  "
              f"AC={bg['ac_val']:.3f}  sig={bg['mean_sig']:.1f}")

        for ci, cls in enumerate(SF_CLASSES + ["PCG"], start=1):
            idx = ci if cls != "PCG" else PCG_IDX
            row = {"species": lab, "assay": entry["assay"],
                   "phylum_code": sp.split("_")[0],
                   "sf_class": cls}
            ok = True
            for m in ("ac_lag","cv","ac_val","mean_sig"):
                vals = df.loc[df["dom_class"]==idx, m].dropna().values
                if len(vals) < MIN_WINDOWS: ok = False; break
                row[f"median_{m}"] = float(np.median(vals))
                row[f"bg_{m}"]     = bg[m]
                row[f"delta_{m}"]  = float(np.median(vals)) - bg[m]
                row[f"n_{m}"]      = len(vals)
            if not ok: continue
            all_rows.append(row)
            print(f"    {cls:12s}: NRL={row['median_ac_lag']:.1f}  "
                  f"ΔNRL={row['delta_ac_lag']:+.1f}  "
                  f"ΔCV={row['delta_cv']:+.3f}  n={row['n_ac_lag']}")

    if not all_rows: print("No data!"); return None
    df_iv = pd.DataFrame(all_rows)
    df_iv.to_csv(os.path.join(OUT,"figF_summary_invivo.tsv"), sep="\t", index=False)
    return df_iv


def _lollipop_panel(ax, items, vals_map, col_map, metric,
                    se_map=None, y_top_to_bot=True,
                    left_labels=True, pcg_shade=False, display=None,
                    points=None, point_colors=None, special_marks=None):
    """
    Draw a single lollipop sub-panel.
    items         : ordered list of row labels
    vals_map      : {item: float}  — mean / representative value
    col_map       : {item: color}
    se_map        : {item: float} or None — SE error bars on the mean
    display       : {internal: display label} rename dict
    points        : {item: [list of floats]} — individual observations
    point_colors  : {item: [list of hex strings]} — one color per observation
    special_marks : {item: float} — extra value drawn as a white diamond
                    (used to overlay PCG per-species on the TE-mean rows)
    """
    display = display or {}
    n   = len(items)
    rng = np.random.default_rng(42)
    yy  = np.arange(n, dtype=float)
    if y_top_to_bot:
        yy = yy[::-1]
    y_map = {item: float(yy[i]) for i, item in enumerate(items)}

    # light band for PCG row
    if pcg_shade and "PCG" in y_map:
        ypcg = y_map["PCG"]
        ax.axhspan(ypcg - 0.45, ypcg + 0.45, alpha=0.06,
                   color="#888888", zorder=0)

    ax.axvline(0, color="#888888", lw=1.0, ls="--", alpha=0.7, zorder=0)

    # ── individual points (behind lollipops) ──────────────────────────
    if points:
        for item in items:
            pts_raw = points.get(item, [])
            pts = [float(v) for v in pts_raw
                   if v is not None and not np.isnan(float(v))]
            if not pts:
                continue
            y   = y_map[item]
            # per-point colors or uniform
            if point_colors and item in point_colors:
                pcols = list(point_colors[item])[:len(pts)]
                if len(pcols) < len(pts):
                    pcols += [col_map.get(item, "#555555")] * (len(pts) - len(pcols))
            else:
                pcols = [col_map.get(item, "#555555")] * len(pts)
            jitter = rng.uniform(-0.30, 0.30, size=len(pts))
            for xi, (px, pc) in enumerate(zip(pts, pcols)):
                ax.scatter([px], [y + jitter[xi]],
                           color=pc, s=44, alpha=0.72,
                           edgecolors="white", linewidths=0.5, zorder=1)

    # ── lollipop stems + heads ────────────────────────────────────────
    for item in items:
        val = vals_map.get(item)
        if val is None or np.isnan(float(val)):
            continue
        val = float(val)
        y   = y_map[item]
        col = col_map.get(item, "#555555")
        se  = float(se_map[item]) if (se_map and item in se_map) else 0.0

        # stem from 0 to dot
        ax.plot([0, val], [y, y], color="#bbbbbb", lw=2.0,
                zorder=2, solid_capstyle="round")
        # SE whisker
        if se > 0:
            ax.plot([val - se, val + se], [y, y], color=col,
                    lw=2.0, zorder=3, solid_capstyle="butt")
            for sign in (-1, 1):
                ax.plot([val + sign * se, val + sign * se],
                        [y - 0.18, y + 0.18], color=col, lw=1.5, zorder=3)
        # dot
        ax.scatter([val], [y], color=col, s=160, zorder=4,
                   edgecolors="white", linewidths=1.5)

    # ── special marks: white diamond (PCG overlay) ───────────────────────
    if special_marks:
        for item in items:
            sv = special_marks.get(item)
            if sv is None:
                continue
            try:
                sv = float(sv)
            except (TypeError, ValueError):
                continue
            if np.isnan(sv):
                continue
            y   = y_map[item]
            col = col_map.get(item, "#555555")
            ax.scatter([sv], [y], marker="D", s=110, zorder=5,
                       facecolors="white", edgecolors=col, linewidths=2.2)

    ax.set_yticks(list(y_map.values()))
    if left_labels:
        ax.set_yticklabels([display.get(it, it) for it in items], fontsize=10.5)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)

    ax.set_ylim(-0.65, n - 0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", labelsize=8.5)
    ax.grid(axis="x", alpha=0.12, lw=0.5, color="gray")


def make_figure2(df_iv, with_points=False):
    """
    Main figure — in vivo chromatin at TE insertions.
    2 rows (A = by species, B = by TE class) × 3 cols (Δ NRL | Δ CV | Δ Phasing).
    Lollipop / Cleveland dot style, white background.
    with_points : overlay individual data points (per-TE-class on species rows;
                  per-species on TE-class rows)
    """
    if df_iv is None or df_iv.empty:
        return

    MKEYS   = [
        ("delta_ac_lag", "Δ NRL (bp)",
         "← more closely packed        further apart →"),
        ("delta_cv",     "Δ CV (rugosity)",
         "← more uniform        more heterogeneous →"),
        ("delta_ac_val", "Δ Phasing",
         "← less periodic        more regular arrays →"),
    ]
    DISPLAY = {"Helitron": "Rolling Circle"}
    TE_SHOW = [c for c in SF_CLASSES if c != "PLE"] + ["PCG"]
    n_m     = len(MKEYS)

    # ── Panel A: per-species TE-pooled ─────────────────────────────────
    sp_vals = {}; sp_cols = {}; sp_labs = {}
    for sp_lab in df_iv["species"].unique():
        sub = df_iv[(df_iv["species"] == sp_lab) & df_iv["sf_class"].isin(SF_CLASSES)]
        if sub.empty:
            continue
        pc    = sub["phylum_code"].iloc[0]
        assay = sub["assay"].iloc[0]
        ww    = sub["n_ac_lag"].clip(lower=1) if "n_ac_lag" in sub.columns else None
        entry = {}
        for m, _, _ in MKEYS:
            if m not in sub.columns:
                continue
            vals = sub[m].dropna()
            if vals.empty:
                continue
            entry[m] = float(np.average(vals, weights=ww[vals.index]) if ww is not None
                             else vals.mean())
        if entry:
            sp_vals[sp_lab] = entry
            sp_cols[sp_lab] = _CP.get(pc, "#888888")
            sp_labs[sp_lab] = f"{sp_lab}  ({assay})"

    sp_list = sorted(sp_vals, key=lambda s: sp_vals[s].get("delta_ac_lag", 0))

    # ── Panel B: per-TE-class (mean ± SE across species) ─────────────
    te_vals = {}; te_se = {}; te_cols = {}
    for cls in TE_SHOW:
        sub = df_iv[df_iv["sf_class"] == cls]
        if sub.empty:
            continue
        entry_v = {}; entry_se = {}
        for m, _, _ in MKEYS:
            if m not in sub.columns:
                continue
            vals = sub[m].dropna().values
            if len(vals) == 0:
                continue
            entry_v[m]  = float(np.mean(vals))
            entry_se[m] = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) \
                          if len(vals) > 1 else 0.0
        if entry_v:
            te_vals[cls] = entry_v
            te_se[cls]   = entry_se
            te_cols[cls] = SF_COLORS.get(cls, "#333333")

    te_list = [c for c in TE_SHOW if c in te_vals]
    te_list.sort(key=lambda c: te_vals[c].get("delta_ac_lag", 0))

    # ── PCG value per species (for diamond overlay in Panel A) ───────────
    sp_pcg = {}
    for sp_lab in sp_list:
        sub_pcg = df_iv[(df_iv["species"] == sp_lab) & (df_iv["sf_class"] == "PCG")]
        sp_pcg[sp_lab] = {}
        for m, _, _ in MKEYS:
            if m in sub_pcg.columns:
                vals = sub_pcg[m].dropna()
                if not vals.empty:
                    sp_pcg[sp_lab][m] = float(vals.iloc[0])

    # ── Individual points for overlay ────────────────────────────────────
    # Panel A — per-species rows: dots colored by TE class
    sp_pts_vals = {}; sp_pts_cols = {}
    for sp_lab in sp_list:
        sub = df_iv[(df_iv["species"] == sp_lab) & df_iv["sf_class"].isin(SF_CLASSES)]
        sp_pts_vals[sp_lab] = {}; sp_pts_cols[sp_lab] = {}
        for m, _, _ in MKEYS:
            if m not in sub.columns:
                continue
            rows_ok = sub[sub[m].notna()]
            sp_pts_vals[sp_lab][m] = rows_ok[m].tolist()
            sp_pts_cols[sp_lab][m] = [SF_COLORS.get(r["sf_class"], "#888888")
                                       for _, r in rows_ok.iterrows()]

    # Panel B — per-TE-class rows: dots colored by species
    te_pts_vals = {}; te_pts_cols = {}
    for cls in te_list:
        sub = df_iv[df_iv["sf_class"] == cls]
        te_pts_vals[cls] = {}; te_pts_cols[cls] = {}
        for m, _, _ in MKEYS:
            if m not in sub.columns:
                continue
            rows_ok = sub[sub[m].notna()]
            te_pts_vals[cls][m] = rows_ok[m].tolist()
            te_pts_cols[cls][m] = [sp_cols.get(r["species"], "#888888")
                                    for _, r in rows_ok.iterrows()]

    # ── Layout ───────────────────────────────────────────────────────────
    n_sp = len(sp_list)
    n_te = len(te_list)
    fig  = plt.figure(figsize=(13, 3.2 + 0.72 * (n_sp + n_te)))
    fig.patch.set_facecolor("white")
    outer = gridspec.GridSpec(
        2, 1, figure=fig,
        height_ratios=[n_sp, n_te],
        hspace=0.62,
        left=0.02, right=0.97, top=0.91, bottom=0.07)

    for pi, (panel_lbl, items, v_map, se_map, c_map, disp_map, pts_v, pts_c) in enumerate([
        ("A   In vivo — by species",
         sp_list,
         {s: {m: sp_vals[s].get(m) for m, _, _ in MKEYS} for s in sp_list},
         None,
         sp_cols,
         sp_labs,
         sp_pts_vals if with_points else None,
         sp_pts_cols if with_points else None),
        ("B   In vivo — by TE class  (mean ± SE across species)",
         te_list,
         {c: {m: te_vals[c].get(m) for m, _, _ in MKEYS} for c in te_list},
         {c: {m: te_se[c].get(m, 0) for m, _, _ in MKEYS} for c in te_list},
         te_cols,
         {c: DISPLAY.get(c, c) for c in te_list},
         te_pts_vals if with_points else None,
         te_pts_cols if with_points else None),
    ]):
        inner = gridspec.GridSpecFromSubplotSpec(
            1, n_m, subplot_spec=outer[pi], wspace=0.09)

        for mi, (m, xlabel, note) in enumerate(MKEYS):
            ax = fig.add_subplot(inner[mi])
            ax.set_facecolor("white")

            _lollipop_panel(
                ax, items,
                vals_map      = {it: v_map[it].get(m) for it in items},
                col_map       = c_map,
                metric        = m,
                se_map        = {it: se_map[it].get(m, 0) for it in items}
                                if se_map else None,
                left_labels   = (mi == 0),
                pcg_shade     = (pi == 1),
                display       = disp_map,
                points        = {it: pts_v[it].get(m, []) for it in items}
                                if pts_v else None,
                point_colors  = {it: pts_c[it].get(m, []) for it in items}
                                if pts_c else None,
                special_marks = {it: sp_pcg[it].get(m) for it in items}
                                if pi == 0 else None,
            )
            ax.set_xlabel(xlabel, fontsize=10, fontweight="bold", labelpad=6)
            ax.text(0.5, -0.18, note, transform=ax.transAxes,
                    ha="center", va="top", fontsize=7, color="#888888",
                    style="italic")
            if mi == 0:
                lbl_suffix = "  (◆ = PCG)" if pi == 0 else ""
                ax.set_title(panel_lbl + lbl_suffix, loc="left", fontsize=11,
                             fontweight="bold", pad=10)

    suffix = "_withpoints" if with_points else ""
    fig.suptitle(
        "Real chromatin organisation at TE insertions vs. genome background\n"
        "MNase / MAINE / ATAC-seq — 10 kb windows"
        + ("  •  dots = individual TE classes / species" if with_points else ""),
        fontsize=12, fontweight="bold", y=0.98)

    for ext in ("pdf", "png"):
        kw = {"bbox_inches": "tight"} if ext == "pdf" \
             else {"dpi": 200, "bbox_inches": "tight"}
        fig.savefig(os.path.join(OUT, f"figF2_invivo_multimetric{suffix}.{ext}"), **kw)
    plt.close(fig)
    print(f"  → figF2_invivo_multimetric{suffix} saved")


def make_figS_nuclpred(df_insil, with_points=False):
    """
    Supplementary — NuclPred in silico chromatin prediction at TEs.
    Same lollipop layout as make_figure2 but:
      Panel A: by fungal clade  (195 species → mean ± SE per clade)
      Panel B: by TE class      (mean ± SE across all species)
    3 metrics: Δ NRL | Δ CV | Δ Phasing
    with_points: overlay individual species dots
    """
    if df_insil is None or df_insil.empty:
        return

    MKEYS   = [
        ("delta_ac_lag", "Δ NRL (bp)",
         "← more closely packed        further apart →"),
        ("delta_cv",     "Δ CV (rugosity)",
         "← more uniform        more heterogeneous →"),
        ("delta_ac_val", "Δ Phasing",
         "← less periodic        more regular arrays →"),
    ]
    DISPLAY_SF = {"Helitron": "Rolling Circle"}
    TE_SHOW = [c for c in SF_CLASSES if c != "PLE"] + ["PCG"]
    n_m     = len(MKEYS)

    # ── Panel A: by clade ─────────────────────────────────────────────
    #   per-species TE-pooled, then aggregate per phylum
    sp_pool_all = {}   # {metric: DataFrame(species, phylum, val)}
    phy_rows = []
    for m, _, _ in MKEYS:
        if m not in df_insil.columns:
            continue
        sp_pool = (df_insil[df_insil["sf_class"].isin(SF_CLASSES)]
                   .groupby(["species", "phylum"])
                   .apply(lambda d: np.average(d[m], weights=d["n_windows"].clip(1)))
                   .reset_index().rename(columns={0: "val"}))
        sp_pool_all[m] = sp_pool
        phy = (sp_pool.groupby("phylum")
               .agg(mean=("val", "mean"), sd=("val", "std"), n=("val", "count"))
               .reset_index())
        phy["se"] = phy["sd"] / np.sqrt(phy["n"].clip(1))
        phy["metric"] = m
        phy_rows.append(phy)
    df_phy = pd.concat(phy_rows, ignore_index=True) if phy_rows else pd.DataFrame()

    # clade order by mean delta_ac_lag
    if not df_phy.empty:
        base = df_phy[df_phy["metric"] == "delta_ac_lag"]
        order_v = {p: float(base.loc[base["phylum"] == p, "mean"].values[0])
                   for p in base["phylum"].unique()
                   if len(base.loc[base["phylum"] == p]) > 0}
        # Use PHYLUM_ORDER where possible, then sort remainder by value
        phy_present = [p for p in PHYLUM_ORDER if p in order_v] + \
                      [p for p in order_v if p not in PHYLUM_ORDER]
        phy_list = sorted(phy_present, key=lambda p: order_v.get(p, 0))
        phy_list = [p for p in phy_list
                    if p not in ("Cryptomycota", "Aphelidiomycota")]
    else:
        phy_list = []

    # ── Panel B: by TE class ──────────────────────────────────────────
    sp_per_cls_all = {}  # {metric: DataFrame(sf_class, species, val)}
    te_rows = []
    for m, _, _ in MKEYS:
        if m not in df_insil.columns:
            continue
        sp_per_cls = (df_insil.groupby(["sf_class", "species"])
                      .apply(lambda d: float(
                          np.average(d[m], weights=d["n_windows"].clip(1))))
                      .reset_index().rename(columns={0: "val"}))
        sp_per_cls_all[m] = sp_per_cls
        cls_agg = (sp_per_cls.groupby("sf_class")
                   .agg(mean=("val", "mean"), sd=("val", "std"), n=("val", "count"))
                   .reset_index())
        cls_agg["se"]     = cls_agg["sd"] / np.sqrt(cls_agg["n"].clip(1))
        cls_agg["metric"] = m
        te_rows.append(cls_agg)
    df_cls = pd.concat(te_rows, ignore_index=True) if te_rows else pd.DataFrame()

    if not df_cls.empty:
        base_cls = df_cls[df_cls["metric"] == "delta_ac_lag"]
        order_cls = {c: float(base_cls.loc[base_cls["sf_class"] == c, "mean"].values[0])
                     for c in base_cls["sf_class"].unique()
                     if len(base_cls.loc[base_cls["sf_class"] == c]) > 0}
        te_list = [c for c in TE_SHOW if c in order_cls and c != "PLE"]
        te_list.sort(key=lambda c: order_cls.get(c, 0))
    else:
        te_list = []

    # ── canonical clade colors from clade_palette ────────────────────
    # TSV phylum names may differ from CLADE_NAMES values; cover both
    name_to_code = {v: k for k, v in CLADE_NAMES.items()}
    # aliases: TSV name → 2-letter code
    _ALIASES = {
        "Chytridiomycota":      "CH",
        "Entomophthoromycota":  "EN",
        "Glomeromycota":        "GL",
        "Kickxellomycota":      "KX",
        "Mortierellomycota":    "MO",
        "Mucoromycota":         "MU",
        "Olpidiomycota":        "OL",
        "Zoopagomycota":        "ZO",
        "Nuclear-other":        "NU",
        "NE":                   "NE",
        # codes passed through directly
        **{c: c for c in _CP},
    }
    # Pretty display names for phylum labels
    _PHY_DISPLAY = {
        "NE":           "Neocallimastigomycetes",
        "Nuclear-other":"Nucleariida + other",
        **{p: p for p in phy_list},   # default: keep as-is
    }
    def _phy_color(p):
        code = name_to_code.get(p) or _ALIASES.get(p, "")
        return _CP.get(code, "#888888")
    phy_col_map = {p: _phy_color(p) for p in phy_list}
    te_col_map    = {c: SF_COLORS.get(c, "#333333") for c in te_list}

    # ── Layout ───────────────────────────────────────────────────────────
    n_phy = len(phy_list)
    n_te  = len(te_list)
    fig   = plt.figure(figsize=(13, 3.2 + 0.58 * (n_phy + n_te)))
    fig.patch.set_facecolor("white")
    outer = gridspec.GridSpec(
        2, 1, figure=fig,
        height_ratios=[n_phy, n_te],
        hspace=0.55,
        left=0.02, right=0.97, top=0.92, bottom=0.05)

    for pi, (panel_lbl, items, col_map) in enumerate([
        ("A   NuclPred (in silico) — by fungal clade  "
         f"(mean ± SE across species, n = {df_insil['species'].nunique()})",
         phy_list, phy_col_map),
        ("B   NuclPred (in silico) — by TE class  "
         "(mean ± SE across species)",
         te_list, te_col_map),
    ]):
        inner = gridspec.GridSpecFromSubplotSpec(
            1, n_m, subplot_spec=outer[pi], wspace=0.09)

        for mi, (m, xlabel, note) in enumerate(MKEYS):
            ax = fig.add_subplot(inner[mi])
            ax.set_facecolor("white")

            # build vals / se and optional per-species points
            if pi == 0 and not df_phy.empty:
                src = df_phy[df_phy["metric"] == m]
                v_map  = {p: float(src.loc[src["phylum"] == p, "mean"].values[0])
                          for p in items if len(src.loc[src["phylum"] == p]) > 0}
                se_map = {p: float(src.loc[src["phylum"] == p, "se"].values[0])
                          for p in items if len(src.loc[src["phylum"] == p]) > 0}
                if with_points and m in sp_pool_all:
                    sp_df = sp_pool_all[m]
                    pts_map  = {p: sp_df.loc[sp_df["phylum"] == p, "val"].dropna().tolist()
                                for p in items}
                    # points colored by clade (same color as lollipop head)
                    pcol_map = {p: [_phy_color(p)] *
                                   len(sp_df.loc[sp_df["phylum"] == p, "val"].dropna())
                                for p in items}
                else:
                    pts_map = None; pcol_map = None
            else:
                src    = df_cls[df_cls["metric"] == m]
                v_map  = {c: float(src.loc[src["sf_class"] == c, "mean"].values[0])
                          for c in items if len(src.loc[src["sf_class"] == c]) > 0}
                se_map = {c: float(src.loc[src["sf_class"] == c, "se"].values[0])
                          for c in items if len(src.loc[src["sf_class"] == c]) > 0}
                if with_points and m in sp_per_cls_all:
                    sp_df = sp_per_cls_all[m]
                    # join phylum info
                    sp_ph = (df_insil[["species","phylum"]]
                             .drop_duplicates(subset="species")
                             .set_index("species")["phylum"])
                    pts_map  = {c: sp_df.loc[sp_df["sf_class"] == c, "val"].dropna().tolist()
                                for c in items}
                    pcol_map = {}
                    for c in items:
                        rows_c = sp_df.loc[sp_df["sf_class"] == c].dropna(subset=["val"])
                        pcol_map[c] = [_phy_color(sp_ph.get(s, ""))
                                       for s in rows_c["species"]]
                else:
                    pts_map = None; pcol_map = None

            disp_map = {c: DISPLAY_SF.get(c, c) for c in items} if pi == 1 \
                       else {p: _PHY_DISPLAY.get(p, p) for p in items}

            _lollipop_panel(
                ax, items,
                vals_map     = v_map,
                col_map      = col_map,
                metric       = m,
                se_map       = se_map,
                left_labels  = (mi == 0),
                pcg_shade    = (pi == 1),
                display      = disp_map,
                points       = pts_map,
                point_colors = pcol_map,
            )
            ax.set_xlabel(xlabel, fontsize=10, fontweight="bold", labelpad=6)
            ax.text(0.5, -0.18, note, transform=ax.transAxes,
                    ha="center", va="top", fontsize=7, color="#888888",
                    style="italic")
            if mi == 0:
                ax.set_title(panel_lbl, loc="left", fontsize=10,
                             fontweight="bold", pad=10)

    suffix = "_withpoints" if with_points else ""
    fig.suptitle(
        "NuclPred-predicted chromatin at TE insertions vs. genome background\n"
        "Sequence-intrinsic nucleosome preference — 10 kb windows — "
        f"{df_insil['species'].nunique()} fungal species"
        + ("  •  dots = individual species" if with_points else ""),
        fontsize=11, fontweight="bold", y=0.985)

    for ext in ("pdf", "png"):
        kw = {"bbox_inches": "tight"} if ext == "pdf" \
             else {"dpi": 200, "bbox_inches": "tight"}
        fig.savefig(os.path.join(OUT, f"figS_nuclpred_chromatin{suffix}.{ext}"), **kw)
    plt.close(fig)
    print(f"  → figS_nuclpred_chromatin{suffix} saved")


# ======================================================================= main
def main():
    df_all = run_all_fungi()
    make_figure1(df_all)

    df_iv = run_invivo()
    make_figure2(df_iv)

    insil_tsv = os.path.join(OUT, "figG1_summary_allfungi_multimetric.tsv")
    df_insil  = pd.read_csv(insil_tsv, sep="\t") if os.path.exists(insil_tsv) else None
    if df_insil is not None:
        print(f"  Loaded in silico TSV: {df_insil['species'].nunique()} species")
        make_figS_nuclpred(df_insil)
        make_figS_nuclpred(df_insil, with_points=True)

    make_figure2(df_iv, with_points=True)

    print(f"\n[DONE] outputs in {OUT}")

if __name__ == "__main__":
    main()
