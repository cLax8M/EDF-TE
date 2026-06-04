#!/usr/bin/env python3
"""
partC_superfamily_invitro.py
============================
Part C: Chromatin metrics by dominant TE superfamily — 4 in vivo species (5 kb windows)
Part D: Same analysis for 5 in vitro NucSeq species
PartCD: Direct comparison of NRL (ac_lag) and CV for in vivo vs in vitro

Windows are classified by their dominant TE superfamily (the class with most bp coverage):
    LTR, LINE, DNA, Helitron, PLE, SINE, Unknown
    "TE-poor": windows where no transposable class covers >2 % of the window

Simple_repeat, Low_complexity, Satellite, tRNA, rRNA, snRNA are excluded from
all TE classifications (not considered mobile elements).

Metrics per window (same as Part B of invivo_nucleosome_landscape.py):
    mean_sig  — mean WIG signal
    cv        — coefficient of variation (rugosity)
    ac_val    — autocorrelation peak 140–240 bp (phasing strength)
    ac_lag    — position of that peak (NRL in bp)

Outputs (all in invivo_landscape/ dir):
    partC_sf_<LABEL>.{pdf,png}          per-species violin plots (in vivo)
    partC_sf_grid.{pdf,png}             grid all in vivo species
    partD_sf_<LABEL>.{pdf,png}          per-species violin plots (in vitro)
    partD_sf_grid.{pdf,png}             grid all in vitro species
    partCD_comparison.{pdf,png}         NRL + CV: in vivo vs in vitro side-by-side
    partCD_mereb1_invivo_vs_invitro.{pdf,png}  BA_Mereb1 ATAC vs NucSeq (same species)
    partCD_summary.tsv
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kruskal, mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics

# ------------------------------------------------------------------ paths ---
MV      = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes"
MV_NUC  = f"{MV}/nuc-seq"
MV_ATC  = f"{MV}/atac-seq"
FAI_FMT = "/legserv/Temp/Carlos/EDF-TE/data/genomes/{sp}.fasta.fai"
TE_DATA = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
TE_MV   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
OUT     = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/invivo_landscape"
os.makedirs(OUT, exist_ok=True)

WIN            = 5_000       # window size
MIN_TE_FRAC    = 0.02        # minimum fraction of a class to call a window as that class
MIN_SF_WINDOWS = 20          # minimum windows in a group to include it in plots

# ----------------------------------------------------------------- panels ---
PANEL_C = [   # in vivo
    dict(label="MU_Rhimi59",     sp="MU_Rhimi59",  assay="MNase", atac=False,
         wig="/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"),
    dict(label="AS_Ncrass2",     sp="AS_Ncrass2",  assay="MNase", atac=False,
         wig=f"{MV_NUC}/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.wig"),
    dict(label="AS_Bcine1",      sp="AS_Bcine1",   assay="MAINE", atac=False,
         wig=f"{MV_NUC}/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.wig"),
    dict(label="BA_Mereb1_atac", sp="BA_Mereb1",   assay="ATAC",  atac=True,
         wig=f"{MV_ATC}/BA_Mereb1/BA_Mereb1_atac.wig"),
]

PANEL_D = [   # in vitro NucSeq
    dict(label="AS_Usnflo1", sp="AS_Usnflo1", assay="NucSeq", atac=False,
         wig=f"{MV_NUC}/AS_Usnflo1/AS_Usnflo1_danpos/pooled/AS_Usnflo1_dedup.smooth.wig"),
    dict(label="BA_Mereb1",  sp="BA_Mereb1",  assay="NucSeq", atac=False,
         wig=f"{MV_NUC}/BA_Mereb1/BA_Mereb1_danpos/pooled/BA_Mereb1_dedup.smooth.wig"),
    dict(label="BA_Mrafri1", sp="BA_Mrafri1", assay="NucSeq", atac=False,
         wig=f"{MV_NUC}/BA_Mrafri1/BA_Mrafri1_danpos/pooled/BA_Mrafri1_dedup.smooth.wig"),
    dict(label="BA_Rhoto",   sp="BA_Rhoto",   assay="NucSeq", atac=False,
         wig=f"{MV_NUC}/BA_Rhoto/BA_Rhoto_danpos/pooled/BA_Rhoto_dedup.smooth.wig"),
    dict(label="BL_Catan2",  sp="BL_Catan2",  assay="NucSeq", atac=False,
         wig=f"{MV_NUC}/BL_Catan2/BL_Catan2_danpos/pooled/BL_Catan2_dedup.smooth.wig"),
]

# --------------------------------------------------------- superfamily defs ---
SF_CLASSES = ["LTR", "LINE", "DNA", "Helitron", "PLE", "SINE", "Unknown"]
SF_ORDER   = ["TE-poor"] + SF_CLASSES
SF_COLORS  = {
    "TE-poor":   "#d0d0d0",
    "LTR":       "#e74c3c",
    "LINE":      "#e67e22",
    "DNA":       "#3498db",
    "Helitron":  "#9b59b6",
    "PLE":       "#1abc9c",
    "SINE":      "#f1c40f",
    "Unknown":   "#7f8c8d",
}

FEATS = [
    ("mean_sig", "Mean signal\n(log)"),
    ("cv",       "CV (rugosity)"),
    ("ac_val",   "AC peak\n(phasing)"),
    ("ac_lag",   "NRL (bp)"),
]

# ---------------------------------------------------------------- helpers ---
def te_bed(sp):
    a = TE_DATA.format(sp=sp); b = TE_MV.format(sp=sp)
    return a if os.path.exists(a) else b

def chrom_sizes(sp):
    fai = FAI_FMT.format(sp=sp)
    df  = pd.read_csv(fai, sep="\t", header=None, names=["c","L","o","lb","lby"])
    return dict(zip(df["c"], df["L"].astype(int)))

def classify_superfamily(col4):
    """
    Return simplified superfamily name from BED col4, or None to exclude.
    Excluded: Simple_repeat, Low_complexity, Satellite, tRNA, rRNA, snRNA, RNA.
    """
    raw = col4.strip()
    c   = raw.split("/")[0].rstrip("?").strip()
    if c.startswith("LTR"):                                  return "LTR"
    if c.startswith("LINE"):                                 return "LINE"
    if c.startswith("DNA"):                                  return "DNA"
    if c.startswith("RC") or c.lower().startswith("heitron"): return "Helitron"
    if c.startswith("PLE") or c.startswith("Penelope"):      return "PLE"
    if c.startswith("SINE"):                                 return "SINE"
    if c in ("Unknown","Retroposon","Unspecified","Other"):   return "Unknown"
    # Exclude repeats that are NOT mobile elements
    return None

def te_superfamily_windows(bed_path, sizes, win):
    """
    Build per-window TE superfamily fractions.
    Returns pd.DataFrame with columns:
        chrom, win_idx, dominant_sf, te_total_frac, te_LTR, te_LINE, te_DNA, ...
    """
    # per-base coverage arrays (one per class)
    bp_class = {c: {sf: np.zeros(L, dtype=np.float32)
                    for sf in SF_CLASSES + ["_total"]}
                for c, L in sizes.items()}

    if os.path.exists(bed_path):
        with open(bed_path) as fh:
            for ln in fh:
                f = ln.split("\t")
                if len(f) < 4: continue
                c = f[0]
                if c not in bp_class: continue
                try:
                    s, e = int(f[1]), int(f[2])
                except ValueError:
                    continue
                sf = classify_superfamily(f[3])
                if sf is None: continue
                L  = len(bp_class[c]["_total"])
                e  = min(e, L); s = max(0, s)
                if s >= e: continue
                bp_class[c]["_total"][s:e] = 1.0
                bp_class[c][sf][s:e]       = 1.0

    rows = []
    for c, L in sizes.items():
        n = L // win
        if n == 0: continue
        tot_arr  = bp_class[c]["_total"][:n*win].reshape(n, win).mean(axis=1)
        sf_arrs  = {sf: bp_class[c][sf][:n*win].reshape(n, win).mean(axis=1)
                    for sf in SF_CLASSES}
        for i in range(n):
            te_tot = float(tot_arr[i])
            row = {"chrom": c, "win_idx": i, "te_total_frac": te_tot}
            for sf in SF_CLASSES:
                row[f"te_{sf}"] = float(sf_arrs[sf][i])
            # Dominant: max-coverage SF class (must beat MIN_TE_FRAC)
            best_sf   = max(SF_CLASSES, key=lambda s: sf_arrs[s][i])
            best_frac = float(sf_arrs[best_sf][i])
            row["dominant_sf"] = best_sf if best_frac >= MIN_TE_FRAC else "TE-poor"
            rows.append(row)

    return pd.DataFrame(rows)

def mwu_p(x, y):
    if len(x) < 3 or len(y) < 3: return np.nan
    try: _, p = mannwhitneyu(x, y, alternative="two-sided"); return p
    except: return np.nan

def cliffs_d(x, y, nmax=5000):
    if len(x) > nmax: x = np.random.choice(x, nmax, replace=False)
    if len(y) > nmax: y = np.random.choice(y, nmax, replace=False)
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x)==0 or len(y)==0: return np.nan
    gt = np.sum(x[:,None] > y[None,:]); lt = np.sum(x[:,None] < y[None,:])
    return float((gt - lt) / (len(x)*len(y)))

def make_violin_ax(ax, df, col, present_sf):
    """Draw violin plot on ax; return (parts, colors)."""
    data   = [df.loc[df["dominant_sf"]==sf, col].dropna().values for sf in present_sf]
    colors = [SF_COLORS[sf] for sf in present_sf]
    valid  = [(d, c) for d, c in zip(data, colors) if len(d) >= 5]
    if not valid: return
    vdata, vcolors = zip(*valid)
    valid_sf = [sf for sf, d in zip(present_sf, data) if len(d) >= 5]
    parts = ax.violinplot(vdata, positions=np.arange(len(valid_sf)),
                          showmedians=True, showextrema=False, widths=0.75)
    for pc, color in zip(parts["bodies"], vcolors):
        pc.set_facecolor(color); pc.set_alpha(0.80); pc.set_edgecolor("black"); pc.set_lw(0.7)
    parts["cmedians"].set_color("black"); parts["cmedians"].set_lw(1.8)
    ax.set_xticks(np.arange(len(valid_sf)))
    n_per = [(df["dominant_sf"]==sf).sum() for sf in valid_sf]
    ax.set_xticklabels([f"{sf}\n(n={n})" for sf, n in zip(valid_sf, n_per)],
                       fontsize=8, rotation=40, ha="right")
    if col == "mean_sig": ax.set_yscale("log")
    if col == "ac_lag":   ax.set_ylim(140, 240)
    ax.grid(axis="y", alpha=0.3, lw=0.6)

# ============================================================= PART C / D ===
def run_superfamily_panel(panel, part_letter):
    title_label = "In vivo" if part_letter == "C" else "In vitro (NucSeq)"
    print(f"\n{'='*60}")
    print(f"PART {part_letter}: {title_label} — chromatin by TE superfamily ({WIN//1000} kb)")
    print("="*60)

    species_results = {}
    summary_rows    = []

    for entry in panel:
        lab = entry["label"]; sp = entry["sp"]; wig = entry["wig"]
        if not os.path.exists(wig):
            print(f"[SKIP] {lab}: {wig}"); continue
        print(f"\n[{lab}]")

        sz   = chrom_sizes(sp)
        arrs = load_nuclpred_wig(wig, sz)

        # TE superfamily window table
        sf_df = te_superfamily_windows(te_bed(sp), sz, WIN)
        print("  Dominant SF counts:")
        vc = sf_df["dominant_sf"].value_counts()
        for sf in SF_ORDER:
            if sf in vc.index: print(f"    {sf:12s}: {vc[sf]:6d}")

        # Chromatin metrics per window
        rows_c = []
        for c, occ in arrs.items():
            n = len(occ) // WIN
            if n == 0: continue
            means, cvs, aclags, acvals, _ = window_metrics(occ, WIN)
            rows_c.append(pd.DataFrame({
                "chrom": c, "win_idx": np.arange(n),
                "mean_sig": means, "cv": cvs, "ac_val": acvals, "ac_lag": aclags,
            }))
        if not rows_c: continue
        chrom_df = pd.concat(rows_c, ignore_index=True)

        df = chrom_df.merge(sf_df[["chrom","win_idx","dominant_sf","te_total_frac"]],
                            on=["chrom","win_idx"], how="inner")
        df = df[df["mean_sig"] > 1e-3].copy()
        species_results[lab] = (df, entry)

        # ---- per-species figure ------------------------------------------
        present_sf = [s for s in SF_ORDER if (df["dominant_sf"]==s).sum() >= MIN_SF_WINDOWS]
        fig, axes = plt.subplots(1, len(FEATS), figsize=(4.5*len(FEATS), 5.0))
        for ax, (col, ttl) in zip(axes, FEATS):
            make_violin_ax(ax, df, col, present_sf)
            ax.set_title(ttl, fontsize=11)
        atac_note = " [ATAC: high=accessible]" if entry.get("atac") else ""
        fig.suptitle(f"{lab} ({entry['assay']}{atac_note}) — {WIN//1000} kb windows by dominant TE superfamily",
                     fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        out_pref = f"part{part_letter}_sf_{lab}"
        fig.savefig(os.path.join(OUT, f"{out_pref}.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(OUT, f"{out_pref}.png"), dpi=180, bbox_inches="tight")
        plt.close(fig)

        # ---- stats vs TE-poor baseline -----------------------------------
        for col, _ in FEATS:
            poor = df.loc[df["dominant_sf"]=="TE-poor", col].dropna().values
            for sf in SF_CLASSES:
                vals = df.loc[df["dominant_sf"]==sf, col].dropna().values
                if len(vals) < MIN_SF_WINDOWS: continue
                summary_rows.append({
                    "species": lab, "assay": entry["assay"],
                    "part": part_letter, "metric": col, "superfamily": sf,
                    "n_windows": len(vals),
                    "median": float(np.median(vals)),
                    "mean":   float(np.mean(vals)),
                    "d_vs_tepoor": cliffs_d(vals, poor),
                    "p_vs_tepoor": mwu_p(vals, poor),
                })

    # ---- grid figure (metrics × species) --------------------------------
    labs = list(species_results.keys())
    if labs:
        fig, axes = plt.subplots(len(FEATS), len(labs),
                                 figsize=(3.8*len(labs), 3.8*len(FEATS)), squeeze=False)
        for j, lab in enumerate(labs):
            df, entry = species_results[lab]
            present_sf = [s for s in SF_ORDER if (df["dominant_sf"]==s).sum() >= MIN_SF_WINDOWS]
            for i, (col, ttl) in enumerate(FEATS):
                ax = axes[i][j]
                make_violin_ax(ax, df, col, present_sf)
                if i == 0:
                    atac_note = "\n[ATAC↑=open]" if entry.get("atac") else ""
                    ax.set_title(f"{lab}\n{entry['assay']}{atac_note}", fontsize=9, fontweight="bold")
                if j == 0: ax.set_ylabel(ttl, fontsize=9)
        fig.suptitle(f"{title_label} chromatin by dominant TE superfamily ({WIN//1000} kb windows)",
                     fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(os.path.join(OUT, f"part{part_letter}_sf_grid.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(OUT, f"part{part_letter}_sf_grid.png"), dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"\n[part {part_letter}] grid saved  ({len(labs)} species)")

    return species_results, summary_rows


# ========================================= COMPARISON in vivo vs in vitro ===
def make_comparison_figure(invivo_res, invitro_res):
    """
    Two panels per metric (NRL + CV):
      left  = in vivo  (colored violins)
      right = in vitro (same layout, hatched)
    One row per metric.
    """
    COMPARE_FEATS = [
        ("ac_lag", "NRL (bp)"),
        ("cv",     "CV (rugosity)"),
    ]
    CAP_PER_SP = 3000  # max windows per species when pooling

    print("\n[partCD] building comparison figure …")

    fig, axes = plt.subplots(len(COMPARE_FEATS), 2,
                             figsize=(15, 5.5*len(COMPARE_FEATS)), squeeze=False)

    for i, (col, ttl) in enumerate(COMPARE_FEATS):
        for j, (panel_label, results) in enumerate([("In vivo", invivo_res),
                                                     ("In vitro (NucSeq)", invitro_res)]):
            ax = axes[i][j]
            pooled = {sf: [] for sf in SF_ORDER}
            for lab, (df, entry) in results.items():
                for sf in SF_ORDER:
                    vals = df.loc[df["dominant_sf"]==sf, col].dropna().values
                    if len(vals) > 0:
                        pooled[sf].extend(
                            vals[np.random.choice(len(vals),
                                 min(CAP_PER_SP, len(vals)), replace=False)].tolist())

            present = [sf for sf in SF_ORDER if len(pooled[sf]) >= MIN_SF_WINDOWS]
            data    = [pooled[sf] for sf in present]
            colors  = [SF_COLORS[sf] for sf in present]

            parts = ax.violinplot(data, positions=np.arange(len(present)),
                                  showmedians=True, showextrema=False, widths=0.80)
            for pc, color in zip(parts["bodies"], colors):
                pc.set_facecolor(color); pc.set_alpha(0.80)
                pc.set_edgecolor("black"); pc.set_lw(0.8)
                if j == 1: pc.set_hatch("//"); pc.set_alpha(0.60)
            parts["cmedians"].set_color("black"); parts["cmedians"].set_lw(2)

            ax.set_xticks(np.arange(len(present)))
            n_per = [len(pooled[sf]) for sf in present]
            ax.set_xticklabels([f"{sf}\n(n={n})" for sf, n in zip(present, n_per)],
                               fontsize=9, rotation=40, ha="right")
            if col == "ac_lag": ax.set_ylim(140, 240)
            ax.set_title(panel_label, fontsize=12, fontweight="bold")
            ax.set_ylabel(ttl, fontsize=11)
            ax.grid(axis="y", alpha=0.3, lw=0.6)
            # Grand median reference line
            all_vals = np.concatenate([d for d in data if d])
            if len(all_vals):
                ax.axhline(np.median(all_vals), color="#555555", ls="--", lw=0.9, alpha=0.5)

    fig.suptitle(f"Chromatin by TE superfamily — in vivo vs in vitro pooled\n"
                 f"(pooled across species, ≤{CAP_PER_SP} windows/species, {WIN//1000} kb)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(OUT, "partCD_comparison.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "partCD_comparison.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  → partCD_comparison saved")


def make_mereb1_comparison(invivo_res, invitro_res):
    """
    BA_Mereb1 special: ATAC (in vivo) vs NucSeq (in vitro) — same genome.
    Side-by-side for all 4 metrics.
    """
    atac_key = "BA_Mereb1_atac" if "BA_Mereb1_atac" in invivo_res else None
    nucs_key = "BA_Mereb1"      if "BA_Mereb1"      in invitro_res else None
    if atac_key is None or nucs_key is None:
        print("[partCD] Mereb1 comparison skipped — one dataset missing"); return

    df_atac, e_atac = invivo_res[atac_key]
    df_nucs, e_nucs = invitro_res[nucs_key]

    fig, axes = plt.subplots(len(FEATS), 2,
                             figsize=(11, 4.0*len(FEATS)), squeeze=False)

    for i, (col, ttl) in enumerate(FEATS):
        for j, (df, entry, subtitle) in enumerate([
                (df_atac, e_atac, "BA_Mereb1 ATAC (in vivo)"),
                (df_nucs, e_nucs, "BA_Mereb1 NucSeq (in vitro)"),
        ]):
            ax = axes[i][j]
            present_sf = [s for s in SF_ORDER if (df["dominant_sf"]==s).sum() >= MIN_SF_WINDOWS]
            make_violin_ax(ax, df, col, present_sf)
            if i == 0: ax.set_title(subtitle, fontsize=10, fontweight="bold")
            if j == 0: ax.set_ylabel(ttl, fontsize=9)

    fig.suptitle("BA_Mereb1: in vivo ATAC vs in vitro NucSeq\n(same genome — tests DNA-sequence encoding)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(OUT, "partCD_mereb1_invivo_vs_invitro.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "partCD_mereb1_invivo_vs_invitro.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  → partCD_mereb1 comparison saved")


# ======================================================================= main
def main():
    np.random.seed(42)

    invivo_res,  rows_C = run_superfamily_panel(PANEL_C, "C")
    invitro_res, rows_D = run_superfamily_panel(PANEL_D, "D")

    # Summary TSV
    all_rows = rows_C + rows_D
    if all_rows:
        summ = pd.DataFrame(all_rows)
        out_tsv = os.path.join(OUT, "partCD_summary.tsv")
        summ.to_csv(out_tsv, sep="\t", index=False)
        print(f"\n[partCD] summary TSV: {out_tsv}")
        # Quick print
        print("\n=== PART C: NRL by superfamily (median, Cliff's d vs TE-poor) ===")
        for row in [r for r in rows_C if r["metric"]=="ac_lag"]:
            print(f"  {row['species']:20s}  {row['superfamily']:10s}  "
                  f"NRL={row['median']:.1f} bp  d={row['d_vs_tepoor']:+.3f}  "
                  f"p={row['p_vs_tepoor']:.2e}")
        print("\n=== PART D: NRL by superfamily (median, Cliff's d vs TE-poor) ===")
        for row in [r for r in rows_D if r["metric"]=="ac_lag"]:
            print(f"  {row['species']:20s}  {row['superfamily']:10s}  "
                  f"NRL={row['median']:.1f} bp  d={row['d_vs_tepoor']:+.3f}  "
                  f"p={row['p_vs_tepoor']:.2e}")

    make_comparison_figure(invivo_res, invitro_res)
    make_mereb1_comparison(invivo_res, invitro_res)

    print(f"\n[DONE] all outputs in {OUT}")

if __name__ == "__main__":
    main()
