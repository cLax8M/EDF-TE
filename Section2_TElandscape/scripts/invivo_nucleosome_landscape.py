#!/usr/bin/env python3
"""
invivo_nucleosome_landscape.py
==============================
Integrated analysis of native chromatin structure in the 4 species with
real in vivo chromatin data:

  MU_Rhimi59      MNase-seq  (legacy WIG only)
  AS_Ncrass2      MNase-seq  (DANPOS positions.xls + WIG)
  AS_Bcine1       MAINE-seq  (DANPOS positions.xls + WIG)
  BA_Mereb1_atac  ATAC-seq   (WIG only; signal = accessibility, inverted)

=== Part A: Per-nucleosome metrics (Ncrass2 + Bcine1 only) ===
Source: DANPOS positions.xls — NOT affected by multi-mapping coverage bias
Groups: TE  /  Gene  /  Intergenic
Metrics:
  fuzziness_score  — SD of dyad position (lower = better positioned)
  smt_value        — peak height / sharpness (higher = more nucleosome-like)
  nrl_bp           — distance to next nucleosome (100–400 bp filter)

=== Part B: Regional analysis (all 4 in vivo) ===
Source: WIG occupancy (or accessibility for ATAC)
Windows: 5 kb AND 10 kb
TE density quartiles Q1 (TE-poor) → Q4 (TE-rich)
Metrics:
  mean_sig  — mean WIG signal per window
  cv        — coefficient of variation (rugosity)
  ac_val    — autocorrelation peak 140-240 bp (phasing strength)
  ac_lag    — position of that peak (NRL)

Outputs:
  .../invivo_landscape/
    partA_nucleosome_groups_grid.{pdf,png}    (fuzz+sharp+NRL, 2 species × 3 metrics)
    partA_<LABEL>_boxplots.{pdf,png}
    partA_summary.tsv
    partB_5kb_grid.{pdf,png}
    partB_10kb_grid.{pdf,png}
    partB_<WIN>kb_<LABEL>.{pdf,png}
    partB_summary_5kb.tsv
    partB_summary_10kb.tsv
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics

# ------------------------------------------------------------------ paths ---
MV_NUC = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes/nuc-seq"
MV_ATC = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes/atac-seq"
GFF_DIR  = "/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff"
EXTRA    = "/legserv/Temp/Carlos/EDF-TE/data/Extra"
FAI_FMT = "/legserv/Temp/Carlos/EDF-TE/data/genomes/{sp}.fasta.fai"
TE_DATA = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
TE_MV   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"

OUT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/invivo_landscape"
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ panel ---
# Part A species (have DANPOS positions.xls)
PANEL_A = [
    dict(label="MU_Rhimi59",  sp="MU_Rhimi59", assay="MNase",
         xls=f"{EXTRA}/Users_carloslax_Desktop_Ficheros_mta1_RmNewAssembly.Fnor.smooth.positions.ref_adjust.xls",
         wig="/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig",
         gff=f"{GFF_DIR}/MU_Rhimi59.gff"),
    dict(label="AS_Ncrass2", sp="AS_Ncrass2", assay="MNase",
         xls=f"{MV_NUC}/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.positions.xls",
         wig=f"{MV_NUC}/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.wig",
         gff=f"{GFF_DIR}/AS_Ncrass2.gff"),
    dict(label="AS_Bcine1",  sp="AS_Bcine1",  assay="MAINE",
         xls=f"{MV_NUC}/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.positions.xls",
         wig=f"{MV_NUC}/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.wig",
         gff=f"{GFF_DIR}/AS_Bcine1.gff3"),
]

# Part B species (WIG only — includes Rhimi59 legacy + ATAC Mereb1)
PANEL_B = [
    dict(label="MU_Rhimi59",     sp="MU_Rhimi59",  assay="MNase", atac=False,
         wig="/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"),
    dict(label="AS_Ncrass2",     sp="AS_Ncrass2",  assay="MNase", atac=False,
         wig=f"{MV_NUC}/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.wig"),
    dict(label="AS_Bcine1",      sp="AS_Bcine1",   assay="MAINE", atac=False,
         wig=f"{MV_NUC}/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.wig"),
    dict(label="BA_Mereb1_atac", sp="BA_Mereb1",   assay="ATAC",  atac=True,
         wig=f"{MV_ATC}/BA_Mereb1/BA_Mereb1_atac.wig"),
]

# ---------------------------------------------------------------- helpers ---
def te_bed(sp):
    a = TE_DATA.format(sp=sp); b = TE_MV.format(sp=sp)
    return a if os.path.exists(a) else b

def chrom_sizes(sp):
    fai = FAI_FMT.format(sp=sp)
    df = pd.read_csv(fai, sep="\t", header=None, names=["c","L","o","lb","lby"])
    return dict(zip(df["c"], df["L"].astype(int)))

def load_intervals(path, is_gff=False):
    """Return dict chrom -> sorted Nx2 int32 array."""
    d = {}
    if not os.path.exists(path):
        return d
    use = None  # for gff: detect feature type
    if is_gff:
        want = set()
        with open(path) as fh:
            for ln in fh:
                if ln.startswith("#"): continue
                f = ln.split("\t")
                if len(f) >= 3: want.add(f[2])
        use = "gene" if "gene" in want else "exon"
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.split("\t")
            if len(f) < 5: continue
            if is_gff:
                if f[2] != use: continue
                c = f[0]
                try: s, e = int(f[3])-1, int(f[4])
                except ValueError: continue
            else:
                c = f[0]
                try: s, e = int(f[1]), int(f[2])
                except ValueError: continue
            d.setdefault(c, []).append((s, e))
    return {c: np.array(sorted(v), dtype=np.int32) for c, v in d.items()}

def overlap_vec(df, ivs_dict):
    """Bool array: True where smt_pos falls inside any interval."""
    result = np.zeros(len(df), dtype=bool)
    for c, grp in df.groupby("chr"):
        ivs = ivs_dict.get(c)
        if ivs is None or len(ivs) == 0: continue
        pos = grp["smt_pos"].values
        idx = np.searchsorted(ivs[:, 0], pos, side="right") - 1
        in_t = np.where(idx >= 0, pos < ivs[idx.clip(min=0), 1], False)
        result[grp.index.values] = in_t
    return result

def compute_nrl(df):
    df = df.sort_values(["chr", "smt_pos"]).copy()
    df["nrl_bp"] = np.nan
    for c, grp in df.groupby("chr"):
        idx = grp.index; pos = grp["smt_pos"].values
        if len(pos) < 2: continue
        d = np.diff(pos).astype(float)
        d[(d < 100) | (d > 400)] = np.nan
        df.loc[idx[:-1], "nrl_bp"] = d
    return df

def cliffs_d(x, y, nmax=5000):
    x = np.asarray(x); y = np.asarray(y)
    if len(x) == 0 or len(y) == 0: return np.nan
    rng = np.random.RandomState(42)
    if len(x) > nmax: x = rng.choice(x, nmax, replace=False)
    if len(y) > nmax: y = rng.choice(y, nmax, replace=False)
    return float((x[:, None] > y[None, :]).sum() - (x[:, None] < y[None, :]).sum()) / (len(x)*len(y))

def mwu_p(x, y):
    if len(x) < 3 or len(y) < 3: return np.nan
    try: _, p = mannwhitneyu(x, y, alternative="two-sided"); return float(p)
    except: return np.nan

# ============================================================= PART A ======
def run_partA():
    print("\n" + "="*60)
    print("PART A: Per-nucleosome metrics (DANPOS positions.xls)")
    print("="*60)

    ORDER  = ["Gene", "Intergenic", "TE"]
    COLORS = {"Gene": "#2CA02C", "Intergenic": "#7F7F7F", "TE": "#D62728"}
    COLS   = [
        ("fuzziness_score", "Fuzziness\n(lower = better positioned)"),
        ("smt_value",       "Sharpness\n(log scale)"),
        ("nrl_bp",          "NRL\n(bp, 100–400 filter)"),
    ]
    summary_rows = []
    species_dfs = {}

    for entry in PANEL_A:
        lab = entry["label"]; sp = entry["sp"]
        print(f"\n[{lab}]")
        df = pd.read_csv(entry["xls"], sep="\t")
        print(f"  {len(df):,} nucleosomes")

        te_ivs   = load_intervals(te_bed(sp))
        gene_ivs = load_intervals(entry["gff"], is_gff=True)
        df["in_te"]   = overlap_vec(df, te_ivs)
        df["in_gene"] = overlap_vec(df, gene_ivs)
        df["group"] = "Intergenic"
        df.loc[df["in_gene"] & ~df["in_te"], "group"] = "Gene"
        df.loc[df["in_te"],                  "group"] = "TE"
        df = compute_nrl(df)

        counts = df["group"].value_counts().to_dict()
        print(f"  {counts}")

        # per-species figure
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        for ax, (col, ttl) in zip(axes, COLS):
            data = [df.loc[df["group"]==g, col].dropna().values for g in ORDER]
            bp = ax.boxplot(data, positions=np.arange(3), widths=0.6,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.5))
            for patch, g in zip(bp["boxes"], ORDER):
                patch.set_facecolor(COLORS[g]); patch.set_edgecolor("black")
            ax.set_xticks(np.arange(3))
            ax.set_xticklabels([f"{g}\n(n={counts.get(g,0):,})" for g in ORDER], fontsize=9)
            ax.set_title(ttl, fontsize=10)
            if col == "smt_value": ax.set_yscale("log")
            if col == "nrl_bp":    ax.set_ylim(100, 400)
        fig.suptitle(f"{lab} ({entry['assay']}, in vivo) — per-nucleosome metrics", fontsize=12)
        fig.tight_layout(rect=[0,0,1,0.94])
        fig.savefig(os.path.join(OUT, f"partA_{lab}.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(OUT, f"partA_{lab}.png"), dpi=180, bbox_inches="tight")
        plt.close(fig)

        # stats
        row = {"species": lab, "assay": entry["assay"]}
        for col, _ in COLS:
            for g1, g2 in [("TE","Gene"), ("TE","Intergenic"), ("Gene","Intergenic")]:
                x = df.loc[df["group"]==g1, col].dropna().values
                y = df.loc[df["group"]==g2, col].dropna().values
                k = f"{col}__{g1}_vs_{g2}"
                row[f"{k}_delta"] = float(np.mean(x)-np.mean(y)) if len(x) and len(y) else np.nan
                row[f"{k}_d"]     = cliffs_d(x, y)
                row[f"{k}_p"]     = mwu_p(x, y)
            for g in ORDER:
                v = df.loc[df["group"]==g, col].dropna().values
                row[f"{col}_mean_{g}"] = float(np.mean(v)) if len(v) else np.nan
        summary_rows.append(row)
        species_dfs[lab] = df

    # combined grid
    fig, axes = plt.subplots(len(COLS), len(PANEL_A),
                             figsize=(5.5*len(PANEL_A), 4*len(COLS)), squeeze=False)
    for j, entry in enumerate(PANEL_A):
        lab = entry["label"]; df = species_dfs[lab]
        counts = df["group"].value_counts().to_dict()
        for i, (col, ttl) in enumerate(COLS):
            ax = axes[i][j]
            data = [df.loc[df["group"]==g, col].dropna().values for g in ORDER]
            bp = ax.boxplot(data, positions=np.arange(3), widths=0.6,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.3))
            for patch, g in zip(bp["boxes"], ORDER):
                patch.set_facecolor(COLORS[g]); patch.set_edgecolor("black")
            ax.set_xticks(np.arange(3))
            ax.set_xticklabels([f"{g}\n({counts.get(g,0):,})" for g in ORDER], fontsize=9)
            if col == "smt_value": ax.set_yscale("log")
            if col == "nrl_bp":    ax.set_ylim(100, 400)
            if i == 0: ax.set_title(f"{lab}\n{entry['assay']} in vivo", fontsize=11, fontweight="bold")
            if j == 0: ax.set_ylabel(ttl, fontsize=10)
    fig.suptitle("Per-nucleosome DANPOS metrics — in vivo MNase/MAINE-seq\nTE vs Protein-coding genes vs Intergenic", fontsize=13)
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(os.path.join(OUT, "partA_nucleosome_groups_grid.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "partA_nucleosome_groups_grid.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)

    # summary TSV
    summ = pd.DataFrame(summary_rows)
    summ.to_csv(os.path.join(OUT, "partA_summary.tsv"), sep="\t", index=False)

    print("\n=== PART A: TE vs Gene (Cliff's d) ===")
    show = [["species","assay"]]
    for col, _ in COLS:
        print(f"\n  {col}:")
        for row in summary_rows:
            print(f"    {row['species']}: d(TE-Gene)={row[f'{col}__TE_vs_Gene_d']:.3f}  "
                  f"Δ={row[f'{col}__TE_vs_Gene_delta']:.2f}  "
                  f"p={row[f'{col}__TE_vs_Gene_p']:.2e}")
    return species_dfs

# ============================================================= PART B ======
def te_density_per_window(bed_path, sizes, win):
    bp = {c: np.zeros(L, dtype=np.uint8) for c, L in sizes.items()}
    if os.path.exists(bed_path):
        with open(bed_path) as fh:
            for ln in fh:
                f = ln.split("\t")
                if len(f) < 3: continue
                c = f[0]
                if c not in bp: continue
                try: s, e = int(f[1]), int(f[2])
                except ValueError: continue
                e = min(e, len(bp[c])); s = max(0, s)
                if s < e: bp[c][s:e] = 1
    out = {}
    for c, L in sizes.items():
        n = L // win
        if n == 0: out[c] = np.zeros(0, dtype=np.float32); continue
        a = bp[c][:n*win].reshape(n, win)
        out[c] = a.sum(axis=1).astype(np.float32) / win
    return out

def run_partB_window(win, label_suffix=""):
    print(f"\n{'='*60}")
    print(f"PART B: Regional analysis — {win//1000} kb windows {label_suffix}")
    print("="*60)

    ORDER = ["Q1","Q2","Q3","Q4"]
    cmap = plt.cm.viridis
    COLORS = {q: cmap(i/3) for i, q in enumerate(ORDER)}
    FEATS = [
        ("mean_sig", "Mean signal"),
        ("cv",       "CV (rugosity)"),
        ("ac_val",   "AC peak (phasing)"),
        ("ac_lag",   "NRL (bp)"),
    ]
    summary_rows = []; species_dfs = {}

    for entry in PANEL_B:
        lab = entry["label"]; sp = entry["sp"]
        wig = entry["wig"]
        if not os.path.exists(wig):
            print(f"[SKIP] {lab}: {wig}"); continue
        print(f"\n[{lab}]")
        sz   = chrom_sizes(sp)
        arrs = load_nuclpred_wig(wig, sz)
        te_d = te_density_per_window(te_bed(sp), sz, win)

        rows = []
        for c, occ in arrs.items():
            n = len(occ) // win
            if n == 0: continue
            means, cvs, aclags, acvals, _ = window_metrics(occ, win)
            tp = te_d.get(c, np.zeros(n, dtype=np.float32))[:n]
            rows.append(pd.DataFrame({
                "chrom": c, "start": np.arange(n)*win, "end": np.arange(n)*win+win,
                "te_pct": tp,
                "mean_sig": means, "cv": cvs, "ac_val": acvals, "ac_lag": aclags,
            }))
        if not rows: print(f"  no windows"); continue
        df = pd.concat(rows, ignore_index=True)
        df = df[df["mean_sig"] > 1e-3].copy()
        df["quartile"] = pd.qcut(df["te_pct"].rank(method="first"), 4,
                                 labels=ORDER).astype(str)
        print(f"  {len(df):,} windows  TE% Q1={df[df['quartile']=='Q1']['te_pct'].mean()*100:.1f}%  "
              f"Q4={df[df['quartile']=='Q4']['te_pct'].mean()*100:.1f}%")
        species_dfs[lab] = df

        # per-species figure — 4 quartile boxplots
        fig, axes = plt.subplots(1, 4, figsize=(15, 3.8))
        for ax, (col, ttl) in zip(axes, FEATS):
            data = [df.loc[df["quartile"]==q, col].values for q in ORDER]
            bp = ax.boxplot(data, positions=np.arange(4), widths=0.65,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.5))
            for patch, q in zip(bp["boxes"], ORDER):
                patch.set_facecolor(COLORS[q]); patch.set_edgecolor("black")
            ax.set_xticks(np.arange(4))
            ax.set_xticklabels([f"{q}\n({(df['quartile']==q).sum():,})" for q in ORDER], fontsize=9)
            ax.set_title(ttl, fontsize=10)
            if col == "mean_sig": ax.set_yscale("log")
            if col == "ac_lag":   ax.set_ylim(140, 240)
        atac_note = " [ATAC: high=accessible]" if entry["atac"] else ""
        fig.suptitle(f"{lab} ({entry['assay']}{atac_note}) — {win//1000} kb windows by TE-density quartile",
                     fontsize=11)
        fig.tight_layout(rect=[0,0,1,0.93])
        fig.savefig(os.path.join(OUT, f"partB_{win//1000}kb_{lab}.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(OUT, f"partB_{win//1000}kb_{lab}.png"), dpi=180, bbox_inches="tight")
        plt.close(fig)

        # stats Q4 vs Q1
        row = {"species": lab, "assay": entry["assay"], "atac": entry["atac"],
               "win_kb": win//1000, "n_windows": len(df),
               "te_pct_Q1": float(df[df["quartile"]=="Q1"]["te_pct"].mean()),
               "te_pct_Q4": float(df[df["quartile"]=="Q4"]["te_pct"].mean())}
        for col, _ in FEATS:
            x = df.loc[df["quartile"]=="Q4", col].dropna().values
            y = df.loc[df["quartile"]=="Q1", col].dropna().values
            row[f"{col}_Q4_mean"]  = float(np.mean(x)) if len(x) else np.nan
            row[f"{col}_Q1_mean"]  = float(np.mean(y)) if len(y) else np.nan
            row[f"{col}_delta"]    = row[f"{col}_Q4_mean"] - row[f"{col}_Q1_mean"]
            row[f"{col}_d"]        = cliffs_d(x, y)
            row[f"{col}_p"]        = mwu_p(x, y)
            # also Kruskal across all 4 quartiles
            try:
                _, row[f"{col}_kruskal_p"] = kruskal(
                    *[df.loc[df["quartile"]==q, col].dropna().values for q in ORDER])
            except: row[f"{col}_kruskal_p"] = np.nan
        summary_rows.append(row)

    # grid figure (4 features × N species)
    labs = list(species_dfs.keys())
    if labs:
        fig, axes = plt.subplots(len(FEATS), len(labs),
                                 figsize=(3.0*len(labs), 3.2*len(FEATS)), squeeze=False)
        for j, lab in enumerate(labs):
            df = species_dfs[lab]
            meta = next(e for e in PANEL_B if e["label"]==lab)
            for i, (col, ttl) in enumerate(FEATS):
                ax = axes[i][j]
                data = [df.loc[df["quartile"]==q, col].values for q in ORDER]
                bp = ax.boxplot(data, positions=np.arange(4), widths=0.65,
                                showfliers=False, patch_artist=True,
                                medianprops=dict(color="black", lw=1.2))
                for patch, q in zip(bp["boxes"], ORDER):
                    patch.set_facecolor(COLORS[q]); patch.set_edgecolor("black")
                ax.set_xticks(np.arange(4)); ax.set_xticklabels(ORDER, fontsize=8)
                if col == "mean_sig": ax.set_yscale("log")
                if col == "ac_lag":   ax.set_ylim(140, 240)
                if i == 0:
                    atac_note = "\n[ATAC↑=open]" if meta["atac"] else ""
                    ax.set_title(f"{lab}\n{meta['assay']} in vivo{atac_note}",
                                 fontsize=9, fontweight="bold")
                if j == 0: ax.set_ylabel(ttl, fontsize=9)
        fig.suptitle(f"In vivo chromatin — {win//1000} kb windows by TE-density quartile\n"
                     "Q1=TE-poor → Q4=TE-rich", fontsize=12)
        fig.tight_layout(rect=[0,0,1,0.95])
        fig.savefig(os.path.join(OUT, f"partB_{win//1000}kb_grid.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(OUT, f"partB_{win//1000}kb_grid.png"), dpi=180, bbox_inches="tight")
        plt.close(fig)

    # summary TSV
    summ = pd.DataFrame(summary_rows)
    summ.to_csv(os.path.join(OUT, f"partB_summary_{win//1000}kb.tsv"), sep="\t", index=False)

    print(f"\n=== PART B {win//1000} kb: Q4 vs Q1 (Cliff's d) ===")
    for col, ttl in FEATS:
        print(f"\n  {col}:")
        for row in summary_rows:
            print(f"    {row['species']:20s}: d={row[f'{col}_d']:.3f}  "
                  f"Δ(Q4-Q1)={row[f'{col}_delta']:+.3g}  p={row[f'{col}_p']:.2e}")
    return summ


# =================================================================== main ---
def main():
    run_partA()
    run_partB_window(5_000,  "(5 kb)")
    run_partB_window(10_000, "(10 kb)")
    print(f"\n[DONE] all outputs in {OUT}")

if __name__ == "__main__":
    main()
