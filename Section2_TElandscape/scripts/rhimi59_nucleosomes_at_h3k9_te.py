#!/usr/bin/env python3
"""
rhimi59_nucleosomes_at_h3k9_te.py
=================================
Descriptive analysis: in Rhizopus microsporus (MU_Rhimi59), look at the
nucleosome features of 500 bp windows split into 4 categories:

  TE+ H3K9+   : window overlaps a TE AND an H3K9me3 broad peak
  TE+ H3K9-   : window overlaps a TE but NOT a H3K9me3 peak
  TE- H3K9+   : H3K9me3 peak without annotated TE
  TE- H3K9-   : neither (background)

Features per window (real DANPOS MNase-seq, no model, no rank-norm):
  mean_occ      mean nucleosome occupancy
  cv            sd / mean (rugosity; high = NDR-rich, low = smooth)
  ac_val        autocorrelation peak in 140-240 bp (phasing strength)
  ac_lag        position of that peak (nucleosome repeat length, NRL)

Outputs:
  rhimi59_h3k9_te_windows.tsv.gz
  rhimi59_h3k9_te_features.pdf   (4 boxplots + occupancy histogram)
  rhimi59_h3k9_te_summary.tsv    (mean / median per category)
"""
import os, gzip, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics

# ----- inputs -----
WIG  = "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"
FAI  = "/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Rhimi59.fasta.fai"
H3K9 = "/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/Epigenetics/RmH3K9me3/H3K9_WT_MACS2_broad_peaks.bed"
TE   = "/data/carlos/EDF-TE/TE_analysis/annotated/MU_Rhimi59_EarlGrey/MU_Rhimi59_summaryFiles/MU_Rhimi59.filteredRepeats.bed"
OUT  = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/rhimi59_h3k9_te"
os.makedirs(OUT, exist_ok=True)
WIN = 500
MIN_OVERLAP = 0.5     # >=50% of the window covered by feature -> "+"

# ----- helpers -----
def chrom_sizes_from_fai(p):
    fai = pd.read_csv(p, sep="\t", header=None,
                      names=["chrom","L","off","lb","lby"])
    return dict(zip(fai["chrom"], fai["L"].astype(int)))

def overlap_array(bed_path, sizes, win):
    """Return dict chrom -> 1D bool array of length n_windows where True
    means window has >=MIN_OVERLAP fraction overlapping any feature."""
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
            e = min(e, len(bp[c]))
            if s < 0: s = 0
            if s < e:
                bp[c][s:e] = 1
    out = {}
    for c, L in sizes.items():
        n = L // win
        if n == 0:
            out[c] = np.zeros(0, dtype=bool); continue
        a = bp[c][:n*win].reshape(n, win)
        out[c] = (a.mean(axis=1) >= MIN_OVERLAP)
    return out

# ----- compute -----
print("[1] sizes")
sizes = chrom_sizes_from_fai(FAI)
print("[2] WIG")
arrs = load_nuclpred_wig(WIG, sizes)
print("[3] window metrics")
rows = []
for c, occ in arrs.items():
    n = len(occ) // WIN
    if n == 0: continue
    means, cvs, aclags, acvals, _ = window_metrics(occ, WIN)
    rows.append(pd.DataFrame({
        "chrom": c,
        "start": np.arange(n) * WIN,
        "end":   np.arange(n) * WIN + WIN,
        "mean_occ": means,
        "cv": cvs,
        "ac_val": acvals,
        "ac_lag": aclags,
    }))
df = pd.concat(rows, ignore_index=True)
print(f"  windows: {len(df)}")

print("[4] H3K9 + TE overlap")
te_ov   = overlap_array(TE,   sizes, WIN)
h3k9_ov = overlap_array(H3K9, sizes, WIN)
te_flag = np.zeros(len(df), dtype=bool)
h3_flag = np.zeros(len(df), dtype=bool)
i0 = 0
for c, sub in df.groupby("chrom", sort=False):
    n = len(sub)
    te_flag[i0:i0+n] = te_ov[c][:n]
    h3_flag[i0:i0+n] = h3k9_ov[c][:n]
    i0 += n
df["te"]   = te_flag
df["h3k9"] = h3_flag
df["category"] = np.where(df["te"] & df["h3k9"],   "TE+ H3K9+",
                  np.where(df["te"] & ~df["h3k9"], "TE+ H3K9-",
                  np.where(~df["te"] & df["h3k9"], "TE- H3K9+",
                                                   "TE- H3K9-")))
# require informative occupancy
df_info = df[df["mean_occ"] > 1e-3].copy()
df_info.to_csv(os.path.join(OUT, "rhimi59_h3k9_te_windows.tsv.gz"),
               sep="\t", index=False, compression="gzip")
print("  by category (informative windows):")
print(df_info["category"].value_counts())

# ----- summary -----
ORDER  = ["TE- H3K9-", "TE+ H3K9-", "TE- H3K9+", "TE+ H3K9+"]
COLORS = {"TE- H3K9-":"#7F7F7F",
          "TE+ H3K9-":"#1F77B4",
          "TE- H3K9+":"#FFA500",
          "TE+ H3K9+":"#D62728"}
summary = (df_info.groupby("category")
                  [["mean_occ","cv","ac_val","ac_lag"]]
                  .agg(["mean","median"])
                  .reindex(ORDER))
summary["n"] = df_info["category"].value_counts().reindex(ORDER)
summary.to_csv(os.path.join(OUT, "rhimi59_h3k9_te_summary.tsv"), sep="\t")
print(summary)

# ----- plot -----
print("[5] plotting")
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
plt.subplots_adjust(wspace=0.35)
features = [
    ("mean_occ", "Nucleosome occupancy (mean per 500bp)", None),
    ("cv",       "CV (sd/mean) — rugosity",               None),
    ("ac_val",   "Autocorrelation peak (140–240 bp)\nphasing strength", None),
    ("ac_lag",   "Autocorrelation peak position (bp)\nnucleosome repeat length", (140, 240)),
]
for ax, (col, ttl, ylim) in zip(axes, features):
    data = [df_info.loc[df_info["category"]==c, col].values for c in ORDER]
    bp = ax.boxplot(data, positions=np.arange(len(ORDER)), widths=0.65,
                    showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=1.5))
    for patch, c in zip(bp["boxes"], ORDER):
        patch.set_facecolor(COLORS[c]); patch.set_edgecolor("black")
    ax.set_xticks(np.arange(len(ORDER)))
    ax.set_xticklabels(ORDER, rotation=20, ha="right", fontsize=9)
    ax.set_title(ttl, fontsize=10)
    if ylim: ax.set_ylim(ylim)
    if col == "mean_occ":
        ax.set_yscale("log")
fig.suptitle("Rhizopus microsporus — nucleosome features by H3K9me3 × TE overlap (500 bp windows)",
             fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(os.path.join(OUT, "rhimi59_h3k9_te_features.pdf"))
fig.savefig(os.path.join(OUT, "rhimi59_h3k9_te_features.png"), dpi=180)
plt.close(fig)

# ---- focus plot: distribution of ac_val and occupancy in H3K9+ subset ----
sub = df_info[df_info["h3k9"]].copy()
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, col, xlabel in [(axes[0], "mean_occ", "Mean occupancy (log)"),
                         (axes[1], "ac_val",   "Autocorrelation peak (140–240 bp)")]:
    for tag, color in [(False, "#FFA500"), (True, "#D62728")]:
        d = sub.loc[sub["te"]==tag, col].values
        if col == "mean_occ":
            d = np.log10(d + 1e-3)
        ax.hist(d, bins=40, alpha=0.55, color=color,
                label=f"{'TE+' if tag else 'TE-'} (n={len(d)})",
                density=True)
    ax.set_xlabel(xlabel); ax.set_ylabel("density")
    ax.legend(fontsize=9)
fig.suptitle("Within H3K9me3 peaks: TE+ vs TE- windows", fontsize=12)
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(os.path.join(OUT, "rhimi59_h3k9_te_within_h3k9.pdf"))
fig.savefig(os.path.join(OUT, "rhimi59_h3k9_te_within_h3k9.png"), dpi=180)
plt.close(fig)
print("[DONE]")
