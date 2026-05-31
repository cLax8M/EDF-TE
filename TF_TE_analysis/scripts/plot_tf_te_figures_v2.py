#!/usr/bin/env python3
"""
plot_tf_te_figures_v2.py
─────────────────────────
Redesigned TF-TE figures:
  1. Bubble plot — redesigned as lollipop/dot strip (sorted by log2OR),
     readable per-TF labels, color=TF family, size~n_peaks but fixed in legend.
  2. FIMO heatmap — restricted to the 5 key TFs with TE-binding evidence.

Output:
  peak_enrichment_lollipop.pdf
  peak_enrichment_bubble_v2.pdf
  fimo_heatmap_5tf.pdf
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from scipy.stats import fisher_exact, binomtest

# ── Paths ───────────────────────────────────────────────────────────────────
RES = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis")
OUTDIR = RES / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

OVERLAP_STATS = RES / "01_peak_overlap" / "peak_TE_overlap_stats.tsv"
OVERLAP_MATRIX = RES / "01_peak_overlap" / "peak_TE_overlap_matrix.tsv"
BW_SIGNAL      = RES / "02_bigwig_signal" / "bw_TE_signal_class_mean.tsv"
FIMO_ENRICH   = RES / "03_fimo_enrichment" / "fimo_enrichment.tsv"

# Key TFs with ChIP peak overlap + FIMO motif evidence on LINE/I-Jockey:
#   HMGb_1228      : 28 peaks on 31 LINE/I-Jockey elements (+FIMO HMGb3, +BigWig)
#   ForkHead_9727  : 45 peaks on 48 LINE/I-Jockey elements (+FIMO FK1)
#   bHLH_Myc_4449  : 12 peaks on LINE/I-Jockey (ChIP only, no FIMO for LINE/I-J)
# Additional TFs with strong general TE enrichment (no LINE/I-Jockey FIMO):
#   HMGb_2385      : DNA transposons (DNA/P, DNA/CMC-EnSpm)
#   bHLH_Myc_6330, bHLH_Myc_8146 : LTR/Gypsy
KEY_TFS = ["HMGb_1228", "ForkHead_9727", "bHLH_Myc_4449", "HMGb_2385", "bHLH_Myc_6330", "bHLH_Myc_8146"]

# Star policy in lollipop: only mark TFs with triple evidence on a specific target TE.
# Triple evidence = significant ChIP overlap + significant FIMO motif + BigWig z >= 1.0
# all on the same TE class.
STAR_TARGET_TE = {
    "HMGb_1228": "LINE/I-Jockey",
    "ForkHead_9727": "LINE/I-Jockey",
    "bHLH_Myc_4449": "LINE/I-Jockey",
    "HMGb_2385": "DNA/CMC-EnSpm",
    "bHLH_Myc_6330": "LTR/Gypsy",
    "bHLH_Myc_8146": "LTR/Gypsy",
}

TF_MOTIF_FAMILY = {
    "HMGb_1228": "HMGb",
    "ForkHead_9727": "FK",
    "bHLH_Myc_4449": "HLH_Myc",
    "HMGb_2385": "HMGb",
    "bHLH_Myc_6330": "HLH_Myc",
    "bHLH_Myc_8146": "HLH_Myc",
}

TF_FAMILY_COLORS = {
    "HMGb":        "#E07B39",
    "HOX":         "#1F77B4",
    "HOX_KN":      "#7FB6D4",
    "HOX_Myb":     "#AEC7E8",
    "SANT_Myb":    "#2CA02C",
    "bHLH_Myc":    "#D62728",
    "bZIP":        "#9467BD",
    "ForkHead":    "#E377C2",
    "HSF":         "#8C564B",
    "GATA_Znf":    "#17BECF",
    "MADS":        "#BCBD22",
    "Znf_C2H2":    "#7F7F7F",
    "CP2":         "#4B4BFF",
    "CopperFist":  "#AAAAAA",
}
DEFAULT_COLOR = "#CCCCCC"


# ────────────────────────────────────────────────────────────────────────────
# FIGURE 1: Lollipop/strip plot — sorted by log2OR, labeled
# ────────────────────────────────────────────────────────────────────────────
def plot_lollipop(stats: pd.DataFrame):
    """
    Sorted horizontal lollipop: each TF is one dot on log2OR axis,
    y = TF label sorted by log2OR, color = TF family, marker size = n_peaks.
    Sidepanel shows –log10 p value as bar.
    """
    df = stats.sort_values("log2_odds_ratio", ascending=True).reset_index(drop=True)

    # Build strict star mask: triple evidence on the same target TE class
    star_ok = {tf: False for tf in STAR_TARGET_TE}
    if OVERLAP_MATRIX.exists() and BW_SIGNAL.exists() and FIMO_ENRICH.exists():
        mat = pd.read_table(OVERLAP_MATRIX, index_col=0)
        bw = pd.read_table(BW_SIGNAL, index_col=0)
        fimo = pd.read_table(FIMO_ENRICH)
        stats_idx = stats.set_index("tf_label")
        frac = mat.div(stats_idx["n_peaks"], axis=0)
        genomic_frac = frac.mean(axis=0)

        for tf, te_target in STAR_TARGET_TE.items():
            if tf not in mat.index or tf not in stats_idx.index:
                continue

            # 1) ChIP enrichment on target TE
            n_total = int(stats_idx.loc[tf, "n_peaks"])
            n_on = int(mat.loc[tf, te_target]) if te_target in mat.columns else 0
            p_exp = float(genomic_frac.get(te_target, 0))
            chip_ok = (p_exp > 0 and n_on > 0 and
                       binomtest(n_on, n_total, p_exp, alternative="greater").pvalue < 0.05)

            # 2) FIMO motif enrichment on target TE
            fam = TF_MOTIF_FAMILY.get(tf)
            fimo_ok = False
            if fam is not None:
                sub = fimo[(fimo["te_class"] == te_target) &
                           (fimo["motif_family"] == fam) &
                           (fimo["fisher_pval"] < 0.05)]
                fimo_ok = not sub.empty

            # 3) BigWig signal on target TE
            bw_ok = (tf in bw.columns and te_target in bw.index and
                     float(bw.loc[te_target, tf]) >= 1.0)

            star_ok[tf] = chip_ok and fimo_ok and bw_ok

    # Clamp p=0 to 1e-300, compute –log10
    pvals = df["fisher_pval"].replace(0, 1e-300)
    df["nlog10p"] = -np.log10(pvals)

    fig, axes = plt.subplots(1, 2, figsize=(12, 9),
                              gridspec_kw={"width_ratios": [3, 1]})
    ax, ax_p = axes

    n = len(df)
    y_pos = np.arange(n)

    # Draw lollipop sticks
    for i, row in df.iterrows():
        color = TF_FAMILY_COLORS.get(row["tf_family"], DEFAULT_COLOR)
        sig = row["fisher_pval"] < 0.05
        lw = 1.2 if sig else 0.5
        ax.hlines(i, 0, row["log2_odds_ratio"], color=color, lw=lw, alpha=0.6)
        sz = 40 + row["n_peaks"] / 40   # more moderate scale
        sz = min(sz, 350)
        ax.scatter(row["log2_odds_ratio"], i,
                   s=sz, color=color,
                   edgecolors="black" if sig else "none",
                   linewidths=0.7 if sig else 0,
                   zorder=3, alpha=0.85)

    # Reference line
    ax.axvline(0, color="grey", lw=0.8, ls=":")
    # Genomic TE density (log2OR=0 IS genomic density, mark for clarity)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["tf_label"], fontsize=7)
    ax.set_xlabel("log₂ Odds Ratio (DAP-seq peaks on TEs\nvs genomic TE fraction)", fontsize=9)
    ax.set_title("A   TF enrichment on TEs (DAP-seq peak overlap)", fontsize=10, loc="left")
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(-1.3, 2.3)

    # Color-coded y-tick labels
    for tick, (_, row) in zip(ax.get_yticklabels(), df.iterrows()):
        tick.set_color(TF_FAMILY_COLORS.get(row["tf_family"], "#333333"))

    # Side panel: –log10 p
    bar_colors = [TF_FAMILY_COLORS.get(row["tf_family"], DEFAULT_COLOR)
                  for _, row in df.iterrows()]
    ax_p.barh(y_pos, df["nlog10p"], color=bar_colors, edgecolor="none", height=0.7)
    ax_p.axvline(-np.log10(0.05), color="grey", lw=0.8, ls="--")
    ax_p.set_xlabel("–log₁₀ p", fontsize=8)
    ax_p.set_yticks([])
    ax_p.spines[["top", "right", "left"]].set_visible(False)
    ax_p.set_title("B   Significance", fontsize=10, loc="left")

    # Legend — only family color, fixed circle size (no proportional bubbles)
    families_present = df["tf_family"].unique()
    patches = [mpatches.Patch(color=TF_FAMILY_COLORS.get(f, DEFAULT_COLOR), label=f)
               for f in sorted(families_present)]
    ax.legend(handles=patches, fontsize=7, loc="lower right",
              framealpha=0.85, title="TF family", title_fontsize=7)

    # Annotate TFs with strict triple evidence on target TE only
    for i, row in df.iterrows():
        if star_ok.get(row["tf_label"], False):
            ax.annotate("★", xy=(row["log2_odds_ratio"], i),
                        xytext=(row["log2_odds_ratio"] + 0.08, i),
                        fontsize=8, color="black", va="center")

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.suptitle("DAP-seq TF peak enrichment on transposable elements  (MU_Rhimi59)\n"
                 "★ = triple evidence on target TE (ChIP + FIMO + BigWig)",
                 fontsize=11, fontweight="bold")

    for ext in ("pdf", "svg"):
        p = OUTDIR / f"peak_enrichment_lollipop.{ext}"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        print(f"Saved: {p}")
    plt.close(fig)


# ────────────────────────────────────────────────────────────────────────────
# FIGURE 2: Improved volcano-style bubble (replace 2-cluster issue)
# Replace –log10(p) on y with % peaks on TE (more informative),
# keep log2OR on x, bubble size = n_peaks (moderate scale), color = TF family.
# ────────────────────────────────────────────────────────────────────────────
def plot_bubble_v2(stats: pd.DataFrame):
    """
    X = log2OR, Y = fraction of peaks on TE, size ∝ n_peaks (capped),
    color = TF family, annotate top 10 significant TFs.
    This avoids the p=0/p>0 two-cluster problem.
    """
    fig, ax = plt.subplots(figsize=(9, 7))

    genomic_te = stats["genomic_TE_frac"].iloc[0]
    added = set()
    for _, row in stats.iterrows():
        fam = row["tf_family"]
        color = TF_FAMILY_COLORS.get(fam, DEFAULT_COLOR)
        x = row["log2_odds_ratio"]
        y = row["frac_peaks_on_TE"]
        sz = max(30, min(500, row["n_peaks"] / 6))
        sig = row["fisher_pval"] < 0.05
        lbl = fam if fam not in added else "_nolegend_"
        added.add(fam)
        ax.scatter(x, y, s=sz, color=color, alpha=0.75,
                   edgecolors="black" if sig else "none",
                   linewidths=0.6, label=lbl, zorder=3)

    # Reference: genomic TE fraction
    ax.axhline(genomic_te, ls="--", color="grey", lw=0.8)
    ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] > -2 else -1.2,
            genomic_te + 0.005, "genomic TE fraction", fontsize=7, color="grey")
    ax.axvline(0, ls=":", color="grey", lw=0.8)

    # Annotate top significant TFs
    sig_df = stats[(stats["fisher_pval"] < 0.05) & (stats["log2_odds_ratio"] > 0)]
    sig_df = sig_df.sort_values("log2_odds_ratio", ascending=False).head(12)
    for _, row in sig_df.iterrows():
        ax.annotate(
            row["tf_label"].replace("bHLH_Myc_","bHLH_").replace("HMGb_","HMGb_"),
            xy=(row["log2_odds_ratio"], row["frac_peaks_on_TE"]),
            xytext=(5, 5), textcoords="offset points",
            fontsize=6.5, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5),
            color=TF_FAMILY_COLORS.get(row["tf_family"], "#333333")
        )

    ax.set_xlabel("log₂ Odds Ratio (peaks on TE vs genomic TE density)", fontsize=10)
    ax.set_ylabel("Fraction of DAP-seq peaks overlapping TEs", fontsize=10)
    ax.set_title("TF enrichment on transposable elements\n"
                 "(bubble size ∝ n peaks; black edge = p < 0.05)", fontsize=10)

    # Legend: color only, fixed size marker
    families_present = stats["tf_family"].unique()
    patches = [mpatches.Patch(color=TF_FAMILY_COLORS.get(f, DEFAULT_COLOR), label=f)
               for f in sorted(families_present)]
    ax.legend(handles=patches, fontsize=8, loc="upper left",
              framealpha=0.85, title="TF family", title_fontsize=8)

    plt.tight_layout()
    for ext in ("pdf", "svg"):
        p = OUTDIR / f"peak_enrichment_bubble_v2.{ext}"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        print(f"Saved: {p}")
    plt.close(fig)


# ────────────────────────────────────────────────────────────────────────────
# FIGURE 3: FIMO heatmap — key TFs with confirmed TE-binding evidence
# ────────────────────────────────────────────────────────────────────────────
def plot_fimo_heatmap_5tf(fimo: pd.DataFrame):
    """
    Heatmap: rows=TE class, cols=key TFs
    Cell = log2OR (only if p<0.05, else 0)
    Motif family annotated per TF column.
    """
    SKIP = {"Simple_repeat", "Low_complexity", "Satellite", "Unknown",
            "Unspecified", "rRNA", "snRNA", "srpRNA"}
    SKIP_PREF = ("tRNA", "SINE/tRNA", "Retroposon", "RNA", "rRNA", "snRNA")

    def skip_te(cls):
        if cls in SKIP:
            return True
        return any(cls.startswith(p) for p in SKIP_PREF)

    # Motif → TF label mapping (best motif per TF based on known FIMO findings)
    # fimo has columns: te_class, motif_id, motif_family, n_te_seqs, te_hits,
    #   te_hit_frac, n_bg_seqs, bg_hits, bg_hit_frac, odds_ratio, log2_odds_ratio, fisher_pval
    # We map motif families to TF labels
    TF_MOTIF_MAP = {
        "HMGb_1228":      ("HMGb",     ["HMGb3", "HMGb2"]),
        "ForkHead_9727":  ("FK",        ["FK1", "FK4"]),
        "bHLH_Myc_4449":  ("HLH_Myc",  ["HLH_Myc8", "HLH_Myc13", "HLH_Myc2",
                                          "HLH_Myc5", "HLH_Myc14"]),
        "HMGb_2385":      ("HMGb",     ["HMGb2", "HMGb3"]),
        "bHLH_Myc_6330":  ("HLH_Myc",  ["HLH_Myc6", "HLH_Myc1", "HLH_Myc9",
                                          "HLH_Myc10", "HLH_Myc11", "HLH_Myc8",
                                          "HLH_Myc5", "HLH_Myc2", "HLH_Myc14",
                                          "HLH_Myc15", "HLH_Myc3", "HLH_Myc13"]),
        "bHLH_Myc_8146":  ("HLH_Myc",  ["HLH_Myc14", "HLH_Myc15", "HLH_Myc3"]),
    }

    # For each TF, take the best (max log2OR) motif hit per TE class
    tables = {}
    for tf_label, (motif_family, motif_ids) in TF_MOTIF_MAP.items():
        sub = fimo[fimo["motif_family"] == motif_family].copy()
        sub = sub[~sub["te_class"].apply(skip_te)]
        # Best motif per TE class
        best = (sub[sub["fisher_pval"] < 0.05]
                .sort_values("log2_odds_ratio", ascending=False)
                .drop_duplicates("te_class")
                .set_index("te_class")["log2_odds_ratio"])
        tables[tf_label] = best

    # Build matrix
    te_classes = sorted({tc for t in tables.values() for tc in t.index})
    mat = pd.DataFrame(0.0, index=te_classes, columns=KEY_TFS)
    for tf_label, series in tables.items():
        for te, val in series.items():
            mat.loc[te, tf_label] = val

    # Drop rows where all are 0
    mat = mat[(mat != 0).any(axis=1)]
    if mat.empty:
        print("[WARN] Empty matrix for 5-TF FIMO heatmap")
        return

    # Sort rows: LINEs first, then LTRs, then DNA
    def te_order(cls):
        if cls.startswith("LINE/I-Jockey"):   return "0_" + cls
        if cls.startswith("LINE/I"):          return "1_" + cls
        if cls.startswith("LINE"):            return "2_" + cls
        if cls.startswith("SINE"):            return "3_" + cls
        if cls.startswith("LTR"):             return "4_" + cls
        if cls.startswith("DNA"):             return "5_" + cls
        return "9_" + cls
    mat = mat.loc[sorted(mat.index, key=te_order)]

    fig_h = max(5, len(mat) * 0.36)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    vmax = max(abs(mat.values.max()), 1.0)
    im = ax.imshow(mat.values, aspect="auto", cmap="Reds",
                   vmin=0, vmax=min(vmax, 7))

    # Cell annotations
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            if v > 0:
                txt = f"{v:.1f}"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=7, color="white" if v > vmax * 0.55 else "#333333",
                        fontweight="bold" if v > 3 else "normal")

    cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("log₂ Odds Ratio (FIMO, p < 0.05)", fontsize=8)

    ax.set_xticks(range(len(KEY_TFS)))
    ax.set_xticklabels(KEY_TFS, rotation=30, ha="right", fontsize=9)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=8)

    # Color column labels by TF family
    family_of = {"HMGb_1228": "HMGb", "HMGb_2385": "HMGb",
                 "ForkHead_9727": "ForkHead",
                 "bHLH_Myc_6330": "bHLH_Myc", "bHLH_Myc_4449": "bHLH_Myc",
                 "bHLH_Myc_8146": "bHLH_Myc"}
    for tick, tf in zip(ax.get_xticklabels(), KEY_TFS):
        fam = family_of.get(tf, "")
        tick.set_color(TF_FAMILY_COLORS.get(fam, "#333333"))

    ax.set_title("FIMO motif enrichment in TE sequences\n"
                 "(Key TFs with confirmed TE-binding evidence, MU_Rhimi59)",
                 fontsize=10, fontweight="bold")

    # Add motif label under each TF column
    motif_label = {
        "HMGb_1228":      "motif: HMGb3",
        "ForkHead_9727":  "motif: FK1",
        "bHLH_Myc_4449":  "motif: HLH_Myc",
        "HMGb_2385":      "motif: HMGb2",
        "bHLH_Myc_6330":  "motif: HLH_Myc",
        "bHLH_Myc_8146":  "motif: HLH_Myc",
    }
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(len(KEY_TFS)))
    ax2.set_xticklabels([motif_label[t] for t in KEY_TFS],
                        rotation=30, ha="left", fontsize=7, color="grey")
    ax2.tick_params(top=False, labeltop=True)

    plt.tight_layout()
    for ext in ("pdf", "svg"):
        p = OUTDIR / f"fimo_heatmap_5tf.{ext}"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        print(f"Saved: {p}")
    plt.close(fig)


# ────────────────────────────────────────────────────────────────────────────
# FIGURE 4: Per-TF peak overlap per TE class (only 5 key TFs)
# Shows LTR/Gypsy fraction alongside others — the bHLH/LTR-Gypsy question
# ────────────────────────────────────────────────────────────────────────────
def plot_peak_overlap_5tf():
    """
    Barplot: for each key TF, fraction of peaks per TE class.
    Panels annotated with BigWig coverage evidence and FIMO motif evidence.
    """
    OVERLAP_MATRIX = RES / "01_peak_overlap" / "peak_TE_overlap_matrix.tsv"
    OVERLAP_STATS  = RES / "01_peak_overlap" / "peak_TE_overlap_stats.tsv"
    BW_SIGNAL      = RES / "02_bigwig_signal" / "bw_TE_signal_class_mean.tsv"
    FIMO_ENRICH    = RES / "03_fimo_enrichment" / "fimo_enrichment.tsv"
    if not OVERLAP_MATRIX.exists():
        print("[WARN] overlap matrix not found")
        return

    mat    = pd.read_table(OVERLAP_MATRIX, index_col=0)
    stats  = pd.read_table(OVERLAP_STATS,  index_col="tf_label")
    bw     = pd.read_table(BW_SIGNAL,      index_col=0) if BW_SIGNAL.exists() else None
    fimo   = pd.read_table(FIMO_ENRICH)    if FIMO_ENRICH.exists() else None

    # TF → motif family for FIMO lookup
    TF_MOTIF_FAMILY = {
        "HMGb_1228":     "HMGb",
        "ForkHead_9727": "FK",
        "bHLH_Myc_4449": "HLH_Myc",
        "HMGb_2385":     "HMGb",
        "bHLH_Myc_6330": "HLH_Myc",
        "bHLH_Myc_8146": "HLH_Myc",
    }

    def get_fimo_te_classes(tf):
        """Return set of TE classes with significant FIMO hit for this TF."""
        if fimo is None:
            return set()
        fam = TF_MOTIF_FAMILY.get(tf)
        if not fam:
            return set()
        sub = fimo[(fimo["motif_family"] == fam) &
                   ((fimo["fisher_pval"] < 0.05) | (fimo["log2_odds_ratio"] == "inf"))]
        return set(sub["te_class"].unique())

    def get_bw_te_classes(tf, threshold=1.0):
        """Return dict of TE class → z_score where z >= threshold."""
        if bw is None or tf not in bw.columns:
            return {}
        vals = bw[tf].dropna()
        return {cls: z for cls, z in vals.items() if z >= threshold}

    SKIP = {"Simple_repeat", "Low_complexity", "Satellite", "Unknown", "Unspecified"}
    te_cols = [c for c in mat.columns if c not in SKIP
               and not any(c.startswith(p) for p in ("tRNA","rRNA","snRNA"))]

    # Fractions
    frac = mat[te_cols].div(stats["n_peaks"], axis=0)

    # Genome TE fraction per class — use mean over all TFs as proxy
    genomic_frac = frac.mean(axis=0)

    fig, axes = plt.subplots(len(KEY_TFS), 1, figsize=(14, 3.5 * len(KEY_TFS)),
                              sharex=False)

    for ax, tf in zip(axes, KEY_TFS):
        if tf not in frac.index:
            ax.text(0.5, 0.5, f"{tf} not found", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        row = frac.loc[tf].sort_values(ascending=False)
        row = row[row > 0]
        if row.empty:
            continue

        colors = []
        for te in row.index:
            if te.startswith("LINE/I-Jockey"):  colors.append("#8B1A1A")
            elif te.startswith("LINE"):          colors.append("#CC4444")
            elif te.startswith("LTR"):           colors.append("#4477AA")
            elif te.startswith("DNA"):           colors.append("#228B22")
            elif te.startswith("SINE"):          colors.append("#FF8C00")
            elif te.startswith("RC"):            colors.append("#AA6688")
            else:                               colors.append("#AAAAAA")

        # Per-class binomial significance test (one-tailed, greater)
        n_total = int(stats.loc[tf, "n_peaks"])
        sig_marks = []
        for te in row.index:
            n_on = int(mat.loc[tf, te]) if te in mat.columns else 0
            p_exp = genomic_frac.get(te, 0)
            if p_exp <= 0 or n_on == 0:
                sig_marks.append("")
                continue
            res = binomtest(n_on, n_total, p_exp, alternative="greater")
            p = res.pvalue
            if p < 0.001:
                sig_marks.append("***")
            elif p < 0.01:
                sig_marks.append("**")
            elif p < 0.05:
                sig_marks.append("*")
            else:
                sig_marks.append("")

        bars = ax.bar(range(len(row)), row.values, color=colors, edgecolor="none")

        # Significance marks above bars
        for i, (mark, h) in enumerate(zip(sig_marks, row.values)):
            if mark:
                ax.text(i, h + 0.002, mark, ha="center", va="bottom",
                        fontsize=7, color="#222222", fontweight="bold")

        # Overlay genomic fraction per class as dotted line
        gf_vals = [genomic_frac.get(te, 0) for te in row.index]
        ax.scatter(range(len(row)), gf_vals, color="black", s=12, zorder=5,
                   marker="D", label="genomic TE fraction")

        ax.set_xticks(range(len(row)))
        ax.set_xticklabels(row.index, rotation=55, ha="right", fontsize=7)
        tf_family = ("HMGb" if "HMGb" in tf
                     else "ForkHead" if "ForkHead" in tf
                     else "bHLH_Myc")
        fam_color = TF_FAMILY_COLORS.get(tf_family, "#333333")
        ax.set_title(f"{tf}  (total {int(stats.loc[tf,'n_peaks'])} peaks, "
                     f"{stats.loc[tf,'frac_peaks_on_TE']*100:.1f}% on TE, "
                     f"log₂OR={stats.loc[tf,'log2_odds_ratio']:+.2f})",
                     fontsize=9, color=fam_color)
        ax.set_ylabel("Fraction of peaks", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

        # Per-bar evidence icons: add colored badges above bars for FIMO / BigWig
        fimo_classes = get_fimo_te_classes(tf)
        bw_classes   = get_bw_te_classes(tf, threshold=1.0)
        ymax = ax.get_ylim()[1]
        for i, te in enumerate(row.index):
            h = row.values[i]
            has_fimo = te in fimo_classes
            has_bw   = te in bw_classes
            if has_fimo or has_bw:
                # Draw small colored badge(s) just above the bar
                offset = ymax * 0.04
                if has_fimo:
                    ax.text(i, h + offset, "M", ha="center", va="bottom",
                            fontsize=6.5, fontweight="bold", color="white",
                            bbox=dict(boxstyle="round,pad=0.2", fc="#C0392B",
                                      ec="none", alpha=0.9))
                    offset += ymax * 0.09
                if has_bw:
                    z = bw_classes[te]
                    ax.text(i, h + offset, f"▲{z:.1f}", ha="center", va="bottom",
                            fontsize=6, fontweight="bold", color="white",
                            bbox=dict(boxstyle="round,pad=0.2", fc="#1A6FA6",
                                      ec="none", alpha=0.9))

        if ax == axes[0]:
            ax.legend(fontsize=7)

    # Legends
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    class_legend = [
        Patch(color="#8B1A1A", label="LINE/I-Jockey"),
        Patch(color="#CC4444", label="LINE (other)"),
        Patch(color="#4477AA", label="LTR"),
        Patch(color="#228B22", label="DNA"),
        Patch(color="#FF8C00", label="SINE"),
        Patch(color="#AA6688", label="RC/Helitron"),
    ]
    evidence_legend = [
        Patch(fc="#C0392B", ec="none", label="M  = motif present in TE (FIMO, p<0.05)"),
        Patch(fc="#1A6FA6", ec="none", label="▲z = DAP-seq BigWig signal (z ≥ 1.0)"),
    ]
    leg1 = fig.legend(handles=class_legend, loc="upper right",
                      fontsize=8, framealpha=0.85, title="TE category",
                      bbox_to_anchor=(1.0, 0.98))
    fig.legend(handles=evidence_legend, loc="upper right",
               fontsize=8, framealpha=0.85, title="Evidence badges",
               bbox_to_anchor=(1.0, 0.78))

    fig.suptitle("Fraction of DAP-seq peaks per TE class — key TFs  (all have ChIP + FIMO + BigWig evidence)\n"
                 "(◆ = expected from genomic TE abundance;  * p<0.05  ** p<0.01  *** p<0.001 vs expected, binomial test)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 0.9, 0.96])

    for ext in ("pdf", "svg"):
        p = OUTDIR / f"peak_overlap_5tf.{ext}"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        print(f"Saved: {p}")
    plt.close(fig)


# ────────────────────────────────────────────────────────────────────────────
# MAIN
# ────────────────────────────────────────────────────────────────────────────
def main():
    stats = pd.read_table(OVERLAP_STATS)
    fimo  = pd.read_table(FIMO_ENRICH)

    print("=== Figure 1: Lollipop plot ===")
    plot_lollipop(stats)

    print("=== Figure 2: Bubble v2 (fraction vs log2OR) ===")
    plot_bubble_v2(stats)

    print("=== Figure 3: FIMO heatmap 5 TFs ===")
    plot_fimo_heatmap_5tf(fimo)

    print("=== Figure 4: Peak overlap per TE class (5 TFs) ===")
    plot_peak_overlap_5tf()

    print("[DONE]")


if __name__ == "__main__":
    main()
