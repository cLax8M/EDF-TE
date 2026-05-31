#!/usr/bin/env python3
"""03_te_genomic_correlations.py

Two figures for the TE × genomic-features correlation section:

  1. scatter_genome_size_vs_TE.pdf
     Scatter of genome size (Mb) vs total TE content (% of genome),
     points coloured by clade (canonical 18-clade palette), with linear
     regression line, Pearson r and p-value annotated.

  2. heatmap_te_vs_features_compact.pdf
     Compact PGLS-style heatmap restricted to:
       rows  = genome_size, intron_density, total_genes
       cols  = TE types (DNA, LTR, LINE, SINE, RC, PLE, Unknown)
     using a single-colour gradient (white → dark red), with significance
     stars (FDR < 0.05) read from D03_pgls/pgls_extended_results.tsv.

Inputs (defaults inside script):
  - merged TE+genomic TSV  (D01_te_genomics_merged_provisional.tsv)
  - PGLS results TSV       (pgls_extended_results.tsv)

Output dir defaults to results/GenomicFeatures/correlations/.
"""
from __future__ import annotations
from pathlib import Path
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

# Canonical clade palette
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "Section2_TElandscape", "scripts"))
from clade_palette import CLADE_COLORS, CLADE_NAMES, CLADE_ORDER, get_clade  # noqa: E402

MERGED = Path(os.environ.get("D01_MERGED",
              "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
              "D01_te_genomics_merged_provisional.tsv"))
PGLS   = Path(os.environ.get("D03_PGLS",
              "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
              "D03_pgls/pgls_extended_results.tsv"))
SYNTENY = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny_guide_corrected/"
               "perspecies_corrected_synteny.tsv")
OUTDIR = Path(os.environ.get("CORR_OUTDIR",
              "/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures/correlations"))


# ── helpers ───────────────────────────────────────────────────────────────────
# get_clade imported from clade_palette (handles NE_SPECIES override)


# ── 1. scatter genome-size vs TE content ──────────────────────────────────────
def scatter_size_vs_te(df: pd.DataFrame, outdir: Path):
    d = df[["species", "Genome_Size_bp", "pct_TE_total"]].dropna().copy()
    d["clade"]   = d["species"].map(get_clade)
    d["size_Mb"] = d["Genome_Size_bp"] / 1e6
    x = d["size_Mb"].values
    y = d["pct_TE_total"].values

    slope, intercept, r, p, _ = stats.linregress(x, y)
    # Also compute Spearman + log-x linear regression (better for skewed sizes)
    rho_s, p_s = stats.spearmanr(x, y)
    log_x = np.log10(x)
    slope_l, intercept_l, r_l, p_l, _ = stats.linregress(log_x, y)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.2))

    for ax, log_scale in zip(axes, [False, True]):
        for c in CLADE_ORDER:
            sub = d[d["clade"] == c]
            if sub.empty:
                continue
            ax.scatter(sub["size_Mb"], sub["pct_TE_total"],
                       s=42, color=CLADE_COLORS.get(c, "#999"),
                       edgecolor="k", linewidth=0.4, alpha=0.85,
                       label=CLADE_NAMES.get(c, c) if not log_scale else None)

        if log_scale:
            xs = np.logspace(np.log10(x.min()), np.log10(x.max()), 200)
            ax.plot(xs, intercept_l + slope_l * np.log10(xs),
                    color="#222", linewidth=1.4, linestyle="--",
                    label=f"y = {slope_l:.2f}·log₁₀(x) + {intercept_l:.2f}")
            txt = (f"Pearson (log x) r = {r_l:.3f}, p = {p_l:.2e}\n"
                   f"Spearman ρ = {rho_s:.3f}, p = {p_s:.2e}\n"
                   f"n = {len(d)}")
            ax.set_xscale("log")
            ax.set_title("Log-x scale", fontsize=11)
        else:
            xs = np.linspace(x.min(), x.max(), 100)
            ax.plot(xs, intercept + slope * xs,
                    color="#222", linewidth=1.4, linestyle="--",
                    label=f"y = {slope:.3g}·x + {intercept:.2f}")
            txt = (f"Pearson r = {r:.3f}, p = {p:.2e}\n"
                   f"Spearman ρ = {rho_s:.3f}, p = {p_s:.2e}\n"
                   f"n = {len(d)}")
            ax.set_title("Linear scale", fontsize=11)

        ax.text(0.03, 0.97, txt, transform=ax.transAxes,
                ha="left", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#999",
                          alpha=0.85))

        ax.set_xlabel("Genome size (Mb)", fontsize=11)
        ax.set_ylabel("Total TE content (% of genome)", fontsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(alpha=0.25, which="both")

    fig.suptitle(f"Genome size vs TE content ({len(d)} species)",
                 fontsize=13, y=1.0)

    # Single legend on the right side of the figure
    handles = []
    labels  = []
    for c in CLADE_ORDER:
        if not d[d["clade"] == c].empty:
            handles.append(plt.Line2D([0], [0], marker="o", linestyle="",
                                      markerfacecolor=CLADE_COLORS.get(c, "#999"),
                                      markeredgecolor="k", markeredgewidth=0.4,
                                      markersize=7))
            labels.append(CLADE_NAMES.get(c, c))
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=8, ncol=1, frameon=False, title="Clade")

    fig.tight_layout()
    out = outdir / "scatter_genome_size_vs_TE.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out}  (r={r:.3f}, p={p:.2e})")


# ── 2. compact heatmap (4 features × TE types) ────────────────────────────────
def compact_pgls_heatmap(outdir: Path, df_merged: pd.DataFrame):
    """Compact heatmap of TE class × genomic feature.

    Cells are coloured by the **Pearson correlation coefficient r derived
    from the PGLS test statistic** (r = t / √(t² + df)) so that the colour
    and the FDR-significance annotation come from the *same* model.
    A `*` on a cell marks PGLS FDR < 0.05.

    Adds a Synteny conservation column (mean genome-wide pairwise synteny r,
    from Synteny_guide_corrected/perspecies_corrected_synteny.tsv).
    """
    if not PGLS.is_file():
        print(f"[WARN] {PGLS} not found — skipping compact heatmap")
        return

    res = pd.read_csv(PGLS, sep="\t")
    sub = res[res["direction"] == "TE~genomic"].copy()
    # Recover signed Pearson r from the t-statistic and df:
    #   r = t / sqrt(t^2 + df),     df = n - 2
    sub["df"]    = sub["n"] - 2
    sub["r_pgls"] = sub["t_value"] / np.sqrt(sub["t_value"]**2 + sub["df"])

    # ── Add gene_density (= Total_Genes / Mb) as fresh tip-Pearson rows
    #    (no PGLS fit available) ----------------------------------------
    if {"Total_Genes", "Genome_Size_bp"}.issubset(df_merged.columns):
        df_merged = df_merged.copy()
        df_merged["gene_density"] = (df_merged["Total_Genes"]
                                       / (df_merged["Genome_Size_bp"] / 1e6))
        for r_te in ["pct_TE_total", "pct_DNA", "pct_LTR", "pct_LINE",
                      "pct_SINE", "pct_RC", "pct_PLE", "pct_Unknown"]:
            if r_te not in df_merged.columns:
                continue
            d = df_merged[[r_te, "gene_density"]].dropna()
            if len(d) < 5:
                continue
            r_v, p_v = stats.pearsonr(d[r_te].values, d["gene_density"].values)
            sub = pd.concat([sub, pd.DataFrame([{
                "response":  r_te, "predictor": "gene_density",
                "direction": "TE~genomic", "estimate": np.nan, "se": np.nan,
                "t_value":   np.nan, "p_raw": p_v, "n": len(d),
                "p_bonferroni": np.nan, "p_fdr": p_v,
                "sig_bonf": bool(p_v < 0.05/len(d)),
                "sig_fdr":  bool(p_v < 0.05),
                "df": len(d) - 2, "r_pgls": r_v,
            }])], ignore_index=True)

    # ── If we have a synteny TSV, run a fresh PGLS-equivalent OLS per
    #    TE class vs mean_r_corrected and merge the rows in.  We use a
    #    plain linear regression here because the synteny scores are
    #    already corrected for fragmentation.
    if SYNTENY.is_file():
        syn = pd.read_csv(SYNTENY, sep="\t")
        df_syn = df_merged.merge(syn[["species", "mean_r_corrected"]],
                                  on="species", how="inner")
        for r_te in ["pct_TE_total", "pct_DNA", "pct_LTR", "pct_LINE",
                      "pct_SINE", "pct_RC", "pct_PLE", "pct_Unknown"]:
            if r_te not in df_syn.columns:
                continue
            d = df_syn[[r_te, "mean_r_corrected"]].dropna()
            if len(d) < 5:
                continue
            r_v, p_v = stats.pearsonr(d[r_te].values,
                                       d["mean_r_corrected"].values)
            sub = pd.concat([sub, pd.DataFrame([{
                "response":   r_te,
                "predictor":  "mean_r_corrected",
                "direction":  "TE~genomic",
                "estimate":   np.nan,
                "se":         np.nan,
                "t_value":    np.nan,
                "p_raw":      p_v,
                "n":          len(d),
                "p_bonferroni": np.nan,
                "p_fdr":      p_v,            # provisional, no joint FDR
                "sig_bonf":   bool(p_v < 0.05/len(d)),
                "sig_fdr":    bool(p_v < 0.05),
                "df":         len(d) - 2,
                "r_pgls":     r_v,
            }])], ignore_index=True)

    # Rows (TE type) and columns (4 genomic features). pct_TE_total = total TE.
    te_rows = ["pct_TE_total",
               "pct_DNA", "pct_LTR", "pct_LINE", "pct_SINE",
               "pct_RC",  "pct_PLE", "pct_Unknown"]
    feat_cols = ["Genome_Size_bp", "Avg_Introns_Multiexon",
                 "gene_density", "mean_r_corrected"]
    feat_labels = {
        "Genome_Size_bp":         "Genome size",
        "Avg_Introns_Multiexon":  "Intron density",
        "gene_density":           "Gene density\n(genes / Mb)",
        "mean_r_corrected":       "Synteny conservation\n(mean corrected r)",
    }
    row_labels = {
        "pct_TE_total": "Total TE",
        "pct_DNA": "DNA",     "pct_LTR": "LTR",     "pct_LINE": "LINE",
        "pct_SINE": "SINE",   "pct_RC": "Rolling Circle",
        "pct_PLE": "Penelope","pct_Unknown": "Unclassified",
    }

    # ── Build matrix of Pearson r (from PGLS t-statistic) ─────────────────
    r_pivot = sub.pivot_table(index="response", columns="predictor",
                               values="r_pgls", aggfunc="first")
    r_mat = r_pivot.reindex(index=te_rows, columns=feat_cols)

    # FDR significance from PGLS (TE~genomic direction)
    sig = sub.pivot_table(index="response", columns="predictor",
                          values="sig_fdr", aggfunc="first")
    sig = sig.reindex(index=te_rows, columns=feat_cols).fillna(False)

    # n per cell for footnote
    n_pivot = sub.pivot_table(index="response", columns="predictor",
                               values="n", aggfunc="first")
    n_mat = n_pivot.reindex(index=te_rows, columns=feat_cols)

    # Diverging colour scale: GREEN = positive r, RED = negative r.
    # Saturated at |r|>=0.30 (user request 2026-04-21).
    vmax = 0.30
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "RedWhiteGreen",
        ["#b2182b", "#ef8a62", "#fddbc7", "#ffffff",
         "#d9f0d3", "#7fbf7b", "#1b7837"], N=256)

    fig, ax = plt.subplots(figsize=(2.4 + 1.4 * len(feat_cols),
                                    1.0 + 0.55 * len(te_rows)))
    im = ax.imshow(r_mat.values.astype(float), aspect="auto", cmap=cmap,
                   vmin=-vmax, vmax=vmax)

    # Tick labels
    ax.set_xticks(range(len(feat_cols)))
    ax.set_xticklabels([feat_labels[c] for c in feat_cols],
                       rotation=20, ha="right", fontsize=9)
    ax.set_yticks(range(len(te_rows)))
    ax.set_yticklabels([row_labels[r] for r in te_rows], fontsize=10)

    # Annotate each cell with r value (and a star if PGLS-FDR significant)
    for i, r in enumerate(te_rows):
        for j, c in enumerate(feat_cols):
            v = r_mat.loc[r, c]
            if not np.isfinite(v):
                ax.text(j, i, "·", ha="center", va="center",
                        fontsize=9, color="#999")
                continue
            star = "*" if sig.loc[r, c] else ""
            txt_color = "white" if abs(v) > 0.20 else "black"
            ax.text(j, i, f"{v:+.2f}{star}", ha="center", va="center",
                    fontsize=9, color=txt_color, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Pearson r derived from PGLS t-statistic\n"
                   "(r = t / √(t² + df))", fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title("TE class vs genomic feature — Pearson r from PGLS\n"
                 "(* = PGLS FDR < 0.05; synteny column uses ordinary Pearson r)",
                 fontsize=9.5)
    ax.set_xlabel("Genomic feature", fontsize=10)
    ax.set_ylabel("TE class", fontsize=10)

    fig.tight_layout()
    out = outdir / "heatmap_te_vs_features_compact.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out}")
    # Also dump the underlying r/p table for the manuscript supplements
    sub.to_csv(outdir / "heatmap_te_vs_features_compact.tsv",
               sep="\t", index=False)


def raw_pearson_heatmap(outdir: Path, df_merged: pd.DataFrame):
    """Sister of compact_pgls_heatmap, but cells coloured by **raw Pearson r
    on tip data** (no phylogenetic correction). Significance star = raw
    p < 0.05 after BH-FDR across all cells of this matrix."""
    te_rows = ["pct_TE_total",
               "pct_DNA", "pct_LTR", "pct_LINE", "pct_SINE",
               "pct_RC",  "pct_PLE", "pct_Unknown"]
    feat_cols = ["Genome_Size_bp", "Avg_Introns_Multiexon",
                 "gene_density", "mean_r_corrected"]
    feat_labels = {
        "Genome_Size_bp":         "Genome size",
        "Avg_Introns_Multiexon":  "Intron density",
        "gene_density":           "Gene density\n(genes / Mb)",
        "mean_r_corrected":       "Synteny conservation\n(mean corrected r)",
    }
    row_labels = {
        "pct_TE_total": "Total TE",
        "pct_DNA": "DNA",     "pct_LTR": "LTR",     "pct_LINE": "LINE",
        "pct_SINE": "SINE",   "pct_RC": "Rolling Circle",
        "pct_PLE": "Penelope","pct_Unknown": "Unclassified",
    }

    df = df_merged.copy()
    # gene_density (genes / Mb)
    if {"Total_Genes", "Genome_Size_bp"}.issubset(df.columns):
        df["gene_density"] = df["Total_Genes"] / (df["Genome_Size_bp"] / 1e6)
    if SYNTENY.is_file():
        syn = pd.read_csv(SYNTENY, sep="\t")
        df = df.merge(syn[["species", "mean_r_corrected"]],
                       on="species", how="left")

    r_mat = pd.DataFrame(index=te_rows, columns=feat_cols, dtype=float)
    p_mat = pd.DataFrame(index=te_rows, columns=feat_cols, dtype=float)
    n_mat = pd.DataFrame(index=te_rows, columns=feat_cols, dtype=float)
    for r_te in te_rows:
        if r_te not in df.columns: continue
        for c_ft in feat_cols:
            if c_ft not in df.columns: continue
            d = df[[r_te, c_ft]].dropna()
            if len(d) < 5: continue
            r_v, p_v = stats.pearsonr(d[r_te].values, d[c_ft].values)
            r_mat.loc[r_te, c_ft] = r_v
            p_mat.loc[r_te, c_ft] = p_v
            n_mat.loc[r_te, c_ft] = len(d)

    # BH-FDR on the flat list of finite p-values
    pflat = p_mat.values.flatten()
    finite = np.isfinite(pflat)
    sig = np.zeros_like(pflat, dtype=bool)
    if finite.any():
        # BH-FDR (manual, no statsmodels dep)
        idx_fin = np.where(finite)[0]
        ps = pflat[idx_fin]
        order = np.argsort(ps)
        m = len(ps)
        ranked = ps[order]
        thresh = (np.arange(1, m + 1) / m) * 0.05
        passed = ranked <= thresh
        if passed.any():
            kmax = np.max(np.where(passed)[0])
            keep = order[:kmax + 1]
            sig[idx_fin[keep]] = True
    sig_mat = pd.DataFrame(sig.reshape(p_mat.shape),
                           index=p_mat.index, columns=p_mat.columns)

    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "RedWhiteGreen",
        ["#b2182b", "#ef8a62", "#fddbc7", "#ffffff",
         "#d9f0d3", "#7fbf7b", "#1b7837"], N=256)

    fig, ax = plt.subplots(figsize=(2.4 + 1.4 * len(feat_cols),
                                    1.0 + 0.55 * len(te_rows)))
    im = ax.imshow(r_mat.values.astype(float), aspect="auto", cmap=cmap,
                   vmin=-0.30, vmax=0.30)

    ax.set_xticks(range(len(feat_cols)))
    ax.set_xticklabels([feat_labels[c] for c in feat_cols],
                       rotation=20, ha="right", fontsize=9)
    ax.set_yticks(range(len(te_rows)))
    ax.set_yticklabels([row_labels[r] for r in te_rows], fontsize=10)

    for i, r in enumerate(te_rows):
        for j, c in enumerate(feat_cols):
            v = r_mat.loc[r, c]
            if not np.isfinite(v):
                ax.text(j, i, "·", ha="center", va="center",
                        fontsize=9, color="#999"); continue
            star = "*" if sig_mat.loc[r, c] else ""
            txt_color = "white" if abs(v) > 0.20 else "black"
            ax.text(j, i, f"{v:+.2f}{star}", ha="center", va="center",
                    fontsize=9, color=txt_color, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Raw Pearson r (no phylogenetic correction)", fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title("TE class vs genomic feature — RAW Pearson r on tip data\n"
                 "(* = BH-FDR < 0.05 across this matrix; "
                 "compare with PGLS-derived heatmap)",
                 fontsize=9.5)
    ax.set_xlabel("Genomic feature", fontsize=10)
    ax.set_ylabel("TE class", fontsize=10)

    fig.tight_layout()
    out = outdir / "heatmap_te_vs_features_raw_pearson.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out}")


EDF_CLADES = {"AP","BL","BS","CH","CR","EN","GL","KX","MI","MO","MU","NE","OL","ZO"}


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(MERGED, sep="\t")
    # Drop control yeasts (VL clade) globally — user request 2026-04-21
    if "species" in df.columns:
        n_before = len(df)
        df = df[~df["species"].str.startswith("VL_")].reset_index(drop=True)
        if len(df) != n_before:
            print(f"[filter] dropped {n_before - len(df)} VL species")
    print(f"Loaded {len(df)} species from {MERGED.name}")
    scatter_size_vs_te(df, OUTDIR)
    compact_pgls_heatmap(OUTDIR, df)
    raw_pearson_heatmap(OUTDIR, df)

    # ── EDF-only versions (early-diverging fungi) ────────────────────
    edf_outdir = OUTDIR / "EDF_only"
    edf_outdir.mkdir(parents=True, exist_ok=True)
    df_edf = df[df["species"].str[:2].isin(EDF_CLADES)].reset_index(drop=True)
    print(f"[EDF] {len(df_edf)} early-diverging fungi species")
    if len(df_edf) >= 10:
        compact_pgls_heatmap(edf_outdir, df_edf)
        raw_pearson_heatmap(edf_outdir, df_edf)


if __name__ == "__main__":
    main()
