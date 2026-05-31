#!/usr/bin/env python3
"""plot_nucleosome_te_v2.py — Nucleosome-TE interaction figure (2-panel).

Matches previous figure style:
  A. Per-clade Δ NuclPred score (TE insertion sites − genome background)
     Horizontal bar chart with whiskers (mean ± SD across species)
  B. Per-TE-class Δ NuclPred score (averaged across all species)
     Horizontal bar chart with whiskers

Data: F3_bias per-species TSVs (205 species).
Six species use real Nuc-seq data instead of NuclPred predictions:
  BA_Mrafri1, BA_Mereb1, AS_Usnflo1, BL_Catan2, BA_Rhoto, MU_Rhimi59
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

F3_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
              "nucleosomes/f3_bias")
OUTDIR = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures")

CLADE_COLORS = {
    "MU": "#E07B39", "BL": "#9C5FA8", "AP": "#59A14F",
    "MI": "#76B7B2", "CH": "#F1CE63", "BA": "#4E79A7",
    "AS": "#E15759", "NU": "#B07AA1", "MO": "#FF9DA7",
    "GL": "#9C755F", "HO": "#80b1d3", "CR": "#a6cee3",
    "OL": "#1f78b4", "ZO": "#b2df8a", "EN": "#e31a1c",
    "KX": "#fdbf6f", "BS": "#cab2d6",
}
CLADE_NAMES = {
    "MU": "Mucoromycota", "BL": "Blastocladiomycota",
    "AP": "Aphelidiomycota", "MI": "Microsporidia",
    "CH": "Chytridiomycota", "BA": "Basidiomycota",
    "AS": "Ascomycota", "NU": "Nucleariida",
    "MO": "Mortierellomycota", "GL": "Glomeromycotina",
    "HO": "Holozoa", "CR": "Cryptomycota",
    "OL": "Olpidiomycota", "ZO": "Zoopagomycota",
    "EN": "Entomophthoromycota", "KX": "Kickxellomycota",
    "BS": "Basidiobolomycota",
}

# TE class grouping: broad classes + keep tRNA, Simple, Low, Satellite separate
TE_BROAD_MAP = {
    "LTR": "LTR", "DNA": "DNA", "LINE": "LINE", "SINE": "SINE",
    "RC": "RC", "Unknown": "Unknown",
    "Satellite": "Satellite",
    # Excluded: tRNA, Simple_repeat, Low_complexity (not real TEs)
}

TE_CLASS_COLORS = {
    "LTR": "#e41a1c", "DNA": "#377eb8", "LINE": "#4daf4a",
    "SINE": "#f1c232", "RC": "#984ea3", "Unknown": "#bdbdbd",
    "Satellite": "#cab2d6",
}


def to_broad(sf):
    """Map superfamily to broad class."""
    if sf in TE_BROAD_MAP:
        return TE_BROAD_MAP[sf]
    top = sf.split("/")[0]
    return TE_BROAD_MAP.get(top, None)  # None → skip (e.g. rRNA)


def load_f3():
    """Load all per-species F3 bias TSVs."""
    frames = []
    for f in sorted(F3_DIR.glob("*_f3_bias.tsv")):
        if f.stat().st_size > 0:
            try:
                frames.append(pd.read_csv(f, sep="\t"))
            except Exception:
                pass
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = load_f3()
    print(f"[INFO] Loaded {len(df)} rows, {df['species'].nunique()} species")

    df["clade"] = df["species"].str[:2]
    df["broad_class"] = df["superfamily"].apply(to_broad)
    df = df.dropna(subset=["broad_class"])
    # Only keep entries with at least 5 TE insertions (avoid noisy estimates)
    df = df[df["n_insertions"] >= 5]

    # Normalize delta: real nuc-seq species have occupancy in read-count scale
    # (background ~200-5000) vs NuclPred [0,1]. Normalize delta by background
    # to get fractional change: delta_norm = delta / |background|
    # This makes both scales comparable.
    REAL_NUC_SP = {"BA_Mrafri1", "BA_Mereb1", "AS_Usnflo1",
                   "BL_Catan2", "BA_Rhoto", "MU_Rhimi59"}
    df["is_real_nuc"] = df["species"].isin(REAL_NUC_SP)
    df["delta_norm"] = df["delta_occupancy"] / df["background_occupancy"].abs().clip(lower=0.01)
    # For NuclPred species, keep original delta (already [0,1] scale)
    # For real nuc species, use normalized delta
    df["delta_plot"] = np.where(df["is_real_nuc"], df["delta_norm"], df["delta_occupancy"])

    n_real = df[df["is_real_nuc"]]["species"].nunique()
    print(f"[INFO] After filtering: {len(df)} rows")
    print(f"[INFO] {n_real} species with real Nuc-seq data (normalized)")
    print(f"[INFO] {df['species'].nunique() - n_real} species with NuclPred predictions")

    # ── Per-species weighted mean delta ───────────────────────────────────
    # For Panel A: aggregate by species (all superfamilies), then group by clade
    sp_delta = (df.groupby("species")
                .apply(lambda g: np.average(g["delta_plot"],
                                            weights=g["n_insertions"]))
                .reset_index(name="delta"))
    sp_delta["clade"] = sp_delta["species"].str[:2]

    # For Panel B: aggregate by species × broad_class, then average across species
    cls_delta = (df.groupby(["species", "broad_class"])
                 .apply(lambda g: pd.Series({
                     "delta": np.average(g["delta_plot"],
                                         weights=g["n_insertions"]),
                     "n": g["n_insertions"].sum()
                 }))
                 .reset_index())

    # ── Figure ────────────────────────────────────────────────────────────
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(16, 7),
                                      gridspec_kw={"width_ratios": [1.1, 1],
                                                    "wspace": 0.35})

    # ── Panel A: Per-clade delta (horizontal bars) ────────────────────────
    # Sort clades by mean delta
    clade_stats = (sp_delta.groupby("clade")["delta"]
                   .agg(["mean", "std", "count"])
                   .sort_values("mean", ascending=True))  # most negative at top

    clades = clade_stats.index.tolist()
    n_c = len(clades)
    y_pos = np.arange(n_c)

    for i, c in enumerate(clades):
        val = clade_stats.loc[c, "mean"]
        sd = clade_stats.loc[c, "std"]
        n = clade_stats.loc[c, "count"]
        color = CLADE_COLORS.get(c, "#aaa")
        # Draw bar
        ax_a.barh(i, val, height=0.7, color=color, edgecolor="none", alpha=0.85)
        # Error bar (SD)
        if not np.isnan(sd) and n > 1:
            ax_a.errorbar(val, i, xerr=sd, fmt="none", ecolor="black",
                         capsize=3, capthick=1, elinewidth=1)

    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels([CLADE_NAMES.get(c, c) for c in clades], fontsize=9)
    ax_a.axvline(0, color="black", linewidth=1, linestyle="--", alpha=0.7)
    ax_a.set_xlabel("Δ NuclPred score\n(TE insertion sites − genome background)", fontsize=10)
    ax_a.set_title("A   Nucleosome Occupancy at TE Insertions\nvs Genome Background",
                   fontsize=11, fontweight="bold", loc="left")
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    # Annotation box
    ax_a.text(0.55, 0.10,
              "Δ < 0 → TEs insert into\nnucleosome-depleted regions (NDRs)\n"
              "Δ > 0 → TEs insert into\nnucleosome-enriched regions",
              transform=ax_a.transAxes, fontsize=7,
              bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                       edgecolor="gray", alpha=0.8),
              va="bottom", ha="left")
    # Limit x-axis to avoid Olpidiomycota (1 species) pulling the scale
    xlim = max(abs(clade_stats["mean"].quantile(0.05)),
               abs(clade_stats["mean"].quantile(0.95))) * 1.5
    xlim = max(xlim, 0.2)
    ax_a.set_xlim(-xlim, xlim)

    # ── Panel B: Per-TE-class delta (horizontal bars) ─────────────────────
    # Average across species for each class
    class_stats = (cls_delta.groupby("broad_class")["delta"]
                   .agg(["mean", "std", "count"])
                   .sort_values("mean", ascending=True))

    classes = class_stats.index.tolist()
    n_cl = len(classes)
    y_pos_b = np.arange(n_cl)

    for i, cl in enumerate(classes):
        val = class_stats.loc[cl, "mean"]
        sd = class_stats.loc[cl, "std"]
        n = class_stats.loc[cl, "count"]
        color = TE_CLASS_COLORS.get(cl, "#aaa")
        ax_b.barh(i, val, height=0.7, color=color, edgecolor="none", alpha=0.85)
        if not np.isnan(sd) and n > 1:
            ax_b.errorbar(val, i, xerr=sd, fmt="none", ecolor="black",
                         capsize=3, capthick=1, elinewidth=1)

    ax_b.set_yticks(y_pos_b)
    ax_b.set_yticklabels(classes, fontsize=9)
    ax_b.axvline(0, color="black", linewidth=1, linestyle="--", alpha=0.7)
    ax_b.set_xlabel("Δ NuclPred score (TE − background)", fontsize=10)
    ax_b.set_title("B   Nucleosome Bias by TE Class\n(averaged across all species)",
                   fontsize=11, fontweight="bold", loc="left")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    fig.tight_layout()
    out = OUTDIR / "nucleosome_te_v2.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "nucleosome_te_v2.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[DONE] {out}")

    # ── Summary ───────────────────────────────────────────────────────────
    real_nuc_species = {"BA_Mrafri1", "BA_Mereb1", "AS_Usnflo1",
                        "BL_Catan2", "BA_Rhoto", "MU_Rhimi59"}
    n_real = sum(1 for sp in sp_delta["species"] if sp in real_nuc_species)
    print(f"\nSpecies with real Nuc-seq data: {n_real}")
    print(f"Species with NuclPred predictions: {sp_delta['species'].nunique() - n_real}")

    print(f"\n{'='*60}")
    print("Panel A — Per-clade Δ NuclPred (sorted):")
    for c in clades:
        row = clade_stats.loc[c]
        print(f"  {CLADE_NAMES.get(c,c):25s}  mean={row['mean']:+.4f}  "
              f"sd={row['std']:.4f}  n={int(row['count'])}")

    print(f"\nPanel B — Per-TE-class Δ NuclPred (sorted):")
    for cl in classes:
        row = class_stats.loc[cl]
        print(f"  {cl:15s}  mean={row['mean']:+.4f}  "
              f"sd={row['std']:.4f}  n={int(row['count'])}")


if __name__ == "__main__":
    main()
