#!/usr/bin/env python3
"""
11_synteny_boxplot_by_clade.py — Boxplot of mean per-species synteny by clade.

Reads corrected (or raw) synteny ratios and produces a panel showing
the distribution of per-species mean synteny scores across clades.

Output: synteny_boxplot_by_clade.pdf
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────────────────
CORRECTED = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny_guide_corrected/"
                 "corrected_synteny_ratios.tsv")
PERSP     = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny_guide_corrected/"
                 "perspecies_corrected_synteny.tsv")
OUTDIR    = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny_guide_corrected/figures_heatmaps")

# Canonical clade palette (single source of truth)
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                   "..", "..", "Section2_TElandscape", "scripts"))
from clade_palette import CLADE_COLORS, CLADE_NAMES, CLADE_ORDER  # noqa: E402

NE_SPECIES = {
    "CH_Anasp1", "CH_Neocon1", "CH_Neolan1",
    "CH_Orpsp1", "CH_Pecora1", "CH_Caecom1", "CH_Pirfi3",
}


def get_clade(sp):
    if sp in NE_SPECIES:
        return "NE"
    return sp[:2]


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── Load pairwise corrected ratios ────────────────────────────────────
    df = pd.read_csv(CORRECTED, sep="\t")
    print(f"Loaded {len(df)} pairwise ratios across {df['clade'].nunique()} clades")

    # Remap NE species
    df["clade_fixed"] = df.apply(
        lambda r: "NE" if r["species_A"] in NE_SPECIES or r["species_B"] in NE_SPECIES
        else r["clade"], axis=1)

    # For pairs where BOTH species are NE, clade is NE
    # For pairs where one is NE and other is CH, clade stays CH (cross-clade)
    # Actually, let's re-derive: if BOTH species in NE_SPECIES -> NE
    # If one NE one CH -> keep as CH (it's a cross-group comparison)
    def fix_clade(row):
        a_ne = row["species_A"] in NE_SPECIES
        b_ne = row["species_B"] in NE_SPECIES
        if a_ne and b_ne:
            return "NE"
        return row["clade"]

    df["clade_fixed"] = df.apply(fix_clade, axis=1)

    # ── Compute per-species mean synteny within each clade ────────────────
    # For each species, compute its mean raw_score with all partners in the same clade
    records = []
    for clade, grp in df.groupby("clade_fixed"):
        species = set(grp["species_A"]) | set(grp["species_B"])
        for sp in species:
            mask = (grp["species_A"] == sp) | (grp["species_B"] == sp)
            sp_rows = grp[mask]
            if len(sp_rows) == 0:
                continue
            mean_raw = sp_rows["raw_score"].mean()
            mean_norm = sp_rows["r_ij_norm"].mean()
            mean_corr = sp_rows["r_corrected_norm"].mean()
            records.append({
                "clade": clade,
                "species": sp,
                "mean_raw": mean_raw,
                "mean_norm": mean_norm,
                "mean_corr": mean_corr,
                "n_pairs": len(sp_rows),
            })

    df_sp = pd.DataFrame(records)
    print(f"Per-species stats: {len(df_sp)} species across {df_sp['clade'].nunique()} clades")

    # ── Use canonical CLADE_ORDER (phylogenetic), keep clades with ≥1 species
    present_clades = set(df_sp["clade"].unique())
    clades_show = [c for c in CLADE_ORDER if c in present_clades and
                   df_sp[df_sp["clade"] == c]["species"].nunique() >= 1]
    extra = sorted(present_clades - set(CLADE_ORDER))
    clades_show += extra
    print(f"Showing boxplots for {len(clades_show)} clades (canonical order)")

    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 6))

    data_for_box = []
    positions = []
    colors = []
    labels = []

    for i, c in enumerate(clades_show):
        vals = df_sp[df_sp["clade"] == c]["mean_corr"].dropna().values
        data_for_box.append(vals)
        positions.append(i)
        colors.append(CLADE_COLORS.get(c, "#999"))
        n_sp = len(vals)
        labels.append(f"{CLADE_NAMES.get(c, c)}\n(n={n_sp})")

    bp = ax.boxplot(data_for_box, positions=positions, widths=0.6,
                    patch_artist=True, showfliers=True,
                    flierprops=dict(marker="o", markersize=3, alpha=0.5),
                    medianprops=dict(color="black", linewidth=1.5))

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Overlay individual points (jittered)
    rng = np.random.default_rng(42)
    for i, (vals, color) in enumerate(zip(data_for_box, colors)):
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(positions[i] + jitter, vals, s=12, alpha=0.6,
                   color=color, edgecolor="k", linewidth=0.3, zorder=3)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("Mean pairwise synteny (corrected score)", fontsize=10)
    ax.set_title("Synteny conservation by clade\n(per-species mean of within-clade corrected synteny scores)",
                 fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(0.5, color="#ccc", linestyle="--", linewidth=0.5, zorder=0)

    fig.tight_layout()

    out = OUTDIR / "synteny_boxplot_by_clade.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(str(out).replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[DONE] {out}")


if __name__ == "__main__":
    main()
