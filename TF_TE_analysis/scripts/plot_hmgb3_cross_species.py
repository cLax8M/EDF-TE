#!/usr/bin/env python3
"""
plot_hmgb3_cross_species.py
───────────────────────────
Figure: Fraction of LINE/I-Jockey elements carrying HMGb3 binding site
        across all fungal species with EarlGrey annotations.

Reads: cross_species_fimo_results.tsv
Output:
  - hmgb3_jockey_cross_species.pdf   (bar chart, species sorted by fraction)
  - hmgb3_jockey_cross_species.svg   (same, vector)
  - hmgb_line_heatmap.pdf            (heatmap: HMGb2/3 × all LINE classes)

Usage
  python plot_hmgb3_cross_species.py \\
      --input /path/to/cross_species_fimo_results.tsv \\
      --outdir /path/to/figures/
"""

import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Clade color map (by 2-letter prefix) ───────────────────────────────────
CLADE_COLORS = {
    "MU": "#E07B39",   # Mucoromycota — orange
    "BL": "#9C5FA8",   # Blastocladiomycota — purple
    "AP": "#59A14F",   # Aphelidiomycota — green
    "MI": "#76B7B2",   # Microsporidia — teal
    "CH": "#F1CE63",   # Chytridiomycota — yellow
    "BA": "#4E79A7",   # Basidiomycota — blue
    "AS": "#E15759",   # Ascomycota — red
    "NU": "#B07AA1",   # Neocallimastigomycota — mauve
    "MO": "#FF9DA7",   # Mortierellomycota — pink
    "GL": "#9C755F",   # Glomeromycotina — brown
}
CLADE_LABELS = {
    "MU": "Mucoromycota", "BL": "Blastocladiomycota",
    "AP": "Aphelidiomycota", "MI": "Microsporidia",
    "CH": "Chytridiomycota", "BA": "Basidiomycota",
    "AS": "Ascomycota", "NU": "Neocallimastigomycota",
    "MO": "Mortierellomycota", "GL": "Glomeromycotina",
}
DEFAULT_COLOR = "#AAAAAA"

# Species to highlight (label on bars)
HIGHLIGHT = {
    "MU_Rhimi59":  "Rhizopus microsporus",
    "MU_Mucci3":   "Mucor circinelloides",
    "MU_Muclu1":   "Mucor lusitanicus",
    "MU_Rhide1":   "Rhizopus delemar",
    "MU_Licmin1":  "Lichtheimia corymbifera",
    "MU_Abscae1":  "Absidia corymbifera",
}


def clade_of(species: str) -> str:
    prefix = species[:2].upper()
    return prefix if prefix in CLADE_COLORS else "XX"


def plot_bar_jockey(df: pd.DataFrame, outdir: Path):
    """Bar chart: fraction of LINE/I-Jockey with HMGb3 per species."""
    sub = df[(df["te_class"] == "LINE/I-Jockey") & (df["motif"] == "HMGb3")].copy()
    if sub.empty:
        print("[WARN] No LINE/I-Jockey × HMGb3 data")
        return

    # Remove species with < 5 elements (fractions from n<5 are unreliable)
    MIN_N = 5
    sub = sub[sub["n_elements"] >= MIN_N].copy()
    sub = sub.sort_values("fraction_with_motif", ascending=True)

    colors = [CLADE_COLORS.get(clade_of(sp), DEFAULT_COLOR)
              for sp in sub["species"]]

    fig_height = max(6, len(sub) * 0.22)
    fig, ax = plt.subplots(figsize=(9, fig_height))

    bars = ax.barh(range(len(sub)), sub["fraction_with_motif"],
                   color=colors, edgecolor="none", height=0.8)

    # Label highlighted species
    for i, (_, row) in enumerate(sub.iterrows()):
        if row["species"] in HIGHLIGHT:
            label = HIGHLIGHT[row["species"]]
            ax.text(row["fraction_with_motif"] + 0.01, i,
                    f" {label}",
                    va="center", ha="left", fontsize=6.5, style="italic",
                    color="#333333")

    ax.set_yticks([])
    ax.set_xlabel("Fraction of LINE/I-Jockey elements with HMGb3 binding site",
                  fontsize=10)
    ax.set_title("HMGb3 motif in LINE/I-Jockey retrotransposons\nacross fungi",
                 fontsize=11, fontweight="bold")
    ax.set_xlim(0, 1.05)
    ax.axvline(0.5, color="grey", linewidth=0.7, linestyle="--", alpha=0.5)

    # Add n_elements annotation for top 5 (by fraction)
    top5 = sub.nlargest(5, "fraction_with_motif")
    for i, (_, row) in enumerate(sub.iterrows()):
        if row["species"] in top5["species"].values:
            ax.text(0.01, i,
                    f"n={int(row['n_elements'])}",
                    va="center", ha="left", fontsize=5.5, color="white",
                    fontweight="bold")

    # Legend (clades present in data)
    clades_present = sorted({clade_of(sp) for sp in sub["species"]
                              if clade_of(sp) in CLADE_COLORS})
    patches = [mpatches.Patch(color=CLADE_COLORS[c],
                               label=CLADE_LABELS.get(c, c))
               for c in clades_present]
    ax.legend(handles=patches, loc="lower right", fontsize=7,
              framealpha=0.9, title="Clade", title_fontsize=7)

    plt.tight_layout()
    for ext in ("pdf", "svg"):
        path = outdir / f"hmgb3_jockey_cross_species.{ext}"
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {path}")
    plt.close(fig)


def plot_hmgb_line_heatmap(df: pd.DataFrame, outdir: Path):
    """
    Heatmap: average fraction of LINE elements with HMGb2/3
    grouped by clade and LINE subclass.
    """
    hmgb = df[df["motif"].isin(["HMGb2", "HMGb3"]) &
              df["te_class"].str.startswith("LINE")].copy()
    if hmgb.empty:
        print("[WARN] No HMGb × LINE data for heatmap")
        return

    hmgb["clade"] = hmgb["species"].apply(lambda s: clade_of(s))
    pivot = (hmgb.groupby(["clade", "te_class", "motif"])["fraction_with_motif"]
             .mean()
             .unstack(["te_class", "motif"])
             .fillna(0))

    if pivot.empty:
        return

    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 0.9),
                                    max(4, pivot.shape[0] * 0.6)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd",
                   vmin=0, vmax=1)
    plt.colorbar(im, ax=ax, label="Mean fraction with motif")

    col_labels = [f"{te}/{mot}" for te, mot in pivot.columns]
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(col_labels, rotation=55, ha="right", fontsize=7)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels([CLADE_LABELS.get(c, c) for c in pivot.index], fontsize=8)
    ax.set_title("HMGb motif enrichment in LINE elements by clade",
                 fontsize=10, fontweight="bold")

    plt.tight_layout()
    path = outdir / "hmgb_line_heatmap.pdf"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    print(f"  Saved: {path}")
    plt.close(fig)


def plot_all_te_classes_heatmap(df: pd.DataFrame, outdir: Path):
    """
    Heatmap: average fraction per (clade × TE class) for the top motif of
    each TE class (max fraction per te_class across all motifs).
    """
    # Best motif per TE class
    best = (df.groupby(["te_class", "motif"])["fraction_with_motif"]
              .mean()
              .reset_index()
              .sort_values("fraction_with_motif", ascending=False)
              .drop_duplicates("te_class"))
    te_motif_map = {row["te_class"]: row["motif"] for _, row in best.iterrows()}

    rows = []
    for te_class, motif in te_motif_map.items():
        sub = df[(df["te_class"] == te_class) & (df["motif"] == motif)].copy()
        sub["clade"] = sub["species"].apply(clade_of)
        grp = sub.groupby("clade")["fraction_with_motif"].mean()
        for clade, frac in grp.items():
            rows.append({"clade": clade, "te_class": te_class,
                         "motif": motif, "fraction": frac})
    if not rows:
        return

    wide = (pd.DataFrame(rows)
            .pivot(index="clade", columns="te_class", values="fraction")
            .fillna(0))

    fig, ax = plt.subplots(figsize=(max(12, wide.shape[1] * 0.8),
                                    max(4, wide.shape[0] * 0.6)))
    im = ax.imshow(wide.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    plt.colorbar(im, ax=ax, label="Mean fraction with best motif")
    ax.set_xticks(range(wide.shape[1]))
    ax.set_xticklabels(wide.columns, rotation=55, ha="right", fontsize=7)
    ax.set_yticks(range(wide.shape[0]))
    ax.set_yticklabels([CLADE_LABELS.get(c, c) for c in wide.index], fontsize=8)
    ax.set_title("Best TF motif enrichment per TE class across clades",
                 fontsize=10, fontweight="bold")
    plt.tight_layout()
    path = outdir / "te_class_motif_clade_heatmap.pdf"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    print(f"  Saved: {path}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  required=True,
                    help="cross_species_fimo_results.tsv")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for figures")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_table(args.input)
    print(f"Loaded {len(df)} rows from {args.input}")
    print(f"Species: {df['species'].nunique()}, "
          f"TE classes: {df['te_class'].nunique()}, "
          f"Motifs: {df['motif'].nunique()}")

    plot_bar_jockey(df, out)
    plot_hmgb_line_heatmap(df, out)
    plot_all_te_classes_heatmap(df, out)

    print("[DONE]")


if __name__ == "__main__":
    main()
