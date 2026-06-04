#!/usr/bin/env python3
"""
fig_chimte_basalfungi.py
ChimeraTE chimera counts for basal fungi species.

Reads merged *_final.tsv files from per-species chimeraTE_out directories.
Produces a stacked-bar figure with one bar per species, grouped by clade,
and coloured by chimera type.

Outputs:
  figChimTE_basalfungi.pdf / .png
"""

from __future__ import annotations

import os
import glob
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ─── Paths ───────────────────────────────────────────────────────────────────
CHIMTE_DIR = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/chimerateTE"
CLADE_FILE = "/legserv/Temp/Carlos/EDF-TE/results/HTT/clade_assignment.tsv"
OUT_DIR    = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"

# Chimera types and display names (order = stacking order)
CHIMERA_TYPES = [
    ("TE-initiated",           "TE-initiated"),
    ("TE-terminated",          "TE-terminated"),
    ("TE-exonized_embedded",   "TE-exonized (embedded)"),
    ("TE-exonized_overlapped", "TE-exonized (overlapped)"),
    ("TE-exonized_intronic",   "TE-exonized (intronic)"),
]

# Colours for chimera types
CHIMERA_COLORS = {
    "TE-initiated":           "#d62728",
    "TE-terminated":          "#ff7f0e",
    "TE-exonized_embedded":   "#2ca02c",
    "TE-exonized_overlapped": "#1f77b4",
    "TE-exonized_intronic":   "#9467bd",
}

CLADE_COLORS = {
    "AP": "#5c9e3a", "AS": "#e8887a", "BA": "#4472c4",
    "BL": "#7030a0", "BS": "#d0b4e0", "CH": "#e0d040",
    "CR": "#6aafe6", "EN": "#e01010", "GL": "#8b5e48",
    "HO": "#b4d0e8", "KX": "#f0a820", "MI": "#3aab98",
    "MO": "#f4a8b0", "MU": "#e07840", "NU": "#1a2850",
    "OL": "#3668b0", "ZO": "#c8dc70",
}

CLADE_LABELS = {
    "AP": "Aphelida",    "AS": "Ascomycota",      "BA": "Basidiomycota",
    "BL": "Blasto-\ncladiomycota", "BS": "Basidiobolo-\nmycota",
    "CH": "Chytridiomycota", "CR": "Cryptomycota", "EN": "Entomophthoromycota",
    "GL": "Glomeromycota", "KX": "Kickxellomycotina", "MI": "Microsporidia",
    "MO": "Mortierellomycota", "MU": "Mucoromycotina",
    "NU": "Nuclearida", "OL": "Olpidiomycota", "ZO": "Zoopagomycota",
}

# Clade display order (basal fungi first)
CLADE_ORDER = ["MU", "MO", "GL", "KX", "ZO", "EN", "BS", "CH", "MI",
               "BL", "OL", "CR", "AP", "AS", "BA", "HO", "NU"]


def load_species_clades(clade_file: str) -> dict[str, str]:
    df = pd.read_csv(clade_file, sep="\t", names=["species", "clade"], header=0)
    return dict(zip(df["species"], df["clade"]))


def count_unique_genes(final_tsv: str) -> int:
    """Return number of unique gene_id rows in a *_final.tsv file."""
    if not os.path.exists(final_tsv) or os.path.getsize(final_tsv) == 0:
        return 0
    try:
        df = pd.read_csv(final_tsv, sep="\t")
        if "gene_id" in df.columns:
            return df["gene_id"].nunique()
        return len(df)
    except Exception:
        return 0


def collect_counts(chimte_dir: str) -> pd.DataFrame:
    """Collect per-species, per-chimera-type unique gene counts."""
    records = []
    for sp_dir in sorted(glob.glob(os.path.join(chimte_dir, "*"))):
        sp = os.path.basename(sp_dir)
        if sp.startswith("@"):
            continue
        out_dir = os.path.join(sp_dir, "projects", "chimeraTE_out")
        has_any = False
        row = {"species": sp}
        for ctype, _ in CHIMERA_TYPES:
            final = os.path.join(out_dir, f"{ctype}_final.tsv")
            n = count_unique_genes(final)
            row[ctype] = n
            if n > 0:
                has_any = True
        if has_any:
            records.append(row)
    return pd.DataFrame(records)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chimte-dir", default=CHIMTE_DIR)
    ap.add_argument("--clade-file", default=CLADE_FILE)
    ap.add_argument("--out-dir", default=OUT_DIR)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    sp2clade = load_species_clades(args.clade_file)
    counts = collect_counts(args.chimte_dir)

    if counts.empty:
        print("[WARN] No species with final ChimeraTE TSVs found.")
        return

    print(f"[INFO] Species with data: {len(counts)}")

    # Assign clade and sort
    counts["clade"] = counts["species"].map(sp2clade).fillna("?")
    counts["clade_order"] = counts["clade"].apply(
        lambda c: CLADE_ORDER.index(c) if c in CLADE_ORDER else 999
    )
    counts = counts.sort_values(["clade_order", "species"]).reset_index(drop=True)

    ctypes = [ct for ct, _ in CHIMERA_TYPES]
    species = counts["species"].tolist()
    clades  = counts["clade"].tolist()

    # ── Normalise to 100% ────────────────────────────────────────────────────
    totals = counts[ctypes].sum(axis=1).values.astype(float)
    totals[totals == 0] = 1  # avoid div-by-zero

    # ── Figure ────────────────────────────────────────────────────────────────
    fig_width = max(8, len(species) * 0.55 + 3)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    x = np.arange(len(species))
    bottom = np.zeros(len(species))

    for ctype, label in CHIMERA_TYPES:
        vals = counts[ctype].values.astype(float) / totals * 100
        bars = ax.bar(x, vals, bottom=bottom, color=CHIMERA_COLORS[ctype],
                      label=label, edgecolor="white", linewidth=0.4)
        bottom += vals

    # ── Total count labels on top of each bar ────────────────────────────────
    for i, tot in enumerate(totals.astype(int)):
        ax.text(x[i], 101, f"n={tot}", ha="center", va="bottom",
                fontsize=5.5, color="#444444", rotation=90)

    # ── Species labels (tilted) ───────────────────────────────────────────────
    ax.set_xticks(x)
    ax.set_xticklabels(species, rotation=45, ha="right", fontsize=7)

    # ── Color the x-tick labels by clade ────────────────────────────────────
    for tick, sp, clade in zip(ax.get_xticklabels(), species, clades):
        tick.set_color(CLADE_COLORS.get(clade, "#333333"))

    # ── Clade separators and labels ──────────────────────────────────────────
    prev_clade = None
    clade_start = {}
    clade_end = {}
    for i, (sp, cl) in enumerate(zip(species, clades)):
        if cl != prev_clade:
            clade_start[cl] = i
            if prev_clade is not None:
                clade_end[prev_clade] = i - 1
                ax.axvline(x=i - 0.5, color="#999999", linewidth=0.7, linestyle="--")
            prev_clade = cl
    if prev_clade:
        clade_end[prev_clade] = len(species) - 1

    ax.set_ylim(0, 100)
    clade_label_y = 119  # above the n= annotations
    for cl, start in clade_start.items():
        end = clade_end.get(cl, start)
        mid = (start + end) / 2.0
        label = CLADE_LABELS.get(cl, cl)
        ax.text(mid, clade_label_y, label, ha="center", va="bottom",
                fontsize=6.5, color=CLADE_COLORS.get(cl, "#333333"),
                rotation=0, clip_on=False)

    ax.set_ylabel("% chimeric genes", fontsize=9)
    ax.set_title("ChimeraTE chimeric genes across basal fungi (proportion)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Legend ───────────────────────────────────────────────────────────────
    handles = [
        mpatches.Patch(facecolor=CHIMERA_COLORS[ct], label=lbl)
        for ct, lbl in CHIMERA_TYPES
    ]
    ax.legend(handles=handles, fontsize=7, frameon=False,
              loc="upper right", bbox_to_anchor=(1.0, 1.0))

    plt.tight_layout()

    for ext in ("pdf", "png"):
        out = os.path.join(args.out_dir, f"figChimTE_basalfungi.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"[OK] Saved {out}")

    plt.close(fig)

    # ── Save TSV summary ─────────────────────────────────────────────────────
    tsv_out = os.path.join(args.out_dir, "figChimTE_basalfungi_counts.tsv")
    out_cols = ["species", "clade"] + ctypes
    counts[out_cols].to_csv(tsv_out, sep="\t", index=False)
    print(f"[OK] Saved {tsv_out}")


if __name__ == "__main__":
    main()
