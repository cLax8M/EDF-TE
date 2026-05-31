#!/usr/bin/env python3
"""
C03_kimura_landscape_grid.py
============================
Expanded Kimura landscape panels — per-clade multi-page grid covering ALL
species with available EarlGrey divergence_summary_table.tsv.

For each species we build a "Kimura landscape" by binning TE families on
their mean_div (0–0.50, 1% bins), weighted by total_bp, stacked by TE class
(DNA / LTR / LINE / SINE / Rolling Circle / PLE / Unclassified).

Outputs (all PDF, multipage):
  C03_kimura_landscape_global.pdf  → all species sorted by recency_index
  C03_kimura_landscape_per_clade.pdf → one page per major clade, all species

Usage:
  python3 C03_kimura_landscape_grid.py [annotated_dir] [out_dir] [kimura_species_tsv]
"""

import sys
import os
import math
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Canonical clade palette + species display
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clade_palette import (CLADE_COLORS, CLADE_NAMES, get_clade,
                            species_display)  # noqa: E402

warnings.filterwarnings("ignore")

ANNOT = sys.argv[1] if len(sys.argv) > 1 else \
    "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures"
KIM_SP_TSV = sys.argv[3] if len(sys.argv) > 3 else \
    "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/C01_kimura_species.tsv"
os.makedirs(OUT_DIR, exist_ok=True)

# Non-TE classes filtered out
NON_TE_SF = frozenset({
    "Simple_repeat", "Low_complexity", "Satellite",
    "rRNA", "tRNA", "snRNA", "scRNA", "sno", "srpRNA",
    "ARTEFACT", "acro", "centr", "subtelo", "DNA_virus", "RNA", "Y",
})

# Canonical TE class colours (match cladogram legend)
CLASS_COLOURS = {
    "DNA":            "#1f77b4",   # blue
    "LTR":            "#d62728",   # red
    "LINE":           "#2ca02c",   # green
    "SINE":           "#ffcc33",   # gold/yellow
    "Rolling Circle": "#ff7f0e",   # orange
    "PLE":            "#5cb8a8",   # teal (Penelope)
    "Unclassified":   "#bdbdbd",   # gray
    "Other":          "#bdbdbd",
}
CLASS_ORDER = ["DNA", "LTR", "LINE", "SINE", "Rolling Circle", "PLE",
               "Unclassified", "Other"]

# Bin edges: 0–0.50 in 1% (0.01) bins → 50 bins
BINS = np.arange(0.0, 0.51, 0.01)
BIN_CENTERS = (BINS[:-1] + BINS[1:]) / 2


def normalise_class(subclass: str, superfamily: str) -> str:
    """Map EarlGrey (subclass, superfamily) → canonical class label."""
    if superfamily in NON_TE_SF:
        return None
    sc = (subclass or "").strip()
    sf = (superfamily or "").strip()
    if sc == "DNA":
        return "DNA"
    if sc == "LTR":
        return "LTR"
    if sc in ("LINE", "Non-LTR"):
        return "LINE"
    if sc == "SINE":
        return "SINE"
    if sc in ("RC", "Rolling Circle", "Helitron"):
        return "Rolling Circle"
    if sc == "PLE" or sf == "Penelope":
        return "PLE"
    if sc in ("Unknown", "Unclassified") or sf in ("Unknown", "Unclassified"):
        return "Unclassified"
    return "Other"


def load_species_landscape(sp: str):
    """Return dict {class: bp_per_bin (np.array len=50)}."""
    fp = os.path.join(ANNOT, f"{sp}_EarlGrey",
                      f"{sp}_summaryFiles",
                      f"{sp}_divergence_summary_table.tsv")
    if not os.path.isfile(fp):
        return None
    try:
        df = pd.read_csv(fp, sep="\t")
    except Exception:
        return None
    if df.empty or "mean_div" not in df.columns or "total_bp" not in df.columns:
        return None
    df["mean_div"] = pd.to_numeric(df["mean_div"], errors="coerce")
    df["total_bp"] = pd.to_numeric(df["total_bp"], errors="coerce")
    df = df.dropna(subset=["mean_div", "total_bp"])
    if df.empty:
        return None
    landscape = {c: np.zeros(len(BIN_CENTERS)) for c in CLASS_ORDER}
    for _, r in df.iterrows():
        cls = normalise_class(str(r.get("subclass", "")),
                              str(r.get("superfamily", "")))
        if cls is None:
            continue
        d = float(r["mean_div"])
        bp = float(r["total_bp"])
        if d < 0 or bp <= 0:
            continue
        # cap at 0.50
        d_capped = min(d, 0.499)
        idx = int(d_capped / 0.01)
        landscape[cls][idx] += bp
    return landscape


def panel_landscape(ax, landscape, title, max_y=None, show_x=True, show_y=True,
                    title_color="black"):
    bottom = np.zeros(len(BIN_CENTERS))
    for cls in CLASS_ORDER:
        if landscape[cls].sum() == 0:
            continue
        ax.bar(BIN_CENTERS, landscape[cls], width=0.01,
               bottom=bottom, color=CLASS_COLOURS[cls],
               edgecolor="none", align="center", label=cls)
        bottom += landscape[cls]
    ax.set_xlim(0, 0.5)
    if max_y:
        ax.set_ylim(0, max_y)
    ax.set_title(title, fontsize=10.5, pad=3, color=title_color,
                 fontweight="bold")
    if show_x:
        ax.set_xlabel("Kimura div.", fontsize=6)
    if show_y:
        ax.set_ylabel("bp", fontsize=6)
    ax.tick_params(axis="both", labelsize=5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def add_legend(fig):
    # Drop "Other" from the visible legend (still used internally for binning).
    legend_classes = [c for c in CLASS_ORDER if c != "Other"]
    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLOURS[c]) for c in legend_classes]
    fig.legend(handles, legend_classes, ncol=len(legend_classes), loc="upper center",
               bbox_to_anchor=(0.5, 0.995), fontsize=8, frameon=False)


def build_grid(species_list, landscapes, title, ncols=6,
               show_title=True, show_legend=True):
    n = len(species_list)
    nrows = math.ceil(n / ncols)
    fig_h = max(4.0, nrows * 1.7)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=(ncols * 2.4, fig_h),
                             squeeze=False)
    if show_title and title:
        fig.suptitle(title, fontsize=12, y=0.995)
    if show_legend:
        add_legend(fig)
    for idx, sp in enumerate(species_list):
        r, c = divmod(idx, ncols)
        ax = axes[r][c]
        ls = landscapes.get(sp)
        if ls is None:
            ax.axis("off")
            continue
        clade = get_clade(sp)
        title_color = CLADE_COLORS.get(clade, "black")
        sp_label = species_display(sp, italic=True)
        panel_landscape(ax, ls, sp_label,
                        show_x=(r == nrows - 1), show_y=(c == 0),
                        title_color=title_color)
    # hide empty
    for k in range(n, nrows * ncols):
        r, c = divmod(k, ncols)
        axes[r][c].axis("off")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


def main():
    # 1. Discover species
    species = sorted([
        d.replace("_EarlGrey", "")
        for d in os.listdir(ANNOT)
        if d.endswith("_EarlGrey")
        and os.path.isdir(os.path.join(ANNOT, d))
        and not d.startswith("VL_")   # drop control yeasts (user request 2026-04-21)
    ])
    print(f"[INFO] Discovered {len(species)} EG species directories")

    # 2. Load landscapes
    landscapes = {}
    skipped = 0
    for sp in species:
        ls = load_species_landscape(sp)
        if ls is None:
            skipped += 1
            continue
        landscapes[sp] = ls
    print(f"[INFO] Built landscapes for {len(landscapes)} species ({skipped} skipped)")

    # 3. Read C01_kimura_species.tsv to sort by recency_index
    sp_order = list(landscapes.keys())
    if os.path.isfile(KIM_SP_TSV):
        kim = pd.read_csv(KIM_SP_TSV, sep="\t")
        kim = kim.dropna(subset=["recency_index"])
        kim = kim.sort_values("recency_index", ascending=False)
        ordered = [s for s in kim["species"].tolist() if s in landscapes]
        # append any landscapes not in kim
        ordered += [s for s in landscapes if s not in ordered]
        sp_order = ordered
        print(f"[INFO] Sorted {len(ordered)} species by recency_index")

    # 4. Global grid (all species sorted by RI), one big multi-page PDF
    out_global = os.path.join(OUT_DIR, "C03_kimura_landscape_global.pdf")
    PER_PAGE = 48  # 6 cols × 8 rows
    with PdfPages(out_global) as pdf:
        for i in range(0, len(sp_order), PER_PAGE):
            chunk = sp_order[i:i + PER_PAGE]
            page_n = i // PER_PAGE + 1
            total_pages = math.ceil(len(sp_order) / PER_PAGE)
            # User request (2026-04-22): drop the "Kimura landscape …"
            # suptitle on every page, keep the class legend only on page 1.
            fig = build_grid(chunk, landscapes, title=None, ncols=6,
                             show_title=False,
                             show_legend=(page_n == 1))
            pdf.savefig(fig)
            plt.close(fig)
    print(f"[OK] {out_global}")

    # 5. Per-clade grid
    by_clade = defaultdict(list)
    for sp in sp_order:
        clade = get_clade(sp)
        by_clade[clade].append(sp)

    out_clade = os.path.join(OUT_DIR, "C03_kimura_landscape_per_clade.pdf")
    with PdfPages(out_clade) as pdf:
        for clade in sorted(by_clade.keys()):
            sps = by_clade[clade]
            ncols = 6 if len(sps) > 12 else (4 if len(sps) > 4 else len(sps))
            ncols = max(1, ncols)
            title = f"Clade {clade} — Kimura landscapes ({len(sps)} species)"
            fig = build_grid(sps, landscapes, title, ncols=ncols)
            pdf.savefig(fig)
            plt.close(fig)
    print(f"[OK] {out_clade}")


if __name__ == "__main__":
    main()
