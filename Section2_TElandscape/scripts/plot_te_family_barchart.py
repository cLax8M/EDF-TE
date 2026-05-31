#!/usr/bin/env python3
"""plot_te_family_barchart.py — TF-style bar chart for TE families.

3-panel figure analogous to the TF families figure:
  A. Collapsed clade view with stacked bars (mean # TE families by class) + Shannon
  B. Cladogram with gain/loss bubbles at key internal nodes
  C. Full species tree with stacked bars (# TE families by class) + Shannon dots

Requires: ete3 for tree manipulation, matplotlib for plotting.
"""
from pathlib import Path
import glob
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from collections import defaultdict, Counter

# ── Paths ─────────────────────────────────────────────────────────────────
PA_MATRIX = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                 "ancestral_count/04_pa_matrix_v2/te_family_pa_matrix.count")
CLUSTER_TSV = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                   "ancestral_count/03_clusters_v2/cluster_cluster.tsv")
GAINS_LOSSES = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                    "ancestral_count/05_count_v2b/gains_losses_summary.tsv")
TREE_FILE = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                 "ancestral_count/05_count/guide_clean.treefile")
SHANNON = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
               "B01_shannon_diversity.tsv")
LIB_DIR = Path("/data/carlos/EDF-TE/TE_analysis/annotated")
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
CLADE_ORDER = ["HO", "CR", "MI", "AP", "OL", "BL", "CH", "NU", "BS",
               "EN", "ZO", "KX", "MO", "GL", "MU", "BA", "AS"]

# TE classes for stacked bars (broad)
TE_CLASSES = ["LTR", "DNA", "LINE", "SINE", "RC", "Unknown"]
TE_CLASS_COLORS = {
    "LTR": "#e41a1c", "DNA": "#377eb8", "LINE": "#4daf4a",
    "SINE": "#f1c232", "RC": "#984ea3", "Unknown": "#bdbdbd",
}


def build_family_class_map():
    """Parse all combined libraries to map (species, family_name_upper) → TE class."""
    family_class = {}  # (species, family_upper) → broad_class
    libs = sorted(LIB_DIR.glob("*_EarlGrey/*_Curated_Library/*_combined_library.fasta"))
    for lib in libs:
        species = lib.parent.name.replace("_Curated_Library", "")
        with open(lib) as f:
            for line in f:
                if line.startswith(">"):
                    # Format: >family_name#CLASS/SUPERFAMILY or >family_name#CLASS
                    header = line[1:].strip().split()[0]
                    if "#" in header:
                        name, cls_sf = header.split("#", 1)
                        broad = cls_sf.split("/")[0]
                        # Normalize
                        name_upper = name.upper()
                        if broad not in TE_CLASSES:
                            broad = "Unknown"
                        family_class[(species, name_upper)] = broad
    return family_class


def classify_pa_families(family_class):
    """For each PA matrix family (cluster rep), assign TE class."""
    # Read PA matrix header (species columns) and first column (family names)
    family_classes = {}
    with open(PA_MATRIX) as f:
        header = f.readline().strip().split("\t")
        species_cols = header[1:]  # first col is "Family"
        for line in f:
            parts = line.strip().split("\t")
            fam_name = parts[0]  # e.g. "CH_Blyhe1|RND-1_FAMILY-195"
            if "|" in fam_name:
                sp, local = fam_name.split("|", 1)
                key = (sp, local.upper())
                cls = family_class.get(key, "Unknown")
            else:
                cls = "Unknown"
            family_classes[fam_name] = cls
    return family_classes, species_cols


def count_per_species(family_classes, species_cols):
    """Count FL TE families per species per TE class from PA matrix."""
    counts = defaultdict(lambda: Counter())  # species → Counter(class → count)

    with open(PA_MATRIX) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            fam_name = parts[0]
            cls = family_classes.get(fam_name, "Unknown")
            for i, val in enumerate(parts[1:]):
                if val == "1":
                    sp = species_cols[i]
                    counts[sp][cls] += 1
    return counts


def read_tree_order():
    """Read tree and return tip order (top to bottom) using ete3."""
    from ete3 import Tree
    t = Tree(str(TREE_FILE), format=1)
    # Ladderize for consistent ordering
    t.ladderize(direction=1)
    tips = [leaf.name for leaf in t.iter_leaves()]
    return tips, t


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Build family → class mapping from libraries
    print("[1/5] Building family → TE class mapping from 238 libraries...")
    family_class = build_family_class_map()
    print(f"       Mapped {len(family_class)} family entries")

    # Step 2: Classify PA matrix families
    print("[2/5] Classifying PA matrix families...")
    family_classes, species_cols = classify_pa_families(family_class)
    class_dist = Counter(family_classes.values())
    print(f"       {len(family_classes)} families classified: {dict(class_dist)}")

    # Step 3: Count per species
    print("[3/5] Counting families per species per class...")
    sp_counts = count_per_species(family_classes, species_cols)
    print(f"       {len(sp_counts)} species with counts")

    # Step 4: Read tree order
    print("[4/5] Reading tree...")
    tip_order, tree = read_tree_order()
    print(f"       {len(tip_order)} tips")

    # Step 5: Read supplementary data
    print("[5/5] Reading Shannon diversity and gains/losses...")
    df_shannon = pd.read_csv(SHANNON, sep="\t")
    shannon_map = dict(zip(df_shannon["species"], df_shannon["H_prime"]))

    df_gl = pd.read_csv(GAINS_LOSSES, sep="\t")
    gl_map = {r["node"]: r for _, r in df_gl.iterrows()}

    # ── Build data matrix ─────────────────────────────────────────────────
    # For each tip in tree order, build class counts
    species_in_tree = [sp for sp in tip_order if sp in sp_counts]
    n_sp = len(species_in_tree)
    print(f"\n{n_sp} species with both tree position and family counts")

    # Build arrays
    data = np.zeros((n_sp, len(TE_CLASSES)))
    shannons = np.zeros(n_sp)
    for i, sp in enumerate(species_in_tree):
        for j, cls in enumerate(TE_CLASSES):
            data[i, j] = sp_counts[sp].get(cls, 0)
        shannons[i] = shannon_map.get(sp, 0)

    # ── Create figure ─────────────────────────────────────────────────────
    # Panel C only (full species tree with stacked bars + Shannon)
    # since 238 species is too large for panels A+B+C on one page
    fig, axes = plt.subplots(1, 4, figsize=(18, 45),
                              gridspec_kw={"width_ratios": [0.3, 4, 4, 1],
                                           "wspace": 0.02})
    ax_clade, ax_tree, ax_bars, ax_shan = axes

    # ── Clade strip ───────────────────────────────────────────────────────
    for i, sp in enumerate(species_in_tree):
        clade = sp[:2]
        color = CLADE_COLORS.get(clade, "#aaa")
        ax_clade.barh(i, 1, height=0.8, color=color, edgecolor="none")
    ax_clade.set_ylim(-0.5, n_sp - 0.5)
    ax_clade.invert_yaxis()
    ax_clade.set_xlim(0, 1)
    ax_clade.axis("off")

    # ── Tree panel (simple cladogram) ─────────────────────────────────────
    # Draw a simplified cladogram using ete3 coordinates
    # We just show species names since the full tree is complex
    for i, sp in enumerate(species_in_tree):
        clade = sp[:2]
        ax_tree.text(0, i, sp.replace("_", " "), fontsize=2.5,
                    va="center", ha="left",
                    color=CLADE_COLORS.get(clade, "#333"))
    ax_tree.set_ylim(-0.5, n_sp - 0.5)
    ax_tree.invert_yaxis()
    ax_tree.axis("off")

    # ── Stacked bar chart ─────────────────────────────────────────────────
    left = np.zeros(n_sp)
    for j, cls in enumerate(TE_CLASSES):
        ax_bars.barh(range(n_sp), data[:, j], left=left, height=0.8,
                    color=TE_CLASS_COLORS.get(cls, "#aaa"), edgecolor="none",
                    label=cls)
        left += data[:, j]
    ax_bars.set_ylim(-0.5, n_sp - 0.5)
    ax_bars.invert_yaxis()
    ax_bars.set_xlabel("Number of FL TE families", fontsize=9)
    ax_bars.tick_params(axis="y", left=False, labelleft=False)
    ax_bars.spines["top"].set_visible(False)
    ax_bars.spines["right"].set_visible(False)
    ax_bars.spines["left"].set_visible(False)

    # Legend for TE classes
    handles = [mpatches.Patch(color=TE_CLASS_COLORS[c], label=c)
               for c in TE_CLASSES]
    ax_bars.legend(handles=handles, loc="lower right", fontsize=7,
                   title="TE class", title_fontsize=8, ncol=2,
                   framealpha=0.8)

    # ── Shannon diversity dots ────────────────────────────────────────────
    ax_shan.scatter(shannons, range(n_sp), s=6, c="steelblue",
                   edgecolor="none", alpha=0.7)
    ax_shan.set_ylim(-0.5, n_sp - 0.5)
    ax_shan.invert_yaxis()
    ax_shan.set_xlabel("Shannon\nH'", fontsize=8)
    ax_shan.tick_params(axis="y", left=False, labelleft=False)
    ax_shan.spines["top"].set_visible(False)
    ax_shan.spines["right"].set_visible(False)
    ax_shan.spines["left"].set_visible(False)
    ax_shan.set_xlim(0, max(shannons) * 1.1 if max(shannons) > 0 else 1)
    ax_shan.tick_params(axis="x", labelsize=6)

    fig.suptitle("FL TE family composition across fungal phylogeny",
                fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.995])

    out = OUTDIR / "te_family_barchart_tree.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    fig.savefig(OUTDIR / "te_family_barchart_tree.png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[DONE] Full tree panel: {out}")

    # ── Panel A: Collapsed clade summary ──────────────────────────────────
    fig2, (ax2_bars, ax2_shan) = plt.subplots(1, 2, figsize=(10, 6),
                                               gridspec_kw={"width_ratios": [3, 1],
                                                             "wspace": 0.08})
    # Aggregate by clade
    clade_data = {}
    clade_shannon = {}
    for c in CLADE_ORDER:
        idx = [i for i, sp in enumerate(species_in_tree) if sp[:2] == c]
        if not idx:
            continue
        clade_data[c] = data[idx].mean(axis=0)
        clade_shannon[c] = shannons[idx].mean()

    present_clades = [c for c in CLADE_ORDER if c in clade_data]
    n_c = len(present_clades)

    left2 = np.zeros(n_c)
    for j, cls in enumerate(TE_CLASSES):
        vals = [clade_data[c][j] for c in present_clades]
        ax2_bars.barh(range(n_c), vals, left=left2, height=0.7,
                     color=TE_CLASS_COLORS.get(cls, "#aaa"),
                     edgecolor="none", label=cls)
        left2 += vals

    ax2_bars.set_yticks(range(n_c))
    ax2_bars.set_yticklabels([f"{CLADE_NAMES.get(c,c)} ({c})" for c in present_clades],
                             fontsize=8)
    ax2_bars.invert_yaxis()
    ax2_bars.set_xlabel("Mean # FL TE families per species", fontsize=9)
    ax2_bars.spines["top"].set_visible(False)
    ax2_bars.spines["right"].set_visible(False)
    handles = [mpatches.Patch(color=TE_CLASS_COLORS[c], label=c) for c in TE_CLASSES]
    ax2_bars.legend(handles=handles, loc="lower right", fontsize=7,
                   title="TE class", ncol=2)

    # Shannon per clade
    shan_vals = [clade_shannon[c] for c in present_clades]
    ax2_shan.barh(range(n_c), shan_vals, height=0.6, color="steelblue",
                  alpha=0.7, edgecolor="none")
    ax2_shan.set_yticks(range(n_c))
    ax2_shan.set_yticklabels([])
    ax2_shan.invert_yaxis()
    ax2_shan.set_xlabel("Mean Shannon H'", fontsize=9)
    ax2_shan.spines["top"].set_visible(False)
    ax2_shan.spines["right"].set_visible(False)

    fig2.suptitle("A. FL TE family composition by clade",
                 fontsize=12, fontweight="bold")
    fig2.tight_layout()
    out2 = OUTDIR / "te_family_barchart_clade.pdf"
    fig2.savefig(out2, dpi=300, bbox_inches="tight")
    fig2.savefig(OUTDIR / "te_family_barchart_clade.png", dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"[DONE] Clade panel: {out2}")

    # ── Panel B: Net changes at key nodes ─────────────────────────────────
    # Summarize: per clade, sum gains and losses for terminal branches
    fig3, ax3 = plt.subplots(figsize=(8, 5))
    clade_gains = defaultdict(int)
    clade_losses = defaultdict(int)
    for _, row in df_gl.iterrows():
        node = row["node"]
        ntype = row["node_type"]
        if ntype == "tip":
            clade = node[:2]
            clade_gains[clade] += row["family_gain"]
            clade_losses[clade] += row["family_loss"]

    present_clades_b = [c for c in CLADE_ORDER if c in clade_gains]
    n_cb = len(present_clades_b)
    x = np.arange(n_cb)
    w = 0.35

    gains = [clade_gains[c] for c in present_clades_b]
    losses = [clade_losses[c] for c in present_clades_b]

    ax3.bar(x - w/2, gains, w, color="#2166ac", label="Gains", alpha=0.8)
    ax3.bar(x + w/2, losses, w, color="#b2182b", label="Losses", alpha=0.8)

    ax3.set_xticks(x)
    ax3.set_xticklabels([CLADE_NAMES.get(c,c)[:12] for c in present_clades_b],
                        fontsize=7, rotation=45, ha="right")
    ax3.set_ylabel("Total TE family changes\n(terminal branches)", fontsize=9)
    ax3.legend(fontsize=8)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.set_title("B. TE family gains and losses by clade",
                 fontsize=11, fontweight="bold", loc="left")

    fig3.tight_layout()
    out3 = OUTDIR / "te_family_gains_losses_clade.pdf"
    fig3.savefig(out3, dpi=300, bbox_inches="tight")
    fig3.savefig(OUTDIR / "te_family_gains_losses_clade.png", dpi=150, bbox_inches="tight")
    plt.close(fig3)
    print(f"[DONE] Gains/losses panel: {out3}")

    # ── Summary stats ─────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Per-species FL TE family counts by class:")
    total_per_class = data.sum(axis=0)
    for j, cls in enumerate(TE_CLASSES):
        print(f"  {cls:10s}  total={int(total_per_class[j]):>7,}  "
              f"mean/sp={data[:,j].mean():.1f}  max={int(data[:,j].max())}")
    print(f"  {'TOTAL':10s}  total={int(data.sum()):>7,}  "
          f"mean/sp={data.sum(axis=1).mean():.1f}")

    print(f"\nTop 10 species by total FL TE families:")
    totals = data.sum(axis=1)
    top_idx = np.argsort(totals)[::-1][:10]
    for i in top_idx:
        print(f"  {species_in_tree[i]:25s}  {int(totals[i]):>5} families")


if __name__ == "__main__":
    main()
