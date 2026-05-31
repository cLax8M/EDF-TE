#!/usr/bin/env python3
"""03_gremlines_tree_figure.py — Enhanced Grem-LINE1 tree with centromere labels.

For each copy on the tree, the tip label shows:
  Cen<N>_<side>_<zone>_<status>  (e.g. Cen1_L_peri_AUT)
where:
  Cen<N>   = centromere identity (Cen1–Cen9, or "NoCen" if no centromere on that scaffold)
  side     = L (left / 5') or R (right / 3') of centromere center
  zone     = core (≤10kb), peri (10-100kb), no_centromere
  status   = AUT (autonomous, clusters 1-3) or REM (remnant)

Also:
  - Colors tips by centromere identity
  - Shapes by AUT/REM
  - Infers putative parent–offspring pairs (sister taxa with one young + one older)
  - Tests: do new insertions go inside (towards centromere core) or outside?

Usage:
  python3 03_gremlines_tree_figure.py
"""
from __future__ import annotations
import re, sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

# ── paths ──────────────────────────────────────────────────────────────────
OUT       = Path("/legserv/Temp/Carlos/EDF-TE/results/Gremlines/MU402")
CEN_BED   = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/MU402_centromere_windows_clipped.bed")
COPY_ZONE = OUT / "copy_zone.tsv"
TREE_FILE = OUT / "gremlin_iqtree.treefile"
FAI_FILE  = Path("/legserv/Temp/Carlos/EDF-TE/data/Extra/mu402GenomeJGI.fa.fai")

# ── centromere scaffold → number mapping ───────────────────────────────────
# MU402 centromeres: scaffolds 1, 2, 3, 4, 5, 6, 8, 9, 11 (9 centromeres)
CEN_SCAFF_ORDER = {
    "scaffold_1": 1, "scaffold_2": 2, "scaffold_3": 3,
    "scaffold_4": 4, "scaffold_5": 5, "scaffold_6": 6,
    "scaffold_8": 7, "scaffold_9": 8, "scaffold_11": 9,
}

def fasta_to_tree_name(fasta_name: str) -> str:
    """Convert copy_zone fasta_name e.g. GremLINE1_1.01_AUT(+) → tree tip name GremLINE1_1.01_AUT___"""
    return fasta_name.replace("(+)", "___").replace("(-)", "_-_")

def tree_to_fasta_name(tree_name: str) -> str:
    """Convert tree tip name back to fasta_name."""
    if tree_name.endswith("___"):
        return tree_name[:-3] + "(+)"
    elif tree_name.endswith("_-_"):
        return tree_name[:-3] + "(-)"
    return tree_name

CEN_COLORS = {
    1: "#e41a1c", 2: "#377eb8", 3: "#4daf4a",
    4: "#984ea3", 5: "#ff7f00", 6: "#a65628",
    7: "#f781bf", 8: "#999999", 9: "#66c2a5",
    0: "#cccccc",  # no centromere
}


def load_centromeres():
    """Return dict scaffold → (start, end, center)."""
    cen = {}
    with open(CEN_BED) as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            cen[chrom] = (start, end, (start + end) // 2)
    return cen


def assign_centromere_label(row, cen_map):
    """For a copy, determine CenN and side (L/R)."""
    chrom = row["chrom"]
    mid = int(row["mid"])
    if chrom not in cen_map:
        return 0, "NoCen", "-"
    cen_start, cen_end, cen_center = cen_map[chrom]
    cen_num = CEN_SCAFF_ORDER.get(chrom, 0)
    side = "L" if mid < cen_center else "R"
    return cen_num, f"Cen{cen_num}", side


def make_tree_label(row, cen_num, cen_name, side):
    """Build informative tree tip label."""
    zone = row["zone"]
    status = row["status"]
    original = row["name"].replace("_AUT", "").replace("_REM", "")
    # Extract cluster number from original name
    m = re.match(r'GremLINE1_(\d+)\.(\d+)', original)
    if m:
        cluster, copy_num = m.group(1), m.group(2)
    else:
        cluster, copy_num = "?", "?"

    return f"{cen_name}_{side}_{zone}_{status}|Cl{cluster}.{copy_num}"


def parse_newick_tips(nwk_str):
    """Extract tip names from a Newick string."""
    tips = re.findall(r'[(,]([^:,()]+):', nwk_str)
    return tips


def find_sister_pairs(nwk_str):
    """Find pairs of taxa that are sister to each other (leaf,leaf) in the tree.
    Returns list of (taxon_a, taxon_b) tuples — using tree tip names."""
    pairs = re.findall(r'\(([^:,()]+):[^,]+,([^:,()]+):[^)]+\)', nwk_str)
    return pairs


def main():
    cen_map = load_centromeres()

    # Load copy data
    df = pd.read_csv(COPY_ZONE, sep="\t")
    print(f"[INFO] {len(df)} copies loaded from copy_zone.tsv")

    # Build tree tip name column
    df["tree_tip"] = df["fasta_name"].apply(fasta_to_tree_name)

    # Assign centromere labels
    cen_nums, cen_names, sides = [], [], []
    for _, row in df.iterrows():
        cn, cname, side = assign_centromere_label(row, cen_map)
        cen_nums.append(cn)
        cen_names.append(cname)
        sides.append(side)
    df["cen_num"] = cen_nums
    df["cen_name"] = cen_names
    df["side"] = sides

    # Build new tree labels (mapping tree_tip → new label)
    label_map = {}          # tree_tip → readable label
    label_map_fasta = {}    # fasta_name → readable label
    for _, row in df.iterrows():
        new_label = make_tree_label(row, row["cen_num"], row["cen_name"], row["side"])
        label_map[row["tree_tip"]] = new_label
        label_map_fasta[row["fasta_name"]] = new_label

    # Save enriched copy table
    df["tree_label"] = [label_map_fasta.get(fn, fn) for fn in df["fasta_name"]]
    df.to_csv(OUT / "copy_zone_enriched.tsv", sep="\t", index=False)
    print(f"[INFO] Enriched table saved to copy_zone_enriched.tsv")

    # ── Build lookup maps keyed by tree_tip name ─────────────────────────
    tree_str = TREE_FILE.read_text().strip()
    tree_tips = parse_newick_tips(tree_str)
    print(f"[INFO] {len(tree_tips)} tips in tree")

    sisters = find_sister_pairs(tree_str)
    print(f"[INFO] {len(sisters)} sister pairs found in tree")

    # Maps keyed by tree_tip name
    k80_map = dict(zip(df["tree_tip"], df["k80"]))
    dist_map = dict(zip(df["tree_tip"], df["dist_centro"]))
    zone_map = dict(zip(df["tree_tip"], df["zone"]))
    cen_map_copy = dict(zip(df["tree_tip"], df["cen_num"]))

    daughter_closer = 0
    daughter_farther = 0
    same_centromere = 0
    diff_centromere = 0
    pair_details = []

    for a, b in sisters:
        # a, b are tree tip names (e.g. GremLINE1_1.01_AUT___)
        if a not in k80_map or b not in k80_map:
            continue
        k_a, k_b = k80_map[a], k80_map[b]
        d_a, d_b = dist_map.get(a, np.nan), dist_map.get(b, np.nan)
        if pd.isna(k_a) or pd.isna(k_b):
            continue

        # Determine parent (older = higher K80) vs daughter (younger)
        if k_a > k_b:
            parent, daughter = a, b
            k_parent, k_daughter = k_a, k_b
            d_parent, d_daughter = d_a, d_b
        else:
            parent, daughter = b, a
            k_parent, k_daughter = k_b, k_a
            d_parent, d_daughter = d_b, d_a

        # Same centromere?
        c_parent = cen_map_copy.get(parent, -1)
        c_daughter = cen_map_copy.get(daughter, -1)
        if c_parent == c_daughter and c_parent > 0:
            same_centromere += 1
        elif c_parent > 0 and c_daughter > 0:
            diff_centromere += 1

        if not pd.isna(d_parent) and not pd.isna(d_daughter):
            if d_daughter < d_parent:
                daughter_closer += 1
            else:
                daughter_farther += 1

        pair_details.append({
            "parent": label_map.get(parent, parent),
            "daughter": label_map.get(daughter, daughter),
            "k80_parent": k_parent,
            "k80_daughter": k_daughter,
            "dist_parent": d_parent,
            "dist_daughter": d_daughter,
            "same_centromere": c_parent == c_daughter and c_parent > 0,
        })

    pairs_df = pd.DataFrame(pair_details)
    if len(pairs_df):
        pairs_df.to_csv(OUT / "sister_pairs_directionality.tsv", sep="\t", index=False)

    # Binomial test: are daughters more likely to insert closer or farther?
    total_directional = daughter_closer + daughter_farther
    if total_directional > 0:
        binom_p = stats.binomtest(daughter_closer, total_directional, 0.5).pvalue
    else:
        binom_p = float("nan")

    # ── Write directionality results ─────────────────────────────────────
    with open(OUT / "directionality_results.txt", "w") as f:
        f.write("Grem-LINE1 insertion directionality analysis (MU402 JGI)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Sister pairs analysed: {len(sisters)}\n")
        f.write(f"  with K80 data for both: {len(pair_details)}\n\n")
        f.write("[Q1] Do new insertions go TOWARDS or AWAY from centromere core?\n")
        f.write(f"  Daughter closer to centro: {daughter_closer}\n")
        f.write(f"  Daughter farther from centro: {daughter_farther}\n")
        f.write(f"  Binomial test p = {binom_p:.4g}\n")
        if daughter_farther > daughter_closer:
            f.write("  → New insertions tend to be FARTHER from centromere (expansion outward)\n")
        elif daughter_closer > daughter_farther:
            f.write("  → New insertions tend to be CLOSER to centromere (insertion inward)\n")
        else:
            f.write("  → No directional bias\n")
        f.write(f"\n[Q2] Do sister copies stay in the SAME centromere?\n")
        f.write(f"  Same centromere: {same_centromere}\n")
        f.write(f"  Different centromere: {diff_centromere}\n")
        if same_centromere + diff_centromere > 0:
            pct = 100 * same_centromere / (same_centromere + diff_centromere)
            f.write(f"  → {pct:.0f}% of sister pairs in same centromere\n")
        f.write("\n[Q3] Per-centromere copy counts:\n")
        for cn in sorted(df["cen_num"].unique()):
            sub = df[df["cen_num"] == cn]
            aut_n = (sub["status"] == "AUT").sum()
            rem_n = (sub["status"] == "REM").sum()
            name = f"Cen{cn}" if cn > 0 else "NoCen"
            f.write(f"  {name:8s}: {len(sub):3d} copies ({aut_n} AUT, {rem_n} REM)\n")

    print("[INFO] Directionality results saved")

    # ── Figure 1: Tree with centromere-colored tips ──────────────────────
    # Since we don't have ete3/toytree necessarily, make a simpler plot:
    # Horizontal layout showing copies grouped by centromere, colored by K80

    fig, axes = plt.subplots(1, 2, figsize=(14, 8), gridspec_kw={"width_ratios": [2, 1]})

    # Panel A: Copies organized by centromere, colored by K80, shape by AUT/REM
    ax = axes[0]
    df_sorted = df.sort_values(["cen_num", "dist_centro"])
    y_positions = np.arange(len(df_sorted))
    vmin, vmax = df_sorted["k80"].quantile(0.05), df_sorted["k80"].quantile(0.95)

    for i, (_, row) in enumerate(df_sorted.iterrows()):
        color = CEN_COLORS.get(row["cen_num"], "#cccccc")
        marker = "o" if row["status"] == "AUT" else "s"
        ax.scatter(row["k80"], i, c=row["k80"], cmap="viridis_r",
                   vmin=vmin, vmax=vmax, s=40, marker=marker,
                   edgecolor=color, linewidth=1.5, zorder=3)
        ax.text(row["k80"] + 0.0005, i, row["tree_label"],
                fontsize=4.5, va="center", ha="left",
                color=color, fontweight="bold" if row["status"] == "AUT" else "normal")

    ax.set_xlabel("K80 distance to consensus (age →)")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_title("A. Grem-LINE1 copies by centromere (MU402 JGI)")
    ax.axvline(df_sorted[df_sorted["status"]=="AUT"]["k80"].median(),
               ls="--", color="green", alpha=0.5, label="AUT median")
    ax.axvline(df_sorted[df_sorted["status"]=="REM"]["k80"].median(),
               ls="--", color="red", alpha=0.5, label="REM median")
    ax.legend(fontsize=7)

    # Panel B: Boxplot K80 per centromere
    ax2 = axes[1]
    cen_groups = []
    cen_labels = []
    for cn in sorted(df["cen_num"].unique()):
        sub = df[df["cen_num"] == cn]["k80"].dropna()
        if len(sub) > 0:
            cen_groups.append(sub.values)
            cen_labels.append(f"Cen{cn}" if cn > 0 else "NoCen")
    bp = ax2.boxplot(cen_groups, labels=cen_labels, patch_artist=True,
                     showfliers=False, vert=True)
    for i, (patch, cn) in enumerate(zip(bp["boxes"], sorted(df["cen_num"].unique()))):
        patch.set_facecolor(CEN_COLORS.get(cn, "#cccccc"))
        patch.set_alpha(0.6)
        data = cen_groups[i]
        ax2.scatter(np.random.normal(i+1, 0.08, len(data)), data,
                    s=12, alpha=0.5, c="k", zorder=3)
    ax2.set_ylabel("K80 distance to consensus")
    ax2.set_title("B. K80 per centromere")
    ax2.tick_params(axis="x", rotation=45)

    fig.tight_layout()
    fig.savefig(OUT / "gremlines_centromere_tree.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 2: Contig map with centromere labels + directionality arrows ──
    fai = pd.read_csv(FAI_FILE, sep="\t", header=None, names=["chrom", "length"],
                       usecols=[0, 1])
    # Only scaffolds with Grem-LINE1 copies
    chroms_with = df["chrom"].unique()
    fai_sub = fai[fai["chrom"].isin(chroms_with)].sort_values("length", ascending=False)

    cen_info = load_centromeres()
    n_chroms = len(fai_sub)
    fig, axes = plt.subplots(n_chroms, 1, figsize=(13, max(4, 0.6 * n_chroms)),
                              sharex=False)
    if n_chroms == 1:
        axes = [axes]

    for ax, (_, frow) in zip(axes, fai_sub.iterrows()):
        ch = frow["chrom"]
        L = frow["length"]

        # Scaffold bar
        ax.barh(0, L / 1e6, height=0.4, color="#e8e8e8", edgecolor="#999", linewidth=0.5)

        # Centromere region
        if ch in cen_info:
            cs, ce, cc = cen_info[ch]
            ax.barh(0, (ce - cs) / 1e6, left=cs / 1e6, height=0.5,
                    color="#ffcccc", edgecolor="#cc0000", linewidth=0.8, alpha=0.6)
            ax.axvline(cc / 1e6, color="#cc0000", ls=":", linewidth=0.8, alpha=0.5)

        # Plot copies
        sub = df[df["chrom"] == ch]
        for _, crow in sub.iterrows():
            color = CEN_COLORS.get(crow["cen_num"], "#cccccc")
            marker = "v" if crow["status"] == "AUT" else "^"
            y_off = 0.35 if crow["strand"] == "+" else -0.35
            ax.plot(crow["mid"] / 1e6, y_off, marker=marker, color=color,
                    markersize=6, markeredgecolor="black", markeredgewidth=0.4)
            # Label
            short = crow["tree_label"].split("|")[1] if "|" in str(crow.get("tree_label","")) else crow["name"][-8:]
            ax.text(crow["mid"] / 1e6, y_off + (0.15 if y_off > 0 else -0.15),
                    short, fontsize=3.5, ha="center",
                    va="bottom" if y_off > 0 else "top", rotation=45)

        cn = CEN_SCAFF_ORDER.get(ch, 0)
        label = f"{ch} (Cen{cn})" if cn > 0 else ch
        ax.set_ylabel(label, fontsize=7, rotation=0, ha="right", va="center")
        ax.set_xlim(0, L / 1e6)
        ax.set_ylim(-0.8, 0.8)
        ax.set_yticks([])
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)

    axes[-1].set_xlabel("Genomic position (Mb)")
    # Legend
    legend_elements = [
        Line2D([0], [0], marker="v", color="w", markerfacecolor="gray",
               markeredgecolor="k", markersize=8, label="AUT (+)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor="gray",
               markeredgecolor="k", markersize=8, label="REM (−)"),
    ]
    for cn in sorted(CEN_SCAFF_ORDER.values()):
        legend_elements.append(
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=CEN_COLORS[cn], markeredgecolor="k",
                   markersize=8, label=f"Cen{cn}")
        )
    fig.legend(handles=legend_elements, loc="upper right", fontsize=6,
               ncol=2, framealpha=0.9)
    fig.suptitle("MU402 Grem-LINE1 — Genomic map with centromere identity", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "gremlines_contig_map_labeled.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 3: Directionality arrow plot ─────────────────────────────
    if len(pairs_df) > 0:
        fig, ax = plt.subplots(figsize=(7, 5))
        valid = pairs_df.dropna(subset=["dist_parent", "dist_daughter"])
        ax.scatter(valid["dist_parent"] / 1e3, valid["dist_daughter"] / 1e3,
                   c=valid["same_centromere"].map({True: "#2ca02c", False: "#d62728"}),
                   s=30, alpha=0.7, edgecolor="k", linewidth=0.3)
        maxval = max(valid["dist_parent"].max(), valid["dist_daughter"].max()) / 1e3
        ax.plot([0, maxval], [0, maxval], "k:", alpha=0.4)
        ax.set_xlabel("Parent distance to centromere center (kb)")
        ax.set_ylabel("Daughter distance to centromere center (kb)")
        ax.set_title(f"Insertion directionality\n"
                     f"Daughter closer: {daughter_closer}, farther: {daughter_farther} "
                     f"(binom p={binom_p:.3g})")
        legend_elements = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#2ca02c",
                   markersize=8, label="Same centromere"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#d62728",
                   markersize=8, label="Different centromere"),
        ]
        ax.legend(handles=legend_elements, fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / "gremlines_directionality.pdf", dpi=300, bbox_inches="tight")
        plt.close(fig)

    print("[DONE] All enhanced Grem-LINE1 figures generated")


if __name__ == "__main__":
    main()
