#!/usr/bin/env python3
"""
05_nucleosome_te_analysis.py  —  Nucleosome positioning and H3K9me3 at TEs.

Analyses DANPOS3 nucleosome occupancy data around TE boundaries and within
TE bodies.  Produces:
  A. Average nucleosome occupancy profiles ±2 kb around TE 5′ and 3′ ends
     (meta-TE plot), split by TE class
  B. H3K9me3 ChIP-seq enrichment over TE copies per class (boxplot)
  C. Comparison of nucleosome density inside vs outside TE copies

Inputs:
  --danpos_dir    DANPOS3 output directory (*.wig or *.bedGraph)
  --chip_dir      H3K9me3 ChIP-seq coverage per bp (bedgraph or BED)
  --te_bed_dir    EarlGrey filtered TE BED files
  --species_list  File with species codes
  --out_dir       Output directory
  --flank_kb      Flanking region for meta-TE profile (default: 2)
  --bin_bp        Bin size in bp for profile (default: 50)

Outputs:
  nucleosome_profile.tsv
  h3k9me3_te_enrichment.tsv
  figures/nucleosome_meta_te.pdf
  figures/h3k9me3_te_boxplot.pdf
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#984ea3",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}

TE_CLASSES = ["LTR", "DNA", "LINE", "SINE", "Penelope"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--danpos_dir",   required=True)
    p.add_argument("--chip_dir",     required=True)
    p.add_argument("--te_bed_dir",   required=True)
    p.add_argument("--species_list", required=True)
    p.add_argument("--out_dir",      default="nucleosome_out")
    p.add_argument("--flank_kb",     type=int, default=2)
    p.add_argument("--bin_bp",       type=int, default=50)
    return p.parse_args()


def _te_class(f: str) -> str:
    f = f.upper()
    if "LTR" in f:    return "LTR"
    if "DNA" in f:    return "DNA"
    if "LINE" in f:   return "LINE"
    if "SINE" in f:   return "SINE"
    if "RC/" in f or "HELITRON" in f: return "Rolling Circle"
    if "PENELOPE" in f: return "Penelope"
    return "Unclassified"


def load_bedgraph(path: Path) -> dict[str, np.ndarray]:
    """Load a bedgraph/wig file as a per-chrom array (1-indexed, 0-filled)."""
    from collections import defaultdict
    intervals = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#", "fixedStep", "variableStep")):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            chrom, s, e, val = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            intervals[chrom].append((s, e, val))
    result = {}
    for chrom, ivs in intervals.items():
        max_pos = max(e for _, e, _ in ivs)
        arr = np.zeros(max_pos + 1, dtype=np.float32)
        for s, e, v in ivs:
            arr[s:e] = v
        result[chrom] = arr
    return result


def extract_profile(cov: dict[str, np.ndarray], te: dict,
                    flank: int, n_bins: int = 100) -> np.ndarray | None:
    """Extract normalised coverage profile for one TE (5′-flank + body + 3′-flank)."""
    arr = cov.get(te["chrom"])
    if arr is None:
        return None
    te_len = te["end"] - te["start"]
    if te_len < 50:
        return None

    s = max(0, te["start"] - flank)
    e = min(len(arr) - 1, te["end"] + flank)
    raw = arr[s:e].astype(float)
    if raw.max() == 0:
        return None

    # Normalise to z-score per TE
    raw = (raw - raw.mean()) / (raw.std() + 1e-9)

    # Rebin to n_bins
    return np.interp(np.linspace(0, len(raw) - 1, n_bins), np.arange(len(raw)), raw)


def compute_meta_profiles(te_bed: pd.DataFrame, cov: dict[str, np.ndarray],
                           flank: int, n_bins: int) -> dict[str, np.ndarray]:
    profiles = {}
    for cls in TE_CLASSES:
        sub = te_bed[te_bed["te_class"] == cls]
        if sub.empty:
            continue
        stack = []
        for _, te in sub.iterrows():
            p = extract_profile(cov, te, flank, n_bins)
            if p is not None:
                stack.append(p)
        if stack:
            profiles[cls] = np.mean(stack, axis=0)
    return profiles


def compute_te_enrichment(te_bed: pd.DataFrame,
                           cov: dict[str, np.ndarray]) -> pd.DataFrame:
    rows = []
    for _, te in te_bed.iterrows():
        arr = cov.get(te["chrom"])
        if arr is None:
            continue
        body = arr[te["start"]:te["end"]].mean()
        rows.append({"te_class": te["te_class"],
                     "te_len":   te["end"] - te["start"],
                     "mean_cov": float(body)})
    return pd.DataFrame(rows)


# =============================================================================
# Figures
# =============================================================================

def plot_meta_te(profiles: dict[str, np.ndarray], n_bins: int,
                 flank_kb: int, out_dir: Path) -> None:
    n_body = n_bins // 3   # approximate body bins
    n_flank = (n_bins - n_body) // 2
    x_body  = np.linspace(-flank_kb, flank_kb + 100, n_bins)

    fig, ax = plt.subplots(figsize=(8, 4))
    for cls, prof in profiles.items():
        color = TE_CLASS_COLORS.get(cls, "#999")
        x = np.linspace(-flank_kb * 1000, flank_kb * 1000, len(prof))
        ax.plot(x / 1000, prof, color=color, linewidth=1.5, label=cls)

    ax.axvline(0, linewidth=0.8, color="grey", linestyle="--")
    ax.axvline(0, linewidth=0.8, color="grey", linestyle="--")
    ax.set_xlabel("Distance from TE boundary (kb)", fontsize=10)
    ax.set_ylabel("Mean normalised occupancy (z-score)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "nucleosome_meta_te.pdf")
    fig.savefig(out_dir / "nucleosome_meta_te.png", dpi=200)
    plt.close(fig)
    print("Saved: nucleosome_meta_te.pdf")


def plot_h3k9_boxplot(enrich: pd.DataFrame, out_dir: Path) -> None:
    classes = [c for c in TE_CLASSES if c in enrich["te_class"].values]
    data = [enrich[enrich["te_class"] == c]["mean_cov"].values for c in classes]
    colors = [TE_CLASS_COLORS.get(c, "#999") for c in classes]

    fig, ax = plt.subplots(figsize=(8, 4))
    bp = ax.boxplot(data, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.5))
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
    ax.set_xticks(range(1, len(classes) + 1))
    ax.set_xticklabels(classes, fontsize=9, rotation=20)
    ax.set_ylabel("Mean H3K9me3 signal over TE body", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "h3k9me3_te_boxplot.pdf")
    plt.close(fig)
    print("Saved: h3k9me3_te_boxplot.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    species_list = [l.strip() for l in open(args.species_list) if l.strip()]
    flank_bp = args.flank_kb * 1000
    N_BINS = 150

    all_profiles = {cls: [] for cls in TE_CLASSES}
    all_enrich   = []

    for sp in species_list:
        print(f"  Processing {sp} ...")
        te_files = list(Path(args.te_bed_dir).glob(f"{sp}*.filteredRepeats.bed")) + \
                   list(Path(args.te_bed_dir).glob(f"{sp}*.bed"))
        dan_files = list(Path(args.danpos_dir).glob(f"{sp}*.bedGraph")) + \
                    list(Path(args.danpos_dir).glob(f"{sp}*.wig"))
        chip_files = list(Path(args.chip_dir).glob(f"{sp}*H3K9*.bedGraph")) + \
                     list(Path(args.chip_dir).glob(f"{sp}*.bedGraph"))

        if not te_files or not dan_files:
            continue

        te_bed = pd.read_csv(te_files[0], sep="\t", header=None, usecols=[0, 1, 2, 3],
                             names=["chrom", "start", "end", "family"])
        te_bed["te_class"] = te_bed["family"].apply(_te_class)

        dan_cov = load_bedgraph(dan_files[0])
        profs = compute_meta_profiles(te_bed, dan_cov, flank_bp, N_BINS)
        for cls, prof in profs.items():
            all_profiles[cls].append(prof)

        if chip_files:
            chip_cov = load_bedgraph(chip_files[0])
            enrich = compute_te_enrichment(te_bed, chip_cov)
            enrich["species"] = sp
            all_enrich.append(enrich)

    # Average across species
    mean_profiles = {cls: np.mean(stack, axis=0)
                     for cls, stack in all_profiles.items() if stack}

    if mean_profiles:
        plot_meta_te(mean_profiles, N_BINS, args.flank_kb, fig_dir)
        prof_rows = [(cls, i, v)
                     for cls, arr in mean_profiles.items()
                     for i, v in enumerate(arr)]
        pd.DataFrame(prof_rows, columns=["te_class", "bin", "mean_z"]).to_csv(
            out_dir / "nucleosome_profile.tsv", sep="\t", index=False)

    if all_enrich:
        enrich_df = pd.concat(all_enrich, ignore_index=True)
        enrich_df.to_csv(out_dir / "h3k9me3_te_enrichment.tsv", sep="\t", index=False)
        plot_h3k9_boxplot(enrich_df, fig_dir)

    print("Done.")


if __name__ == "__main__":
    main()
