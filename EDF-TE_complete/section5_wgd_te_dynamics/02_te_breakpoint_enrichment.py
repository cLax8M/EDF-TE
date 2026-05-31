#!/usr/bin/env python3
"""
02_te_breakpoint_enrichment.py  —  TE enrichment at synteny breakpoints.

Tests whether TE copies are enriched at genomic positions corresponding to
synteny block boundaries, consistent with TE-mediated chromosomal rearrangement.

Analysis:
  1. Define breakpoint regions: ±5 kb from synteny block endpoints
  2. Compute TE density (bp TE / bp window) in breakpoint windows
  3. Compare to matched randomly sampled non-breakpoint intergenic windows
     (same size, same chromosome, 100 per breakpoint)
  4. Stratify by TE class: which classes drive breakpoint enrichment?
  5. Wilcoxon signed-rank test per species; meta-analysis across species

Inputs:
  --synteny_dir   Directory with pairwise synteny block BED files
                  (format: sp1, sp2 subdirectory with blocks.bed)
  --te_bed_dir    EarlGrey filtered TE BED files
  --species_list  Mucoromycotina species with synteny data
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory
  --window_kb     Half-window around breakpoint (default: 5)
  --n_random      Random comparison windows per breakpoint (default: 100)

Outputs:
  breakpoint_te_enrichment.tsv
  figures/breakpoint_te_density.pdf
  figures/breakpoint_te_class_composition.pdf
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#984ea3",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--synteny_dir",  required=True)
    p.add_argument("--te_bed_dir",   required=True)
    p.add_argument("--species_list", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="breakpoint_out")
    p.add_argument("--window_kb",    type=int, default=5)
    p.add_argument("--n_random",     type=int, default=100)
    p.add_argument("--seed",         type=int, default=42)
    return p.parse_args()


def _te_class(f: str) -> str:
    f = f.upper()
    if "LTR" in f:    return "LTR"
    if "DNA" in f:    return "DNA"
    if "LINE" in f:   return "LINE"
    if "SINE" in f:   return "SINE"
    if "RC/" in f:    return "Rolling Circle"
    if "PENELOPE" in f: return "Penelope"
    return "Unclassified"


def load_synteny_breakpoints(synteny_dir: str, sp: str,
                               window_bp: int) -> pd.DataFrame:
    """Load synteny blocks and extract breakpoint windows."""
    block_files = list(Path(synteny_dir).glob(f"{sp}*blocks*.bed")) + \
                  list(Path(synteny_dir).glob(f"**/{sp}*blocks*.bed")) + \
                  list(Path(synteny_dir).glob(f"{sp}*/*.bed"))
    if not block_files:
        return pd.DataFrame()

    blocks = pd.read_csv(block_files[0], sep="\t", header=None, usecols=[0, 1, 2],
                         names=["chrom", "start", "end"])
    # Breakpoints at block starts and ends
    bps = []
    for _, row in blocks.iterrows():
        for pos in [row["start"], row["end"]]:
            bps.append({
                "chrom": row["chrom"],
                "bp_pos": pos,
                "win_start": max(0, pos - window_bp),
                "win_end":   pos + window_bp,
            })
    return pd.DataFrame(bps).drop_duplicates()


def te_density_in_window(te_bed: pd.DataFrame, chrom: str,
                          start: int, end: int) -> float:
    win_len = end - start
    if win_len <= 0:
        return 0.0
    sub = te_bed[(te_bed["chrom"] == chrom) &
                 (te_bed["te_end"] > start) & (te_bed["te_start"] < end)]
    if sub.empty:
        return 0.0
    ov = (sub["te_end"].clip(upper=end) - sub["te_start"].clip(lower=start)).sum()
    return ov / win_len


def te_class_bp_in_window(te_bed: pd.DataFrame, chrom: str,
                            start: int, end: int) -> dict[str, float]:
    win_len = end - start
    if win_len <= 0:
        return {}
    sub = te_bed[(te_bed["chrom"] == chrom) &
                 (te_bed["te_end"] > start) & (te_bed["te_start"] < end)]
    if sub.empty:
        return {}
    result = {}
    for cls, grp in sub.groupby("te_class"):
        ov = (grp["te_end"].clip(upper=end) - grp["te_start"].clip(lower=start)).sum()
        result[cls] = ov / win_len
    return result


def random_windows_for_chrom(chrom: str, chrom_len: int,
                              win_len: int, n: int,
                              rng: np.random.Generator) -> list[tuple[int, int]]:
    out = []
    for _ in range(n * 10):
        s = int(rng.integers(0, max(1, chrom_len - win_len)))
        out.append((s, s + win_len))
        if len(out) >= n:
            break
    return out[:n]


def estimate_chrom_lengths(te_bed: pd.DataFrame) -> dict[str, int]:
    return te_bed.groupby("chrom")["te_end"].max().to_dict()


# =============================================================================
# Figures
# =============================================================================

def plot_density_comparison(df: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))

    sp_list = sorted(df["species"].unique())
    bp_vals   = [df[df["species"] == sp]["bp_te_density"].values    for sp in sp_list]
    rand_vals = [df[df["species"] == sp]["rand_te_density"].values  for sp in sp_list]

    combined_bp   = np.concatenate(bp_vals)
    combined_rand = np.concatenate(rand_vals)

    ax.boxplot([combined_rand, combined_bp], patch_artist=True,
               medianprops=dict(color="black", linewidth=1.5),
               labels=["Random windows", "Breakpoint windows"])
    ax.get_children()[0].set_facecolor("#a6cee3")
    ax.get_children()[2].set_facecolor("#e31a1c")

    if len(combined_rand) >= 5:
        stat, pval = wilcoxon(combined_bp[:len(combined_rand)],
                              combined_rand[:len(combined_bp)])
        ax.text(0.5, 0.96, f"Wilcoxon p={pval:.3g}",
                transform=ax.transAxes, ha="center", va="top", fontsize=9)

    ax.set_ylabel("TE density (bp TE / bp window)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "breakpoint_te_density.pdf")
    fig.savefig(out_dir / "breakpoint_te_density.png", dpi=200)
    plt.close(fig)
    print("Saved: breakpoint_te_density.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))
    species_list = [l.strip() for l in open(args.species_list) if l.strip()]
    rng = np.random.default_rng(args.seed)
    window_bp = args.window_kb * 1000

    all_rows = []
    for sp in species_list:
        te_files = list(Path(args.te_bed_dir).glob(f"{sp}*.filteredRepeats.bed")) + \
                   list(Path(args.te_bed_dir).glob(f"{sp}*.bed"))
        if not te_files:
            continue
        te_bed = pd.read_csv(te_files[0], sep="\t", header=None, usecols=[0, 1, 2, 3],
                             names=["chrom", "te_start", "te_end", "family"])
        te_bed["te_class"] = te_bed["family"].apply(_te_class)

        bps = load_synteny_breakpoints(args.synteny_dir, sp, window_bp)
        if bps.empty:
            continue

        chrom_lens = estimate_chrom_lengths(te_bed)
        print(f"  {sp}: {len(bps)} breakpoint windows ...")

        for _, bp in bps.iterrows():
            bp_dens = te_density_in_window(te_bed, bp["chrom"],
                                           bp["win_start"], bp["win_end"])
            win_len = bp["win_end"] - bp["win_start"]
            c_len   = chrom_lens.get(bp["chrom"], bp["win_end"] + 1)
            rw      = random_windows_for_chrom(bp["chrom"], c_len,
                                               win_len, args.n_random, rng)
            rand_dens = np.mean([te_density_in_window(te_bed, r[0] if isinstance(r, tuple) else bp["chrom"],
                                                       r[0] if isinstance(r, tuple) else r[0],
                                                       r[1] if isinstance(r, tuple) else r[1])
                                 for r in rw]) if rw else 0.0

            # Correct the call above — rw is list of (start, end) on same chrom
            rand_dens = np.mean([te_density_in_window(te_bed, bp["chrom"], s, e)
                                  for s, e in rw]) if rw else 0.0

            all_rows.append({
                "species":        sp,
                "clade":          clade_map.get(sp, "?"),
                "chrom":          bp["chrom"],
                "bp_pos":         bp["bp_pos"],
                "bp_te_density":  bp_dens,
                "rand_te_density": rand_dens,
            })

    df = pd.DataFrame(all_rows)
    if df.empty:
        print("No breakpoint data generated.", file=sys.stderr)
        sys.exit(1)

    df.to_csv(out_dir / "breakpoint_te_enrichment.tsv", sep="\t", index=False)
    plot_density_comparison(df, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
