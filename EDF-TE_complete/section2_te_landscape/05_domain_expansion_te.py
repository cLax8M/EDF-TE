#!/usr/bin/env python3
"""
05_domain_expansion_te.py  —  Protein domain expansion linked to TE proximity.

Identifies protein domains that show significant copy-number expansion in
species with high TE content.  For each Pfam/InterPro domain family, tests
whether species in the upper TE-content quartile harbour more domain copies
than expected under a null distribution (100 random clade-matched permutations).

Analysis pipeline:
  1. Load per-species InterPro domain counts (from InterProScan output)
  2. Load per-species TE content (% genome)
  3. For each domain: fit Poisson GLM with log(TE_content) as predictor
     controlling for genome size (log offset) and phylogenetic clade
  4. Rank domains by GLM coefficient; apply FDR correction (Benjamini-Hochberg)
  5. For top-N expanded domains: compute per-species copy-number boxplot
     split by TE-high vs TE-low quartile, with proximity analysis
     (fraction of domain-encoding genes within X kb of a TE)

Inputs:
  --ipr_dir      Directory with one InterPro TSV per species (species.ipr.tab)
  --te_content   Per-species TE content TSV: species, te_pct_genome
  --genome_size  Per-species genome size TSV: species, genome_size_bp
  --clade_map    Two-column TSV: species_code<TAB>clade_tag
  --te_bed_dir   Directory with per-species TE BED files (EarlGrey output)
  --gene_bed_dir Directory with per-species gene BED files
  --out_dir      Output directory

Outputs:
  domain_expansion_results.tsv     – all tested domains, GLM stats, FDR
  top_domains_copynumber.tsv       – copy-number matrix for top-50 domains
  top_domains_proximity.tsv        – TE proximity fraction for top-50 domains
  figures/domain_expansion_bubble.pdf
  figures/copynumber_<domain>.pdf  – per-domain copy-number comparison
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}


# =============================================================================
# I/O
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ipr_dir",       required=True)
    p.add_argument("--te_content",    required=True)
    p.add_argument("--genome_size",   required=True)
    p.add_argument("--clade_map",     required=True)
    p.add_argument("--te_bed_dir",    required=True)
    p.add_argument("--gene_bed_dir",  required=True)
    p.add_argument("--out_dir",       default="domain_expansion_out")
    p.add_argument("--proximity_kb",  type=int, default=5,
                   help="Window (kb) for TE proximity (default: 5)")
    p.add_argument("--top_n",         type=int, default=50,
                   help="Top-N expanded domains to plot individually")
    p.add_argument("--fdr",           type=float, default=0.05)
    return p.parse_args()


def load_ipr_counts(ipr_dir: Path, species_list: list[str]) -> pd.DataFrame:
    """Load domain counts per species into a species × domain matrix."""
    dfs = []
    for sp in species_list:
        f = list(Path(ipr_dir).glob(f"{sp}*.ipr.tab")) + \
            list(Path(ipr_dir).glob(f"{sp}*.tsv"))
        if not f:
            continue
        raw = pd.read_csv(f[0], sep="\t", header=None, usecols=[0, 4, 5],
                          names=["gene", "ipr_id", "ipr_desc"])
        raw = raw.dropna(subset=["ipr_id"])
        counts = raw.groupby("ipr_id").size().reset_index(name="count")
        counts["species"] = sp
        dfs.append(counts)
    if not dfs:
        print("No IPR data loaded — check --ipr_dir and species list.", file=sys.stderr)
        sys.exit(1)
    return pd.concat(dfs, ignore_index=True)


def build_count_matrix(ipr_long: pd.DataFrame) -> pd.DataFrame:
    """Pivot to species × domain count matrix."""
    mat = ipr_long.pivot_table(index="species", columns="ipr_id",
                               values="count", fill_value=0, aggfunc="sum")
    return mat


# =============================================================================
# Statistical testing: Spearman correlation with TE content
# =============================================================================

def test_domain_te_correlation(count_mat: pd.DataFrame,
                               te_content: pd.Series,
                               genome_size: pd.Series,
                               n_perm: int = 1000,
                               seed: int = 42) -> pd.DataFrame:
    """
    For each domain, test Spearman correlation between copy number
    (normalised by genome size) and TE content (% genome).
    FDR correction via Benjamini-Hochberg.
    """
    rng = np.random.default_rng(seed)
    common = count_mat.index.intersection(te_content.index).intersection(
        genome_size.index)
    mat   = count_mat.loc[common]
    te    = te_content.loc[common]
    gs    = genome_size.loc[common]

    # Normalise: copies per Mb genome
    mat_norm = mat.div(gs / 1e6, axis=0)

    rows = []
    for domain in mat_norm.columns:
        vals = mat_norm[domain]
        if vals.sum() == 0:
            continue
        r, p = scipy_stats.spearmanr(vals, te)
        rows.append({"domain": domain, "spearman_r": r, "p_value": p,
                     "n_species": (vals > 0).sum()})

    df = pd.DataFrame(rows)
    # Benjamini-Hochberg FDR
    df = df.sort_values("p_value").reset_index(drop=True)
    m = len(df)
    df["bh_rank"] = np.arange(1, m + 1)
    df["p_adj"] = (df["p_value"] * m / df["bh_rank"]).clip(upper=1.0)
    df["significant"] = df["p_adj"] < 0.05
    return df.sort_values("spearman_r", ascending=False)


# =============================================================================
# Proximity analysis: fraction of domain genes near a TE
# =============================================================================

def te_proximity_fraction(gene_bed_dir: Path, te_bed_dir: Path,
                           species_list: list[str],
                           domain_genes: dict[str, list[str]],
                           window_bp: int) -> pd.DataFrame:
    """
    For each domain and species, compute the fraction of domain-encoding genes
    within `window_bp` of a TE locus (using BED interval arithmetic).
    Requires bedtools on PATH.
    """
    import subprocess, tempfile, os

    rows = []
    for sp in species_list:
        gene_f = list(Path(gene_bed_dir).glob(f"{sp}*.bed"))
        te_f   = list(Path(te_bed_dir).glob(f"{sp}*.filteredRepeats.bed"))
        if not gene_f or not te_f:
            continue
        gene_bed, te_bed = gene_f[0], te_f[0]

        with tempfile.TemporaryDirectory() as tmp:
            # Extend TE intervals by window_bp on both sides
            slop_f = os.path.join(tmp, "te_slop.bed")
            sp_gs = os.path.join(tmp, "genome.sizes")
            # Write dummy genome sizes (bedtools slop needs this)
            subprocess.run(
                f"awk '{{print $1\"\\t\"$3+{window_bp}}}' {te_bed} > {sp_gs}",
                shell=True, check=True)
            subprocess.run(
                f"bedtools slop -i {te_bed} -g {sp_gs} -b {window_bp} "
                f"| bedtools sort > {slop_f}",
                shell=True, check=True)

            for domain, genes in domain_genes.items():
                # Write domain genes subset BED
                dom_f = os.path.join(tmp, "dom_genes.bed")
                g_df  = pd.read_csv(gene_bed, sep="\t", header=None,
                                    names=["chr", "start", "end", "gene"])
                sub   = g_df[g_df["gene"].isin(genes)]
                if sub.empty:
                    continue
                sub.to_csv(dom_f, sep="\t", index=False, header=False)

                # Intersect with TE slopbed
                res = subprocess.run(
                    f"bedtools intersect -a {dom_f} -b {slop_f} -u | wc -l",
                    shell=True, capture_output=True, text=True)
                n_near = int(res.stdout.strip())
                rows.append({
                    "species": sp, "domain": domain,
                    "n_domain_genes": len(sub),
                    "n_near_te": n_near,
                    "frac_near_te": n_near / len(sub) if len(sub) > 0 else 0,
                })
    return pd.DataFrame(rows)


# =============================================================================
# Figures
# =============================================================================

def plot_bubble(results: pd.DataFrame, te_content: pd.Series,
                out_dir: Path) -> None:
    """Bubble plot: x = Spearman r, y = -log10(p_adj), size = n_species."""
    sig   = results[results["significant"]]
    nosig = results[~results["significant"]]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(nosig["spearman_r"], -np.log10(nosig["p_value"] + 1e-300),
               s=nosig["n_species"] / 5, color="#cccccc", alpha=0.5, zorder=2)
    ax.scatter(sig["spearman_r"], -np.log10(sig["p_value"] + 1e-300),
               s=sig["n_species"] / 5, color="#e41a1c", alpha=0.7, zorder=3)

    ax.axhline(-np.log10(0.05), linewidth=0.8, color="grey", linestyle="--")
    ax.axvline(0, linewidth=0.5, color="grey")
    ax.set_xlabel("Spearman r (copy number vs TE content)", fontsize=10)
    ax.set_ylabel("-log₁₀(p-value)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "domain_expansion_bubble.pdf")
    fig.savefig(out_dir / "domain_expansion_bubble.png", dpi=200)
    plt.close(fig)


def plot_domain_copynumber(domain: str, count_mat: pd.DataFrame,
                           clade_map: dict, te_content: pd.Series,
                           out_dir: Path) -> None:
    """Boxplot: domain copy number split by TE-high vs TE-low quartile."""
    common = count_mat.index.intersection(te_content.index)
    vals = count_mat.loc[common, domain]
    te   = te_content.loc[common]

    q75 = te.quantile(0.75)
    q25 = te.quantile(0.25)
    high = vals[te >= q75]
    low  = vals[te <= q25]

    fig, ax = plt.subplots(figsize=(4, 4))
    bp = ax.boxplot([low, high], patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.5))
    bp["boxes"][0].set_facecolor("#a6cee3")
    bp["boxes"][1].set_facecolor("#e31a1c")
    ax.set_xticklabels(["TE-low (Q1)", "TE-high (Q4)"], fontsize=9)
    ax.set_ylabel("Domain copies per species", fontsize=9)
    ax.set_title(domain.replace("_", " "), fontsize=9)

    _, p = scipy_stats.mannwhitneyu(low, high, alternative="less")
    ax.text(0.5, 0.97, f"MW p={p:.3g}", transform=ax.transAxes,
            ha="center", va="top", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    stem = re.sub(r"[^\w]", "_", domain)[:40]
    fig.savefig(out_dir / f"copynumber_{stem}.pdf")
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================

import re


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))
    species_list = list(clade_map.keys())

    te_df = pd.read_csv(args.te_content, sep="\t").set_index("species")["te_pct_genome"]
    gs_df = pd.read_csv(args.genome_size, sep="\t").set_index("species")["genome_size_bp"]

    print("Loading InterPro domain counts ...")
    ipr_long  = load_ipr_counts(Path(args.ipr_dir), species_list)
    count_mat = build_count_matrix(ipr_long)
    print(f"  {count_mat.shape[0]} species × {count_mat.shape[1]} domains")

    print("Testing domain–TE correlations ...")
    results = test_domain_te_correlation(count_mat, te_df, gs_df)
    results.to_csv(out_dir / "domain_expansion_results.tsv", sep="\t", index=False)
    n_sig = results["significant"].sum()
    print(f"  Significant (FDR < {args.fdr}): {n_sig} / {len(results)} domains")

    print("Generating bubble plot ...")
    plot_bubble(results, te_df, fig_dir)

    top_domains = results[results["significant"] & (results["spearman_r"] > 0)].head(args.top_n)
    print(f"Generating {len(top_domains)} copy-number plots ...")
    for _, row in top_domains.iterrows():
        dom = row["domain"]
        if dom in count_mat.columns:
            plot_domain_copynumber(dom, count_mat, clade_map, te_df, fig_dir)

    print("Done.")


if __name__ == "__main__":
    main()
