#!/usr/bin/env python3
"""
D03_pgls.py — Phylogenetic Generalised Least Squares (extended)
================================================================
Extends D02_pgls_regression.R with:
  • Additional TE metrics: Shannon diversity (H', J), Kimura (recency,
    burst, mean divergence), TE composition proportions.
  • Additional genomic predictors: Pct_Unspliced_Genes.
  • Both directions: TE metric ~ genomic feature  AND  genomic ~ TE metric.
  • BH-FDR correction alongside Bonferroni.
  • Significance heatmap in PDF.

Requires: pandas, numpy, scipy, statsmodels, matplotlib, seaborn.
Tree is parsed from Newick without ete3/dendropy using a built-in parser.

Usage
-----
  python D03_pgls.py [options]

  --merged   D01_te_genomics_merged_provisional.tsv
  --shannon  B01_shannon_diversity.tsv
  --kimura   C01_kimura_species.tsv
  --summary  ALL_SPECIES_ALL_SUMMARY_V4.txt   (only for Pct_Unspliced_Genes)
  --tree     SpeciesTree_rooted.txt   (Newick, tips may have .aa suffix)
  --outdir   output directory

Outputs
-------
  <outdir>/pgls_extended_results.tsv   — all pairwise PGLS results
  <outdir>/pgls_heatmap.pdf            — signed-log10(p) heatmap
"""

import argparse
import os
import sys
import re
import warnings
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

warnings.filterwarnings("ignore", category=FutureWarning)

# ── paths ─────────────────────────────────────────────────────────────────────
_DATA = "/legserv/Temp/Carlos/EDF-TE"
_SEC2 = f"{_DATA}/results/Section2_TElandscape"
_GF   = f"{_DATA}/results/GenomicFeatures_V3"
_TREE = f"{_DATA}/results/Orthofinder_run/Species_Tree/SpeciesTree_rooted.txt"


def parse_args():
    p = argparse.ArgumentParser(description="Extended PGLS: TE landscape × genomic architecture")
    p.add_argument("--merged",  default=f"{_SEC2}/D01_te_genomics_merged_provisional.tsv")
    p.add_argument("--shannon", default=f"{_SEC2}/B01_shannon_diversity.tsv")
    p.add_argument("--kimura",  default=f"{_SEC2}/C01_kimura_species.tsv")
    p.add_argument("--summary", default=f"{_GF}/ALL_SPECIES_ALL_SUMMARY_V4.txt")
    p.add_argument("--tree",    default=_TREE)
    p.add_argument("--outdir",  default=f"{_SEC2}/D03_pgls")
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
# MINIMAL NEWICK PARSER → VCV MATRIX (Brownian Motion)
# ═══════════════════════════════════════════════════════════════════════════════

class _Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children: list["_Node"] = []
        self.parent: Optional["_Node"] = None

    def is_leaf(self):
        return len(self.children) == 0

    def add_child(self, child: "_Node"):
        child.parent = self
        self.children.append(child)


def _parse_newick(nwk: str) -> _Node:
    """Recursive-descent parser for a Newick string.  Returns root node."""
    nwk = nwk.strip().rstrip(";").strip()
    pos = [0]  # mutable index

    def _read_node() -> _Node:
        node = _Node()
        if pos[0] < len(nwk) and nwk[pos[0]] == "(":
            pos[0] += 1  # skip '('
            child = _read_node()
            node.add_child(child)
            while pos[0] < len(nwk) and nwk[pos[0]] == ",":
                pos[0] += 1
                child = _read_node()
                node.add_child(child)
            if pos[0] < len(nwk) and nwk[pos[0]] == ")":
                pos[0] += 1
        # read optional label:branch_length
        label = []
        while pos[0] < len(nwk) and nwk[pos[0]] not in (",", ")", ";", "("):
            label.append(nwk[pos[0]])
            pos[0] += 1
        label_str = "".join(label).strip()
        if ":" in label_str:
            parts = label_str.rsplit(":", 1)
            node.name = parts[0].strip()
            try:
                node.length = float(parts[1])
            except ValueError:
                node.length = 0.0
        else:
            node.name = label_str
        return node

    root = _read_node()
    return root


def _get_leaves(node: _Node) -> list[_Node]:
    if node.is_leaf():
        return [node]
    leaves = []
    for c in node.children:
        leaves.extend(_get_leaves(c))
    return leaves


def _root_to_node_dist(node: _Node) -> float:
    """Distance from the root to this node (sum of branch lengths up)."""
    d = 0.0
    n = node
    while n.parent is not None:
        d += n.length
        n = n.parent
    return d


def _compute_vcv(root: _Node, tip_names: list[str]) -> np.ndarray:
    """Compute the BM variance–covariance matrix.

    V[i,j] = shared path length from the root to the MRCA of tips i and j.
    """
    n = len(tip_names)
    tip_idx = {name: i for i, name in enumerate(tip_names)}

    # For each node, collect the set of tip indices below it
    def _tip_set(node: _Node) -> set[int]:
        if node.is_leaf():
            idx = tip_idx.get(node.name)
            return {idx} if idx is not None else set()
        s = set()
        for c in node.children:
            s |= _tip_set(c)
        return s

    # Build VCV by traversing and accumulating shared branch time
    V = np.zeros((n, n), dtype=np.float64)

    def _add_bl(node: _Node):
        tips_below = _tip_set(node)
        bl = node.length if node.parent is not None else 0.0
        if bl > 0 and tips_below:
            idx_list = sorted(tips_below)
            for i in idx_list:
                for j in idx_list:
                    V[i, j] += bl
        for c in node.children:
            _add_bl(c)

    _add_bl(root)
    return V


def load_tree(path: str) -> tuple[_Node, list[str]]:
    """Read Newick, strip .aa suffix, return (root, tip_names)."""
    with open(path) as f:
        raw = f.read().strip()
    root = _parse_newick(raw)
    # strip .aa / .faa suffix
    for leaf in _get_leaves(root):
        leaf.name = re.sub(r"\.(aa|faa)$", "", leaf.name)
    tip_names = [l.name for l in _get_leaves(root)]
    return root, tip_names


# ═══════════════════════════════════════════════════════════════════════════════
# PGLS via statsmodels GLS
# ═══════════════════════════════════════════════════════════════════════════════

def pgls(y: np.ndarray, x: np.ndarray, sigma: np.ndarray):
    """Run GLS:  y = a + b*x + ε,  Cov(ε) = σ² Σ.

    Returns (estimate, se, t, p, n).
    """
    import statsmodels.api as sm

    X = sm.add_constant(x)
    try:
        model = sm.GLS(y, X, sigma=sigma)
        result = model.fit()
        # slope is index 1
        est = result.params[1]
        se  = result.bse[1]
        t   = result.tvalues[1]
        p   = result.pvalues[1]
        return est, se, t, p, len(y)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, len(y)


# ═══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════════

def load_data(args) -> pd.DataFrame:
    """Merge all data sources into one table indexed by species."""
    # 1. D01 merged (TE% + genomic features)
    merged = pd.read_csv(args.merged, sep="\t")

    # 2. Shannon diversity
    if os.path.isfile(args.shannon):
        shan = pd.read_csv(args.shannon, sep="\t")
        shan = shan[["species", "H_prime", "J_evenness"]].copy()
        merged = merged.merge(shan, on="species", how="left")

    # 3. Kimura species-level metrics
    if os.path.isfile(args.kimura):
        kim = pd.read_csv(args.kimura, sep="\t")
        kim = kim[["species", "recency_index", "burst_score", "mean_div_weighted"]].copy()
        merged = merged.merge(kim, on="species", how="left")

    # 4. Pct_Unspliced_Genes from V4 summary
    if os.path.isfile(args.summary):
        v4 = pd.read_csv(args.summary, sep="\t")
        if "Pct_Unspliced_Genes" in v4.columns:
            v4 = v4[["Species", "Pct_Unspliced_Genes"]].rename(
                columns={"Species": "species"})
            merged = merged.merge(v4, on="species", how="left")

    merged.set_index("species", inplace=True)
    print(f"[DATA] {merged.shape[0]} species × {merged.shape[1]} variables")
    return merged


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

TE_VARS = [
    "pct_TE_total", "pct_DNA", "pct_LTR", "pct_LINE", "pct_SINE",
    "pct_RC", "pct_PLE", "pct_Unknown", "pct_Other",
    "H_prime", "J_evenness",
    "recency_index", "burst_score", "mean_div_weighted",
]

GENOMIC_VARS = [
    "Genome_Size_bp", "Total_Genes", "Avg_Exons_Per_Gene",
    "Avg_Introns_Multiexon", "Pct_Genic", "Pct_Exonic",
    "Pct_Intronic", "Pct_Intergenic", "Pct_Unspliced_Genes",
]


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ── load data ─────────────────────────────────────────────────────────────
    df = load_data(args)
    root, tip_names = load_tree(args.tree)
    tip_set = set(tip_names)

    # intersect species
    common = sorted(set(df.index) & tip_set)
    if len(common) < 20:
        sys.exit(f"[ERROR] Only {len(common)} species in common between data and tree.")
    print(f"[TREE] {len(common)} species in common")

    df = df.loc[common]

    # prune tree: rebuild VCV for common species only
    # (keep only leaves in common; reconstruct)
    # To avoid re-implementing tree pruning, compute full VCV and subselect
    full_tip_idx = {name: i for i, name in enumerate(tip_names)}
    V_full = _compute_vcv(root, tip_names)
    common_idx = [full_tip_idx[s] for s in common]
    V = V_full[np.ix_(common_idx, common_idx)]

    # regularise: small ridge to ensure positive definiteness
    V += np.eye(V.shape[0]) * 1e-6 * np.mean(np.diag(V))

    # ── filter available variables ────────────────────────────────────────────
    te_cols = [c for c in TE_VARS if c in df.columns]
    gen_cols = [c for c in GENOMIC_VARS if c in df.columns]
    print(f"[VARS] TE variables:      {len(te_cols)}")
    print(f"[VARS] Genomic variables: {len(gen_cols)}")

    # ── run all pairwise PGLS ─────────────────────────────────────────────────
    records = []

    # Direction 1: TE ~ genomic
    for response in te_cols:
        for predictor in gen_cols:
            y_raw = df[response].values.astype(float)
            x_raw = df[predictor].values.astype(float)
            keep = np.isfinite(y_raw) & np.isfinite(x_raw)
            if keep.sum() < 15:
                continue
            V_sub = V[np.ix_(keep, keep)]
            est, se, t, p, n = pgls(y_raw[keep], x_raw[keep], V_sub)
            records.append({
                "response": response, "predictor": predictor,
                "direction": "TE~genomic",
                "estimate": est, "se": se, "t_value": t, "p_raw": p, "n": n,
            })

    # Direction 2: genomic ~ TE
    for response in gen_cols:
        for predictor in te_cols:
            y_raw = df[response].values.astype(float)
            x_raw = df[predictor].values.astype(float)
            keep = np.isfinite(y_raw) & np.isfinite(x_raw)
            if keep.sum() < 15:
                continue
            V_sub = V[np.ix_(keep, keep)]
            est, se, t, p, n = pgls(y_raw[keep], x_raw[keep], V_sub)
            records.append({
                "response": response, "predictor": predictor,
                "direction": "genomic~TE",
                "estimate": est, "se": se, "t_value": t, "p_raw": p, "n": n,
            })

    if not records:
        sys.exit("[ERROR] No PGLS tests completed.")

    out = pd.DataFrame(records)

    # ── multiple testing correction ───────────────────────────────────────────
    valid = out["p_raw"].notna()
    n_tests = valid.sum()
    out.loc[valid, "p_bonferroni"] = np.minimum(out.loc[valid, "p_raw"] * n_tests, 1.0)

    # BH-FDR
    pvals = out.loc[valid, "p_raw"].values
    order = np.argsort(pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(pvals) + 1)
    fdr = np.minimum(pvals * len(pvals) / ranks, 1.0)
    # enforce monotonicity (decreasing in rank order)
    fdr_sorted = fdr[order]
    for i in range(len(fdr_sorted) - 2, -1, -1):
        fdr_sorted[i + 1] = min(fdr_sorted[i + 1], fdr_sorted[i]) if i + 1 < len(fdr_sorted) else fdr_sorted[i + 1]
    fdr[order] = fdr_sorted
    out.loc[valid, "p_fdr"] = fdr

    out["sig_bonf"] = out["p_bonferroni"] < 0.05
    out["sig_fdr"]  = out["p_fdr"] < 0.05
    out.sort_values("p_raw", inplace=True)

    # ── write results ─────────────────────────────────────────────────────────
    tsv_path = os.path.join(args.outdir, "pgls_extended_results.tsv")
    out.to_csv(tsv_path, sep="\t", index=False, float_format="%.6g")
    n_bonf = out["sig_bonf"].sum()
    n_fdr  = out["sig_fdr"].sum()
    print(f"\n[RESULTS] {len(out)} tests | {n_bonf} significant (Bonferroni) | {n_fdr} significant (FDR)")
    print(f"[SAVED]   {tsv_path}")

    # ── heatmap ───────────────────────────────────────────────────────────────
    _plot_heatmap(out, args.outdir, te_cols, gen_cols)
    print("[DONE]")


# ═══════════════════════════════════════════════════════════════════════════════
# HEATMAP
# ═══════════════════════════════════════════════════════════════════════════════

def _plot_heatmap(results: pd.DataFrame, outdir: str,
                  te_cols: list[str], gen_cols: list[str]):
    """Signed −log10(p) heatmap for direction 'TE~genomic'."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sub = results[results["direction"] == "TE~genomic"].copy()
    if sub.empty:
        return

    # signed -log10(p_fdr)
    sub["signed_logp"] = -np.log10(sub["p_fdr"].clip(lower=1e-300)) * np.sign(sub["estimate"])

    pivot = sub.pivot_table(index="response", columns="predictor",
                            values="signed_logp", aggfunc="first")
    # reorder
    row_order = [c for c in te_cols if c in pivot.index]
    col_order = [c for c in gen_cols if c in pivot.columns]
    pivot = pivot.reindex(index=row_order, columns=col_order)

    # stars for significance
    annot_df = sub.pivot_table(index="response", columns="predictor",
                               values="sig_fdr", aggfunc="first")
    annot_df = annot_df.reindex(index=row_order, columns=col_order)
    annot = annot_df.map(lambda x: "*" if x else "")

    vmax = min(pivot.abs().max().max(), 50)

    fig, ax = plt.subplots(figsize=(1.0 + 0.8 * len(col_order),
                                    0.8 + 0.6 * len(row_order)))
    sns.heatmap(pivot, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
                annot=annot, fmt="", linewidths=0.5, ax=ax,
                cbar_kws={"label": "signed −log₁₀(FDR p)"})
    ax.set_title("PGLS: TE metric ~ genomic feature", fontsize=11)
    ax.set_ylabel("TE response")
    ax.set_xlabel("Genomic predictor")
    plt.tight_layout()

    pdf_path = os.path.join(outdir, "pgls_heatmap.pdf")
    fig.savefig(pdf_path, dpi=150)
    plt.close(fig)
    print(f"[PLOT]    {pdf_path}")


if __name__ == "__main__":
    main()
