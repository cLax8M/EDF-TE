#!/usr/bin/env python3
"""
13_duf_te_fusions.py
DUF domain × TE domain co-occurrence analysis.

Identifies Domain of Unknown Function (DUF/UPF) domains that co-occur
with TE-derived domains in the same gene (potential neo-functionalization /
TE domestication events with unknown host functional output).

Logic:
  For each gene that contains ≥1 TE domain AND ≥1 DUF/UPF domain:
    → record the (DUF_id, TE_superfamily) pair
  Count across all species; rank by frequency; compute enrichment vs
  null expectation (DUF freq × TE freq).

Inputs:
  - Per-species domains.tsv   (238/238 — runs immediately)
  - config/TE_domains_master.tsv

Outputs:
  - figures/13_duf_te_fusions.tsv     (ranked co-occurrence table)
  - figures/13_duf_te_top_pairs.pdf   (lollipop chart of top 30 pairs)

Usage:
  python3 13_duf_te_fusions.py
"""

import os
import sys
import collections
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ─────────────────────────────────────────────────────────────────────
COOPTION_DIR   = "/legserv/Temp/Carlos/EDF-TE/results/Co-option"
DOMAINS_DIR    = os.path.join(COOPTION_DIR, "per_species", "domains")
TE_MASTER      = os.path.join(COOPTION_DIR, "config", "TE_domains_master.tsv")
OUT_DIR        = os.path.join(COOPTION_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_TSV = os.path.join(OUT_DIR, "13_duf_te_fusions.tsv")
OUT_FIG = os.path.join(OUT_DIR, "13_duf_te_top_pairs.pdf")

# ── helpers ───────────────────────────────────────────────────────────────────

def is_duf(domain_name, domain_id):
    """Return True if this domain is a DUF (domain of unknown function)."""
    name = str(domain_name).upper()
    did  = str(domain_id).upper()
    return (
        name.startswith("DUF")
        or name.startswith("UPF")
        or "DOMAIN OF UNKNOWN" in name
        or "UNCHARACTERISED" in name
        or "UNCHARACTERIZED" in name
    )


def load_te_master(path):
    """Return dict: domain_id → TE_superfamily, and set of TE domain IDs."""
    df = pd.read_csv(path, sep="\t")
    te_superfam = {}
    for _, row in df.iterrows():
        did = str(row["domain_id"]).strip().upper()
        sf  = str(row.get("TE_superfamily", "Unknown")).strip()
        te_superfam[did] = sf if sf and sf != "nan" else "Unknown"
    return te_superfam


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    print("[INFO] Loading TE domain master...")
    te_superfam = load_te_master(TE_MASTER)
    te_ids = set(te_superfam.keys())

    domain_files = sorted(
        f for f in os.listdir(DOMAINS_DIR) if f.endswith(".domains.tsv")
    )
    print(f"[INFO] Processing {len(domain_files)} species...")

    # Counters
    pair_counts   = collections.Counter()   # (duf_id, te_superf) → n_genes
    pair_species  = collections.defaultdict(set)   # pair → {species}
    duf_total     = collections.Counter()   # duf_id → n_genes (background)
    te_total_     = collections.Counter()   # te_superfam → n_genes (background)
    total_genes   = 0
    total_te_duf  = 0
    skipped       = []

    for fname in domain_files:
        sp = fname.replace(".domains.tsv", "")
        fpath = os.path.join(DOMAINS_DIR, fname)
        try:
            dom = pd.read_csv(fpath, sep="\t")
        except Exception as e:
            skipped.append(f"{sp}: {e}")
            continue

        # Group by gene
        for gene_id, grp in dom.groupby("gene_id"):
            domains_in_gene = grp[["domain_id", "domain_name"]].drop_duplicates()

            duf_ids_gene = set()
            te_sfs_gene  = set()

            for _, row in domains_in_gene.iterrows():
                did  = str(row["domain_id"]).strip().upper()
                dname = str(row["domain_name"]).strip()
                if did in te_ids:
                    te_sfs_gene.add(te_superfam[did])
                elif is_duf(dname, did):
                    duf_ids_gene.add(f"{did}|{dname}")

            # Background counts
            for d in duf_ids_gene:
                duf_total[d] += 1
            for t in te_sfs_gene:
                te_total_[t] += 1
            total_genes += 1

            # Co-occurrence pairs
            if duf_ids_gene and te_sfs_gene:
                total_te_duf += 1
                for d in duf_ids_gene:
                    for t in te_sfs_gene:
                        pair_counts[(d, t)] += 1
                        pair_species[(d, t)].add(sp)

    print(f"  {total_genes} total genes processed")
    print(f"  {total_te_duf} genes have BOTH DUF and TE domains")
    print(f"  {len(pair_counts)} unique DUF × TE-superfamily pairs")

    if not pair_counts:
        print("[WARN] No DUF+TE co-occurrences found. Exiting.")
        sys.exit(0)

    # Build results table with enrichment stats
    rows = []
    duf_bg_total = sum(duf_total.values())
    te_bg_total  = sum(te_total_.values())

    for (duf_did, te_sf), n_pair in pair_counts.items():
        dname_part = duf_did.split("|", 1)[1] if "|" in duf_did else duf_did
        did_part   = duf_did.split("|", 1)[0]

        n_duf_total = duf_total.get(duf_did, 0)
        n_te_total  = te_total_.get(te_sf, 0)

        # Fisher test: co-occurrence vs background
        # 2×2: [genes with BOTH, genes with DUF only, genes with TE only, genes with neither]
        n_both      = n_pair
        n_duf_only  = n_duf_total - n_pair
        n_te_only   = n_te_total  - n_pair
        n_neither   = total_genes - n_duf_total - n_te_total + n_pair

        if n_duf_only < 0:  n_duf_only = 0
        if n_te_only  < 0:  n_te_only  = 0
        if n_neither  < 0:  n_neither  = 0

        tbl = [[n_both, n_duf_only], [n_te_only, n_neither]]
        try:
            or_, pv = fisher_exact(tbl, alternative="greater")
        except Exception:
            or_, pv = float("nan"), 1.0

        rows.append({
            "duf_domain_id":    did_part,
            "duf_domain_name":  dname_part,
            "te_superfamily":   te_sf,
            "n_genes_copresent": n_both,
            "n_species":        len(pair_species[(duf_did, te_sf)]),
            "n_genes_with_duf": n_duf_total,
            "n_genes_with_te":  n_te_total,
            "odds_ratio":       round(or_, 4),
            "pval_fisher":      pv,
        })

    df = pd.DataFrame(rows)
    if len(df) > 0:
        _, padj, _, _ = multipletests(df["pval_fisher"].values, method="fdr_bh")
        df["padj_fisher"] = padj
    df = df.sort_values(["n_genes_copresent", "n_species"], ascending=False).reset_index(drop=True)

    df.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")
    print(f"[OUT] DUF×TE co-occurrence table → {OUT_TSV}")
    sig = (df["padj_fisher"] < 0.05) if "padj_fisher" in df.columns else pd.Series([False])
    print(f"      {int(sig.sum())} pairs significant after FDR<0.05 correction")

    # ── figure ────────────────────────────────────────────────────────────────
    top_n = min(30, len(df))
    df_top = df.head(top_n)

    fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.35)))
    labels = [
        f"{row['duf_domain_name']} + {row['te_superfamily']}"
        for _, row in df_top.iterrows()
    ]
    values = df_top["n_genes_copresent"].values
    y_pos  = np.arange(len(labels))

    # Color bars by TE superfamily
    superfams = df_top["te_superfamily"].unique()
    cmap = plt.get_cmap("tab20", len(superfams))
    sf_colors = {sf: cmap(i) for i, sf in enumerate(superfams)}
    bar_colors = [sf_colors[row["te_superfamily"]] for _, row in df_top.iterrows()]

    ax.barh(y_pos, values, color=bar_colors, height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Number of genes with DUF+TE co-occurrence")
    ax.set_title(f"Top {top_n} DUF × TE-superfamily co-occurrence pairs\n"
                 f"({len(domain_files)} species, {total_te_duf} fusion genes total)")
    # Legend for TE superfamily colors
    from matplotlib.patches import Patch
    legend_els = [Patch(facecolor=sf_colors[sf], label=sf) for sf in superfams]
    ax.legend(handles=legend_els, loc="lower right", fontsize=7, title="TE superfamily")

    plt.tight_layout()
    plt.savefig(OUT_FIG, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"[OUT] Figure → {OUT_FIG}")
    print("[DONE]")


if __name__ == "__main__":
    main()
