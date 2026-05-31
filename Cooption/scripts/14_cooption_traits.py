#!/usr/bin/env python3
"""
14_cooption_traits.py
Per-species co-option burden table for tree visualization.

Computes per-species metrics that can be mapped as tip traits on the species
tree (e.g. in iTOL, FigTree, or R ggtree plots).

Metrics per species:
  - n_coopted            : number of TE-coopted candidate genes
  - n_total_genes        : total annotated genes (from gene_classes.tsv)
  - pct_coopted          : n_coopted / n_total_genes * 100
  - n_te_core            : genes with ≥1 TE_core domain
  - n_te_assoc           : genes with ONLY TE_assoc domains (no TE_core)
  - dominant_te_superfam : most common TE superfamily in coopted genes
  - dominant_host_domain : most common non-TE domain in coopted genes
  - top_architecture     : most common domain architecture (TE + host)
  - genome_size_mb       : from GenomicFeatures summary (if available)
  - pct_genome_TE        : % genome covered by TEs (from summary, if available)

Inputs (238/238 — runs immediately):
  - per_species/classification/{SP}.true_TE_coopted_genes.tsv
  - per_species/classification/{SP}.gene_classes.tsv
  - per_species/domains/{SP}.domains.tsv
  - config/TE_domains_master.tsv
  - results/GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY.txt  (or V2)

Outputs:
  - figures/14_cooption_traits.tsv     (main output for tree annotation)
  - figures/14_cooption_traits.pdf     (overview scatter plots)

Usage:
  python3 14_cooption_traits.py
"""

import os
import sys
import re
import collections
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ─────────────────────────────────────────────────────────────────────
DATA_BASE      = "/legserv/Temp/Carlos/EDF-TE"
COOPTION_DIR   = os.path.join(DATA_BASE, "results", "Co-option")
CLASS_DIR      = os.path.join(COOPTION_DIR, "per_species", "classification")
DOMAINS_DIR    = os.path.join(COOPTION_DIR, "per_species", "domains")
TE_MASTER      = os.path.join(COOPTION_DIR, "config", "TE_domains_master.tsv")

# GenomicFeatures: try V3 first, fall back to V2
GENOMIC_FEAT_PATHS = [
    os.path.join(DATA_BASE, "results", "GenomicFeatures_V3", "ALL_SPECIES_ALL_SUMMARY.txt"),
    os.path.join(DATA_BASE, "results", "GenomicFeatures_V2", "ALL_SPECIES_ALL_SUMMARY.txt"),
]

OUT_DIR = os.path.join(COOPTION_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_TSV = os.path.join(OUT_DIR, "14_cooption_traits.tsv")
OUT_FIG = os.path.join(OUT_DIR, "14_cooption_traits.pdf")

# ── helpers ───────────────────────────────────────────────────────────────────

def load_te_master(path):
    """Return dict: domain_id → TE_superfamily; set of TE_core ids; set of TE_assoc ids."""
    df = pd.read_csv(path, sep="\t")
    superfam   = {}
    te_core_ids  = set()
    te_assoc_ids = set()
    for _, row in df.iterrows():
        did = str(row["domain_id"]).strip().upper()
        sf  = str(row.get("TE_superfamily", "Unknown")).strip()
        superfam[did] = sf if sf and sf != "nan" else "Unknown"
        status = str(row.get("TE_status", "")).strip().upper()
        if "CORE" in status:
            te_core_ids.add(did)
        elif "ASSOC" in status:
            te_assoc_ids.add(did)
    return superfam, te_core_ids, te_assoc_ids


def load_genomic_features(paths):
    """Load ALL_SPECIES_ALL_SUMMARY.txt; returns dict: sp → {col: value}."""
    for p in paths:
        if os.path.exists(p):
            try:
                df = pd.read_csv(p, sep="\t")
                # normalise species column name
                id_col = next((c for c in df.columns
                               if c.lower() in ("species", "species_id", "sp")), None)
                if id_col:
                    df = df.rename(columns={id_col: "species"})
                    return df.set_index("species").to_dict(orient="index")
            except Exception:
                pass
    return {}


def most_common(lst, default="None"):
    if not lst:
        return default
    c = collections.Counter(lst)
    return c.most_common(1)[0][0]


def architecture_string(te_sfs, host_domains):
    """Build a compact label: TE1+TE2_HOST1+HOST2"""
    if not te_sfs and not host_domains:
        return "unknown"
    te_part   = "+".join(sorted(set(te_sfs)))   if te_sfs     else "none"
    host_part = "+".join(sorted(set(host_domains))) if host_domains else "none"
    return f"{te_part}_{host_part}"


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    print("[INFO] Loading TE domain master...")
    superfam, te_core_ids, te_assoc_ids = load_te_master(TE_MASTER)
    all_te_ids = te_core_ids | te_assoc_ids

    print("[INFO] Loading genomic features...")
    gf_data = load_genomic_features(GENOMIC_FEAT_PATHS)
    print(f"  {len(gf_data)} species in GenomicFeatures summary")

    # Discover species from classification dir
    all_files = os.listdir(CLASS_DIR)
    species_list = sorted(set(
        f.replace(".true_TE_coopted_genes.tsv", "")
        for f in all_files
        if f.endswith(".true_TE_coopted_genes.tsv")
    ))
    print(f"[INFO] {len(species_list)} species with classification data")

    rows = []
    skipped = []

    for sp in species_list:
        coopted_path = os.path.join(CLASS_DIR, f"{sp}.true_TE_coopted_genes.tsv")
        classes_path = os.path.join(CLASS_DIR, f"{sp}.gene_classes.tsv")
        domains_path = os.path.join(DOMAINS_DIR, f"{sp}.domains.tsv")

        # Gene-class counts
        if not os.path.exists(classes_path):
            skipped.append(f"{sp} (no gene_classes.tsv)")
            continue
        try:
            cls_df = pd.read_csv(classes_path, sep="\t")
            n_total = len(cls_df)
        except Exception as e:
            skipped.append(f"{sp} (read error: {e})")
            continue

        # Coopted genes
        if not os.path.exists(coopted_path):
            n_coopted = 0
            df_coopted = pd.DataFrame()
        else:
            try:
                df_coopted = pd.read_csv(coopted_path, sep="\t")
                n_coopted = len(df_coopted)
            except Exception:
                n_coopted = 0
                df_coopted = pd.DataFrame()

        pct_coopted = round(100.0 * n_coopted / n_total, 4) if n_total > 0 else 0

        # TE_core vs TE_assoc counts
        n_te_core_genes  = 0
        n_te_assoc_genes = 0
        if not df_coopted.empty and "n_TE_core" in df_coopted.columns:
            n_te_core_genes  = int((df_coopted["n_TE_core"]  > 0).sum())
            n_te_assoc_genes = int(
                ((df_coopted["n_TE_assoc"] > 0) & (df_coopted["n_TE_core"] == 0)).sum()
            )

        # Dominant TE superfamily and host domain (from domains.tsv)
        dominant_te_sf      = "None"
        dominant_host_dom   = "None"
        top_architecture    = "None"
        n_unique_te_sfs     = 0

        if os.path.exists(domains_path) and not df_coopted.empty:
            try:
                dom_df = pd.read_csv(domains_path, sep="\t")
                coopted_genes = set(df_coopted["gene_id"].astype(str))
                dom_coopted = dom_df[dom_df["gene_id"].astype(str).isin(coopted_genes)]

                te_sfs_all      = []
                host_doms_all   = []
                arch_per_gene   = []

                for gene_id, grp in dom_coopted.groupby("gene_id"):
                    gene_te_sfs   = []
                    gene_host_dom = []
                    for _, row in grp.iterrows():
                        did = str(row["domain_id"]).strip().upper()
                        if did in all_te_ids:
                            sf = superfam.get(did, "Unknown")
                            gene_te_sfs.append(sf)
                            te_sfs_all.append(sf)
                        else:
                            # non-TE domain → host domain
                            dname = str(row.get("domain_name", "")).strip()
                            gene_host_dom.append(dname)
                            host_doms_all.append(dname)
                    arch_per_gene.append(architecture_string(gene_te_sfs, gene_host_dom))

                dominant_te_sf    = most_common(te_sfs_all)
                dominant_host_dom = most_common(host_doms_all)
                top_architecture  = most_common(arch_per_gene)
                n_unique_te_sfs   = len(set(te_sfs_all))

            except Exception:
                pass

        # Genomic features
        gf = gf_data.get(sp, {})
        genome_size_mb = gf.get("genome_size_mb", gf.get("Genome_Size_Mb",
                          gf.get("genome_size", np.nan)))
        pct_te_genome  = gf.get("pct_TE", gf.get("pct_genome_TE",
                          gf.get("TE_percent", np.nan)))

        rows.append({
            "species":             sp,
            "n_total_genes":       n_total,
            "n_coopted":           n_coopted,
            "pct_coopted":         pct_coopted,
            "n_te_core_genes":     n_te_core_genes,
            "n_te_assoc_genes":    n_te_assoc_genes,
            "n_unique_te_sfs":     n_unique_te_sfs,
            "dominant_te_superfam":dominant_te_sf,
            "dominant_host_domain":dominant_host_dom,
            "top_architecture":    top_architecture,
            "genome_size_mb":      genome_size_mb,
            "pct_genome_TE":       pct_te_genome,
        })

    df_out = pd.DataFrame(rows)
    df_out = df_out.sort_values("pct_coopted", ascending=False)
    df_out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")
    print(f"[OUT] Co-option traits table → {OUT_TSV}")
    print(f"      {len(rows)} species | {df_out['n_coopted'].sum()} total coopted genes")
    print(f"      pct_coopted range: {df_out['pct_coopted'].min():.3f}%–{df_out['pct_coopted'].max():.3f}%")

    if skipped:
        print(f"\n[WARN] Skipped {len(skipped)} species:")
        for s in skipped:
            print(f"  {s}")

    # ── figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # Panel A: histogram of pct_coopted across species
    ax = axes[0]
    ax.hist(df_out["pct_coopted"].dropna(), bins=30, color="#ff7f0e", edgecolor="white")
    ax.set_xlabel("% TE-coopted genes per species")
    ax.set_ylabel("Number of species")
    ax.set_title(f"Co-option burden distribution\n(n={len(df_out)} species)")

    # Panel B: pct_coopted vs pct_genome_TE (if available)
    ax2 = axes[1]
    sub = df_out.dropna(subset=["pct_coopted", "pct_genome_TE"])
    if len(sub) > 5:
        ax2.scatter(sub["pct_genome_TE"], sub["pct_coopted"],
                    alpha=0.7, s=30, color="#d62728")
        # Add regression line
        m, b = np.polyfit(sub["pct_genome_TE"], sub["pct_coopted"], 1)
        xs = np.linspace(sub["pct_genome_TE"].min(), sub["pct_genome_TE"].max(), 100)
        ax2.plot(xs, m * xs + b, "k--", linewidth=0.9)
        corr = sub[["pct_genome_TE", "pct_coopted"]].corr().iloc[0, 1]
        ax2.set_title(f"Co-option % vs TE genome %\n(r={corr:.3f}, n={len(sub)})")
    else:
        ax2.text(0.5, 0.5, "Insufficient genomic\nfeature data",
                 ha="center", va="center", transform=ax2.transAxes)
        ax2.set_title("Co-option % vs TE genome %")
    ax2.set_xlabel("% genome covered by TEs")
    ax2.set_ylabel("% TE-coopted genes")

    # Panel C: top dominant TE superfamilies
    ax3 = axes[2]
    sf_counts = df_out["dominant_te_superfam"] \
        .value_counts() \
        .drop("None", errors="ignore") \
        .head(15)
    if len(sf_counts) > 0:
        sf_counts.plot.barh(ax=ax3, color="#1f77b4")
        ax3.set_xlabel("Number of species (dominant TE SF)")
        ax3.set_title("Dominant TE superfamily\nin co-opted genes")
        ax3.invert_yaxis()
    else:
        ax3.text(0.5, 0.5, "No data", ha="center", transform=ax3.transAxes)

    plt.tight_layout()
    plt.savefig(OUT_FIG, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"[OUT] Figure → {OUT_FIG}")
    print("[DONE]")


if __name__ == "__main__":
    main()
