#!/usr/bin/env python3
"""
11_orphan_gene_TE_domains.py
Orphan gene × TE domain enrichment analysis.

Tests whether orphan genes (genes not assigned to any orthogroup by OrthoFinder)
are enriched for TE-derived domains compared to non-orphan genes.

Inputs (read automatically from standard paths):
  - OrthoFinder Orthogroups_UnassignedGenes.tsv  → orphan gene IDs per species
  - Per-species domains.tsv                       → domain annotations
  - TE_domains_master.tsv                         → TE domain classification
  - results/unknown_function/unknown_function_summary.tsv → unknown% context

Outputs:
  - results/Co-option/figures/11_orphan_te_enrichment.tsv   (per-species stats)
  - results/Co-option/figures/11_orphan_global_summary.tsv  (global Fisher test)
  - results/Co-option/figures/11_orphan_te_enrichment.pdf   (figure)

Design notes:
  - Skips species with missing domains.tsv (not yet processed).
  - Re-runnable: when all 238 EarlGrey finish and domains.tsv are present,
    just re-run to include previously missing species.
  - "Orphan gene" = gene not in any orthogroup in OrthoFinder output.
  - For each species: Fisher test (TE domain presence vs orphan status).
  - Global test: genes pooled across all species.

Usage:
  python3 11_orphan_gene_TE_domains.py

All paths are set below and auto-resolved from COOPTION_DIR and ORTHOFINDER_DIR.
"""

import os
import sys
import warnings
import collections
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── paths ─────────────────────────────────────────────────────────────────────
COOPTION_DIR   = "/legserv/Temp/Carlos/EDF-TE/results/Co-option"
DOMAINS_DIR    = os.path.join(COOPTION_DIR, "per_species", "domains")
TE_MASTER      = os.path.join(COOPTION_DIR, "config", "TE_domains_master.tsv")
OG_UNASSIGNED  = ("/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/"
                  "Orthogroups/Orthogroups_UnassignedGenes.tsv")
UNKNOWN_FN_TSV = "/legserv/Temp/Carlos/EDF-TE/results/unknown_function/unknown_function_summary.tsv"
OUT_DIR        = os.path.join(COOPTION_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PER_SP  = os.path.join(OUT_DIR, "11_orphan_te_enrichment.tsv")
OUT_GLOBAL  = os.path.join(OUT_DIR, "11_orphan_global_summary.tsv")
OUT_FIG     = os.path.join(OUT_DIR, "11_orphan_te_enrichment.pdf")

# ── helpers ───────────────────────────────────────────────────────────────────

def load_te_master(path):
    """Return set of TE domain IDs and set of TE domain names (both case-insensitive)."""
    df = pd.read_csv(path, sep="\t")
    te_ids   = set(df["domain_id"].str.strip().str.upper())
    te_names = set(df["domain_name"].str.strip().str.upper())
    return te_ids, te_names


def load_orphan_genes(path):
    """
    Parse Orthogroups_UnassignedGenes.tsv.
    Returns dict: species_code → set of orphan gene IDs.
    Species column headers look like 'AP_Amoocci.aa'.
    """
    warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    orphans = {}
    for col in df.columns:
        sp = col.replace(".aa", "")
        genes = set(df[col].dropna().astype(str).str.strip())
        if genes:
            orphans[sp] = genes
    return orphans


def has_te_domain(domain_ids, domain_names, te_ids, te_names):
    """Return True if any domain in lists is TE-derived."""
    for d in domain_ids:
        if str(d).upper() in te_ids:
            return True
    for n in domain_names:
        if str(n).upper() in te_names:
            return True
    return False


def process_species(sp, domains_file, orphan_set, te_ids, te_names):
    """
    For one species, classify each gene as orphan/non-orphan × TE-domain / not.
    Returns 2×2 contingency table and summary dict.
    """
    try:
        dom = pd.read_csv(domains_file, sep="\t")
    except Exception as e:
        return None, None

    # Group domains per gene
    grouped = dom.groupby("gene_id").agg(
        domain_ids=("domain_id", list),
        domain_names=("domain_name", list)
    ).reset_index()

    n_orphan_te   = 0
    n_orphan_note = 0
    n_other_te    = 0
    n_other_note  = 0

    for _, row in grouped.iterrows():
        gene  = row["gene_id"]
        te    = has_te_domain(row["domain_ids"], row["domain_names"], te_ids, te_names)
        is_orphan = gene in orphan_set
        if is_orphan:
            if te:   n_orphan_te   += 1
            else:    n_orphan_note += 1
        else:
            if te:   n_other_te    += 1
            else:    n_other_note  += 1

    total_genes  = len(grouped)
    n_orphan     = n_orphan_te + n_orphan_note
    n_te_total   = n_orphan_te + n_other_te

    if n_orphan == 0 or total_genes == 0:
        return None, None

    # Fisher exact: rows = orphan/other, cols = TE/not-TE
    table = [[n_orphan_te, n_orphan_note],
             [n_other_te,  n_other_note]]
    oddsratio, pval = fisher_exact(table, alternative="two-sided")

    pct_orphan_te = 100.0 * n_orphan_te / n_orphan if n_orphan > 0 else 0
    pct_other_te  = 100.0 * n_other_te  / (total_genes - n_orphan) if (total_genes - n_orphan) > 0 else 0

    summary = {
        "species":         sp,
        "total_genes_dom": total_genes,
        "n_orphan":        n_orphan,
        "n_orphan_te":     n_orphan_te,
        "pct_orphan_te":   round(pct_orphan_te, 3),
        "n_other_te":      n_other_te,
        "pct_other_te":    round(pct_other_te, 3),
        "odds_ratio":      round(oddsratio, 4),
        "pval_fisher":     pval,
    }
    return summary, table


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    print("[INFO] Loading TE domain master list...")
    te_ids, te_names = load_te_master(TE_MASTER)
    print(f"  {len(te_ids)} TE domain IDs, {len(te_names)} TE domain names")

    print("[INFO] Loading orphan genes from OrthoFinder...")
    orphans = load_orphan_genes(OG_UNASSIGNED)
    print(f"  {len(orphans)} species with orphan gene data")

    # Collect per-species results
    results   = []
    skipped   = []
    global_table = [0, 0, 0, 0]  # orphan_te, orphan_note, other_te, other_note

    domain_files = {
        os.path.basename(f).replace(".domains.tsv", ""): os.path.join(DOMAINS_DIR, f)
        for f in os.listdir(DOMAINS_DIR) if f.endswith(".domains.tsv")
    }

    print(f"[INFO] Processing {len(domain_files)} species with domains.tsv...")
    for sp, dom_path in sorted(domain_files.items()):
        orphan_set = orphans.get(sp, set())
        summary, table = process_species(sp, dom_path, orphan_set, te_ids, te_names)
        if summary is None:
            skipped.append(sp)
            continue
        results.append(summary)
        global_table[0] += table[0][0]
        global_table[1] += table[0][1]
        global_table[2] += table[1][0]
        global_table[3] += table[1][1]

    print(f"  {len(results)} species processed, {len(skipped)} skipped")

    if not results:
        print("[ERROR] No results. Exiting.")
        sys.exit(1)

    df_res = pd.DataFrame(results)

    # FDR correction
    pvals = df_res["pval_fisher"].values
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")
    df_res["padj_fisher"] = padj
    df_res["significant"] = padj < 0.05
    df_res["enriched_in_orphans"] = (df_res["odds_ratio"] > 1) & (padj < 0.05)
    df_res["depleted_in_orphans"] = (df_res["odds_ratio"] < 1) & (padj < 0.05)

    df_res = df_res.sort_values("padj_fisher")
    df_res.to_csv(OUT_PER_SP, sep="\t", index=False, float_format="%.6g")
    print(f"[OUT] Per-species results → {OUT_PER_SP}")

    # Global Fisher test across all species
    global_tbl = [[global_table[0], global_table[1]],
                  [global_table[2], global_table[3]]]
    g_or, g_pval = fisher_exact(global_tbl, alternative="two-sided")

    total_orphans = global_table[0] + global_table[1]
    total_others  = global_table[2] + global_table[3]
    pct_orphan_te = 100.0 * global_table[0] / total_orphans if total_orphans > 0 else 0
    pct_other_te  = 100.0 * global_table[2] / total_others  if total_others  > 0 else 0

    df_global = pd.DataFrame([{
        "metric":            "global_fisher",
        "n_orphan_genes":    total_orphans,
        "n_other_genes":     total_others,
        "n_orphan_with_TE":  global_table[0],
        "n_other_with_TE":   global_table[2],
        "pct_orphan_TE":     round(pct_orphan_te, 3),
        "pct_other_TE":      round(pct_other_te, 3),
        "odds_ratio":        round(g_or, 4),
        "pval":              g_pval,
        "n_species":         len(results),
        "n_species_sig_enrich": int(df_res["enriched_in_orphans"].sum()),
        "n_species_sig_deplete": int(df_res["depleted_in_orphans"].sum()),
    }])
    df_global.to_csv(OUT_GLOBAL, sep="\t", index=False, float_format="%.6g")
    print(f"[OUT] Global summary → {OUT_GLOBAL}")
    print(f"      Global OR={g_or:.3f}, p={g_pval:.2e}")
    print(f"      Orphan genes with TE domain: {pct_orphan_te:.2f}%")
    print(f"      Other genes with TE domain:  {pct_other_te:.2f}%")
    print(f"      Sig enriched in orphans: {int(df_res['enriched_in_orphans'].sum())} species")
    print(f"      Sig depleted in orphans: {int(df_res['depleted_in_orphans'].sum())} species")

    # ── figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: odds ratio distribution across species
    ax = axes[0]
    colors = []
    for _, row in df_res.iterrows():
        if row["enriched_in_orphans"]:
            colors.append("#d62728")
        elif row["depleted_in_orphans"]:
            colors.append("#1f77b4")
        else:
            colors.append("#aec7e8")
    log_or = np.log2(df_res["odds_ratio"].clip(lower=1e-3, upper=1e3))
    ax.barh(range(len(df_res)), log_or, color=colors, height=0.7)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("log₂ Odds Ratio (TE domain: orphan vs non-orphan genes)")
    ax.set_title("Per-species log₂(OR)")
    ax.set_yticks([])
    red_p   = mpatches.Patch(color="#d62728", label="Enriched (FDR<0.05)")
    blue_p  = mpatches.Patch(color="#1f77b4", label="Depleted (FDR<0.05)")
    grey_p  = mpatches.Patch(color="#aec7e8", label="Not significant")
    ax.legend(handles=[red_p, blue_p, grey_p], fontsize=9)

    # Panel B: global % comparison bar
    ax2 = axes[1]
    cats = ["Orphan genes", "Non-orphan genes"]
    pcts = [pct_orphan_te, pct_other_te]
    bars = ax2.bar(cats, pcts, color=["#ff7f0e", "#1f77b4"], width=0.5)
    ax2.set_ylabel("% genes with ≥1 TE domain")
    ax2.set_title(
        f"Global: OR={g_or:.2f}, p={g_pval:.2e}\n"
        f"n_orphan={total_orphans:,}  n_other={total_others:,}  "
        f"n_species={len(results)}"
    )
    for bar, pct in zip(bars, pcts):
        ax2.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + 0.05, f"{pct:.2f}%",
                 ha="center", va="bottom", fontsize=11)
    ax2.set_ylim(0, max(pcts) * 1.3 + 0.5)

    plt.tight_layout()
    plt.savefig(OUT_FIG, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"[OUT] Figure → {OUT_FIG}")
    print("[DONE]")


if __name__ == "__main__":
    main()
