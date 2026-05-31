#!/usr/bin/env python3
"""
E04_outlier_genes.py — Orthogroup enrichment in TE-landscape outliers
=====================================================================
Compares gene content (OrthoFinder orthogroups) of clade outliers
(from E01_clade_outliers.py) vs non-outlier species within each clade.

For each clade with ≥2 outlier and ≥2 non-outlier species:
  1. Binarise orthogroup presence/absence from Orthogroups.GeneCount.tsv
  2. Fisher's exact test per OG (enriched/depleted in outliers)
  3. BH-FDR correction
  4. Annotate significant OGs with InterPro domain descriptions

Outputs
-------
  <outdir>/E04_per_clade/        — per-clade TSVs of significant OGs
  <outdir>/E04_combined.tsv      — all clades combined
  <outdir>/E04_summary.tsv       — per-clade summary counts

Usage
-----
  python E04_outlier_genes.py [options]

  --outliers     clade_outliers.tsv from E01
  --orthogroups  Orthogroups.GeneCount.tsv
  --og_proteins  Orthogroups.tsv  (for OG→protein annotation)
  --ipr_dir      directory with *_IPR.tab files
  --outdir       output directory
  --alpha        FDR significance threshold (default: 0.1)
  --min_outliers minimum outlier species per clade (default: 2)

Dependencies: pandas, numpy, scipy.
"""

import argparse
import os
import sys
import glob
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# ── paths ─────────────────────────────────────────────────────────────────────
_DATA = "/legserv/Temp/Carlos/EDF-TE"
_SEC2 = f"{_DATA}/results/Section2_TElandscape"


def parse_args():
    p = argparse.ArgumentParser(
        description="Orthogroup enrichment in TE-landscape outlier species")
    p.add_argument("--outliers",
                   default=f"{_SEC2}/clade_outliers/clade_outliers.tsv")
    p.add_argument("--orthogroups",
                   default=f"{_DATA}/results/Orthofinder_run/Orthogroups/Orthogroups.GeneCount.tsv")
    p.add_argument("--og_proteins",
                   default=f"{_DATA}/results/Orthofinder_run/Orthogroups/Orthogroups.tsv")
    p.add_argument("--ipr_dir",
                   default=f"{_DATA}/data/functional_interpro")
    p.add_argument("--outdir", default=f"{_SEC2}/E04_outlier_genes")
    p.add_argument("--alpha", type=float, default=0.1,
                   help="FDR significance threshold")
    p.add_argument("--min_outliers", type=int, default=2,
                   help="Minimum outlier species per clade to test")
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════════

def load_outliers(path: str) -> pd.DataFrame:
    """Load E01 outlier table.  Returns df with 'species', 'clade',
    'outlier_univariate', 'outlier_composition'."""
    df = pd.read_csv(path, sep="\t")
    # Mark as outlier if flagged by either method
    df["is_outlier"] = df["outlier_univariate"] | df["outlier_composition"]
    return df


def load_og_binary(path: str, species_set: set) -> pd.DataFrame:
    """Load Orthogroups.GeneCount.tsv, binarise to presence/absence.

    Returns DataFrame: index=OrthoGroup, columns=species (no .aa suffix).
    Reads in chunks to manage memory on large files.
    """
    print(f"[LOAD] Reading OrthoGroup gene counts (chunked) …")
    chunks = []
    reader = pd.read_csv(path, sep="\t", index_col=0, chunksize=5000)
    for chunk in reader:
        # strip .aa suffix from columns
        chunk.columns = [c.replace(".aa", "").replace(".faa", "")
                         for c in chunk.columns]
        cols_keep = [c for c in chunk.columns if c in species_set]
        if not cols_keep:
            continue
        sub = chunk[cols_keep]
        # binarise
        sub = (sub >= 1).astype("int8")
        # drop OGs absent in all species
        sub = sub[sub.sum(axis=1) > 0]
        if not sub.empty:
            chunks.append(sub)

    if not chunks:
        sys.exit("[ERROR] No orthogroups found for the given species set.")

    og = pd.concat(chunks)
    # remove duplicates (shouldn't happen, but safe)
    og = og[~og.index.duplicated(keep="first")]
    print(f"[LOAD] OG binary matrix: {og.shape[0]} OGs × {og.shape[1]} species")
    return og


# ═══════════════════════════════════════════════════════════════════════════════
# PER-CLADE ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════════

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    fdr = np.minimum(pvals * n / ranks, 1.0)
    # enforce monotonicity
    fdr_sorted = fdr[order].copy()
    running_min = fdr_sorted[-1]
    for i in range(n - 2, -1, -1):
        running_min = min(running_min, fdr_sorted[i])
        fdr_sorted[i] = running_min
    fdr[order] = fdr_sorted
    return fdr


def test_clade(clade: str, outlier_species: list[str],
               non_outlier_species: list[str],
               og_binary: pd.DataFrame) -> pd.DataFrame:
    """Fisher's exact test for each OG: enriched/depleted in outliers."""
    n_out = len(outlier_species)
    n_non = len(non_outlier_species)

    # select columns present in OG matrix
    out_cols = [s for s in outlier_species if s in og_binary.columns]
    non_cols = [s for s in non_outlier_species if s in og_binary.columns]
    if len(out_cols) < 2 or len(non_cols) < 2:
        return pd.DataFrame()

    out_mat = og_binary[out_cols]
    non_mat = og_binary[non_cols]

    # per-OG counts
    out_present = out_mat.sum(axis=1)
    non_present = non_mat.sum(axis=1)
    out_absent  = len(out_cols) - out_present
    non_absent  = len(non_cols) - non_present

    # filter: skip OGs present in ALL or NONE
    total_present = out_present + non_present
    total = len(out_cols) + len(non_cols)
    variable = (total_present > 0) & (total_present < total)
    ogs = og_binary.index[variable]

    records = []
    for og in ogs:
        a = int(out_present[og])   # outlier + OG present
        b = int(out_absent[og])    # outlier + OG absent
        c = int(non_present[og])   # non-outlier + OG present
        d = int(non_absent[og])    # non-outlier + OG absent

        # Fisher's exact test (two-sided)
        _, p = fisher_exact([[a, b], [c, d]])

        # odds ratio (with Haldane correction for zeros)
        OR = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

        records.append({
            "orthogroup": og,
            "clade": clade,
            "n_outlier": len(out_cols),
            "n_non_outlier": len(non_cols),
            "outlier_present": a,
            "outlier_absent": b,
            "nonoutlier_present": c,
            "nonoutlier_absent": d,
            "odds_ratio": round(OR, 4),
            "p_fisher": p,
        })

    if not records:
        return pd.DataFrame()

    result = pd.DataFrame(records)
    result["p_fdr"] = bh_fdr(result["p_fisher"].values)
    result["direction"] = np.where(result["odds_ratio"] > 1, "enriched_in_outliers",
                                   "depleted_in_outliers")
    result.sort_values("p_fisher", inplace=True)
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# ANNOTATION (InterPro domains for significant OGs)
# ═══════════════════════════════════════════════════════════════════════════════

def annotate_ogs(results: pd.DataFrame, og_proteins_path: str,
                 ipr_dir: str) -> pd.DataFrame:
    """Add InterPro domain descriptions to significant OGs.

    Strategy: for each significant OG, look up which proteins belong to it
    (from Orthogroups.tsv), then find their InterPro hits.
    Returns the results DataFrame with added 'top_ipr_domains' column.
    """
    sig_ogs = set(results["orthogroup"].unique())
    if not sig_ogs:
        results["top_ipr_domains"] = ""
        return results

    # ── read OG → protein mapping (only for significant OGs) ─────────────
    og_proteins = {}
    if os.path.isfile(og_proteins_path):
        print(f"[ANNOT] Reading protein lists for {len(sig_ogs)} significant OGs …")
        with open(og_proteins_path) as f:
            header = f.readline().strip().split("\t")
            species_cols = header[1:]  # species names
            for line in f:
                parts = line.strip().split("\t")
                og = parts[0]
                if og not in sig_ogs:
                    continue
                proteins = set()
                for cell in parts[1:]:
                    if cell.strip():
                        for prot in cell.split(","):
                            prot = prot.strip()
                            if prot:
                                proteins.add(prot)
                og_proteins[og] = proteins
    else:
        print(f"[WARN] Orthogroups.tsv not found at {og_proteins_path}. Skipping annotation.")
        results["top_ipr_domains"] = ""
        return results

    # ── build protein → IPR description from IPR files ────────────────────
    # Only load proteins that are in significant OGs
    all_prots = set()
    for prots in og_proteins.values():
        all_prots |= prots

    print(f"[ANNOT] Loading InterPro annotations for {len(all_prots)} proteins …")
    prot_to_iprs = {}
    ipr_files = glob.glob(os.path.join(ipr_dir, "*_IPR.tab"))
    for fpath in ipr_files:
        try:
            df = pd.read_csv(fpath, sep="\t", usecols=[0, 1, 2],
                             names=["proteinId", "iprId", "iprDesc"],
                             header=0, dtype=str)
            for _, row in df.iterrows():
                prot = str(row["proteinId"])
                if prot in all_prots and pd.notna(row.get("iprId")):
                    desc = str(row["iprDesc"]) if pd.notna(row.get("iprDesc")) else row["iprId"]
                    prot_to_iprs.setdefault(prot, []).append(desc)
        except Exception:
            continue

    # ── aggregate: OG → most common InterPro descriptions ─────────────────
    og_annot = {}
    for og, prots in og_proteins.items():
        domain_counts = Counter()
        for prot in prots:
            for desc in prot_to_iprs.get(prot, []):
                domain_counts[desc] += 1
        if domain_counts:
            top3 = domain_counts.most_common(3)
            og_annot[og] = "; ".join(f"{d} ({n})" for d, n in top3)
        else:
            og_annot[og] = ""

    results["top_ipr_domains"] = results["orthogroup"].map(og_annot).fillna("")
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    per_clade_dir = os.path.join(args.outdir, "E04_per_clade")
    os.makedirs(per_clade_dir, exist_ok=True)

    # ── load outlier data ─────────────────────────────────────────────────
    outliers = load_outliers(args.outliers)
    all_species = set(outliers["species"])
    print(f"[DATA] {len(outliers)} species, {outliers['is_outlier'].sum()} outliers across {outliers['clade'].nunique()} clades")

    # ── load OG binary matrix ─────────────────────────────────────────────
    og_binary = load_og_binary(args.orthogroups, all_species)

    # ── test each clade ───────────────────────────────────────────────────
    all_results = []
    summary_rows = []

    for clade, grp in outliers.groupby("clade"):
        out_sp  = grp.loc[grp["is_outlier"], "species"].tolist()
        non_sp  = grp.loc[~grp["is_outlier"], "species"].tolist()

        if len(out_sp) < args.min_outliers:
            continue
        if len(non_sp) < 2:
            continue

        print(f"\n[CLADE] {clade}: {len(out_sp)} outliers vs {len(non_sp)} non-outliers")

        result = test_clade(clade, out_sp, non_sp, og_binary)
        if result.empty:
            print(f"  No variable OGs to test.")
            continue

        n_tested = len(result)
        sig = result[result["p_fdr"] < args.alpha]
        n_sig = len(sig)
        n_enriched = (sig["direction"] == "enriched_in_outliers").sum()
        n_depleted = (sig["direction"] == "depleted_in_outliers").sum()

        print(f"  Tested: {n_tested} OGs | Significant (FDR<{args.alpha}): {n_sig} "
              f"(↑{n_enriched} enriched, ↓{n_depleted} depleted)")

        # save per-clade
        clade_path = os.path.join(per_clade_dir, f"{clade}_enrichment.tsv")
        sig.to_csv(clade_path, sep="\t", index=False)

        all_results.append(sig)
        summary_rows.append({
            "clade": clade,
            "n_outlier": len(out_sp),
            "n_non_outlier": len(non_sp),
            "n_tested": n_tested,
            "n_significant": n_sig,
            "n_enriched": n_enriched,
            "n_depleted": n_depleted,
            "outlier_species": ", ".join(sorted(out_sp)),
        })

    # ── combine results ───────────────────────────────────────────────────
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        print(f"\n[TOTAL] {len(combined)} significant OGs across all clades")

        # ── annotate ──────────────────────────────────────────────────────
        combined = annotate_ogs(combined, args.og_proteins, args.ipr_dir)

        combined_path = os.path.join(args.outdir, "E04_combined.tsv")
        combined.to_csv(combined_path, sep="\t", index=False)
        print(f"[SAVED] {combined_path}")
    else:
        print("\n[TOTAL] No significant enrichments found at FDR < {:.2f}".format(args.alpha))

    # ── summary ───────────────────────────────────────────────────────────
    if summary_rows:
        summary = pd.DataFrame(summary_rows)
        summary_path = os.path.join(args.outdir, "E04_summary.tsv")
        summary.to_csv(summary_path, sep="\t", index=False)
        print(f"[SAVED] {summary_path}")
        print("\n" + summary.to_string(index=False))

    print("\n[DONE]")


if __name__ == "__main__":
    main()
