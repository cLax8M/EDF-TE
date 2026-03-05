#!/usr/bin/env python3
"""
02_annotate_results.py
----------------------
Collect treeWAS results across all 16 runs, annotate features with
descriptions and TE flags, and write a single consolidated table.

Usage:
  python 02_annotate_results.py [--resultsdir PATH] [--ipr_dir PATH]

Output:
  <resultsdir>/treewas_annotated_hits.tsv

Columns:
  feature, phenotype, matrix, n_tests_significant, tests_significant,
  mean_score, robust, description, is_TE_pfam
"""

import argparse
import os
import glob

import pandas as pd

# ── TE-related Pfam IDs ───────────────────────────────────────────────────────
TE_PFAM = {
    "PF00078",  # Reverse transcriptase (LTR retrotransposons)
    "PF03184",  # DDE transposase (class II TEs)
    "PF13456",  # Reverse transcriptase-like
    "PF14214",  # Helitron helicase
    "PF13683",  # Gypsy integrase
    "PF01497",  # Mariner transposase
    "PF13952",  # AAA+ATPase (TE mobilization)
}


def parse_args():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_base = "/legserv/Temp/Carlos/EDF-TE"
    p = argparse.ArgumentParser(description="Annotate treeWAS hits")
    p.add_argument("--resultsdir",
                   default=os.path.join(base, "results"),
                   help="Directory containing per-run subdirectories")
    p.add_argument("--ipr_dir",
                   default=f"{data_base}/data/functional_interpro",
                   help="Directory with *_IPR.tab files (for IPR descriptions)")
    p.add_argument("--orthogroups",
                   default=f"{data_base}/results/Orthofinder_run/Orthogroups/Orthogroups.GeneCount.tsv",
                   help="OrthoFinder GeneCount TSV (to list species per OG)")
    return p.parse_args()


# ── 1. COLLECT SUMMARY TSVs ──────────────────────────────────────────────────

def collect_results(resultsdir: str) -> pd.DataFrame:
    pattern = os.path.join(resultsdir, "*_OG", "summary.tsv")
    pattern2 = os.path.join(resultsdir, "*_IPR", "summary.tsv")
    files = glob.glob(pattern) + glob.glob(pattern2)

    if not files:
        raise FileNotFoundError(
            f"No summary.tsv files found under {resultsdir}. "
            "Run 01_submit_treewas.sh first."
        )

    dfs = []
    for f in sorted(files):
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            continue
        dfs.append(df)

    if not dfs:
        print("WARNING: All summary files are empty — no significant hits found.")
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)
    print(f"Collected {len(combined)} feature-phenotype associations from {len(files)} runs")
    robust_n = combined["robust"].sum() if "robust" in combined.columns else 0
    print(f"  Robust hits (≥2 tests): {robust_n}")
    return combined


# ── 2. BUILD IPR DESCRIPTION LOOKUP ─────────────────────────────────────────

def build_ipr_lookup(ipr_dir: str) -> dict:
    """Read a sample of IPR files to map iprId → iprDesc."""
    print("Building IPR description lookup …")
    lookup = {}
    files = sorted(glob.glob(os.path.join(ipr_dir, "*_IPR.tab")))
    # scan all files to maximise coverage of descriptions
    for fpath in files:
        try:
            df = pd.read_csv(fpath, sep="\t", header=None,
                             usecols=[1, 2],
                             names=["iprId", "iprDesc"])
            for _, row in df.dropna(subset=["iprId"]).iterrows():
                ipr_id = str(row["iprId"]).strip()
                if ipr_id and ipr_id not in lookup:
                    lookup[ipr_id] = str(row["iprDesc"]).strip()
        except Exception:
            continue
    print(f"  {len(lookup)} unique IPR terms with descriptions")
    return lookup


# ── 3. ANNOTATE ──────────────────────────────────────────────────────────────

def annotate(df: pd.DataFrame, ipr_lookup: dict) -> pd.DataFrame:
    if df.empty:
        return df

    descriptions = []
    is_te_list = []

    for _, row in df.iterrows():
        feature = str(row["feature"])
        matrix  = str(row.get("matrix", "")).upper()

        if matrix == "IPR":
            desc   = ipr_lookup.get(feature, "")
            # flag if any TE-related Pfam sub-ID maps to this IPR
            # InterPro IDs themselves don't carry PF prefix; flag based on known IPR↔Pfam
            # We also flag if the description mentions transposable element keywords
            te_keywords = ["transposase", "reverse transcriptase", "retroviral", "integrase",
                           "helitron", "mariner", "gypsy", "transposon", "transposable",
                           "AAA+ ATPase", "retrotransposon"]
            is_te = any(kw.lower() in desc.lower() for kw in te_keywords)
        elif matrix == "OG":
            desc   = ""  # OG descriptions not available from OrthoFinder directly
            is_te  = False
        else:
            desc  = ""
            is_te = False

        descriptions.append(desc)
        is_te_list.append(is_te)

    df = df.copy()
    df["description"] = descriptions
    df["is_TE_pfam"]  = is_te_list
    return df


# ── 4. FLAG OG HITS WITH TE PFAM OVERLAP ────────────────────────────────────

def flag_og_te_overlap(df: pd.DataFrame, ipr_dir: str, resultsdir: str) -> pd.DataFrame:
    """
    Check if OG hits correspond to features that overlap known TE Pfam domains.
    We do this by cross-referencing OG names with IPR hits from the IPR matrix:
    if an OG feature ID appears in the IPR results with a TE-related description,
    flag it.

    Simpler approach: check if any Pfam-mapped IPR (PF00078, PF03184, etc.) is
    known to appear in TE proteins. We can't do a perfect OG→Pfam mapping without
    functional annotations per OG, so we rely on the IPR matrix + description scan.
    For now, set is_TE_pfam=False for all OG hits and note that validation requires
    manual inspection.
    """
    return df


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    combined = collect_results(args.resultsdir)
    if combined.empty:
        print("No hits to annotate. Exiting.")
        return

    ipr_lookup = build_ipr_lookup(args.ipr_dir)
    annotated  = annotate(combined, ipr_lookup)
    annotated  = flag_og_te_overlap(annotated, args.ipr_dir, args.resultsdir)

    # sort: robust first, then by n_tests_significant desc, then by mean_score desc
    sort_cols = []
    if "robust" in annotated.columns:
        sort_cols.append(("robust", False))
    if "n_tests_significant" in annotated.columns:
        sort_cols.append(("n_tests_significant", False))
    if "mean_score" in annotated.columns:
        sort_cols.append(("mean_score", False))

    if sort_cols:
        annotated = annotated.sort_values(
            [c for c, _ in sort_cols],
            ascending=[asc for _, asc in sort_cols]
        )

    out_path = os.path.join(args.resultsdir, "treewas_annotated_hits.tsv")
    annotated.to_csv(out_path, sep="\t", index=False)
    print(f"\nAnnotated hits written to: {out_path}")
    print(f"  Total rows      : {len(annotated)}")
    if "robust" in annotated.columns:
        print(f"  Robust hits     : {annotated['robust'].sum()}")
    if "is_TE_pfam" in annotated.columns:
        print(f"  TE-flagged IPRs : {annotated['is_TE_pfam'].sum()}")

    # print top 20 robust hits
    if "robust" in annotated.columns:
        top = annotated[annotated["robust"]].head(20)
        if not top.empty:
            print("\nTop 20 robust hits:")
            display_cols = [c for c in ["feature", "phenotype", "matrix",
                                         "n_tests_significant", "mean_score",
                                         "description", "is_TE_pfam"]
                            if c in top.columns]
            print(top[display_cols].to_string(index=False))


if __name__ == "__main__":
    main()
