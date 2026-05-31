#!/usr/bin/env python3
"""
C01_kimura_metrics.py
Compute Kimura divergence metrics per species from EarlGrey
divergence_summary_table.tsv files.

Metrics per species:
  - recency_index (RI): fraction of TE bp with mean_div < 0.1  (recent activity)
  - burst_score:        fraction of total TE bp in the dominant superfamily
  - mean_div_weighted:  genome-wide weighted mean Kimura divergence
  - n_superfamilies:    number of TE superfamilies detected

Metrics per species × superfamily (separate output):
  - RI per superfamily
  - mean_div per superfamily (weighted)

Usage:
    python C01_kimura_metrics.py <earlgrey_annotated_dir> \
        <output_species_tsv> <output_superfamily_tsv>
"""

import sys
import os
import math
import pandas as pd
from collections import defaultdict


# Non-TE repeat classes within subclass="Other" that must be excluded:
# structural RNAs, simple/satellite repeats, artefacts, chromosomal features.
NON_TE_SF = frozenset({
    "Simple_repeat", "Low_complexity", "Satellite",
    "rRNA", "tRNA", "snRNA", "scRNA", "sno", "srpRNA",
    "ARTEFACT", "acro", "centr", "subtelo", "DNA_virus", "RNA", "Y",
})


def load_divergence(divergence_file):
    """Load *_divergence_summary_table.tsv, return DataFrame of true TE rows."""
    df = pd.read_csv(divergence_file, sep="\t")
    # Drop rows with no divergence data
    df = df.dropna(subset=["mean_div", "total_bp"])
    df["mean_div"] = pd.to_numeric(df["mean_div"], errors="coerce")
    df["total_bp"] = pd.to_numeric(df["total_bp"], errors="coerce")
    df = df.dropna(subset=["mean_div", "total_bp"])
    # Exclude non-TE repeats (structural RNAs, simple/satellite repeats, etc.)
    # that RepeatMasker classifies under subclass="Other"
    non_te_mask = (df["subclass"] == "Other") & (df["superfamily"].isin(NON_TE_SF))
    df = df[~non_te_mask].copy()
    return df


def species_metrics(df):
    total_bp = df["total_bp"].sum()
    if total_bp == 0:
        return None

    # Recency index: bp fraction with mean_div < 0.1
    recent_bp = df.loc[df["mean_div"] < 0.1, "total_bp"].sum()
    ri = recent_bp / total_bp

    # Burst score: dominant superfamily bp fraction
    sf_bp = df.groupby("superfamily")["total_bp"].sum()
    burst = sf_bp.max() / total_bp
    dominant_sf = sf_bp.idxmax()

    # Weighted mean divergence
    mean_div_w = (df["mean_div"] * df["total_bp"]).sum() / total_bp

    return {
        "recency_index": round(ri, 4),
        "burst_score": round(burst, 4),
        "dominant_superfamily": dominant_sf,
        "mean_div_weighted": round(mean_div_w, 4),
        "n_superfamilies": df["superfamily"].nunique(),
        "total_te_bp": int(total_bp),
    }


def superfamily_metrics(df, sp):
    rows = []
    total_bp = df["total_bp"].sum()
    for sf, grp in df.groupby("superfamily"):
        sf_total = grp["total_bp"].sum()
        recent = grp.loc[grp["mean_div"] < 0.1, "total_bp"].sum()
        ri = recent / sf_total if sf_total > 0 else float("nan")
        mean_div_w = (grp["mean_div"] * grp["total_bp"]).sum() / sf_total if sf_total > 0 else float("nan")
        rows.append({
            "species": sp,
            "superfamily": sf,
            "total_bp": int(sf_total),
            "pct_of_TE": round(sf_total / total_bp * 100, 3) if total_bp > 0 else 0,
            "recency_index": round(ri, 4) if not math.isnan(ri) else "NA",
            "mean_div_weighted": round(mean_div_w, 4) if not math.isnan(mean_div_w) else "NA",
            "n_families": grp["family"].nunique(),
        })
    return rows


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: C01_kimura_metrics.py <annotated_dir> <species_tsv> <superfamily_tsv>")

    annotated_dir = sys.argv[1]
    out_species = sys.argv[2]
    out_sf = sys.argv[3]

    sp_records = []
    sf_records = []
    missing = []

    for entry in sorted(os.listdir(annotated_dir)):
        sp_dir = os.path.join(annotated_dir, entry)
        if not os.path.isdir(sp_dir) or not entry.endswith("_EarlGrey"):
            continue
        sp = entry.replace("_EarlGrey", "")
        summary_dir = os.path.join(sp_dir, f"{sp}_summaryFiles")
        div_file = os.path.join(summary_dir, f"{sp}_divergence_summary_table.tsv")
        if not os.path.isfile(div_file):
            missing.append(sp)
            continue

        df = load_divergence(div_file)
        if df.empty:
            missing.append(sp)
            continue

        m = species_metrics(df)
        if m:
            m["species"] = sp
            sp_records.append(m)
            sf_records.extend(superfamily_metrics(df, sp))

    os.makedirs(os.path.dirname(out_species) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(out_sf) or ".", exist_ok=True)

    pd.DataFrame(sp_records).to_csv(out_species, sep="\t", index=False)
    pd.DataFrame(sf_records).to_csv(out_sf, sep="\t", index=False)

    print(f"[OK] Species Kimura metrics: {out_species} ({len(sp_records)} species)")
    print(f"[OK] Superfamily Kimura metrics: {out_sf} ({len(sf_records)} rows)")
    if missing:
        print(f"[WARN] divergence_summary_table not found for {len(missing)} species")


if __name__ == "__main__":
    main()
