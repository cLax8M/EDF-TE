#!/usr/bin/env python3
"""
prepare_figure_data.py
======================
Consolidates all per-species data tracks needed for the circular phylogenetic
figure:
  1. % unknown genes
  2. % paralog / multicopy genes
  3. % OG families above fungal average
  4. Genome size (Mb)
  5. Gene density (genes / Mb)
  6. Average introns per gene
  7. % real TEs (Class I + Class II, no Unclassified / Simple Repeats)
  8. % Class I TEs  (of genome)
  9. % Class II TEs (of genome)
 10. % genes with ≥2 isoforms
 11. Clade annotation (for arc coloring)

Inputs  (all paths relative to project roots):
  - TE + genomic features : D01_te_genomics_merged_provisional.tsv
  - Unknown function       : unknown_function_summary.tsv
  - Duplication            : DuplicationCG/orthogroup_stats.tsv
  - Isoforms               : GenomicFeatures_V3/{SP}_ALL_isoform_counts.txt
  - Tree tip labels        : guide_rooted_Naegrub.treefile  (for taxon list)

Output:
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figure_species_data.tsv

Usage:
  python3 prepare_figure_data.py
"""

import re
import os
import sys
import pandas as pd
from pathlib import Path

# ─── Paths ────────────────────────────────────────────────────────────────────
TREE        = Path("/data/carlos/EDF-TE/SpeciesTree/08_guide_tree/guide_rooted_Naegrub.treefile")
D01         = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/D01_te_genomics_merged_provisional.tsv")
UNKNOWN     = Path("/legserv/Temp/Carlos/EDF-TE/results/unknown_function/unknown_function_summary.tsv")
DUPLICATION = Path("/legserv/Temp/Carlos/EDF-TE/results/DuplicationCG/orthogroup_stats.tsv")
GENFEAT_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures_V3")
OUT         = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figure_species_data.tsv")

# ─── Clade mapping (2-3 letter prefix → group name for arc coloring) ─────────
CLADE_MAP = {
    "AP": "Aphelida",
    "AS": "Ascomycota",
    "BA": "Basidiomycota",
    "BL": "Blastocladiomycota",
    "BS": "Basidiobolales",          # Entomophthoromycotina
    "CH": "Chytridiomycota",
    "CR": "Cryptomycota",
    "EN": "Entomophthoromycota",
    "GL": "Glomeromycota",
    "HO": "Holozoa",                 # outgroup
    "KX": "Kickxellomycotina",
    "MI": "Microsporidia",
    "MO": "Mortierellomycota",
    "MU": "Mucoromycota",
    "NU": "Nucleariida",
    "OL": "Olpidiomycota",
    "ZO": "Zoopagomycota",
}

# ─── 1. Read taxon list from tree ─────────────────────────────────────────────
tree_text = TREE.read_text()
# Tree tips look like: AP_Amoocci.aa  (with optional .aa suffix)
raw_tips = sorted(set(re.findall(r'[A-Z]{2,3}_[A-Za-z0-9_]+(?:\.aa)?', tree_text)))
# Normalise: strip .aa suffix for matching against data files
tips_clean = [t.replace(".aa", "") for t in raw_tips]
# Keep the .aa version for the tree label
tip_map = {t.replace(".aa", ""): t for t in raw_tips}   # clean → tree_label

print(f"[1] Tree tips: {len(tips_clean)}")

# ─── 2. D01 base table (TE + genomic features already merged) ─────────────────
d01 = pd.read_csv(D01, sep="\t")
d01 = d01.rename(columns={"species": "taxon"})

# Fix known naming discrepancies between D01 and tree tip labels
d01["taxon"] = d01["taxon"].replace({
    "AS_AspfumUD":  "AS_AspfumUD1",
    "NU_Parvatl":   "NU_Parvatl_c",
})

# Real TE % = Class I + Class II (exclude Unclassified and Other)
# Class I  = LTR + LINE + SINE + PLE (Penelope)
# Class II = DNA + RC  (Rolling Circle)
d01["pct_ClassI"]  = d01[["pct_LTR","pct_LINE","pct_SINE","pct_PLE"]].sum(axis=1)
d01["pct_ClassII"] = d01[["pct_DNA","pct_RC"]].sum(axis=1)
d01["pct_realTE"]  = d01["pct_ClassI"] + d01["pct_ClassII"]
d01["genome_size_mb"] = d01["Genome_Size_bp"] / 1e6
d01["gene_density"]   = d01["Total_Genes"] / d01["genome_size_mb"]

print(f"[2] D01 rows: {len(d01)}")

# ─── 3. Unknown function ──────────────────────────────────────────────────────
unk = pd.read_csv(UNKNOWN, sep="\t")[["Species","Percent_unknown"]]
unk = unk.rename(columns={"Species":"taxon","Percent_unknown":"pct_unknown"})
print(f"[3] Unknown function rows: {len(unk)}")

# ─── 4. Duplication / multicopy ───────────────────────────────────────────────
dup = pd.read_csv(DUPLICATION, sep="\t")[[
    "Species","Pct_paralog_genes","Pct_OGs_above_fungal_avg"
]]
dup = dup.rename(columns={
    "Species":                  "taxon",
    "Pct_paralog_genes":        "pct_paralog_genes",
    "Pct_OGs_above_fungal_avg": "pct_OGs_above_fungal_avg",
})
print(f"[4] Duplication rows: {len(dup)}")

# ─── 5. % genes with ≥2 isoforms (from isoform_counts files) ─────────────────
# If ALL genes in a file have Isoform_count = 1, the annotation is a single-
# isoform-per-gene gene-prediction (JGI ab-initio with no RNA-seq isoforms).
# These are indistinguishable from "no data" and are reported as NaN.
# Only files that contain at least one gene with ≥2 isoforms are treated as
# real transcript-diversity data.
isoform_records = []
n_real = 0
n_placeholder = 0
for sp in tips_clean:
    iso_file = GENFEAT_DIR / f"{sp}_ALL_isoform_counts.txt"
    if not iso_file.exists():
        # Some species have only PROTEIN_CODING version
        iso_file = GENFEAT_DIR / f"{sp}_PROTEIN_CODING_isoform_counts.txt"
    if iso_file.exists():
        try:
            df_iso = pd.read_csv(iso_file, sep="\t")
            total   = len(df_iso)
            multi   = int((df_iso["Isoform_count"] >= 2).sum())
            max_iso = int(df_iso["Isoform_count"].max()) if total > 0 else 1
            # If no gene has >1 isoform, the annotation lacks splicing evidence
            if max_iso <= 1:
                isoform_records.append({"taxon": sp, "pct_multiisoform": float("nan")})
                n_placeholder += 1
            else:
                pct_multi = 100.0 * multi / total if total > 0 else 0.0
                isoform_records.append({"taxon": sp, "pct_multiisoform": pct_multi})
                n_real += 1
        except Exception as e:
            print(f"  WARN isoform {sp}: {e}", file=sys.stderr)
            isoform_records.append({"taxon": sp, "pct_multiisoform": float("nan")})
    else:
        isoform_records.append({"taxon": sp, "pct_multiisoform": float("nan")})

iso_df = pd.DataFrame(isoform_records)
print(f"[5] Isoform data: real={n_real}, placeholder_NaN={n_placeholder}, missing={iso_df['pct_multiisoform'].isna().sum() - n_placeholder}")

# ─── 6. Merge all ─────────────────────────────────────────────────────────────
df = pd.DataFrame({"taxon": tips_clean})
df = df.merge(d01[[
    "taxon","pct_realTE","pct_ClassI","pct_ClassII","pct_TE_total",
    "pct_Unknown","genome_size_mb","gene_density",
    "Avg_Introns_All","Total_Genes","Avg_Isoforms_Per_Gene",
]].rename(columns={"pct_Unknown":"pct_TE_unknown"}),
    on="taxon", how="left")
df = df.merge(unk,    on="taxon", how="left")
df = df.merge(dup,    on="taxon", how="left")
df = df.merge(iso_df, on="taxon", how="left")

# ─── 7. Add clade annotation ──────────────────────────────────────────────────
df["prefix"] = df["taxon"].str.extract(r'^([A-Z]{2,3})_')
df["clade"]  = df["prefix"].map(CLADE_MAP).fillna("Unknown")

# Add tree label (with .aa suffix) for joining in R
df["tree_label"] = df["taxon"].map(tip_map)

# ─── 8. Round and save ────────────────────────────────────────────────────────
float_cols = [c for c in df.columns if df[c].dtype == float]
df[float_cols] = df[float_cols].round(4)

OUT.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT, sep="\t", index=False)

print(f"\n[OK] Written {len(df)} rows × {len(df.columns)} columns → {OUT}")
print(f"Columns: {df.columns.tolist()}")

# Quick coverage report
print("\nCoverage (non-NA counts per track):")
key_cols = ["pct_realTE","pct_unknown","pct_paralog_genes",
            "pct_OGs_above_fungal_avg","genome_size_mb",
            "gene_density","Avg_Introns_All","pct_multiisoform"]
for c in key_cols:
    if c in df.columns:
        n = df[c].notna().sum()
        print(f"  {c:35s}: {n}/{len(df)}")
