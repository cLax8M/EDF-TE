#!/usr/bin/env python3

import pandas as pd
import sys
import os

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: python 01_annotate_domains.py <interpro.tsv> <TE_domains_master.tsv> <output.tsv>\n"
    )
    sys.exit(1)

ipr_file = sys.argv[1]
te_domains_file = sys.argv[2]
output_file = sys.argv[3]

# --------------------------------------------------
# Load InterProScan TSV (robust to column number)
# --------------------------------------------------
ipr = pd.read_csv(ipr_file, sep="\t", header=None)

n_cols = ipr.shape[1]
if n_cols < 10:
    raise ValueError(
        f"InterProScan file has {n_cols} columns, expected at least 10"
    )

# Standard InterProScan core columns (first 10)
ipr.columns = [
    "protein_id",
    "sequence_md5",
    "seq_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "end",
    "score",
    "status",
] + [f"extra_{i}" for i in range(n_cols - 10)]

# Keep only required fields
ipr = ipr[[
    "protein_id",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "end",
]]

ipr.rename(
    columns={
        "analysis": "domain_db",
        "signature_accession": "domain_id",
        "signature_description": "domain_name",
    },
    inplace=True,
)

# --------------------------------------------------
# Map InterProScan analyses to Pfam / CDD
# --------------------------------------------------
ipr["domain_db"] = ipr["domain_db"].str.strip()

analysis_to_db = {
    "HMMPfam": "Pfam",
    "Pfam-A": "Pfam",
    "pfam": "Pfam",
    "CDD": "CDD",
    "cd": "CDD",
    "cdd": "CDD"
}

ipr["domain_db"] = ipr["domain_db"].map(analysis_to_db)

# Keep only Pfam and CDD
ipr = ipr[ipr["domain_db"].isin(["Pfam", "CDD"])]

# --------------------------------------------------
# Load TE domain master table
# --------------------------------------------------
te = pd.read_csv(te_domains_file, sep="\t")

# Normalise fields
te["domain_db"] = te["domain_db"].str.strip()
te["domain_id"] = te["domain_id"].astype(str).str.strip()

# --------------------------------------------------
# Annotate domains with TE status
# --------------------------------------------------
annot = ipr.merge(
    te,
    on=["domain_db", "domain_id"],
    how="left",
)

# Domains not in TE list → non-TE
annot["TE_status"] = annot["TE_status"].fillna("non_TE")
annot["TE_superfamily"] = annot["TE_superfamily"].fillna("non_TE")

# --------------------------------------------------
# Sort and export
# --------------------------------------------------
annot = annot.sort_values(
    by=["protein_id", "start", "end"]
)

os.makedirs(os.path.dirname(output_file), exist_ok=True)

annot.to_csv(output_file, sep="\t", index=False)

print(f"[OK] Annotated domains written to: {output_file}")
print(f"[INFO] Total domain hits: {len(annot)}")
print(
    "[INFO] TE-associated hits:",
    (annot["TE_status"] != "non_TE").sum(),
)