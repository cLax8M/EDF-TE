#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) != 4:
    sys.exit(
        "Usage: 03_split_TE_host_domains.py "
        "<protein_classes.tsv> <domains.tsv> <output.tsv>"
    )

protein_classes_file = sys.argv[1]
domains_file = sys.argv[2]
output_file = sys.argv[3]

# Load files
pc = pd.read_csv(protein_classes_file, sep="\t")
domains = pd.read_csv(domains_file, sep="\t")

# Required columns
pc_required = {
    "gene_id", "architecture",
    "n_TE_domains", "n_nonTE_domains"
}
dom_required = {
    "gene_id", "domain_name", "domain_id"
}

missing_pc = pc_required - set(pc.columns)
missing_dom = dom_required - set(domains.columns)

if missing_pc:
    raise ValueError(f"Missing columns in protein_classes: {missing_pc}")
if missing_dom:
    raise ValueError(f"Missing columns in domains file: {missing_dom}")

# Keep only TE-host fusion proteins
fusions = pc[pc["architecture"] == "TE_host_fusion"].copy()

# Merge domains
merged = fusions.merge(domains, on="gene_id", how="left")

# Assign domain type based on architecture counts
def classify_domain(row):
    if row["n_TE_domains"] > 0 and row["n_nonTE_domains"] > 0:
        return "mixed"
    elif row["n_TE_domains"] > 0:
        return "TE_domain"
    else:
        return "host_domain"

merged["domain_class"] = merged.apply(classify_domain, axis=1)

# Select output columns
out = merged[
    [
        "gene_id",
        "domain_id",
        "domain_name",
        "domain_class",
        "n_TE_domains",
        "n_nonTE_domains"
    ]
]

out.to_csv(output_file, sep="\t", index=False)
print(f"[OK] Written {output_file}")