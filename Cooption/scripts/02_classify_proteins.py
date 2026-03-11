#!/usr/bin/env python3
import sys
import pandas as pd

TE_KEYWORDS = [
    "retrotrans", "integrase", "reverse transcriptase",
    "gag", "pol", "endonuclease", "ribonuclease h",
    "helicase", "protease", "transposase"
]

def is_TE_domain(name):
    name = str(name).lower()
    return any(k in name for k in TE_KEYWORDS)

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile, sep="\t")

records = []

for gene, sub in df.groupby("gene_id"):
    te = sub["domain_name"].apply(is_TE_domain)
    n_te = te.sum()
    n_non = (~te).sum()

    if n_te > 0 and n_non > 0:
        arch = "TE_host_fusion"
    elif n_te > 0:
        arch = "TE_only"
    elif n_non > 0:
        arch = "host_only"
    else:
        arch = "no_domains"

    records.append([
        gene, n_te > 0, n_non > 0,
        n_te, n_non, arch
    ])

out = pd.DataFrame(records, columns=[
    "gene_id",
    "has_TE_domain",
    "has_nonTE_domain",
    "n_TE_domains",
    "n_nonTE_domains",
    "architecture"
])

out.to_csv(outfile, sep="\t", index=False)
print(f"[OK] Written {outfile}")