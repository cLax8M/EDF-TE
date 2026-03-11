#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: 00_normalize_interpro.py <input_IPR.tab> <output.domains.tsv>")
    sys.exit(1)

inp, out = sys.argv[1], sys.argv[2]

# --------------------------------------------------
# Inspect first line to infer format
# --------------------------------------------------
with open(inp) as f:
    first = f.readline()

# --------------------------------------------------
# JGI-like format
# --------------------------------------------------
if "domainStarts" in first and "domainEnds" in first:
    df = pd.read_csv(
        inp,
        sep="\t",
        comment="#",
        header=None,
        engine="python",
        on_bad_lines="skip"
    )

    df_out = pd.DataFrame({
        "gene_id": df.iloc[:, 0],
        "domain_id": df.iloc[:, 4],
        "domain_name": df.iloc[:, 5].fillna("NA"),
        "source_db": df.iloc[:, 3],
        "start": df.iloc[:, 7].astype(str).str.rstrip(","),
        "end": df.iloc[:, 8].astype(str).str.rstrip(",")
    })

# --------------------------------------------------
# InterProScan / Galaxy format
# --------------------------------------------------
else:
    df = pd.read_csv(
        inp,
        sep="\t",
        header=None,
        engine="python",
        on_bad_lines="skip"
    )

    if df.shape[1] < 8:
        raise ValueError(f"Unrecognized InterPro format in {inp}")

    df_out = pd.DataFrame({
        "gene_id": df.iloc[:, 0],
        "domain_id": df.iloc[:, 4],
        "domain_name": df.iloc[:, 5].fillna("NA"),
        "source_db": df.iloc[:, 3],
        "start": df.iloc[:, 6],
        "end": df.iloc[:, 7]
    })

# --------------------------------------------------
# Clean & export
# --------------------------------------------------
df_out = df_out.dropna(subset=["domain_id"])
df_out.to_csv(out, sep="\t", index=False)

print(f"[OK] Normalized domains written to: {out}")
print(f"[INFO] Total domain hits: {len(df_out)}")