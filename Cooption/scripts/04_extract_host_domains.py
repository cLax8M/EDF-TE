#!/usr/bin/env python3

import sys
import pandas as pd

fusion_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(fusion_file, sep="\t")

# Accept both historical and current fusion schemas.
# Historical schema: protein_id + host_domain_names
# Current schema: gene_id + domain_name (+ optional domain_class)
id_col = None
for cand in ["protein_id", "gene_id"]:
    if cand in df.columns:
        id_col = cand
        break

if id_col is None:
    raise ValueError("Input fusion file must include either 'protein_id' or 'gene_id'.")

rows = []

if "host_domain_names" in df.columns:
    for _, r in df.iterrows():
        host_raw = "" if pd.isna(r["host_domain_names"]) else str(r["host_domain_names"]).strip()
        if host_raw == "":
            rows.append({"protein_id": r[id_col], "host_domain": "NONE"})
        else:
            for d in host_raw.split(","):
                dom = d.strip()
                if dom:
                    rows.append({"protein_id": r[id_col], "host_domain": dom})
else:
    # Build host-domain list from per-domain rows.
    # If domain_class exists, avoid TE-only rows when possible.
    use_df = df.copy()
    if "domain_class" in use_df.columns:
        keep = use_df["domain_class"].astype(str).str.lower().isin(["host_domain", "mixed"])
        if keep.any():
            use_df = use_df.loc[keep].copy()

    if "domain_name" in use_df.columns:
        grouped = use_df.groupby(id_col, dropna=False)["domain_name"]
        for pid, vals in grouped:
            uniq = sorted({str(v).strip() for v in vals if pd.notna(v) and str(v).strip()})
            if not uniq:
                rows.append({"protein_id": pid, "host_domain": "NONE"})
            else:
                for dom in uniq:
                    rows.append({"protein_id": pid, "host_domain": dom})
    else:
        for pid in sorted(set(use_df[id_col].dropna().astype(str))):
            rows.append({"protein_id": pid, "host_domain": "NONE"})

out = pd.DataFrame(rows, columns=["protein_id", "host_domain"])
out.to_csv(output_file, sep="\t", index=False)

print(f"[OK] Host domains extracted: {output_file}")
print(f"[INFO] Proteins with host domain annotations: "
      f"{(out.host_domain != 'NONE').sum()}")