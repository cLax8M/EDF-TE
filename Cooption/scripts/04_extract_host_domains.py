#!/usr/bin/env python3

import sys
import pandas as pd

fusion_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(fusion_file, sep="\t")

# Expand host domains
rows = []
for _, r in df.iterrows():
    if pd.isna(r["host_domain_names"]) or r["host_domain_names"] == "":
        rows.append({
            "protein_id": r["protein_id"],
            "host_domain": "NONE"
        })
    else:
        for d in r["host_domain_names"].split(","):
            rows.append({
                "protein_id": r["protein_id"],
                "host_domain": d.strip()
            })

out = pd.DataFrame(rows)
out.to_csv(output_file, sep="\t", index=False)

print(f"[OK] Host domains extracted: {output_file}")
print(f"[INFO] Proteins with host domain annotations: "
      f"{(out.host_domain != 'NONE').sum()}")