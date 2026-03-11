#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: 06_filter_true_cooption.py <gene_classes.tsv> <output.tsv>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep="\t")

# --------------------------------------------------
# Filtrar solo candidatos
# --------------------------------------------------

candidates = df[df["class"] == "TE_coopted_candidate"].copy()

def parse_domains(s):
    if pd.isna(s) or s.strip() == "":
        return []
    return [d.strip() for d in s.split(",")]

candidates["domains_list"] = candidates["domains"].apply(parse_domains)

# --------------------------------------------------
# Criterio de co-option real
# --------------------------------------------------
# Aceptar si:
# - tiene al menos un dominio non-TE
# - no todos los dominios son TE-like

def is_true_cooption(domains):
    if not domains:
        return False

    has_non_te = False
    all_te = True

    for d in domains:
        # heurística simple: dominios TE conocidos
        if d in {
            "RVT_1", "rve", "RNase_H", "Integrase",
            "Reverse_transcriptase", "Retrotrans_gag",
            "zf-CCHC", "Chromo", "Protease", "DDE_1"
        }:
            continue
        else:
            has_non_te = True
            all_te = False

    return has_non_te and not all_te

candidates["true_cooption"] = candidates["domains_list"].apply(is_true_cooption)

true_coopted = candidates[candidates["true_cooption"]].copy()

# --------------------------------------------------
# Output
# --------------------------------------------------

true_coopted.drop(columns=["domains_list", "true_cooption"]) \
            .to_csv(output_file, sep="\t", index=False)

print(f"[OK] True TE co-option genes written to: {output_file}")
print(f"[INFO] Initial candidates: {len(candidates)}")
print(f"[INFO] High-confidence co-opted genes: {len(true_coopted)}")