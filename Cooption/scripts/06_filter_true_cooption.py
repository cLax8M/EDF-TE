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

if candidates.empty:
    # Keep schema consistent even when there are no candidate genes.
    candidates = candidates.copy()
    true_coopted = candidates.copy()
    true_coopted.to_csv(output_file, sep="\t", index=False)
    print(f"[OK] True TE co-option genes written to: {output_file}")
    print(f"[INFO] Initial candidates: 0")
    print(f"[INFO] High-confidence co-opted genes: 0")
    sys.exit(0)

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
        # TE-derived domain identifiers — both free-text names and Pfam accessions
        # (Pfam IDs are stored by 05_classify_genes.py when matched via TE_CORE_IDS)
        if d in {
            # Legacy name matches
            "RVT_1", "rve", "RNase_H", "RNase_H2", "Integrase",
            "Reverse_transcriptase", "Retrotrans_gag", "zf-CCHC",
            "Chromo", "Chromo_shadow", "Protease", "DDE_1",
            "DDE_Tnp_IS240", "DDE_Tnp_4", "DDE_Tnp_Tc5", "DDE_3",
            "hAT_N", "MULE", "Tnp_P_element", "CENP-B_N",
            "Helitron_like_N", "Helicase_C", "PIF1",
            "Exo_endo_phos", "Gag_MA", "Gag_p24", "Gag_p30",
            "Peptidase_A17", "Peptidase_A2B", "TLV_coat",
            # Pfam accessions (stored by step 05 when matched by ID)
            "PF00078", "PF07727", "PF13846", "PF00075",
            "PF00665", "PF02022", "PF13683",
            "PF00540", "PF00607", "PF02093", "PF03732",
            "PF01184", "PF08284", "PF00207",
            "PF03184", "PF01609", "PF01359", "PF05699",
            "PF10551", "PF13359", "PF13843", "PF13358",
            "PF03221", "PF18107",
            "PF14214", "PF00271", "PF05970",
            "PF03372", "PF00098",
            "PF00385", "PF01772", "PF13714",
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

true_coopted.drop(columns=["domains_list", "true_cooption"], errors="ignore") \
            .to_csv(output_file, sep="\t", index=False)

print(f"[OK] True TE co-option genes written to: {output_file}")
print(f"[INFO] Initial candidates: {len(candidates)}")
print(f"[INFO] High-confidence co-opted genes: {len(true_coopted)}")