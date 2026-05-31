#!/usr/bin/env python3

import sys
import pandas as pd
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: step3_classify_genes.py <domains.tsv> <output.tsv>")
    sys.exit(1)

domains_file = sys.argv[1]
output_file = sys.argv[2]

# Domain name heuristics (legacy, for backward compatibility)
TE_CORE_NAMES = {
    "RT", "RVT", "RNaseH", "RNase_H",
    "Integrase", "INT", "Gag", "Protease"
}

# Pfam accessions for TE core machinery (Oggenfuss 2024 + literature)
# When matched by ID, the accession string is stored as the domain label.
TE_CORE_IDS = {
    # Reverse transcriptases (LTR retrotransposons, LINEs)
    "PF00078",  # RVT_1 — RNA-dependent DNA polymerase (RT)
    "PF07727",  # RVT_2 — RT thumb domain
    "PF13846",  # RVT_3 — RT-like
    "PF00075",  # RNase_H — ribonuclease H
    # Integrases (LTR retrotransposons)
    "PF00665",  # rve — integrase core domain
    "PF02022",  # Integrase_H2C — integrase zinc-binding
    "PF13683",  # rve_3 — integrase C-terminal
    # Gag proteins (LTR)
    "PF00540",  # Gag_MA — matrix domain
    "PF00607",  # Gag_p24 — capsid domain
    "PF02093",  # Gag_p30
    "PF03732",  # Retrotrans_gag — Copia/Gypsy gag
    # Retroviral proteases
    "PF01184",  # Peptidase_A17 — retroviral protease
    "PF08284",  # Peptidase_A2B — retroviral protease-like
    "PF00207",  # TLV_coat — retroviral coat protein
    # TIR / DDE transposases (DNA transposons)
    "PF03184",  # DDE_1 — Tc1/Mariner transposase
    "PF01609",  # DDE_Tnp_IS240 — Tc1/Mariner (SETMAR domain)
    "PF01359",  # Transposase_1 — IS3 family
    "PF05699",  # hAT_N — hAT transposase N-terminal
    "PF10551",  # MULE — Mutator-like transposase
    "PF13359",  # DDE_Tnp_4 — DDE transposase (IS family)
    "PF13843",  # DDE_Tnp_Tc5 — Tc5B transposase
    "PF13358",  # DDE_3 — DDE transposase superfamily
    "PF03221",  # Tnp_P_element — P-element transposase
    "PF18107",  # CENP-B_N — pogo-like / CENP-B transposase
    # Helitron
    "PF14214",  # Helitron_like_N — Helitron N-terminal
    "PF00271",  # Helicase_C — superfamily-II helicase (Helitron)
    "PF05970",  # PIF1 — helicase associated with Helitrons
    # LINE-specific
    "PF03372",  # Exo_endo_phos — apurinic endonuclease (LINE ORF2)
    "PF00098",  # zf-CCHH (zf-CCHC) — zinc finger in LINEs/LTRs
}

# Domain name heuristics (legacy)
TE_ASSOC_NAMES = {
    "Chromo", "Chromo_shadow"
}

# Pfam accessions for TE-associated chromatin readers
TE_ASSOC_IDS = {
    "PF00385",  # Chromo — chromodomain (Gypsy-associated; chromatin reader)
    "PF01772",  # Chromo_shadow
    "PF13714",  # Chromo_2
}

genes = defaultdict(lambda: {
    "te_core": set(),
    "te_assoc": set(),
    "non_te": set()
})

df = pd.read_csv(domains_file, sep="\t", dtype=str).fillna("")

if "gene_id" not in df.columns:
    raise ValueError(f"Missing required column 'gene_id' in {domains_file}")
if "domain_name" not in df.columns:
    raise ValueError(f"Missing required column 'domain_name' in {domains_file}")

for _, r in df.iterrows():
    gene_id = str(r["gene_id"]).strip()
    domain_name = str(r["domain_name"]).strip()
    domain_id   = str(r["domain_id"]).strip() if "domain_id" in df.columns else ""
    domain_class = str(r["domain_class"]).strip() if "domain_class" in df.columns else ""

    # Prefer whichever field looks like a Pfam accession as the display label.
    # In JGI-legacy files the columns are swapped: domain_id = DB name ("HMMPfam")
    # and domain_name = accession ("PF00078"), so check both.
    _looks_pfam = lambda s: bool(s and s.startswith("PF") and s[2:].isdigit())
    if _looks_pfam(domain_name):
        label = domain_name
    elif domain_id and domain_id not in ("", "nan"):
        label = domain_id
    else:
        label = domain_name

    # Match against both columns so JGI-swapped files are handled correctly
    if (domain_name in TE_CORE_NAMES or domain_id in TE_CORE_IDS
            or domain_name in TE_CORE_IDS or domain_class == "TE_core"):
        genes[gene_id]["te_core"].add(label)
    elif (domain_name in TE_ASSOC_NAMES or domain_id in TE_ASSOC_IDS
            or domain_name in TE_ASSOC_IDS or domain_class == "TE_associated"):
        genes[gene_id]["te_assoc"].add(label)
    else:
        genes[gene_id]["non_te"].add(label)

with open(output_file, "w") as out:
    out.write("gene_id\tclass\tn_TE_core\tn_TE_assoc\tn_non_TE\tdomains\n")

    for gene, d in sorted(genes.items()):
        # --- Context-dependent rule: RNase_H (PF00075) ---
        # PF00075 is a cellular enzyme that also appears in LINE/LTR TEs.
        # It only counts as TE_core when co-occurring with a genuine RT domain
        # (PF00078 RVT_1). Alone it is treated as a host/non-TE domain to avoid
        # false positives in species with many cellular RNase H genes (e.g.
        # some Glomeromycota where PF00075 with ProSite/SUPERFAMILY redundant
        # hits accounts for the majority of spurious co-option calls).
        RT_DOMAINS = {"PF00078", "PF07727", "PF13846"}   # RVT_1/2/3
        if "PF00075" in d["te_core"] and not d["te_core"] & RT_DOMAINS:
            d["te_core"].discard("PF00075")
            d["non_te"].add("PF00075")

        n_core = len(d["te_core"])
        n_assoc = len(d["te_assoc"])
        n_non = len(d["non_te"])

        if n_core >= 2:
            gene_class = "TE_complete"
        elif n_core == 1 and n_non == 0:
            gene_class = "TE_fragment"
        elif n_core >= 1 and n_non >= 1:
            gene_class = "TE_coopted_candidate"
        elif n_core == 0 and n_assoc >= 1:
            gene_class = "TE_coopted_candidate"
        else:
            gene_class = "non_TE"

        all_domains = ",".join(
            sorted(d["te_core"] | d["te_assoc"] | d["non_te"])
        )

        out.write(
            f"{gene}\t{gene_class}\t{n_core}\t{n_assoc}\t{n_non}\t{all_domains}\n"
        )

print(f"[OK] Classification written to {output_file}")