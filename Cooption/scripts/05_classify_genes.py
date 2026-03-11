#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: step3_classify_genes.py <domains.tsv> <output.tsv>")
    sys.exit(1)

domains_file = sys.argv[1]
output_file = sys.argv[2]

TE_CORE = {
    "RT", "RVT", "RNaseH", "RNase_H",
    "Integrase", "INT", "Gag", "Protease"
}

TE_ASSOC = {
    "Chromo", "Chromo_shadow"
}

genes = defaultdict(lambda: {
    "te_core": set(),
    "te_assoc": set(),
    "non_te": set()
})

with open(domains_file) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip().split("\t")
        gene_id = fields[0]
        domain_name = fields[3]
        domain_class = fields[7]

        if domain_name in TE_CORE or domain_class == "TE_core":
            genes[gene_id]["te_core"].add(domain_name)
        elif domain_name in TE_ASSOC or domain_class == "TE_associated":
            genes[gene_id]["te_assoc"].add(domain_name)
        else:
            genes[gene_id]["non_te"].add(domain_name)

with open(output_file, "w") as out:
    out.write("gene_id\tclass\tn_TE_core\tn_TE_assoc\tn_non_TE\tdomains\n")

    for gene, d in sorted(genes.items()):
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