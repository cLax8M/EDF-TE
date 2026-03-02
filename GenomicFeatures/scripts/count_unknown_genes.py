#!/usr/bin/env python3

import os
from collections import defaultdict

IPR_DIR = "/Users/carloslaxmolina/Library/Mobile Documents/com~apple~CloudDocs/Lab_Postdoctorado/UniNE-Croll Lab/Holomycota_TE_Evolution/data/functional_interpro"
PROT_DIR = "/Users/carloslaxmolina/Library/Mobile Documents/com~apple~CloudDocs/Lab_Postdoctorado/UniNE-Croll Lab/Holomycota_TE_Evolution/data/proteins_annot_faa"
OUT_DIR  = "/Users/carloslaxmolina/Library/Mobile Documents/com~apple~CloudDocs/Lab_Postdoctorado/UniNE-Croll Lab/Holomycota_TE_Evolution/results/unknown_function"

os.makedirs(OUT_DIR, exist_ok=True)

# --------------------------------------------------
# 1. Leer proteomas
# --------------------------------------------------
proteins = defaultdict(set)

for fname in os.listdir(PROT_DIR):
    if not fname.endswith(".fasta") and not fname.endswith(".aa.fasta"):
        continue

    species = fname.split(".")[0]
    with open(os.path.join(PROT_DIR, fname)) as f:
        for line in f:
            if line.startswith(">"):
                pid = line[1:].split()[0]
                proteins[species].add(pid)

# --------------------------------------------------
# 2. Leer anotaciones InterPro
# --------------------------------------------------
annotations = defaultdict(lambda: defaultdict(list))

for fname in os.listdir(IPR_DIR):
    if not fname.endswith("_IPR.tab"):
        continue

    species = fname.replace("_IPR.tab", "")
    with open(os.path.join(IPR_DIR, fname)) as f:
        for line in f:
            if line.strip() == "":
                continue

            cols = line.rstrip("\n").split("\t")
            pid = cols[0]
            desc = cols[5] if len(cols) > 5 else ""
            go   = cols[11] if len(cols) > 11 else ""

            annotations[species][pid].append((desc, go))

# --------------------------------------------------
# 3. Clasificar unknown function genes
# --------------------------------------------------
results = []

for species in sorted(proteins):
    total = len(proteins[species])
    unknown = 0

    for pid in proteins[species]:
        annots = annotations.get(species, {}).get(pid, [])

        # Caso 1: sin ninguna anotación
        if not annots:
            unknown += 1
            continue

        # Caso 2: solo DUF y sin GO
        only_duf = True
        has_go = False

        for desc, go in annots:
            if go and go != "-":
                has_go = True
            if desc and "DUF" not in desc and "Domain of Unknown Function" not in desc:
                only_duf = False

        if only_duf and not has_go:
            unknown += 1

    annotated = total - unknown
    perc_unknown = (unknown / total * 100) if total > 0 else 0

    results.append((
        species,
        total,
        annotated,
        unknown,
        round(perc_unknown, 2)
    ))

# --------------------------------------------------
# 4. Guardar resultados
# --------------------------------------------------
header = (
    "Species\tTotal_proteins\tAnnotated_proteins\t"
    "Unknown_function_proteins\tPercent_unknown\n"
)

txt_out = os.path.join(OUT_DIR, "unknown_function_summary.txt")
tsv_out = os.path.join(OUT_DIR, "unknown_function_summary.tsv")

for out in [txt_out, tsv_out]:
    with open(out, "w") as f:
        f.write(header)
        for r in results:
            f.write("\t".join(map(str, r)) + "\n")

print("✔ Analysis finished correctly")
print(f"Results written to:\n - {txt_out}\n - {tsv_out}")