#!/usr/bin/env python3
"""Rename sequences to HoloTE_nonfungal_TE_{i}_{tag}, preserving classification."""
import re
from Bio import SeqIO

WORK  = "/data/carlos/EDF-TE/TE_analysis/Library_curation"
IN_FA = f"{WORK}/10_tagged/consensus_tagged.fa"
OUT   = f"{WORK}/HoloTE_nonfungal_v1.fa"

i = 1
with open(OUT, "w") as out:
    for rec in SeqIO.parse(IN_FA, "fasta"):
        # Extract classification (#LTR/Gypsy etc) from original header
        cls_match = re.search(r'(#\S+)', rec.description)
        cls = cls_match.group(1) if cls_match else "#Unknown"
        # Extract evidence tag (_PE, _DA, _NE) from id
        tag_match = re.search(r'(_PE|_DA|_NE)$', rec.id)
        tag = tag_match.group(1) if tag_match else "_NE"
        new_id = f"HoloTE_nonfungal_TE_{i}{tag}{cls}"
        out.write(f">{new_id}\n{rec.seq}\n")
        i += 1

print(f"Final library: {OUT}")
print(f"Total sequences: {i-1}")
