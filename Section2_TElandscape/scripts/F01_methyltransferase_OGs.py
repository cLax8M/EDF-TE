#!/usr/bin/env python3
"""
F01_methyltransferase_OGs.py
=============================
Map C-5 cytosine methyltransferases (IPR001525/IPR018117) to OrthoFinder
orthogroups.  For each OG containing at least one methyltransferase protein,
report which species are present/absent.

Usage:
    python F01_methyltransferase_OGs.py
"""

import os, re, glob, sys
from collections import defaultdict

IPR_DIR   = "/legserv/Temp/Carlos/EDF-TE/data/functional_interpro"
OG_FILE   = "/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups/Orthogroups.tsv"
TARGET_IPRS = {"IPR001525", "IPR018117"}

# ── Step 1: collect methyltransferase protein IDs per species ────────────
print("Step 1: Collecting C-5 methyltransferase proteins from IPR files...", flush=True)
species_meth = defaultdict(set)          # species -> set of protein IDs (str)

for fpath in sorted(glob.glob(os.path.join(IPR_DIR, "*_IPR.tab"))):
    sp = os.path.basename(fpath).replace("_IPR.tab", "")
    with open(fpath) as fh:
        hdr = fh.readline().strip().split("\t")
        hdr = [c.lstrip("#") for c in hdr]
        if "iprId" not in hdr or "proteinId" not in hdr:
            continue
        ii = hdr.index("iprId")
        pi = hdr.index("proteinId")
        for line in fh:
            fs = line.split("\t")
            if len(fs) > max(ii, pi) and fs[ii] in TARGET_IPRS:
                species_meth[sp].add(fs[pi])

total = sum(len(v) for v in species_meth.values())
print(f"  {len(species_meth)} species, {total} proteins", flush=True)

# reverse map: (species_short, protein_id) -> full_species
# species_meth keys are like "MU_Phybl4", fasta headers are "jgi|Phybl4|ID|..."
# strip clade prefix: "MU_Phybl4" -> "Phybl4"
sp_pid_lookup = set()    # (species_short, protein_id)
full_sp_map = {}         # species_short -> full species name
for sp, pids in species_meth.items():
    sp_short = sp.split("_", 1)[1] if "_" in sp else sp
    full_sp_map[sp_short] = sp
    for p in pids:
        sp_pid_lookup.add((sp_short, p))

# ── Step 2: scan Orthogroups.tsv ─────────────────────────────────────────
print("Step 2: Scanning Orthogroups.tsv ...", flush=True)
og_hits = defaultdict(lambda: defaultdict(set))   # OG -> species -> prot IDs

with open(OG_FILE) as fh:
    fh.readline()                                  # skip header
    for ln, line in enumerate(fh, 2):
        if ln % 30000 == 0:
            print(f"  ... {ln:,} OGs", flush=True)
        for m in re.finditer(r'jgi\|(\w+)\|(\d+)\|', line):
            sp_short, pid = m.group(1), m.group(2)
            if (sp_short, pid) in sp_pid_lookup:
                og = line[:line.index('\t')]
                full_sp = full_sp_map[sp_short]
                og_hits[og][full_sp].add(pid)

print(f"\n  {len(og_hits)} OGs contain C-5 methyltransferases\n", flush=True)

# ── Step 3: report ───────────────────────────────────────────────────────
all_mu = sorted(s for s in species_meth if s.startswith("MU_"))

print("Top 20 OGs by species count:")
for og in sorted(og_hits, key=lambda x: len(og_hits[x]), reverse=True)[:20]:
    n_sp = len(og_hits[og])
    n_mu = sum(1 for s in og_hits[og] if s.startswith("MU_"))
    n_p  = sum(len(p) for p in og_hits[og].values())
    print(f"  {og}: {n_sp} sp ({n_mu} MU), {n_p} prots")

print("\nNeurospora crassa (AS_Ncrass2) OGs:")
for og in sorted(og_hits):
    if "AS_Ncrass2" in og_hits[og]:
        n_sp = len(og_hits[og])
        n_mu = sum(1 for s in og_hits[og] if s.startswith("MU_"))
        pids = og_hits[og]["AS_Ncrass2"]
        print(f"  {og}: Ncrass2 prots={pids} | {n_sp} sp ({n_mu} MU)")

print("\nMU presence per large OG:")
for og in sorted(og_hits, key=lambda x: len(og_hits[x]), reverse=True)[:10]:
    mu_in  = sorted(s for s in og_hits[og] if s.startswith("MU_"))
    mu_out = [s for s in all_mu if s not in og_hits[og]]
    print(f"\n  {og}: {len(mu_in)}/{len(all_mu)} MU present")
    if mu_out and len(mu_out) <= 15:
        print(f"    absent: {', '.join(mu_out)}")
    elif mu_out:
        print(f"    absent: {len(mu_out)} MU species")
    if mu_in and len(mu_in) <= 15:
        print(f"    present: {', '.join(mu_in)}")
