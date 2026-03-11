#!/usr/bin/env python3
"""
Annotate treeWAS robust genome-size OGs with IPR domains.

For each robust OG from log10_Genome_Size_OG, looks up all member proteins
in the per-species IPR annotation files and identifies TE-related functions.

Usage:
    python annotate_og_ipr.py
"""

import os
import csv
import sys
from collections import defaultdict

# ─── Paths ───────────────────────────────────────────────────────────────────
SUMMARY_TSV    = "/data/carlos/EDF-TE/EDF-TE/TreeWAS/results/log10_Genome_Size_OG/summary.tsv"
ORTHOGROUPS_TSV = "/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups/Orthogroups.tsv"
IPR_DIR        = "/legserv/Temp/Carlos/EDF-TE/data/functional_interpro"
OUTPUT_TSV     = "/data/carlos/EDF-TE/EDF-TE/TreeWAS/results/log10_Genome_Size_OG/og_ipr_annotation.tsv"

# ─── TE-related keywords ─────────────────────────────────────────────────────
TE_KEYWORDS = [
    "transposase", "reverse transcriptase", "retroviral", "integrase",
    "helitron", "mariner", "gypsy", "transposon", "transposable",
    "aaa+ atpase", "aaa atpase", "aaa-atpase", "retrotransposon",
    "retrovirus", "cchc", "zinc knuckle", "rnase h", "mobile element",
    "dna transposon", "ty1", "ty3", "copia", "chromovirus",
    "zinc finger, cchc", "znknc", "gag", "pol polyprotein",
]

def is_te_related(ipr_desc: str) -> bool:
    desc_lower = ipr_desc.lower()
    return any(kw in desc_lower for kw in TE_KEYWORDS)

def parse_protein_id(prot_id: str) -> str:
    """Extract IPR lookup key from a protein ID string."""
    prot_id = prot_id.strip()
    if prot_id.startswith("jgi|"):
        parts = prot_id.split("|")
        if len(parts) >= 3:
            return parts[2]   # numeric JGI protein ID
    return prot_id            # NCBI accession, used as-is

# ─── Step 1: Load robust OGs ─────────────────────────────────────────────────
print("Loading robust OGs ...", file=sys.stderr)
robust_ogs = set()
robust_meta = {}
with open(SUMMARY_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["robust"].strip() == "TRUE":
            og = row["feature"].strip()
            robust_ogs.add(og)
            robust_meta[og] = {
                "n_tests":  row["n_tests_significant"].strip(),
                "tests":    row["tests_significant"].strip(),
                "score":    row["mean_score"].strip(),
            }
print(f"  {len(robust_ogs)} robust OGs found", file=sys.stderr)

# ─── Step 2: Build species → IPR file mapping ────────────────────────────────
# IPR filenames: PREFIX_SPECIESCODE_IPR.tab → key = PREFIX_SPECIESCODE
ipr_files = {}
for fname in os.listdir(IPR_DIR):
    if fname.endswith("_IPR.tab"):
        key = fname[: -len("_IPR.tab")]
        ipr_files[key] = os.path.join(IPR_DIR, fname)
print(f"  {len(ipr_files)} IPR files indexed", file=sys.stderr)

# ─── Step 3: Map Orthogroups.tsv column names → IPR keys ─────────────────────
# Column names look like "AP_Amoocci.aa" or "BA_Hetan2" (no .aa suffix)
with open(ORTHOGROUPS_TSV) as f:
    header_cols = f.readline().rstrip("\n").split("\t")
species_cols = header_cols[1:]  # drop "Orthogroup"

col_to_ipr_key = {}
for col in species_cols:
    # Strip common FASTA-extension suffixes
    for suffix in (".aa", ".faa", ".prot", ".fa", ".fasta"):
        if col.endswith(suffix):
            key = col[: -len(suffix)]
            if key in ipr_files:
                col_to_ipr_key[col] = key
                break
    else:
        if col in ipr_files:
            col_to_ipr_key[col] = col

mapped   = len(col_to_ipr_key)
unmapped  = [c for c in species_cols if c not in col_to_ipr_key]
print(f"  {mapped}/{len(species_cols)} species columns mapped to IPR files", file=sys.stderr)
if unmapped:
    print(f"  Unmapped: {unmapped[:10]}", file=sys.stderr)

# ─── Step 4: Parse Orthogroups.tsv — collect proteins for robust OGs ─────────
print("Scanning Orthogroups.tsv for robust OG members ...", file=sys.stderr)
og_by_species = defaultdict(lambda: defaultdict(list))  # og → species_col → [prot_ids]

col_idx = {col: i + 1 for i, col in enumerate(species_cols)}

with open(ORTHOGROUPS_TSV) as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for row in reader:
        og_id = row[0]
        if og_id not in robust_ogs:
            continue
        for sp_col, idx in col_idx.items():
            if idx >= len(row) or not row[idx].strip():
                continue
            prots = [p.strip() for p in row[idx].split(",") if p.strip()]
            if prots:
                og_by_species[og_id][sp_col].extend(prots)

print(f"  Proteins collected for {len(og_by_species)} OGs", file=sys.stderr)

# ─── Step 5: Load IPR data for needed species ————————————————————————————
print("Loading IPR annotation files ...", file=sys.stderr)
needed_species = set()
for sp_dict in og_by_species.values():
    needed_species.update(sp_dict.keys())

species_ipr = {}  # species_col → {lookup_key → set[(ipr_id, ipr_desc)]}
for sp_col in sorted(needed_species):
    if sp_col not in col_to_ipr_key:
        continue
    ipr_path = ipr_files[col_to_ipr_key[sp_col]]
    prot_annots = defaultdict(set)
    with open(ipr_path, errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            pid     = parts[0].strip()
            ipr_id  = parts[1].strip()
            ipr_desc = parts[2].strip()
            if ipr_id and ipr_id != "-":
                prot_annots[pid].add((ipr_id, ipr_desc))
    species_ipr[sp_col] = prot_annots

print(f"  IPR data loaded for {len(species_ipr)} species", file=sys.stderr)

# ─── Step 6: Aggregate IPR counts per OG ─────────────────────────────────────
print("Aggregating IPR terms per robust OG ...", file=sys.stderr)
og_summary = {}  # og → {total_prots, annotated_prots, ipr: {id → {desc, count, te}}}

for og_id in robust_ogs:
    ipr_counts = defaultdict(lambda: {"desc": "", "count": 0, "te": False})
    total_prots = 0
    annotated_prots = 0

    for sp_col, prots in og_by_species.get(og_id, {}).items():
        sp_ipr = species_ipr.get(sp_col, {})
        for prot in prots:
            total_prots += 1
            key  = parse_protein_id(prot)
            hits = sp_ipr.get(key, set())
            if hits:
                annotated_prots += 1
            for ipr_id, ipr_desc in hits:
                ipr_counts[ipr_id]["desc"]  = ipr_desc
                ipr_counts[ipr_id]["count"] += 1
                if is_te_related(ipr_desc):
                    ipr_counts[ipr_id]["te"] = True

    og_summary[og_id] = {
        "total_prots":     total_prots,
        "annotated_prots": annotated_prots,
        "ipr":             dict(ipr_counts),
    }

# ─── Step 7: Write full output table ─────────────────────────────────────────
print(f"Writing {OUTPUT_TSV} ...", file=sys.stderr)
os.makedirs(os.path.dirname(OUTPUT_TSV), exist_ok=True)

with open(OUTPUT_TSV, "w") as out:
    out.write("og_id\tn_tests\ttests_significant\tmean_score\t"
              "total_prots\tannotated_prots\tipr_id\tipr_desc\tprotein_count\tte_related\n")
    for og_id in sorted(og_summary.keys()):
        data = og_summary[og_id]
        meta = robust_meta.get(og_id, {})
        ipr  = data["ipr"]
        if not ipr:
            out.write(f"{og_id}\t{meta.get('n_tests','')}\t{meta.get('tests','')}\t"
                      f"{meta.get('score','')}\t{data['total_prots']}\t"
                      f"{data['annotated_prots']}\tNA\tNA\t0\tFALSE\n")
        else:
            for ipr_id, info in sorted(ipr.items(), key=lambda x: -x[1]["count"]):
                te = "TRUE" if info["te"] else "FALSE"
                out.write(f"{og_id}\t{meta.get('n_tests','')}\t{meta.get('tests','')}\t"
                          f"{meta.get('score','')}\t{data['total_prots']}\t"
                          f"{data['annotated_prots']}\t{ipr_id}\t{info['desc']}\t"
                          f"{info['count']}\t{te}\n")

print("Done.", file=sys.stderr)

# ─── Step 8: Print TE summary to stdout ──────────────────────────────────────
te_ogs = {
    og: [(k, v) for k, v in d["ipr"].items() if v["te"]]
    for og, d in og_summary.items()
}
te_ogs = {og: hits for og, hits in te_ogs.items() if hits}

print(f"\n{'='*60}")
print(f"TE-RELATED IPR HITS IN ROBUST log10_Genome_Size OGs")
print(f"{'='*60}")
if not te_ogs:
    print("No TE-related IPR domains found in any robust genome-size OG.")
else:
    for og_id in sorted(te_ogs.keys()):
        meta = robust_meta[og_id]
        hits = te_ogs[og_id]
        data = og_summary[og_id]
        print(f"\n{og_id}  (tests: {meta['tests']}, score: {meta['score']}, "
              f"{data['total_prots']} proteins):")
        for ipr_id, info in sorted(hits, key=lambda x: -x[1]["count"]):
            print(f"  {ipr_id}  {info['desc']}  (n={info['count']})")

# ─── Print top-5 IPR per OG (for context) ────────────────────────────────────
print(f"\n{'='*60}")
print("TOP-3 IPR DOMAINS PER ROBUST OG (by protein count)")
print(f"{'='*60}")
for og_id in sorted(og_summary.keys()):
    data = og_summary[og_id]
    ipr  = data["ipr"]
    meta = robust_meta[og_id]
    top3 = sorted(ipr.items(), key=lambda x: -x[1]["count"])[:3]
    te_flag = "** TE **" if og_id in te_ogs else ""
    print(f"\n{og_id}  {te_flag}  [n_prots={data['total_prots']}, "
          f"annotated={data['annotated_prots']}, tests={meta['tests']}]")
    for ipr_id, info in top3:
        te_marker = " [TE]" if info["te"] else ""
        print(f"  {ipr_id}: {info['desc']} (n={info['count']}){te_marker}")
