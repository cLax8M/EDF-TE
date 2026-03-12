#!/usr/bin/env python3
"""
Two analyses:

1. Sirtuin (OG0007438) copy number per species vs. Pct_OGs_multicopy,
   Pct_paralog_genes, and TE% — testing whether sirtuin loss correlates
   with more duplication (TE-silencing hypothesis).

2. TE content (total_masked_pct, LTR%, LINE%, DNA%) vs. all 8 phenotypes
   — Spearman correlations across the 181 species with EarlGrey data.

Outputs:
  TreeWAS/data/sirtuin_per_species.tsv    — OG0007438 copy number per species
  TreeWAS/data/te_phenotype_correlations.tsv — Spearman r + p-value table
"""

import os
import csv
import sys
from collections import defaultdict
from scipy import stats
import numpy as np

# ─── Paths ───────────────────────────────────────────────────────────────────
BASE          = "/data/carlos/EDF-TE/EDF-TE/TreeWAS"
ORTHOGROUPS   = "/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups/Orthogroups.tsv"
PHENOTYPES    = os.path.join(BASE, "results/phenotypes.tsv")
TE_CONTENT    = os.path.join(BASE, "data/te_content_per_species.tsv")
OUT_SIRTUIN   = os.path.join(BASE, "data/sirtuin_per_species.tsv")
OUT_CORR      = os.path.join(BASE, "data/te_phenotype_correlations.tsv")

SIRTUIN_OG    = "OG0007438"
EDF_PREFIXES  = {"BL","BS","CH","CR","EN","GL","KX","MI","MO","MU","ZO"}

TE_COLS   = ["total_masked_pct", "bases_masked_pct",
             "LTR_pct", "LINE_pct", "SINE_pct", "DNA_pct",
             "Helitron_pct", "Unclassified_pct", "genome_Mb"]
PHENO_COLS = ["log10_Genome_Size", "Pct_Intergenic", "Avg_Introns_Multiexon",
              "Avg_Exons_Per_Gene", "Pct_paralog_genes", "Pct_OGs_multicopy",
              "Percent_unknown", "mean_r_intra"]


# ═══════════════════════════════════════════════════════════════════════════
# PART 1 — Sirtuin copy number per species
# ═══════════════════════════════════════════════════════════════════════════
print("="*60)
print("PART 1: Sirtuin (OG0007438) copy number per species")
print("="*60)

sirtuin_counts = {}   # species_code → n_proteins

with open(ORTHOGROUPS) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    sp_cols = header[1:]
    for row in reader:
        if row[0] != SIRTUIN_OG:
            continue
        for i, col in enumerate(sp_cols):
            # col looks like "BL_Allma1.aa" or "BL_Allma1"
            sp = col.split(".")[0]
            prefix = sp.split("_")[0] if "_" in sp else sp[:2].upper()
            if prefix not in EDF_PREFIXES:
                continue
            val = row[i+1].strip() if i+1 < len(row) else ""
            prots = [p.strip() for p in val.split(",") if p.strip()]
            sirtuin_counts[sp] = len(prots)
        break

print(f"  {SIRTUIN_OG} found in {sum(1 for v in sirtuin_counts.values() if v>0)} EDF species")

# Load phenotype data
pheno = {}
with open(PHENOTYPES) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sp = row["species"].strip()
        prefix = sp.split("_")[0] if "_" in sp else sp[:2].upper()
        if prefix not in EDF_PREFIXES:
            continue
        pheno[sp] = {k: float(row[k]) for k in PHENO_COLS if row.get(k, "").strip() not in ("", "NA", "nan")}

# Load TE content
te = {}
with open(TE_CONTENT) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sp = row["species"].strip()
        te[sp] = {k: float(row[k]) for k in TE_COLS if row.get(k,"").strip() not in ("","NA","nan")}

# Merge: species that have ALL three
common = set(sirtuin_counts) & set(pheno) & set(te)
print(f"  Species with sirtuin + phenotype + TE data: {len(common)}")

# Write sirtuin table
with open(OUT_SIRTUIN, "w") as out:
    out.write("species\tprefix\tsirtuin_copies\tPct_OGs_multicopy\tPct_paralog_genes\ttotal_masked_pct\tLTR_pct\n")
    for sp in sorted(common):
        prefix = sp.split("_")[0]
        n = sirtuin_counts.get(sp, 0)
        mc  = pheno[sp].get("Pct_OGs_multicopy", float("nan"))
        par = pheno[sp].get("Pct_paralog_genes", float("nan"))
        tot = te[sp].get("total_masked_pct", float("nan"))
        ltr = te[sp].get("LTR_pct", float("nan"))
        out.write(f"{sp}\t{prefix}\t{n}\t{mc:.4f}\t{par:.4f}\t{tot:.2f}\t{ltr:.2f}\n")

# Summary by clade
print("\nSirtuin copies per clade  (mean ± range):")
by_clade = defaultdict(list)
for sp in sorted(sirtuin_counts):
    prefix = sp.split("_")[0]
    if prefix in EDF_PREFIXES:
        by_clade[prefix].append(sirtuin_counts[sp])

print(f"{'Clade':<8} {'n':>4}  {'mean':>6}  {'min':>4}  {'max':>4}  {'0-copy':>6}")
print("-" * 40)
for cl in sorted(by_clade):
    vals = by_clade[cl]
    n0 = sum(1 for v in vals if v == 0)
    print(f"{cl:<8} {len(vals):>4}  {np.mean(vals):>6.2f}  {min(vals):>4}  {max(vals):>4}  {n0:>6}")

# Correlations: sirtuin copies vs. key phenotypes + TE%
print("\nSpearman correlations: sirtuin copies vs. metrics:")
sirtuin_vec = []
mc_vec, par_vec, te_vec, ltr_vec, gs_vec = [], [], [], [], []
for sp in sorted(common):
    n = sirtuin_counts.get(sp, 0)
    if n == 0 and sirtuin_counts.get(sp) == 0:
        pass  # keep zeros
    sirtuin_vec.append(n)
    mc_vec.append(pheno[sp].get("Pct_OGs_multicopy", float("nan")))
    par_vec.append(pheno[sp].get("Pct_paralog_genes", float("nan")))
    te_vec.append(te[sp].get("total_masked_pct", float("nan")))
    ltr_vec.append(te[sp].get("LTR_pct", float("nan")))
    gs_vec.append(pheno[sp].get("log10_Genome_Size", float("nan")))

def spearman_clean(a, b):
    pairs = [(x, y) for x, y in zip(a, b) if not (np.isnan(x) or np.isnan(y))]
    if len(pairs) < 5:
        return float("nan"), float("nan"), 0
    xs, ys = zip(*pairs)
    r, p = stats.spearmanr(xs, ys)
    return r, p, len(pairs)

for label, vec in [("Pct_OGs_multicopy", mc_vec),
                   ("Pct_paralog_genes", par_vec),
                   ("total_masked_pct (TE%)", te_vec),
                   ("LTR_pct", ltr_vec),
                   ("log10_Genome_Size", gs_vec)]:
    r, p, n = spearman_clean(sirtuin_vec, vec)
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
    print(f"  {label:<30}  r={r:+.3f}  p={p:.3e}  n={n}  {sig}")


# ═══════════════════════════════════════════════════════════════════════════
# PART 2 — TE content vs. all phenotypes (Spearman)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("PART 2: TE content × phenotype correlations")
print("="*60)

# Build merged by species
te_pheno_rows = []
for sp in sorted(set(pheno) & set(te)):
    row = {"species": sp, "prefix": sp.split("_")[0]}
    row.update({k: te[sp].get(k, float("nan")) for k in TE_COLS})
    row.update({k: pheno[sp].get(k, float("nan")) for k in PHENO_COLS})
    te_pheno_rows.append(row)

n_merged = len(te_pheno_rows)
print(f"  {n_merged} species with both TE and phenotype data\n")

results = []
print(f"{'TE metric':<25} {'Phenotype':<28} {'r':>7} {'p':>10} {'n':>5}")
print("-" * 75)
for te_col in TE_COLS:
    for ph_col in PHENO_COLS:
        a = [r[te_col] for r in te_pheno_rows]
        b = [r[ph_col] for r in te_pheno_rows]
        r, p, n = spearman_clean(a, b)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        results.append((te_col, ph_col, r, p, n, sig))
        if abs(r) >= 0.20 or p < 0.05:   # print only meaningful ones
            print(f"  {te_col:<23} {ph_col:<28} {r:+.3f} {p:10.3e}  {n:>4}  {sig}")

# Write full table
with open(OUT_CORR, "w") as out:
    out.write("te_metric\tphenotype\tspearman_r\tp_value\tn\tsignificant\n")
    for te_col, ph_col, r, p, n, sig in results:
        out.write(f"{te_col}\t{ph_col}\t{r:.4f}\t{p:.4e}\t{n}\t{sig}\n")

print(f"\nFull table → {OUT_CORR}")
print(f"Sirtuin table → {OUT_SIRTUIN}")

# ─── Summary: strongest TE × phenotype correlations ─────────────────────────
print("\n" + "="*60)
print("TOP SIGNIFICANT TE × PHENOTYPE CORRELATIONS")
print("="*60)
sig_results = [(r, p, te_col, ph_col, sig, n)
               for te_col, ph_col, r, p, n, sig in results
               if p < 0.05]
sig_results.sort(key=lambda x: -abs(x[0]))
for r, p, te_col, ph_col, sig, n in sig_results[:25]:
    print(f"  {te_col:<23} × {ph_col:<28}  r={r:+.3f}  p={p:.3e}  {sig}")
