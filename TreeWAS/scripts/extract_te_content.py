#!/usr/bin/env python3
"""
Extract TE content (% masked) per species from EarlGrey RepeatMasker .tbl files.

Input:  /data/carlos/EDF-TE/TE_analysis/annotated/*_EarlGrey/
Output: /data/carlos/EDF-TE/EDF-TE/TreeWAS/data/te_content_per_species.tsv

Only EDF species are kept (prefixes: BL, BS, CH, CR, EN, GL, KX, MI, MO, MU, ZO).
"""

import os
import re
import sys
from collections import defaultdict

# ─── Paths ───────────────────────────────────────────────────────────────────
ANNOTATED_DIR = "/data/carlos/EDF-TE/TE_analysis/annotated"
OUTPUT_TSV    = "/data/carlos/EDF-TE/EDF-TE/TreeWAS/data/te_content_per_species.tsv"

EDF_PREFIXES = {"BL", "BS", "CH", "CR", "EN", "GL", "KX", "MI", "MO", "MU", "ZO"}

# ─── Regexes for parsing RepeatMasker .tbl ───────────────────────────────────
# Format: "total length:   57060573 bp"
RE_GENOME   = re.compile(r"total length:\s+([\d,]+)\s+bp", re.IGNORECASE)
# Format: "bases masked:    5585041 bp ( 9.79 %)"
RE_BASES_MASKED = re.compile(r"bases masked:\s+[\d,]+\s+bp\s*\(\s*([\d.]+)\s*%\)", re.IGNORECASE)
# Format: "Total interspersed repeats:     4070064 bp    7.13 %"
RE_MASKED   = re.compile(r"Total interspersed repeats:\s+[\d,]+\s+bp\s+([\d.]+)\s*%", re.IGNORECASE)
# Format: "   SINEs:              298       158010 bp    0.28 %"
RE_SINELINES = re.compile(r"^\s*SINEs:\s+\d+\s+[\d,]+\s+bp\s+([\d.]+)\s*%", re.MULTILINE)
RE_LINES    = re.compile(r"^\s*LINEs:\s+\d+\s+[\d,]+\s+bp\s+([\d.]+)\s*%", re.MULTILINE)
RE_LTRS     = re.compile(r"^\s*LTR elements:\s+\d+\s+[\d,]+\s+bp\s+([\d.]+)\s*%", re.MULTILINE)
RE_DNA      = re.compile(r"^\s*DNA transposons\s+\d+\s+[\d,]+\s+bp\s+([\d.]+)\s*%", re.MULTILINE)
RE_HELITRON = re.compile(r"^\s*Rolling-circles\s+\d+\s+[\d,]+\s+bp\s+([\d.]+)\s*%", re.MULTILINE)
RE_UNCLASS  = re.compile(r"^\s*Unclassified:\s+\d+\s+[\d,]+\s+bp\s+([\d.]+)\s*%", re.MULTILINE)


def parse_tbl(tbl_path: str) -> dict:
    """Parse a RepeatMasker .tbl file and return a dict of TE metrics."""
    text = open(tbl_path, errors="replace").read()

    def _get(rx, default=0.0):
        m = rx.search(text)
        if m:
            return float(m.group(1).replace(",", ""))
        return default

    genome_bp = _get(RE_GENOME)
    return {
        "genome_bp":            int(genome_bp),
        "genome_Mb":            round(genome_bp / 1e6, 2),
        "bases_masked_pct":     _get(RE_BASES_MASKED),   # all repeats incl. simple
        "total_masked_pct":     _get(RE_MASKED),          # interspersed TEs only
        "LTR_pct":              _get(RE_LTRS),
        "LINE_pct":             _get(RE_LINES),
        "SINE_pct":             _get(RE_SINELINES),
        "DNA_pct":              _get(RE_DNA),
        "Helitron_pct":         _get(RE_HELITRON),
        "Unclassified_pct":     _get(RE_UNCLASS),
    }


def find_tbl(eg_dir: str) -> str | None:
    """Find the RepeatMasker .tbl file inside an EarlGrey directory tree.

    Prefers the _Against_Custom_Library directory (final EarlGrey run);
    falls back to any other .tbl found.
    """
    custom_lib = None
    fallback   = None
    for root, _dirs, files in os.walk(eg_dir):
        for f in files:
            if f.endswith(".tbl"):
                path = os.path.join(root, f)
                if "Against_Custom_Library" in root:
                    custom_lib = path
                elif fallback is None:
                    fallback = path
    return custom_lib or fallback


# ─── Main ─────────────────────────────────────────────────────────────────────
rows = []

if not os.path.isdir(ANNOTATED_DIR):
    sys.exit(f"ERROR: {ANNOTATED_DIR} does not exist")

for entry in sorted(os.listdir(ANNOTATED_DIR)):
    if not entry.endswith("_EarlGrey"):
        continue
    # Species code = entry minus the "_EarlGrey" suffix
    species = entry[: -len("_EarlGrey")]
    prefix  = species.split("_")[0] if "_" in species else species[:2].upper()

    if prefix not in EDF_PREFIXES:
        print(f"  SKIP (non-EDF): {species}", file=sys.stderr)
        continue

    eg_dir = os.path.join(ANNOTATED_DIR, entry)
    tbl    = find_tbl(eg_dir)
    if tbl is None:
        print(f"  WARNING: no .tbl file found in {eg_dir}", file=sys.stderr)
        continue

    try:
        metrics = parse_tbl(tbl)
    except Exception as e:
        print(f"  ERROR parsing {tbl}: {e}", file=sys.stderr)
        continue

    rows.append({"species": species, "prefix": prefix, **metrics})
    print(f"  {species}  {metrics['genome_Mb']} Mb  "
          f"masked={metrics['total_masked_pct']:.2f}%  "
          f"LTR={metrics['LTR_pct']:.2f}%  LINE={metrics['LINE_pct']:.2f}%  "
          f"DNA={metrics['DNA_pct']:.2f}%", file=sys.stderr)

if not rows:
    sys.exit("No EDF species with EarlGrey results found.")

# ─── Write TSV ────────────────────────────────────────────────────────────────
os.makedirs(os.path.dirname(OUTPUT_TSV), exist_ok=True)
cols = ["species", "prefix", "genome_bp", "genome_Mb",
        "bases_masked_pct", "total_masked_pct", "LTR_pct", "LINE_pct", "SINE_pct",
        "DNA_pct", "Helitron_pct", "Unclassified_pct"]

with open(OUTPUT_TSV, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in rows:
        out.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")

print(f"\nWrote {len(rows)} species → {OUTPUT_TSV}", file=sys.stderr)

# ─── Summary by clade ────────────────────────────────────────────────────────
print("\n" + "="*60)
print("TE CONTENT BY CLADE (mean % masked)")
print("="*60)

by_clade = defaultdict(list)
for r in rows:
    by_clade[r["prefix"]].append(r)

all_clades = sorted(by_clade.keys())
print(f"\n{'Clade':<8} {'n':>4} {'Genome(Mb)':>12} {'Total%':>8} "
      f"{'LTR%':>7} {'LINE%':>7} {'DNA%':>7} {'Uncl%':>7}")
print("-" * 65)
for clade in all_clades:
    sp_list = by_clade[clade]
    n = len(sp_list)

    def mean(col):
        vals = [r[col] for r in sp_list if r[col] > 0]
        return sum(vals) / len(vals) if vals else 0.0

    print(f"{clade:<8} {n:>4} {mean('genome_Mb'):>12.1f} {mean('total_masked_pct'):>8.2f} "
          f"{mean('LTR_pct'):>7.2f} {mean('LINE_pct'):>7.2f} "
          f"{mean('DNA_pct'):>7.2f} {mean('Unclassified_pct'):>7.2f}")

print("\nAll clades:")
all_sp = rows
n = len(all_sp)
print(f"  n={n}  mean total masked={(sum(r['total_masked_pct'] for r in all_sp)/n):.2f}%  "
      f"LTR={(sum(r['LTR_pct'] for r in all_sp)/n):.2f}%  "
      f"LINE={(sum(r['LINE_pct'] for r in all_sp)/n):.2f}%  "
      f"DNA={(sum(r['DNA_pct'] for r in all_sp)/n):.2f}%")
