#!/usr/bin/env python3
"""
B01_shannon_diversity.py
Compute Shannon diversity (H') and Pielou's evenness (J) of TE superfamilies
per species from EarlGrey familyLevelCount.txt files.

Usage:
    python B01_shannon_diversity.py <earlgrey_annotated_dir> <output_tsv>

Input:
    Directory containing {SP}_EarlGrey/{SP}_summaryFiles/{SP}.familyLevelCount.txt
    Family names follow the pattern: NAME#Class/Superfamily or NAME#Class

Output:
    TSV with columns: species, H_prime, J_evenness, n_superfamilies,
                      dominant_superfamily, dominant_pct, total_te_bp
"""

import sys
import os
import re
import math
from collections import defaultdict
import pandas as pd

# Non-TE repeat classes to exclude from diversity calculations.
# These appear as NAME#Simple_repeat, NAME#Low_complexity, etc. in familyLevelCount.txt.
NON_TE_CLASSES = frozenset({
    "Simple_repeat", "Low_complexity", "Satellite",
    "rRNA", "tRNA", "snRNA", "scRNA", "sno", "srpRNA",
    "ARTEFACT",
})


def parse_superfamily(family_name):
    """Extract (te_class, superfamily) from EarlGrey family label (NAME#Class/Superfamily)."""
    m = re.search(r"#([^/]+)/(.+)$", family_name)
    if m:
        return m.group(1), m.group(1) + "/" + m.group(2)
    m = re.search(r"#(.+)$", family_name)
    if m:
        return m.group(1), m.group(1)
    return "Unknown", "Unknown"


def shannon(counts):
    """Shannon H' from a list of counts."""
    total = sum(counts)
    if total == 0:
        return 0.0
    probs = [c / total for c in counts if c > 0]
    return -sum(p * math.log(p) for p in probs)


def pielou(h, n):
    """Pielou's evenness J = H' / ln(S)."""
    if n <= 1:
        return float("nan")
    return h / math.log(n)


def process_species(family_count_file):
    """Read familyLevelCount.txt and return superfamily-level coverage dict (TEs only)."""
    sf_bp = defaultdict(float)
    with open(family_count_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("TE Family"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            family = parts[0]
            try:
                bp = float(parts[1])
            except ValueError:
                continue
            te_class, sf = parse_superfamily(family)
            # Skip non-TE repeats (structural RNAs, simple repeats, artefacts)
            if te_class in NON_TE_CLASSES:
                continue
            sf_bp[sf] += bp
    return sf_bp


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: B01_shannon_diversity.py <annotated_dir> <output_tsv>")

    annotated_dir = sys.argv[1]
    output_tsv = sys.argv[2]

    records = []
    missing = []

    for entry in sorted(os.listdir(annotated_dir)):
        sp_dir = os.path.join(annotated_dir, entry)
        if not os.path.isdir(sp_dir) or not entry.endswith("_EarlGrey"):
            continue
        sp = entry.replace("_EarlGrey", "")
        summary_dir = os.path.join(sp_dir, f"{sp}_summaryFiles")
        fcf = os.path.join(summary_dir, f"{sp}.familyLevelCount.txt")
        if not os.path.isfile(fcf):
            missing.append(sp)
            continue

        sf_bp = process_species(fcf)
        if not sf_bp:
            missing.append(sp)
            continue

        counts = list(sf_bp.values())
        total_bp = sum(counts)
        h = shannon(counts)
        n = len(counts)
        j = pielou(h, n)
        dominant = max(sf_bp, key=sf_bp.get)
        dom_pct = sf_bp[dominant] / total_bp * 100 if total_bp > 0 else 0.0

        records.append({
            "species": sp,
            "H_prime": round(h, 4),
            "J_evenness": round(j, 4) if not math.isnan(j) else "NA",
            "n_superfamilies": n,
            "dominant_superfamily": dominant,
            "dominant_pct": round(dom_pct, 2),
            "total_te_bp": int(total_bp),
        })

    df = pd.DataFrame(records)
    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)
    df.to_csv(output_tsv, sep="\t", index=False)

    print(f"[OK] Shannon diversity written: {output_tsv}")
    print(f"[INFO] Species processed: {len(records)}")
    if missing:
        print(f"[WARN] familyLevelCount.txt not found for {len(missing)} species: {missing[:5]}{'...' if len(missing) > 5 else ''}")


if __name__ == "__main__":
    main()
