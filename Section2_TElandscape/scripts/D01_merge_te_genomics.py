#!/usr/bin/env python3
"""
D01_merge_te_genomics.py
Join EarlGrey TE summary (highLevelCount.txt) with GenomicFeatures V3
to prepare the merged table for PGLS regression (D02).

Usage:
    python D01_merge_te_genomics.py \
        <earlgrey_annotated_dir> \
        <genomicfeatures_summary_tsv> \
        <output_merged_tsv>

Input:
    annotated_dir: contains {SP}_EarlGrey/{SP}_summaryFiles/{SP}.highLevelCount.txt
    genomicfeatures_summary_tsv: ALL_SPECIES_ALL_SUMMARY.txt
        cols: Species, Genome_Size_bp, Total_Genes, Total_Isoforms,
              Avg_Isoforms_Per_Gene, Total_Exons, Total_Introns,
              Avg_Exons_Per_Gene, Avg_Introns_All, Avg_Introns_Multiexon,
              Pct_Genic, Pct_Exonic, Pct_Intronic, Pct_Intergenic

Output:
    TSV with TE class percentages + genomic features per species
"""

import sys
import os
import re
import pandas as pd


# Map EarlGrey high-level class labels → short column names
CLASS_MAP = {
    "DNA": "pct_DNA",
    "LTR": "pct_LTR",
    "LINE": "pct_LINE",
    "SINE": "pct_SINE",
    "Rolling Circle": "pct_RC",
    "Unknown": "pct_Unknown",
    "Satellite": "pct_Satellite",
    "Simple Repeat": "pct_SimpleRepeat",
    "Low Complexity": "pct_LowComplexity",
}

NESTED_RE = re.compile(r"-nested$", re.IGNORECASE)


def load_highlevel(path):
    """
    Parse {SP}.highLevelCount.txt.
    Columns: TE Classification, Coverage (bp), Copy Number,
             % Genome Coverage, Genome Size, TE Family Count
    Returns dict: class_label → pct_genome_coverage
    Nested classes (e.g. DNA-nested) are added to their parent.
    """
    data = {}
    genome_size = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("TE Classification"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            te_class = parts[0].strip()
            try:
                pct = float(parts[3])
                gs = int(parts[4]) if len(parts) > 4 else None
            except ValueError:
                continue
            if genome_size is None and gs:
                genome_size = gs
            # Merge nested into parent class
            base_class = NESTED_RE.sub("", te_class)
            data[base_class] = data.get(base_class, 0.0) + pct

    total_te = sum(data.values())
    return data, total_te, genome_size


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: D01_merge_te_genomics.py <annotated_dir> <genomicfeatures_tsv> <output_tsv>")

    annotated_dir = sys.argv[1]
    gf_tsv = sys.argv[2]
    output_tsv = sys.argv[3]

    # ── load genomic features ─────────────────────────────────
    gf = pd.read_csv(gf_tsv, sep="\t")
    gf["Species"] = gf["Species"].str.strip()
    gf_indexed = gf.set_index("Species")

    # ── collect TE data per species ───────────────────────────
    te_records = []
    missing_gf = []

    for entry in sorted(os.listdir(annotated_dir)):
        sp_dir = os.path.join(annotated_dir, entry)
        if not os.path.isdir(sp_dir) or not entry.endswith("_EarlGrey"):
            continue
        sp = entry.replace("_EarlGrey", "")
        summary_dir = os.path.join(sp_dir, f"{sp}_summaryFiles")
        hlc = os.path.join(summary_dir, f"{sp}.highLevelCount.txt")
        if not os.path.isfile(hlc):
            continue

        class_data, total_te, genome_size = load_highlevel(hlc)

        row = {"species": sp, "pct_TE_total": round(total_te, 4)}
        for label, col in CLASS_MAP.items():
            row[col] = round(class_data.get(label, 0.0), 4)
        if genome_size:
            row["genome_size_TE"] = genome_size

        te_records.append(row)

    te_df = pd.DataFrame(te_records)

    # ── merge with genomic features ───────────────────────────
    merged_rows = []
    for _, row in te_df.iterrows():
        sp = row["species"]
        if sp not in gf_indexed.index:
            missing_gf.append(sp)
            # Keep TE data even without genomic features
            merged_rows.append(row.to_dict())
        else:
            gf_row = gf_indexed.loc[sp]
            merged = {**row.to_dict(), **gf_row.to_dict()}
            merged_rows.append(merged)

    merged_df = pd.DataFrame(merged_rows)

    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)
    merged_df.to_csv(output_tsv, sep="\t", index=False)

    print(f"[OK] Merged TE + genomics table: {output_tsv}")
    print(f"[INFO] Species with TE data: {len(te_records)}")
    print(f"[INFO] Species with both TE + genomic features: {len(te_records) - len(missing_gf)}")
    if missing_gf:
        print(f"[WARN] No genomic features for {len(missing_gf)} species: {missing_gf[:5]}")


if __name__ == "__main__":
    main()
