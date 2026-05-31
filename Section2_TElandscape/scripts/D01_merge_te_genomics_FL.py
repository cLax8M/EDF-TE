#!/usr/bin/env python3
"""
D01_merge_te_genomics_FL.py — Same as D01_merge_te_genomics.py but using
ONLY full-length copies (copy_len >= FL_FRAC * consensus_len, default 0.8).

For each species we re-walk the EarlGrey filteredRepeats.gff together with
the combined_library.fasta to (a) classify each family by high-level class
from its FASTA header (FAMILY#CLASS/SUBCLASS), (b) compute consensus length,
(c) keep only FL copies, and (d) sum bp per class. Percentages are then
expressed as a fraction of genome size (taken from highLevelCount.txt to be
consistent with the non-FL pipeline).

Usage:
    python D01_merge_te_genomics_FL.py <annotated_dir> <gf_tsv> <output_tsv> [--fl_frac 0.8]
"""

import sys, os, re, argparse
import pandas as pd
from collections import defaultdict


# Map EarlGrey high-level class → short column names (same as D01).
CLASS_MAP = {
    "DNA":            "pct_DNA",
    "LTR":            "pct_LTR",
    "LINE":           "pct_LINE",
    "SINE":           "pct_SINE",
    "RC":             "pct_RC",          # rolling-circle / Helitron
    "Helitron":       "pct_RC",          # alias
    "PLE":            "pct_PLE",
    "Penelope":       "pct_PLE",
    "Unknown":        "pct_Unknown",
    "Unclassified":   "pct_Unknown",
}
TE_CLASSES = frozenset({"pct_DNA", "pct_LTR", "pct_LINE", "pct_SINE",
                        "pct_RC", "pct_PLE", "pct_Unknown"})

LEGACY_TO_CANONICAL = {
    "AS_AspfumUD": "AS_AspfumUD1",
    "NU_Parvatl":  "NU_Parvatl_c",
}

_NAME_RE = re.compile(r'Name=([^;]+)')


def load_consensus(lib_fa):
    """Return dict family_name -> (consensus_length, class_short)."""
    out = {}
    name = cls = None
    L = 0
    with open(lib_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = (L, cls)
                hdr = line[1:].rstrip().split()[0]
                # FAMILY#CLASS/SUBCLASS
                if "#" in hdr:
                    fam_raw, cls_full = hdr.split("#", 1)
                    cls_top = cls_full.split("/")[0]
                else:
                    fam_raw, cls_top = hdr, "Unknown"
                name = fam_raw.upper()
                cls  = CLASS_MAP.get(cls_top, "pct_Unknown")
                L = 0
            else:
                L += len(line.rstrip())
        if name is not None:
            out[name] = (L, cls)
    return out


def parse_genome_size(hlc_path):
    gs = None
    with open(hlc_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 4:
                try:
                    gs = int(parts[4]); break
                except ValueError:
                    continue
    return gs


def process_species(eg_dir, fl_frac):
    base = os.path.basename(eg_dir)
    sp = base.replace("_EarlGrey", "")
    summary = os.path.join(eg_dir, f"{sp}_summaryFiles")
    gff = os.path.join(summary, f"{sp}.filteredRepeats.gff")
    lib = os.path.join(summary, f"{sp}_combined_library.fasta")
    hlc = os.path.join(summary, f"{sp}.highLevelCount.txt")
    if not (os.path.isfile(gff) and os.path.isfile(lib) and os.path.isfile(hlc)):
        return None

    cons = load_consensus(lib)
    gs = parse_genome_size(hlc)
    if not gs:
        return None

    bp_by_class = defaultdict(int)
    with open(gff) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            try:
                start = int(cols[3]); end = int(cols[4])
            except ValueError:
                continue
            m = _NAME_RE.search(cols[8])
            if not m:
                continue
            fam = m.group(1).split("#")[0].upper()
            info = cons.get(fam)
            if not info:
                continue
            cl, cls = info
            if cl == 0:
                continue
            copy_len = end - start + 1
            if copy_len >= fl_frac * cl:
                bp_by_class[cls] += copy_len

    row = {"species": LEGACY_TO_CANONICAL.get(sp, sp)}
    for col in ("pct_DNA", "pct_LTR", "pct_LINE", "pct_SINE",
                "pct_RC", "pct_PLE", "pct_Unknown"):
        row[col] = round(100.0 * bp_by_class.get(col, 0) / gs, 4)
    row["pct_TE_total"] = round(sum(row[c] for c in TE_CLASSES), 4)
    row["genome_size_TE"] = gs
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("annotated_dir")
    ap.add_argument("gf_tsv")
    ap.add_argument("out_tsv")
    ap.add_argument("--fl_frac", type=float, default=0.8)
    args = ap.parse_args()

    gf = pd.read_csv(args.gf_tsv, sep="\t")
    gf["Species"] = gf["Species"].str.strip()
    gf_idx = gf.set_index("Species")

    rows = []
    for entry in sorted(os.listdir(args.annotated_dir)):
        if not entry.endswith("_EarlGrey"):
            continue
        full = os.path.join(args.annotated_dir, entry)
        if not os.path.isdir(full):
            continue
        r = process_species(full, args.fl_frac)
        if r:
            rows.append(r)

    te_df = pd.DataFrame(rows)
    merged = []
    for _, r in te_df.iterrows():
        sp = r["species"]
        d = r.to_dict()
        if sp in gf_idx.index:
            d.update(gf_idx.loc[sp].to_dict())
        merged.append(d)

    out = pd.DataFrame(merged)
    # Drop VL (control yeasts) globally
    if "species" in out.columns:
        out = out[~out["species"].str.startswith("VL_")].reset_index(drop=True)

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[OK] Wrote FL-only merged table: {args.out_tsv}  ({len(out)} species)")


if __name__ == "__main__":
    main()
