#!/usr/bin/env python3
"""
12_te_proximity_genes.py
Distance of genes to nearest TE insertion, by gene category.

For each species with a filteredRepeats.bed file, computes the distance from
each annotated gene to its nearest TE using a pure-Python bisect approach
(no external bedtools module required).

Gene categories:
  - orphan          : OrthoFinder unassigned genes (from script 11)
  - coopted         : TE-coopted candidate genes (true_TE_coopted_genes.tsv)
  - unknown_fn      : genes encoding only Unknown-Function / DUF domains
  - other           : all remaining genes

Inputs (auto-resolved):
  - GFF annotations   : /legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff/{SP}.gff[3]
  - filteredRepeats   : TE_analysis/annotated/{SP}_EarlGrey/{SP}_summaryFiles/{SP}.filteredRepeats.bed
  - Orthogroups_UnassignedGenes.tsv (orphan list)
  - per_species/classification/{SP}.true_TE_coopted_genes.tsv (coopted genes)
  - per_species/domains/{SP}.domains.tsv + TE_domains_master.tsv

Outputs:
  - figures/12_te_proximity_per_gene.tsv.gz   (per-gene distances)
  - figures/12_te_proximity_summary.tsv        (per-category stats per species)
  - figures/12_te_proximity_boxplot.pdf

Design notes:
  - Species missing BED are **skipped** (logged) and added when EarlGrey finishes.
  - Re-runnable; simply re-run after all EarlGrey jobs complete.
  - Uses bisect_left on sorted TE midpoint lists per scaffold for O(log n) lookups.
  - Gene midpoint = (gene_start + gene_end) / 2 (0-based BED coordinates).
  - Minimum distance = 0 when gene overlaps a TE.

Usage:
  python3 12_te_proximity_genes.py
"""

import os
import sys
import gzip
import bisect
import collections
import warnings
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ─────────────────────────────────────────────────────────────────────
DATA_BASE      = "/legserv/Temp/Carlos/EDF-TE"
COOPTION_DIR   = os.path.join(DATA_BASE, "results", "Co-option")
DOMAINS_DIR    = os.path.join(COOPTION_DIR, "per_species", "domains")
CLASS_DIR      = os.path.join(COOPTION_DIR, "per_species", "classification")
TE_MASTER      = os.path.join(COOPTION_DIR, "config", "TE_domains_master.tsv")
OG_UNASSIGNED  = os.path.join(DATA_BASE, "results", "Orthofinder_run",
                               "Orthogroups", "Orthogroups_UnassignedGenes.tsv")
GFF_DIR        = os.path.join(DATA_BASE, "data", "genes_annot_gff")
# EarlGrey BED paths are resolved per-species
TE_ANNOT_BASE  = os.path.join(DATA_BASE, "data", "De_novo_TE_annotation_non_fungal")

# Alternative BED discovery (EarlGrey output layout may vary):
EARLGREY_BASE  = "/data/carlos/EDF-TE/TE_analysis/annotated"

OUT_DIR        = os.path.join(COOPTION_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PERGENE  = os.path.join(OUT_DIR, "12_te_proximity_per_gene.tsv.gz")
OUT_SUMMARY  = os.path.join(OUT_DIR, "12_te_proximity_summary.tsv")
OUT_FIG      = os.path.join(OUT_DIR, "12_te_proximity_boxplot.pdf")

CATEGORIES = ["orphan", "coopted", "unknown_fn", "other"]
CAT_COLORS  = {
    "orphan":     "#d62728",
    "coopted":    "#ff7f0e",
    "unknown_fn": "#9467bd",
    "other":      "#aec7e8",
}

# ── helpers ───────────────────────────────────────────────────────────────────

def find_bed_file(sp):
    """
    Locate filteredRepeats.bed for a species.
    Tries multiple path patterns from EarlGrey output layout.
    """
    patterns = [
        os.path.join(EARLGREY_BASE, f"{sp}_EarlGrey",
                     f"{sp}_summaryFiles", f"{sp}.filteredRepeats.bed"),
        os.path.join(EARLGREY_BASE, f"{sp}_EarlGrey",
                     f"{sp}_summaryFiles", f"{sp}.filteredRepeats.bed.gz"),
        os.path.join(TE_ANNOT_BASE, sp, f"{sp}.filteredRepeats.bed"),
    ]
    for p in patterns:
        if os.path.exists(p):
            return p
    return None


def find_gff_file(sp):
    """Locate GFF/GFF3 annotation. Returns path or None."""
    for ext in [".gff", ".gff3", ".gff.gz", ".gff3.gz"]:
        p = os.path.join(GFF_DIR, sp + ext)
        if os.path.exists(p):
            return p
    return None


def open_file(path):
    """Open plain or gzipped file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_te_domains_set(master_path):
    df = pd.read_csv(master_path, sep="\t")
    return set(df["domain_id"].str.strip().str.upper())


def load_orphans(og_path):
    """Returns dict: sp → set of orphan gene IDs."""
    warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)
    df = pd.read_csv(og_path, sep="\t", index_col=0, low_memory=False)
    out = {}
    for col in df.columns:
        sp = col.replace(".aa", "")
        genes = set(df[col].dropna().astype(str).str.strip())
        if genes:
            out[sp] = genes
    return out


def load_coopted(class_dir, sp):
    """Returns set of coopted gene IDs for species sp."""
    fpath = os.path.join(class_dir, f"{sp}.true_TE_coopted_genes.tsv")
    if not os.path.exists(fpath):
        return set()
    df = pd.read_csv(fpath, sep="\t")
    return set(df["gene_id"].astype(str))


def load_unknown_fn_genes(domains_dir, te_ids, sp):
    """
    Returns set of gene IDs that ONLY have unknown-function domains
    (no TE domains, no Pfam/SUPERFAMILY functional domains).
    A gene is 'unknown_fn' if it has NO recognized non-TE, non-DUF domains.
    Operationally: all its domain names start with 'DUF', 'UPF' or contain
    'domain of unknown' (case-insensitive), AND no TE domains.
    """
    dom_path = os.path.join(domains_dir, f"{sp}.domains.tsv")
    if not os.path.exists(dom_path):
        return set()
    df = pd.read_csv(dom_path, sep="\t")
    result = set()
    for gene_id, grp in df.groupby("gene_id"):
        ids   = grp["domain_id"].fillna("").str.upper().tolist()
        names = grp["domain_name"].fillna("").str.lower().tolist()
        has_te = any(i in te_ids for i in ids)
        if has_te:
            continue
        # All domains must qualify as "unknown"
        all_unk = all(
            n.startswith("duf") or n.startswith("upf") or "unknown" in n
            for n in names
        )
        if all_unk:
            result.add(str(gene_id))
    return result


def load_te_bed(bed_path):
    """
    Load filteredRepeats.bed.
    Returns dict: scaffold → sorted list of (start, end) tuples.
    Coordinates are 0-based.
    """
    te_by_scaffold = collections.defaultdict(list)
    with open_file(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            te_by_scaffold[chrom].append((start, end))
    # Sort by start
    for chrom in te_by_scaffold:
        te_by_scaffold[chrom].sort()
    return dict(te_by_scaffold)


def load_genes_gff(gff_path):
    """
    Parse GFF3/GFF and return dict: gene_id → (chrom, start, end) with 0-based coords.

    Two-pass strategy:
    1. Parse 'gene' features → locus_tag → (chrom, start, end)
    2. Parse 'CDS' features → protein_id → parent locus_tag (via Parent chain)
    3. Build protein_id → (chrom, start, end) using gene coords

    Also records protein_id directly as keys when protein_id == gene_id (simple genomes).
    Falls back to using 'gene' feature ID directly if no CDS protein_id mapping is found.
    """
    # Pass 1: collect gene-level coords indexed by all possible IDs from that feature
    gene_coords   = {}   # raw ID (e.g. gene-MP638_005444) → (chrom, s, e)
    locus_to_id   = {}   # locus_tag → raw gene ID

    # Pass 2: CDS protein_id → set of parent mRNA IDs
    protein_to_gene = {}  # protein_id → raw gene ID (via parent chain)

    # mRNA parent mapping: mRNA ID → gene ID
    rna_to_gene = {}

    with open_file(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat  = parts[2]
            chrom = parts[0]
            start = int(parts[3]) - 1
            end   = int(parts[4])
            attrs_str = parts[8]
            attrs = {}
            for a in attrs_str.split(";"):
                a = a.strip()
                if "=" in a:
                    k, v = a.split("=", 1)
                    attrs[k.strip()] = v.strip()

            if feat == "gene":
                raw_id = attrs.get("ID", "")
                gene_coords[raw_id] = (chrom, start, end)
                lt = attrs.get("locus_tag", "")
                if lt:
                    locus_to_id[lt] = raw_id
                # Also index by Name
                nm = attrs.get("Name", "")
                if nm and nm != raw_id:
                    locus_to_id[nm] = raw_id

            elif feat in ("mRNA", "transcript"):
                rna_id  = attrs.get("ID", "")
                parent  = attrs.get("Parent", "")
                if rna_id and parent:
                    rna_to_gene[rna_id] = parent
                # Also by locus_tag
                lt = attrs.get("locus_tag", "")
                if lt and parent:
                    rna_to_gene[lt] = parent

            elif feat == "CDS":
                pid    = attrs.get("protein_id", "")
                parent = attrs.get("Parent", "")
                if pid and parent:
                    # parent is a mRNA/transcript ID; we'll resolve it post-loop
                    protein_to_gene[pid] = parent

    # Resolve protein_id → gene coords through the rna_to_gene chain
    genes = {}

    for pid, rna_id in protein_to_gene.items():
        gene_raw_id = rna_to_gene.get(rna_id, rna_id)  # may be gene ID directly
        coords = gene_coords.get(gene_raw_id)
        if coords is None:
            # Try stripping prefixes
            for prefix in ("gene-", "gene:", "rna-"):
                if gene_raw_id.startswith(prefix):
                    coords = gene_coords.get(gene_raw_id)
                    break
        if coords:
            genes[pid] = coords

    # Fallback: if no protein→gene mapping found, use gene feature IDs directly
    # (handles simple genomes where gene ID IS the protein/domain ID)
    if not genes:
        for raw_id, coords in gene_coords.items():
            clean = raw_id
            for prefix in ("gene-", "gene:", "rna-", "transcript-"):
                if clean.startswith(prefix):
                    clean = clean[len(prefix):]
            genes[clean] = coords

    return genes


def dist_to_nearest_te(chrom, gstart, gend, te_by_scaffold):
    """
    Compute minimum distance from gene [gstart, gend) to nearest TE on same scaffold.
    Uses bisect on sorted TE start positions.
    Returns:
      - int distance (0 = overlap)
      - None if scaffold not in te_by_scaffold
    """
    tes = te_by_scaffold.get(chrom)
    if not tes:
        return None

    gene_mid = (gstart + gend) // 2

    # Bisect on TE starts
    starts = [t[0] for t in tes]
    idx = bisect.bisect_left(starts, gene_mid)

    best = float("inf")
    for i in range(max(0, idx - 2), min(len(tes), idx + 3)):
        ts, te = tes[i]
        if gstart >= ts and gend <= te + 1:
            return 0   # gene inside TE
        if ts > gend:
            best = min(best, ts - gend)
        elif te < gstart:
            best = min(best, gstart - te)
        else:
            return 0   # overlap
    return int(best) if best < float("inf") else None


def categorize(gene_id, orphan_set, coopted_set, unknown_fn_set):
    if gene_id in coopted_set:
        return "coopted"
    if gene_id in orphan_set:
        return "orphan"
    if gene_id in unknown_fn_set:
        return "unknown_fn"
    return "other"


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    print("[INFO] Loading shared data...")
    te_ids    = load_te_domains_set(TE_MASTER)
    all_orphans = load_orphans(OG_UNASSIGNED)

    # Detect species with domains.tsv
    all_species = sorted(
        f.replace(".domains.tsv", "")
        for f in os.listdir(DOMAINS_DIR)
        if f.endswith(".domains.tsv")
    )
    print(f"[INFO] {len(all_species)} species with domain data")

    per_gene_rows = []
    summary_rows  = []
    skipped  = []
    no_bed   = []

    for sp in all_species:
        bed_path = find_bed_file(sp)
        if not bed_path:
            no_bed.append(sp)
            continue

        gff_path = find_gff_file(sp)
        if not gff_path:
            skipped.append(f"{sp}  (no GFF)")
            continue

        # Load TE BED
        try:
            te_by_scaffold = load_te_bed(bed_path)
        except Exception as e:
            skipped.append(f"{sp}  (BED error: {e})")
            continue

        # Load gene coords
        try:
            genes = load_genes_gff(gff_path)
        except Exception as e:
            skipped.append(f"{sp}  (GFF error: {e})")
            continue

        if not genes:
            skipped.append(f"{sp}  (0 genes in GFF)")
            continue

        # Category sets
        orphan_set     = all_orphans.get(sp, set())
        coopted_set    = load_coopted(CLASS_DIR, sp)
        unknown_fn_set = load_unknown_fn_genes(DOMAINS_DIR, te_ids, sp)

        # Compute distances
        cat_dists = collections.defaultdict(list)
        for gene_id, (chrom, gstart, gend) in genes.items():
            d = dist_to_nearest_te(chrom, gstart, gend, te_by_scaffold)
            if d is None:
                continue
            cat = categorize(gene_id, orphan_set, coopted_set, unknown_fn_set)
            cat_dists[cat].append(d)
            per_gene_rows.append({
                "species":     sp,
                "gene_id":     gene_id,
                "category":    cat,
                "distance_bp": d,
            })

        # Per-species summary stats per category
        cats_used = set(cat_dists.keys())
        for cat in CATEGORIES:
            d_arr = cat_dists.get(cat, [])
            n = len(d_arr)
            if n == 0:
                continue
            summary_rows.append({
                "species":    sp,
                "category":   cat,
                "n_genes":    n,
                "median_dist":  int(np.median(d_arr)),
                "mean_dist":    int(np.mean(d_arr)),
                "pct_overlap":  round(100.0 * sum(1 for x in d_arr if x == 0) / n, 3),
            })

        # Mann-Whitney U between orphan and other (if both present)
        if "orphan" in cats_used and "other" in cats_used:
            u, p = mannwhitneyu(cat_dists["orphan"], cat_dists["other"],
                                alternative="two-sided")
        sys.stderr.write(f"  {sp}  genes:{len(genes)}  "
                         f"orphan:{len(cat_dists.get('orphan',[]))}  "
                         f"coopted:{len(cat_dists.get('coopted',[]))}\n")

    # ── results ──────────────────────────────────────────────────────────────
    df_gene = pd.DataFrame(per_gene_rows)
    df_gene.to_csv(OUT_PERGENE, sep="\t", index=False, compression="gzip")
    print(f"\n[OUT] Per-gene distances → {OUT_PERGENE}")

    df_sum = pd.DataFrame(summary_rows)
    df_sum.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"[OUT] Per-category summary → {OUT_SUMMARY}")

    if no_bed:
        print(f"\n[INFO] Skipped {len(no_bed)} species (no BED file yet):")
        for s in sorted(no_bed):
            print(f"  {s}")

    if skipped:
        print(f"\n[WARN] Skipped {len(skipped)} species (other error):")
        for s in skipped:
            print(f"  {s}")

    # ── figure ────────────────────────────────────────────────────────────────
    if df_gene.empty:
        print("[WARN] No data for figure. Exiting.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: global distance distribution per category (log-scale violin)
    ax = axes[0]
    data_by_cat = []
    labels_used = []
    for cat in CATEGORIES:
        sub = df_gene[df_gene["category"] == cat]["distance_bp"].values
        if len(sub) > 0:
            data_by_cat.append(np.log10(sub + 1))
            labels_used.append(cat)
    vp = ax.violinplot(data_by_cat, positions=range(len(labels_used)),
                       showmedians=True, widths=0.6)
    for i, cat in enumerate(labels_used):
        pc = vp["bodies"][i]
        pc.set_facecolor(CAT_COLORS[cat])
        pc.set_alpha(0.7)
    ax.set_xticks(range(len(labels_used)))
    ax.set_xticklabels(labels_used, fontsize=10)
    ax.set_ylabel("log₁₀(distance to nearest TE + 1) [bp]")
    ax.set_title(f"TE proximity by gene category\n"
                 f"(all species pooled, n_genes={len(df_gene):,})")

    # Panel B: per-species median distance — orphan vs other
    ax2 = axes[1]
    pivot = df_sum.pivot_table(index="species", columns="category",
                               values="median_dist", aggfunc="first")
    if "orphan" in pivot.columns and "other" in pivot.columns:
        common = pivot.dropna(subset=["orphan", "other"])
        ax2.scatter(common["other"], common["orphan"],
                    alpha=0.7, s=30, color="#d62728", label="species")
        lim = max(common[["orphan", "other"]].max())
        ax2.plot([0, lim], [0, lim], "k--", linewidth=0.8)
        ax2.set_xlabel("Median dist. non-orphan genes (bp)")
        ax2.set_ylabel("Median dist. orphan genes (bp)")
        ax2.set_title("Orphan vs non-orphan TE proximity\n(per species median)")
    else:
        ax2.text(0.5, 0.5, "Insufficient data for scatter", ha="center",
                 transform=ax2.transAxes)

    plt.tight_layout()
    plt.savefig(OUT_FIG, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"[OUT] Figure → {OUT_FIG}")
    print("[DONE]")


if __name__ == "__main__":
    main()
