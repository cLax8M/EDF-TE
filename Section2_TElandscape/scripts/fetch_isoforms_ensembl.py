#!/usr/bin/env python3
"""
fetch_isoforms_ensembl.py
=========================
Downloads GFF3 annotations from EnsemblFungi (via FTP) and FungiDB for
species in our dataset that currently have 0% multi-isoform genes, computes
isoform statistics, and writes updated *_ALL_isoform_counts.txt files.

Only uses curated, high-confidence species↔database mappings.
Strain mismatches are noted: isoform fractions are species-level properties
and are transferable across strains of the same species.

Usage:
    python fetch_isoforms_ensembl.py [--dry-run] [--force]

Requirements: none beyond stdlib + urllib (no pandas needed)
"""

import sys, os, re, gzip, urllib.request, collections, argparse
from pathlib import Path

# ─── Config ──────────────────────────────────────────────────────────────────
GENFEAT_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures_V3")
LOG_FILE    = GENFEAT_DIR / "_fetch_isoforms.log"

# EnsemblFungi FTP base (current release symlink)
EF_FTP = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/current/gff3"

# FungiDB FTP base
FUNGIDB_FTP = "https://fungidb.org/common/downloads/Current_Release"

# ─── Curated species mapping ─────────────────────────────────────────────────
# Format:
#   taxon_id : {
#     "db"      : "ensemblfungi" | "fungidb",
#     "path"    : FTP path suffix (species dir for ensembl, org abbrev for fungidb),
#     "species" : common/publication name for logging,
#     "note"    : strain note if applicable,
#   }
#
# Only species where:
#  a) public annotation is clearly the same species (strains may differ)
#  b) the database has RNA-seq-supported isoform annotations
#  c) we are confident the gene model set is comparable in quality

SPECIES_MAP = {
    # ── Ascomycota ────────────────────────────────────────────────────────────
    "AS_Ncrass2": {
        "db": "ensemblfungi",
        "path": "neurospora_crassa",
        "species": "Neurospora crassa OR74A",
        "note": "JGI v2 vs EnsemblFungi OR74A — same strain",
    },
    "AS_Sacce1": {
        "db": "ensemblfungi",
        "path": "saccharomyces_cerevisiae",
        "species": "Saccharomyces cerevisiae S288C",
        "note": "JGI Sacce1 is S288C — same strain",
    },
    "AS_Schpo1": {
        "db": "ensemblfungi",
        "path": "schizosaccharomyces_pombe",
        "species": "Schizosaccharomyces pombe 972h-",
        "note": "JGI Schpo1 is 972h — same strain",
    },
    "AS_Bcine1": {
        "db": "ensemblfungi",
        "path": "botrytis_cinerea",
        "species": "Botrytis cinerea B05.10",
        "note": "JGI Bcine1 is B05.10 — same strain",
    },
    "AS_Cocim1": {
        "db": "ensemblfungi",
        "path": "fungi_ascomycota1_collection/coccidioides_immitis_rs_gca_000149335",
        "species": "Coccidioides immitis RS",
        "note": "JGI Cocim1 is RS — same strain (in ascomycota1 collection)",
    },
    "AS_YarliYB419": {
        "db": "ensemblfungi",
        "path": "yarrowia_lipolytica",
        "species": "Yarrowia lipolytica CLIB89",
        "note": "JGI strain YB419 differs from EnsemblFungi CLIB89; isoform fraction transferable at species level",
    },
    "AS_AspfumUD1": {
        "db": "ensemblfungi",
        "path": "aspergillus_fumigatus",
        "species": "Aspergillus fumigatus Af293",
        "note": "JGI strain UD1 differs from EnsemblFungi Af293; both well-annotated clinical isolates",
    },
    "AS_Metac1": {
        "db": "ensemblfungi",
        "path": "fungi_ascomycota1_collection/metarhizium_acridum_cqma_102_gca_000187405",
        "species": "Metarhizium acridum CQMa102",
        "note": "Same strain (in ascomycota1 collection)",
    },
    "AS_Dekbr2": {
        "db": "ensemblfungi",
        "path": "fungi_ascomycota1_collection/brettanomyces_bruxellensis_awri1499_gca_000259595",
        "species": "Brettanomyces (Dekkera) bruxellensis AWRI1499",
        "note": "Dekkera=Brettanomyces anamorph; JGI Dekbr2 is CBS 2499, EnsemblFungi AWRI1499 — different strains, same species",
    },
    # ── Basidiomycota ─────────────────────────────────────────────────────────
    "BA_Ustma2": {
        "db": "ensemblfungi",
        "path": "ustilago_maydis",
        "species": "Ustilago maydis 521",
        "note": "Same strain as JGI",
    },
    "BA_Lacbi2": {
        "db": "ensemblfungi",
        "path": "fungi_basidiomycota1_collection/laccaria_bicolor_s238n_h82_gca_000143565",
        "species": "Laccaria bicolor S238N-H82",
        "note": "Same strain as JGI v2 (in basidiomycota1 collection)",
    },
    "BA_Rhoto": {
        "db": "ensemblfungi",
        "path": "fungi_basidiomycota1_collection/rhodotorula_toruloides_np11_gca_000320785",
        "species": "Rhodotorula toruloides NP11",
        "note": "JGI Rhoto is NP11 (in basidiomycota1 collection)",
    },
    "BA_Treme1": {
        "db": "ensemblfungi",
        "path": "fungi_basidiomycota1_collection/tremella_mesenterica_dsm_1558_gca_000271645",
        "species": "Tremella mesenterica DSM 1558",
        "note": "JGI Treme1 is DSM 1558 (in basidiomycota1 collection)",
    },
    # ── Microsporidia ─────────────────────────────────────────────────────────
    "MI_Enccu1": {
        "db": "ensemblfungi",
        "path": "fungi_microsporidia1_collection/encephalitozoon_cuniculi_gb_m1_gca_000091225",
        "species": "Encephalitozoon cuniculi GB-M1",
        "note": "Same strain (in microsporidia1 collection)",
    },
    "MI_Enche1": {
        "db": "ensemblfungi",
        "path": "fungi_microsporidia1_collection/encephalitozoon_hellem_atcc_50504_gca_000277815",
        "species": "Encephalitozoon hellem ATCC 50504",
        "note": "Same strain (in microsporidia1 collection)",
    },
    "MI_Encin1": {
        "db": "ensemblfungi",
        "path": "fungi_microsporidia1_collection/encephalitozoon_intestinalis_atcc_50506_gca_000146465",
        "species": "Encephalitozoon intestinalis ATCC 50506",
        "note": "Same strain (in microsporidia1 collection)",
    },
    "MI_Encro1": {
        "db": "ensemblfungi",
        "path": "fungi_microsporidia1_collection/encephalitozoon_romaleae_sj_2008_gca_000280035",
        "species": "Encephalitozoon romaleae SJ-2008",
        "note": "Same strain (in microsporidia1 collection)",
    },
    "MI_Edhae1": {
        "db": "ensemblfungi",
        "path": "fungi_microsporidia1_collection/edhazardia_aedis_usnm_41457_gca_000230595",
        "species": "Edhazardia aedis USNM 41457",
        "note": "Same strain (in microsporidia1 collection)",
    },
    "MI_AmpWSBS2006": {
        "db": "ensemblfungi",
        "path": "fungi_microsporidia1_collection/amphiamblys_sp_wsbs2006_gca_001875675",
        "species": "Amphiamblys sp. WSBS2006",
        "note": "Same strain (in microsporidia1 collection)",
    },
    # ── Chytridiomycota ───────────────────────────────────────────────────────
    "CH_Batde5": {
        "db": "ensemblfungi",
        "path": "fungi_chytridiomycota1_collection/batrachochytrium_dendrobatidis_jel423_gca_000149865",
        "species": "Batrachochytrium dendrobatidis JEL423",
        "note": "JGI Batde5 is JEL423 — same strain (in chytridiomycota1 collection)",
    },
}


# ─── Helpers ─────────────────────────────────────────────────────────────────
def log(msg, file=None):
    print(msg)
    if file:
        print(msg, file=file)


def find_ensemblfungi_gff(species_path):
    """Return the GFF3 URL for a given EnsemblFungi species path component."""
    listing_url = f"{EF_FTP}/{species_path}/"
    try:
        with urllib.request.urlopen(listing_url, timeout=30) as r:
            html = r.read().decode()
    except Exception as e:
        return None, str(e)

    # Find .gff3.gz file (not chromosome, not abinitio)
    matches = re.findall(r'href="([^"]+\.gff3\.gz)"', html)
    # Prefer the "top-level" or main annotation file
    candidates = [m for m in matches if "abinitio" not in m and "chromosome" not in m.lower()]
    if not candidates:
        candidates = matches
    if not candidates:
        return None, "no .gff3.gz found"
    # Pick smallest suffix (usually the primary annotation)
    candidates.sort(key=len)
    fname = candidates[0]
    if fname.startswith("http"):
        return fname, None
    return f"{EF_FTP}/{species_path}/{fname}", None


def compute_isoforms_from_gff3_gz(gff_path):
    """
    Read a (possibly gzipped) GFF3 file and return (gene_mrna_dict, error_str).
    gene_mrna_dict: {gene_id: n_mRNA}
    """
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    gene_mrna = collections.defaultdict(int)
    try:
        with opener(gff_path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.strip().split("\t")
                if len(p) < 9:
                    continue
                feat = p[2]
                if feat in ("mRNA", "transcript", "mrna"):
                    m = re.search(r"Parent=([^;,\s]+)", p[8])
                    if m:
                        gene_mrna[m.group(1)] += 1
    except Exception as e:
        return None, str(e)
    return gene_mrna, None


def write_isoform_counts(out_path, gene_mrna):
    with open(out_path, "w") as f:
        f.write("Gene\tIsoform_count\n")
        for gene, cnt in sorted(gene_mrna.items()):
            f.write(f"{gene}\t{cnt}\n")


def download_file(url, dest):
    """Download url to dest, return error string or None."""
    try:
        urllib.request.urlretrieve(url, dest)
        return None
    except Exception as e:
        return str(e)


# ─── Main ────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="Only resolve URLs, don't download")
    ap.add_argument("--force", action="store_true", help="Re-download even if isoform file already non-zero")
    args = ap.parse_args()

    GENFEAT_DIR.mkdir(parents=True, exist_ok=True)
    tmp_dir = GENFEAT_DIR / "_gff3_cache"
    tmp_dir.mkdir(exist_ok=True)

    results = {}

    with open(LOG_FILE, "w") as logf:
        log(f"=== fetch_isoforms_ensembl.py ===", logf)
        log(f"Dry-run: {args.dry_run}  Force: {args.force}", logf)
        log("", logf)

        for taxon, info in SPECIES_MAP.items():
            sp       = info["species"]
            note     = info["note"]
            out_file = GENFEAT_DIR / f"{taxon}_ALL_isoform_counts.txt"

            # Skip if already has non-zero data (unless --force)
            if not args.force and out_file.exists():
                with open(out_file) as f:
                    next(f)
                    counts = [int(l.strip().split("\t")[1]) for l in f]
                pct = 100 * sum(1 for c in counts if c >= 2) / len(counts) if counts else 0
                if pct > 0:
                    log(f"[SKIP] {taxon} — already has {pct:.1f}% multi-isoform", logf)
                    results[taxon] = ("skipped", pct)
                    continue

            log(f"\n[PROCESS] {taxon} → {sp}", logf)
            log(f"  Note: {note}", logf)

            # ── Resolve GFF URL ───────────────────────────────────────────────
            if info["db"] == "ensemblfungi":
                gff_url, err = find_ensemblfungi_gff(info["path"])
                if err:
                    log(f"  ERROR resolving URL: {err}", logf)
                    results[taxon] = ("url_error", err)
                    continue
            else:
                log(f"  SKIP: FungiDB not yet implemented", logf)
                results[taxon] = ("skipped", "fungidb not implemented")
                continue

            log(f"  URL: {gff_url}", logf)

            if args.dry_run:
                results[taxon] = ("dry-run", gff_url)
                continue

            # ── Download ─────────────────────────────────────────────────────
            gff_local = tmp_dir / f"{taxon}.gff3.gz"
            if not gff_local.exists():
                log(f"  Downloading...", logf)
                err = download_file(gff_url, gff_local)
                if err:
                    log(f"  ERROR downloading: {err}", logf)
                    results[taxon] = ("download_error", err)
                    continue
                log(f"  Downloaded: {gff_local.stat().st_size/1e6:.1f} MB", logf)
            else:
                log(f"  Using cached: {gff_local}", logf)

            # ── Parse ─────────────────────────────────────────────────────────
            gene_mrna, err = compute_isoforms_from_gff3_gz(gff_local)
            if err:
                log(f"  ERROR parsing: {err}", logf)
                results[taxon] = ("parse_error", err)
                continue

            if not gene_mrna:
                log(f"  WARNING: no mRNA/transcript features found", logf)
                results[taxon] = ("empty", "no mRNA features")
                continue

            n_genes = len(gene_mrna)
            n_multi = sum(1 for v in gene_mrna.values() if v >= 2)
            pct_multi = 100 * n_multi / n_genes
            avg_iso = sum(gene_mrna.values()) / n_genes

            log(f"  Genes: {n_genes}  multi-isoform: {n_multi} ({pct_multi:.1f}%)  avg: {avg_iso:.3f}", logf)

            # Safety check: reject if the EnsemblFungi annotation is suspiciously
            # sparse compared to the existing file (< 30% of existing gene count).
            # This guards against partially-imported EnsemblFungi records.
            if out_file.exists():
                with open(out_file) as f:
                    existing_genes = sum(1 for _ in f) - 1  # subtract header
                if existing_genes > 0 and n_genes < 0.3 * existing_genes:
                    log(f"  WARN: downloaded {n_genes} genes vs existing {existing_genes} "
                        f"— likely incomplete EnsemblFungi annotation; skipping", logf)
                    results[taxon] = ("skipped_sparse", n_genes)
                    continue

            write_isoform_counts(out_file, gene_mrna)
            log(f"  Written: {out_file}", logf)
            results[taxon] = ("updated", pct_multi)

        # ── Summary ───────────────────────────────────────────────────────────
        log("\n\n=== SUMMARY ===", logf)
        updated  = [(t,v) for t,(s,v) in results.items() if s == "updated"]
        skipped  = [(t,v) for t,(s,v) in results.items() if s == "skipped"]
        errors   = [(t,v) for t,(s,v) in results.items() if "error" in s]

        log(f"\nUpdated ({len(updated)}):", logf)
        for t, v in sorted(updated):
            log(f"  {t:<30}  {v:.1f}%", logf)
        log(f"\nSkipped ({len(skipped)}):", logf)
        for t, v in skipped:
            log(f"  {t:<30}  {v}", logf)
        if errors:
            log(f"\nErrors ({len(errors)}):", logf)
            for t, v in errors:
                log(f"  {t:<30}  {v}", logf)

    print(f"\nLog written to {LOG_FILE}")
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
