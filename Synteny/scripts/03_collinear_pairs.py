#!/usr/bin/env python3
"""
03_collinear_pairs.py — Identify collinear (syntenic) ortholog pairs for all
                        within-clade species comparisons.

Definition (following Sebé-Pedrós et al. 2013 / Putnam et al. 2007):
  Two ortholog pairs (A1-B1) and (A2-B2) are collinear if A1 and A2 are
  ADJACENT on the same scaffold in species A AND B1 and B2 are ADJACENT
  on the same scaffold in species B, irrespective of gene orientation.

Output per pair:
  {COLLINEAR_DIR}/{CLADE}/{SpA}__{SpB}.collinear.tsv
  Columns: scaffold_A  gene_A1  gene_A2  scaffold_B  gene_B1  gene_B2

Summary TSV (all pairs):
  {COLLINEAR_DIR}/{CLADE}/{CLADE}_collinear_summary.tsv
  Columns: species_A  species_B  c_ij  N_ij

  c_ij = number of collinear ortholog pairs observed
  N_ij = number of ortholog pairs that could potentially be adjacent
         (= sum over scaffolds of (n_orthologs_on_scaffold - 1) for each species,
          averaged between SpA and SpB)

Usage:
  python 03_collinear_pairs.py --clade CLADE [--clade ...] | --all
                               --bed_dir BED_DIR
                               --pairs_dir ORTHOLOG_PAIRS_DIR
                               --out_dir COLLINEAR_DIR

Author: EDF-TE project
Last updated: 2026-02-20
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from itertools import combinations


# Clade definitions (same as 02_extract_orthologs.py)
CLADES = {
    "AP": ["AP_Amoocci", "AP_Amoprot"],
    "AS": ["AS_AspfumUD1", "AS_Bcine1", "AS_Ccin1", "AS_ChetHAW225", "AS_Cocim1",
           "AS_Dekbr2", "AS_Metac1", "AS_Ncrass2", "AS_Sacce1", "AS_Schpo1",
           "AS_Talathe1", "AS_Tuber1", "AS_Zymtr1"],
    "BA": ["BA_Amamu1", "BA_Copci2", "BA_Fomme1", "BA_Lacbi2", "BA_Monro1",
           "BA_Pucgr1", "BA_Schco1", "BA_Serpul2", "BA_Treme1", "BA_Ustma2"],
    "BL": ["BL_Allma1", "BL_Blabri1", "BL_Catan2", "BL_Parsed1"],
    "BS": ["BS_Basme", "BS_Basme92"],
    "CH": ["CH_Batde5", "CH_Caecom1", "CH_Gobsp1", "CH_Lobvas1", "CH_Newsp1",
           "CH_Pirfi3", "CH_Podang1", "CH_Rhirose", "CH_Spipu1", "CH_Triarc1"],
    "CR": ["CR_Parsa1", "CR_Rozal1"],
    "EN": ["EN_Conco1", "EN_Conth1", "EN_Entmai1", "EN_Mascica", "EN_Zoorad1"],
    "GL": ["GL_Geopyr1", "GL_Gigmar1", "GL_Paroc1", "GL_Rhice1", "GL_Rhicl1",
           "GL_RhiirA4"],
    "HO": ["HO_Acancaste", "HO_Callivi", "HO_Capowcza", "HO_Ciointe",
           "HO_Dicdisco", "HO_Drome", "HO_Homosap", "HO_Hydravulg",
           "HO_Lucicu", "HO_Luciser", "HO_Magisepte", "HO_Mneleidy",
           "HO_Monbrevi", "HO_Musdom", "HO_Naegrub", "HO_Pluxyl",
           "HO_Salpirose", "HO_Sphaearti", "HO_Thecatra", "HO_Trichadha"],
    "KX": ["KX_Coemoj1", "KX_DimcrSC1", "KX_Kicala1", "KX_SmicuMNP",
           "KX_Smimuc2", "KX_Zancul2"],
    "MI": ["MI_Enccu1", "MI_Encin1", "MI_MicG_I_3", "MI_Nosap1",
           "MI_Traho1", "MI_Vitco1"],
    "MO": ["MO_Moramb1", "MO_Morel2", "MO_Morgam1", "MO_Morcap1",
           "MO_Disdec1", "MO_Lobtra1"],
    "MU": ["MU_Mucci3", "MU_Mucmuc1", "MU_Phybl4", "MU_Rhimi59",
           "MU_Rhipu1", "MU_Rhisto1", "MU_Zycmex1"],
    "NU": ["NU_Fontalba", "NU_Parvatl"],
    "OL": ["OL_Olpbor1"],
    "ZO": ["ZO_Pipcy3", "ZO_Synfus1", "ZO_Synplu1", "ZO_Thasp1"],
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed_dir",   required=True, help="Directory with *.bed files (from step 01)")
    p.add_argument("--pairs_dir", required=True, help="Directory with ortholog pair TSVs (from step 02)")
    p.add_argument("--out_dir",   required=True, help="Output directory for collinear results")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--clade", nargs="+", choices=list(CLADES.keys()))
    group.add_argument("--all", action="store_true")
    return p.parse_args()


def load_bed(bed_path):
    """
    Load a BED file → dict: gene_id -> (chrom, start, rank_on_chrom)
    Rank is the 0-based position of the gene on its scaffold (sorted by start).
    Also returns: chrom -> [gene_id_in_order]
    """
    chrom_genes = defaultdict(list)  # chrom -> [(start, gene_id)]

    with open(bed_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            chrom, start, end, gene_id = cols[0], int(cols[1]), int(cols[2]), cols[3]
            chrom_genes[chrom].append((start, gene_id))

    # Sort each chrom by start position
    gene_info = {}   # gene_id -> (chrom, rank)
    chrom_order = {} # chrom -> [gene_id in genome order]

    for chrom, entries in chrom_genes.items():
        entries_sorted = sorted(entries, key=lambda x: x[0])
        chrom_order[chrom] = [gene_id for _, gene_id in entries_sorted]
        for rank, (start, gene_id) in enumerate(entries_sorted):
            gene_info[gene_id] = (chrom, rank)

    return gene_info, chrom_order


def load_ortholog_pairs(pairs_path):
    """Load ortholog pairs TSV → list of (gene_a, gene_b)."""
    pairs = []
    with open(pairs_path) as fh:
        next(fh)  # skip header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                pairs.append((cols[1], cols[2]))
    return pairs


def find_collinear_pairs(gene_info_a, chrom_order_a,
                         gene_info_b, chrom_order_b,
                         ortholog_pairs):
    """
    Identify collinear ortholog pairs.

    Two ortholog pairs (a1,b1) and (a2,b2) are collinear if:
      - a1 and a2 are adjacent on scaffold X in species A  (rank difference = 1)
      - b1 and b2 are adjacent on scaffold Y in species B  (rank difference = 1)

    Strategy:
      1. Build a set of ortholog genes in each species from the pair list.
      2. For each scaffold in species A, collect orthologs in that scaffold's order.
      3. For each consecutive pair of orthologs on the scaffold [oa1, oa2]:
         check if their partners [ob1, ob2] are also consecutive on the same
         scaffold in species B.

    Returns:
      collinear_pairs: list of (chrom_a, gene_a1, gene_a2, chrom_b, gene_b1, gene_b2)
      c_ij: count of collinear pairs
      N_ij: number of comparable ortholog-adjacent pairs (potential syntenic pairs)
    """
    # Build lookup: gene_in_A -> gene_in_B
    a_to_b = {}
    b_to_a = {}
    for gene_a, gene_b in ortholog_pairs:
        # Only include genes present in BED files
        if gene_a in gene_info_a and gene_b in gene_info_b:
            a_to_b[gene_a] = gene_b
            b_to_a[gene_b] = gene_a

    # For each scaffold in species A, get only the orthologs in order
    chrom_orthologs_a = defaultdict(list)
    for gene_a in a_to_b:
        if gene_a in gene_info_a:
            chrom, rank = gene_info_a[gene_a]
            chrom_orthologs_a[chrom].append((rank, gene_a))

    # Sort ortholog lists by rank
    for chrom in chrom_orthologs_a:
        chrom_orthologs_a[chrom].sort(key=lambda x: x[0])

    # Also build the adjacency map for species B orthologs on each scaffold
    # gene_b -> (chrom_b, rank_among_orthologs_on_that_scaffold)
    # We need consecutive ortholog adjacency in B, not just genomic adjacency
    chrom_orthologs_b = defaultdict(list)
    for gene_b in b_to_a:
        if gene_b in gene_info_b:
            chrom, rank = gene_info_b[gene_b]
            chrom_orthologs_b[chrom].append((rank, gene_b))
    for chrom in chrom_orthologs_b:
        chrom_orthologs_b[chrom].sort(key=lambda x: x[0])

    # Build lookup: gene_b -> (chrom_b, rank_among_orthologs)
    gene_b_ortholog_rank = {}
    for chrom, entries in chrom_orthologs_b.items():
        for orth_rank, (genomic_rank, gene_b) in enumerate(entries):
            gene_b_ortholog_rank[gene_b] = (chrom, orth_rank)

    collinear_pairs = []
    N_ij = 0  # potential adjacencies

    for chrom_a, entries in chrom_orthologs_a.items():
        if len(entries) < 2:
            continue
        N_ij += len(entries) - 1  # n adjacent pairs possible on this scaffold in A

        for i in range(len(entries) - 1):
            _, gene_a1 = entries[i]
            _, gene_a2 = entries[i + 1]

            gene_b1 = a_to_b[gene_a1]
            gene_b2 = a_to_b[gene_a2]

            if gene_b1 not in gene_b_ortholog_rank or gene_b2 not in gene_b_ortholog_rank:
                continue

            chrom_b1, orth_rank_b1 = gene_b_ortholog_rank[gene_b1]
            chrom_b2, orth_rank_b2 = gene_b_ortholog_rank[gene_b2]

            # Same scaffold in B AND adjacent (rank difference = 1, direction-agnostic)
            if chrom_b1 == chrom_b2 and abs(orth_rank_b1 - orth_rank_b2) == 1:
                collinear_pairs.append((chrom_a, gene_a1, gene_a2,
                                        chrom_b1, gene_b1, gene_b2))

    return collinear_pairs, len(collinear_pairs), N_ij


def process_clade(clade_code, bed_dir, pairs_dir, out_dir):
    species_list = CLADES.get(clade_code, [])
    if len(species_list) < 2:
        print(f"  [{clade_code}] Only {len(species_list)} species — skipping.")
        return []

    clade_out = Path(out_dir) / clade_code
    clade_out.mkdir(parents=True, exist_ok=True)
    pairs_clade_dir = Path(pairs_dir) / clade_code
    bed_dir = Path(bed_dir)

    summary_rows = []
    all_pairs = list(combinations(sorted(species_list), 2))
    print(f"  [{clade_code}] Processing {len(all_pairs)} species pairs ...")

    # Cache BED data to avoid re-reading
    bed_cache = {}
    def get_bed(sp):
        if sp not in bed_cache:
            bed_path = bed_dir / f"{sp}.bed"
            if not bed_path.exists():
                print(f"    WARNING: BED not found for {sp} — skipping pairs involving this species")
                bed_cache[sp] = (None, None)
            else:
                bed_cache[sp] = load_bed(bed_path)
        return bed_cache[sp]

    for sp_a, sp_b in all_pairs:
        gene_info_a, chrom_order_a = get_bed(sp_a)
        gene_info_b, chrom_order_b = get_bed(sp_b)
        if gene_info_a is None or gene_info_b is None:
            continue

        pairs_file = pairs_clade_dir / f"{sp_a}__{sp_b}.pairs.tsv"
        if not pairs_file.exists():
            print(f"    WARNING: pairs file missing for {sp_a} × {sp_b}")
            continue

        ortholog_pairs = load_ortholog_pairs(pairs_file)
        if not ortholog_pairs:
            continue

        collinear, c_ij, N_ij = find_collinear_pairs(
            gene_info_a, chrom_order_a,
            gene_info_b, chrom_order_b,
            ortholog_pairs
        )

        # Write collinear pairs file
        out_file = clade_out / f"{sp_a}__{sp_b}.collinear.tsv"
        with open(out_file, "w") as fh:
            fh.write("scaffold_A\tgene_A1\tgene_A2\tscaffold_B\tgene_B1\tgene_B2\n")
            for row in collinear:
                fh.write("\t".join(str(x) for x in row) + "\n")

        summary_rows.append({
            "species_A": sp_a, "species_B": sp_b,
            "n_orthologs": len(ortholog_pairs),
            "c_ij": c_ij, "N_ij": N_ij
        })

    # Write clade summary
    summary_file = clade_out / f"{clade_code}_collinear_summary.tsv"
    with open(summary_file, "w") as fh:
        fh.write("species_A\tspecies_B\tn_orthologs\tc_ij\tN_ij\n")
        for row in summary_rows:
            fh.write(f"{row['species_A']}\t{row['species_B']}\t"
                     f"{row['n_orthologs']}\t{row['c_ij']}\t{row['N_ij']}\n")

    print(f"    → summary written to {summary_file}")
    return summary_rows


def main():
    args = parse_args()
    target_clades = list(CLADES.keys()) if args.all else args.clade
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for clade_code in target_clades:
        process_clade(clade_code, args.bed_dir, args.pairs_dir, out_dir)

    print("Step 03 complete.")


if __name__ == "__main__":
    main()
