#!/usr/bin/env python3
"""
02_extract_orthologs.py — Extract 1:1 ortholog pairs for all within-clade
                          species comparisons from OrthoFinder Orthologues/ tables.

For each ordered pair (SpA, SpB) within a clade, reads:
  {ORTHOLOGUES_DIR}/Orthologues_{SpA}.aa/{SpA}.aa__v__{SpB}.aa.tsv
  (OrthoFinder uses .aa suffix in directory/file names)

Outputs a TSV per pair:
  {ORTHOLOG_PAIRS_DIR}/{CLADE}/{SpA}__{SpB}.pairs.tsv
  Columns: orthogroup  gene_A  gene_B

Only rows where both columns contain a single gene (1:1) are retained.

Usage:
  python 02_extract_orthologs.py --clade CLADE [--clade CLADE2 ...]
  python 02_extract_orthologs.py --all

Author: EDF-TE project
Last updated: 2026-02-20
"""

import argparse
import os
import sys
from pathlib import Path
from itertools import combinations


# ---------------------------------------------------------------------------
# Clade definitions (mirrors config.sh)
# ---------------------------------------------------------------------------

CLADES = {
    "AP": ["AP_Amoocci", "AP_Amoprot"],
    "AS": ["AS_AspfumUD1", "AS_Bcine1", "AS_ChetHAW225", "AS_Cocim1", "AS_Dekbr2",
           "AS_Metac1", "AS_Ncrass2", "AS_Sacce1", "AS_Schpo1", "AS_Tapde1",
           "AS_Usnflo1", "AS_YarliYB419", "AS_Zymtr1"],
    "BA": ["BA_Amamu1", "BA_Aurde3", "BA_Fomme1", "BA_Hetan2", "BA_Lacbi2",
           "BA_Maglo", "BA_Mereb1", "BA_Mrafri1", "BA_PleosDSM11191", "BA_Rhisol1",
           "BA_Rhoto", "BA_SerlaS7", "BA_Treme1", "BA_Ustma2"],
    "BL": ["BL_Allma1", "BL_Allomarbu", "BL_Blabri1", "BL_Catan2", "BL_Parsed1"],
    "BS": ["BS_Basme", "BS_Basme92"],
    "CH": ["CH_Anasp1", "CH_Batde5", "CH_Batsa1", "CH_Blyhe1", "CH_Caecom1",
           "CH_Caupr1", "CH_Chyhya1", "CH_Chylag1", "CH_Clapol1", "CH_Clarep1",
           "CH_Dinkinn", "CH_Enthel1", "CH_Entlut1", "CH_Fimjon1", "CH_Gaesem1",
           "CH_Ganpr1", "CH_Germichi", "CH_Gervar1", "CH_Glopol1", "CH_Gorhay1",
           "CH_Hompol1", "CH_Hyacur1", "CH_Lobangu", "CH_Neocon1", "CH_Neolan1",
           "CH_Obemuc1", "CH_Orpsp1", "CH_Pecora1", "CH_Phybulla", "CH_Pirfi3",
           "CH_Polagg1", "CH_Powhir1", "CH_Rhihy1", "CH_Rhirose", "CH_Spipu1",
           "CH_Terramyc", "CH_Triarc1"],
    "CR": ["CR_Parsa1", "CR_Rozal1"],
    "EN": ["EN_Conco1", "EN_Conth1", "EN_Entmai1", "EN_Mascica", "EN_Zoorad1"],
    "GL": ["GL_Geopyr1", "GL_Gigmar1", "GL_Gigro1", "GL_Paroc1", "GL_Rhice1",
           "GL_Rhicl1", "GL_RhiirA4"],
    "HO": ["HO_Acancaste", "HO_Callivi", "HO_Capowcza", "HO_Ciointe",
           "HO_Dicdisco", "HO_Drome", "HO_Homosap", "HO_Hydravulg",
           "HO_Lucicu", "HO_Luciser", "HO_Magisepte_c", "HO_Mneleidy",
           "HO_Monbrevi", "HO_Musdom", "HO_Naegrub", "HO_Pluxyl",
           "HO_Salpirose", "HO_Sphaearti", "HO_Thecatra", "HO_Trichadha"],
    "KX": ["KX_Coemoj1", "KX_Coespi1", "KX_DimcrSC1", "KX_Dispsimple",
           "KX_Kicala1", "KX_Linpe1", "KX_Marpt1", "KX_Pincor1", "KX_Rambr1",
           "KX_SmicuMNP", "KX_Smimuc2", "KX_Spiroaspi", "KX_Zancul2"],
    "MI": ["MI_AmpWSBS2006", "MI_Antlo1", "MI_Edhae1", "MI_Enccu1", "MI_Enche1",
           "MI_Encin1", "MI_Encro1", "MI_Entbi1", "MI_Entca1", "MI_Heper1",
           "MI_MicG_I_3", "MI_Mitdap1", "MI_Nemdi1", "MI_Nempa1", "MI_Nosap1",
           "MI_Nosbo1", "MI_Ordco1", "MI_Pseneu1", "MI_Sprlo1", "MI_Traho1",
           "MI_Vavcu1", "MI_Vitco1"],
    "MO": ["MO_Disdec1", "MO_Disorn1", "MO_Entech1", "MO_HabiNA12553",
           "MO_Lobtra1", "MO_MoalAD071", "MO_MoamBC1065", "MO_MobeAD1040",
           "MO_MocaN2610", "MO_MochAD033", "MO_MoclAM1000", "MO_Mocy1230",
           "MO_MoelGBA34", "MO_MoepAD058", "MO_MoepN5512", "MO_MomuN6456",
           "MO_Moramb1", "MO_Morcap1", "MO_MoreMES2147", "MO_Morel2",
           "MO_Morgam1", "MO_Mormul1", "MO_Morsel1", "MO_Morwol1",
           "MO_MowoN66265", "MO_MozyN2592"],
    "MU": ["MU_Abscae1", "MU_Absrep1", "MU_Amyrou1", "MU_Apoossi", "MU_Bacci1",
           "MU_Baclam1", "MU_Benpoi1", "MU_Bifad1", "MU_Blatri", "MU_Calth1",
           "MU_Chabre1", "MU_Chlapade", "MU_Chocucu1", "MU_Cirumb1", "MU_Cokrec1",
           "MU_Cunech1", "MU_Dicele1", "MU_Dicrob1", "MU_Ellano2", "MU_Endsp1",
           "MU_Fenlin", "MU_Gilper1", "MU_Gonbut1", "MU_Halrad1", "MU_Helpul1",
           "MU_Hesve2", "MU_Jimfl", "MU_Jimlac1", "MU_Kircor1", "MU_Lenpar1",
           "MU_Lenzyc1", "MU_Liccor1", "MU_Licra1", "MU_Mucam1", "MU_Mucci3",
           "MU_Muccirc1", "MU_Mucend1", "MU_Muclan1", "MU_Mucmuc1", "MU_Mucpir1",
           "MU_Mucrac1", "MU_Mycafr1", "MU_Mycind1", "MU_Parpar1", "MU_Phaart1",
           "MU_Phybl4", "MU_Pilano1", "MU_Pilumb1", "MU_PnitS607", "MU_Radspe1",
           "MU_Rhimi59", "MU_Rhior3", "MU_Rhipu1", "MU_Rhisto1", "MU_Sakvas1",
           "MU_Spifus1", "MU_Spoumb1", "MU_Synrac1", "MU_Syzme1", "MU_Thaele1",
           "MU_Umban1", "MU_Umbra1", "MU_Zycmex1"],
    "NU": ["NU_Fontalba", "NU_Parvatl_c"],
    "OL": ["OL_Olpbor1"],
    "ZO": ["ZO_Pipcy3", "ZO_Synfus1", "ZO_Synplu1", "ZO_Thasp1"],
}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orthologues_dir", required=True,
                   help="Path to OrthoFinder Orthologues/ directory")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for .pairs.tsv files")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--clade", action="append", choices=list(CLADES.keys()),
                       help="Clade code(s) to process")
    group.add_argument("--all", action="store_true",
                       help="Process all clades")
    p.add_argument("--dry_run", action="store_true",
                   help="Print expected file paths without reading them")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------

def find_orthologues_file(orth_dir, sp_a, sp_b):
    """
    Locate the OrthoFinder pairwise ortholog TSV for (sp_a, sp_b).
    OrthoFinder creates files as:
      Orthologues_{SpA}.aa/{SpA}.aa__v__{SpB}.aa.tsv  (or reversed)
    Returns (path, reversed) where reversed=True if sp_b is the "source".
    """
    orth_dir = Path(orth_dir)

    # Try A vs B
    candidate = orth_dir / f"Orthologues_{sp_a}.aa" / f"{sp_a}.aa__v__{sp_b}.aa.tsv"
    if candidate.exists():
        return candidate, False

    # Try B vs A
    candidate_rev = orth_dir / f"Orthologues_{sp_b}.aa" / f"{sp_b}.aa__v__{sp_a}.aa.tsv"
    if candidate_rev.exists():
        return candidate_rev, True

    # Fallback: species names without .aa suffix
    candidate2 = orth_dir / f"Orthologues_{sp_a}" / f"{sp_a}__v__{sp_b}.tsv"
    if candidate2.exists():
        return candidate2, False

    candidate2_rev = orth_dir / f"Orthologues_{sp_b}" / f"{sp_b}__v__{sp_a}.tsv"
    if candidate2_rev.exists():
        return candidate2_rev, True

    # Last resort: glob for mixed .aa / no-.aa naming (e.g. BA_Treme1, CH_Batde5)
    for sp_x, sp_y, is_rev in [(sp_a, sp_b, False), (sp_b, sp_a, True)]:
        for dir_suffix in [".aa", ""]:
            d = orth_dir / f"Orthologues_{sp_x}{dir_suffix}"
            if d.exists():
                matches = list(d.glob(f"{sp_x}*__v__{sp_y}*.tsv"))
                if matches:
                    return matches[0], is_rev

    return None, None


def extract_one_to_one(tsv_path, reversed_order=False):
    """
    Read an OrthoFinder __v__.tsv file and return only 1:1 ortholog rows.
    OrthoFinder format:
      Orthogroup  SpA_gene(s)  SpB_gene(s)
    Multi-copy entries have comma-separated gene IDs — we skip those.

    Returns list of (orthogroup, gene_a, gene_b).
    If reversed_order=True, columns are (SpB, SpA) → swap back.
    """
    pairs = []
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if len(header) < 3:
            return pairs

        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            orthogroup = cols[0]
            genes_col1 = cols[1].strip()
            genes_col2 = cols[2].strip()

            # Skip multi-copy (comma-separated)
            if "," in genes_col1 or "," in genes_col2:
                continue
            if not genes_col1 or not genes_col2:
                continue

            if reversed_order:
                gene_a, gene_b = genes_col2, genes_col1
            else:
                gene_a, gene_b = genes_col1, genes_col2

            pairs.append((orthogroup, gene_a, gene_b))

    return pairs


def process_clade(clade_code, species_list, orth_dir, out_dir, dry_run=False):
    """Process all pairs within a clade."""
    if len(species_list) < 2:
        print(f"  [{clade_code}] Only {len(species_list)} species — skipping (need ≥2).")
        return

    clade_out = Path(out_dir) / clade_code
    clade_out.mkdir(parents=True, exist_ok=True)

    pairs_list = list(combinations(sorted(species_list), 2))
    print(f"  [{clade_code}] {len(species_list)} species → {len(pairs_list)} pairs")

    found = 0
    missing = 0
    total_pairs = 0

    for sp_a, sp_b in pairs_list:
        out_file = clade_out / f"{sp_a}__{sp_b}.pairs.tsv"

        tsv_path, reversed_order = find_orthologues_file(orth_dir, sp_a, sp_b)

        if dry_run:
            status = "OK" if tsv_path else "MISSING"
            print(f"    {sp_a} × {sp_b}  [{status}]")
            continue

        if tsv_path is None:
            print(f"    WARNING: no OrthoFinder file for {sp_a} × {sp_b} — skipping")
            missing += 1
            continue

        pairs = extract_one_to_one(tsv_path, reversed_order=reversed_order)

        with open(out_file, "w") as fh:
            fh.write("orthogroup\tgene_A\tgene_B\n")
            for og, ga, gb in pairs:
                fh.write(f"{og}\t{ga}\t{gb}\n")

        found += 1
        total_pairs += len(pairs)

    if not dry_run:
        print(f"    → {found} pair files written, {missing} missing, "
              f"{total_pairs} total 1:1 ortholog pairs")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if args.all:
        target_clades = list(CLADES.keys())
    else:
        target_clades = args.clade

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for clade_code in target_clades:
        species_list = CLADES.get(clade_code, [])
        process_clade(clade_code, species_list, args.orthologues_dir,
                      out_dir, dry_run=args.dry_run)

    print("Step 02 complete.")


if __name__ == "__main__":
    main()
