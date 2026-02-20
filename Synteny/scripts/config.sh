#!/usr/bin/env bash
# =============================================================================
# config.sh — paths and clade definitions for the EDF-TE synteny analysis
# Source this file in every pipeline step: source config.sh
# =============================================================================

# ---------------------------------------------------------------------------
# Data directories (read-only — not in the repo)
# ---------------------------------------------------------------------------
export GFF_DIR="/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff"
export GENOME_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
export ORTHOLOGUES_DIR="/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthologues"
export ORTHOGROUPS_DIR="/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups"

# ---------------------------------------------------------------------------
# Output directories (generated locally, gitignored)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SYNTENY_DIR="${SCRIPT_DIR}/.."
export OUT_DIR="${SYNTENY_DIR}/data"
export FIG_DIR="${SYNTENY_DIR}/figures"
export BED_DIR="${OUT_DIR}/beds"
export ORTHOLOG_PAIRS_DIR="${OUT_DIR}/ortholog_pairs"
export COLLINEAR_DIR="${OUT_DIR}/collinear"
export SHUFFLE_DIR="${OUT_DIR}/shuffle"
export RATIO_DIR="${OUT_DIR}/ratios"

mkdir -p "${BED_DIR}" "${ORTHOLOG_PAIRS_DIR}" "${COLLINEAR_DIR}" \
         "${SHUFFLE_DIR}" "${RATIO_DIR}" "${FIG_DIR}"

# ---------------------------------------------------------------------------
# Analysis parameters
# ---------------------------------------------------------------------------
export N_SHUFFLES=100          # bedtools shuffle replicates
export MIN_GENES_PER_SCAFFOLD=3  # minimum genes on a scaffold to be included
export SEED_BASE=42            # base random seed; replicate i uses SEED_BASE+i

# ---------------------------------------------------------------------------
# Clade definitions  (associative array: CODE -> "space-separated species prefixes")
# Species prefix = {CODE}_{Name} as used in GFF/genome filenames
# ---------------------------------------------------------------------------
declare -A CLADES

# AP — Aphelida (2 spp.)
CLADES[AP]="AP_Amoocci AP_Amoprot"

# AS — Ascomycota (15 spp.)
CLADES[AS]="AS_AspfumUD1 AS_Bcine1 AS_Ccin1 AS_ChetHAW225 AS_Cocim1 AS_Dekbr2 AS_Metac1 AS_Ncrass2 AS_Sacce1 AS_Schpo1 AS_Talathe1 AS_Tuber1 AS_Zymtr1"

# BA — Basidiomycota (14 spp.)
CLADES[BA]="BA_Amamu1 BA_Copci2 BA_Fomme1 BA_Lacbi2 BA_Monro1 BA_Pucgr1 BA_Schco1 BA_Serpul2 BA_Treme1 BA_Ustma2"

# BL — Blastocladiomycota (5 spp.)
CLADES[BL]="BL_Allma1 BL_Blabri1 BL_Catan2 BL_Parsed1"

# BS — Basidiobolomycota (2 spp.)
CLADES[BS]="BS_Basme BS_Basme92"

# CH — Chytridiomycota (37 spp.)
CLADES[CH]="CH_Batde5 CH_Caecom1 CH_Gobsp1 CH_Lobvas1 CH_Newsp1 CH_Pirfi3 CH_Podang1 CH_Rhirose CH_Spipu1 CH_Triarc1"

# CR — Cryptomycota / Rozellida (2 spp.)
CLADES[CR]="CR_Parsa1 CR_Rozal1"

# EN — Entomophthorales (5 spp.)
CLADES[EN]="EN_Conco1 EN_Conth1 EN_Entmai1 EN_Mascica EN_Zoorad1"

# GL — Glomeromycota (7 spp.)
CLADES[GL]="GL_Geopyr1 GL_Gigmar1 GL_Paroc1 GL_Rhice1 GL_Rhicl1 GL_RhiirA4"

# HO — Holozoa outgroups (20 spp.)
CLADES[HO]="HO_Acancaste HO_Callivi HO_Capowcza HO_Ciointe HO_Dicdisco HO_Drome HO_Homosap HO_Hydravulg HO_Lucicu HO_Luciser HO_Magisepte HO_Mneleidy HO_Monbrevi HO_Musdom HO_Naegrub HO_Pluxyl HO_Salpirose HO_Sphaearti HO_Thecatra HO_Trichadha"

# KX — Kickxellomycotina (13 spp.)
CLADES[KX]="KX_Coemoj1 KX_DimcrSC1 KX_Kicala1 KX_SmicuMNP KX_Smimuc2 KX_Zancul2"

# MI — Microsporidia (22 spp.)
CLADES[MI]="MI_Enccu1 MI_Encin1 MI_MicG_I_3 MI_Nosap1 MI_Traho1 MI_Vitco1"

# MO — Mortierellomycota (26 spp.)
CLADES[MO]="MO_Moramb1 MO_Morel2 MO_Morgam1 MO_Morcap1 MO_Disdec1 MO_Lobtra1"

# MU — Mucoromycota / Mucorales (64 spp.)
CLADES[MU]="MU_Mucci3 MU_Mucmuc1 MU_Phybl4 MU_Rhimi59 MU_Rhipu1 MU_Rhisto1 MU_Zycmex1"

# NU — Nucleariida (2 spp.)
CLADES[NU]="NU_Fontalba NU_Parvatl"

# OL — Olpidium (1 sp. — degenerate, skipped in pairwise analysis)
CLADES[OL]="OL_Olpbor1"

# ZO — Zoopagomycota (4 spp.)
CLADES[ZO]="ZO_Pipcy3 ZO_Synfus1 ZO_Synplu1 ZO_Thasp1"

export CLADES

# Ordered list of all clade codes (for looping)
export ALL_CLADES=(AP AS BA BL BS CH CR EN GL HO KX MI MO MU NU OL ZO)

echo "[config] EDF-TE Synteny configuration loaded."
echo "         GFF_DIR   : ${GFF_DIR}"
echo "         OUT_DIR   : ${OUT_DIR}"
