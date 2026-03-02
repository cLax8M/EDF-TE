#!/usr/bin/env bash
# =============================================================================
# config.sh — paths and clade definitions for the EDF-TE synteny analysis
# Source this file in every pipeline step: source config.sh
# =============================================================================

# ---------------------------------------------------------------------------
# Data directories (read-only — not in the repo)
# ---------------------------------------------------------------------------
export GFF_DIR="/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff"
export FAA_DIR="/legserv/Temp/Carlos/EDF-TE/data/proteins_annot_faa"
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

# AS — Ascomycota (13 spp.)
CLADES[AS]="AS_AspfumUD1 AS_Bcine1 AS_ChetHAW225 AS_Cocim1 AS_Dekbr2 AS_Metac1 AS_Ncrass2 AS_Sacce1 AS_Schpo1 AS_Tapde1 AS_Usnflo1 AS_YarliYB419 AS_Zymtr1"

# BA — Basidiomycota (14 spp.)
CLADES[BA]="BA_Amamu1 BA_Aurde3 BA_Fomme1 BA_Hetan2 BA_Lacbi2 BA_Maglo BA_Mereb1 BA_Mrafri1 BA_PleosDSM11191 BA_Rhisol1 BA_Rhoto BA_SerlaS7 BA_Treme1 BA_Ustma2"

# BL — Blastocladiomycota (5 spp.)
CLADES[BL]="BL_Allma1 BL_Allomarbu BL_Blabri1 BL_Catan2 BL_Parsed1"

# BS — Basidiobolomycota (2 spp.)
CLADES[BS]="BS_Basme BS_Basme92"

# CH — Chytridiomycota (37 spp.)
CLADES[CH]="CH_Anasp1 CH_Batde5 CH_Batsa1 CH_Blyhe1 CH_Caecom1 CH_Caupr1 CH_Chyhya1 CH_Chylag1 CH_Clapol1 CH_Clarep1 CH_Dinkinn CH_Enthel1 CH_Entlut1 CH_Fimjon1 CH_Gaesem1 CH_Ganpr1 CH_Germichi CH_Gervar1 CH_Glopol1 CH_Gorhay1 CH_Hompol1 CH_Hyacur1 CH_Lobangu CH_Neocon1 CH_Neolan1 CH_Obemuc1 CH_Orpsp1 CH_Pecora1 CH_Phybulla CH_Pirfi3 CH_Polagg1 CH_Powhir1 CH_Rhihy1 CH_Rhirose CH_Spipu1 CH_Terramyc CH_Triarc1"

# CR — Cryptomycota / Rozellida (2 spp.)
CLADES[CR]="CR_Parsa1 CR_Rozal1"

# EN — Entomophthorales (5 spp.)
CLADES[EN]="EN_Conco1 EN_Conth1 EN_Entmai1 EN_Mascica EN_Zoorad1"

# GL — Glomeromycota (7 spp.)
CLADES[GL]="GL_Geopyr1 GL_Gigmar1 GL_Gigro1 GL_Paroc1 GL_Rhice1 GL_Rhicl1 GL_RhiirA4"

# HO — Holozoa outgroups (20 spp.)
CLADES[HO]="HO_Acancaste HO_Callivi HO_Capowcza HO_Ciointe HO_Dicdisco HO_Drome HO_Homosap HO_Hydravulg HO_Lucicu HO_Luciser HO_Magisepte_c HO_Mneleidy HO_Monbrevi HO_Musdom HO_Naegrub HO_Pluxyl HO_Salpirose HO_Sphaearti HO_Thecatra HO_Trichadha"

# KX — Kickxellomycotina (13 spp.)
CLADES[KX]="KX_Coemoj1 KX_Coespi1 KX_DimcrSC1 KX_Dispsimple KX_Kicala1 KX_Linpe1 KX_Marpt1 KX_Pincor1 KX_Rambr1 KX_SmicuMNP KX_Smimuc2 KX_Spiroaspi KX_Zancul2"

# MI — Microsporidia (22 spp.)
CLADES[MI]="MI_AmpWSBS2006 MI_Antlo1 MI_Edhae1 MI_Enccu1 MI_Enche1 MI_Encin1 MI_Encro1 MI_Entbi1 MI_Entca1 MI_Heper1 MI_MicG_I_3 MI_Mitdap1 MI_Nemdi1 MI_Nempa1 MI_Nosap1 MI_Nosbo1 MI_Ordco1 MI_Pseneu1 MI_Sprlo1 MI_Traho1 MI_Vavcu1 MI_Vitco1"

# MO — Mortierellomycota (26 spp.)
CLADES[MO]="MO_Disdec1 MO_Disorn1 MO_Entech1 MO_HabiNA12553 MO_Lobtra1 MO_MoalAD071 MO_MoamBC1065 MO_MobeAD1040 MO_MocaN2610 MO_MochAD033 MO_MoclAM1000 MO_Mocy1230 MO_MoelGBA34 MO_MoepAD058 MO_MoepN5512 MO_MomuN6456 MO_Moramb1 MO_Morcap1 MO_MoreMES2147 MO_Morel2 MO_Morgam1 MO_Mormul1 MO_Morsel1 MO_Morwol1 MO_MowoN66265 MO_MozyN2592"

# MU — Mucoromycota / Mucorales (63 spp.)
CLADES[MU]="MU_Abscae1 MU_Absrep1 MU_Amyrou1 MU_Apoossi MU_Bacci1 MU_Baclam1 MU_Benpoi1 MU_Bifad1 MU_Blatri MU_Calth1 MU_Chabre1 MU_Chlapade MU_Chocucu1 MU_Cirumb1 MU_Cokrec1 MU_Cunech1 MU_Dicele1 MU_Dicrob1 MU_Ellano2 MU_Endsp1 MU_Fenlin MU_Gilper1 MU_Gonbut1 MU_Halrad1 MU_Helpul1 MU_Hesve2 MU_Jimfl MU_Jimlac1 MU_Kircor1 MU_Lenpar1 MU_Lenzyc1 MU_Liccor1 MU_Licra1 MU_Mucam1 MU_Mucci3 MU_Muccirc1 MU_Mucend1 MU_Muclan1 MU_Mucmuc1 MU_Mucpir1 MU_Mucrac1 MU_Mycafr1 MU_Mycind1 MU_Parpar1 MU_Phaart1 MU_Phybl4 MU_Pilano1 MU_Pilumb1 MU_PnitS607 MU_Radspe1 MU_Rhimi59 MU_Rhior3 MU_Rhipu1 MU_Rhisto1 MU_Sakvas1 MU_Spifus1 MU_Spoumb1 MU_Synrac1 MU_Syzme1 MU_Thaele1 MU_Umban1 MU_Umbra1 MU_Zycmex1"

# NU — Nucleariida (2 spp.)
CLADES[NU]="NU_Fontalba NU_Parvatl_c"

# OL — Olpidium (1 sp. — degenerate, skipped in pairwise analysis)
CLADES[OL]="OL_Olpbor1"

# ZO — Zoopagomycota (4 spp.)
CLADES[ZO]="ZO_Pipcy3 ZO_Synfus1 ZO_Synplu1 ZO_Thasp1"

export CLADES

# Ordered list of all clade codes (for looping)
export ALL_CLADES=(AP AS BA BL BS CH CR EN GL HO KX MI MO MU NU OL ZO)

# Tool paths — syntenyEDFTE micromamba environment
export MAMBA_ENV="syntenyEDFTE"
# If the env is already activated (e.g. inside a Slurm job), use direct binaries
# to avoid concurrent micromamba lock contention under set -euo pipefail.
if command -v bedtools &>/dev/null && command -v parallel &>/dev/null; then
    export BEDTOOLS="bedtools"
    export PARALLEL="parallel"
    export PYTHON="python3"
else
    export BEDTOOLS="micromamba run -n ${MAMBA_ENV} bedtools"
    export PARALLEL="micromamba run -n ${MAMBA_ENV} parallel"
    export PYTHON="micromamba run -n ${MAMBA_ENV} python3"
fi

echo "[config] EDF-TE Synteny configuration loaded."
echo "         GFF_DIR   : ${GFF_DIR}"
echo "         FAA_DIR   : ${FAA_DIR}"
echo "         OUT_DIR   : ${OUT_DIR}"
