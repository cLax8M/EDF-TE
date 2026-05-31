#!/usr/bin/env bash
# =============================================================================
# config.sh  —  shared configuration for the EDF-TE analysis pipeline
#
# Edit BASE_DIR and RESULTS_DIR to match your system before running any
# section script.  All other paths are derived from these two variables.
# =============================================================================

# --- Root directories --------------------------------------------------------
BASE_DIR="${BASE_DIR:-/path/to/edf-te}"          # code + small data
RESULTS_DIR="${RESULTS_DIR:-/path/to/results}"   # large outputs

# --- Data sub-directories ----------------------------------------------------
DATA_DIR="${BASE_DIR}/data"
GENOMES_DIR="${DATA_DIR}/genomes"
GFF_DIR="${DATA_DIR}/genes_annot_gff"
FAA_DIR="${DATA_DIR}/proteins_annot_faa"
IPR_DIR="${DATA_DIR}/functional_interpro"
TE_ANNOT_DIR="${RESULTS_DIR}/TE_annotation"      # EarlGrey per-species folders

# --- Reference species (three focal Mucoromycotina genomes) ------------------
REF_MUCCI="${GENOMES_DIR}/MU_Mucci3.fasta"
REF_RHIMI="${GENOMES_DIR}/MU_Rhimi59.fasta"
REF_PHYBL="${GENOMES_DIR}/MU_Phybl4.fasta"

# --- Species / taxonomy metadata ---------------------------------------------
SPECIES_LIST="${DATA_DIR}/species_list.txt"
CLADE_MAP="${DATA_DIR}/clade_codes.tsv"
PHYLO_TREE="${RESULTS_DIR}/species_tree/species_tree_rooted.nwk"

# --- Results sub-directories (created by each section) -----------------------
S1_DIR="${RESULTS_DIR}/section1_genome_phylogenetics"
S2_DIR="${RESULTS_DIR}/section2_te_landscape"
S3_DIR="${RESULTS_DIR}/section3_te_epigenetics"
S4_DIR="${RESULTS_DIR}/section4_te_cooption"
S5_DIR="${RESULTS_DIR}/section5_wgd_te_dynamics"
S6_DIR="${RESULTS_DIR}/section6_htt"

# --- Software (conda/micromamba environments) --------------------------------
ENV_BASE="base"
ENV_EARLGREY="earlgrey"
ENV_PHYLO="phylo"
ENV_PYTHON="edf_py"
ENV_R="edf_r"

# --- Compute resources (SLURM defaults, override per-job as needed) ----------
SLURM_PARTITION="normal.72h"
SLURM_CPUS=8
SLURM_MEM="32G"
LOGS_DIR="${RESULTS_DIR}/logs"

mkdir -p "${LOGS_DIR}"
