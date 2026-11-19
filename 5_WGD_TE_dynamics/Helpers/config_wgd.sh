#!/usr/bin/env bash

# Global configuration for WGD pipeline
# Adjust these paths to match the LEGcompute server layout.

set -euo pipefail

# Root of the code (this repo)
export WGD_CODE_ROOT="/home/ghizlane/WGD"

# Root of the data on LEGcompute (OrthoFinder outputs, proteins, etc.)
export WGD_DATA_ROOT="/data/ghizlane/project/WGD"

# OrthoFinder main directory (where Gene_Duplication_Events, Species_Tree, etc. live)
export ORTHOFINDER_DIR="${WGD_DATA_ROOT}"

# Genomes and annotations (adjusted to actual LEGcompute layout)
#   - Genomes:          /data/ghizlane/project/WGD/genomes/<species>.fasta
#   - Annotations GFF:  /data/ghizlane/project/WGD/genes_annot_gff/<species>.gff[3]
export GENOME_DIR="${WGD_DATA_ROOT}/genomes"           # FASTA of genomes
export GFF_DIR="${WGD_DATA_ROOT}/genes_annot_gff"      # GFF/GFF3 files
export PROTEIN_DIR="${WGD_DATA_ROOT}/proteins_annot_faa"     # filtered proteomes (.aa.fasta)

# Results root (relative to code root when running on server)
export RESULTS_DIR="${WGD_CODE_ROOT}/results"

# Threads and resources
export WGD_THREADS="8"

# Example: activate conda/mamba env if needed (uncomment and adapt)
# if [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
#   source "$HOME/mambaforge/etc/profile.d/conda.sh"
#   conda activate wgd_env
# fi

# Create results subdirs expected by early modules
mkdir -p "${RESULTS_DIR}/duplications" "${RESULTS_DIR}/ks" "${RESULTS_DIR}/gene_loss" "${RESULTS_DIR}/synteny" "${RESULTS_DIR}/summary"