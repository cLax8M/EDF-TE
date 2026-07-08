#!/usr/bin/env bash
#SBATCH --job-name=WGD_Ks_PAML
#SBATCH --output=/home/ghizlane/WGD/results/ks/paml_yn00/slurm-%A_%a.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00

set -euo pipefail

# Load MAFFT module and activate micromamba environment with PAL2NAL + PAML
module load MAFFT
if [ -f "$HOME/micromamba/etc/profile.d/micromamba.sh" ]; then
  # shellcheck disable=SC1090
  source "$HOME/micromamba/etc/profile.d/micromamba.sh"
  micromamba activate wgd
fi

# Use absolute paths because within SLURM the script is copied to a spool directory
SCRIPT_DIR="/home/ghizlane/WGD/scripts/ks_wgd"
source "/home/ghizlane/WGD/scripts/core/config_wgd.sh"

SPECIES_LIST_FILE="${RESULTS_DIR}/ks/species_for_paml_ks.txt"
PARALOG_TSV="${RESULTS_DIR}/ks/paralog_pairs_for_ks.filtered.tsv"
CDS_DIR="${RESULTS_DIR}/ks/cds"

# Build species list (species with both pairs and CDS) if it does not exist
if [ ! -f "${SPECIES_LIST_FILE}" ]; then
  echo "[WGD][KsPAML] Building species list at ${SPECIES_LIST_FILE}" >&2
  tmp_sp_pairs=$(mktemp)
  tmp_sp_cds=$(mktemp)
  cut -f1 "${PARALOG_TSV}" | tail -n +2 | sort -u > "${tmp_sp_pairs}"
  ls "${CDS_DIR}"/*.cds.fa 2>/dev/null | xargs -n1 basename | sed 's/\.cds\.fa$//' | sort -u > "${tmp_sp_cds}"
  comm -12 "${tmp_sp_pairs}" "${tmp_sp_cds}" > "${SPECIES_LIST_FILE}"
  rm -f "${tmp_sp_pairs}" "${tmp_sp_cds}"
fi

if [ ! -s "${SPECIES_LIST_FILE}" ]; then
  echo "[WGD][KsPAML] ERROR: species list ${SPECIES_LIST_FILE} is empty; nothing to run." >&2
  exit 1
fi

# Derive species for this array task
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
LINE_NO=$((TASK_ID + 1))
SP=$(sed -n "${LINE_NO}p" "${SPECIES_LIST_FILE}")

if [ -z "${SP:-}" ]; then
  echo "[WGD][KsPAML] WARNING: no species for array index ${TASK_ID}; exiting." >&2
  exit 0
fi

# Optionally activate a micromamba/conda env here if needed
# source "$HOME/micromamba/etc/profile.d/micromamba.sh"
# micromamba activate wgd_tools

cd "${WGD_CODE_ROOT}"

bash "${SCRIPT_DIR}/03_compute_ks_paml.sh" "${SP}"
