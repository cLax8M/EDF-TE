#!/bin/bash
# ============================================================
# run_fimo_cross_species.sh
# SLURM array: FIMO motif enrichment for key TF × TE class
# pairs across ALL fungal species with EarlGrey annotations.
#
# Key associations (from MU_Rhimi59 enrichment analysis):
#   LINE/I-Jockey × HMGb3  (100% hit in Rhimi59 — KEY FINDING)
#   LINE/R1        × HMGb3  (56% in Rhimi59)
#   LINE           × HMGb2
#   LTR/Gypsy      × GATA_Znf9
#   DNA/CMC-EnSpm  × HLH_Myc6
#   DNA/Merlin     × bZip3
#   ... and more (see 03b_fimo_cross_species_per_sp.py)
#
# Output:
#   ${RESULTS_DIR}/cross_species_fimo_results.tsv
#   Columns: species, te_class, motif, n_elements, n_with_motif, fraction
#
# Usage:
#   sbatch run_fimo_cross_species.sh
# ============================================================

#SBATCH --job-name=fimo_cross_sp
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-236
#SBATCH --output=/data/carlos/EDF-TE/logs/fimo_cross/fimo_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/fimo_cross/fimo_%A_%a.err

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────
SCRIPTS_DIR="/data/carlos/EDF-TE/EDF-TE/TF_TE_analysis/scripts"
GENOME_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
TE_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
MOTIF_DIR="/legserv/Temp/Carlos/EDF-TE/data/TFsDAP/motifs"
RESULTS_DIR="/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/cross_species_fimo"
WORK_DIR="${RESULTS_DIR}/work"
SPECIES_LIST="${SCRIPTS_DIR}/species_all_te.txt"
OUT_TSV="${RESULTS_DIR}/cross_species_fimo_results.tsv"

mkdir -p /data/carlos/EDF-TE/logs/fimo_cross
mkdir -p "${RESULTS_DIR}"
mkdir -p "${WORK_DIR}"

# Activate micromamba env with pandas/scipy (needed if running figure later)
# The per-sp script uses only stdlib + subprocess, so base python3 is fine.
source /etc/profile.d/lmod.sh 2>/dev/null || true
module load BEDTools/2.31.1-GCC-13.3.0
module load MEME/5.5.7-gompi-2023b

# ── Select species for this task ──────────────────────────────
SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SPECIES_LIST}")
if [[ -z "${SP}" ]]; then
    echo "[SKIP] No species at index ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# ── Resolve BED path (handle AS_AspfumUD1 symlink case) ───────
BED="${TE_DIR}/${SP}_EarlGrey/${SP}_summaryFiles/${SP}.filteredRepeats.bed"
if [[ ! -s "${BED}" ]]; then
    echo "[ERROR] No BED for ${SP}: ${BED}"
    exit 1
fi

# ── Resolve genome FASTA ──────────────────────────────────────
GENOME="${GENOME_DIR}/${SP}.fasta"
if [[ ! -s "${GENOME}" ]]; then
    echo "[ERROR] No genome for ${SP}: ${GENOME}"
    exit 1
fi

# ── Skip if this species already completed ────────────────────
# Check if output TSV already has rows for this species
if [[ -s "${OUT_TSV}" ]] && grep -q "^${SP}	" "${OUT_TSV}"; then
    echo "[SKIP] ${SP} already in ${OUT_TSV}"
    exit 0
fi

echo "[$(date)] START ${SP}"

# ── Run per-species FIMO cross-species script ─────────────────
python3 "${SCRIPTS_DIR}/03b_fimo_cross_species_per_sp.py" \
    --species  "${SP}" \
    --bed      "${BED}" \
    --genome   "${GENOME}" \
    --motif_dir "${MOTIF_DIR}" \
    --outdir   "${WORK_DIR}" \
    --out_tsv  "${OUT_TSV}"

echo "[$(date)] DONE ${SP}"
