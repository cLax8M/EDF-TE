#!/bin/bash
#SBATCH --job-name=fimo_TE
#SBATCH --output=/data/carlos/EDF-TE/logs/fimo_TE_%A_%a.out
#SBATCH --error=/data/carlos/EDF-TE/logs/fimo_TE_%A_%a.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-34%6     # 35 TE classes (0-34), max 6 simultaneous

#
# SLURM array script: run 03_fimo_motif_enrichment.py for one TE class per task.
#
# To submit:
#   1. Run: python 03_fimo_motif_enrichment.py --list-classes
#      → prints the list of TE classes and the correct --array upper bound
#   2. Adjust the --array line above accordingly
#   3. sbatch run_fimo_enrichment.sh
#
# After all tasks complete:
#   python 03_fimo_motif_enrichment.py --merge-only
#   → merges all per-class TSVs and generates the final figures
#

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────
PYTHON=/data/carlos/EDF-TE/PanGenome/.venv/bin/python
SCRIPT=/data/carlos/EDF-TE/EDF-TE/TF_TE_analysis/scripts/03_fimo_motif_enrichment.py
OUTDIR=/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/03_fimo_enrichment

mkdir -p /data/carlos/EDF-TE/logs

# ── Load modules ────────────────────────────────────────────────────────────
module purge
module load BEDTools/2.31.1-GCC-13.3.0
module load MEME/5.5.7-gompi-2023b

# ── TE class list (generated once by --list-classes mode) ──────────────────
# If the file doesn't exist yet, create it from the BED
TE_CLASS_LIST=${OUTDIR}/te_class_list.txt
if [[ ! -f "${TE_CLASS_LIST}" ]]; then
    "${PYTHON}" "${SCRIPT}" --list-classes --outdir "${OUTDIR}"
fi

# ── Pick this task's TE class ───────────────────────────────────────────────
TE_CLASS=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${TE_CLASS_LIST}")
if [[ -z "${TE_CLASS}" ]]; then
    echo "No TE class found for task ${SLURM_ARRAY_TASK_ID}, exiting"
    exit 0
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: processing TE class '${TE_CLASS}'"
echo "Start: $(date)"

"${PYTHON}" "${SCRIPT}" \
    --te-class "${TE_CLASS}" \
    --outdir "${OUTDIR}"

echo "End: $(date)"
