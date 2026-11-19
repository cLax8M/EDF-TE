#!/usr/bin/env bash

# calculates TE annotation per genome using EarlGrey 


# imports / strict mode
set -euo pipefail
# define script location and load shared configuration 
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# define EarlGrey parameters
EG_THREADS="${SLURM_CPUS:-16}"
EG_SPECIES="Fungi"          # RepeatMasker species for library lookup
EG_GENOME_SIZE="200000000"  # upper bound for softmasking

# parse  arguments 
DRY_RUN=false; RESUME=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --resume)  RESUME=true;  shift ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

# define log function 
log() { echo "[$(date '+%H:%M:%S')] $*"; }

# create output dirs
mkdir -p "${TE_ANNOT_DIR}" "${LOGS_DIR}"

# Build job list 
JOB_LIST="${TE_ANNOT_DIR}/job_list.txt"
> "${JOB_LIST}"
while IFS= read -r sp; do
    [[ -z "${sp}" || "${sp}" == \#* ]] && continue
    GENOME="${GENOMES_DIR}/${sp}.fasta"
    DONE_FLAG="${TE_ANNOT_DIR}/${sp}_EarlGrey/${sp}_summaryFiles/${sp}.highLevelCount.txt"
    if [[ "${RESUME}" == true && -f "${DONE_FLAG}" ]]; then
        log "  Skipping ${sp} (already complete)"
        continue
    fi
    if [[ ! -f "${GENOME}" ]]; then
        log "  WARNING: genome not found for ${sp}: ${GENOME}"
        continue
    fi
    echo "${sp}" >> "${JOB_LIST}"
done < "${SPECIES_LIST}"

# count jobs
N=$(wc -l < "${JOB_LIST}")
log "Submitting EarlGrey for ${N} species"

# dry run exit
if [[ "${DRY_RUN}" == true ]]; then
    log "DRY-RUN: would submit array 1-${N}"
    exit 0
fi

# exit if nothing to do
if [[ "${N}" -eq 0 ]]; then
    log "Nothing to do."
    exit 0
fi

# submit SLURM array
sbatch \
    --partition="${SLURM_PARTITION}" \
    --cpus-per-task="${EG_THREADS}" \
    --mem="${SLURM_MEM}" \
    --array="1-${N}%20" \
    --output="${LOGS_DIR}/earlgrey_%a_%j.log" \
    --error="${LOGS_DIR}/earlgrey_%a_%j.err" \
    --wrap='
SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" '"${JOB_LIST}"')
GENOME='"${GENOMES_DIR}"'/${SP}.fasta
OUT_DIR='"${TE_ANNOT_DIR}"'/${SP}_EarlGrey

micromamba run -n '"${ENV_EARLGREY}"' earlGrey \
    -g "${GENOME}" \
    -s "${SP}" \
    -o "${OUT_DIR}" \
    -t '"${EG_THREADS}"' \
    -r '"${EG_SPECIES}"' \
    -l "${OUT_DIR}/${SP}EarlGrey.log"
'


# Library curation (run after all EarlGrey jobs complete)
# Submit separately with --dependency=afterok:<array_jobid>

log "To run cross-species library curation after annotation completes:"
log "  bash helpers/library_curation.sh"
