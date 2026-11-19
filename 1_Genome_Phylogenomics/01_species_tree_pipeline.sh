#!/usr/bin/env bash

# builds a phylogenomic species tree using a concatenated supermatrix approach from single-copy orthologous markers inferred by OrthoFinder.


set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# Per-section paths 
ORTHO_DIR="${RESULTS_DIR}/orthofinder"
WORK_DIR="${S1_DIR}/species_tree"
LOGS="${WORK_DIR}/logs"

DIR_MARKERS="${WORK_DIR}/01_markers"
DIR_ALIGNED="${WORK_DIR}/02_aligned"
DIR_TRIMMED="${WORK_DIR}/03_trimmed"
DIR_FILTERED="${WORK_DIR}/04_filtered"
DIR_GENE_TREES="${WORK_DIR}/05_gene_trees"
DIR_FINAL_GENES="${WORK_DIR}/06_final_genes"
DIR_CONCAT="${WORK_DIR}/07_concatenated"
DIR_GUIDE_TREE="${WORK_DIR}/08_guide_tree"
DIR_FINAL_TREE="${WORK_DIR}/09_final_tree"

# Filtering thresholds 
MIN_SPECIES=95      # min genomes per marker (≥40% of 238)
MAX_COPIES=2        # max gene copies per species per marker
MIN_ALN_LEN=60      # min alignment length (AA columns) after trimming
MIN_PIS=25          # min parsimony-informative sites
MAX_RCV=0.5         # max relative composition variability (PhyKIT)
PARALOG_ZSCORE=3.5  # branch-length z-score for paralog flagging
THREADS="${SLURM_CPUS:-8}"

# Argument parsing 
START_STEP=1; DRY_RUN=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --start)   START_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --help) grep '^#' "$0" | head -30 | sed 's/^# \?//'; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

log() { echo "[$(date '+%H:%M:%S')] $*"; }
run() {
    local label="$1"; shift
    log "Running: $label"
    if [[ "${DRY_RUN}" == true ]]; then echo "  DRY-RUN: $*"; return; fi
    "$@"
}

mkdir -p "${LOGS}" "${DIR_MARKERS}" "${DIR_ALIGNED}" "${DIR_TRIMMED}" \
         "${DIR_FILTERED}" "${DIR_GENE_TREES}" "${DIR_FINAL_GENES}" \
         "${DIR_CONCAT}" "${DIR_GUIDE_TREE}" "${DIR_FINAL_TREE}"


# STEP 1 — Extract marker sequences

if [[ "${START_STEP}" -le 1 ]]; then
    log "=== STEP 1: Extract marker sequences ==="
    run "extract_markers" python3 "${SCRIPT_DIR}/helpers/extract_markers.py" \
        --ortho_dir  "${ORTHO_DIR}" \
        --out_dir    "${DIR_MARKERS}" \
        --min_sp     "${MIN_SPECIES}" \
        --max_copies "${MAX_COPIES}"
fi


# STEP 2 — Align with MAFFT (SLURM array)

if [[ "${START_STEP}" -le 2 ]]; then
    log "=== STEP 2: MAFFT alignment ==="
    ls "${DIR_MARKERS}"/*.faa > "${WORK_DIR}/markers_list.txt"
    N=$(wc -l < "${WORK_DIR}/markers_list.txt")
    log "  Aligning ${N} markers"
    if [[ "${DRY_RUN}" == false ]]; then
        sbatch --partition="${SLURM_PARTITION}" \
               --cpus-per-task="${THREADS}" --mem="${SLURM_MEM}" \
               --array="1-${N}%50" \
               --output="${LOGS}/mafft_%a.log" \
               --wrap='FA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" '"${WORK_DIR}/markers_list.txt"')
                       BASE=$(basename "${FA}" .faa)
                       mafft --auto --thread '"${THREADS}"' "${FA}" \
                           > '"${DIR_ALIGNED}/${BASE}"'.faa'
    fi
fi


# STEP 3 — Trim with trimAl

if [[ "${START_STEP}" -le 3 ]]; then
    log "=== STEP 3: trimAl ==="
    ls "${DIR_ALIGNED}"/*.faa > "${WORK_DIR}/aligned_list.txt"
    N=$(wc -l < "${WORK_DIR}/aligned_list.txt")
    if [[ "${DRY_RUN}" == false ]]; then
        sbatch --partition="${SLURM_PARTITION}" \
               --cpus-per-task=1 --mem="4G" \
               --array="1-${N}%100" \
               --output="${LOGS}/trimal_%a.log" \
               --wrap='FA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" '"${WORK_DIR}/aligned_list.txt"')
                       BASE=$(basename "${FA}" .faa)
                       trimal -in "${FA}" -out '"${DIR_TRIMMED}/${BASE}"'.faa -automated1'
    fi
fi


# STEP 4 — Quality filter (PhyKIT)

if [[ "${START_STEP}" -le 4 ]]; then
    log "=== STEP 4: PhyKIT quality filter ==="
    run "quality_filter" python3 "${SCRIPT_DIR}/helpers/quality_filter.py" \
        --in_dir     "${DIR_TRIMMED}" \
        --out_dir    "${DIR_FILTERED}" \
        --min_len    "${MIN_ALN_LEN}" \
        --min_pis    "${MIN_PIS}" \
        --max_rcv    "${MAX_RCV}"
fi


# STEP 5 — Per-gene trees (IQ-TREE, LG+G4)

if [[ "${START_STEP}" -le 5 ]]; then
    log "=== STEP 5: Per-gene trees ==="
    ls "${DIR_FILTERED}"/*.faa > "${WORK_DIR}/filtered_list.txt"
    N=$(wc -l < "${WORK_DIR}/filtered_list.txt")
    if [[ "${DRY_RUN}" == false ]]; then
        sbatch --partition="${SLURM_PARTITION}" \
               --cpus-per-task="${THREADS}" --mem="8G" \
               --array="1-${N}%30" \
               --output="${LOGS}/iqtree_gene_%a.log" \
               --wrap='FA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" '"${WORK_DIR}/filtered_list.txt"')
                       BASE=$(basename "${FA}" .faa)
                       iqtree2 -s "${FA}" -m LG+G4 -T '"${THREADS}"' \
                           --prefix '"${DIR_GENE_TREES}/${BASE}"' -redo -quiet'
    fi
fi


# STEP 6 — Paralog removal (branch-length z-score)

if [[ "${START_STEP}" -le 6 ]]; then
    log "=== STEP 6: Paralog removal ==="
    run "filter_paralogs" python3 "${SCRIPT_DIR}/helpers/filter_paralogs.py" \
        --tree_dir   "${DIR_GENE_TREES}" \
        --seq_dir    "${DIR_FILTERED}" \
        --out_dir    "${DIR_FINAL_GENES}" \
        --zscore     "${PARALOG_ZSCORE}"
fi


# STEP 7 — Concatenation into supermatrix

if [[ "${START_STEP}" -le 7 ]]; then
    log "=== STEP 7: Concatenation ==="
    run "concatenate" python3 "${SCRIPT_DIR}/helpers/concatenate.py" \
        --in_dir  "${DIR_FINAL_GENES}" \
        --out_dir "${DIR_CONCAT}"
fi


# STEP 8 — Guide tree (LG+F+G4)

if [[ "${START_STEP}" -le 8 ]]; then
    log "=== STEP 8: Guide tree ==="
    if [[ "${DRY_RUN}" == false ]]; then
        sbatch --partition="${SLURM_PARTITION}" \
               --cpus-per-task="${THREADS}" --mem="${SLURM_MEM}" \
               --output="${LOGS}/iqtree_guide.log" \
               --wrap="iqtree2 -s ${DIR_CONCAT}/supermatrix.faa \
                          -p ${DIR_CONCAT}/partitions.nex \
                          -m LG+F+G4 -T ${THREADS} \
                          --prefix ${DIR_GUIDE_TREE}/guide -redo"
    fi
fi


# STEP 9 — Final tree (LG+C60+F+G4, 1000 ultrafast bootstraps)

if [[ "${START_STEP}" -le 9 ]]; then
    log "=== STEP 9: Final species tree ==="
    if [[ "${DRY_RUN}" == false ]]; then
        sbatch --partition="${SLURM_PARTITION}" \
               --cpus-per-task="${THREADS}" --mem="64G" \
               --output="${LOGS}/iqtree_final.log" \
               --wrap="iqtree2 -s ${DIR_CONCAT}/supermatrix.faa \
                          -p ${DIR_CONCAT}/partitions.nex \
                          -m LG+C60+F+G4 \
                          -t ${DIR_GUIDE_TREE}/guide.treefile \
                          -T ${THREADS} -B 1000 \
                          --prefix ${DIR_FINAL_TREE}/species_tree -redo"
    fi
fi


# STEP 10 — Root and annotate

if [[ "${START_STEP}" -le 10 ]]; then
    log "=== STEP 10: Root and annotate ==="
    run "root_tree" python3 "${SCRIPT_DIR}/helpers/root_and_annotate.py" \
        --tree   "${DIR_FINAL_TREE}/species_tree.treefile" \
        --outgroup "HO" \
        --output "${S1_DIR}/species_tree_rooted.nwk"
    log "Final tree written to: ${S1_DIR}/species_tree_rooted.nwk"
fi

log "Pipeline complete."
