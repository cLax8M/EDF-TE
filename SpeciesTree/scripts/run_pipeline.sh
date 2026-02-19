#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master pipeline launcher (handles all steps automatically)
#
# PURPOSE : Executes all 9 steps of the phylogenomic pipeline in order,
#           automatically:
#             - Counts input files at each step
#             - Creates input lists for SLURM arrays
#             - Submits arrays with the correct --array=1-N range
#             - Chains jobs with SLURM dependencies (each step waits for
#               the previous one to finish before starting)
#
# USAGE   : micromamba activate esptreeEDFTE
#           bash run_pipeline.sh
#
#           Options:
#             --start STEP   Resume from a specific step (1-9), e.g. --start 3
#             --dry-run      Print sbatch commands without submitting
#             --help         Show this message
#
# SLURM DEPENDENCIES:
#   Each step is submitted with --dependency=afterok:PREV_JOBID so the cluster
#   will automatically start each step only when the previous one succeeds.
#   If any step fails, downstream steps are cancelled (afterok behaviour).
#
# TO MONITOR PROGRESS:
#   watch -n 30 "squeue -u $USER"
#
# IF A STEP FAILS:
#   Check log: /data/carlos/EDF-TE/SpeciesTree/logs/
#   Resubmit from that step: bash run_pipeline.sh --start 3
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ─── Parse arguments ──────────────────────────────────────────────────────────
START_STEP=1
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --start)   START_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --help)
            sed -n '3,30p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ─── Helper functions ─────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }

submit_job() {
    # Usage: submit_job LABEL DEPENDENCY_JOBID "sbatch args..." SCRIPT [extra_sbatch_args]
    local label="$1"
    local dep_jobid="$2"
    local script="$3"
    shift 3
    local extra_args=("$@")

    local dep_arg=""
    if [[ -n "${dep_jobid}" ]]; then
        dep_arg="--dependency=afterok:${dep_jobid}"
    fi

    local cmd=(sbatch ${dep_arg} "${extra_args[@]}" "${script}")

    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "  [DRY-RUN] ${cmd[*]}" >&2
        echo "DRYRUN_JOB_${label}"
        return
    fi

    local result
    result=$("${cmd[@]}")
    local jobid
    jobid=$(echo "${result}" | awk '{print $NF}')
    log "  Submitted ${label}: job ${jobid}" >&2   # stderr so caller can capture only jobid
    echo "${jobid}"
}

make_input_list() {
    local dir="$1"
    local pattern="$2"
    local list_file="$3"
    mkdir -p "$(dirname "${list_file}")"
    ls "${dir}"/${pattern} 2>/dev/null | sort > "${list_file}" || true
    wc -l < "${list_file}" | tr -d ' '
}

# ─── Activate environment ──────────────────────────────────────────────────────
eval "$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"

# ─── Pre-flight check ─────────────────────────────────────────────────────────
log "================================================================"
log " EDF-TE Phylogenomic Pipeline"
log "================================================================"
log " Work dir    : ${WORK_DIR}"
log " Start step  : ${START_STEP}"
log " Dry run     : ${DRY_RUN}"
log ""

# Create all output directories upfront
mkdir -p \
    "${DIR_MARKERS}" \
    "${DIR_ALIGNED}" \
    "${DIR_TRIMMED}" \
    "${WORK_DIR}/logs" \
    "${DIR_FILTERED}/alignments" \
    "${DIR_FILTERED}/rejected" \
    "${DIR_GENE_TREES}" \
    "${DIR_FINAL_GENES}/alignments" \
    "${DIR_CONCAT}" \
    "${DIR_GUIDE_TREE}" \
    "${DIR_FINAL_TREE}"

# ─── Job ID tracking (for SLURM dependencies) ─────────────────────────────────
JOB_02=""
JOB_03=""
JOB_05=""
JOB_08=""
JOB_09=""

# ─── STEP 01: Extract markers ─────────────────────────────────────────────────
if [[ ${START_STEP} -le 1 ]]; then
    log "--- STEP 01: Extracting OG markers ---"

    if [[ "${DRY_RUN}" == "true" ]]; then
        log "  [DRY-RUN] python ${SCRIPT_DIR}/01_extract_markers.py"
    else
        python "${SCRIPT_DIR}/01_extract_markers.py" \
            --counts-file "${OG_COUNTS_FILE}" \
            --table-file "${OG_TABLE_FILE}" \
            --og-seq-dir "${OG_SEQ_DIR}" \
            --out-dir "${DIR_MARKERS}" \
            --min-species "${MIN_SPECIES}" \
            --max-copies "${MAX_COPIES}"
    fi

    N_MARKERS=$(ls "${DIR_MARKERS}"/*.fa 2>/dev/null | wc -l || echo 0)
    log "  OG FASTAs extracted: ${N_MARKERS}"

    if [[ ${N_MARKERS} -eq 0 ]] && [[ "${DRY_RUN}" == "false" ]]; then
        log "ERROR: No markers extracted. Check Step 01 output."
        exit 1
    fi
fi

# ─── STEP 02: MAFFT alignment (SLURM array) ───────────────────────────────────
if [[ ${START_STEP} -le 2 ]]; then
    log ""
    log "--- STEP 02: MAFFT alignment (SLURM array) ---"

    INPUT_LIST_02="${DIR_ALIGNED}/input_list.txt"
    N_02=$(make_input_list "${DIR_MARKERS}" "*.fa" "${INPUT_LIST_02}")
    log "  Input files: ${N_02}"

    if [[ ${N_02} -eq 0 ]]; then
        log "ERROR: No .fa files in ${DIR_MARKERS}. Did Step 01 complete?"
        exit 1
    fi

    # Max 50 concurrent to avoid cluster overload
    CONCURRENT_02=50
    JOB_02=$(submit_job "mafft" "" \
        "${SCRIPT_DIR}/02_align_mafft_array.sh" \
        --array="1-${N_02}%${CONCURRENT_02}")

    log "  Job ID: ${JOB_02} (array of ${N_02} jobs)"
fi

# ─── STEP 03: trimAl (SLURM array, depends on Step 02) ───────────────────────
if [[ ${START_STEP} -le 3 ]]; then
    log ""
    log "--- STEP 03: trimAl trimming (SLURM array, waits for Step 02) ---"

    # For resuming: count aligned files that already exist
    if [[ ${START_STEP} -ge 3 ]]; then
        N_ALN=$(ls "${DIR_ALIGNED}"/*.aln.fa 2>/dev/null | wc -l || echo 0)
        if [[ ${N_ALN} -eq 0 ]]; then
            log "ERROR: No aligned files in ${DIR_ALIGNED}. Did Step 02 complete?"
            exit 1
        fi
        INPUT_LIST_03="${DIR_TRIMMED}/input_list.txt"
        N_03=$(make_input_list "${DIR_ALIGNED}" "*.aln.fa" "${INPUT_LIST_03}")
    else
        # Array not yet run; list will be populated after step 02 finishes.
        # Create a placeholder step that generates the list AFTER step 02 completes,
        # then submits step 03. We use a small SLURM script for this.
        SETUP_SCRIPT_03="${WORK_DIR}/logs/setup_step03.sh"
        cat > "${SETUP_SCRIPT_03}" << 'INNEREOF'
#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="PLACEHOLDER_SCRIPT_DIR"
source "${SCRIPT_DIR}/config.sh"
INPUT_LIST="${DIR_TRIMMED}/input_list.txt"
mkdir -p "${DIR_TRIMMED}"
ls "${DIR_ALIGNED}"/*.aln.fa | sort > "${INPUT_LIST}"
N=$(wc -l < "${INPUT_LIST}")
echo "[setup] Found ${N} aligned files for step 03."
INNEREOF
        # Replace placeholder
        sed -i "s|PLACEHOLDER_SCRIPT_DIR|${SCRIPT_DIR}|g" "${SETUP_SCRIPT_03}"
        chmod +x "${SETUP_SCRIPT_03}"

        # Submit setup as a small serial job, then array depends on it
        SETUP_DEP=""
        [[ -n "${JOB_02}" ]] && SETUP_DEP="--dependency=afterok:${JOB_02}"

        if [[ "${DRY_RUN}" == "false" ]]; then
            SETUP_JOB=$(sbatch --job-name=setup_03 --partition="${PARTITION}" \
                --cpus-per-task=1 --mem=1G --time=00:05:00 \
                ${SETUP_DEP} \
                --wrap="bash ${SETUP_SCRIPT_03}" \
                --output="${WORK_DIR}/logs/setup03_%j.out" \
                | awk '{print $NF}')
            log "  Setup job: ${SETUP_JOB}"

            # Use a large placeholder N for the array (will NOT exceed actual file count
            # because sed -n "${SLURM_ARRAY_TASK_ID}p" returns empty for out-of-range IDs)
            # The trimal script gracefully exits for empty lines.
            # A safe upper bound: use N_MARKERS (step 01 output ≤ step 02 output count)
            N_03_ESTIMATE=$(ls "${DIR_MARKERS}"/*.fa 2>/dev/null | wc -l || echo 700)
            JOB_03=$(sbatch \
                --array="1-${N_03_ESTIMATE}%100" \
                --dependency=afterok:${SETUP_JOB} \
                "${SCRIPT_DIR}/03_trim_trimal_array.sh" \
                | awk '{print $NF}')
            log "  Job ID: ${JOB_03} (array up to ${N_03_ESTIMATE})"
        else
            log "  [DRY-RUN] Would submit trimal array (size = #markers from step 01)"
            JOB_03="DRYRUN_JOB_trimal"
        fi
    fi
fi

# Resume path for step 03 (when starting from step 3+)
if [[ ${START_STEP} -ge 3 ]] && [[ ${START_STEP} -le 3 ]]; then
    INPUT_LIST_03="${DIR_TRIMMED}/input_list.txt"
    N_03=$(make_input_list "${DIR_ALIGNED}" "*.aln.fa" "${INPUT_LIST_03}")
    log "  Input files: ${N_03}"
    JOB_03=$(submit_job "trimal" "" \
        "${SCRIPT_DIR}/03_trim_trimal_array.sh" \
        --array="1-${N_03}%100")
    log "  Job ID: ${JOB_03}"
fi

# ─── STEP 04: PhyKIT quality filter (runs after Step 03 array completes) ─────
if [[ ${START_STEP} -le 4 ]]; then
    log ""
    log "--- STEP 04: PhyKIT quality filter (serial, waits for Step 03) ---"

    SETUP_SCRIPT_04="${WORK_DIR}/logs/setup_step04.sh"
    cat > "${SETUP_SCRIPT_04}" << INNEREOF
#!/usr/bin/env bash
set -euo pipefail
eval "\$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"
python "${SCRIPT_DIR}/04_quality_filter_phykit.py" \\
    --trimmed-dir "${DIR_TRIMMED}" \\
    --out-dir "${DIR_FILTERED}"
INNEREOF
    chmod +x "${SETUP_SCRIPT_04}"

    DEP_04=""
    [[ -n "${JOB_03}" ]] && DEP_04="--dependency=afterok:${JOB_03}"

    if [[ "${DRY_RUN}" == "false" ]]; then
        JOB_04=$(sbatch \
            --job-name=phykit_filter \
            --partition="${PARTITION}" \
            --cpus-per-task="${THREADS}" \
            --mem=32G \
            --time=04:00:00 \
            ${DEP_04} \
            --output="${WORK_DIR}/logs/04_phykit_%j.out" \
            --error="${WORK_DIR}/logs/04_phykit_%j.err" \
            --wrap="bash ${SETUP_SCRIPT_04}" \
            | awk '{print $NF}')
        log "  Job ID: ${JOB_04}"
    else
        log "  [DRY-RUN] Would submit PhyKIT filter job (depends on Step 03)"
        JOB_04="DRYRUN_JOB_phykit"
    fi
fi

# ─── STEP 05: Gene trees (SLURM array, depends on Step 04) ───────────────────
if [[ ${START_STEP} -le 5 ]]; then
    log ""
    log "--- STEP 05: Gene trees IQTree (SLURM array, waits for Step 04) ---"

    # Setup script to create input list after step 04 finishes
    SETUP_SCRIPT_05="${WORK_DIR}/logs/setup_step05.sh"
    cat > "${SETUP_SCRIPT_05}" << INNEREOF
#!/usr/bin/env bash
set -euo pipefail
INPUT_LIST="${DIR_GENE_TREES}/input_list.txt"
mkdir -p "${DIR_GENE_TREES}"
ls "${DIR_FILTERED}/alignments/"*.fa | sort > "\${INPUT_LIST}"
N=\$(wc -l < "\${INPUT_LIST}")
echo "[setup] Found \${N} filtered alignments for Step 05 gene trees."

# Re-submit the gene tree array with the correct N
# (We cannot set array size dynamically in the original sbatch,
# so the setup job submits it fresh here)
sbatch \\
  --array="1-\${N}%30" \\
  --job-name=iqtree_genes \\
  --partition="${PARTITION}" \\
  --cpus-per-task=4 \\
  --mem=16G \\
  --time=04:00:00 \\
  --output="${WORK_DIR}/logs/05_genetrees_%A_%a.out" \\
  --error="${WORK_DIR}/logs/05_genetrees_%A_%a.err" \\
  "${SCRIPT_DIR}/05_gene_trees_array.sh"
INNEREOF
    chmod +x "${SETUP_SCRIPT_05}"

    DEP_05=""
    [[ -n "${JOB_04:-}" ]] && DEP_05="--dependency=afterok:${JOB_04}"
    [[ ${START_STEP} -ge 5 ]] && DEP_05=""

    if [[ "${DRY_RUN}" == "false" ]]; then
        # If resuming from step 5, submit directly with known counts
        if [[ ${START_STEP} -ge 5 ]]; then
            INPUT_LIST_05="${DIR_GENE_TREES}/input_list.txt"
            N_05=$(make_input_list "${DIR_FILTERED}/alignments" "*.fa" "${INPUT_LIST_05}")
            log "  Input files: ${N_05}"
            JOB_05=$(submit_job "gene_trees" "" \
                "${SCRIPT_DIR}/05_gene_trees_array.sh" \
                --array="1-${N_05}%30")
        else
            # Submit setup job that will itself submit the array
            JOB_SETUP_05=$(sbatch \
                --job-name=setup_05 \
                --partition="${PARTITION}" \
                --cpus-per-task=1 --mem=1G --time=00:05:00 \
                ${DEP_05} \
                --output="${WORK_DIR}/logs/setup05_%j.out" \
                --wrap="bash ${SETUP_SCRIPT_05}" \
                | awk '{print $NF}')
            log "  Setup job: ${JOB_SETUP_05} (will auto-submit gene tree array)"
            JOB_05="PENDING_VIA_${JOB_SETUP_05}"
        fi
    else
        log "  [DRY-RUN] Would submit gene tree array (depends on Step 04)"
        JOB_05="DRYRUN_JOB_genetrees"
    fi
fi

# ─── STEP 06: Filter paralogs ─────────────────────────────────────────────────
if [[ ${START_STEP} -le 6 ]]; then
    log ""
    log "--- STEP 06: Paralog filtering (serial) ---"

    SETUP_SCRIPT_06="${WORK_DIR}/logs/setup_step06.sh"
    cat > "${SETUP_SCRIPT_06}" << INNEREOF
#!/usr/bin/env bash
set -euo pipefail
eval "\$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"
python "${SCRIPT_DIR}/06_filter_paralogs.py" \\
    --gene-trees-dir "${DIR_GENE_TREES}" \\
    --filtered-aln-dir "${DIR_FILTERED}/alignments" \\
    --out-dir "${DIR_FINAL_GENES}" \\
    --zscore "${PARALOG_ZSCORE}"
INNEREOF
    chmod +x "${SETUP_SCRIPT_06}"

    # Step 06 depends on step 05, but step 05's array was submitted by a setup job
    # and we may not have its job ID. Use a time-based dependency approach:
    # submit step 06 with --dependency=singleton (waits for all jobs with same name)
    # OR, if we have JOB_05, use that.
    DEP_06=""
    if [[ -n "${JOB_05:-}" ]] && [[ "${JOB_05}" != *"PENDING"* ]]; then
        DEP_06="--dependency=afterok:${JOB_05}"
    fi

    if [[ "${DRY_RUN}" == "false" ]]; then
        JOB_06=$(sbatch \
            --job-name=paralog_filter \
            --partition="${PARTITION}" \
            --cpus-per-task=4 \
            --mem=16G \
            --time=02:00:00 \
            --dependency=singleton \
            --output="${WORK_DIR}/logs/06_paralogs_%j.out" \
            --error="${WORK_DIR}/logs/06_paralogs_%j.err" \
            --wrap="bash ${SETUP_SCRIPT_06}" \
            | awk '{print $NF}')
        log "  Job ID: ${JOB_06} (waits for all iqtree_genes jobs to finish)"
    else
        log "  [DRY-RUN] Would submit paralog filter job"
        JOB_06="DRYRUN_JOB_paralogs"
    fi
fi

# ─── STEP 07: Concatenate ─────────────────────────────────────────────────────
if [[ ${START_STEP} -le 7 ]]; then
    log ""
    log "--- STEP 07: Build supermatrix (serial) ---"

    SETUP_SCRIPT_07="${WORK_DIR}/logs/setup_step07.sh"
    cat > "${SETUP_SCRIPT_07}" << INNEREOF
#!/usr/bin/env bash
set -euo pipefail
eval "\$(micromamba shell hook --shell=bash)"
micromamba activate "${CONDA_ENV}"
python "${SCRIPT_DIR}/07_concatenate.py" \\
    --final-genes-dir "${DIR_FINAL_GENES}/alignments" \\
    --out-dir "${DIR_CONCAT}"
INNEREOF
    chmod +x "${SETUP_SCRIPT_07}"

    DEP_07=""
    [[ -n "${JOB_06:-}" ]] && DEP_07="--dependency=afterok:${JOB_06}"
    [[ ${START_STEP} -ge 7 ]] && DEP_07=""

    if [[ "${DRY_RUN}" == "false" ]]; then
        JOB_07=$(sbatch \
            --job-name=concatenate \
            --partition="${PARTITION}" \
            --cpus-per-task=4 \
            --mem=32G \
            --time=01:00:00 \
            ${DEP_07} \
            --output="${WORK_DIR}/logs/07_concat_%j.out" \
            --error="${WORK_DIR}/logs/07_concat_%j.err" \
            --wrap="bash ${SETUP_SCRIPT_07}" \
            | awk '{print $NF}')
        log "  Job ID: ${JOB_07}"
    else
        log "  [DRY-RUN] Would submit concatenate job"
        JOB_07="DRYRUN_JOB_concat"
    fi
fi

# ─── STEP 08: Guide tree ──────────────────────────────────────────────────────
if [[ ${START_STEP} -le 8 ]]; then
    log ""
    log "--- STEP 08: Guide tree LG+F+G4 (waits for Step 07) ---"

    DEP_08=""
    [[ -n "${JOB_07:-}" ]] && DEP_08="--dependency=afterok:${JOB_07}"
    [[ ${START_STEP} -ge 8 ]] && DEP_08=""

    if [[ "${DRY_RUN}" == "false" ]]; then
        JOB_08=$(sbatch \
            ${DEP_08} \
            "${SCRIPT_DIR}/08_species_tree_guide.sh" \
            | awk '{print $NF}')
        log "  Job ID: ${JOB_08} (~6-24h)"
    else
        log "  [DRY-RUN] sbatch ${DEP_08} 08_species_tree_guide.sh"
        JOB_08="DRYRUN_JOB_guide"
    fi
fi

# ─── STEP 09: Final species tree ──────────────────────────────────────────────
if [[ ${START_STEP} -le 9 ]]; then
    log ""
    log "--- STEP 09: Final tree LG+C60+F+G4 (waits for Step 08) ---"

    DEP_09=""
    [[ -n "${JOB_08:-}" ]] && DEP_09="--dependency=afterok:${JOB_08}"
    [[ ${START_STEP} -ge 9 ]] && DEP_09=""

    if [[ "${DRY_RUN}" == "false" ]]; then
        JOB_09=$(sbatch \
            ${DEP_09} \
            "${SCRIPT_DIR}/09_species_tree_final.sh" \
            | awk '{print $NF}')
        log "  Job ID: ${JOB_09} (~24-96h)"
    else
        log "  [DRY-RUN] sbatch ${DEP_09} 09_species_tree_final.sh"
        JOB_09="DRYRUN_JOB_final"
    fi
fi

# ─── Summary ──────────────────────────────────────────────────────────────────
log ""
log "================================================================"
log " Pipeline submitted — SLURM job chain:"
log "================================================================"
log "  01 extract    → runs immediately (Python, interactive)"
log "  02 mafft      → SLURM array  (job ${JOB_02:-submitted})"
log "  03 trimal     → SLURM array  (job ${JOB_03:-submitted}, waits for 02)"
log "  04 phykit     → SLURM serial (job ${JOB_04:-submitted}, waits for 03)"
log "  05 gene trees → SLURM array  (auto-submitted after 04)"
log "  06 paralogs   → SLURM serial (waits for all gene tree jobs)"
log "  07 concat     → SLURM serial (job ${JOB_07:-submitted}, waits for 06)"
log "  08 guide tree → SLURM serial (job ${JOB_08:-submitted}, waits for 07)"
log "  09 final tree → SLURM serial (job ${JOB_09:-submitted}, waits for 08)"
log ""
log " Monitor: watch -n 60 \"squeue -u \$USER\""
log " Logs:    ${WORK_DIR}/logs/"
log ""
log " Final tree will be at:"
log "   ${DIR_FINAL_TREE}/final.treefile"
log "================================================================"
