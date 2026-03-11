#!/bin/bash
# Launch full EarlGrey annotation for all genomes using HoloTE_combined library
# Usage: bash launch_earlgrey_annotate.sh [--dry-run]
#
# What it does:
#   - Iterates all real genome files (skips symlinks) in GENOME_DIR
#   - Submits one SLURM job per genome running:
#       earlGrey -g GENOME -s SPECIES -o OUTPUT_DIR -t 16 -r eukarya -l LIBRARY
#   - Skips genomes that already have a completed output directory
#
# Prerequisites:
#   - HoloTE_combined_v1.fa must exist at COMBINED_LIB path
#   - EarlGrey module available via: module use /data/toby/tools/modulefiles && module load earlgrey

GENOME_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
OUTPUT_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
LOG_DIR="/data/carlos/EDF-TE/logs/annotate"
SCRIPT_DIR="/data/carlos/EDF-TE/logs/sbatch_annotate"
COMBINED_LIB="/data/carlos/EDF-TE/TE_analysis/Library_curation/HoloTE_combined_v1.fa"
THREADS=16
MEM="80G"
PARTITION="normal.1000h"

DRY_RUN=0
[[ "$1" == "--dry-run" ]] && DRY_RUN=1

# ── Sanity checks ─────────────────────────────────────────────────────────────
if [[ ! -f "$COMBINED_LIB" ]]; then
    if [[ $DRY_RUN -eq 0 ]]; then
        echo "ERROR: Combined library not found: $COMBINED_LIB"
        echo "       Wait for HoloTE_nonfungal pipeline to finish, then run:"
        echo "       bash combine_libs.sh"
        exit 1
    else
        echo "[WARN] Combined library not found yet (dry-run only — OK)"
    fi
fi

mkdir -p "$LOG_DIR" "$SCRIPT_DIR" "$OUTPUT_DIR"

echo "================================================================"
echo "EarlGrey annotation launcher"
echo "  Genome dir : $GENOME_DIR"
echo "  Library    : $COMBINED_LIB"
echo "  Output     : $OUTPUT_DIR"
echo "  Dry run    : $DRY_RUN"
echo "================================================================"

SUBMITTED=0
SKIPPED_DONE=0
SKIPPED_MISSING=0

for GENOME_PATH in "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fna; do
    # Skip if glob matched nothing
    [[ -e "$GENOME_PATH" ]] || continue

    # Skip symlinks (duplicates: _c and _1 variants)
    [[ -L "$GENOME_PATH" ]] && continue

    GENOME_FILE=$(basename "$GENOME_PATH")
    # Strip extension to get species code (e.g. MU_Mucci3)
    SPECIES="${GENOME_FILE%.*}"

    # Skip if output already exists (has a _EarlGrey dir with final output)
    OUT_SPECIES_DIR="${OUTPUT_DIR}/${SPECIES}_EarlGrey"
    if [[ -d "$OUT_SPECIES_DIR" ]]; then
        STRAINED=$(find "$OUT_SPECIES_DIR" -name "*.strained" ! -path "*/TS_*" 2>/dev/null | head -1)
        if [[ -n "$STRAINED" ]]; then
            echo "SKIP (done): $SPECIES"
            ((SKIPPED_DONE++))
            continue
        fi
    fi

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[DRY-RUN] $SPECIES  →  $(du -sh "$GENOME_PATH" 2>/dev/null | cut -f1)"
        ((SUBMITTED++))
        continue
    fi

    # Generate sbatch script
    JOB_NAME="${SPECIES}_EGA"
    SCRIPT="${SCRIPT_DIR}/${SPECIES}.sh"

    cat > "$SCRIPT" <<SBATCH
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --partition=${PARTITION}
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=${MEM}
#SBATCH --output=${LOG_DIR}/${SPECIES}_%j.log

module use /data/toby/tools/modulefiles
module load earlgrey

earlGrey \\
    -g ${GENOME_PATH} \\
    -s ${SPECIES} \\
    -o ${OUTPUT_DIR} \\
    -t ${THREADS} \\
    -r eukarya \\
    -l ${COMBINED_LIB}
SBATCH

    JOB_ID=$(sbatch "$SCRIPT" | awk '{print $4}')
    echo "Submitted $SPECIES → job $JOB_ID"
    ((SUBMITTED++))

    # Small delay to avoid overwhelming the scheduler
    sleep 0.2
done

echo "================================================================"
if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY-RUN complete: $SUBMITTED genomes would be submitted"
else
    echo "Done: $SUBMITTED submitted, $SKIPPED_DONE already complete"
fi
echo "Monitor with: squeue -u \$USER | grep EGA"
echo "================================================================"
