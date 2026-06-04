#!/usr/bin/env bash
# =============================================================================
# setup_chimeraTE.sh
# Install ChimeraTE and prepare per-species analysis scripts for EDF-TE.
#
# ChimeraTE detects chimeric transcripts (gene + TE fusion transcripts)
# from RNA-seq alignments.
# GitHub: https://github.com/OliveiraDS-hub/ChimeraTE
#
# Usage:
#   1. Install:   bash setup_chimeraTE.sh --install
#   2. Prepare:   bash setup_chimeraTE.sh --prepare --species MU_Mlus402
#   3. Run all:   bash setup_chimeraTE.sh --submit-all
# =============================================================================

set -euo pipefail

INSTALL_DIR=/data/carlos/EDF-TE/tools/ChimeraTE
CONDA_ENV=chimeraTE
OUT_BASE=/legserv/Temp/Carlos/EDF-TE/results/Co-option/chimerateTE
RNASEQ_DIR=/legserv/Temp/Carlos/EDF-TE/data/rnaseq   # deposit RNA-seq reads here
GENOMES_DIR=/legserv/Temp/Carlos/EDF-TE/data/genomes
GENE_GFF_DIR=/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff
TE_BED_BASE=/data/carlos/EDF-TE/TE_analysis/annotated
LOGS_DIR=/data/carlos/EDF-TE/logs/chimeraTE
SCRIPT_DIR=/data/carlos/EDF-TE/EDF-TE/Cooption/scripts/chimeraTE_jobs

mkdir -p "$OUT_BASE" "$LOGS_DIR" "$SCRIPT_DIR"

# ─── Target species for RNA-seq (EDF clades covered) ─────────────────────────
# User-provided reads (place in $RNASEQ_DIR/<species>/<R1.fastq.gz> <R2.fastq.gz>)
TARGET_SPECIES=(
  MU_Mlus402      # Mucor lusitanicus 402 (Mucoromycota)
  MU_Rblksl1      # Rhizopus microsporus (Mucoromycota)
  MU_Phybl4       # Phycomyces blakesleeanus (Mucoromycota)
  # Add more species as RNA-seq data becomes available:
  # AS_Afum1       # Aspergillus fumigatus (Ascomycota)
  # BA_Cneofor     # Cryptococcus neoformans (Basidiomycota)
  # CH_Alloma1     # Allomyces macrogynus (Chytridiomycota)
)

# ─── INSTALL ─────────────────────────────────────────────────────────────────
do_install () {
  echo "[INSTALL] Cloning ChimeraTE…"
  git clone https://github.com/OliveiraDS-hub/ChimeraTE.git "$INSTALL_DIR" || {
    echo "Already cloned, updating…"
    git -C "$INSTALL_DIR" pull
  }

  echo "[INSTALL] Creating conda environment ${CONDA_ENV}…"
  eval "$(micromamba shell hook --shell bash)"
  micromamba create -n "$CONDA_ENV" \
    -c bioconda -c conda-forge -c defaults \
    star=2.7.11a \
    samtools=1.19 \
    bedtools=2.31 \
    python=3.10 \
    r-base=4.3 \
    r-ggplot2 r-dplyr r-data.table \
    bioconductor-deseq2 \
    hisat2 \
    -y || echo "Environment already exists"

  echo "[INSTALL] Done. Activate with: micromamba activate ${CONDA_ENV}"
  echo "ChimeraTE binary: ${INSTALL_DIR}/ChimeraTE.py"
}

# ─── PREPARE per-species job ──────────────────────────────────────────────────
prepare_species () {
  local sp=$1
  local RNASEQ="${RNASEQ_DIR}/${sp}"
  local GENOME="${GENOMES_DIR}/${sp}.fasta"
  local GFF="${GENE_GFF_DIR}/${sp}.gff"
  local TE_BED="${TE_BED_BASE}/${sp}_EarlGrey/${sp}_summaryFiles/${sp}.filteredRepeats.bed"
  local OUT="${OUT_BASE}/${sp}"
  local SCRIPT="${SCRIPT_DIR}/run_${sp}.sh"
  mkdir -p "$OUT"

  cat > "$SCRIPT" << SLURM_EOF
#!/bin/bash
#SBATCH --job-name=chimeraTE_${sp}
#SBATCH --partition=normal.10
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=${LOGS_DIR}/${sp}_%j.out
#SBATCH --error=${LOGS_DIR}/${sp}_%j.err

set -euo pipefail

eval "\$(micromamba shell hook --shell bash)"
micromamba activate ${CONDA_ENV}

RNASEQ=${RNASEQ}
GENOME=${GENOME}
GFF=${GFF}
TE_BED=${TE_BED}
OUT=${OUT}

# Sanity checks
[[ -d "\$RNASEQ" ]] || { echo "ERROR: RNA-seq dir missing: \$RNASEQ"; exit 1; }
[[ -f "\$GENOME" ]] || { echo "ERROR: Genome missing: \$GENOME"; exit 1; }
[[ -f "\$TE_BED" ]] || { echo "ERROR: TE BED missing: \$TE_BED"; exit 1; }

# Detect paired-end reads
R1=\$(ls "\$RNASEQ"/*_R1*.fastq.gz 2>/dev/null | head -1 || ls "\$RNASEQ"/*_1.fastq.gz 2>/dev/null | head -1)
R2=\$(ls "\$RNASEQ"/*_R2*.fastq.gz 2>/dev/null | head -1 || ls "\$RNASEQ"/*_2.fastq.gz 2>/dev/null | head -1)
[[ -n "\$R1" ]] || { echo "ERROR: No R1 reads in \$RNASEQ"; exit 1; }

echo "Date: \$(date)"
echo "Species: ${sp}"
echo "R1: \$R1"
echo "R2: \$R2"

# Step 1: Build STAR index
STAR_IDX="\${OUT}/star_index"
mkdir -p "\$STAR_IDX"
STAR --runMode genomeGenerate \
     --genomeDir "\$STAR_IDX" \
     --genomeFastaFiles "\$GENOME" \
     --sjdbGFFfile "\$GFF" \
     --runThreadN 8 \
     --genomeSAindexNbases 11   # adjust for small genomes; default 14

# Step 2: Align with STAR
BAM="\${OUT}/${sp}_aligned.bam"
STAR --runMode alignReads \
     --genomeDir "\$STAR_IDX" \
     --readFilesIn "\$R1" \$([[ -n "\$R2" ]] && echo "\$R2" || echo "") \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix "\${OUT}/${sp}_" \
     --outSAMattributes NH HI AS NM \
     --runThreadN 8 \
     --limitBAMsortRAM 24000000000
samtools index "\${OUT}/${sp}_Aligned.sortedByCoord.out.bam"

# Step 3: Run ChimeraTE
python3 ${INSTALL_DIR}/ChimeraTE.py \
  --bam "\${OUT}/${sp}_Aligned.sortedByCoord.out.bam" \
  --gtf "\$GFF" \
  --te_bed "\$TE_BED" \
  --genome "\$GENOME" \
  --outdir "\$OUT" \
  --threads 8

echo "Done: \$(date)"
SLURM_EOF

  chmod +x "$SCRIPT"
  echo "[PREPARED] ${SCRIPT}"
}

# ─── SUBMIT ALL ───────────────────────────────────────────────────────────────
submit_all () {
  for sp in "${TARGET_SPECIES[@]}"; do
    prepare_species "$sp"
    if [[ -d "${RNASEQ_DIR}/${sp}" ]]; then
      echo "[SUBMIT] ${sp}"
      sbatch "${SCRIPT_DIR}/run_${sp}.sh"
    else
      echo "[SKIP] ${sp}: RNA-seq not yet deposited at ${RNASEQ_DIR}/${sp}/"
    fi
  done
}

# ─── DISPATCH ────────────────────────────────────────────────────────────────
case "${1:-}" in
  --install)    do_install ;;
  --prepare)    prepare_species "${2:?Usage: --prepare SPECIES}" ;;
  --submit-all) submit_all ;;
  *)
    echo "Usage: bash setup_chimeraTE.sh --install | --prepare SP | --submit-all"
    echo ""
    echo "Next steps:"
    echo "  1. bash setup_chimeraTE.sh --install"
    echo "  2. Deposit RNA-seq reads at: ${RNASEQ_DIR}/<species>/"
    echo "  3. bash setup_chimeraTE.sh --submit-all"
    echo ""
    echo "Target species:"
    for sp in "${TARGET_SPECIES[@]}"; do
      rnaseq_ready=$([[ -d "${RNASEQ_DIR}/${sp}" ]] && echo "READY" || echo "waiting for reads")
      te_ready=$([[ -f "${TE_BED_BASE}/${sp}_EarlGrey/${sp}_summaryFiles/${sp}.filteredRepeats.bed" ]] && echo "OK" || echo "missing")
      echo "  ${sp}  RNA-seq:${rnaseq_ready}  TE_annotation:${te_ready}"
    done
    ;;
esac
