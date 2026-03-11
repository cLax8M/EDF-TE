#!/bin/bash
# ============================================================
# run_nucseq_array.sh
# SLURM array: download SRA + process Nuc-seq for 5 species
# Pipeline: fasterq-dump → BBDuk QC → bowtie1 align →
#           Picard MarkDuplicates → DANPOS dpos
#
# Requires conda env: nucseq_pipeline
#   (sra-tools, bowtie, samtools, picard, danpos2, bbduk)
#   Install: conda create -n nucseq_pipeline -c bioconda -c conda-forge \
#              sra-tools bowtie samtools picard bedtools bbtools python=3.9
#            pip install danpos  # or conda install -c bioconda danpos
#
# Usage:
#   sbatch run_nucseq_array.sh
#   bash run_nucseq_array.sh  # local test (SLURM_ARRAY_TASK_ID must be set)
# ============================================================

#SBATCH --job-name=nucseq
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --array=1-5
#SBATCH --output=/data/carlos/EDF-TE/logs/nucseq/nucseq_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nucseq/nucseq_%A_%a.err

set -euo pipefail

# ── environment ───────────────────────────────────────────────
eval "$(micromamba shell hook --shell bash)"
micromamba activate nucseq_pipeline

# ── paths ────────────────────────────────────────────────────
GENOMES_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
OUT_DIR="/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq"
TMPDIR_BASE="/scratch/nucseq_$$"
mkdir -p "${OUT_DIR}" "${TMPDIR_BASE}"
TMPDIR="${TMPDIR_BASE}"

# ── species manifest ──────────────────────────────────────────
# FORMAT: SRR_ID  SPECIES_CODE  GENOME_FASTA
declare -A SRR=(
  [1]="SRR28069066"
  [2]="SRR28069068"
  [3]="SRR7476903"
  [4]="SRR7472154"
  [5]="SRR7472065"
)
declare -A SP=(
  [1]="BA_Mereb1"
  [2]="BA_Rhoto"
  [3]="BL_Catan2"
  [4]="AS_Usnflo1"
  [5]="BA_Mrafri1"
)
declare -A GENOME=(
  [1]="${GENOMES_DIR}/BA_Mereb1.fasta"
  [2]="${GENOMES_DIR}/BA_Rhoto.fasta"
  [3]="${GENOMES_DIR}/BL_Catan2.fasta"
  [4]="${GENOMES_DIR}/AS_Usnflo1.fasta"
  [5]="${GENOMES_DIR}/BA_Mrafri1.fasta"
)

IDX="${SLURM_ARRAY_TASK_ID}"
SRR_ID="${SRR[$IDX]}"
SPECIES="${SP[$IDX]}"
GENOME_FA="${GENOME[$IDX]}"
SP_OUT="${OUT_DIR}/${SPECIES}"
mkdir -p "${SP_OUT}"

echo "[START] ${SPECIES} — SRR: ${SRR_ID}"

# ── skip if already complete ──────────────────────────────────
if [[ -f "${SP_OUT}/${SPECIES}_danpos/pooled/${SPECIES}.wig" ]]; then
    echo "[SKIP] Already complete: ${SPECIES}"
    exit 0
fi

# ── Step 1: download SRA ──────────────────────────────────────
echo "[1/5] Downloading ${SRR_ID}..."
prefetch "${SRR_ID}" --output-directory "${TMPDIR}"
fasterq-dump "${TMPDIR}/${SRR_ID}" \
    --outdir "${TMPDIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --split-files

R1="${TMPDIR}/${SRR_ID}_1.fastq"
R2="${TMPDIR}/${SRR_ID}_2.fastq"
# Handle single-end if no R2
if [[ ! -f "${R2}" ]]; then
    R2=""
fi

# ── Step 2: QC with BBDuk ─────────────────────────────────────
echo "[2/5] QC with BBDuk..."
ADAPTERS_REF="$(micromamba info --base)/pkgs/bbtools-*/resources/adapters.fa"
ADAPTERS_REF=$(ls ${ADAPTERS_REF} 2>/dev/null | head -1 || echo "adapters")

if [[ -n "${R2}" ]]; then
    bbduk.sh \
        in1="${R1}" in2="${R2}" \
        out1="${TMPDIR}/${SPECIES}_R1_clean.fastq" \
        out2="${TMPDIR}/${SPECIES}_R2_clean.fastq" \
        ref="${ADAPTERS_REF}" \
        ktrim=r k=25 mink=11 hdist=1 \
        qtrim=r trimq=6 \
        minlen=41 minavgquality=13 maxns=0 \
        threads="${SLURM_CPUS_PER_TASK}"
    CLEAN_R1="${TMPDIR}/${SPECIES}_R1_clean.fastq"
    CLEAN_R2="${TMPDIR}/${SPECIES}_R2_clean.fastq"
else
    bbduk.sh \
        in="${R1}" \
        out="${TMPDIR}/${SPECIES}_clean.fastq" \
        ref="${ADAPTERS_REF}" \
        ktrim=r k=25 mink=11 hdist=1 \
        qtrim=r trimq=6 \
        minlen=41 minavgquality=13 maxns=0 \
        threads="${SLURM_CPUS_PER_TASK}"
    CLEAN_R1="${TMPDIR}/${SPECIES}_clean.fastq"
    CLEAN_R2=""
fi

# ── Step 3: bowtie1 align ─────────────────────────────────────
echo "[3/5] Building bowtie1 index and aligning..."
INDEX="${TMPDIR}/${SPECIES}_bt1"
bowtie-build --threads "${SLURM_CPUS_PER_TASK}" "${GENOME_FA}" "${INDEX}"

if [[ -n "${CLEAN_R2}" ]]; then
    bowtie -S -t -m 1 -v 2 \
        -p "${SLURM_CPUS_PER_TASK}" \
        "${INDEX}" \
        -1 "${CLEAN_R1}" -2 "${CLEAN_R2}" \
        | samtools sort -@ "${SLURM_CPUS_PER_TASK}" \
            -o "${TMPDIR}/${SPECIES}.bam"
else
    bowtie -S -t -m 1 -v 2 \
        -p "${SLURM_CPUS_PER_TASK}" \
        "${INDEX}" "${CLEAN_R1}" \
        | samtools sort -@ "${SLURM_CPUS_PER_TASK}" \
            -o "${TMPDIR}/${SPECIES}.bam"
fi
samtools index "${TMPDIR}/${SPECIES}.bam"

# ── Step 4: Picard MarkDuplicates ─────────────────────────────
echo "[4/5] Removing PCR duplicates..."
picard MarkDuplicates \
    I="${TMPDIR}/${SPECIES}.bam" \
    O="${TMPDIR}/${SPECIES}_dedup.bam" \
    M="${SP_OUT}/${SPECIES}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT
samtools index "${TMPDIR}/${SPECIES}_dedup.bam"

# ── Step 5: DANPOS nucleosome calling ─────────────────────────
echo "[5/5] DANPOS nucleosome occupancy calling..."
mkdir -p "${SP_OUT}/${SPECIES}_danpos"
python $(which danpos.py) dpos \
    "${TMPDIR}/${SPECIES}_dedup.bam" \
    -m 1 -a 1 \
    -o "${SP_OUT}/${SPECIES}_danpos"

echo "[DONE] ${SPECIES} — WIG at ${SP_OUT}/${SPECIES}_danpos/"

# ── cleanup tmp ───────────────────────────────────────────────
rm -rf "${TMPDIR_BASE}"
