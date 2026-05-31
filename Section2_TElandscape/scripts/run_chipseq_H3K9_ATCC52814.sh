#!/bin/bash
# ============================================================
# run_chipseq_H3K9_ATCC52814.sh
# H3K9me3 ChIP-seq pipeline for R. microsporus ATCC52814
# Input reads from ATCC11559 (same species, cross-mapping)
#
# Steps:
#   1. Index genome (bowtie2-build + samtools faidx)
#   2. Quality trim (fastp, paired-end)
#   3. Align IP + Input to ATCC52814 (bowtie2 --very-sensitive)
#   4. Sort, mark duplicates (picard), index BAMs
#   5. bigwig coverage (deeptools bamCoverage, CPM, bin=10bp)
#   6. bigwig log2 IP/Input ratio (bamCompare, bin=50bp)
#   7. Peak calling (MACS2 broad, BAMPE mode)
# ============================================================

#SBATCH --job-name=chipseq_H3K9_ATCC52814
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/data/carlos/EDF-TE/logs/chipseq/chipseq_H3K9_ATCC52814_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/chipseq/chipseq_H3K9_ATCC52814_%j.err

set -euo pipefail

# ── Activate micromamba environment ──────────────────────────
eval "$(micromamba shell hook --shell bash)"
micromamba activate chipseq_env

# ── Paths ────────────────────────────────────────────────────
READS_DIR="/legserv/Temp/Carlos/EDF-TE/data/Extra/ChIPRMH3K9me3"
GENOME_FA="/legserv/Temp/Carlos/EDF-TE/data/Extra/ATCC52814.fa"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/Epigenetics/RmH3K9me3_ATCC52814"
SCRATCH="/scratch/chipseq_H3K9_${SLURM_JOB_ID}"
THREADS="${SLURM_CPUS_PER_TASK}"
GENOME_SIZE="2.7e7"   # Rh. microsporus ~27Mb

mkdir -p "${OUT_DIR}" "${SCRATCH}"

# ── Input files ───────────────────────────────────────────────
IP_R1="${READS_DIR}/H3K9_IP_WT_1.fq.gz"
IP_R2="${READS_DIR}/H3K9_IP_WT_2.fq.gz"
IN_R1="${READS_DIR}/H3K9InputWT_1.fq.gz"
IN_R2="${READS_DIR}/H3K9InputWT_2.fq.gz"

echo "[$(date)] === H3K9me3 ChIP-seq pipeline started ==="
echo "[$(date)] Genome: ${GENOME_FA}"
echo "[$(date)] IP reads: ${IP_R1} / ${IP_R2}"
echo "[$(date)] Input reads: ${IN_R1} / ${IN_R2}"

# ── Step 1: Index genome ──────────────────────────────────────
IDX="${SCRATCH}/ATCC52814"
if [[ ! -f "${IDX}.1.bt2" ]]; then
    echo "[$(date)] Indexing genome with bowtie2-build..."
    bowtie2-build --threads "${THREADS}" "${GENOME_FA}" "${IDX}"
fi
if [[ ! -f "${GENOME_FA}.fai" ]]; then
    echo "[$(date)] Indexing genome with samtools faidx..."
    samtools faidx "${GENOME_FA}"
fi
echo "[$(date)] Genome index ready."

# ── Step 2: Trim ──────────────────────────────────────────────
echo "[$(date)] Trimming IP reads with fastp..."
fastp \
  --in1 "${IP_R1}" --in2 "${IP_R2}" \
  --out1 "${SCRATCH}/IP_R1_clean.fq.gz" \
  --out2 "${SCRATCH}/IP_R2_clean.fq.gz" \
  --thread "${THREADS}" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 36 \
  --json "${OUT_DIR}/fastp_IP.json" \
  --html "${OUT_DIR}/fastp_IP.html"

echo "[$(date)] Trimming Input reads with fastp..."
fastp \
  --in1 "${IN_R1}" --in2 "${IN_R2}" \
  --out1 "${SCRATCH}/Input_R1_clean.fq.gz" \
  --out2 "${SCRATCH}/Input_R2_clean.fq.gz" \
  --thread "${THREADS}" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 36 \
  --json "${OUT_DIR}/fastp_Input.json" \
  --html "${OUT_DIR}/fastp_Input.html"

# ── Step 3: Align ─────────────────────────────────────────────
echo "[$(date)] Aligning IP reads..."
bowtie2 \
  -x "${IDX}" \
  -1 "${SCRATCH}/IP_R1_clean.fq.gz" \
  -2 "${SCRATCH}/IP_R2_clean.fq.gz" \
  --very-sensitive \
  --no-mixed --no-discordant \
  -p "${THREADS}" \
  2>"${OUT_DIR}/bowtie2_IP_stats.txt" \
| samtools sort -@ "${THREADS}" -m 2G -T "${SCRATCH}/sort_IP" -o "${SCRATCH}/IP.bam"
samtools index "${SCRATCH}/IP.bam"

echo "[$(date)] Aligning Input reads..."
bowtie2 \
  -x "${IDX}" \
  -1 "${SCRATCH}/Input_R1_clean.fq.gz" \
  -2 "${SCRATCH}/Input_R2_clean.fq.gz" \
  --very-sensitive \
  --no-mixed --no-discordant \
  -p "${THREADS}" \
  2>"${OUT_DIR}/bowtie2_Input_stats.txt" \
| samtools sort -@ "${THREADS}" -m 2G -T "${SCRATCH}/sort_Input" -o "${SCRATCH}/Input.bam"
samtools index "${SCRATCH}/Input.bam"

# ── Step 3b: Add read groups (required by Picard MarkDuplicates) ──────────
echo "[$(date)] Adding read groups to IP/Input BAMs..."
samtools addreplacerg \
  -@ "${THREADS}" \
  -r "ID:IP" -r "SM:H3K9me3_IP_ATCC52814" -r "LB:H3K9me3" -r "PL:ILLUMINA" -r "PU:IP_WT" \
  -o "${SCRATCH}/IP.rg.bam" \
  "${SCRATCH}/IP.bam"
samtools index "${SCRATCH}/IP.rg.bam"

samtools addreplacerg \
  -@ "${THREADS}" \
  -r "ID:INPUT" -r "SM:H3K9me3_Input_ATCC52814" -r "LB:H3K9me3" -r "PL:ILLUMINA" -r "PU:Input_WT" \
  -o "${SCRATCH}/Input.rg.bam" \
  "${SCRATCH}/Input.bam"
samtools index "${SCRATCH}/Input.rg.bam"

# ── Step 4: Mark duplicates ───────────────────────────────────
echo "[$(date)] Marking duplicates for IP..."
picard MarkDuplicates \
  I="${SCRATCH}/IP.rg.bam" \
  O="${OUT_DIR}/IP_dedup.bam" \
  M="${OUT_DIR}/IP_dup_metrics.txt" \
  REMOVE_DUPLICATES=true \
  VALIDATION_STRINGENCY=LENIENT
samtools index "${OUT_DIR}/IP_dedup.bam"

echo "[$(date)] Marking duplicates for Input..."
picard MarkDuplicates \
  I="${SCRATCH}/Input.rg.bam" \
  O="${OUT_DIR}/Input_dedup.bam" \
  M="${OUT_DIR}/Input_dup_metrics.txt" \
  REMOVE_DUPLICATES=true \
  VALIDATION_STRINGENCY=LENIENT
samtools index "${OUT_DIR}/Input_dedup.bam"

# ── Step 5: Alignment stats ───────────────────────────────────
echo "[$(date)] Flagstat IP:"
samtools flagstat -@ "${THREADS}" "${OUT_DIR}/IP_dedup.bam" | tee "${OUT_DIR}/IP_dedup_flagstat.txt"
echo "[$(date)] Flagstat Input:"
samtools flagstat -@ "${THREADS}" "${OUT_DIR}/Input_dedup.bam" | tee "${OUT_DIR}/Input_dedup_flagstat.txt"

# ── Step 6: bigwig coverage (CPM, bin 10bp) ───────────────────
echo "[$(date)] Generating bigwig for IP (CPM, bin=10bp)..."
bamCoverage \
  -b "${OUT_DIR}/IP_dedup.bam" \
  -o "${OUT_DIR}/H3K9me3_IP_CPM10.bigwig" \
  --normalizeUsing CPM \
  --binSize 10 \
  --smoothLength 50 \
  --extendReads \
  -p "${THREADS}"

echo "[$(date)] Generating bigwig for Input (CPM, bin=10bp)..."
bamCoverage \
  -b "${OUT_DIR}/Input_dedup.bam" \
  -o "${OUT_DIR}/H3K9me3_Input_CPM10.bigwig" \
  --normalizeUsing CPM \
  --binSize 10 \
  --smoothLength 50 \
  --extendReads \
  -p "${THREADS}"

# ── Step 7: log2 IP/Input ratio bigwig ───────────────────────
echo "[$(date)] Generating log2 IP/Input ratio bigwig (bin=50bp)..."
bamCompare \
  -b1 "${OUT_DIR}/IP_dedup.bam" \
  -b2 "${OUT_DIR}/Input_dedup.bam" \
  -o "${OUT_DIR}/H3K9me3_log2ratio.bigwig" \
  --operation log2 \
  --normalizeUsing CPM \
  --binSize 50 \
  --smoothLength 150 \
  --extendReads \
  --pseudocount 1 \
  -p "${THREADS}"

# ── Step 8: Peak calling (MACS2 broad) ────────────────────────
echo "[$(date)] Calling peaks with MACS2 (broad mode)..."
macs2 callpeak \
  -t "${OUT_DIR}/IP_dedup.bam" \
  -c "${OUT_DIR}/Input_dedup.bam" \
  -f BAMPE \
  --broad \
  -g "${GENOME_SIZE}" \
  --broad-cutoff 0.1 \
  -n H3K9me3_ATCC52814 \
  --outdir "${OUT_DIR}" \
  --nomodel \
  --extsize 200 \
  2>"${OUT_DIR}/macs2_broad.log"

# Also call narrow peaks for comparison
echo "[$(date)] Calling narrow peaks with MACS2..."
macs2 callpeak \
  -t "${OUT_DIR}/IP_dedup.bam" \
  -c "${OUT_DIR}/Input_dedup.bam" \
  -f BAMPE \
  -g "${GENOME_SIZE}" \
  -n H3K9me3_ATCC52814_narrow \
  --outdir "${OUT_DIR}" \
  --nomodel \
  --extsize 200 \
  2>"${OUT_DIR}/macs2_narrow.log"

# ── Step 9: Cleanup scratch ───────────────────────────────────
echo "[$(date)] Cleaning up scratch..."
rm -rf "${SCRATCH}"

echo "[$(date)] === Pipeline complete ==="
echo "Outputs in: ${OUT_DIR}"
echo "  bigwigs:   H3K9me3_IP_CPM10.bigwig, H3K9me3_Input_CPM10.bigwig, H3K9me3_log2ratio.bigwig"
echo "  peaks:     H3K9me3_ATCC52814_peaks.broadPeak, H3K9me3_ATCC52814_narrow_peaks.narrowPeak"
echo "  BAMs:      IP_dedup.bam, Input_dedup.bam"
