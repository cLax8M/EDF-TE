#!/bin/bash
# ============================================================
# run_invivo_mnase_array.sh
# Multi-replicate MNase-seq pipeline for in vivo nucleosome data:
#   - AS_Ncrass2  (N. crassa OR74A; SRR13067301-303)
#   - AS_Bcine1   (B. cinerea B05.10 MAINE-seq; SRR13224644/47/48)
#
# Pipeline per species: download all reps -> BBDuk QC -> bowtie1 align ->
#   Picard MarkDuplicates -> samtools merge dedup BAMs -> DANPOS dpos.
#
# Output: /data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/<SP>/<SP>_danpos/pooled/<SP>_dedup.smooth.wig
# (Same path/format as the previous 5 species, so downstream scripts work unchanged.)
#
# Run interactively:
#   SLURM_ARRAY_TASK_ID=1 bash run_invivo_mnase_array.sh
#   SLURM_ARRAY_TASK_ID=2 bash run_invivo_mnase_array.sh
# ============================================================
#SBATCH --job-name=invivo_mnase
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=1-2
#SBATCH --output=/data/carlos/EDF-TE/logs/nucseq/invivo_mnase_%A_%a.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nucseq/invivo_mnase_%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate nucseq_pipeline
# wget+pigz come from earlgrey env (APPENDED so nucseq_pipeline python wins for danpos)
export PATH="${PATH}:/home/carlos/micromamba/envs/earlgrey/bin"

GENOMES_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
OUT_DIR="/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq"
TMPDIR_BASE="/scratch/invivo_mnase_$$"
mkdir -p "${OUT_DIR}" "${TMPDIR_BASE}"
TMPDIR="${TMPDIR_BASE}"
NCPU="${SLURM_CPUS_PER_TASK:-8}"

# --- species manifest (multiple SRRs per species, comma-separated) ----
declare -A SRRS=(
  [1]="SRR13067301,SRR13067302,SRR13067303"
  [2]="SRR13224644,SRR13224647,SRR13224648"
)
declare -A SP=(
  [1]="AS_Ncrass2"
  [2]="AS_Bcine1"
)

IDX="${SLURM_ARRAY_TASK_ID:?need SLURM_ARRAY_TASK_ID}"
SPECIES="${SP[$IDX]}"
SRR_LIST="${SRRS[$IDX]}"
GENOME_FA="${GENOMES_DIR}/${SPECIES}.fasta"
SP_OUT="${OUT_DIR}/${SPECIES}"
mkdir -p "${SP_OUT}"

OUT_WIG="${SP_OUT}/${SPECIES}_danpos/pooled/${SPECIES}_dedup.smooth.wig"
if [[ -s "${OUT_WIG}" ]]; then
  echo "[SKIP] ${SPECIES} already complete: ${OUT_WIG}"; exit 0
fi

echo "[START] ${SPECIES}  SRRs: ${SRR_LIST}"

# --- helper: ENA path for an SRR ---
ena_url() {
  local s="$1"
  local pre="${s:0:6}"
  local len=${#s}
  if [[ $len -le 9 ]]; then
    echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${pre}/${s}"
  else
    local extra=$(( len - 9 ))
    local suf
    suf=$(printf '%03d' "${s: -${extra}}")
    echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${pre}/${suf}/${s}"
  fi
}

# --- bowtie1 index (built once per species) ---
INDEX="${TMPDIR}/${SPECIES}_bt1"
echo "[idx] building bowtie1 index"
bowtie-build --threads "${NCPU}" "${GENOME_FA}" "${INDEX}"

DEDUPED_BAMS=()
IFS=',' read -ra SRR_ARR <<<"${SRR_LIST}"
for SRR_ID in "${SRR_ARR[@]}"; do
  echo "==== ${SPECIES}  rep ${SRR_ID} ===="
  REP_TMP="${TMPDIR}/${SRR_ID}"
  mkdir -p "${REP_TMP}"

  # 1) download from ENA
  base=$(ena_url "${SRR_ID}")
  ok=false
  if wget -q -P "${REP_TMP}" "${base}/${SRR_ID}_1.fastq.gz" && [[ -s "${REP_TMP}/${SRR_ID}_1.fastq.gz" ]]; then
    wget -q -P "${REP_TMP}" "${base}/${SRR_ID}_2.fastq.gz" || true
    ok=true
  elif wget -q -P "${REP_TMP}" "${base}/${SRR_ID}.fastq.gz" && [[ -s "${REP_TMP}/${SRR_ID}.fastq.gz" ]]; then
    ok=true
  fi
  if [[ "${ok}" != "true" ]]; then
    echo "[ENA fail] falling back to prefetch+fasterq-dump for ${SRR_ID}"
    prefetch -O "${REP_TMP}" "${SRR_ID}"
    fasterq-dump --threads "${NCPU}" -O "${REP_TMP}" -t "${REP_TMP}" "${REP_TMP}/${SRR_ID}/${SRR_ID}.sra"
    [[ -f "${REP_TMP}/${SRR_ID}_1.fastq" ]] && pigz "${REP_TMP}/${SRR_ID}_1.fastq" && pigz "${REP_TMP}/${SRR_ID}_2.fastq"
    [[ -f "${REP_TMP}/${SRR_ID}.fastq"   ]] && pigz "${REP_TMP}/${SRR_ID}.fastq"
  fi

  # decompress
  for fq in "${REP_TMP}"/*.fastq.gz; do
    [[ -e "${fq}" ]] || continue
    pigz -d "${fq}" 2>/dev/null || gunzip "${fq}"
  done
  R1="${REP_TMP}/${SRR_ID}_1.fastq"
  R2="${REP_TMP}/${SRR_ID}_2.fastq"
  SE="${REP_TMP}/${SRR_ID}.fastq"
  PAIRED=false
  if [[ -f "${R1}" && -f "${R2}" ]]; then PAIRED=true
  elif [[ -f "${SE}" ]]; then PAIRED=false
  else echo "[ERR] no fastq for ${SRR_ID}"; exit 1
  fi

  # 2) BBDuk QC
  ADAPTERS_REF=$(ls "${CONDA_PREFIX}/opt/bbmap"*/resources/adapters.fa 2>/dev/null | head -1)
  [[ -z "${ADAPTERS_REF}" ]] && ADAPTERS_REF="adapters"
  # auto minlen: 25 if reads are short (<=50bp), else 41
  RLEN=$(awk 'NR==2 {print length($0); exit}' "${R1}" 2>/dev/null || echo 100)
  if [[ "${RLEN}" -le 50 ]]; then MINLEN=25; else MINLEN=41; fi
  echo "[bbduk] read length=${RLEN}, using minlen=${MINLEN}"
  if [[ "${PAIRED}" == "true" ]]; then
    bbduk.sh -Xmx16g \
      in1="${R1}" in2="${R2}" \
      out1="${REP_TMP}/clean_R1.fastq" out2="${REP_TMP}/clean_R2.fastq" \
      ref="${ADAPTERS_REF}" ktrim=r k=25 mink=11 hdist=1 \
      qtrim=r trimq=6 minlen=${MINLEN} minavgquality=13 maxns=0 \
      threads="${NCPU}"
  else
    bbduk.sh -Xmx16g \
      in="${SE}" out="${REP_TMP}/clean.fastq" \
      ref="${ADAPTERS_REF}" ktrim=r k=25 mink=11 hdist=1 \
      qtrim=r trimq=6 minlen=${MINLEN} minavgquality=13 maxns=0 \
      threads="${NCPU}"
  fi

  # 3) bowtie1 align
  if [[ "${PAIRED}" == "true" ]]; then
    bowtie -S -t -m 1 -v 2 -p "${NCPU}" "${INDEX}" \
      -1 "${REP_TMP}/clean_R1.fastq" -2 "${REP_TMP}/clean_R2.fastq" \
      | samtools sort -@ "${NCPU}" -o "${REP_TMP}/${SRR_ID}.bam"
  else
    bowtie -S -t -m 1 -v 2 -p "${NCPU}" "${INDEX}" "${REP_TMP}/clean.fastq" \
      | samtools sort -@ "${NCPU}" -o "${REP_TMP}/${SRR_ID}.bam"
  fi
  samtools index "${REP_TMP}/${SRR_ID}.bam"

  # 4) Picard MarkDuplicates
  picard MarkDuplicates \
    I="${REP_TMP}/${SRR_ID}.bam" \
    O="${REP_TMP}/${SRR_ID}_dedup.bam" \
    M="${SP_OUT}/${SRR_ID}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
  samtools index "${REP_TMP}/${SRR_ID}_dedup.bam"
  DEDUPED_BAMS+=("${REP_TMP}/${SRR_ID}_dedup.bam")

  # cleanup intermediate fastqs to save space
  rm -f "${REP_TMP}"/*.fastq
done

# --- merge replicate BAMs ---
echo "[merge] pooling ${#DEDUPED_BAMS[@]} replicate BAMs"
MERGED="${TMPDIR}/${SPECIES}_dedup.bam"
samtools merge -@ "${NCPU}" -f "${MERGED}" "${DEDUPED_BAMS[@]}"
samtools index "${MERGED}"

# --- DANPOS3 ---
echo "[danpos] calling nucleosomes"
DANPOS3="/data/carlos/EDF-TE/TE_analysis/DANPOS3/danpos.py"
mkdir -p "${SP_OUT}/${SPECIES}_danpos"
python "${DANPOS3}" dpos "${MERGED}" -m 1 -a 1 -o "${SP_OUT}/${SPECIES}_danpos"

echo "[DONE] ${SPECIES} -> ${OUT_WIG}"
rm -rf "${TMPDIR_BASE}"
