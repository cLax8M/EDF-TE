#!/bin/bash
# ============================================================
# run_atac_mereb1.sh
# ATAC-seq pipeline for M. eburnea (BA_Mereb1):
#   SRR28069093 (heavy grind)
#   SRR28069094 (heavy grind + 5 micron filter)
#   SRR28069106 (light grind)
#
# Pipeline per replicate:
#   ENA/SRA fastq -> fastp (adapter+QC, 2x150 PE) ->
#   bowtie2 (--very-sensitive -X 2000 --no-mixed --no-discordant) ->
#   samtools sort + filter (q>=30, properly paired, no chrM) ->
#   Picard MarkDuplicates -> samtools merge across replicates ->
#   MACS2 callpeak (BAMPE, --shift -75 --extsize 150) ->
#   bamCoverage (binSize=10, normalizeUsing=RPGC, extendReads) -> BigWig
#
# Output:
#   /data/carlos/EDF-TE/TE_analysis/nucleosomes/atac-seq/BA_Mereb1/
#       BA_Mereb1_atac.bam              (pooled, dedup, filtered)
#       BA_Mereb1_atac_peaks.narrowPeak
#       BA_Mereb1_atac.bw               (RPGC accessibility signal)
#       BA_Mereb1_atac.wig              (converted from .bw for our pipeline)
# ============================================================
set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate chipseq_env
# wget+pigz come from earlgrey env
export PATH="/home/carlos/micromamba/envs/earlgrey/bin:${PATH}"

NCPU=16
SP="BA_Mereb1"
SRRS=(SRR28069093 SRR28069094 SRR28069106)

GENOMES_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
GENOME_FA="${GENOMES_DIR}/${SP}.fasta"
FAI="${GENOME_FA}.fai"
OUT_DIR="/data/carlos/EDF-TE/TE_analysis/nucleosomes/atac-seq/${SP}"
TMP="/scratch/atac_${SP}_$$"
mkdir -p "${OUT_DIR}" "${TMP}"

# effective genome size
EGS=$(awk '{s+=$2} END {print s}' "${FAI}")
echo "[info] genome size: ${EGS}"

OUT_BW="${OUT_DIR}/${SP}_atac.bw"
if [[ -s "${OUT_BW}" ]]; then
  echo "[SKIP] already complete: ${OUT_BW}"; exit 0
fi

# bowtie2 index
INDEX="${TMP}/${SP}_bt2"
echo "[idx] bowtie2-build"
bowtie2-build --threads "${NCPU}" "${GENOME_FA}" "${INDEX}"

ena_url() {
  local s="$1"; local pre="${s:0:6}"; local len=${#s}
  if [[ $len -le 9 ]]; then
    echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${pre}/${s}"
  else
    local extra=$(( len - 9 ))
    local suf; suf=$(printf '%03d' "${s: -${extra}}")
    echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${pre}/${suf}/${s}"
  fi
}

DEDUP_BAMS=()
for SRR in "${SRRS[@]}"; do
  echo "==== ${SRR} ===="
  REP="${TMP}/${SRR}"; mkdir -p "${REP}"

  # --- 1) download from ENA ---
  base=$(ena_url "${SRR}")
  ok=false
  if wget -q -P "${REP}" "${base}/${SRR}_1.fastq.gz" && [[ -s "${REP}/${SRR}_1.fastq.gz" ]]; then
    wget -q -P "${REP}" "${base}/${SRR}_2.fastq.gz" || true
    ok=true
  elif wget -q -P "${REP}" "${base}/${SRR}.fastq.gz" && [[ -s "${REP}/${SRR}.fastq.gz" ]]; then
    ok=true
  fi
  if [[ "${ok}" != "true" ]]; then
    echo "[ENA fail] using prefetch"
    micromamba run -n nucseq_pipeline prefetch -O "${REP}" "${SRR}"
    micromamba run -n nucseq_pipeline fasterq-dump --threads "${NCPU}" -O "${REP}" -t "${REP}" "${REP}/${SRR}/${SRR}.sra"
    pigz "${REP}"/*.fastq 2>/dev/null || true
  fi

  R1="${REP}/${SRR}_1.fastq.gz"
  R2="${REP}/${SRR}_2.fastq.gz"
  PAIRED=true
  if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
    echo "[ERR] expected PE for ATAC ${SRR}"; ls "${REP}"; exit 1
  fi

  # --- 2) fastp QC + adapter trim ---
  fastp -i "${R1}" -I "${R2}" \
        -o "${REP}/clean_R1.fq.gz" -O "${REP}/clean_R2.fq.gz" \
        -w "${NCPU}" --detect_adapter_for_pe \
        -h "${OUT_DIR}/${SRR}_fastp.html" -j "${OUT_DIR}/${SRR}_fastp.json"

  # --- 3) bowtie2 align (ATAC params) ---
  bowtie2 -x "${INDEX}" -1 "${REP}/clean_R1.fq.gz" -2 "${REP}/clean_R2.fq.gz" \
          -p "${NCPU}" --very-sensitive -X 2000 \
          --no-mixed --no-discordant 2> "${OUT_DIR}/${SRR}_bowtie2.log" \
    | samtools sort -@ "${NCPU}" -o "${REP}/${SRR}.sorted.bam"
  samtools index "${REP}/${SRR}.sorted.bam"

  # --- 4) filter: properly paired, MAPQ>=30 ---
  samtools view -@ "${NCPU}" -b -f 2 -q 30 \
          "${REP}/${SRR}.sorted.bam" \
    > "${REP}/${SRR}.filt.bam"
  samtools index "${REP}/${SRR}.filt.bam"

  # --- 5) markdup (needs MC tags from fixmate; requires name-sorted input) ---
  samtools sort -n -@ "${NCPU}" -o "${REP}/${SRR}.nsort.bam" "${REP}/${SRR}.filt.bam"
  samtools fixmate -m -@ "${NCPU}" "${REP}/${SRR}.nsort.bam" "${REP}/${SRR}.fix.bam"
  samtools sort -@ "${NCPU}" -o "${REP}/${SRR}.fixsort.bam" "${REP}/${SRR}.fix.bam"
  samtools markdup -@ "${NCPU}" -r "${REP}/${SRR}.fixsort.bam" "${REP}/${SRR}.dedup.bam"
  samtools index "${REP}/${SRR}.dedup.bam"
  rm -f "${REP}/${SRR}.nsort.bam" "${REP}/${SRR}.fix.bam" "${REP}/${SRR}.fixsort.bam"
  DEDUP_BAMS+=("${REP}/${SRR}.dedup.bam")

  rm -f "${REP}"/*.fastq.gz "${REP}/clean_"* "${REP}/${SRR}.sorted.bam"* "${REP}/${SRR}.filt.bam"*
done

# --- 6) merge replicates ---
MERGED="${OUT_DIR}/${SP}_atac.bam"
samtools merge -@ "${NCPU}" -f "${MERGED}" "${DEDUP_BAMS[@]}"
samtools index "${MERGED}"

# --- 7) MACS2 peaks (ATAC mode) ---
echo "[macs2]"
macs2 callpeak -t "${MERGED}" -f BAMPE -g "${EGS}" \
   -n "${SP}_atac" --outdir "${OUT_DIR}" \
   --nomodel --shift -75 --extsize 150 --keep-dup all -q 0.05

# --- 8) bamCoverage -> bigWig ---
echo "[deeptools] bigwig"
bamCoverage -b "${MERGED}" -o "${OUT_BW}" \
            --binSize 10 --normalizeUsing RPGC \
            --effectiveGenomeSize "${EGS}" \
            --extendReads --numberOfProcessors "${NCPU}"

# --- 9) bigWig -> wig (variableStep, span=10) for our pipeline ---
echo "[bw->wig]"
WIG="${OUT_DIR}/${SP}_atac.wig"
python - "${OUT_BW}" "${WIG}" <<'PY'
import sys, pyBigWig
bw = pyBigWig.open(sys.argv[1])
out = open(sys.argv[2], "w")
for chrom, length in bw.chroms().items():
    out.write(f"variableStep chrom={chrom} span=10\n")
    vals = bw.values(chrom, 0, length, numpy=True)
    for i in range(0, length, 10):
        v = vals[i]
        if v != v: continue
        out.write(f"{i+1}\t{v:.4f}\n")
out.close()
bw.close()
PY

echo "[DONE] ${SP}"
ls -lh "${OUT_DIR}"
rm -rf "${TMP}"
