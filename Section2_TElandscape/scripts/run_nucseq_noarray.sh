#!/bin/bash
# ============================================================
# run_nucseq_noarray.sh
# Sequential (non-array) run of the 5-species Nuc-seq pipeline.
# ============================================================

#SBATCH --job-name=nucseq_noarray
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --output=/data/carlos/EDF-TE/logs/nucseq/nucseq_noarray_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/nucseq/nucseq_noarray_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate nucseq_pipeline

GENOMES_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
OUT_DIR="/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq"
mkdir -p "${OUT_DIR}" /data/carlos/EDF-TE/logs/nucseq

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

run_one() {
    local idx="$1"
    local srr_id="${SRR[$idx]}"
    local species="${SP[$idx]}"
    local genome_fa="${GENOME[$idx]}"
    local sp_out="${OUT_DIR}/${species}"
    local tmpbase="/scratch/nucseq_${SLURM_JOB_ID}_${idx}"
    local tmpdir="${tmpbase}"

    mkdir -p "${sp_out}" "${tmpbase}"

    if [[ -f "${sp_out}/${species}_danpos/pooled/${species}_dedup.Frag.wig" ]]; then
        echo "[SKIP] ${species} already complete"
        rm -rf "${tmpbase}"
        return 0
    fi

    echo "[START] ${species} (${srr_id})"
    prefetch "${srr_id}" --output-directory "${tmpdir}"
    fasterq-dump "${tmpdir}/${srr_id}" --outdir "${tmpdir}" --threads "${SLURM_CPUS_PER_TASK}" --split-files

    local r1="${tmpdir}/${srr_id}_1.fastq"
    local r2="${tmpdir}/${srr_id}_2.fastq"
    [[ ! -f "${r2}" ]] && r2=""

    local adapters_ref
    adapters_ref=$(ls "${CONDA_PREFIX}"/opt/bbmap*/resources/adapters.fa 2>/dev/null | head -1)
    [[ -z "${adapters_ref}" ]] && adapters_ref="adapters"

    local clean_r1 clean_r2
    if [[ -n "${r2}" ]]; then
        bbduk.sh in1="${r1}" in2="${r2}" out1="${tmpdir}/${species}_R1_clean.fastq" out2="${tmpdir}/${species}_R2_clean.fastq" ref="${adapters_ref}" ktrim=r k=25 mink=11 hdist=1 qtrim=r trimq=6 minlen=41 minavgquality=13 maxns=0 threads="${SLURM_CPUS_PER_TASK}"
        clean_r1="${tmpdir}/${species}_R1_clean.fastq"
        clean_r2="${tmpdir}/${species}_R2_clean.fastq"
    else
        bbduk.sh in="${r1}" out="${tmpdir}/${species}_clean.fastq" ref="${adapters_ref}" ktrim=r k=25 mink=11 hdist=1 qtrim=r trimq=6 minlen=41 minavgquality=13 maxns=0 threads="${SLURM_CPUS_PER_TASK}"
        clean_r1="${tmpdir}/${species}_clean.fastq"
        clean_r2=""
    fi

    local index="${tmpdir}/${species}_bt1"
    bowtie-build --threads "${SLURM_CPUS_PER_TASK}" "${genome_fa}" "${index}"

    if [[ -n "${clean_r2}" ]]; then
        bowtie -S -t -m 1 -v 2 -p "${SLURM_CPUS_PER_TASK}" "${index}" -1 "${clean_r1}" -2 "${clean_r2}" | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${tmpdir}/${species}.bam"
    else
        bowtie -S -t -m 1 -v 2 -p "${SLURM_CPUS_PER_TASK}" "${index}" "${clean_r1}" | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${tmpdir}/${species}.bam"
    fi
    samtools index "${tmpdir}/${species}.bam"

    picard MarkDuplicates I="${tmpdir}/${species}.bam" O="${tmpdir}/${species}_dedup.bam" M="${sp_out}/${species}_dup_metrics.txt" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
    samtools index "${tmpdir}/${species}_dedup.bam"

    local danpos3="/data/carlos/EDF-TE/TE_analysis/DANPOS3/danpos.py"
    mkdir -p "${sp_out}/${species}_danpos"
    python "${danpos3}" dpos "${tmpdir}/${species}_dedup.bam" -m 1 -a 1 -o "${sp_out}/${species}_danpos"

    echo "[DONE] ${species}"
    rm -rf "${tmpbase}"
}

for idx in 1 2 3 4 5; do
    run_one "$idx"
done

echo "[SUMMARY] nucseq_noarray finished"
