#!/bin/bash
# ============================================================
# run_hmgb3_jockey_cross_species.sh
# Cross-species conservation of HMGb3 binding site in
# LINE/I-Jockey elements across MU/BL species.
#
# Rhimi59 reference: 4/4 elements (100%) have HMGb3 site (p=0.00274)
# ============================================================
#SBATCH --job-name=hmgb3_jockey
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/hmgb3_jockey_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/hmgb3_jockey_%j.err
#SBATCH --exclude=LEGcompute1

set -euo pipefail

module load MEME/5.5.7-gompi-2023b
module load BEDTools/2.31.1-GCC-13.3.0

ANNOT_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
GENOMES_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
MOTIF="/legserv/Temp/Carlos/EDF-TE/data/TFsDAP/motifs/HMGb3.meme"
PARSER="/data/carlos/EDF-TE/EDF-TE/TF_TE_analysis/scripts/hmgb3_jockey_fimo_parser.py"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/hmgb3_jockey_cross_species"
SCRATCH="/scratch/hmgb3_jockey_${SLURM_JOB_ID}"

mkdir -p "${OUT_DIR}" "${SCRATCH}"
mkdir -p /data/carlos/EDF-TE/logs

SUMMARY="${OUT_DIR}/hmgb3_jockey_summary.tsv"
echo -e "species\tgenome_mb\tn_jockey_elements\tn_with_hmgb3\tfraction_with_hmgb3\tbest_score\tmean_score\tmin_qval" \
    > "${SUMMARY}"

echo "[$(date)] Starting HMGb3 x LINE/I-Jockey cross-species analysis"

for BED_PATH in "${ANNOT_DIR}"/MU_*_EarlGrey/*_summaryFiles/*.filteredRepeats.bed \
                "${ANNOT_DIR}"/BL_*_EarlGrey/*_summaryFiles/*.filteredRepeats.bed; do

    SP=$(basename "${BED_PATH}" .filteredRepeats.bed)
    GENOME="${GENOMES_DIR}/${SP}.fasta"
    [[ ! -f "${GENOME}" ]] && { echo "[SKIP] ${SP}: no genome"; continue; }

    JOCKEY_BED="${SCRATCH}/${SP}_jockey.bed"
    grep "LINE/I-Jockey" "${BED_PATH}" > "${JOCKEY_BED}" || true
    N_ELEM=$(wc -l < "${JOCKEY_BED}")
    [[ "${N_ELEM}" -eq 0 ]] && { rm -f "${JOCKEY_BED}"; continue; }

    echo "[$(date)] ${SP}: ${N_ELEM} LINE/I-Jockey elements"
    GENOME_MB=$(awk '/^>/{next}{s+=length($0)} END{printf "%.2f", s/1e6}' "${GENOME}")

    JOCKEY_FA="${SCRATCH}/${SP}_jockey.fa"
    bedtools getfasta -fi "${GENOME}" -bed "${JOCKEY_BED}" -fo "${JOCKEY_FA}" -name 2>/dev/null

    N_SEQS=$(grep -c "^>" "${JOCKEY_FA}" 2>/dev/null || echo 0)
    if [[ "${N_SEQS}" -eq 0 ]]; then
        echo "[WARN] ${SP}: no sequences extracted"
        rm -f "${JOCKEY_BED}" "${JOCKEY_FA}"
        continue
    fi

    FIMO_DIR="${OUT_DIR}/fimo_${SP}"
    mkdir -p "${FIMO_DIR}"
    fimo --thresh 1e-3 --oc "${FIMO_DIR}" --verbosity 1 \
        "${MOTIF}" "${JOCKEY_FA}" \
        2>"${FIMO_DIR}/fimo.stderr" || true

    FIMO_TSV="${FIMO_DIR}/fimo.tsv"
    if [[ -f "${FIMO_TSV}" ]]; then
        python3 "${PARSER}" "${SP}" "${GENOME_MB}" "${N_ELEM}" "${FIMO_TSV}" "${SUMMARY}"
    else
        echo "  ${SP}: no FIMO hits"
        printf '%s\t%s\t%d\t0\t0.0000\t0.000\t0.000\t1.00e+00\n' \
            "${SP}" "${GENOME_MB}" "${N_ELEM}" >> "${SUMMARY}"
    fi

    rm -f "${JOCKEY_BED}" "${JOCKEY_FA}"
done

echo ""
echo "[$(date)] === Done. Summary: ${SUMMARY} ($(wc -l < "${SUMMARY}") rows) ==="
echo ""
echo "=== Species with >=50% LINE/I-Jockey elements having HMGb3 site ==="
awk -F'\t' 'NR>1 && $5+0>=0.5 {printf "%s\t%s elements\t%s with site\t%.1f%%\n",$1,$3,$4,$5*100}' "${SUMMARY}"

echo ""
echo "=== Full results ==="
column -t "${SUMMARY}"

rm -rf "${SCRATCH}"
