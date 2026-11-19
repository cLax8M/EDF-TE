#!/usr/bin/env bash

# prepares per-clade TE libraries and orthogroups for HTT analyses
set -euo pipefail

# load config
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# define dirs
HTT_DIR="${RESULTS_DIR}/htt"
LIB_DIR="${HTT_DIR}/clade_libraries"
OG_DIR="${HTT_DIR}/orthogroups"
LOG_DIR="${HTT_DIR}/logs"
mkdir -p "${LIB_DIR}" "${OG_DIR}" "${LOG_DIR}"

# define programs
MMSEQS="${TOOLS_DIR}/mmseqs/bin/mmseqs"
REPEATCLASSIFIER="${TOOLS_DIR}/RepeatModeler/RepeatClassifier"
ORTHOFINDER="${TOOLS_DIR}/OrthoFinder/orthofinder"


# generate per-clade TE libraries

echo "generating per-clade TE libraries"

# load species per clade
declare -A CLADE_SPECIES
while IFS=$'\t' read -r species clade; do
    CLADE_SPECIES["${clade}"]+=" ${species}"
done < "${CLADE_MAP}"

# concatenate TE libraries
for clade in "${!CLADE_SPECIES[@]}"; do
    out_fa="${LIB_DIR}/${clade}.raw.fa"
    : > "${out_fa}"
    for sp in ${CLADE_SPECIES["${clade}"]}; do
        lib="${TE_LIB_DIR}/${sp}/${sp}-families.fa"
        if [[ -f "${lib}" ]]; then
            # add species prefix to sequence headers
            awk -v sp="${sp}" '/^>/{print ">" sp "__" substr($0,2)} !/^>/{print}' \
                "${lib}" >> "${out_fa}"
        fi
    done
    n_seqs="$(grep -c '^>' "${out_fa}" 2>/dev/null || echo 0)"
    echo "  ${clade}: ${n_seqs} raw sequences"
done


# cluster te libraries

echo "clustering TE libraries"

for clade in "${!CLADE_SPECIES[@]}"; do
    raw_fa="${LIB_DIR}/${clade}.raw.fa"
    [[ -s "${raw_fa}" ]] || continue
    out_fa="${LIB_DIR}/${clade}.fa"
    tmp_clust="${LIB_DIR}/tmp_${clade}"
    mkdir -p "${tmp_clust}"

    sbatch \
        --job-name="mmseqs_${clade}" \
        --partition="${SLURM_PARTITION}" \
        --cpus-per-task=8 \
        --mem=16G \
        --time=04:00:00 \
        --output="${LOG_DIR}/mmseqs_${clade}.out" \
        --wrap="
            micromamba run -n edf_py \
                ${MMSEQS} easy-cluster \
                    '${raw_fa}' \
                    '${tmp_clust}/cluster' \
                    '${tmp_clust}' \
                    --min-seq-id 0.95 \
                    -c 0.80 \
                    --cov-mode 1 \
                    --threads 8 \
                    --cluster-reassign 1
            cp '${tmp_clust}/cluster_rep_seq.fasta' '${out_fa}'
        "
done

echo "submitted clustering jobs"


# classify TE libraries

echo "classifying TE libraries"

for clade in "${!CLADE_SPECIES[@]}"; do
    clust_fa="${LIB_DIR}/${clade}.fa"
    [[ -s "${clust_fa}" ]] || continue
    sbatch \
        --job-name="rclass_${clade}" \
        --partition="${SLURM_PARTITION}" \
        --cpus-per-task=4 \
        --mem=8G \
        --time=06:00:00 \
        --dependency=singleton \
        --output="${LOG_DIR}/rclass_${clade}.out" \
        --wrap="
            cd '${LIB_DIR}'
            micromamba run -n earlgrey \
                ${REPEATCLASSIFIER} \
                    -consensi '${clust_fa}' \
                    -engine rmblast
        "
done


# run orthofinder per clade

echo "running orthofinder"

for clade in "${!CLADE_SPECIES[@]}"; do
    prot_dir="${OG_DIR}/${clade}/proteins"
    mkdir -p "${prot_dir}"
    for sp in ${CLADE_SPECIES["${clade}"]}; do
        prot_fa="${PROTEOME_DIR}/${sp}.faa"
        [[ -f "${prot_fa}" ]] && ln -sf "${prot_fa}" "${prot_dir}/${sp}.faa"
    done

    n_sp="$(ls "${prot_dir}"/*.faa 2>/dev/null | wc -l)"
    [[ "${n_sp}" -lt 3 ]] && continue

    sbatch \
        --job-name="of_${clade}" \
        --partition="${SLURM_PARTITION_LONG}" \
        --cpus-per-task=16 \
        --mem=32G \
        --time=48:00:00 \
        --output="${LOG_DIR}/orthofinder_${clade}.out" \
        --wrap="
            micromamba run -n phylo \
                ${ORTHOFINDER} \
                    -f '${prot_dir}' \
                    -t 16 \
                    -a 4 \
                    -o '${OG_DIR}/${clade}/of_output'
        "
done

echo "all jobs submitted"
echo "monitor with: squeue -u ${USER}"
