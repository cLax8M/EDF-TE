#!/usr/bin/env bash
# =============================================================================
# 01_htt_phylo_prep.sh  —  Prepare per-clade TE libraries for HTT detection
#
# Constructs per-clade TE libraries and per-species OrthoFinder single-copy
# orthogroup alignments required for horizontal TE transfer (HTT) candidate
# detection.
#
# Steps:
#   1. Concatenate EarlGrey TE consensus sequences per clade
#   2. Cluster within-clade sequences (MMseqs2, 95% identity) to remove
#      redundancy arising from closely related species
#   3. Annotate remaining sequences with RepeatClassifier
#   4. Run OrthoFinder on single-copy orthogroups per target clade for
#      per-node Ks dating (used downstream in 03_htt_ks_dating.py)
#
# Usage:
#   source config.sh
#   bash 01_htt_phylo_prep.sh
#
# Outputs:
#   ${HTT_DIR}/clade_libraries/<clade>.fa        – clustered TE library
#   ${HTT_DIR}/clade_libraries/<clade>.classify  – RepeatClassifier annotations
#   ${HTT_DIR}/orthogroups/<clade>/              – OrthoFinder output
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

HTT_DIR="${RESULTS_DIR}/htt"
LIB_DIR="${HTT_DIR}/clade_libraries"
OG_DIR="${HTT_DIR}/orthogroups"
LOG_DIR="${HTT_DIR}/logs"
mkdir -p "${LIB_DIR}" "${OG_DIR}" "${LOG_DIR}"

MMSEQS="${TOOLS_DIR}/mmseqs/bin/mmseqs"
REPEATCLASSIFIER="${TOOLS_DIR}/RepeatModeler/RepeatClassifier"
ORTHOFINDER="${TOOLS_DIR}/OrthoFinder/orthofinder"

# ---------------------------------------------------------------------------
# Step 1: per-clade TE library concatenation
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step 1: Concatenating per-clade TE libraries"

declare -A CLADE_SPECIES
while IFS=$'\t' read -r species clade; do
    CLADE_SPECIES["${clade}"]+=" ${species}"
done < "${CLADE_MAP}"

for clade in "${!CLADE_SPECIES[@]}"; do
    out_fa="${LIB_DIR}/${clade}.raw.fa"
    : > "${out_fa}"
    for sp in ${CLADE_SPECIES["${clade}"]}; do
        lib="${TE_LIB_DIR}/${sp}/${sp}-families.fa"
        if [[ -f "${lib}" ]]; then
            # Prefix sequence headers with species code
            awk -v sp="${sp}" '/^>/{print ">" sp "__" substr($0,2)} !/^>/{print}' \
                "${lib}" >> "${out_fa}"
        fi
    done
    n_seqs="$(grep -c '^>' "${out_fa}" 2>/dev/null || echo 0)"
    echo "  ${clade}: ${n_seqs} raw sequences"
done

# ---------------------------------------------------------------------------
# Step 2: Within-clade clustering (MMseqs2, 95% id, 80% coverage)
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step 2: Within-clade clustering"

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

echo "[$(date '+%H:%M:%S')] Clustering jobs submitted"

# ---------------------------------------------------------------------------
# Step 3: RepeatClassifier annotation
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step 3: RepeatClassifier"

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

# ---------------------------------------------------------------------------
# Step 4: Per-clade OrthoFinder for downstream Ks dating
# ---------------------------------------------------------------------------
echo "[$(date '+%H:%M:%S')] Step 4: OrthoFinder per clade"

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

echo "[$(date '+%H:%M:%S')] All jobs submitted.  Monitor with: squeue -u \${USER}"
