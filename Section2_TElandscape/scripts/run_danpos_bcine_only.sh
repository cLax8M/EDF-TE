#!/bin/bash
# Just the DANPOS step on the already-pooled Bcine1 BAM.
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate nucseq_pipeline

SP="AS_Bcine1"
MERGED="/scratch/invivo_mnase_3812637/${SP}_dedup.bam"
SP_OUT="/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/${SP}"
mkdir -p "${SP_OUT}/${SP}_danpos"

if [[ ! -s "${MERGED}" ]]; then
  echo "[ERR] missing pooled BAM: ${MERGED}"; exit 1
fi

echo "[danpos] dpos on ${MERGED}"
DANPOS3="/data/carlos/EDF-TE/TE_analysis/DANPOS3/danpos.py"
python "${DANPOS3}" dpos "${MERGED}" -m 1 -a 1 -o "${SP_OUT}/${SP}_danpos"

echo "[DONE] ${SP}"
ls -lh "${SP_OUT}/${SP}_danpos"
