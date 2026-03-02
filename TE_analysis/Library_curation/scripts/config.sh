#!/bin/bash
# Central config for TE library curation pipeline

WORK_DIR="/data/carlos/EDF-TE/TE_analysis/Library_curation"
EARlgrey_OUT="/data/carlos/EDF-TE/TE_analysis"
DB_DIR="${WORK_DIR}/databases"

REFSEQ_DB="${DB_DIR}/refseq_eukaryota"          # Diamond .dmnd (no extension)
REPEATPEPS="${DB_DIR}/RepeatPeps.lib"
REPET_HMM="${DB_DIR}/ProfilesBankForREPET_Pfam35.0_GypsyDB.hmm"

THREADS=16
MIN_LEN=300
MAX_LEN=15000
MAX_N_FRAC=0.05                                  # 5% N content threshold
