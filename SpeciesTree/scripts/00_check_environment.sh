#!/usr/bin/env bash
# =============================================================================
# 00_check_environment.sh — Verify all required tools and their versions
#
# PURPOSE : Before running the pipeline, confirm that the conda environment
#           'esptreeEDFTE' is active and all required tools are available
#           at the expected versions. Exits with error if any check fails.
#
# USAGE   : micromamba activate esptreeEDFTE
#           bash 00_check_environment.sh
#
# PROJECT : EDF-TE — TE Evolution in Early-Diverging Fungi
# AUTHORS : Carlos Lax-Molina
# CREATED : 2026-02-19
# =============================================================================

set -euo pipefail

SCRIPT_DIR="/data/carlos/EDF-TE/EDF-TE/SpeciesTree/scripts"
source "${SCRIPT_DIR}/config.sh"

PASS=0
FAIL=0

check_tool() {
    local tool="$1"
    local version_flag="$2"
    local expected_pattern="$3"

    if ! command -v "${tool}" &>/dev/null; then
        echo "[FAIL] ${tool}: NOT FOUND in PATH"
        FAIL=$((FAIL + 1))
        return
    fi

    local version_output
    version_output=$("${tool}" ${version_flag} 2>&1 | head -3) || true

    if echo "${version_output}" | grep -qE "${expected_pattern}"; then
        echo "[PASS] ${tool}: $(echo "${version_output}" | grep -oE "${expected_pattern}" | head -1)"
    else
        echo "[WARN] ${tool}: found but could not verify expected version"
        echo "       Output: $(echo "${version_output}" | head -1)"
        PASS=$((PASS + 1))
    fi
    PASS=$((PASS + 1))
}

echo "================================================================"
echo " Checking environment: ${CONDA_ENV}"
echo "================================================================"

# Check micromamba environment
if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    if [[ "${CONDA_DEFAULT_ENV}" == "${CONDA_ENV}" ]]; then
        echo "[PASS] Micromamba env: ${CONDA_DEFAULT_ENV} (active)"
        PASS=$((PASS + 1))
    else
        echo "[WARN] Micromamba env: '${CONDA_DEFAULT_ENV}' is active, expected '${CONDA_ENV}'"
        echo "       Run: micromamba activate ${CONDA_ENV}"
        FAIL=$((FAIL + 1))
    fi
else
    echo "[WARN] No micromamba environment detected. Activate '${CONDA_ENV}' before running pipeline."
fi

echo ""
echo "--- Required tools ---"

# Python
check_tool "python" "--version" "Python 3\."

# MAFFT
check_tool "mafft" "--version" "v7\.[0-9]+"

# trimAl
check_tool "trimal" "--version" "[0-9]+\.[0-9]+"

# IQTree (must be 'iqtree', not 'iqtree2')
check_tool "iqtree" "--version" "IQ-TREE"

# PhyKIT
check_tool "phykit" "--version" "[0-9]+\.[0-9]+"

# SeqKit (for counting sequences and lengths)
check_tool "seqkit" "--version" "seqkit"

echo ""
echo "--- Required Python libraries ---"

python -c "import pandas; print('[PASS] pandas:', pandas.__version__)" 2>/dev/null \
    || { echo "[FAIL] pandas: not importable"; FAIL=$((FAIL + 1)); }

python -c "import numpy; print('[PASS] numpy:', numpy.__version__)" 2>/dev/null \
    || { echo "[FAIL] numpy: not importable"; FAIL=$((FAIL + 1)); }

python -c "import Bio; print('[PASS] Biopython:', Bio.__version__)" 2>/dev/null \
    || { echo "[FAIL] Biopython: not importable"; FAIL=$((FAIL + 1)); }

python -c "from ete3 import Tree; print('[PASS] ete3: importable')" 2>/dev/null \
    || { echo "[FAIL] ete3: not importable (needed for 06_filter_paralogs.py)"; FAIL=$((FAIL + 1)); }

echo ""
echo "--- Required input files ---"

for f in \
    "${HOG_FILE}" \
    "${OG_SEQ_DIR}" \
    "${ORTHO_SPECIES_TREE}"; do
    if [[ -e "${f}" ]]; then
        echo "[PASS] Exists: ${f}"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] Missing: ${f}"
        FAIL=$((FAIL + 1))
    fi
done

# Count sequences in OG seq dir
N_OG=$(ls "${OG_SEQ_DIR}"/*.fa 2>/dev/null | wc -l)
echo "       OG sequence files found: ${N_OG}"

# Count rows in HOG file
N_HOGS=$(tail -n +2 "${HOG_FILE}" | wc -l)
echo "       HOGs in N1.tsv: ${N_HOGS}"

echo ""
echo "--- SLURM availability ---"
if command -v sbatch &>/dev/null; then
    echo "[PASS] sbatch found"
    PASS=$((PASS + 1))
    sinfo --noheader -o "%P %a" 2>/dev/null | while read part state; do
        if [[ "${part}" == "${PARTITION}" || "${part}" == "${PARTITION}*" ]]; then
            echo "[PASS] Partition '${PARTITION}' found (state: ${state})"
        fi
    done
else
    echo "[WARN] sbatch not found — are you on a login node?"
fi

echo ""
echo "================================================================"
echo " Summary: ${PASS} checks passed, ${FAIL} checks failed"
echo "================================================================"

if [[ ${FAIL} -gt 0 ]]; then
    echo ""
    echo "ERROR: ${FAIL} check(s) failed. Fix before running the pipeline."
    exit 1
else
    echo ""
    echo "All critical checks passed. Ready to run the pipeline."
    echo "Next step: bash 01_extract_markers.py  (actually: python 01_extract_markers.py)"
fi
