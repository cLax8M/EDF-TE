#!/usr/bin/env bash
# setup_htt_envs.sh — create micromamba environments required by HTT pipeline steps 02 and 03
# Run once before submitting SLURM jobs:  bash setup_htt_envs.sh
set -euo pipefail

MICROMAMBA_BIN=$(which micromamba)
echo "Using micromamba: $MICROMAMBA_BIN"
echo "Creating HTT tool environments…"
echo ""

# ── HTT_python: standard Python (steps 02 step1) ──────────────────────────────
if ! micromamba env list | grep -q "^  HTT_python "; then
  echo "[1/5] Creating HTT_python …"
  micromamba create -n HTT_python -c conda-forge python=3.10 -y
else
  echo "[1/5] HTT_python already exists"
fi

# ── HTT_phylorank: phylorank tool (step 02 step2) ─────────────────────────────
if ! micromamba env list | grep -q "^  HTT_phylorank "; then
  echo "[2/5] Creating HTT_phylorank …"
  micromamba create -n HTT_phylorank -c conda-forge python=3.10 -y
  micromamba run -n HTT_phylorank pip install phylorank --quiet
else
  echo "[2/5] HTT_phylorank already exists"
fi

# ── HTT_python_ete3: Python + ete3 (step 02 step3; step 03 steps 1-2) ─────────
if ! micromamba env list | grep -q "^  HTT_python_ete3 "; then
  echo "[3/5] Creating HTT_python_ete3 (ete3 + numpy) …"
  micromamba create -n HTT_python_ete3 -c conda-forge -c etetoolkit ete3 python=3.10 -y
else
  echo "[3/5] HTT_python_ete3 already exists"
fi

# ── HTT_mafft: mafft aligner (step 03 step3) ──────────────────────────────────
# esptreeEDFTE already has mafft — just alias it if HTT_mafft doesn't exist
if ! micromamba env list | grep -q "^  HTT_mafft "; then
  echo "[4/5] Creating HTT_mafft …"
  micromamba create -n HTT_mafft -c bioconda -c conda-forge mafft python=3.10 -y
else
  echo "[4/5] HTT_mafft already exists"
fi

# ── HTT_pal2nal: pal2nal + R + seqinr + jsonlite (step 03 steps 4-5) ─────────
if ! micromamba env list | grep -q "^  HTT_pal2nal "; then
  echo "[5/5] Creating HTT_pal2nal …"
  micromamba create -n HTT_pal2nal -c bioconda -c conda-forge pal2nal python=3.10 -y
else
  echo "[5/5] HTT_pal2nal already exists"
fi

# ── HTT_r: R + seqinr + jsonlite (step 03 step5) ─────────────────────────────
if ! micromamba env list | grep -q "^  HTT_r "; then
  echo "[6/6] Creating HTT_r …"
  micromamba create -n HTT_r -c conda-forge -c r r-base r-seqinr r-jsonlite -y
else
  echo "[6/6] HTT_r already exists"
fi

echo ""
echo "All HTT environments ready:"
micromamba env list | grep -E "HTT_"
