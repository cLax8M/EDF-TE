#!/bin/bash
#SBATCH --job-name=44_rna_6mA_extremes
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/logs/44_rna_6mA_extremes_%j.log
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/logs/44_rna_6mA_extremes_%j.log
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
set -euo pipefail
source /home/carlos/micromamba/etc/profile.d/mamba.sh
micromamba activate esptreeEDFTE
export PATH="/home/carlos/micromamba/envs/earlgrey/bin:$PATH"
cd /legserv/Temp/Carlos/EDF-TE/results/TE_dynamics
Rscript scripts/44_rna_6mA_extremes.R
