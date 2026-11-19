#!/bin/bash
#SBATCH --job-name=39_6mA_scatter
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/logs/39_6mA_scatter_%j.log
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/logs/39_6mA_scatter_%j.log
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:30:00

set -euo pipefail
source /home/carlos/micromamba/etc/profile.d/mamba.sh
micromamba activate esptreeEDFTE

cd /legserv/Temp/Carlos/EDF-TE/results/TE_dynamics
Rscript scripts/39_6mA_te_scatter.R
