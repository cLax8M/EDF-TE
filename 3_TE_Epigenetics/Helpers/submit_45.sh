#!/bin/bash
#SBATCH --job-name=45_mac_tss_expression
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/logs/45_mac_tss_expression_%j.log
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/logs/45_mac_tss_expression_%j.log
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0:30:00
set -euo pipefail
source /home/carlos/micromamba/etc/profile.d/mamba.sh
micromamba activate esptreeEDFTE
export PATH="/home/carlos/micromamba/envs/earlgrey/bin:$PATH"
cd /legserv/Temp/Carlos/EDF-TE/results/TE_dynamics
Rscript scripts/45_mac_tss_expression.R
