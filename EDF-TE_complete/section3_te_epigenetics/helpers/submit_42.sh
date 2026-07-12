#!/bin/bash
#SBATCH --job-name=42_mac_tss
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/logs/42_mac_tss_%j.log
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/logs/42_mac_tss_%j.log
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00

set -euo pipefail
source /home/carlos/micromamba/etc/profile.d/mamba.sh
micromamba activate esptreeEDFTE

cd /legserv/Temp/Carlos/EDF-TE/results/TE_dynamics
Rscript scripts/42_6mA_mac_at_tss.R
