#!/bin/bash
#SBATCH --job-name=37b_6mA_expr
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/logs/37b_6mA_expr_%j.log
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/logs/37b_6mA_expr_%j.log
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00

set -euo pipefail
source /home/carlos/micromamba/etc/profile.d/mamba.sh
micromamba activate esptreeEDFTE

cd /legserv/Temp/Carlos/EDF-TE/results/TE_dynamics
Rscript scripts/37b_6mA_expression_combined.R
