#!/bin/bash
#SBATCH --job-name=nuc_profiles
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/data/carlos/EDF-TE/logs/nuc_profiles_%j.log

eval "$(micromamba shell hook --shell bash)"
micromamba activate syntenyEDFTE

cd /data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts
python nucleosome_profiles.py
