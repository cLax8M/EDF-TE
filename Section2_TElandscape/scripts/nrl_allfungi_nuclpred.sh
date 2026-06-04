#!/bin/bash
#SBATCH --job-name=nrl_fungi
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/nrl_allfungi_%j.log

eval "$(micromamba shell hook --shell bash)"
micromamba activate syntenyEDFTE

cd /data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts
python nrl_allfungi_nuclpred.py
