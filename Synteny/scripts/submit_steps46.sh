#!/bin/bash
#SBATCH --job-name=synteny46
#SBATCH --partition=normal.1000h
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --output=/data/carlos/EDF-TE/EDF-TE/Synteny/data/logs/slurm_steps46_%j.log

cd /data/carlos/EDF-TE/EDF-TE/Synteny

eval "$(micromamba shell hook --shell bash)"
micromamba activate syntenyEDFTE

bash scripts/run_pipeline.sh --all --threads 16 --steps 4-6
