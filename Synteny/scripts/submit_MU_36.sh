#!/bin/bash
#SBATCH --job-name=syntenyMU
#SBATCH --partition=normal.1000h
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --output=/data/carlos/EDF-TE/EDF-TE/Synteny/data/logs/slurm_MU36_%j.log

cd /data/carlos/EDF-TE/EDF-TE/Synteny

eval "$(micromamba shell hook --shell bash)"
micromamba activate syntenyEDFTE

bash scripts/run_pipeline.sh --clade MU --steps 3-6 --threads 16
