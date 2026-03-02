#!/bin/bash
#SBATCH --job-name=AcancasteEG
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=16
#SBATCH --mem=45G
#SBATCH --output=/data/carlos/EDF-TE/logs/AcancasteEG_%j.log

module use /data/toby/tools/modulefiles
module load earlgrey

earlGreyLibConstruct \
  -g /legserv/Temp/Carlos/EDF-TE/data/genomes/HO_Acancaste.fasta \
  -s AcantamoebaCastellanii \
  -o /data/carlos/EDF-TE/TE_analysis \
  -t 16
