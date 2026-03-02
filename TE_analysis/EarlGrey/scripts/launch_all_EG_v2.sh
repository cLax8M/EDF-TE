#!/bin/bash
# Launch EarlGrey v2 (Dfam 3.9) for all non-fungal species
# Excludes DictyosteliumDiscoideum and MonosigaBrevicollis (already running as validation)
# Usage: bash launch_all_EG_v2.sh [--dry-run]

GENOME_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
DRY_RUN=0
[[ "$1" == "--dry-run" ]] && DRY_RUN=1

# Format: "SPECIES_NAME|GENOME_FILE|SHORT_NAME"
declare -a SPECIES=(
  "AcantamoebaCastellanii|HO_Acancaste.fasta|AcancasteEG"
  "Calliphoravicina|HO_Callivi.fasta|CalliviEG"
  "CapsasporaOwczarzaki|HO_Capowcza.fasta|CapowczaEG"
  "CionaIntestinalis|HO_Ciointe.fasta|CiointeEG"
  "DrosophilaMelanogaster|HO_Drome.fasta|DromeEG"
  "FonticulaAlba|NU_Fontalba.fasta|FontalbaEG"
  "HydraVulgaris|HO_Hydravulg.fna|HydravulgEG"
  "Luciliacuprina|HO_Lucicu.fasta|LucicuEG"
  "MnemiopsisLeidy|HO_Mneleidy.fasta|MneleidyEG"
  "NaegleriaGruberi|HO_Naegrub.fasta|NaegrubEG"
  "ParvulariaAtlantis|NU_Parvatl.fasta|ParvatlEG"
  "SalpingoecaRosetta|HO_Salpirose.fasta|SalpiroseEG"
  "SphaeroformaArtica|HO_Sphaearti.fasta|SphaerartiEG"
  "ThecamonasTrahens|HO_Thecatra.fasta|ThecatraEG"
  "TrichoplaxAdhaerens|HO_Trichadha.fasta|TrichadhaEG"
)

echo "Submitting EarlGrey v2 jobs..."
echo "================================"

TMPDIR_SCRIPTS="/data/carlos/EDF-TE/logs/sbatch_scripts"
mkdir -p "$TMPDIR_SCRIPTS"

for entry in "${SPECIES[@]}"; do
  SPECIES_NAME=$(echo "$entry" | cut -d'|' -f1)
  GENOME_FILE=$(echo "$entry" | cut -d'|' -f2)
  SHORT_NAME=$(echo "$entry" | cut -d'|' -f3)
  GENOME_PATH="${GENOME_DIR}/${GENOME_FILE}"

  if [[ ! -f "$GENOME_PATH" ]]; then
    echo "WARNING: genome not found for $SPECIES_NAME ($GENOME_PATH) — skipping"
    continue
  fi

  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[DRY-RUN] Would submit: $SPECIES_NAME  genome=$GENOME_FILE  job=$SHORT_NAME"
    continue
  fi

  SCRIPT="${TMPDIR_SCRIPTS}/${SHORT_NAME}.sh"
  cat > "$SCRIPT" <<SBATCH
#!/bin/bash
#SBATCH --job-name=${SHORT_NAME}
#SBATCH --partition=normal.1000h
#SBATCH --cpus-per-task=16
#SBATCH --mem=45G
#SBATCH --output=/data/carlos/EDF-TE/logs/${SHORT_NAME}_%j.log

module use /data/toby/tools/modulefiles
module load earlgrey

earlGreyLibConstruct \\
  -g ${GENOME_PATH} \\
  -s ${SPECIES_NAME} \\
  -o /data/carlos/EDF-TE/TE_analysis \\
  -t 16
SBATCH

  JOB_ID=$(sbatch "$SCRIPT" | awk '{print $4}')
  echo "Submitted $SPECIES_NAME → job $JOB_ID  (log: logs/${SHORT_NAME}_${JOB_ID}.log)"
done

echo "================================"
echo "Done. Check queue with: squeue -u \$USER | grep EG"
