#!/bin/bash
#SBATCH --job-name=identify_gremlins
#SBATCH --partition=normal.168h
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/logs/identify_gremlins_%j.log
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/logs/identify_gremlins_%j.err

set -euo pipefail

module load BLAST+/2.16.0-gompi-2023b

ANNOT_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
QUERY="/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_centromeric.fasta"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/data/gremlin/blast_results"
mkdir -p "$OUT_DIR"

# For each species EarlGrey directory
for SPECIES_DIR in "$ANNOT_DIR"/MU_*_EarlGrey "$ANNOT_DIR"/BL_*_EarlGrey; do
    SPECIES=$(basename "$SPECIES_DIR" _EarlGrey)
    SUMDIR="$SPECIES_DIR/${SPECIES}_summaryFiles"
    LIB="$SUMDIR/${SPECIES}_combined_library.fasta"

    [[ -f "$LIB" ]] || { echo "WARNING: no library for $SPECIES"; continue; }

    OUT="$OUT_DIR/${SPECIES}_gremlin_blast.tsv"
    [[ -f "$OUT" ]] && { echo "SKIP $SPECIES (exists)"; continue; }

    echo "Processing $SPECIES..."

    # Make BLAST db (temporary)
    DB_TMP=$(mktemp -d)
    makeblastdb -in "$LIB" -dbtype nucl -out "$DB_TMP/db" -logfile /dev/null 2>/dev/null

    # BLAST: >=30% identity, >=50% query coverage (as Navarro-Mendoza 2019)
    blastn -query "$QUERY" \
           -db "$DB_TMP/db" \
           -outfmt "6 qseqid sseqid pident qcovs length evalue" \
           -perc_identity 30 \
           -qcov_hsp_perc 50 \
           -num_threads 8 \
           -evalue 1e-5 \
           -out "$OUT" 2>/dev/null || true

    rm -rf "$DB_TMP"

    NHIT=$(wc -l < "$OUT" 2>/dev/null || echo 0)
    echo "  $SPECIES: $NHIT hits"
done

echo "All done. Results in $OUT_DIR"
echo "Summarizing hit families per species..."

# Produce summary: species -> list of unique TE family names (from sseqid, strip after #)
SUMMARY="$OUT_DIR/gremlin_families_per_species.tsv"
echo -e "species\tfamily_name\tnum_queries_hit" > "$SUMMARY"

for TSV in "$OUT_DIR"/*_gremlin_blast.tsv; do
    SPECIES=$(basename "$TSV" _gremlin_blast.tsv)
    if [[ -s "$TSV" ]]; then
        # sseqid format: "rnd-N_family-N#CLASS" or "KNOWN_NAME#CLASS @species"
        # Extract the part before #, strip leading >
        awk -v sp="$SPECIES" '{
            fam = $2
            sub(/#.*/, "", fam)
            fams[fam] = fams[fam]; queries[fam][NR] = $1
        }
        END {
            for (f in fams) {
                n = 0
                for (q in queries[f]) n++
                print sp "\t" f "\t" n
            }
        }' "$TSV" >> "$SUMMARY"
    fi
done

echo "Summary written to $SUMMARY"
wc -l "$SUMMARY"
