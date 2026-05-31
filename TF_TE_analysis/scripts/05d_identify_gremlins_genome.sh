#!/bin/bash
#SBATCH --job-name=gremlin_genome
#SBATCH --partition=normal.168h
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/logs/gremlin_genome_%A_%a.log
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/logs/gremlin_genome_%A_%a.err
#SBATCH --array=0-26

# tBLASTn of GremLINE1 ORF2 proteins vs full genome assemblies for the 27 MU
# species where tBLASTn vs TE library yielded no rnd-* de novo family hits.
# After the BLAST, genomic hit coords are intersected with EarlGrey GFF to
# recover the de novo family name (Name=RND-1_FAMILY-N) of each hit element.
#
# Array indices 0-26 map to the SPECIES_LIST below.

set -euo pipefail

module load BLAST+/2.16.0-gompi-2023b
module load BEDTools/2.31.1-GCC-13.3.0

# ── Species list (27 MU species with 0 rnd-* hits from library tBLASTn) ────
SPECIES_LIST=(
    MU_Apoossi  MU_Bacci1   MU_Baclam1  MU_Bifad1   MU_Calth1
    MU_Chlapade MU_Cunech1  MU_Dicele1  MU_Endsp1   MU_Fenlin
    MU_Halrad1  MU_Helpul1  MU_Jimfl    MU_Lenzyc1  MU_Liccor1
    MU_Licra1   MU_Mucend1  MU_Muclan1  MU_Pilano1  MU_Radspe1
    MU_Rhimi59  MU_Rhipu1   MU_Spoumb1  MU_Synrac1  MU_Syzme1
    MU_Thaele1  MU_Umbra1
)

SPECIES="${SPECIES_LIST[$SLURM_ARRAY_TASK_ID]}"
echo "=== Processing $SPECIES (array task $SLURM_ARRAY_TASK_ID) ==="

# ── Paths ────────────────────────────────────────────────────────────────────
GENOME_DIR="/legserv/Temp/Carlos/EDF-TE/data/genomes"
ANNOT_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
ORF2_QUERY="/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_orf2.faa"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/data/gremlin/genome_tblastn_results"
mkdir -p "$OUT_DIR"

GENOME="${GENOME_DIR}/${SPECIES}.fasta"
GFF=$(ls "${ANNOT_DIR}/${SPECIES}_EarlGrey/${SPECIES}_summaryFiles/${SPECIES}"*.filteredRepeats.gff \
      2>/dev/null | head -1 || true)

if [[ ! -f "$GENOME" ]]; then
    echo "ERROR: genome not found for $SPECIES: $GENOME"
    exit 1
fi

# ── Step 1: tBLASTn ORF2 proteins vs genome ──────────────────────────────────
BLAST_OUT="${OUT_DIR}/${SPECIES}_genome_tblastn.tsv"

if [[ -s "$BLAST_OUT" ]]; then
    echo "SKIP tBLASTn (exists, $(wc -l < "$BLAST_OUT") hits)"
else
    DB_TMP=$(mktemp -d /legserv/Temp/Carlos/tmp_gblast_XXXXXX)
    trap "rm -rf '$DB_TMP'" EXIT

    echo "Building BLAST db from genome $(du -sh "$GENOME" | cut -f1)..."
    makeblastdb -in "$GENOME" -dbtype nucl -out "$DB_TMP/db" \
                -logfile /dev/null 2>/dev/null

    echo "Running tBLASTn (33 ORF2 proteins × 8 threads)..."
    tblastn -query  "$ORF2_QUERY" \
            -db     "$DB_TMP/db" \
            -evalue 1e-5 \
            -qcov_hsp_perc 10 \
            -outfmt "6 qseqid sseqid pident qcovhsp sstart send evalue sstrand" \
            -num_threads 8 \
            -max_target_seqs 20 \
            -out "$BLAST_OUT" 2>/dev/null || true

    rm -rf "$DB_TMP"
    trap - EXIT
fi
NHIT=$(wc -l < "$BLAST_OUT")
echo "$SPECIES: $NHIT raw tBLASTn hits (evalue<=1e-5, qcov>=10%)"

# ── Step 2: convert BLAST hits to BED (0-based) ──────────────────────────────
# BLAST sstart > send for minus-strand hits; normalize to start < end
BED_OUT="${OUT_DIR}/${SPECIES}_genome_hits.bed"
awk 'BEGIN{OFS="\t"} {
    chr    = $2
    strand = ($5 < $6) ? "+" : "-"
    s      = ($5 < $6) ? $5-1 : $6-1
    e      = ($5 < $6) ? $6   : $5
    print chr, s, e, $1, $3, strand
}' "$BLAST_OUT" | sort -k1,1 -k2,2n > "$BED_OUT"
echo "  BED: $(wc -l < "$BED_OUT") entries"

# ── Step 3: intersect BLAST hits with EarlGrey GFF ───────────────────────────
# Only retain intersections with de novo (rnd-*) OR Unknown de novo entries
INTERSECT_OUT="${OUT_DIR}/${SPECIES}_families.tsv"

if [[ -z "$GFF" || ! -f "$GFF" ]]; then
    echo "WARNING: no EarlGrey GFF found for $SPECIES; skipping intersect"
    touch "$INTERSECT_OUT"
elif [[ ! -s "$BED_OUT" ]]; then
    echo "  No BLAST hits; skipping intersect"
    touch "$INTERSECT_OUT"
else
    # GFF has 9 tab-separated columns; bedtools appends them after the 6 BED cols
    # Result col layout: chr(1) start(2) end(3) qid(4) pident(5) str(6) |
    #   gff_chr(7) src(8) feat(9) gstart(10) gend(11) sc(12) gstr(13) fr(14) attrs(15)
    echo "  Intersecting with EarlGrey GFF..."
    bedtools intersect -a "$BED_OUT" -b "$GFF" -wb 2>/dev/null | \
        awk 'BEGIN{OFS="\t"} {
            attr = $15
            # Extract Name= value (uppercase in GFF, e.g. RND-1_FAMILY-83)
            if (match(attr, /Name=([^;]+)/, arr)) {
                fam = arr[1]
                # Normalize to lowercase rnd- convention used by FASTA headers
                gsub(/^RND-/, "rnd-", fam)
                gsub(/^FAMILY-/, "family-", fam)
                # Keep only de novo rnd-* families (exclude known Repbase entries)
                if (fam ~ /^rnd-/) print fam
            }
        }' | sort | uniq -c | sort -rn > "$INTERSECT_OUT"

    N_FAM=$(wc -l < "$INTERSECT_OUT")
    echo "  GFF intersect: $N_FAM de novo families with hits"
    [[ $N_FAM -gt 0 ]] && head -5 "$INTERSECT_OUT" | awk '{printf "    %s (n=%s)\n", $2, $1}'
fi

echo "=== Done: $SPECIES ==="
