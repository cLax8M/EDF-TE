#!/bin/bash
# Extract 60kb windows and run pairwise BLASTN for v3 figure
# Groups: Rhizopus | Mucor | Phycomyces
set -euo pipefail

SAMTOOLS=/opt/easybuild/software/SAMtools/1.21-GCC-13.3.0/bin/samtools
BLASTN=/opt/easybuild/software/BLAST+/2.16.0-gompi-2024a/bin/blastn
MAKEDB=/opt/easybuild/software/BLAST+/2.16.0-gompi-2024a/bin/makeblastdb

GDIR=/legserv/Temp/Carlos/EDF-TE/data/genomes
EXTRA=/legserv/Temp/Carlos/EDF-TE/data/Extra
PHYBL_FA=/data/carlos/EDF-TE/TE_analysis/sex_locus/extra_tblastn/PhyblU21.fasta
OUTDIR=/data/carlos/EDF-TE/TE_analysis/sex_locus/synteny_figure
WINDIR=${OUTDIR}/windows60
mkdir -p "$WINDIR"

echo "=== Extracting 60kb windows ==="

# species -> fasta, scaffold, winStart-winEnd
declare -A FA=(
  ["Blatri"]="${GDIR}/MU_Blatri.fasta"
  ["Rhisto1"]="${GDIR}/MU_Rhisto1.fasta"
  ["Rhimi59"]="${GDIR}/MU_Rhimi59.fasta"
  ["Rhior3"]="${GDIR}/MU_Rhior3.fasta"
  ["Rhipu1"]="${GDIR}/MU_Rhipu1.fasta"
  ["Spoumb1"]="${GDIR}/MU_Spoumb1.fasta"

  ["Mucci3"]="${GDIR}/MU_Mucci3.fasta"
  ["Muccirc1"]="${GDIR}/MU_Muccirc1.fasta"
  ["ATCC52813"]="${EXTRA}/ATCC52813.fa"
  ["ATCC52814"]="${EXTRA}/ATCC52814.fa"
  ["PS15m"]="${EXTRA}/Mcir_PS15m.fna"
  ["PS15p"]="${EXTRA}/Mcir_PS15p.fna"
  ["Mucmuc1"]="${GDIR}/MU_Mucmuc1.fasta"
  ["Mucpir1"]="${GDIR}/MU_Mucpir1.fasta"
  ["Mycind1"]="${GDIR}/MU_Mycind1.fasta"
  ["Cokrec1"]="${GDIR}/MU_Cokrec1.fasta"

  ["PhyblU21"]="${PHYBL_FA}"
  ["Phybl4"]="${GDIR}/MU_Phybl4.fasta"
  ["PnitS607"]="${GDIR}/MU_PnitS607.fasta"
  ["PnitSJTU"]="${EXTRA}/PnitSJTU.fna"
)

declare -A REGION=(
  ["Blatri"]="scaffold_61:19836-79836"
  ["Rhisto1"]="scaffold_35:171724-231724"
  ["Rhimi59"]="scaffold_1:338056-398056"
  ["Rhior3"]="supercontig_3.1:3279818-3339818"
  ["Rhipu1"]="scaffold_19:114991-174991"
  ["Spoumb1"]="scaffold_13:20947-80947"

  ["Mucci3"]="scaffold_1:2183094-2243094"
  ["Muccirc1"]="KE_123975:43879-103879"
  ["ATCC52813"]="scaffold_2:3186275-3246275"
  ["ATCC52814"]="scaffold_1:3983605-4043605"
  ["PS15m"]="JBOZQX010000005.1:1669161-1729161"
  ["PS15p"]="JBNRQL010000006.1:1709974-1769974"
  ["Mucmuc1"]="scaffold_63:76533-136533"
  ["Mucpir1"]="scaffold_29:149236-209236"
  ["Mycind1"]="scaffold_10:602379-662379"
  ["Cokrec1"]="scaffold_3:440980-500980"

  ["PhyblU21"]="scaffold_17:302700-362700"
  ["Phybl4"]="scaffold_10:696371-756371"
  ["PnitS607"]="scaffold_54:12868-72868"
  ["PnitSJTU"]="JAHRBC010000002.1:6816920-6876920"
)

for SP in "${!REGION[@]}"; do
    OUT="${WINDIR}/${SP}_60kb.fasta"
    if [ -f "$OUT" ]; then
        echo "  $SP: already extracted"
        continue
    fi
    FASTA="${FA[$SP]}"
    REG="${REGION[$SP]}"
    if [ ! -f "$FASTA" ]; then
        echo "  WARNING: missing fasta for $SP: $FASTA"
        continue
    fi
    # Index if needed
    [ -f "${FASTA}.fai" ] || $SAMTOOLS faidx "$FASTA"
    $SAMTOOLS faidx "$FASTA" "$REG" > "$OUT" 2>/dev/null
    if [ ! -s "$OUT" ]; then
        echo "  WARNING: empty window for $SP ($REG)"
        rm -f "$OUT"
    else
        echo "  $SP: extracted $(wc -c < "$OUT") bytes"
    fi
done

echo ""
echo "=== Running pairwise BLASTN ==="

run_blastn() {
    local QSP="$1"
    local TSP="$2"
    local OUTF="${OUTDIR}/${QSP}_vs_${TSP}.paf"
    if [ -f "$OUTF" ]; then
        echo "  ${QSP} vs ${TSP}: already done"
        return
    fi
    local QFA="${WINDIR}/${QSP}_60kb.fasta"
    local TFA="${WINDIR}/${TSP}_60kb.fasta"
    [ -f "$QFA" ] || { echo "  SKIP ${QSP} vs ${TSP}: missing query fasta"; return; }
    [ -f "$TFA" ] || { echo "  SKIP ${QSP} vs ${TSP}: missing target fasta"; return; }

    # Make DB from target
    local DB="${WINDIR}/${TSP}_blastdb"
    $MAKEDB -in "$TFA" -dbtype nucl -out "$DB" -parse_seqids 2>/dev/null

    echo -n "  ${QSP} vs ${TSP}: running BLASTN..."
    $BLASTN -query "$QFA" -db "$DB" \
        -out "$OUTF" \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
        -evalue 1e-10 \
        -perc_identity 60 \
        -word_size 11 \
        -num_threads 8 2>/dev/null
    echo " $(wc -l < "$OUTF") hits"
}

# Rhizopus group
run_blastn Blatri   Rhisto1   # existing but re-run for new 60kb
run_blastn Rhisto1  Rhimi59
run_blastn Rhimi59  Rhior3
run_blastn Rhior3   Rhipu1
run_blastn Rhipu1   Spoumb1

# Mucor group
run_blastn Mucci3   Muccirc1
run_blastn Muccirc1 ATCC52813
run_blastn ATCC52813 ATCC52814
run_blastn ATCC52814 PS15m
run_blastn PS15m    PS15p
run_blastn PS15p    Mucmuc1
run_blastn Mucmuc1  Mucpir1
run_blastn Mucpir1  Mycind1
run_blastn Mycind1  Cokrec1

# Phycomyces group
run_blastn PhyblU21 Phybl4
run_blastn Phybl4   PnitS607
run_blastn PnitS607 PnitSJTU

echo ""
echo "=== DONE ==="
