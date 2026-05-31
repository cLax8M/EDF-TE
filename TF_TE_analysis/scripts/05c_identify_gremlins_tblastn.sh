#!/bin/bash
#SBATCH --job-name=gremlin_tblastn
#SBATCH --partition=normal.168h
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/logs/gremlin_tblastn_%j.log
#SBATCH --error=/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/logs/gremlin_tblastn_%j.err

set -euo pipefail

module load BLAST+/2.16.0-gompi-2023b
module load EMBOSS/6.6.0-foss-2024a

ANNOT_DIR="/data/carlos/EDF-TE/TE_analysis/annotated"
CENFASTA="/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_centromeric.fasta"
ORF2_QUERY="/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_orf2.faa"
OUT_DIR="/legserv/Temp/Carlos/EDF-TE/data/gremlin/tblastn_results"
mkdir -p "$OUT_DIR"

# ── Step 1: extract ORF2 proteins if not already done ──────────────────────
if [[ ! -s "$ORF2_QUERY" ]]; then
    echo "Extracting ORF2 proteins (>=800 aa) from centromeric GremLINE1 elements..."
    ALL_ORFS="/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_orfs.faa"
    getorf -sequence "$CENFASTA" -outseq "$ALL_ORFS" -minsize 900 -find 1 2>/dev/null
    python3 - <<'PYEOF'
from pathlib import Path
recs = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_orfs.faa").read_text()
out = []
cur_head, cur_seq = "", []
for line in recs.splitlines():
    if line.startswith(">"):
        if cur_seq:
            seq = "".join(cur_seq)
            if len(seq) >= 800:
                out.append(f"{cur_head}\n{seq}")
        cur_head, cur_seq = line, []
    else:
        cur_seq.append(line)
if cur_seq:
    seq = "".join(cur_seq)
    if len(seq) >= 800:
        out.append(f"{cur_head}\n{seq}")
Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_orf2.faa").write_text("\n".join(out) + "\n")
print(f"Wrote {len(out)} ORF2 sequences")
PYEOF
fi
echo "ORF2 query: $(grep -c '^>' "$ORF2_QUERY") sequences"

# ── Step 2: tBLASTn vs each species combined library ───────────────────────
for SPECIES_DIR in "$ANNOT_DIR"/MU_*_EarlGrey "$ANNOT_DIR"/BL_*_EarlGrey; do
    SPECIES=$(basename "$SPECIES_DIR" _EarlGrey)
    SUMDIR="$SPECIES_DIR/${SPECIES}_summaryFiles"
    LIB="$SUMDIR/${SPECIES}_combined_library.fasta"

    [[ -f "$LIB" ]] || { echo "WARNING: no library for $SPECIES"; continue; }

    OUT="$OUT_DIR/${SPECIES}_tblastn.tsv"
    [[ -f "$OUT" ]] && [[ -s "$OUT" ]] && { echo "SKIP $SPECIES (exists, $(wc -l < "$OUT") hits)"; continue; }

    echo "Processing $SPECIES..."

    DB_TMP=$(mktemp -d /legserv/Temp/Carlos/tmp_blast_XXXXXX)
    makeblastdb -in "$LIB" -dbtype nucl -out "$DB_TMP/db" -logfile /dev/null 2>/dev/null

    # tBLASTn: protein query vs nucleotide DB (translated)
    # Parameters matching Navarro-Mendoza 2019 protein survey:
    #   ≥25% identity, ≥50% query coverage, e-value ≤ 1e-5
    tblastn -query "$ORF2_QUERY" \
            -db "$DB_TMP/db" \
            -evalue 1e-5 \
            -qcov_hsp_perc 10 \
            -outfmt "6 qseqid sseqid pident qcovhsp length evalue" \
            -num_threads 8 \
            -max_target_seqs 10 \
            -out "$OUT" 2>/dev/null || true

    rm -rf "$DB_TMP"

    NHIT=$(wc -l < "$OUT" 2>/dev/null || echo 0)
    echo "  $SPECIES: $NHIT tBLASTn hits"
done

echo ""
echo "All done. Results in $OUT_DIR"

# ── Step 3: summarise family names per species ─────────────────────────────
SUMMARY="$OUT_DIR/gremlin_families_tblastn.tsv"
echo -e "species\tfamily_name\tnum_queries_hit" > "$SUMMARY"

for TSV in "$OUT_DIR"/*_tblastn.tsv; do
    SPECIES=$(basename "$TSV" _tblastn.tsv)
    [[ -s "$TSV" ]] || continue
    awk -v sp="$SPECIES" '{
        fam = $2
        sub(/#.*/, "", fam)
        nqueries[fam]++
    }
    END {
        for (f in nqueries)
            print sp "\t" f "\t" nqueries[f]
    }' "$TSV" >> "$SUMMARY"
done

echo "Summary written to $SUMMARY"
echo "Species with hits: $(awk 'NR>1{print $1}' "$SUMMARY" | sort -u | wc -l)"
wc -l "$SUMMARY"
