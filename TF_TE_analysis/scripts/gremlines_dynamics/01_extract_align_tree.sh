#!/bin/bash
#SBATCH --job-name=grem_dyn_align
#SBATCH --partition=normal.168h
#SBATCH --exclude=LEGcompute1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/grem_dyn_align_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/grem_dyn_align_%j.err
set -euo pipefail

module load BEDTools/2.31.1-GCC-13.3.0
module load MAFFT/7.526-GCC-13.3.0-with-extensions
module load trimAl/1.5.0-GCCcore-13.3.0
module load IQ-TREE/3.0.1-gompi-2024a

# === Use the EXACT JGI Muccir1_3 (MU402) assembly used by Navarro-Mendoza 2019 ===
# GCA_010203745.1 = JGI Muccir1_3 = MU402, 24 scaffolds, 36.82 Mb, PacBio.
# Table S3 fits 56/56 entries (vs 52/56 in the newer ASM4025672v1=MU_Mucci3).
GENOME=/legserv/Temp/Carlos/EDF-TE/data/Extra/mu402GenomeJGI.fa
S3_BED=/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_table_S3.bed
EG_GFF=/data/carlos/EDF-TE/TE_analysis/annotated/MU_Mucci3_EarlGrey/MU_Mucci3_summaryFiles/MU_Mucci3.filteredRepeats.gff
OUT=/legserv/Temp/Carlos/EDF-TE/results/Gremlines/MU402
mkdir -p "$OUT"
cd "$OUT"

# Ensure .fai exists
[ -f "${GENOME}.fai" ] || samtools faidx "$GENOME"

# Build chromosome BED (sanity intersection)
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' "${GENOME}.fai" > chroms.bed

# === Ground truth from Navarro-Mendoza 2019 Table S3 ===
# Clusters 1-3 = 31 AUTONOMOUS Grem-LINE1 (full-length 5.5-6 kbp, ORF1+ORF2)
# Clusters 4-9 + GremLINE1_10 + NoCen = 25 REMNANTS (truncated)

# 1a. AUTONOMOUS = clusters 1-3 (used for the phylogeny)
awk '$4 ~ /^GremLINE1_[123]\./ {print}' "$S3_BED" \
    | bedtools intersect -a - -b chroms.bed -f 1.0 -u \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_AUT",$3-$2,$6}' \
    > gremlin_autonomous.bed

# 1b. REMNANTS = clusters 4-9 + 10 + NoCen (for K80 / zone analyses only)
awk '$4 ~ /^GremLINE1_([4-9]\.|10\.|NoCen)/ {print}' "$S3_BED" \
    | bedtools intersect -a - -b chroms.bed -f 1.0 -u \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_REM",$3-$2,$6}' \
    > gremlin_remnants.bed

cat gremlin_autonomous.bed gremlin_remnants.bed > gremlin_copies_all.bed

# 1c. Outgroup: 5 longest family-218 copies (>=1kb) from EarlGrey on MU_Mucci3
#     (close MU402 relative; gremlin family-218 is the diverged sister L1 ~30%).
#     EG annotation isn't done on the JGI assembly but the consensus library is the
#     same and family-218 is found in all Mucor strains. We extract them here from the
#     CBSv3 fasta (same species, slightly different scaffolding).
CBSV3=/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Mucci3.fasta

grep -iE 'family-218' "$EG_GFF" \
  | awk -F'\t' 'BEGIN{OFS="\t"} { match($9, /Name=([^;]+)/, a);
       fam=a[1]; sub(/.*family-/,"f",fam); sub(/#.*/,"",fam);
       print $1, $4-1, $5, fam, $6, $7 }' \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d 1000 -s -c 4,6 -o distinct,distinct \
  | awk -F'\t' 'BEGIN{OFS="\t"} ($3-$2)>=1000 {print $1,$2,$3,$4"_OUT_"NR,$3-$2,$5}' \
  | sort -k5,5 -nr | head -5 > gremlin_outgroup.bed

# Outgroup sequences (from CBSv3 — same species, family-218 conserved)
bedtools getfasta -fi "$CBSV3" -bed gremlin_outgroup.bed -nameOnly -s > gremlin_outgroup.fasta

NAUT=$(wc -l < gremlin_autonomous.bed)
NREM=$(wc -l < gremlin_remnants.bed)
NOG=$(wc  -l < gremlin_outgroup.bed)
echo "[INFO] Genome: MU402 JGI Muccir1_3 (GCA_010203745.1) — 24 scaffolds, 36.82 Mb"
echo "[INFO] AUTONOMOUS Grem-LINE1 (Navarro Table S3 clusters 1-3): $NAUT (target: 31)"
echo "[INFO] REMNANT Grem-LINE1 (clusters 4-9 + 10 + NoCen):        $NREM (target: 25)"
echo "[INFO] OUTGROUP family-218 (diverged L1, from CBSv3):         $NOG"

# 2. Extract sequences from JGI MU402
bedtools getfasta -fi "$GENOME" -bed gremlin_copies_all.bed   -nameOnly -s > gremlin_copies_all.fasta
bedtools getfasta -fi "$GENOME" -bed gremlin_autonomous.bed   -nameOnly -s > gremlin_autonomous.fasta
cat gremlin_autonomous.fasta gremlin_outgroup.fasta > gremlin_for_tree.fasta

# 3. MAFFT alignments
mafft --auto --thread 8 gremlin_copies_all.fasta > gremlin_aln_all.fasta
mafft --auto --thread 8 gremlin_for_tree.fasta   > gremlin_aln.fasta

# 4. trimAl
trimal -in gremlin_aln_all.fasta -out gremlin_aln_all.trim.fasta -gt 0.4
trimal -in gremlin_aln.fasta     -out gremlin_aln.trim.fasta     -gt 0.5

# 5. IQTree3 ModelFinder + UFBoot1000 (root via family-218 outgroup if possible)
OG=$(awk '/_OUT_/ {print $4}' gremlin_outgroup.bed | head -1)
if [ -n "$OG" ]; then
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 \
            -pre gremlin_iqtree --quiet -redo -o "$OG" || \
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 -pre gremlin_iqtree --quiet -redo
else
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 -pre gremlin_iqtree --quiet -redo
fi

echo "[DONE] alignment + tree at $OUT"
#!/bin/bash
#SBATCH --job-name=grem_dyn_align
#SBATCH --partition=normal.168h
#SBATCH --exclude=LEGcompute1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/grem_dyn_align_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/grem_dyn_align_%j.err
set -euo pipefail

module load BEDTools/2.31.1-GCC-13.3.0
module load MAFFT/7.526-GCC-13.3.0-with-extensions
module load trimAl/1.5.0-GCCcore-13.3.0
module load IQ-TREE/3.0.1-gompi-2024a

GENOME=/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Mucci3.fasta
S3_BED=/legserv/Temp/Carlos/EDF-TE/data/gremlin/GremLINE1_table_S3.bed
OUT=/legserv/Temp/Carlos/EDF-TE/results/Gremlines/Mucci3
mkdir -p "$OUT"
cd "$OUT"

# === Ground truth from Navarro-Mendoza 2019 Table S3 ===
# Clusters 1-3 = 31 AUTONOMOUS Grem-LINE1 (full-length 5.5-6 kbp, ORF1+ORF2)
# Clusters 4-9 + NoCen = 25 REMNANTS (truncated)
# Of the 31 autonomous, 27 fall within MU_Mucci3 (CBSv3) scaffolds (4 in scaffold_15/18 missing).

# Build chromosome BED for intersection (filter out S3 entries on missing scaffolds)
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' "${GENOME}.fai" > chroms.bed

# 1a. AUTONOMOUS = clusters 1-3 (use these for tree)
awk '$4 ~ /^GremLINE1_[123]\./ {print}' "$S3_BED" \
    | bedtools intersect -a - -b chroms.bed -f 1.0 -u \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_AUT",$3-$2,$6}' \
    > gremlin_autonomous.bed

# 1b. REMNANTS = clusters 4-9 + NoCen (for K80/zone)
awk '$4 ~ /^GremLINE1_[4-9]\.|^GremLINE1_NoCen|^GremLINE1_10/ {print}' "$S3_BED" \
    | bedtools intersect -a - -b chroms.bed -f 1.0 -u \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_REM",$3-$2,$6}' \
    > gremlin_remnants.bed

cat gremlin_autonomous.bed gremlin_remnants.bed > gremlin_copies_all.bed

# 1c. Outgroup: 5 longest family-218 copies (>=1kb) from EarlGrey GFF (diverged L1)
EG_GFF=/data/carlos/EDF-TE/TE_analysis/annotated/MU_Mucci3_EarlGrey/MU_Mucci3_summaryFiles/MU_Mucci3.filteredRepeats.gff
grep -iE 'family-218' "$EG_GFF" \
  | awk -F'\t' 'BEGIN{OFS="\t"} { match($9, /Name=([^;]+)/, a);
       fam=a[1]; sub(/.*family-/,"f",fam); sub(/#.*/,"",fam);
       print $1, $4-1, $5, fam, $6, $7 }' \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d 1000 -s -c 4,6 -o distinct,distinct \
  | awk -F'\t' 'BEGIN{OFS="\t"} ($3-$2)>=1000 {print $1,$2,$3,$4"_OUT_"NR,$3-$2,$5}' \
  | sort -k5,5 -nr | head -5 > gremlin_outgroup.bed

cat gremlin_autonomous.bed gremlin_outgroup.bed > gremlin_for_tree.bed
NAUT=$(wc -l < gremlin_autonomous.bed)
NREM=$(wc -l < gremlin_remnants.bed)
NOG=$(wc -l < gremlin_outgroup.bed)
echo "[INFO] AUTONOMOUS Grem-LINE1 (Navarro Table S3 clusters 1-3 in MU_Mucci3): $NAUT"
echo "[INFO] REMNANT Grem-LINE1 (clusters 4-9 + NoCen + 10):                     $NREM"
echo "[INFO] OUTGROUP family-218 (diverged L1):                                  $NOG"

# 2. Extract sequences
bedtools getfasta -fi "$GENOME" -bed gremlin_copies_all.bed -nameOnly -s > gremlin_copies_all.fasta
bedtools getfasta -fi "$GENOME" -bed gremlin_for_tree.bed   -nameOnly -s > gremlin_for_tree.fasta

# 3. MAFFT alignments
mafft --auto --thread 8 gremlin_copies_all.fasta > gremlin_aln_all.fasta
mafft --auto --thread 8 gremlin_for_tree.fasta   > gremlin_aln.fasta

# 4. trimAl
trimal -in gremlin_aln_all.fasta -out gremlin_aln_all.trim.fasta -gt 0.4
trimal -in gremlin_aln.fasta     -out gremlin_aln.trim.fasta     -gt 0.5

# 5. IQTree3 ModelFinder + UFBoot1000 (root via family-218 outgroup if possible)
OG=$(awk '/_OUT_/ {print $4"("$6")"}' gremlin_for_tree.bed | head -1)
if [ -n "$OG" ]; then
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 \
            -pre gremlin_iqtree --quiet -redo -o "$OG" || \
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 -pre gremlin_iqtree --quiet -redo
else
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 -pre gremlin_iqtree --quiet -redo
fi

echo "[DONE] alignment + tree at $OUT"
#!/bin/bash
#SBATCH --job-name=grem_dyn_align
#SBATCH --partition=normal.168h
#SBATCH --exclude=LEGcompute1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/grem_dyn_align_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/grem_dyn_align_%j.err
set -euo pipefail

module load BEDTools/2.31.1-GCC-13.3.0
module load MAFFT/7.526-GCC-13.3.0-with-extensions
module load trimAl/1.5.0-GCCcore-13.3.0
module load IQ-TREE/3.0.1-gompi-2024a

SP=MU_Mucci3
GENOME=/legserv/Temp/Carlos/EDF-TE/data/genomes/${SP}.fasta
EG_GFF=/data/carlos/EDF-TE/TE_analysis/annotated/${SP}_EarlGrey/${SP}_summaryFiles/${SP}.filteredRepeats.gff
OUT=/legserv/Temp/Carlos/EDF-TE/results/Gremlines/Mucci3
mkdir -p "$OUT"
cd "$OUT"

# REAL Grem-LINE1 families validated by tBLASTN of GremLINE1_orf2.faa (Navarro-Mendoza 2019)
# vs MU_Mucci3 EarlGrey library: rnd-1_family-198 / 206 / 180 / 106 hit at >70% identity.
# family-218 is more diverged (~30%) — used as outgroup to root the tree.
#
# Per Navarro-Mendoza 2019: AUTONOMOUS elements = Clusters 1-3 = full-length 5.4-6 kbp
# with both ORF1 (Zn fingers) + ORF2 (APE+RT). Remnants = truncated.
# EarlGrey fragments adjacent hits → we MERGE within 1 kb on same strand to recover
# full-length elements (target ~31 autonomous per Navarro2019).

GREM_FAMS='family-(198|206|180|106)'

# 1a. Raw gremlin hits as 6-col BED, then sort
grep -iE "$GREM_FAMS" "$EG_GFF" \
  | awk -F'\t' 'BEGIN{OFS="\t"} { match($9, /Name=([^;]+)/, a);
       fam=a[1]; sub(/.*family-/,"f",fam); sub(/#.*/,"",fam);
       print $1, $4-1, $5, fam, $6, $7 }' \
  | sort -k1,1 -k2,2n > gremlin_hits_raw.bed

# 1b. Merge adjacent fragments same strand within 1 kb → recover full-length
bedtools merge -i gremlin_hits_raw.bed -d 1000 -s -c 4,6 -o distinct,distinct \
  | awk -F'\t' 'BEGIN{OFS="\t"} {
       len=$3-$2;
       cls=(len>=5000)?"AUT":"REM";
       print $1, $2, $3, $4"_"cls"_"NR, len, $5
     }' > gremlin_copies_all.bed

awk '($3-$2)>=300' gremlin_copies_all.bed > gremlin_copies_filtered.bed
mv gremlin_copies_filtered.bed gremlin_copies_all.bed

awk '/_AUT_/' gremlin_copies_all.bed > gremlin_autonomous.bed

# 1c. Outgroup: 5 longest family-218 copies (>=1kb, also merged)
grep -iE 'family-218' "$EG_GFF" \
  | awk -F'\t' 'BEGIN{OFS="\t"} { match($9, /Name=([^;]+)/, a);
       fam=a[1]; sub(/.*family-/,"f",fam); sub(/#.*/,"",fam);
       print $1, $4-1, $5, fam, $6, $7 }' \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - -d 1000 -s -c 4,6 -o distinct,distinct \
  | awk -F'\t' 'BEGIN{OFS="\t"} ($3-$2)>=1000 {print $1,$2,$3,$4"_OUT_"NR,$3-$2,$5}' \
  | sort -k5,5 -nr | head -5 > gremlin_outgroup.bed

cat gremlin_autonomous.bed gremlin_outgroup.bed > gremlin_for_tree.bed
NALL=$(wc -l < gremlin_copies_all.bed)
NAUT=$(wc -l < gremlin_autonomous.bed)
NOG=$(wc -l < gremlin_outgroup.bed)
echo "[INFO] Total real Grem-LINE1 copies (>=300bp, MERGED): $NALL"
echo "[INFO] AUTONOMOUS copies (>=5kb merged): $NAUT (target ~31 per Navarro2019)"
echo "[INFO] Outgroup family-218 copies: $NOG"

# 2. Extract sequences
bedtools getfasta -fi "$GENOME" -bed gremlin_copies_all.bed -nameOnly -s > gremlin_copies_all.fasta
bedtools getfasta -fi "$GENOME" -bed gremlin_for_tree.bed -nameOnly -s > gremlin_for_tree.fasta

# 3. MAFFT alignments
mafft --auto --thread 8 gremlin_copies_all.fasta > gremlin_aln_all.fasta
mafft --auto --thread 8 gremlin_for_tree.fasta > gremlin_aln.fasta

# 4. trimAl
trimal -in gremlin_aln_all.fasta -out gremlin_aln_all.trim.fasta -gt 0.4
trimal -in gremlin_aln.fasta     -out gremlin_aln.trim.fasta     -gt 0.5

# 5. IQTree3 on autonomous + outgroup only
OG=$(awk '/_OUT_/ {print $4"("$6")"}' gremlin_for_tree.bed | head -1)
if [ -n "$OG" ]; then
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 \
            -pre gremlin_iqtree --quiet -redo -o "$OG" || \
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 -pre gremlin_iqtree --quiet -redo
else
    iqtree3 -s gremlin_aln.trim.fasta -m MFP -B 1000 -T 8 -pre gremlin_iqtree --quiet -redo
fi

echo "[DONE] alignment + tree at $OUT"
#!/bin/bash
#SBATCH --job-name=grem_dyn_align
#SBATCH --partition=normal.168h
#SBATCH --exclude=LEGcompute1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/grem_dyn_align_%j.log
#SBATCH --error=/data/carlos/EDF-TE/logs/grem_dyn_align_%j.err
set -euo pipefail

module load BEDTools/2.31.1-GCC-13.3.0
module load MAFFT/7.526-GCC-13.3.0-with-extensions
module load trimAl/1.5.0-GCCcore-13.3.0
module load IQ-TREE/3.0.1-gompi-2024a

SP=MU_Mucci3
GENOME=/legserv/Temp/Carlos/EDF-TE/data/genomes/${SP}.fasta
EG_BED=/data/carlos/EDF-TE/TE_analysis/annotated/${SP}_EarlGrey/${SP}_summaryFiles/${SP}.filteredRepeats.bed
OUT=/legserv/Temp/Carlos/EDF-TE/results/Gremlines/Mucci3
mkdir -p "$OUT"
cd "$OUT"

# 1. Extract LINE/L1-Tx1 copies (Grem-LINE1 family) ≥500bp
awk -F'\t' 'BEGIN{OFS="\t"} $4 ~ /L1-Tx1/ && ($3-$2)>=500 {print $1,$2,$3,$4"|"NR,$5,$6}' "$EG_BED" > L1Tx1_copies.bed
N=$(wc -l < L1Tx1_copies.bed)
echo "[INFO] L1-Tx1 copies (>=500bp): $N"

# 2. Extract sequences
bedtools getfasta -fi "$GENOME" -bed L1Tx1_copies.bed -nameOnly -s > L1Tx1_copies.fasta
grep -c '>' L1Tx1_copies.fasta

# 3. MAFFT alignment
mafft --auto --thread 8 L1Tx1_copies.fasta > L1Tx1_aln.fasta

# 4. trimAl
trimal -in L1Tx1_aln.fasta -out L1Tx1_aln.trim.fasta -gt 0.5

# 5. IQTree (nucleotide ModelFinder + UFBoot1000)
iqtree3 -s L1Tx1_aln.trim.fasta -m MFP -B 1000 -T 8 -pre L1Tx1_iqtree --quiet -redo

echo "[DONE] alignment + tree at $OUT"
