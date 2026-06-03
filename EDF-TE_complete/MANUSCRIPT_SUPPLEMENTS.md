# EDF-TE — Materials and Methods, Figure Legends, and Source Data

> **Notation.**  Species codes: CLADE_SpeciesID (e.g. MU_Mucci3 = *Mucor circinelloides* v3).
> `ln` = natural log; `log₂` = base-2 log.  Micromamba environments are defined in `config.sh`.
> Database release dates are given at first use.

---

## PART I — MATERIALS AND METHODS

### 1. Genome dataset and taxonomic sampling

We assembled a dataset of 320 eukaryotic genomes spanning approximately 1 billion
years of divergence across 18 broad lineages: Holozoa (HO), Nucleariida (NU),
Aphelida (AP), Cryptomycetes (CR), Microsporidia (MI), Blastocladiomycetes (BL),
Olpidiales (OL), Chytridiomycetes (CH), Neocallimastigomycota (NE), Bacillariophyta
(BS), Entomophthorales (EN), Zoopagomycota (ZO), Kickxellomycota (KX), Glomeromycota
(GL), Mortierellomycota (MO), Mucoromycota (MU), Basidiomycota (BA), and Ascomycota
(AS).  Seven Neocallimastigomycota species (CH_Anasp1, CH_Neocon1, CH_Neolan1,
CH_Orpsp1, CH_Pecora1, CH_Caecom1, CH_Pirfi3) carry the JGI "CH" portal prefix but
were classified in the NE clade for all biological analyses.

Genome assemblies were downloaded from JGI MycoCosm (accessed January 2024), NCBI
RefSeq (release 222, January 2024), and the European Nucleotide Archive (ENA).
Scaffolds shorter than 200 bp were discarded before all downstream analyses.
Species were included only when a gene annotation (GFF/GFF3) and a predicted
proteome (FASTA of longest isoforms) were available.  Functional domain annotations
were either downloaded from JGI or recomputed with InterProScan 5.67-99.0 using
the following databases: Pfam 36.0, PRINTS, PANTHER 18.0, TIGRFAM 15.0, and CDD.

Three Mucoromycota assemblies (MU_Mucci3, MU_Rhimi59, MU_Phybl4) were generated de
novo.  PacBio HiFi CCS reads (mean read quality ≥ Q20; mean read length ≥ 10 kb)
were assembled with Flye 2.9 (`--genome-size` estimated from flow cytometry;
`--min-overlap 5000`; `--threads 32`).  Assembly completeness was assessed with
BUSCO 5.8.2 against the fungi_odb10 dataset (lineage_dataset v10, 758 orthologs).
Only assemblies with BUSCO completeness ≥ 70 % were retained (all three new assemblies
exceeded 90 %).  Repeat-masked scaffolds (soft-masking only) were generated with
RepeatMasker 4.1.5 using species-specific libraries (see Section 5) for all downstream
gene prediction.

---

### 2. Species tree inference

#### 2.1 Orthogroup inference

Single-copy orthogroups were identified with OrthoFinder 3.1 run on the 238 annotated
proteomes with the following parameters:

```
orthofinder -f <proteomes_dir> -M msa -A mafft -T fasttree -S diamond -t 32 -a 8
```

The `-M msa` flag instructs OrthoFinder to generate multiple sequence alignments for
each orthogroup; `-A mafft` and `-T fasttree` specify MAFFT v7.520 and FastTree 2.1.11
as the alignment and tree-estimation engines for the intermediate gene-tree step;
`-S diamond` uses DIAMOND v2.1.8 for the all-against-all protein similarity search
(`--more-sensitive`, e-value ≤ 10⁻⁵).

#### 2.2 Marker selection

From the complete OrthoFinder gene families, markers were retained only if they
simultaneously satisfied all four criteria:

1. **Taxon occupancy** ≥ 40 % (≥ 95 of 238 genomes).
2. **Copy number** ≤ 2 copies in every sampled genome (to limit hidden paralogy).
3. **Post-trimming alignment length** ≥ 60 amino acid columns.
4. **Parsimony-informative sites** (PIS, counted from trimmed alignment) ≥ 25.
5. **Relative composition variability** (RCV; PhyKIT v2.0) ≤ 0.5, computed as:
   $\text{RCV} = \frac{1}{4L}\sum_{i=1}^{20}|f_i - \bar{f}|$,
   where $f_i$ is the frequency of amino acid $i$ and $L$ is the alignment length.

Per-gene trees under LG+G4 (IQ-TREE 2.3.6, `--fast` mode) were used to detect
long-branch outlier sequences: any taxon whose terminal branch length exceeded the
per-gene mean by more than 3.5 standard deviations was flagged as a potential
paralog, removed, and the gene was re-evaluated against criteria 1–5.  A total of
682 markers passed all filters, yielding a concatenated supermatrix of 238 taxa ×
219,476 amino acid positions (58.3 % completeness after accounting for missing data).

#### 2.3 Phylogenetic inference

Each retained marker was individually aligned with MAFFT v7.520 (`--auto`; for
alignments ≥ 200 sequences, `--retree 2 --reorder` was applied to reduce
computational load).  Alignments were trimmed with trimAl v1.4.1 (`-automated1`),
removing columns with more than 20 % gaps.  The trimmed alignments were concatenated
and partitioned by gene.

A first-pass guide tree under the LG+F+G4 model (IQ-TREE 2.3.6, `-B 1000`) was
computed to seed the final analysis.  The definitive tree was inferred under the
LG+C60+F+G4 mixture model (`-m LG+C60+F+G4 -B 1000 -bnni --seed 12345`), which
uses 60 empirical amino acid frequency profiles (C60) estimated from a large database
of protein families to account for across-site compositional heterogeneity.  The
`-bnni` flag enables branch-length optimisation of each bootstrap tree.  The tree was
rooted on the Holozoa clade using midpoint rooting followed by manual inspection.

---

### 3. Genomic feature characterisation

Six per-genome metrics were computed from the assembly FASTA and the GFF3 annotation
file (`section1_genome_phylogenetics/02_genomic_features.py`):

| Metric | Operational definition |
|--------|----------------------|
| Genome size (Mb) | Sum of scaffold sequence lengths ÷ 10⁶; N bases counted but not excluded from the total |
| GC content (%) | (G + C) / (A + T + G + C) × 100; N and ambiguous bases excluded from denominator |
| Gene count | Number of records with type = "gene" or type = "mRNA" in the GFF3; only protein-coding features retained |
| Mean gene length (bp) | Mean of (chromEnd − chromStart + 1) for the single longest isoform per gene locus |
| Intron density | (total intron features) / (total gene features) from the GFF3; 0 for intron-poor genomes with no "intron" GFF3 lines |
| N50 (Mb) | Scaffold length at which the cumulative sum of scaffold sizes ≥ 50 % of total assembly size |

Pseudogenes (GFF3 attribute `pseudo=true`) and transposable-element gene models
(any gene whose product description matched the terms "transposon", "retrotransposon",
"reverse transcriptase", "integrase", or "transposase") were excluded from gene count
and mean gene length calculations.  For genomes where the GFF3 did not contain
"intron" records (i.e. introns were implicit between consecutive exon records), intron
density was computed as: (total exon count − total gene count) / total gene count.
From the per-species genomic feature table, a scatter plot of exon content (% exonic, x-axis)
vs. intron content (% intronic, y-axis) was produced with species coloured by clade
(Fig. 1C; `GenomicFeatures/scripts/02_plot_exon_intron_scatter.py`).  Average intron length
per gene in multi-exon genes was visualised as a box plot per clade (Fig. S1A;
`GenomicFeatures/scripts/01_plot_genomic_features.py`).  Pairwise Pearson correlations between
TE content metrics and genomic features across 238 species were visualised as a heatmap
(Fig. S3B; `GenomicFeatures/scripts/03_te_genomic_correlations.py` combined with
`Section2_TElandscape/scripts/D01_merge_te_genomics.py`).

Pan-genome accumulation curves were generated by 100 random permutations of species-addition
order applied to the OrthoFinder gene presence/absence matrix
(`PanGenome/scripts/01_pangenome_compute.py`).  At each depth (number of species sampled),
orthogroups were classified as core (present in all sampled species), shell (present in ≥ 2 but
< n genomes), or unique/cloud (present in exactly one sampled species).  The mean and standard
deviation across 100 random orderings were calculated per depth for each of the 11 taxonomic
groups and for the combined dataset (Fig. 1B; `PanGenome/scripts/02_plot_pangenome.py`).

---

### 4. Microsynteny analysis by ortholog pair collinearity

Microsynteny conservation across 17 clades of the EDF dataset was estimated from the
frequency of collinear ortholog pairs, following the approach of Sebé-Pedrós *et al.*
(2013) adapted for this broader phylogenetic scope (`Synteny/scripts/`).

#### 4.1 Gene-order extraction

Gene positions were extracted from GFF/GFF3 annotation files for each species.  Per
scaffold, protein-coding `gene` features were sorted by chromosomal position using the
`locus_tag` attribute (NCBI-sourced genomes) or the `Name`/`ID` attribute (JGI-sourced
genomes) as the gene identifier (`Synteny/scripts/01_gff_to_bed.py`).  Scaffolds with
fewer than three annotated genes were excluded.  Before running synteny computations,
species were filtered for assembly quality: assemblies passing N50 ≥ 50 kb OR BUSCO
completeness ≥ 80 % were included (`Synteny/scripts/00_genome_qc_filter.py`).

#### 4.2 Ortholog identification

Single-copy ortholog pairs for each within-clade pairwise comparison were retrieved from
the OrthoFinder v2 pairwise ortholog tables (Section 2).  Only 1:1 ortholog relationships
(one gene per species per orthogroup) were retained; orthogroups with multiple copies in
either species were discarded to reduce the influence of recent gene duplications.

#### 4.3 Collinear pair identification

For each species pair (i, j) within the same clade, collinear ortholog pairs were
identified as follows.  On each scaffold of species i, the ordered sequence of orthologs
shared with species j was extracted.  For each consecutive pair of orthologs on the same
scaffold of species i, we tested whether their corresponding orthologs in species j were
also adjacent on the same scaffold, irrespective of gene orientation (to account for local
inversions; Putnam *et al.*, 2007).  The count of such collinear pairs is $c_{ij}$; the
total number of adjacent-ortholog slots across all scaffolds of species i is $N_{ij}$.

#### 4.4 Null distribution and spurious synteny correction

To account for spurious gene-order conservation from random co-localisation, random
genomic positions were assigned to each gene using `bedtools shuffle` v2.24 (Quinlan &
Hall, 2010) in 100 independent rounds per species, with scaffolds preserved as location
constraints (`-chrom` flag).  Each replicate used a fixed per-replicate random seed.
The mean over 100 replicates ($\bar{s}_{ij}$) was computed for each species pair.

#### 4.5 Corrected synteny ratio

The raw synteny score was computed as:

$$\text{score}_{ij} = \frac{c_{ij} - \bar{s}_{ij}}{N_{ij}}$$

Negative values were set to zero.  Scores were normalised to [0, 1] within each clade
using the maximum observed score in that clade as reference:

$$r_{ij} = \frac{\text{score}_{ij}}{\max_{kl}(\text{score}_{kl})}$$

Normalisation was performed independently per clade because absolute synteny scores are
not comparable across clades with different divergence scales.  A per-species ranking by
distance-corrected synteny ratio is shown in Fig. S1B.  Per-clade boxplots of synteny
ratios are shown in Fig. S1C (`Synteny_guide_corrected/figures/synteny_ratio_by_clade_box.pdf`).

#### 4.6 Circos visualisation of selected pairs

For four species pairs with the highest synteny ratios, syntenic gene-pair links were
visualised on the top-30 scaffolds per species using Circos v0.67 (Krzywinski *et al.*,
2009).  The four Circos panels are shown in Fig. 1D
(`Synteny_guide_corrected/figures_circos_v6_4pairs/`).

---

### 5. TE annotation with EarlGrey

Transposable elements were annotated in all 320 genomes with EarlGrey v1.6
(`section2_te_landscape/01_te_annotation.sh`), which orchestrates the following
sequential sub-modules:

#### 5.1 De novo consensus discovery — RepeatModeler2 v2.0.5

RepeatModeler2 was run in RECON + RepeatScout mode with NCBI BLAST+ as the search
engine (`-engine ncbi`):

```
RepeatModeler -database <species_db> -engine ncbi -pa 16 -LTRStruct
```

`-LTRStruct` activates the LTRharvest + LTR_retriever structural LTR discovery
pipeline.  LTRharvest (GenomeTools v1.6.2) was invoked with the following parameters:
`--similar 85 --vic 10 --seed 20 --seqids yes`.  LTR_retriever v2.9.0 filtered
LTRharvest candidates requiring: LTR length 100–7,000 bp; TSD length 4–6 bp; PBS
and PPT motifs when detectable; and identity between 5′ and 3′ LTR ≥ 60 %.
RepeatScout v1.0.6 was used with default parameters; tandem arrays detected by
Tandem Repeat Finder v4.09 (TRF; `2 7 7 80 10 50 500 -f -d -m -h`) were excluded
from the RepeatScout output.  RECON v1.08 multi-locus assembly was applied to
RepeatScout seeds with minimum element count ≥ 8 (RECON `minnoels = 8`).

#### 5.2 Structural and quality filtering

After initial consensus generation, EarlGrey applied the following sequential filters:
- Consensi shorter than 80 bp or longer than 40,000 bp were discarded.
- For putative DNA transposons: terminal inverted repeats (TIRs) were verified with
  IRScan (EarlGrey internal); consensi lacking TIRs of length ≥ 10 bp with identity
  ≥ 70 % were reclassified as Unknown.
- For putative LTR retrotransposons: LTR boundary integrity was verified by
  LTR_retriever v2.9.0; consensi failing the structural filter were reclassified.
- Consensi with > 50 % simple or low-complexity sequence (identified by RepeatMasker
  `-nolow` scoring) were removed.

#### 5.3 Phylogenetic classification — Dfam nhmmer search

Each consensus was classified against the Dfam 3.7 database (released July 2023;
20,462 HMM profiles) using nhmmer (HMMER 3.3.2):

```
nhmmer --cpu 8 -E 1e-5 --tblout <out.tbl> Dfam3.7.hmm <consensi.fa>
```

The best-scoring Dfam hit determined the superfamily-level classification.  Consensi
with no Dfam hit (E > 10⁻⁵) were passed to RepeatClassifier (RepeatModeler2 v2.0.5)
for automated kmer-based assignment.  Consensi that remained unclassified after both
steps were labelled `Unknown`.

#### 5.4 RepeatMasker genome annotation

The curated per-species consensus library was used to annotate the repeat content of
the corresponding genome assembly with RepeatMasker 4.1.5:

```
RepeatMasker -lib <species_library.fa> -pa 16 -s -dir <out_dir> \
  -gff -inv -a -nolow -norna <genome.fa>
```

`-s` applies the slow/sensitive search mode; `-a` produces the alignment (`.align`)
output file required for Kimura divergence computation (Section 6); `-nolow` excludes
the scoring of simple repeats and low-complexity sequences from masking; `-norna`
excludes small RNA genes.  The output `.filteredRepeats.bed` file records one BED
line per annotated TE copy (chromosome, start, end, family, score, strand).

#### 5.5 Cross-species library construction and host-gene filtering

After per-species annotation, all consensus sequences from all 320 species were
concatenated into a joint FASTA file (headers format: `SPECIES__family#class`).
Within-library redundancy was removed by MMseqs2 v14.7564:

```
mmseqs easy-cluster <all_consensi.fa> <clusterDB> /tmp \
  --min-seq-id 0.80 -c 0.80 --cov-mode 1 -s 7 --threads 16
```

Cluster representatives were retained.  Each representative was then compared to
the RefSeq non-redundant protein database (NCBI NR, release 221, December 2023)
by DIAMOND BLASTX v2.1.8 (`--more-sensitive -e 1e-5 --id 40 -k 1`).  Sequences
with a RefSeq hit (identity ≥ 40 %, e-value ≤ 10⁻⁵) not matching any known TE
protein (i.e. the RefSeq hit annotation did not contain the terms "transposon",
"transposase", "integrase", "reverse transcriptase", "capsid", or "gag") were
classified as host-derived and permanently removed from the TE library.  TE
superfamily nomenclature follows Wicker et al. (2007).

All EarlGrey SLURM jobs were submitted to partition `normal.168h` with 16 CPUs and
64 GB RAM per species.  Two Microsporidia genomes (MI_Encro1, MI_Ordco1) returned
zero TE content after filtering and were excluded from ordination analyses but
retained in all other comparative analyses.

---

### 6. Kimura divergence landscapes and the Recency Index

#### 6.1 K2P distance computation

RepeatMasker produces a per-copy pairwise alignment in the `.align` output file.
Each alignment block records the substitution matrix decomposed into matches (m),
transitions (i) and transversions (v) accumulated between the inserted copy and the
consensus sequence.  For every annotated copy we derived the Kimura two-parameter
(K2P) nucleotide substitution distance:

$$K = -\frac{1}{2}\ln(1 - 2p - q) - \frac{1}{4}\ln(1 - 2q)$$

where $p = i / L_{\text{aln}}$ is the transition fraction and
$q = v / L_{\text{aln}}$ is the transversion fraction, and $L_{\text{aln}}$ is
the number of aligned sites after stripping terminal gaps from both the copy and
the consensus.  CpG-context transitions (C→T or G→A at a CpG dinucleotide) are
reported separately by RepeatMasker; we used the CpG-adjusted substitution matrix
(columns "CpG" in the `.align` statistics line) to obtain a corrected K value that
does not inflate divergence estimates in CpG-rich elements.  Copies with fewer than
10 aligned sites after gap-stripping were discarded.

(`section2_te_landscape/04_kimura_divergence.py`)

#### 6.2 Divergence landscape histograms

Kimura divergence landscapes were constructed as stacked histograms per species.
Each histogram spans K ∈ [0, 1] with a bin width of 0.005 (200 bins).  The y-value
per bin is the sum of genomic footprints (end − start in bp from the `.filteredRepeats.bed`)
of all TE copies whose K falls within that bin, stacked by TE class using the
canonical Wicker et al. class hierarchy (Class I LTR, Class I non-LTR, Class II,
Helitrons, Unknown).

#### 6.3 Recency Index

The Recency Index for species $s$ quantifies the proportion of the TE landscape
that has been recently active:

$$RI_s = \frac{\displaystyle\sum_{\substack{i \in s,\; K_i < 0.10}} bp_i}{\displaystyle\sum_{i \in s} bp_i}$$

where $i$ indexes all annotated TE copies in species $s$, $bp_i$ is the copy
genomic footprint, and $K_i$ is its K2P distance.  The threshold $K < 0.10$ was
chosen empirically to capture copies that have accumulated fewer than ~10 % of the
substitutions expected under neutral decay since the origin of the host lineage.
At typical fungal synonymous substitution rates (~2 × 10⁻⁹ per site per year), a
K value of 0.10 corresponds roughly to insertion events occurring within the last
~50 Myr of the analysed divergence time window.

---

### 7. TE landscape composition analysis

#### 7.1 PCA of TE high-level class composition

TE class composition per species was summarised from the EarlGrey `highLevelCount.txt`
output (eight high-level classes: DNA, LTR, LINE, SINE, Helitron, RC, PLE, Unknown) as
proportions summing to 1.  Principal component analysis (PCA) was applied using
`sklearn.decomposition.PCA` on the centred proportion matrix.  Outlier species were
flagged by univariate Z-score (total TE%) and multivariate Mahalanobis distance (TE class
composition vector) within each clade; the Mahalanobis approach falls back to Euclidean
Z-score when clade n < 4 (`Section2_TElandscape/scripts/E01_clade_outliers.py`).  The
first two principal components are shown in Fig. 2B
(`Section2_TElandscape/clade_outliers/clade_composition_pca_wide.pdf`), coloured by
clade.  Per-clade boxplots of total TE content (Fig. S2A), stacked bars for EDF-only
clades (Fig. S2B), and Shannon diversity boxplots (Fig. S2C) are in
`Section2_TElandscape/clade_outliers/`.

#### 7.2 Shannon diversity and Recency Index

Shannon diversity (H′, bits) and Pielou's evenness (J) per species were computed from the
EarlGrey `familyLevelCount.txt` file, excluding structural RNA and repeat classes
(Simple_repeat, Low_complexity, Satellite, rRNA, tRNA, snRNA, scRNA, sno, srpRNA):

$$H_i = -\sum_j p_{ij} \log_2 p_{ij}, \qquad J_i = \frac{H_i}{\log_2 n_i}$$

where $p_{ij}$ is the proportional abundance (bp) of superfamily $j$ in species $i$,
and $n_i$ is the number of detected superfamilies
(`Section2_TElandscape/scripts/B01_shannon_diversity.py`).

The Recency Index (RI) was computed from EarlGrey Kimura divergence data as the fraction
of TE bp with K2P distance < 0.10 (Section 6.3;
`Section2_TElandscape/scripts/C01_kimura_metrics.py`).

Both H′ and RI are displayed as annotation tracks on the cladogram in Fig. 2A
(`Section2_TElandscape/figures/cladogram_te_types_shannon_ri.pdf`).  The RI is shown as
a boxplot per clade in Fig. 2C
(`Section2_TElandscape/figures/kimura_ri_by_clade_boxplot.pdf`) and as a clustered
summary across superfamilies and species in Fig. 2D
(`Section2_TElandscape/figures/kimura_compact_summary_clustered.pdf`).
Shannon boxplots per clade are shown in Fig. S2C
(`Section2_TElandscape/figures/shannon_boxplot_per_clade.pdf`).

---

### 8. DNMT methyltransferase machinery conservation

#### 8.1 Domain-based DNMT classification

DNMT1-type, DIM2-type, and RID/DNMT2-type DNA methyltransferases were identified in all
320 predicted proteomes using InterProScan 5.67-99.0 (Pfam 36.0 database).  For each
protein, the presence of the following Pfam domains was recorded: C5MT catalytic domain
(PF00145), RFTS replication foci targeting sequence, BAH bromo-adjacent homology,
DNAmet_Nter N-terminal methyltransferase extension, Helicase_C, SNF2_N, and PHD finger.
Per-species copy numbers and domain architectures were compiled in `dnmt_matrix_final.tsv`
(`TE_dynamics/dnmt_machinery/data/`).  The DNMT1 vs. DIM2 vs. RID classification was
based on combinations of these domains following Iyer *et al.* (2011).

The species tree was annotated with per-species DNMT copy numbers and TE content
(% genome) from `dnmt_matrix_final.tsv` and `te_content_wide.tsv` (Fig. 3A;
`TE_dynamics/dnmt_machinery/figures/dnmt_te_figure_rotated.pdf`).  The scatter of
DNMT1/DIM2 copy number vs. total TE content is shown for Mucoromycota species in
Fig. S5A (`TE_dynamics/dnmt_machinery/dnmt1dim2_vs_te_content_MU.pdf`).

---

### 9. 5-methylcytosine (5mC) analysis

#### 9.1 WGBS library processing

Paired-end WGBS reads (Illumina NovaSeq 6000, 2 × 150 bp) were quality-trimmed with
Trim Galore v0.6.10 (Cutadapt v4.6 backend):

```
trim_galore --paired --quality 20 --phred33 \
  --clip_r1 10 --clip_r2 10 --three_prime_clip_r1 10 --three_prime_clip_r2 10 \
  --fastqc <R1.fq.gz> <R2.fq.gz>
```

The additional 10 bp clipping from each read end removes fill-in artefacts common
in WGBS libraries.  Trimmed reads were aligned to the reference genome with Bismark
v0.24.2, which performs simultaneous in silico bisulphite conversion and alignment:

```
bismark --genome <genome_dir> --non_directional --score_min L,0,-0.6 \
  -p 4 --no_mixed --no_discordant \
  -1 <R1_trimmed.fq.gz> -2 <R2_trimmed.fq.gz>
```

`--non_directional` is applied because our libraries were prepared with the random
priming PBAT (post-bisulphite adapter tagging) protocol, which produces reads from
all four DNA strands.  `--score_min L,0,-0.6` sets the minimum alignment score as
a linear function of read length (intercept = 0, slope = −0.6), permitting ≤ 6 %
mismatch in a 100 bp read.  `--no_mixed --no_discordant` requires both mates of
each pair to align concordantly.  Only samples with bisulphite conversion efficiency
≥ 99 % (estimated from lambda phage spike-in reads mapping to the non-converted
lambda genome) and unique mapping rate ≥ 40 % were retained for downstream analyses.

Duplicate reads were removed with `deduplicate_bismark --paired`.  CpG-context
methylation calls were extracted with:

```
bismark_methylation_extractor --paired-end --comprehensive \
  --CpG --cytosine_report --genome_folder <genome_dir> \
  --buffer_size 10G <deduplicated.bam>
```

The `--cytosine_report` flag produces a six-column file (chromosome, position,
strand, count_methylated, count_unmethylated, context) for every cytosine in the
genome, regardless of coverage.

#### 9.2 Per-copy methylation quantification

The cytosine report was intersected with the TE copy BED file using:

```
bedtools intersect -wa -wb -sorted \
  -a <te_copies.bed> -b <cytosine_report.bed>
```

For each TE copy, the methylation level was computed from covered CpG positions only:

$$m_i = \frac{\sum_{k \in i,\; \text{cov}_k \geq 5} C_k}{\sum_{k \in i,\; \text{cov}_k \geq 5} (C_k + U_k)}$$

where $C_k$ and $U_k$ are the methylated and unmethylated call counts at CpG position
$k$ within copy $i$, and the sum runs over positions with coverage ≥ 5×.  Copies
with fewer than 3 such covered CpGs were excluded.  A copy was classified as
methylated if $m_i \geq 0.50$.

#### 9.3 5mC landscape by TE family and Kimura age

5mC levels per TE family stratified by Kimura divergence age class were computed as a
coverage-weighted mean per family per age bin, extending the approach of Chiapello *et al.*
to additional Mucoromycota species
(`TE_dynamics/20_chiapello_extension/`; Fig. 3B = `F4_wm_by_family_age.pdf`).  Spearman's
rank correlation between $m_i$ and $K_i$ (per copy) was computed per species.  A separate
analysis for dark-mycelium conditions is shown in Fig. S5B
(`TE_dynamics/12_te_age_5mc/plots/01_wm_by_age_class.pdf`).

#### 9.4 Gene proximity and expression effects

The distance from each gene to the nearest full-length TE insertion was computed and
stratified by methylation quartile (Q1–Q4 of per-copy $m_i$; Fig. 3C:
`TE_dynamics/31_meth_bins_gene_expr/gene_dist_from_TE_edge_FL.pdf`).  Adjacent-gene
expression (RNA-seq RPKM; see Section 10.4) per methylation quartile is shown in
Fig. 3D (top panels: *Phycomyces blakesleeanus*;
`TE_dynamics/31_meth_bins_gene_expr/adj_gene_expr_by_meth_quartile.pdf`).

#### 9.5 Methylation metaplot across TE bodies

5mC profiles across TcMar-Tc1 and MULE-MuDR superfamily insertions were computed as
metaplots stratified by methylation quartile bins (Q1 vs. Q4), showing methylation
signal across TE bodies and flanking regions
(`TE_dynamics/F1_meth_q1q4_metaplot/F1_meth_q1q4_metaplot.pdf`; Fig. 3H).

---

### 10. H3K9me3 ChIP-seq and chromatin-category analysis

#### 10.1 ChIP-seq read processing

ChIP-seq reads (Illumina NovaSeq 6000, 2 × 75 bp) for H3K9me3 (immunoprecipitation)
and the corresponding total-chromatin input control were independently quality-trimmed
with Trim Galore v0.6.10 (`--quality 20 --paired --fastqc`) and aligned to the
reference genome with Bowtie2 v2.5.1:

```
bowtie2 --very-sensitive -X 2000 -p 16 \
  --no-mixed --no-discordant \
  -x <genome_bt2_index> \
  -1 <R1_trimmed.fq.gz> -2 <R2_trimmed.fq.gz> \
  | samtools sort -@ 8 -o <sorted.bam>
samtools index <sorted.bam>
```

`-X 2000` extends the maximum insert size to 2,000 bp to accommodate mononucleosome-
protection libraries (150–200 bp) as well as occasional di-nucleosome fragments.
`--no-mixed --no-discordant` retains only concordantly aligned pairs.  PCR duplicates
were marked and removed with `samtools markdup -r`.  Only samples with ≥ 10 million
uniquely mapped read pairs and Fraction of Reads in Peaks (FRiP) ≥ 0.01 (after
peak calling; see below) were retained.

#### 10.2 Peak calling

H3K9me3 peaks were called with MACS2 v2.2.9.1, using the input-control BAM for
background correction:

```
macs2 callpeak \
  -t <ip.bam> -c <input.bam> \
  --nomodel --extsize 200 \
  -g <effective_genome_size> \
  --bdg -B -q 0.05 \
  --keep-dup auto \
  -n <sample_prefix> --outdir <out_dir>
```

`--nomodel --extsize 200` bypasses the fragment size estimation step and fixes the
fragment extension to 200 bp, which is appropriate for H3K9me3 ChIP-seq with low
signal-to-noise.  `-q 0.05` sets the minimum FDR (Benjamini–Hochberg corrected
q-value) for peak calling.  The `-B` and `--bdg` flags produce per-base coverage
bedgraph files for both the treatment and control channels.

#### 10.3 Per-TE-copy H3K9me3 enrichment

H3K9me3 enrichment per TE copy was computed as the mean bedgraph signal over the copy
BED coordinates, normalised to RPKM:

$$\text{RPKM}_i = \frac{\text{reads in copy } i \times 10^9}{\text{copy length (bp)} \times \text{total uniquely mapped reads}}$$

Per-species per-modification enrichment analyses are shown in Fig. S5C
(`TE_dynamics/36b_enrichment/`), and the top-6 most enriched families by H3K9me3 or
5mC are shown in Fig. S5D
(`TE_dynamics/38_h3k9_5mC_top6_family/38_h3k9_5mC_top6_family.pdf`).
The Venn summary of epigenetic marks across species is shown in Fig. S5E
(`TE_dynamics/36_epigenetic_summary/plots/36d_venn_summary.pdf`).

#### 10.4 RNA-seq for nearest-gene expression

Total RNA-seq reads (Illumina, 2 × 150 bp, rRNA-depleted) were trimmed with Trim
Galore v0.6.10 (`--paired --quality 20`) and mapped with STAR v2.7.11a:

```
STAR --runThreadN 16 \
  --genomeDir <star_index> \
  --readFilesIn <R1.fq.gz> <R2.fq.gz> \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS NM MD \
  --outFilterMismatchNoverLmax 0.05 \
  --alignIntronMax 50000 \
  --outFilterMultimapNmax 1
```

`--outFilterMismatchNoverLmax 0.05` permits at most 5 % mismatches per aligned read;
`--alignIntronMax 50000` accommodates introns up to 50 kb; `--outFilterMultimapNmax 1`
retains only uniquely mapping reads.  Gene-level read counts were obtained with
featureCounts (Subread v2.0.6):

```
featureCounts -p -B -s 0 -T 8 \
  -a <annotation.gtf> -o <counts.txt> <sorted.bam>
```

`-s 0` specifies strand-agnostic counting; `-B` requires both mates of a pair to be
counted.  RPKM was computed as: reads × 10⁹ / (gene length in bp × total reads).

Gene distance from nearest TE insertion stratified by H3K9me3 enrichment level is
shown in Fig. 3E (`TE_dynamics/33_h3k9_bins_gene_expr/gene_dist_from_TE_edge_H3K9.pdf`).
Adjacent gene expression by H3K9me3 peak presence/absence is shown in Fig. 3F
(`TE_dynamics/33_h3k9_bins_gene_expr/adj_gene_expr_by_h3k9_peak.pdf`).
The methylation-state and H3K9me3 summary of TE families at MAC transcription start
sites is shown in Figs. 3G and 3I
(`TE_dynamics/45_mac_tss_expression/` and `TE_dynamics/44b_mac_tss_families/`).
Adjacent gene expression and methylation scatter for individual TE families is shown
in Fig. S6 (`TE_dynamics/48_adjacent_genes/48_family_summary_scatter_v2.pdf`).

---

### 11. Computational nucleosome occupancy prediction

#### 11.1 nuclpred genome-wide prediction

Nucleosome occupancy per genome was predicted computationally using nuclpred
(Minary & Levitt, 2014) applied to the repeat-masked genome FASTA of each species.
nuclpred generates a per-base predicted occupancy WIG file.  The analysis was applied
to 217 species with sufficiently contiguous assemblies
(`Section2_TElandscape/nucleosomes/nuclpred/{SP}_nuclpred.wig`).

#### 11.2 Chromatin state classification

For each species, the genome was tiled into 500 bp windows and characterised by two
metrics from the nuclpred WIG: (i) mean occupancy and (ii) autocorrelation value
(ac_val, indicating nucleosome phasing/regularity).  Four mutually exclusive chromatin
states were assigned per window
(`Section2_TElandscape/scripts/chromatin_states_from_nuclpred.py`):

| State | Criteria |
|---|---|
| Open | mean_occ < q₃₃ (genome-wide) |
| Phased | mean_occ ≥ q₃₃ AND ac_val ≥ median |
| Fuzzy_mid | q₃₃ ≤ mean_occ < q₆₇ AND ac_val < median |
| Fuzzy_high | mean_occ ≥ q₆₇ AND ac_val < median |

`Fuzzy_high` serves as a heterochromatin-like proxy: in *Rhizopus microsporus* (where
in vivo H3K9me3 data is available), H3K9me3-enriched regions are characterised by high
mean occupancy and low phasing — consistent with the `Fuzzy_high` state definition.

#### 11.3 TE occupancy relative to background

For each species, the mean predicted nucleosome occupancy of TE-overlapping 500 bp
windows was compared to background (non-TE, non-genic) windows.  Delta occupancy
(Δocc = TE − background) per TE superfamily was summarised per species
(`Section2_TElandscape/nucleosomes/f3_bias/{SP}_f3_bias.tsv`; 205 species with
columns: species, superfamily, mean_te_occupancy, background_occupancy, delta_occupancy,
n_insertions).  Cross-species and cross-clade summaries are shown in Fig. S4B
(`figS_nuclpred_chromatin_withpoints.pdf`).  Comparison with in-vivo chromatin metrics
from four experimentally characterised species is shown in Fig. S4A
(`figF2_invivo_multimetric_withpoints.pdf`).

---

### 12. N6-methyladenine detection and small RNA silencing

#### 12.1 6mA detection via PacBio kinetics

PacBio CLR (Continuous Long Reads) were mapped to the reference genome with
pbmm2 v1.9.0:

```
pbmm2 align <reference.fa> <reads.bam> <aligned.bam> \
  --preset SUBREAD --sort --min-length 1000 \
  -j 16
```

IPD (inter-pulse duration) ratios relative to an in-silico unmethylated kinetic
model were computed with ipdSummary (SMRTLink v12.0):

```
ipdSummary <aligned.bam> --reference <reference.fa> \
  --identify m6A \
  --minCoverage 25 \
  --pvalue 0.001 \
  --outfile <kinetics.gff>
```

Genomic positions were called as 6mA if the IPD ratio ≥ 5.0 AND per-base sequencing
coverage ≥ 25× AND the ipdSummary p-value ≤ 0.001.  TE copies with ≥ 1 called 6mA
position were classified as 6mA-positive.  The nearest-gene RNA-seq RPKM (from
Section 9.4) was compared between 6mA-positive and 6mA-negative TE copies using a
two-sided Wilcoxon rank-sum test (`section3_te_epigenetics/04_6mA_expression_analysis.py`).

#### 12.2 Small RNA library processing

Small RNA (sRNA) sequencing reads (Illumina HiSeq 2500, 1 × 50 bp single-end)
were adapter-trimmed with Cutadapt v4.6:

```
cutadapt -a TGGAATTCTCGGGTGCCAAGG \
  --minimum-length 18 --maximum-length 32 \
  --discard-untrimmed \
  -o <trimmed.fq.gz> <raw.fq.gz>
```

The adapter sequence (`-a`) corresponds to the Illumina small RNA 3′ adapter.
`--discard-untrimmed` removes reads without a detectable adapter (i.e. reads whose
insert is longer than 50 bp).  Trimmed reads were size-selected to 20–25 nt
(the canonical piRNA/siRNA range in fungi) by additional filtering with
`awk 'NR%4==2{len=length($0)} NR%4==0{if(len>=20 && len<=25) print prev3 "\n" prev2 "\n" prev1 "\n" $0} {prev3=prev2; prev2=prev1; prev1=$0}'`.

Size-selected reads were mapped to the cross-species TE consensus library with
Bowtie v1.3.1:

```
bowtie -v 2 -M 1 --best --strata -p 8 \
  -x <te_library_index> \
  -q <reads_20-25nt.fq> \
  -S <mapped.sam>
```

`-v 2` permits at most 2 mismatches; `-M 1` reports one randomly chosen alignment
for reads mapping to multiple locations (multimappers), setting the `XM` tag; `--best
--strata` ensures the alignment reported is among the best-scoring ones.  Reads
per million (RPM) per TE family:

$$\text{RPM}_j = \frac{N_j}{N_{\text{total mapped}}} \times 10^6$$

where $N_j$ is the count of aligned reads mapped to family $j$ and $N_{\text{total mapped}}$
is the total count of reads with at least one alignment (including multimappers).

Spearman's rank correlation between per-family RPM and mean H3K9me3 RPKM (over all
copies of that family) was computed per species
(`section3_te_epigenetics/06_srna_silencing_pipeline.py`).

---

### 13. TE co-option analysis

#### 13.1 Domain-based co-option screen

Predicted proteomes from all 238 analysed species were annotated with InterProScan
(Pfam database) using a 16-step pipeline (`Cooption/scripts/00_normalize_interpro.py`
through `15_clade_cooption_analysis.py`).  Each gene was classified on the basis of
its Pfam domain composition into one of three categories using a curated list of
TE-machinery Pfam accessions (`Cooption/scripts/05_classify_genes.py`):

- **te_core**: gene carries ≥ 1 Pfam domain corresponding to a TE structural or
  catalytic component.  The te_core domain set includes (non-exhaustively) PF00078
  (RVT_1), PF07727 (RVT_2), PF00075 (RNase_H), PF00665 (rve integrase), PF02022
  (Integrase_H2C), PF03732 (Retrotrans_gag), PF03184 (DDE_1 Tc1/Mariner),
  PF01609 (DDE/IS240), PF10551 (MULE), PF05699 (hAT_N), PF14214 (Helitron_like_N),
  PF03372 (Exo_endo_phos LINE endonuclease), and PF00098 (zf-CCHC).
- **te_assoc**: gene carries the chromodomain (PF00385 Chromo or PF01772
  Chromo_shadow) — a chromatin-reader domain found in Gypsy-clade Ty3 elements —
  but no te_core domain.
- **non_te**: all other genes.

A gene was designated **co-opted** if it satisfied both of the following:
(i) at least one te_core or te_assoc domain is present; and
(ii) the gene is not classified as a complete TE ORF — i.e. it carries at least
one Pfam domain not belonging to the TE domain set
(`Cooption/scripts/06_filter_true_cooption.py`).

A total of **51,799 co-opted genes** were identified across 238 species (17 clades).
[*NOTE — confirm total and per-clade breakdown numbers with the RESULTS_SUMMARY.txt*]

#### 13.2 ChimeraTE chimeric transcript detection

Chimeric gene–TE transcripts were detected with ChimeraTE (Mode 1, RNA-based)
(`Cooption/scripts/chimte_bedtools_perrep.py`, `chimte_merge.py`,
`prep_chimeraTE_data.py`).

[*NOTE — ChimeraTE analysis was running for MU species at time of writing; EDF
species pending.  Please fill in: (1) exact ChimeraTE version and command-line
parameters used; (2) RNA-seq datasets and STAR alignment parameters; (3) final
number of species analysed and chimeric junctions detected.*]

#### 13.3 TF binding site enrichment in TE sequences

DAP-seq peak data for **58 transcription factors** from *Rhizopus microsporus*
ATCC11559 were analysed in two complementary ways
(`TF_TE_analysis/scripts/01_peak_TE_overlap.py` and
`TF_TE_analysis/scripts/03_fimo_motif_enrichment.py`).

**Peak–TE overlap (script 01).**  DAP-seq peak BED files (filtered at
q-value ≥ 10⁻² ) were intersected with the EarlGrey TE annotation of *R. microsporus*
ATCC11559 using bedtools intersect.  For each TF, the fraction of peaks overlapping
any TE copy was computed; enrichment relative to the genomic TE coverage was tested
with a two-sided Fisher's exact test.  Significance was assessed after
Benjamini–Hochberg correction (FDR 5 %).

**FIMO motif enrichment in TE sequences (script 03).**  FIMO v5.5.5 was run for
each of the 58 TF motif files (MEME format; `data/TFsDAP/motifs/`) against the TE
sequences of each TE class, and separately against 3× size-matched genomic
shuffled sequences generated with `bedtools shuffle`.  The p-value threshold for
FIMO hit reporting was 1 × 10⁻⁴.  Enrichment in TE sequences relative to the
shuffled background was tested with a two-sided Fisher's exact test on the
2 × 2 contingency table (TE sequences with ≥ 1 FIMO hit vs. background sequences
with ≥ 1 FIMO hit); p-values were corrected across all TF × TE class combinations
using the Benjamini–Hochberg procedure (FDR 5 %).

A focused cross-species analysis of the HMGb motif (HMGb_1228) enrichment in
LINE/I-Jockey elements was subsequently performed for all available species with
annotated Jockey elements
(`TF_TE_analysis/scripts/run_hmgb3_jockey_cross_species.sh`).

[*NOTE — please confirm: (1) TF family names and DAP-seq experiment details;
(2) the specific HMGb result statistics (p-value, OR, n species);
(3) whether any cross-species analyses beyond R. microsporus were finalised.*]

#### 13.4 TE–host domain fusion proteins

Proteins carrying at least one TE-derived Pfam domain AND at least one host-function
domain in the same gene were identified across all 238 species
(`Cooption/scripts/02_classify_proteins.py`, architecture class `TE_host_fusion`).
This includes proteins in which a TE catalytic component (te_core) has been
co-recruited with a cellular function domain (non_te), and proteins where DUF
(Domain of Unknown Function) domains co-occur with TE domains
(`Cooption/scripts/13_duf_te_fusions.py`).

For DUF × TE pairs, enrichment relative to the null expectation (marginal DUF
frequency × marginal TE frequency) was tested with Fisher's exact test (BH-FDR 5 %).
The most prevalent co-occurring pairs recovered include DUF4217 (PF13959) with
Helitron and Helicase TE domains (present in 156 species), and DUF1605 (PF07717)
with Helitron elements (44 species).

[*NOTE — please provide: (1) the four canonical fusion classes you highlighted in
the paper (domain pair names + Pfam IDs); (2) any additional enrichment statistics
you want quoted.*]

#### 13.5 DAM methyltransferases embedded within TEs

DAM N6-adenine methyltransferase domain proteins (Pfam PF05869) were identified in
all 320 predicted proteomes by HMMER v3.3.2 hmmsearch:

```
hmmsearch --cpu 8 -E 1e-5 --noali \
  --tblout <out.tbl> PF05869.hmm <proteome.faa>
```

`--noali` suppresses the verbose alignment output.  The E-value threshold 10⁻⁵ was
applied to the full-sequence score (not domain score) to exclude fragmentary matches.

Each DAM-domain gene was cross-referenced to the EarlGrey BED coordinates with:

```
bedtools closest -d -a <dam_genes.bed> -b <te_copies.bed>
```

A gene was classified as **TE-embedded** when the distance field returned by
`bedtools closest -d` equalled 0, indicating that the gene model fully overlapped
with at least one annotated TE copy on the same strand.

AMT1-class N6-adenine methyltransferases were defined by IPR014885 (Methyltransferase,
N6-adenine, N-terminal domain).  The 2×2 Fisher's exact test for co-occurrence of
DAM-TE and AMT1 across 320 species yielded **OR = 5.19, p = 0.028** (two-tailed).
Fifteen species harboured confirmed DAM-TE genes: 14 associated with DIRS-class LTR
retrotransposons, 1 with a Gypsy LTR element.

---

### 14. Whole-genome duplication analysis

WGD detection combined three independent lines of evidence computed for each
species: (i) synonymous substitution rate (Ks) distributions of the whole-paranome
(wgd v2), (ii) collinearity-based detection of multiplicons (wgd v2 `syn` module),
and (iii) gene tree reconciliation–based duplication counts from OrthoFinder.
Evidence was integrated into a four-category classification, described in §13.6.
All WGD analyses were restricted to the 131 species with ≤ 50 assembly scaffolds
(to ensure reliable collinearity detection); the full 320-species dataset was used
only for the OrthoFinder reconciliation arm (§13.3).

#### 14.1 Genome pre-filtering for collinearity-based analyses

Before running the wgd v2 synteny and Ks pipelines, assemblies were filtered to
retain only those with ≤ 50 scaffolds/contigs (`scripts/select_genomes_le50_contigs.py`).
Scaffold count was obtained from the FASTA index (`.fai`) generated with `samtools faidx`.
This threshold was chosen because wgd v2 `syn` is computationally intractable and
statistically unreliable on highly fragmented assemblies where true collinear blocks
are interrupted by assembly gaps.  The 131 assemblies passing this threshold cover
all 18 clades and include all species previously described as WGD candidates in the
literature.

#### 14.2 Whole-paranome Ks distribution — wgd v2 `dmd` + Ks pipeline

Whole-paranome all-vs-all protein similarity was computed for each of the 131
selected species with `wgd dmd` (wgd v2; Zwaenepoel & Van de Peer 2019 updated
2023):

```
wgd dmd <proteome.faa> --out-dir <dmd_dir> \
  --threads 16 --eval 1e-5 --inflation 1.5 \
  --alignmethod mafft --tree fasttree
```

`--inflation 1.5` is the MCL inflation parameter controlling granularity of
paralogous family clustering (higher values = smaller, more cohesive families).
`wgd dmd` performs: (1) DIAMOND BLASTP all-vs-all with E ≤ 10⁻⁵; (2) MCL graph
clustering on bit-score normalised similarities; (3) multiple sequence alignment of
each family with MAFFT v7.520 (`--auto`); (4) FastTree 2.1.11 (`-gamma -lg`)
maximum-likelihood tree per family; (5) extraction of all paralogous gene pairs
within each family.

For each paralogous pair extracted from the MCL-clustered families, the CDS
sequences were codon-aligned with PAL2NAL v14:

```
pal2nal <protein_pair_aln.faa> <cds_pair.fa> -output paml -nogap
```

`-nogap` removes any codon triplet containing a gap character.  Ks and Ka values
were estimated from the resulting codon alignment with yn00 (PAML v4.10.10):

```
yn00   # control file: icode = 0 (standard genetic code)
       #               weighting = 0 (unweighted)
       #               commonf3x4 = 0 (pair-specific base frequencies)
```

yn00 implements the Yang & Nielsen (2000) maximum-likelihood codon model, which
estimates Ks and Ka by simultaneously optimising the transition/transversion ratio,
codon frequencies, and substitution parameters.  Gene pairs were excluded if:
Ks = 0 (CDS identity; likely assembly duplicates), Ks > 3.0 (saturation under any
standard substitution rate), Ka/Ks > 5 (non-synonymous model convergence failure),
or the codon alignment contained fewer than 100 codons after gap removal.
The resulting per-species paranome Ks distributions were used for GMM fitting
(§13.4) and age-correlation analyses (§13.5).

#### 14.3 OrthoFinder gene tree reconciliation — duplication mapping

All 320 species (including fragmented assemblies) were included in the OrthoFinder
3.1 run described in Section 2.  The `Gene_Duplication_Events/Duplications.tsv`
file from OrthoFinder records every gene duplication event inferred by gene tree /
species tree reconciliation (using the STRIDE algorithm; Emms & Kelly 2018), with
the column `Species Tree Node` indicating the node on the ML species tree (Section 2)
at which the duplication was mapped.

Duplication rates per species-tree node were computed as normalised counts
(`scripts/compute_paper_like_duplication_rates.py`):

$$R_n = \frac{D_n}{T_n}$$

where $D_n$ is the count of duplication events mapped to node $n$ in
`Duplications.tsv`, and $T_n$ is the number of **eligible orthogroups** for node
$n$ — defined as orthogroups in which species from **both** descending clades of
node $n$ are represented, after excluding orthogroups where any single species
contributes > 5 copies (the `--max-species-copies 5` filter, which removes
lineage-specific gene family expansions that would otherwise inflate node-specific
duplication counts).

Species-tree nodes at which the normalised rate $R_n$ was classified as an outlier
(Robust Z-score ≥ 1.5) relative to all nodes of the same type (terminal or
internal) were tagged as **duplication-peak nodes** (`scripts/detect_duplication_peaks.py`).
The Robust Z-score was computed as:

$$Z_n = \frac{0.6745\,(R_n - \tilde{R})}{\text{MAD}(R)}$$

where $\tilde{R}$ is the median of all per-node rates across the species tree, and
$\text{MAD}(R) = \text{median}(|R_n - \tilde{R}|)$.  The constant 0.6745 rescales
MAD to be consistent with the standard deviation under a Gaussian distribution.
Per-species signals were summarised as: `deep_peak_nodes` (count of internal
tree nodes at depth ≥ 2 from the tip classified as duplication-peak), and
`recent_peak_nodes` (terminal/depth-1 nodes classified as duplication-peak), and
`max_peak_z` (maximum $Z_n$ over all nodes of that species).

To complement the node-level analysis, a **triplet topology test** was applied to
each orthogroup containing gene copies from at least one focal species and two
outgroup species forming a triplet in the species tree.  For each triplet
(focal, outgroup_1, outgroup_2), the gene tree topology was extracted from the
OrthoFinder resolved gene trees (`Resolved_Gene_Trees/`).  A topology was
classified as `topology_ok` if the gene subtree of the focal species was
monophyletic with respect to both outgroups — i.e. an unambiguous in-paralogue
group could be distinguished from the outgroup sequences.  The count of `topology_ok`
orthogroups per focal species (`wgd_candidates_by_species.tsv`) was used as a
cross-validation signal for the duplication-peak calls.

#### 14.4 GMM fitting and WGD peak detection in the paranome Ks

The per-species paranome Ks distribution (from §13.2; maximum 5,000 randomly
sub-sampled pairs when more pairs existed, to reduce compute time while preserving
distribution shape) was modelled with a Gaussian Mixture Model (GMM) using
scikit-learn GaussianMixture:

```python
from sklearn.mixture import GaussianMixture
gm = GaussianMixture(
    n_components=k,
    covariance_type='full',
    n_init=5,
    max_iter=200,
    random_state=42
).fit(Ks_values.reshape(-1, 1))
```

The number of components $k$ was varied from 1 to 5.  Model selection used the
Bayesian Information Criterion:

$$\text{BIC}(k) = k \ln n - 2 \ln \hat{L}$$

where $n$ is the number of Ks data points and $\hat{L}$ is the maximum likelihood.
The $k$ minimising BIC was selected.  A GMM component was reported as a candidate
WGD Ks peak if its mixture weight $w_k \geq 0.10$ AND its mean
$\mu_k \in [0.05,\, 3.0]$, thereby excluding near-zero components (gene
conversion / very recent duplicates) and saturated components above the reliable
estimation range.  The event timing was classified as `ancient_deep` if all
WGD-peak components had $\mu_k > 0.50$, `recent_terminal` if $\mu_k \leq 0.50$,
or `mixed_deep_and_recent` if both conditions occurred in different components.

#### 14.5 Collinearity analysis — wgd v2 `syn` and multiplicons

For the 131 filtered species (§13.1), intra-genomic collinearity was assessed using
`wgd syn` (wgd v2):

```
wgd syn <dmd_families_tsv> <gff3> \
  --out-dir <syn_dir> \
  --ks-range 0.05 3.0 \
  --anchor-ks \
  --threads 16
```

`wgd syn` builds a collinearity graph using the MCScan-derived algorithm: it first
identifies **anchor gene pairs** — paralogous pairs within MCL families whose Ks
value falls in the specified range — and then clusters anchors into collinear blocks
by the tandem-chaining heuristic.  Each collinear block is reported as a
**multiplicon** (equivalent to a paralogon): a pair of genomic segments on the same
or different chromosomes sharing ≥ 3 anchor gene pairs in the same relative order
and orientation.  The output files comprise:
- `*.multiplicons.tsv`: one row per multiplicon (genomic coordinates, n_anchors,
  n_segments);
- `*.anchorpoints.tsv`: one row per anchor gene pair within a multiplicon;
- `*.segments.tsv`: one row per genomic segment participating in a multiplicon.

Per-species synteny statistics were aggregated as `n_multiplicons`,
`n_anchorpoints`, and `n_segments` (`scripts/summarize_synteny_outputs.py`).
The synteny evidence score was defined as the maximum of the multiplicon and
segment Robust Z-scores over the 131-species distribution:

$$S_i = \max\!\left(Z_{M,i},\; Z_{S,i}\right)$$

where $Z_{M,i}$ and $Z_{S,i}$ are computed by the same MAD-based Robust Z-score
formula (§13.3) applied to `n_multiplicons` and `n_segments`, respectively.
A species was considered to have **synteny evidence** if `syn_completed = 1` AND
`n_multiplicons > 0`.

The anchor-point Ks distribution (Ks values of anchor gene pairs within multiplicons)
provides a more precise WGD age estimate than the full paranome, because it filters
out background small-scale duplications.  The modal Ks of the anchor-Ks distribution,
extracted as the GMM peak component (§13.4 applied to anchor pairs only), was used
as the primary WGD age estimate.

#### 14.6 Multi-evidence WGD classification

Species were classified into four WGD categories by combining the three evidence
streams (`scripts/classify_wgd_evidence.py`; `scripts/integrate_wgd_evidence_paper_like.py`):

| Classification | Criterion |
|---|---|
| **Strong WGD** | $R_n \geq 0.30$ (duplication signal) AND synteny completed AND $n_\text{multiplicons} > 0$ |
| **Uncertain — duplication only** | $R_n \geq 0.30$ AND ($n_\text{multiplicons} = 0$ OR `syn_completed = 0`) |
| **Weak / ancient WGD** | $R_n < 0.30$ AND $n_\text{multiplicons} > 0$ AND `syn_completed = 1` |
| **No WGD** | $R_n < 0.30$ AND ($n_\text{multiplicons} = 0$ OR `syn_completed = 0`) |

The duplication rate threshold $R_n \geq 0.30$ was calibrated against 12 fungal
species with previously published WGD evidence (*Saccharomyces cerevisiae*, the
Rhizopus microsporus complex, *Phytophthora infestans*, etc.); all 12 scored
$R_n > 0.40$ in our normalised OrthoFinder analysis.  Species with
`event_timing = mixed_deep_and_recent` (§13.4) were flagged for manual review,
as this pattern can arise from two independent rounds of duplication or from
allopolyploidy followed by diploidisation.

#### 14.7 TE breakpoint enrichment at WGD synteny boundaries

For each WGD-positive species (Strong WGD or Weak/ancient), the boundaries of
intra-genomic collinear blocks (wgd v2 multiplicons, §13.5) were used to define
breakpoint windows — i.e. the regions of the genome where multiplicon adjacency
is interrupted.  Concretely, a breakpoint was defined as the chromosomal interval
between the last anchor-point coordinate of one multiplicon and the first
anchor-point coordinate of the next multiplicon on the same scaffold, sorted by
genomic position.  Each breakpoint midpoint was expanded by ±5,000 bp.  One hundred
randomly sampled windows of the same length distribution (sampling with replacement
from the scaffold) were drawn per breakpoint as a background.

TE density in each window $W$:

$$d_W = \frac{\displaystyle\sum_{i \in W} \text{overlap}(TE_i, W)}{|W|}$$

A two-sided Wilcoxon signed-rank test (scipy.stats.wilcoxon) compared matched
breakpoint vs. random TE densities (one observation per breakpoint, paired).
Per-superfamily enrichment: same Wilcoxon test restricted to superfamilies with
≥ 10 windows containing ≥ 1 bp of that superfamily; BH-corrected FDR 5 %.

#### 14.8 TE burst — WGD age correlation

WGD ages in Mya were estimated from the modal Ks of the **anchor-pair** GMM peak
component (§13.5), converted using a per-lineage synonymous substitution rate:

$$t_{\text{WGD}} = \frac{\mu_k^{\text{anchor}}}{2 \times r_{\text{syn}}}$$

where $r_{\text{syn}}$ (substitutions per site per year) was derived from the
time-calibrated divergence of the node connecting the WGD lineage to its closest
non-WGD sister in the ML species tree, using the fossil-constrained nodes of the
companion study.

For each WGD-positive species, the TE burst age was summarised as the KDE mode of
the per-copy Kimura K distribution (§6.2) for each superfamily, using
`scipy.stats.gaussian_kde` with Scott's bandwidth.  Spearman's $\rho$ between
$t_{\text{WGD}}$ and the TE burst KDE mode was computed per superfamily across all
WGD-positive species; BH-corrected FDR 5 % across superfamilies.

#### 14.9 Normalised gene duplication rates per phylogenetic node

The normalised rate $R_n$ (§13.3) was plotted on the species tree as a heat-map
column, using the same colour scale across all nodes.  For visualisation, ancestral
node rates were projected onto the LG+C60+F+G4 ML species tree with ete3 v3.1.3
(`NodeStyle` colouring, `TreeStyle` rectangular layout).  Duplication-peak nodes
(Robust Z ≥ 1.5, §13.3) are highlighted with a filled circle; the size is
proportional to $Z_n$.

---

### 15. Horizontal transposon transfer (HTT)

#### 15.1 HTT candidate identification

Per-clade TE consensus libraries (from Section 5.5, clustered at 95 % identity by
MMseqs2 `--min-seq-id 0.95 -c 0.80 --cov-mode 0`) were constructed and used as
queries in an all-vs.-all BLASTN search within each of the 18 clades:

```
blastn -query <clade_library.fa> -db <clade_genome_db> \
  -perc_identity 85 -word_size 11 -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
  -num_threads 16
```

A TE family was flagged as an HTT candidate if at least one cross-species hit
(query and subject from different species) satisfied:
- nucleotide identity ≥ 85 %
- alignment length ≥ 100 bp
- E-value ≤ 10⁻⁵

**Phylogenetic anomaly score**: For each BLASTN hit pair (species A, species B,
TE family F), the anomaly score was:

$$a = d_{\text{expected}} - d_{\text{observed}}$$

where $d_{\text{expected}}$ is the species-tree patristic distance (from the
LG+C60+F+G4 ML tree, in expected substitutions per site) between A and B, and
$d_{\text{observed}} = 1 - p_{\text{ident}}/100$.  HTT events were retained if
$a \geq 0.10$ (i.e. the TE copies are at least 10 % more similar than expected
from neutral host-genome divergence).

#### 15.2 HTT dating

For each confirmed HTT event, the transfer date was estimated as the median
synonymous Ks among all orthologue pairs between the donor and recipient species
at the relevant OrthoFinder orthogroup, computed with yn00 (PAML v4.10.10) as in
Section 14.1.  Ks was converted to Mya using the per-lineage $r_{\text{syn}}$
from Section 14.4.

#### 15.3 Flanking sequence locus conservation

For each HTT candidate TE locus, 10,000 bp of flanking sequence on each side was
extracted.  Both flanks were queried against the proteome and genome of the putative
donor species via BLASTN (≥ 60 % identity, ≥ 50 bp alignment length, E ≤ 10⁻⁵,
`-task blastn-short`):

$$f_{\text{cons}} = \frac{\sum_{\text{hits}} \text{alignment\_length}_{\text{hit}}}{\text{total flank length (bp)}}$$

Loci were classified as:
- **Conserved**: $f_{\text{cons}} \geq 0.30$
- **Partial**: $0.10 \leq f_{\text{cons}} < 0.30$
- **Non-conserved**: $f_{\text{cons}} < 0.10$

#### 15.4 HTT event counting — Dollo parsimony

HTT events were counted under a Dollo parsimony model (an element can be gained only
once, losses are free) using the ete3 v3.1.3 `TreeNode.get_descendant_evol_events()`
framework.  The full candidate set comprised 224 BLASTN-flagged TE family × clade-pair
combinations.  After applying the anomaly-score filter ($a \geq 0.10$) and the
locus-conservation classification (≥ partial), the parsimony-minimum number of
independent transfer events across the species tree was **11**.

---

### 16. Statistical framework

All statistical tests were two-tailed.  Multiple comparisons were corrected by the
Benjamini–Hochberg (BH) procedure at FDR ≤ 5 % whenever ≥ 5 hypotheses were tested
simultaneously (scipy.stats.false_discovery_control(method='bh')).  Effect sizes
are reported as: Spearman $\rho$ or Pearson $r$ (correlations), odds ratio with
95 % Wilson confidence intervals (Fisher's exact), and rank-biserial correlation $r$
(Wilcoxon / Mann–Whitney tests).

Unless otherwise stated, all analyses were implemented in Python 3.11.3 with the
following library versions: numpy 1.26.2, pandas 2.1.4, scipy 1.11.4, scikit-learn
1.3.2, matplotlib 3.8.2, seaborn 0.13.2, pyBigWig 0.3.18, ete3 3.1.3,
networkx 3.2.1.  R 4.4.0 was used for visualisation with ggplot2 4.0.0, ggtree 3.10,
circlize 0.4.16, and ComplexHeatmap 2.18.0.  Network visualisation was performed with
Cytoscape 3.10.2 with the yFiles Organic layout.


---

## Part II — Figure Legends

### Figure 1. Genomic diversity and synteny across 320 EDF genomes

**(A)** Circular maximum-likelihood phylogenetic tree of all 320 genomes (Section 2;
LG+C60+F+G4 model, IQ-TREE 2).  Tips are coloured by clade following the 18-clade
colour scheme; bootstrap support ≥ 95 is indicated at internal nodes.  The 18 clade
abbreviations are labelled on the outer ring.
Source data: Fig1A_tree.

**(B)** Pan-genome accumulation curves (Section 3) computed from OrthoFinder
orthogroups via 100 random species-addition permutations.  Each line shows the mean
(± 1 s.d.) number of orthogroups classified as pangenome (≥ 1 genome), core (all
genomes), shell (≥ 2 but < n genomes), or unique (1 genome) as a function of genome
count, shown per major taxonomic group.
Source data: Fig1B_pangenome.

**(C)** Scatter plot of exon content (% exonic, x-axis) vs. intron content (%
intronic, y-axis) per genome (Section 3).  Each dot is one species; colour indicates
phylogenetic clade.  Selected outlier species are labelled.
Source data: Fig1C_genomic_features.

**(D)** Circos microsynteny diagrams for four selected species pairs with the highest
corrected synteny ratios (Section 4).  Ribbons link collinear ortholog pairs on the
top-30 scaffolds per species; ribbon colour corresponds to the scaffold in the left
(reference) species.  Species pair labels are shown at the start of each arc.
Source data: Fig1D_synteny.

---

### Figure 2. TE landscape diversity across 18 clades

**(A)** Cladogram of 238 species with TE annotation, with four coloured annotation
tracks: (i) TE high-level class composition (stacked bar, Wicker et al. classes);
(ii) total TE content (% genome, bar); (iii) Shannon diversity H′ (Section 7.2);
(iv) Recency Index RI (Section 6.3).  Tips coloured by clade.
Source data: Fig2A_te_content.

**(B)** PCA of TE high-level class composition across 238 species (Section 7.1).
Each dot is one species; colour indicates clade.  Ellipses: 95% confidence intervals
per clade.  Outlier species (within-clade Z-score or Mahalanobis distance above
threshold) are labelled.  Axis labels include variance explained.
Source data: Fig2B_pca_composition.

**(C)** Boxplot of the Recency Index (RI) per clade (Section 6.3).  Clades sorted
by median RI; each box shows the interquartile range; whiskers: 1.5 × IQR; outliers:
individual points.
Source data: Fig2C_kimura_ri.

**(D)** Clustered Kimura divergence summary across 238 species, showing the mean
Kimura divergence per TE superfamily × species combination as a compact heat-map
(Section 6.2).  Species and superfamilies are ordered by hierarchical clustering.
Source data: Fig2D_kimura_summary.

**(E)** TE proximity and mechanism analysis — panel D from
`expansion_analysis/results/figures/fig2_te_proximity_mechanism.pdf` (Section 14).
[ASK: brief description of what panel D shows]
Source data: Fig2EF_proximity.

**(F)** TE proximity and mechanism analysis — panel C from the same file (Section 14).
[ASK: brief description of what panel C shows]
Source data: Fig2EF_proximity.

**(G)** Sex locus (MAT) analysis — panels A and C of
`MAT_locus/figures/sex_locus_figure.pdf`.
Panel A: genomic architecture of the sex locus in representative species.
Panel C: [ASK: description of panel C].
Source data: Fig2G_sex_locus.

**(H)** Sex locus microsynteny — linear synteny comparison of the sex locus region
between *Fenestella* sp. and *Phycomyces blakesleeanus*
(`sex_locus_synteny_v3_Fenlin.pdf` and `sex_locus_synteny_v3_Phycomyces.pdf`).
Coloured links connect collinear gene pairs flanking the MAT locus.
Source data: Fig2H_sex_locus_synteny.

---

### Figure 3. Epigenetic regulation of transposable elements

**(A)** Conservation of DNA methyltransferase machinery across 320 species (Section 8).
Left: species tree with per-species DNMT1, DIM2, and RID/DNMT2 copy number tracks
(colour-coded circles; absent = white).  Right: total TE content (% genome) bar track
per species.
Source data: Fig3A_dnmt.

**(B)** 5-methylcytosine (5mC) per-copy weighted mean by TE family and Kimura
divergence age class (Section 9.3).  Heat-map rows: TE families with ≥ 50 covered
copies; columns: age bins.  Colour: coverage-weighted mean 5mC.
Source data: Fig3B_5mc_family_age.

**(C)** Cumulative distribution of gene distance from the nearest full-length TE
insertion, stratified by per-copy 5mC methylation quartile (Q1–Q4; Section 9.4).
Each curve is one quartile; shaded bands: ± 1 s.e.m.  Kolmogorov–Smirnov tests
between Q1 and Q4 shown.
Source data: Fig3C_gene_dist_5mc.

**(D)** Expression of genes adjacent to TE insertions, stratified by methylation
quartile (Section 9.4).  Top panels: *Phycomyces blakesleeanus*.  Box plots of
adjacent-gene expression (log₂ RPKM) per quartile; Benjamini–Hochberg-corrected
Mann–Whitney U comparisons shown above brackets.
Source data: Fig3D_adj_gene_expr_meth.

**(E)** Gene distance from nearest TE insertion stratified by H3K9me3 ChIP
enrichment level (Section 10.3).  Format as in (C).
Source data: Fig3E_gene_dist_h3k9.

**(F)** Adjacent gene expression by H3K9me3 peak presence/absence (Section 10.4).
Format as in (D).
Source data: Fig3F_adj_gene_expr_h3k9.

**(G)** Analysis of TE-associated transcription start sites (TSS) across species
(Section 10.4).  [ASK: describe what "MAC" means and exact panel content of
45_mac_tss_expression.pdf for mycelium light condition, Rhizopus + Phycomyces.]
Source data: Fig3G_mac_tss.

**(H)** Metaplot of 5mC levels (Q1 vs. Q4 methylation bins) across TcMar-Tc1 and
MULE-MuDR superfamily insertions (Section 9.5).  Each line represents one
methylation quartile bin; shaded region: TE body; bands: ± 1 s.e.m.
Source data: Fig3H_metaplot.

**(I)** TE families enriched at methylated transcription start sites (Section 10.4).
Bar chart of TE superfamilies most frequently overlapping methylated TSS windows,
per species.
Source data: Fig3I_mac_tss_families.

**(J)** RNAi pathway overview — pie charts showing the proportion of TE superfamilies
targeted by the RNAi/dicer pathway per species (Section 12.2).
Source data: Fig3J_rnai_pies.

**(K)** Small RNA (sRNA) targeting of TE superfamilies — panel D scatter plot of
sRNA abundance (RPM, y-axis) vs. TE copy expression (RPKM, x-axis) per superfamily;
co-regulated TE class counts summarised alongside (Section 12.2).
Source data: Fig3K_srna.

**(L)** IGV genome browser snapshots of representative TE loci with overlaid sRNA
coverage tracks (Section 12.2).  Species and loci selected to illustrate concordant
sRNA targeting of active TE insertions.
Source data: Fig3L_igv.

---

### Figure 4. TE co-option and transcription factor binding

**(A)** Co-option overview tree (Section 13.1).  Species tree of all 238 analysed
genomes with tip symbols indicating presence of co-opted TE-derived genes per clade.
[ASK: clarify exact panel layout of figA_tree_cooption_v2.pdf]
Source data: Fig4A_cooption_tree.

**(B)** Domain composition of co-opted genes in basal fungi (Section 13.1).  Stacked
bar chart of the most prevalent host-function Pfam domains (excluding transposase /
integrase / RT domains) in co-opted proteins per clade, restricted to basal fungal
clades.
Source data: Fig4B_domain_barchart.

**(C)** TE class proportions of co-opted genes per clade (Section 13.1).  Stacked
bars show the fraction of co-opted genes derived from each high-level TE class
(Class I LTR, non-LTR; Class II; Helitron; unknown).
Source data: Fig4C_te_class_proportions.

**(D)** Fusion protein co-option network visualised in Cytoscape (Section 13.4).
Nodes: Pfam domain families — diamonds = host functions, circles = TE domains.
Edges: ≥ 1 protein with both domain types co-recruited.  Edge width: number of
species.  Layout: yFiles Organic.
Source data: Fig4D_fusion_network.

**(E)** ChimeraTE chimeric gene–TE transcript detection in basal fungi (Section 13.2).
Bar chart of chimeric junction counts per species, coloured by TE class.
Source data: Fig4E_chimte.

**(F)** Phylogenetic trees of FLYWCH × TE-RT fusion proteins, restricted to
Neocallimastigomycota (NE clade; Section 13.4).  Trees show the evolutionary
diversification of FLYWCH domain fusions.
Source data: Fig4F_flywch_trees.

**(G)** Phylogenetic trees of CotH × TE fusion proteins across clades (Section 13.4).
Source data: Fig4G_coth_trees.

**(H)** Full circular phylogeny of DAM methyltransferase + reverse transcriptase (RT)
domain proteins (Section 13.5).  Tips coloured by TE superfamily; outer ring: host
species clade.
Source data: Fig4H_dam_rt_tree.

**(I)** Peak enrichment of TF binding (DAP-seq, 58 TFs) in TE classes — bubble chart
(Section 13.3).  Bubble size: number of TF × TE combinations significant after
BH correction; colour: median −log₁₀ FDR.
Source data: Fig4I_tf_enrichment.

**(J)** deepTools metaplot of HMGb3 (HMGb_1228) DAP-seq signal over LINE/I-Jockey
TE insertions (Section 13.3).  Solid line: mean signal; shaded band: ± 1 s.e.m.;
vertical dashed lines: TE body boundaries.
Source data: Fig4J_hmgb3_metaplot.

**(K)** MEME motif logo for the HMGb3 binding motif enriched within LINE/I-Jockey
TE sequences (Section 13.3).  Nucleotide height proportional to information content.

**(L)** Clade-level heat-map of HMGb3 occupancy across Jockey superfamily insertions
(Section 13.3).  Rows: EDF clades; columns: TE subfamilies within Jockey; colour:
mean normalised DAP-seq signal.
Source data: Fig4L_hmgb3_heatmap.

---

### Figure 5. Whole-genome duplication and TE dynamics in Mucoromycota

**(A)** Gene duplication rate tree for Mucoromycota (Section 14.9; ≤ 200-scaffold
filter) annotated with normalised duplication rates per phylogenetic node, with an
adjacent bar chart of co-linear paralog counts per species (Section 14.3).
Source data: Fig5A_paralogs.

**(B)** Syntenic block structure and TE content (Section 14.5).  Ribbon diagram
or dot-plot showing syntenic block coverage relative to total genome length, with
TE density per block overlaid as a colour track, per WGD species.
Source data: Fig5B_syntenic_blocks.

**(C)** TE enrichment at WGD synteny breakpoints — log₂ fold-change over randomised
windows, per TE superfamily, per species (Section 14.7).  Colour: BH-corrected
significance; dot size: absolute log₂FC.
Source data: Fig5C_te_breakpoints.

**(D)** TE dynamics combined figure — panel C: [ASK: describe what panel C of
fig_te_dynamics_combined.pdf shows] (Section 14).
Source data: Fig5D_te_dynamics.

**(E)** TE–Ks age correlation — panel A: Spearman correlation between anchor-pair Ks
(WGD age proxy) and TE superfamily burst mode, for species with reliable Ks estimates
(Section 14.8).  Each point: one species × superfamily; colour: TE class.
Source data: Fig5E_te_ks_corr.

---

### Figure 6. Horizontal transposon transfer (HTT)

**(A)** HTT chord diagrams.  Large panel: all HTT candidate links between EDF clades
(`HTT/figures/htt_chord_labeled.pdf`).  Small inset: validated HTT events only
(`HTT/figures/figHTT_chord_validated.pdf`).  Arc width ∝ number of candidates.
Clade colours follow the 18-clade scheme.
Source data: Fig6A_htt_chord.

**(B)** HTT summary — panel A: Ks distribution of 223 HTT candidates, stacked by
TE class (Section 15.1).
Source data: Fig6B_htt_ks.

**(C)** HTT summary — panel B: heat-map of HTT candidate counts per clade pair
(triangular matrix; Section 15.1).  Cell colour: log₁₀(n + 1) candidates.
Source data: Fig6C_htt_heatmap.

**(D)** HTT LINE synteny — synteny plot showing shared locus architecture (flanking
gene conservation) for validated LINE-class HTT events (Section 15.3).
Source data: Fig6D_htt_line_synteny.


## Part III — Source Data Index

The Excel workbook `EDF-TE_source_data.xlsx` contains one sheet per figure panel.
Each sheet provides the numerical data underlying one or more figure panels as
described below.

| Sheet name | Figure panel(s) | Description | Key columns |
|---|---|---|---|
| Fig1A_tree | Fig. 1A | Species tree tip metadata | species, clade, clade_colour, bootstrap_support |
| Fig1B_pangenome | Fig. 1B | Pangenome accumulation curve data | group, n_genomes, n_core, n_shell, n_unique, n_pangenome, mean, sd, permutation_id |
| Fig1C_genomic_features | Fig. 1C | Per-species genomic feature metrics | species, clade, Pct_Exonic, Pct_Intronic, Avg_Introns_Multiexon, Total_Genes |
| Fig1D_synteny | Fig. 1D | Corrected synteny ratio per species pair | species_A, species_B, clade, c_ij, s_ij_mean, N_ij, score_ij, r_ij |
| Fig2A_te_content | Fig. 2A | TE class composition + Shannon + RI per species | species, clade, total_te_pct, DNA_pct, LTR_pct, LINE_pct, shannon_H, RI |
| Fig2B_pca_composition | Fig. 2B | PCA of TE class composition | species, clade, PC1, PC2, var_explained_PC1, var_explained_PC2 |
| Fig2C_kimura_ri | Fig. 2C | Recency Index per species | species, clade, RI, n_te_bp_total |
| Fig2D_kimura_summary | Fig. 2D | Mean Kimura divergence per superfamily × species | species, clade, superfamily, mean_kimura, n_copies |
| Fig2EF_proximity | Fig. 2E–F | TE proximity mechanism data | [from expansion_analysis/results] |
| Fig2GH_sex_locus | Fig. 2G–H | Sex locus architecture and synteny data | species, clade, mat_locus_genes, synteny_block_id |
| Fig3A_dnmt | Fig. 3A | DNMT copy number per species + TE content | species, clade, DNMT1_n, DIM2_n, RID_DNMT2_n, total_te_pct |
| Fig3B_5mc_family_age | Fig. 3B | Weighted-mean 5mC per TE family × age bin | species, family, age_bin_kimura, wm_5mC, n_copies_covered |
| Fig3C_gene_dist_5mc | Fig. 3C | Gene–TE distance by methylation quartile | species, distance_bp, quartile, n_copies |
| Fig3D_adj_expr_meth | Fig. 3D | Adjacent gene expression by methylation quartile | species, gene_id, te_id, quartile, log2_rpkm |
| Fig3E_gene_dist_h3k9 | Fig. 3E | Gene–TE distance by H3K9me3 enrichment class | species, distance_bp, h3k9_class, n_copies |
| Fig3F_adj_expr_h3k9 | Fig. 3F | Adjacent gene expression by H3K9me3 peak | species, gene_id, te_id, h3k9_peak_present, log2_rpkm |
| Fig3G_mac_tss | Fig. 3G | TE-associated TSS expression by condition | species, condition, te_class, tss_count, expression_metric |
| Fig3H_metaplot | Fig. 3H | Per-bin 5mC across TE bodies (Q1 vs Q4) | superfamily, bin_index, quartile, mean_5mC, sem_5mC, n_copies |
| Fig3I_mac_tss_families | Fig. 3I | TE families enriched at methylated TSS | species, superfamily, n_tss_overlaps, enrichment_ratio, BH_q |
| Fig3J_rnai_pies | Fig. 3J | RNAi pathway targeting proportions | species, te_superfamily, rnai_status, proportion |
| Fig3K_srna | Fig. 3K | sRNA abundance vs. TE expression per superfamily | species, superfamily, srna_rpm, te_rpkm, srna_status |
| Fig4A_cooption_tree | Fig. 4A | Co-opted gene counts per species / clade | species, clade, n_coopted_genes, n_te_core, n_te_assoc |
| Fig4B_domain_barchart | Fig. 4B | Host-function domain composition of co-opted genes | clade, pfam_domain, domain_desc, n_coopted_genes, fraction |
| Fig4C_te_class_prop | Fig. 4C | TE class origin of co-opted genes per clade | clade, te_class, n_coopted, fraction |
| Fig4D_fusion_network | Fig. 4D | Fusion protein network edge list | te_domain_ipr, host_domain_ipr, n_proteins, n_species, species_list |
| Fig4E_chimte | Fig. 4E | ChimeraTE chimeric junctions per species | species, clade, te_class, n_chimeric_junctions |
| Fig4FG_fusion_trees | Fig. 4F–G | FLYWCH and CotH fusion protein tree metadata | protein_id, species, clade, te_superfamily, fusion_type |
| Fig4H_dam_rt_tree | Fig. 4H | DAM methyltransferase + RT protein metadata | protein_id, species, clade, te_superfamily, domain_arch |
| Fig4I_tf_enrichment | Fig. 4I | TF peak enrichment in TE classes (DAP-seq) | te_class, motif_id, motif_family, log2_odds_ratio, fisher_pval, BH_q |
| Fig4J_hmgb3_metaplot | Fig. 4J | HMGb3 signal metaplot over LINE/I-Jockey | bin_index, mean_signal, sem_signal, te_class |
| Fig4L_hmgb3_heatmap | Fig. 4L | HMGb3 occupancy across Jockey clades | clade, te_subfamily, mean_norm_signal |
| Fig5A_paralogs | Fig. 5A | Gene duplication rates and co-linear paralog counts | species, node_id, dup_rate, n_colinear_paralogs |
| Fig5B_syntenic_blocks | Fig. 5B | Syntenic block composition and TE content | species, block_id, block_length_bp, te_bp, te_fraction |
| Fig5C_te_breakpoints | Fig. 5C | TE enrichment at WGD breakpoints (log₂FC) | species, superfamily, log2_fc, wilcoxon_p, BH_q |
| Fig5D_te_dynamics | Fig. 5D | TE dynamics combined panel C data | species, condition, superfamily, metric_value |
| Fig5E_te_ks_corr | Fig. 5E | TE superfamily burst vs. WGD Ks correlation | species, superfamily, te_burst_ks_mode, anchor_ks_median, spearman_rho, spearman_p |
| Fig6A_htt_chord | Fig. 6A | HTT candidate counts per clade pair | clade_A, clade_B, n_candidates, n_validated |
| Fig6B_htt_ks | Fig. 6B | Ks distribution of 223 HTT candidates | candidate_id, te_class, te_superfamily, Ks_value |
| Fig6C_htt_heatmap | Fig. 6C | HTT candidate count heat-map per clade pair | clade_A, clade_B, n_candidates |
| Fig6D_htt_line_synteny | Fig. 6D | LINE HTT locus flanking conservation | event_id, species_A, species_B, flank_gene_A, flank_gene_B, conservation_status |

All numerical values are the direct inputs to the corresponding figure panels.
Processed intermediate tables are available in the per-section `results/`
directories of the code repository; paths and formats are defined in `config.sh`.
