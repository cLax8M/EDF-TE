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
5. **Relative composition variability** (RCV; PhyKIT v2.0) ≤ 0.5.  RCV measures whether
   the amino acid composition of the alignment is consistent across all sampled taxa.
   For each of the 20 amino acids, the mean frequency across all taxa in the alignment
   is computed; RCV is then the average absolute deviation of those frequencies from a
   perfectly uniform distribution (i.e. 1/20 = 5 % per amino acid), normalised by the
   alignment length: `RCV = mean(|f_i − 0.05|) / (4 × alignment_length)`.  A low RCV
   (< 0.5) indicates that amino acid frequencies are relatively consistent across
   sampled lineages, a prerequisite for reliable phylogenetic inference under standard
   substitution models; high RCV suggests that some taxa have strongly divergent amino
   acid composition, which can lead to long-branch attraction artefacts by grouping
   unrelated taxa with convergently biased sequences.  The threshold of 0.5 is standard
   in the field (Mullen *et al.*, 2023).

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
computed to seed the final analysis.  The definitive species tree was inferred under
the **LG+G4 model** (IQ-TREE 2.3.6, `-B 1000` UFBoot; 1,000 ultrafast bootstrap
replicates; gamma shape parameter α = 0.824565), which provides a well-calibrated
and computationally tractable fit for concatenated supermatrices of this size.  The
tree was rooted on the Holozoa clade using midpoint rooting followed by manual
inspection.

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
`GenomicFeatures/scripts/01_plot_genomic_features.py`).  Genome assembly size (Mb) vs. total
TE content (%) was plotted per species, coloured by clade (Fig. S3A;
`GenomicFeatures/correlations/scatter_genome_size_vs_TE.pdf`).  Pairwise Pearson correlations
between TE content metrics and genomic features across 238 species were visualised as a
heatmap (Fig. S3B; `GenomicFeatures/correlations/heatmap_te_vs_features_raw_pearson.pdf`;
`GenomicFeatures/scripts/03_te_genomic_correlations.py` combined with
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
inversions; Putnam *et al.*, 2007).  The count of such collinear pairs is noted
as `observed_collinear_pairs`; the total number of adjacent-ortholog slots across
all scaffolds of species i is `total_adjacent_slots`.

#### 4.4 Null distribution and spurious synteny correction

To account for spurious gene-order conservation from random co-localisation, random
genomic positions were assigned to each gene using `bedtools shuffle` v2.24 (Quinlan &
Hall, 2010) in 100 independent rounds per species, with scaffolds preserved as location
constraints (`-chrom` flag).  Each replicate used a fixed per-replicate random seed.
The mean over 100 replicates (the null expectation, `null_mean_collinear`) was computed for each species pair.

#### 4.5 Patristic-distance-corrected synteny ratio

Even after subtracting the randomisation baseline (§4.4), raw pairwise synteny scores
still reflect evolutionary divergence: species pairs that diverged longer ago share fewer
collinear gene arrangements simply from neutral drift, independently of any TE-driven
rearrangement.  To remove this divergence-time confound, synteny scores were corrected
against patristic distances extracted from the species guide tree (Section 2), following
a three-step procedure implemented in `Synteny/scripts/09_patristic_correction.py`.

**Step 1 — raw score per pair.**  For each pairwise comparison (i, j), the baseline
null-corrected score is:

`raw_score = (observed_collinear_pairs − null_mean_collinear) / total_adjacent_slots`

This is the number of collinear ortholog pairs in the real genome minus the mean number
recovered from 100 randomised gene-position replicates, divided by the total number of
consecutive ortholog-pair slots tested in that comparison.  Subtracting the null mean
removes the random co-localisation baseline; dividing by total slots normalises for
assembly fragmentation (highly fragmented assemblies expose fewer adjacent-gene slots
and would otherwise score artificially low).  Negative values — which arise by chance
in very fragmented assemblies — were set to zero.

**Step 2 — patristic distance and polynomial regression.**  The patristic distance
d_ij between species i and j was computed as the sum of branch lengths along the path
connecting the two tips in the rooted guide tree (Bio.Phylo, Cock et al., 2009).  Within
each clade, a second-degree polynomial was fitted to the scatter of raw_score against
d_ij across all species pairs in that clade:

`expected_raw(d_ij) = a × d_ij² + b × d_ij + c`

The polynomial captures the empirical decay of synteny with divergence time within the
clade.  The corrected score is then the residual — how much more or less synteny the
pair retains than expected at its particular divergence level:

`corrected_score = raw_score − expected_raw(d_ij)`

A positive residual means the pair retains more collinear gene order than typical for
two species that diverged as long ago as d_ij suggests; a negative residual means gene
order has eroded faster than expected, consistent with increased rearrangement or TE
activity.

**Step 3 — within-clade normalisation.**  Residuals were linearly rescaled to [0, 1]
within each clade:

`normalised_score = (corrected_score − min_clade) / (max_clade − min_clade)`

This makes scores comparable across species within a clade while acknowledging that
clades with very different divergence ranges cannot be directly compared on the same
axis.  A normalised score of 1.0 corresponds to the most synteny-conserved pair in the
clade after accounting for divergence; a value near 0 indicates retention close to or
below the divergence-adjusted expectation.

A per-species ranking by patristic-distance-corrected synteny ratio is shown in
Fig. S1B.  Per-clade boxplots of corrected synteny ratios are shown in Fig. S1C
(`Synteny_guide_corrected/figures/synteny_ratio_by_clade_box.pdf`).  Column names in
the source data table (Fig1D_synteny) reflect the four quantities computed above:
`observed_collinear_pairs`, `null_mean_collinear`, `total_adjacent_slots`,
`corrected_score`, and `normalised_score`.

#### 4.6 Circos visualisation of selected pairs

For four species pairs with the highest synteny ratios, syntenic gene-pair links were
visualised on the top-30 scaffolds per species using Circos v0.67 (Krzywinski *et al.*,
2009).  The four Circos panels are shown in Fig. 1D
(`Synteny_guide_corrected/figures_circos_v6_4pairs/`).

---

### 5. TE annotation with EarlGrey

Each of the 320 genomes was annotated for transposable elements using EarlGrey v1.6.
EarlGrey orchestrates a modular pipeline combining de novo repeat discovery,
structural filtering, database-guided classification, and genome masking; full details
of each step are given in the subsections below.

#### 5.1 De novo consensus library construction

For each genome, putative TE consensus sequences were generated using the
`earlGreyLibConstruct` module of EarlGrey v1.6, configured with Dfam 3.7 curated
elements as the reference database, using default settings.  Briefly, RepeatModeler2
v2.0.5 was used to identify repetitive elements across each genome assembly, combining
the RECON + RepeatScout pipeline for tandem-masked sequences with structural LTR
discovery via LTRharvest (GenomeTools v1.6.2) and LTR_retriever v2.9.0.

All per-species consensus sequences were combined into a single multi-FASTA file.
To generate a clustered non-redundant library, we employed a scalable cascaded
clustering approach using MMseqs2 v14.7564 `easy-cluster` with
`--min-seq-id 0.8 -c 0.8 --cov-mode 1 --cluster-reassign`, retaining the
representative sequence for each cluster.

#### 5.2 Structural and quality filtering

After consensus generation, EarlGrey applied the following sequential filters:
- Consensi shorter than 80 bp or longer than 40,000 bp were discarded.
- For putative DNA transposons: terminal inverted repeats (TIRs) were verified with
  IRScan (EarlGrey internal); consensi lacking TIRs of length ≥ 10 bp with identity
  ≥ 70 % were reclassified as Unknown.
- For putative LTR retrotransposons: LTR boundary integrity was verified by
  LTR_retriever v2.9.0; consensi failing the structural filter were reclassified.
- Consensi with > 50 % simple or low-complexity sequence were removed.

#### 5.3 Host-gene filtering

To identify and remove putative host genes that may have been captured in the TE
library, open reading frames were identified in all consensus sequences and matched
against the RefSeq non-redundant protein database (NCBI NR, release 221, December
2023) using DIAMOND BLASTX v2.1.8 (`--more-sensitive -e 1e-5 --id 40 -k 1`).
Sequences with a RefSeq protein hit (identity ≥ 40 %, e-value ≤ 10⁻⁵) that did not
also match known TE proteins (i.e. the best RefSeq hit annotation did not contain
the terms *transposon*, *transposase*, *integrase*, *reverse transcriptase*, *capsid*,
or *gag*) were classified as host-derived and removed from the TE library.

#### 5.4 Phylogenetic classification — Dfam nhmmer search

Each consensus was classified against the Dfam 3.7 database (20,462 HMM profiles)
using nhmmer (HMMER 3.3.2):

```
nhmmer --cpu 8 -E 1e-5 --tblout <out.tbl> Dfam3.7.hmm <consensi.fa>
```

The best-scoring Dfam hit determined the superfamily-level classification.  Consensi
with no Dfam hit (E > 10⁻⁵) were passed to RepeatClassifier (RepeatModeler2 v2.0.5)
for automated kmer-based assignment.  Consensi that remained unclassified after both
steps were labelled `Unknown`.  TE superfamily nomenclature follows Wicker et al.
(2007).

#### 5.5 Genome annotation with RepeatMasker

The curated per-species consensus library was used to annotate the repeat content of
the corresponding genome assembly with RepeatMasker 4.1.5:

```
RepeatMasker -lib <species_library.fa> -pa 16 -s -dir <out_dir> \
  -gff -inv -a -nolow -norna <genome.fa>
```

`-s` applies the slow/sensitive search mode; `-a` produces the alignment (`.align`)
output file required for Kimura divergence computation (Section 6); `-nolow` excludes
simple repeats and low-complexity sequences from masking; `-norna` excludes small RNA
genes.  The output `.filteredRepeats.bed` file records one BED line per annotated TE
copy (chromosome, start, end, family, score, strand).

All EarlGrey SLURM jobs were submitted to partition `normal.168h` with 16 CPUs and
64 GB RAM per species.  Two Microsporidia genomes (MI_Encro1, MI_Ordco1) returned
zero TE content after filtering and were excluded from ordination analyses but
retained in all other comparative analyses.

---

### 6. Kimura divergence landscapes and the Recency Index

#### 6.1 K2P distance computation

RepeatMasker produces a per-copy pairwise alignment in the `.align` output file.
Each alignment block records, for the region of the genome occupied by a given TE copy,
how many nucleotide positions match the consensus sequence, how many differ by a
transition (purine–purine or pyrimidine–pyrimidine substitution: A↔G or C↔T), and how
many differ by a transversion (purine–pyrimidine substitution: A/G ↔ C/T).  These
counts are used to compute the Kimura two-parameter (K2P) divergence, which corrects for
the fact that (i) the same position may have been substituted more than once since the TE
inserted (multiple-hit correction), and (ii) transitions occur at a higher rate than
transversions and would otherwise be under-counted.  The K2P formula is:

`K = −0.5 × ln(1 − 2p − q) − 0.25 × ln(1 − 2q)`

where p is the proportion of sites differing by a transition (transitions / aligned sites)
and q is the proportion of sites differing by a transversion (transversions / aligned
sites).  Aligned sites are counted after trimming terminal gaps from both the copy and the
consensus.  Because K increases monotonically with divergence time since insertion, it
serves as a molecular clock proxy for TE insertion age: a copy with K = 0 was inserted
very recently (zero divergence from its consensus), whereas a copy with K > 0.3 has
accumulated so many substitutions that its ancestral boundaries are difficult to recover
reliably.

CpG-context transitions (C→T or G→A at CpG dinucleotides) accumulate approximately
ten-fold faster than non-CpG transitions and would artificially inflate the divergence of
CpG-rich elements; these were corrected using the CpG-adjusted substitution counts
reported separately by RepeatMasker in the `.align` file.  Copies with fewer than 10
aligned sites after gap-stripping were discarded as unreliable.

(`section2_te_landscape/04_kimura_divergence.py`)

#### 6.2 Divergence landscape histograms

Kimura divergence landscapes were constructed as stacked histograms per species.
Each histogram spans K values from 0 to 1, divided into bins of width 0.005 (200 bins
total).  For each bin, the y-axis shows the total genomic space occupied by TE copies
whose K2P distance falls within that bin, expressed in base pairs.  The genomic space
occupied by a copy (hereafter "copy length") is simply the length of the annotated
region on the chromosome in base pairs, taken directly from the start and end coordinates
in the `.filteredRepeats.bed` output file.  Stacking by TE class follows the canonical
Wicker et al. (2007) class hierarchy: Class I LTR, Class I non-LTR, Class II DNA
transposons, Helitrons, and Unknown.

#### 6.3 Recency Index

The Recency Index (RI) for a given species is the proportion of the total TE content
(measured in base pairs) that consists of recently inserted copies.  Specifically, for
each species we summed the lengths of all annotated TE copies with K2P distance < 0.10,
and divided that sum by the total length of all annotated TE copies in the genome:

`RI = (total bp from copies with K < 0.10) / (total TE bp in the genome)`

A copy’s length is defined as the number of base pairs it occupies on the chromosome,
read from the start and end coordinates in the RepeatMasker annotation file.  The RI
therefore ranges from 0 (no recently inserted copies) to 1 (all TE content is recent).

The threshold K < 0.10 was chosen empirically to capture copies that have accumulated
fewer than ~10 % of the substitutions expected under neutral decay.  At typical fungal
synonymous substitution rates (~2 × 10⁻⁹ per site per year), a K value of 0.10
corresponds roughly to insertion events occurring within the last ~50 Myr of the
analysed divergence time window.

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

To measure how diverse the TE repertoire of each species is, we treated each TE
superfamily as analogous to a biological species in a community ecology framework,
and the number of base pairs contributed by that superfamily to the genome as its
"abundance".  Structural RNA and simple repeat classes were excluded before computation
(Simple_repeat, Low_complexity, Satellite, rRNA, tRNA, snRNA, scRNA, sno, srpRNA).
All values were computed from the EarlGrey `familyLevelCount.txt` output file.

The **Shannon diversity index** (H′, bits) measures the overall complexity of the TE
community.  It is calculated as: for each superfamily, multiply its proportional
contribution to total TE base pairs by the log₂ of that proportion, sum these products
across all superfamilies, and negate the result: `H′ = −Σ p_j × log₂(p_j)`.  A value
of H′ = 0 bits means the entire TE content is dominated by a single superfamily;
higher values indicate a more diverse community.  The theoretical maximum,
log₂(number of superfamilies), is reached when all superfamilies contribute equally.

Because H′ increases automatically with the number of detected superfamilies, we also
computed **Pielou’s evenness** (J), which rescales H′ to the range [0, 1] by dividing
by this theoretical maximum: `J = H′ / log₂(number of superfamilies)`.  A value of
J = 1 means all superfamilies contribute equally to the TE landscape; a value near 0
means a handful of superfamilies account for virtually all the TE content.  A high H′
combined with a low J indicates a diverse but uneven community — many superfamilies
present but a handful responsible for most of the genomic mass — whereas high H′ and
high J indicate a genuinely diverse and balanced TE repertoire
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

DNMT1-type, DIM2-type, and RID/DNMT2-type DNA methyltransferases were identified
by a multi-step approach combining domain screening and phylogenetic classification.
First, all proteins containing the C5-cytosine methyltransferase catalytic domain
(PF00145) were retrieved from InterProScan 5.67-99.0 (Pfam 36.0) annotations of
all 320 predicted proteomes.  These candidate proteins were then searched by
reciprocal BLASTP (BLAST+ 2.12.0; E ≤ 10⁻⁵) against a curated set of reference
DNMT1, DIM-2, and RID/DNMT2 sequences from characterised fungal species
(`TE_dynamics/dnmt_machinery/scripts/31_extract_all_C5MT.py`).  The resulting
hits were aligned with MAFFT v7.520 (`--localpair --maxiterate 1000`), trimmed
with trimAl v1.4 (`-gappyout`), and used to build a maximum-likelihood phylogeny
(IQ-TREE2, LG+G4 model, 1,000 UFBoot).  Orthology was assigned based on
monophyletic clustering with the reference sequences.  For each protein, the
presence of additional diagnostic domains (RFTS, BAH, DNAmet_Nter, Helicase_C,
SNF2_N, PHD, Chromo) was recorded.  Per-species copy numbers and domain
architectures were compiled in `dnmt_matrix_final.tsv`
(`TE_dynamics/dnmt_machinery/data/`).  The DNMT1 vs. DIM2 vs. RID classification
followed Iyer *et al.* (2011) domain combination rules.

The species tree was annotated with per-species DNMT copy numbers and TE content
(% genome) from `dnmt_matrix_final.tsv` and `te_content_wide.tsv` (Fig. 3A;
`TE_dynamics/dnmt_machinery/figures/dnmt_te_figure_rotated.pdf`).  The scatter of
DNMT1/DIM2 copy number vs. total TE content is shown for Mucoromycota species in
Fig. S5A (`TE_dynamics/dnmt_machinery/dnmt1dim2_vs_te_content_MU.pdf`).

---

### 9. 5-methylcytosine (5mC) analysis

Whole-genome bisulphite sequencing was performed for *Phycomyces blakesleeanus*.

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

For each TE copy, the methylation level was computed from covered CpG positions only.
For each CpG position within the copy with sequencing coverage ≥ 5×, the number of
methylated calls (C) and unmethylated calls (U) was read from the Bismark cytosine
report.  The per-copy methylation fraction is simply:

`m = sum of methylated calls / sum of (methylated + unmethylated calls)`

where the sum runs only over positions with coverage ≥ 5×.  Copies with fewer than
3 covered CpG positions were excluded.  A copy was classified as methylated if
m ≥ 0.50.

#### 9.3 5mC landscape by TE family and Kimura age

5mC levels per TE family stratified by Kimura divergence age class were computed as a
coverage-weighted mean per family per age bin, extending the approach of Chiapello *et al.*
to additional Mucoromycota species
(`TE_dynamics/20_chiapello_extension/`; Fig. 3B = `F4_wm_by_family_age.pdf`).  Spearman's
rank correlation between the per-copy methylation fraction (m) and K2P divergence (K) was computed per species.  A separate
analysis for dark-mycelium conditions is shown in Fig. S5B
(`TE_dynamics/12_te_age_5mc/plots/01_wm_by_age_class.pdf`).

#### 9.4 Gene proximity and expression effects

The distance from each gene to the nearest full-length TE insertion was computed and
stratified by methylation quartile (Q1–Q4 of the per-copy methylation fraction m; Fig. 3C:
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

H3K9me3 ChIP-seq was performed for *Rhizopus microsporus* and *Mucor lusitanicus*.

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

H3K9me3 enrichment per TE copy was computed as the mean bedgraph signal over the
copy coordinates, normalised to RPKM.  RPKM (reads per kilobase per million mapped
reads) standardises read depth for both library size and feature length:

`RPKM = (reads overlapping copy × 10⁹) / (copy length in bp × total uniquely mapped reads)`

This allows direct comparison of enrichment levels across copies of different
lengths and between samples with different sequencing depths.

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
in Fig. S8 (`TE_dynamics/48_adjacent_genes/48_family_summary_scatter_v2.pdf`).

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

`Fuzzy_high` serves as a heterochromatin-like proxy: in *Rhizopus microsporus* and
*Mucor lusitanicus* (where in vivo H3K9me3 data is available), H3K9me3-enriched
regions are characterised by high mean occupancy and low phasing — consistent with
the `Fuzzy_high` state definition.

#### 11.3 TE occupancy relative to background

For each species, the mean predicted nucleosome occupancy of TE-overlapping 500 bp
windows was compared to background (non-TE, non-genic) windows.  Delta occupancy
(Δocc = TE − background) per TE superfamily was summarised per species
(`Section2_TElandscape/nucleosomes/f3_bias/{SP}_f3_bias.tsv`; 205 species with
columns: species, superfamily, mean_te_occupancy, background_occupancy, delta_occupancy,
n_insertions).  Cross-species and cross-clade summaries are shown in Fig. S4B
(`Section2_TElandscape/supervised_chromatin/invivo_landscape/figS_nuclpred_chromatin_withpoints.pdf`).

#### 11.4 In vivo nucleosome and chromatin accessibility analysis

To validate the computational nuclPred predictions with native chromatin data, in
vivo chromatin accessibility and nucleosome positioning were measured in four
species using three experimental approaches.

**MNase-seq (***Rhizopus microsporus* **ATCC11559, legacy dataset).**  A
previously generated MNase-seq WIG track for *R. microsporus* (MU_Rhimi59) was
used directly as input; raw reads were not reprocessed.  This dataset provides
the single-nucleosome resolution occupancy signal used as an in vivo reference
throughout the study.

**MNase-seq and MAINE-seq (***Neurospora crassa*** OR74A and ***Botrytis cinerea***
B05.10).**  Raw reads were downloaded from the ENA (N. crassa SRR13067301,
SRR13067302, SRR13067303; B. cinerea SRR13224644, SRR13224647, SRR13224648).
For each species, replicate FASTQ files were quality-trimmed with BBDuk v38.87:

```
bbduk.sh in1=<R1.fq.gz> in2=<R2.fq.gz> \
  out1=<trimmed_R1.fq.gz> out2=<trimmed_R2.fq.gz> \
  ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 minlen=50
```

Trimmed reads were aligned to the respective reference genomes with Bowtie v1.3.1
(single-end mode, using the mono-nucleosomal insert range):

```
bowtie -v 2 -m 1 --best --strata -p 16 \
  -x <genome_index> -q <trimmed.fq> -S <mapped.sam>
```

`-m 1` retains only uniquely mapping reads to avoid multi-mapping artefacts in
repetitive regions.  BAMs were coordinate-sorted with samtools v1.17 and PCR
duplicates were removed with Picard MarkDuplicates v3.1 (`REMOVE_DUPLICATES=true`).
Per-replicate duplicate-removed BAMs were merged into a single pooled BAM per
species with `samtools merge`.

DANPOS v3 (`dpos` subcommand) was applied to the pooled BAM to call
per-nucleosome dyad positions and compute per-nucleosome metrics:

```
python danpos.py dpos <pooled_dedup.bam> -m 1 -a 1 -o <danpos_dir>
```

`-m 1` enables background model subtraction; `-a 1` enables DANPOS
normalisation.  The output `positions.xls` file provides, for every called
nucleosome, the smoothed dyad position (`smt_pos`), peak height (`smt_value`,
proportional to nucleosome sharpness/occupancy), and fuzziness score
(standard deviation of the read-end distribution around the dyad, in bp;
lower values indicate better-positioned nucleosomes).  A smoothed WIG occupancy
file (`*.smooth.wig`) was also produced for regional analyses.

**ATAC-seq (***Microbotryum eburneum***).**  Three paired-end ATAC-seq replicates
(SRR28069093, SRR28069094, SRR28069106) were downloaded from the SRA.  Reads were
quality-trimmed with fastp v0.23.4 (`--detect_adapter_for_pe`), aligned to the
*M. eburneum* genome with Bowtie2 v2.5.2 (`--no-mixed --no-discordant --maxins 2000`),
and duplicates were removed with samtools markdup.  A per-base accessibility WIG
was computed with DANPOS `dpos` (same parameters as above); note that ATAC-seq
signal represents chromatin accessibility (high signal = open chromatin), which
is directionally opposite to MNase occupancy.

**Per-nucleosome feature analysis (***N. crassa*** and ***B. cinerea*** only).**
Three features were extracted from the DANPOS `positions.xls` for each called
nucleosome and compared across three genomic contexts — TE loci (dyad inside an
EarlGrey-annotated TE interval), protein-coding gene bodies (dyad inside a GFF
`gene` feature), and intergenic regions (neither):
(i) `fuzziness_score` — positional SD of the read population (lower = more precisely
positioned);
(ii) `smt_value` — smoothed peak height (higher = more nucleosome-like, better
protected);
(iii) `nrl_bp` — distance to the next dyad on the same chromosome (nucleosome
repeat length proxy; filtered to 100–400 bp).
Group comparisons used two-sided Mann–Whitney U tests (BH-corrected FDR 5 %).

**Regional TE-density analysis (all four in vivo datasets).**  The genome was
partitioned into 5 kb and 10 kb non-overlapping windows.  Windows were ranked
into four TE-density quartiles (Q1 = TE-poorest, Q4 = TE-richest) based on
fractional overlap with EarlGrey annotations.  Four features were computed from
the smoothed WIG per window: mean signal (`mean_sig`), coefficient of variation
(`cv`), autocorrelation peak value in the 140–240 bp lag range (`ac_val`, nucleosome
phasing strength), and the lag of that peak (`ac_lag`, nucleosome repeat length).
For ATAC-seq data, `mean_sig` was sign-inverted before comparison with MNase
occupancy to align the directionality.  Quartile comparisons used Kruskal–Wallis
tests followed by pairwise Mann–Whitney U tests (BH-corrected FDR 5 %)
(`Section2_TElandscape/scripts/invivo_nucleosome_landscape.py`,
`mnase_invivo_nucleosomes_at_te.py`, `danpos_fuzziness_at_te.py`).

The combined in vivo dataset (MNase, MAINE and ATAC) is shown in Fig. S4A
(`Section2_TElandscape/supervised_chromatin/invivo_landscape/figF2_invivo_multimetric_withpoints.pdf`).

---

### 12. N6-methyladenine detection and small RNA silencing

N6-methyladenine (6mA) PacBio kinetics data were generated for *Phycomyces
blakesleeanus* and *Rhizopus microsporus*.  Small RNA sequencing was performed
for *Mucor lusitanicus*.

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
coverage ≥ 25× AND the ipdSummary p-value ≤ 0.001.  Methylated adenine clusters
(MACs) were defined as dense runs of 6mA-positive positions within TE-associated
regions; MACs overlapping the TE transcription start site window (−150 to +400 bp)
were used for the TSS analysis (Figs. 3G, 3I).  The nearest-gene RNA-seq RPKM (from
Section 10.4) was compared between 6mA-positive and 6mA-negative TE copies using a
two-sided Wilcoxon rank-sum test (`TE_dynamics/35_srna_mucor/` and
`TE_dynamics/45_mac_tss_expression/`).

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
--strata` ensures the alignment reported is among the best-scoring ones.  Reads per
million (RPM) per TE family were computed as: the number of reads mapped to that
family divided by the total number of mapped reads (including multimappers), multiplied
by 10⁶.

`RPM = (reads mapped to family / total mapped reads) × 1,000,000`

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

#### 13.2 ChimeraTE chimeric transcript detection

Chimeric gene–TE transcripts were detected with ChimeraTE v1.2 (Mode 1, genome-guided
RNA-based mode; Oliveira *et al.*, 2022; https://github.com/OliveiraDS-hub/ChimeraTE)
across 24 species spanning 15 clades
(`Cooption/scripts/chimte_bedtools_perrep.py`, `chimte_merge.py`,
`prep_chimeraTE_data.py`).

**RNA-seq datasets.**  Three species were analysed using previously published
RNA-seq data from the MycoMobilome project (MU_Mucci3, MU_Rhimi59, MU_Phybl4;
read format: interleaved FASTQ, PE).  For the remaining 21 EDF species, publicly
available RNA-seq runs were downloaded from the SRA and processed in paired-end
mode (species listed in Source Data; see §17 for accession IDs):
AP_Amoprot, BL_Blabri1, BS_Basme92, CH_Batde5, CH_Batsa1, CH_Pirfi3,
EN_Entmai1, EN_Zoorad1, GL_Geopyr1, GL_Rhicl1, GL_RhiirA4,
KX_Coemoj1, KX_Kicala1, MI_Enccu1, MI_Sprlo1, MI_Vavcu1,
MO_Disdec1, MO_Morel2, NU_Fontalba, ZO_Pipcy3.

**Alignment.**  Reads were aligned to each species' reference genome (repeat-soft-masked)
with STAR v2.7.11b in 2-pass mode:

```
STAR --runMode alignReads \
  --genomeDir <STAR_index> \
  --readFilesIn <R1.fastq.gz> <R2.fastq.gz> \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 100 \
  --winAnchorMultimapNmax 200 \
  --twopassMode Basic \
  --runThreadN 16
```

`--outFilterMultimapNmax 100` and `--winAnchorMultimapNmax 200` allow reads to
multi-map to repetitive TE copies, which is necessary to detect chimeric transcripts
spanning gene–TE junctions; the more permissive alignment policy is standard for
ChimeraTE.

**Chimeric junction detection.**  Per-replicate chimeric read counts were computed
with a BEDtools-based reimplementation (`chimte_bedtools_perrep.py`) that reproduces
the ChimeraTE Mode 1 logic: reads spanning the boundary between a gene body and an
overlapping or adjacent TE copy are classified by junction type (TE-initiated,
TE-terminated, or TE-internal).  Results were merged across replicates per species
(`chimte_merge.py`).  TE classes excluded from ChimeraTE analysis: `Simple_repeat`,
`Satellite`, `Low_complexity`, `rRNA`, `snRNA`, `srpRNA`, `tRNA`, `scRNA`, `RNA`,
`ARTEFACT`.  A chimeric event was reported if ≥ 2 independent reads supported the
junction in ≥ 1 replicate.

#### 13.3 TF binding site enrichment in TE sequences

DAP-seq peak data for **58 transcription factors** from *Rhizopus microsporus*
ATCC11559 were analysed in two complementary ways
(`TF_TE_analysis/scripts/01_peak_TE_overlap.py` and
`TF_TE_analysis/scripts/03_fimo_motif_enrichment.py`).

**Peak–TE overlap (script 01).**  DAP-seq peak BED files (filtered at
q-value ≥ 10⁻²) were intersected with the EarlGrey TE annotation of *R. microsporus*
ATCC11559 using `bedtools intersect`.  For each TF, the fraction of peaks overlapping
any TE copy was computed; enrichment relative to the genomic TE coverage was tested
with a two-sided Fisher's exact test (BH-FDR 5 %).

**FIMO motif enrichment in TE sequences (script 03).**  FIMO v5.5.5 was run for
each of the 58 TF motif files (MEME format; `data/TFsDAP/motifs/`) against the TE
sequences of each TE class, and separately against 3× size-matched genomic shuffled
sequences generated with `bedtools shuffle`.  The p-value threshold for FIMO hit
reporting was 1 × 10⁻⁴.  Enrichment in TE sequences relative to the shuffled
background was tested with a two-sided Fisher's exact test on the 2 × 2 contingency
table (TE sequences with ≥ 1 FIMO hit vs. background sequences with ≥ 1 FIMO hit);
p-values were corrected across all TF × TE class combinations using the
Benjamini–Hochberg procedure (FDR 5 %).

The strongest enrichment identified in the *R. microsporus* analysis was for the
HMGb3 transcription factor (DAP-seq motif HMGb_1228) in LINE/I-Jockey elements
(100 % of Jockey copies contained at least one FIMO hit; OR = 18.3, BH-q < 10⁻⁴).
To determine whether this association is phylogenetically widespread, the FIMO
cross-species analysis was extended to all 236 EDF species with EarlGrey TE
annotations and available HMGb3 motif data
(`TF_TE_analysis/scripts/run_fimo_cross_species.sh`; SLURM array, one job per species).
The fraction of LINE/I-Jockey copies carrying an HMGb3 binding site was computed per
species and plotted sorted by clade (`TF_TE_analysis/scripts/plot_hmgb3_cross_species.py`;
`plot_hmgb3_clade_heatmap.py`).  Species without annotated I-Jockey elements were
excluded from the cross-species panel.

**Gypsy LTR × bHLH_Myc motif association.**  A second focused enrichment was the
association between LTR/Gypsy elements and helix–loop–helix transcription factors
(HLH_Myc class; motifs HLH_Myc14, HLH_Myc5, HLH_Myc15).  To determine whether the
HLH_Myc binding affinity co-segregates with a specific Gypsy clade, the RT domain
sequences of all 103 Gypsy elements in *R. microsporus* ATCC11559 carrying an
identifiable RT-domain were extracted and aligned.  A maximum-likelihood tree was
inferred with IQ-TREE2 (`-m MFP -bb 1000 --redo`).  Gypsy elements were classified
into RT-domain clades using TEsorter v1.4.7 with the REXdb/Metazoa database
(`--db REXdb -p 4`), which resolved four clades: Chlamyvir, Reina, chromo-outgroup,
and non-chromo-outgroup.  Each leaf of the RT tree was annotated with: (i) its TEsorter
clade assignment; (ii) presence of DAP-seq peak overlap for each of the three bHLH_Myc
TFs; and (iii) the number of FIMO hits per HLH_Myc motif
(`TF_TE_analysis/scripts/plot_gypsy_tree_hlh_heatmap.py`).  The Chlamyvir clade was
specifically enriched for HLH_Myc binding sites (Fisher's exact, BH-FDR 5 %).

**GremLINE1 retrotransposon distribution.**  GremLINE1 is a LINE/I clade element
originally characterised as a centromeric repeat in Mucoromycota (Navarro-Mendoza
et al. 2019).  To determine its distribution across the EDF dataset, the published
GremLINE1 consensus sequence was used as a BLASTN query against the EarlGrey
consensus TE library of each Mucoromycota and Blastocladiomycota species
(`TF_TE_analysis/scripts/05b_identify_gremlins.sh`; BLAST parameters: ≥ 30 %
identity, ≥ 50 % query coverage, E ≤ 10⁻⁵, following Navarro-Mendoza et al. 2019).
For each species, the TE family names of all libraries with BLAST hits were
extracted (`gremlin_families_per_species.tsv`).  Downstream FIMO motif enrichment
and DAP-seq overlap analyses were performed on the identified GremLINE1-related
families using the same HMGb3 motif pipeline described above
(`TF_TE_analysis/scripts/05c_identify_gremlins_tblastn.sh` cross-validates with
tBLASTn on translated copies; `05d_identify_gremlins_genome.sh` scans genome copies
directly for flanking gene conservation).

#### 13.4 TE–host domain fusion proteins

Proteins carrying at least one TE-derived Pfam domain AND at least one host-function
domain in the **same polypeptide** were identified across all 238 species using the
InterProScan output (IPR.tab files; one per species, generated with InterProScan 5.65
as described in §13.1).  The screen was implemented in
`Cooption/scripts/02_classify_proteins.py` (architecture class `TE_host_fusion`) and
`Cooption/scripts/13_duf_te_fusions.py`.

**ORF-level fusion detection.**  For each protein, all Pfam domain accessions were
collected from the IPR.tab file (stripping version suffixes).  A protein was classified
as a TE–host fusion if its Pfam set contained at least one member of the TE-machinery
domain set (e.g. PF00078 RVT_1, PF00665 rve integrase, PF01609 DDE transposase,
PF10551 MULE, PF02994 L1 transposase — see full list in §13.1) AND at least one
domain not belonging to the TE domain set.  A symmetric check was performed for the
fusion screen (`Cooption/scripts/02_fusion_screen.py`): for each candidate expansion
domain (target Pfam), every protein carrying that domain was tested for co-occurrence
with any TE-domain in the same IPR.tab entry.  This approach provides direct protein-level
evidence (single ORF) that the host functional domain has been incorporated into a TE
coding frame or vice versa, independently of co-expression or proximity.

**DUF × TE enrichment analysis.**  Domain of Unknown Function (DUF) accessions
co-occurring with TE domains in the same protein were compiled across all species.
For each DUF × TE domain pair, enrichment relative to the null expectation (marginal
DUF frequency × marginal TE frequency across all proteins) was tested with Fisher's
exact test (BH-FDR 5 %).  The most prevalent co-occurring pairs include DUF4217
(PF13959) with Helitron/Helicase TE domains (156 species) and DUF1605 (PF07717)
with Helitron elements (44 species).

**Highlighted fusion types.**

*(i) FLYWCH domain (PF04500) × TE RT domains* — the FLYWCH domain is a zinc-finger
motif present in C2H2 transcription factors.  FLYWCH × reverse-transcriptase fusions
were identified in **Neocallimastigomycota** (clade NE) where copy numbers of
FLYWCH-domain proteins are markedly elevated relative to all other clades
(Mann–Whitney, BH-q < 0.05, fold ≥ 5; `expansion_analysis/01_copynumber_screen.py`).
Protein-level fusions carrying both PF04500 and a TE RT domain (PF00078/PF07727) in
a single ORF were detected in multiple NE species, suggesting domestication of TE
reverse-transcriptase activity into a chimeric FLYWCH effector.

*(ii) CotH kinase (PF08757) × TE domains* — CotH is the eponymous kinase of
coat-protein biosynthesis in Mucorales.  CotH genes are expanded in the CH clade
(Chytridiomycota) and show protein-level fusions with TE integrase (PF00665) and
transposase (PF01609) domains in multiple species.  The CotH expansion was the first
domain identified in the original seven-domain targeted screen
(`expansion_analysis/00b_discover_expansions.py`).

*(iii) NACHT NTPase domain (PF05729) × TE RT/integrase domains* — NACHT domains
are central to innate-immunity NLR proteins in animals.  In fungi, NACHT × TE
fusions were found scattered across multiple clades, with the TE partner being
predominantly LINE-class RT domains, suggesting the NACHT signalling cassette may
recruit TE-derived reverse transcriptase or endonuclease activity as a sensing
component.

*(iv) DAM N6-adenine methyltransferase (PF05869/IPR003356) fused within DIRS LTR
retrotransposons* — described in §13.5 below.  Briefly, 14 out of 15 confirmed
DAM-TE fusions involve DIRS-class elements, and the DAM domain is found within the
TE ORF rather than as an independently transcribed host gene, indicating that the
DIRS element has captured and carries the methyltransferase as a functional passenger.

**Phylogenetic trees for fusion domain families.**  For each of the 31 candidate
domains entering the master expansion table, phylogenetic trees were built to
visualise the distribution of TE-fused proteins relative to their non-fused
paralogues (`expansion_analysis/05_phylo_pipeline.sh`, SLURM array of 31 jobs).
For each domain: (1) all proteins carrying the target Pfam domain were extracted
from the 238 proteomes; (2) sequences were aligned with MAFFT v7.526 (`--auto`);
(3) the alignment was trimmed with trimAl v1.5 (`-automated1`); (4) a maximum-
likelihood tree was inferred with FastTree 2.2 (`-lg -gamma`).  Each leaf was
annotated with the clade identity of its host species and flagged if the protein
also carried any TE-machinery domain (from the TE-core domain set; see §13.1).
TE-fused proteins are displayed as large coloured symbols in the circular ggtree
visualisation (`fig_te_fusion_trees.R`), colour-coded by associated TE domain
type (RVT_1: blue; rve integrase: pink; Chromo: green; MULE transposase:
orange-red; Tc5 transposase: yellow-orange).  Trees for the highlighted domains
(FLYWCH, CotH, NACHT) are provided in Supplementary Figures.

For the Gypsy/DIRS RT domain context, a broader reverse-transcriptase phylogeny
was built using all proteins carrying both a chromodomain (PF00385 Chromo) and an
RT domain (PF00078 RVT_1 or PF07727 RVT_2) across all 238 proteomes
(`Cooption/scripts/fig_H_gypsy_chromo_phylotree.py`).  Sequences were aligned
with MAFFT v7.526 (`--auto`) and a maximum-likelihood tree was inferred with
IQ-TREE2 (`-m MFP -bb 1000`), allowing model selection.  This tree was used to
classify Gypsy elements and is shown in Fig. 4.

For the DAM + RT phylogeny specifically, the RT domain of each DAM-containing ORF
was extracted by profile-guided alignment (`hmmalign --trim`) to the PF00078 HMM;
the DAM domain was aligned separately with MAFFT after trimming to the PF05869
profile.  Both trees were inferred with FastTree 2.2 (`-lg -gamma`), with seven
confirmed bacterial DAM orthologues (UniProt: P0AEE8, P43870, Q9KNA4, Q8ZI12,
P26402, O34634, A1KSH9) included as outgroup in the DAM tree
(`DAM_search/06_DAM_RT_phylogeny.py`).

#### 13.5 DAM methyltransferases embedded within TEs

**Systematic methylation domain screening in TE ORFs.**  To assess the full
breadth of methyltransferase domains encoded within TE consensus sequences, the
EarlGrey TE libraries for all 238 species were first translated into all six
reading frames with EMBOSS `getorf` (minimum ORF length 100 aa).  The resulting
ORF databases were then screened by `hmmscan` against a curated set of
methyltransferase HMMs comprising three functional classes:

- **C5-cytosine DNA methyltransferase** (DNMT): PF00145 (DNA_methylase), PF17980
  (ADD_DNMT3) — the canonical eukaryotic 5mC maintenance and de-novo enzymes.
- **MT-A70 / m6A RNA/DNA methyltransferase** (MTA70): PF05063 (MT-A70) — the
  catalytic subunit of the METTL3/METTL14 m6A writer complex and its TE-encoded
  homologues.
- **DAM-class N6-adenine DNA methyltransferase**: PF05869 (Dam), PF01555
  (N6_N4_Mtase) — the bacterial-type Dam domain found in DIRS and Gypsy elements.

A TE ORF was classified as carrying a given methyltransferase domain when ≥ 1
hmmscan hit at E-value ≤ 10⁻⁵ (domain-level score) was detected.  Of the three
classes, only DAM-domain ORFs also carried TE-structural domains in the same ORF
(co-occurrence with RVT_1, rve, Chromo, or gag_pre-integrase at E ≤ 10⁻⁵).
No TE ORF carried a DNMT (C5-methyltransferase PF00145/PF17980) or MTA70 domain
together with a TE-structural domain, indicating that C5-DNMT and MTA70 fusions
within TE coding frames are either absent from or extremely rare in fungi
(`te_dnmt_search/parse_te_dnmt.py`; `te_dnmt_search/te_with_dnmt.tsv`).

**DAM-class methyltransferases in annotated proteomes.**  DAM N6-adenine
methyltransferase domain proteins (Pfam PF05869) were also identified in all 320
predicted proteomes by HMMER v3.3.2 hmmsearch:

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

**DAM domain co-fusion within DIRS ORFs.**  The phase-2 screening of EarlGrey TE
ORFs searched each species' ORF database with a combined HMM profile set including
Dam (PF05869), N6_N4_Mtase (PF01555), and the full TE-structural domain set (RVT_1,
RVT_2, rve, Chromo, gag_pre-integrase).  A DAM-TE co-fusion was defined as a single
ORF carrying ≥ 1 methylase hit AND ≥ 1 TE-context domain hit, both at E ≤ 10⁻⁵
(`DAM_search/02_hmmsearch_array.sh`, `DAM_search/03_parse_results.py`).  The
`same_orf = YES` criterion in the phase-2 master table ensures that both domains
are encoded within the same open reading frame rather than in adjacent ORFs on the
same TE consensus.  This is a key distinction: the DAM methyltransferase is a
functional passenger embedded within the retrotransposon's coding region, not a
neighbouring host gene co-annotated with the TE.

Fifteen species harboured confirmed DAM-TE genes: 14 associated with DIRS-class LTR
retrotransposons, 1 with a Gypsy LTR element.

#### 13.6 Clade-specific domain expansion and TE proximity analysis

To identify protein domains that have undergone clade-specific expansions in
parallel with TE activity, a systematic scan of all Pfam domains was performed
across the 238-species IPR.tab dataset
(`TE_analysis/expansion_analysis/00b_discover_expansions.py`).

**Discovery of expanded domains.**  For each Pfam domain with ≥ 5 total protein
hits across the dataset, two expansion tests were applied:
*(i) Cladal expansion*: a one-sided Mann–Whitney U test comparing per-species copy
counts in each of the 18 clades against all other species (rest-of-dataset).  A
domain was flagged as cladally expanded if the corrected p-value was < 10⁻⁴ (BH-FDR
across all clade × domain combinations), the fold-change over the rest-of-dataset
median was ≥ 5, and the within-clade median copy count was ≥ 3 copies per species.
*(ii) Species-specific expansion*: a domain was flagged if any single species
harboured ≥ 10 copies while the median of all other species was ≤ 1 copy.
Ubiquitous housekeeping domains (present in > 40 % of species with ≥ 2 copies),
ribosomal and proteasomal components, and the seven previously characterised
domains (CotH, TLD, Pil1, SUR7, SMG1_N, Trypsin, LPMO_10) were excluded.

The final expansion candidate table (`results/clade_expansion_table.tsv`) was
produced by `expansion_analysis/10_clade_expansion_table.py`, ranking candidates by
fold-change relative to the rest of the dataset.

**Per-species copy number screen.**  For each of the 31 domains entering the master
expansion table, per-species raw copy counts and protein IDs were extracted from all
238 IPR.tab files
(`expansion_analysis/01_copynumber_screen.py`; output: `results/expansion_copynumber.tsv`
and `results/expansion_proteins.tsv`).

**TE proximity analysis (Fig. 2E, Supp. 3C).**  For each species with at least one
gene carrying an expanded domain, the distance from every such gene to the nearest
annotated TE copy was computed with `bedtools closest -d`
(`expansion_analysis/04_proximity_analysis.py`).  Three comparisons were made per
domain × species: (A) expanded-clade species vs. non-expanded-clade species for the
same domain; (B) target-domain genes vs. a size-matched random background of other
genes in the same genome; (C) high-copy species (top 25 % by domain copy number)
vs. low-copy species (bottom 75 %) within the expanded clade.  Background genes
were sampled without replacement from the same genome.  Two-sided Mann–Whitney U
tests with BH correction at FDR 5 % were applied; results are summarised in
`results/proximity/proximity_summary.tsv`.

**TE correlation screen (Fig. 2F, Supp. 3D).**  For each expanded domain × clade
pair, the Spearman rank correlation between per-species domain copy number and
genome TE content (% genome covered by EarlGrey-annotated TEs) was computed across
all species in the dataset (`expansion_analysis/09_te_correlation_screen.py`).
Significant positive correlations (ρ > 0, BH-FDR 5 %) indicate that species with
more TE content tend to carry more copies of the expanded domain, consistent with
TE-mediated amplification.  The scatter plot and per-domain summary are shown in
Fig. 2E and Supplementary Fig. 3C–D.

**Structural analysis of pack-element candidates (Supp. 3C–D).**  For domains with
significant TE proximity and TE-content correlation, structural analysis was
performed to identify putative pack-TE (Pack-MULE or Pack-HELITRON) elements —
non-autonomous TEs that carry host gene fragments
(`expansion_analysis/06_structural_analysis.py`).  Each candidate gene was
intersected with the EarlGrey annotation; genes fully embedded within a TE interval
(overlap ≥ 90 % of gene length by a single TE copy) and shorter than the median
domain-family gene length by > 40 % were flagged as putative pack elements.

**Evidence integration and mechanism classification.**  Each domain in the master
expansion table was assigned to one of three mutually exclusive mechanism
categories based on two binary evidence flags (`expansion_analysis/08_summary_figures.R`):

- **TE-proximity evidence** (flag `has_prox`): met if pack-element candidates
  ≥ 30 OR the genomic-proximity Mann–Whitney p-value (§13.6 TE proximity analysis)
  was < 0.05 after BH correction.  This indicates that gene copies are
  non-randomly concentrated near TE insertions, consistent with TE-mediated
  amplification.

- **TE co-option/fusion evidence** (flag `has_fus`): met if ≥ 3 proteins in the
  dataset carry both the target domain and a TE-machinery domain in the same ORF
  (§13.4 fusion screen), OR if ≥ 1 TE consensus sequence in the EarlGrey library
  encodes the target domain (`n_te_consensus_with_domain > 0`).  This indicates
  direct incorporation of the host domain into a TE coding frame.

Domains satisfying both flags were classified as **"1+2: TE-proximity AND TE
co-option"**.  Domains satisfying `has_prox` only were classified as
**"1: TE-proximity mediated expansion"**.  Domains satisfying `has_fus` only
were classified as **"2: TE co-option / domestication"**.  All remaining
candidate expansions lacking either flag were classified as **"Reference: copy
expansion, no clear TE link"** and were excluded from the main TE-expansion
figure panels but retained in Supplementary Table S1.

---

### 14. Whole-genome duplication analysis

WGD detection combined three independent lines of evidence computed for each
species: (i) synonymous substitution rate (Ks) distributions of the whole-paranome
(wgd v2), (ii) collinearity-based detection of multiplicons (wgd v2 `syn` module),
and (iii) gene tree reconciliation–based duplication counts from OrthoFinder.
Evidence was integrated into a four-category classification, described in §14.6.
All WGD analyses were restricted to the 131 species with ≤ 50 assembly scaffolds
(to ensure reliable collinearity detection); the full 320-species dataset was used
only for the OrthoFinder reconciliation arm (§14.3).

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
(§14.4) and age-correlation analyses (§14.5).

#### 14.3 OrthoFinder gene tree reconciliation — duplication mapping

All 320 species (including fragmented assemblies) were included in the OrthoFinder
3.1 run described in Section 2.  The `Gene_Duplication_Events/Duplications.tsv`
file from OrthoFinder records every gene duplication event inferred by gene tree /
species tree reconciliation (using the STRIDE algorithm; Emms & Kelly 2018), with
the column `Species Tree Node` indicating the node on the ML species tree (Section 2)
at which the duplication was mapped.

Duplication rates per species-tree node were computed as normalised counts
(`scripts/compute_paper_like_duplication_rates.py`).  For each node n in the tree,
the rate R_n is the number of duplication events mapped to that node (from
`Duplications.tsv`) divided by the number of eligible orthogroups for that node:

`R_n = duplication events at node n / eligible orthogroups at node n`

Eligible orthogroups are those in which species from both descending clades of the
node are represented, excluding any orthogroup where a single species contributes more
than 5 copies (the `--max-species-copies 5` filter, which removes lineage-specific
gene family expansions that would otherwise inflate node-specific duplication counts).

Species-tree nodes at which the normalised rate R_n was classified as an outlier
(Robust Z-score ≥ 1.5) relative to all nodes of the same type (terminal or
internal) were tagged as **duplication-peak nodes** (`scripts/detect_duplication_peaks.py`).
The Robust Z-score is a median-based measure of how many "typical deviations" a value
lies from the median of the distribution, computed as:

`Z = 0.6745 × (value − median) / MAD`

where MAD is the median absolute deviation (MAD = median of |R_n − median(R)|) across
all nodes of the same type, and the constant 0.6745 rescales MAD to be equivalent to
the standard deviation under a normal distribution.  Per-species signals were summarised
as: `deep_peak_nodes` (count of internal tree nodes at depth ≥ 2 from the tip classified
as duplication-peak), `recent_peak_nodes` (terminal/depth-1 nodes classified as
duplication-peak), and `max_peak_z` (maximum Robust Z-score over all nodes of that
species).

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

The per-species paranome Ks distribution (from §14.2; maximum 5,000 randomly
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

The number of components k was varied from 1 to 5.  Model selection used the
Bayesian Information Criterion (BIC), which penalises model complexity to avoid
overfitting:

`BIC(k) = k × ln(n) − 2 × ln(maximum likelihood)`

where n is the number of Ks data points.  The k minimising BIC was selected.  A GMM component was reported as a candidate WGD Ks peak if its mixture weight
(proportion of data points assigned to it) was ≥ 0.10 AND its mean Ks fell within
[0.05, 3.0], thereby excluding near-zero components (gene conversion / very recent
duplicates) and saturated components above the reliable estimation range.  The event
timing was classified as `ancient_deep` if all WGD-peak components had a mean Ks
> 0.50, `recent_terminal` if the mean Ks was ≤ 0.50, or `mixed_deep_and_recent` if
both conditions occurred in different components.

#### 14.5 Collinearity analysis — wgd v2 `syn` and multiplicons

For the 131 filtered species (§14.1), intra-genomic collinearity was assessed using
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
The synteny evidence score for each species was defined as the larger of the two
Robust Z-scores (described in §14.3) computed for that species’ multiplicon count
and segment count, respectively, across the 131-species distribution:

`S = max(Robust Z-score of multiplicon count, Robust Z-score of segment count)`

A species was considered to have **synteny evidence** if `syn_completed = 1` AND
`n_multiplicons > 0`.

The anchor-point Ks distribution (Ks values of anchor gene pairs within multiplicons)
provides a more precise WGD age estimate than the full paranome, because it filters
out background small-scale duplications.  The modal Ks of the anchor-Ks distribution,
extracted as the GMM peak component (§14.4 applied to anchor pairs only), was used
as the primary WGD age estimate.

#### 14.6 Multi-evidence WGD classification

Species were classified into four WGD categories by combining the three evidence
streams (`scripts/classify_wgd_evidence.py`; `scripts/integrate_wgd_evidence_paper_like.py`):

| Classification | Criterion |
|---|---|
| **Strong WGD** | duplication rate R >= 0.30 AND synteny completed AND n_multiplicons > 0 |
| **Uncertain — duplication only** | duplication rate R >= 0.30 AND (n_multiplicons = 0 OR syn_completed = 0) |
| **Weak / ancient WGD** | duplication rate R < 0.30 AND n_multiplicons > 0 AND syn_completed = 1 |
| **No WGD** | duplication rate R < 0.30 AND (n_multiplicons = 0 OR syn_completed = 0) |

The duplication rate threshold of 0.30 was calibrated against 12 fungal
species with previously published WGD evidence (*Saccharomyces cerevisiae*, the
Rhizopus microsporus complex, *Phytophthora infestans*, etc.); all 12 scored
R > 0.40 in our normalised OrthoFinder analysis.  Species with
`event_timing = mixed_deep_and_recent` (§14.4) were flagged for manual review,
as this pattern can arise from two independent rounds of duplication or from
allopolyploidy followed by diploidisation.

#### 14.7 TE breakpoint enrichment at WGD synteny boundaries

For each WGD-positive species (Strong WGD or Weak/ancient), the boundaries of
intra-genomic collinear blocks (wgd v2 multiplicons, §14.5) were used to define
breakpoint windows — i.e. the regions of the genome where multiplicon adjacency
is interrupted.  Concretely, a breakpoint was defined as the chromosomal interval
between the last anchor-point coordinate of one multiplicon and the first
anchor-point coordinate of the next multiplicon on the same scaffold, sorted by
genomic position.  Each breakpoint midpoint was expanded by ±5,000 bp.  One hundred
randomly sampled windows of the same length distribution (sampling with replacement
from the scaffold) were drawn per breakpoint as a background.

TE density in each window W was computed as the total number of base pairs covered
by TE annotations within the window, divided by the window length:

`TE density = total TE bp in window / window length in bp`

A two-sided Wilcoxon signed-rank test (scipy.stats.wilcoxon) compared matched
breakpoint vs. random TE densities (one observation per breakpoint, paired).
Per-superfamily enrichment: same Wilcoxon test restricted to superfamilies with
≥ 10 windows containing ≥ 1 bp of that superfamily; BH-corrected FDR 5 %.

#### 14.8 TE burst — WGD age correlation

**WGD age estimation from anchor-Ks peaks and a time-calibrated molecular clock.**
WGD ages in Mya were estimated from the modal Ks of the anchor-pair GMM peak
component (§14.5), converted to absolute time using a synonymous substitution rate
calibrated against the fossil-constrained timetree of Szánthó *et al.* (2025).
The conversion follows the standard molecular clock formula:

`WGD age (Mya) = anchor_Ks_peak / (2 × r_syn)`

where `r_syn` is the synonymous substitution rate per site per year.

**Branch-specific rate calibration.**  For each WGD-positive lineage, `r_syn` was
derived by matching the lineage to the nearest anchored node in the Szánthó *et al.*
(2025) time-calibrated ultrametric chronogram (Supplementary Material, NEXUS tree).
Node ages were extracted as cumulative branch lengths from each internal node to
the tips; in an ultrametric tree all tips are at time 0, so each node's age equals
the sum of branch lengths connecting it to any terminal.  The key calibration ages
used are listed below (all in Ma):

| Szánthó node | Fungal group | Age (Ma) |
|---|---|---|
| Sz_199 | Mucor s.s. crown | 100 |
| Sz_198 | Mucor + Rhizopus | 143.9 |
| Sz_197 | Crown Mucorineae (excl. Phycomyces) | 260.9 |
| Sz_196 | Crown Mucorales | 283.5 |
| Sz_195 | Broader Mucorales | 307.3 |
| Sz_194 | Mucoromycotina (excl. Endogonales + Glomeromycota) | 514.0 |
| Sz_212 | Glomeromycota crown | 471.7 |
| Sz_217 | Mortierella crown | 74.5 |
| Sz_141 | Crown Ascomycota | 839.7 |
| Sz_184 | Saccharomycotina (Scer + Dbru) | 342.4 |
| Sz_59  | Crown Basidiomycota | 790.1 |

For each WGD lineage, the patristic distance between the WGD clade and its closest
non-WGD sister clade in the IQ-TREE ML species tree was measured.  This patristic
distance was then compared to the corresponding divergence time from the Szánthó
chronogram (interpolating linearly between bracket nodes where a direct age
assignment was not available, using branch-length-proportional interpolation).
The per-lineage `r_syn` was then computed as:

`r_syn = patristic_distance / (2 × node_age_Myr × 10⁶)`

(the factor 10⁶ converts Myr to years).  For the Mucoromycota WGD nodes
(WGD-1, WGD-2, WGD-3 mapped to OrthoFinder nodes N141, N153, N163), a global
calibration of µ = 1.5 × 10⁻⁹ substitutions/site/year was used (calibrated by
matching the WGD-1 Ks peak of 0.785 to the 261 Ma Szánthó age).  WGD-4 (node N172,
Backusella; anchor Ks = 0.616) was NOT assigned an age via Ks because its
topological placement within the WGD-3 clade makes a direct Ks-to-age conversion
contradictory; instead, its age was estimated by branch-length interpolation between
the nearest calibrated ancestor (N163, 104 Ma) and the species-tree tip, giving a
rate-adjusted estimate of ~82 Ma (`WGD/scripts/figures/fig_gene_dup_rate_ma.py`).

The resulting three confirmed Mucoromycota WGD ages are:
WGD-1 ~261 Ma, WGD-2 ~157 Ma, WGD-3 ~104 Ma.

**TE burst — WGD temporal correlation.**  For each WGD-positive species, the
TE burst age was summarised as the kernel density estimation (KDE) mode of the
per-copy Kimura divergence K distribution (§6.2) for each TE superfamily, using
`scipy.stats.gaussian_kde` with Scott's bandwidth rule.  The KDE mode approximates
the most recent period of elevated TE activity (lowest K = most recently active).
The KDE mode was converted to Mya using the same per-lineage `r_syn`:

`burst_age (Mya) = K_KDE_mode / (2 × r_syn × 10⁶)`

Spearman's rank correlation between WGD age (Mya) and TE burst KDE mode (Mya)
was computed per superfamily across all WGD-positive species; BH-corrected FDR 5 %
across superfamilies.  Superfamilies with ≥ 5 WGD-positive species passing the KDE
estimation threshold (≥ 20 TE copies per species) were included.
Results are shown in Fig. 5C.

#### 14.9 Normalised gene duplication rates per phylogenetic node

The normalised duplication rate per node (Section 14.3) was plotted on the species
tree as a heat-map column, using the same colour scale across all nodes.  For
visualisation, ancestral node rates were projected onto the LG+G4 ML species tree
with ete3 v3.1.3 (`NodeStyle` colouring, `TreeStyle` rectangular layout).
Duplication-peak nodes (Robust Z ≥ 1.5, Section 14.3) are highlighted with a filled
circle; circle size is proportional to the Robust Z-score.

---

### 15. Horizontal transposon transfer (HTT)

Horizontal transposon transfer (HTT) detection followed a nine-step computational
pipeline (`HTT/scripts/` and `legserv/.../HTT_hpc_scripts/`), closely inspired by
the approach of Romeijn *et al.* (2025) for fungi.  Each step is described below
in pipeline order.

#### 15.1 Species phylogeny and evolutionary divergence calibration

A reference phylogenetic tree for all 320 EDF species was inferred as described in
Section 2 (LG+G4 model, 1000 UFBoot replicates; `HTT/scripts/01_phylo_tree.sh`).
The tree was rooted on the outgroup (*Naegruberia gruberi*, clade HO).

To calibrate evolutionary divergence between each pair of clades at the node level,
the synonymous substitution rate (Ks) of core single-copy genes was computed for
each internal node of the species tree.  For each node, single-copy OrthoFinder
orthogroups (Section 2) common to both descendant sub-clades were used as the
gene set — this replaces BUSCO-based selection and ensures coverage even for
lineages where BUSCO databases are absent or incomplete.  For species lacking
nucleotide CDS files in the BUSCO output, CDS sequences were extracted directly
from the GFF + genome FASTA using the OrthoFinder protein assignments
(`HTT/scripts/03a_extract_cds_for_ogs.py`, `03a_extract_cds_for_missing_busco.py`).

For each orthogroup at a node, one representative CDS per descendant clade was
selected by maximum sequence length (or randomly for ties).  The two CDS were
aligned at the protein level with MAFFT v7.520 (`--auto`), then a codon alignment
was built with PAL2NAL v14 (`-nogap`).  Ks values were estimated with the
`kaks` function in seqinr (v4.2-36) in R (`HTT/scripts/05_ks_seqinr.R`).  The
0.5th percentile of the Ks distribution per node was used as a **node Ks threshold**;
node pairs with a threshold < 0.01 were collapsed (treated as a single taxonomic
unit) to avoid false-positives from near-identical species
(`HTT/scripts/02_div_times.sh`).

#### 15.2 TE annotation and protein-domain classification

TEs were annotated in all 320 genomes with EarlGrey as described in Section 5.
For the HTT pipeline, TE copies < 300 bp were discarded
(`HTT/scripts/04_te_annotation.sh`).  All retained TE sequences were extracted
from genome FASTAs with BEDtools v2.31.1 `getfasta` and compiled into a combined
FASTA (`HTT/scripts/05_prep_te_data.py`).

To restrict HTT analysis to copies likely to encode functional TE proteins,
TE-related protein domains were identified in each copy.  TE sequences were
translated in all six frames with `transeq` (EMBOSS v6.6.0, `-frame 6`).
Translated sequences were searched against:
(i) 142 Conserved Domain Database (CDD) profiles with RPSTBLASTN
(BLAST+ v2.12.0, `-evalue 0.001 -max_target_seqs 1`);
(ii) 110 TE-related Pfam HMM profiles (same domain set as §13.1) with HMMER v3.4
(`hmmsearch -E 0.001`).
Only TE copies with at least one domain hit passing the e-value threshold were
retained for the HTT analysis
(`HTT/scripts/05_te_classification.sh`).

#### 15.3 All-vs-all BLASTN and hit chaining

Retained domain-containing TE sequences were searched against each other in an
all-vs-all BLASTN (BLAST+ v2.12.0):

```
blastn -dbsize 10000000000 -evalue 10 \
  -max_target_seqs 100000 \
  -perc_identity 75 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send bitscore"
```

Cross-clade hits (i.e. query and subject from different clades) were retained if
percent identity ≥ 75 %, alignment length ≥ 300 bp, and bitscore ≥ 200.
Fragmented collinear hits between the same pair of TE sequences were chained into
single extended alignments allowing a maximum gap and maximum overlap of 600 bp.
Pairs were further filtered to require that the chained alignment covered ≥ 60 %
of both the query and subject TE length (to exclude partial matches that could
reflect conserved domain fragments).  This yielded **50,143 chained cross-clade
hit pairs** as input to Ks estimation.

#### 15.4 Ks-based HTT candidate selection

For each chained hit pair, Ks was estimated using the same codon-alignment +
seqinr pipeline as in §15.1, restricted to the aligned domain region.  Protein
domains present in both TEs had to overlap precisely the same CDD/Pfam coordinate
ranges on both sides of the alignment; no frameshift-introducing gaps were
permitted.  A candidate HTT pair was retained if its estimated Ks was lower than
the 0.5th-percentile Ks threshold of the most recent common ancestor (MRCA) node
of the two species' clades (i.e. the TE copies are more similar than expected for
any gene evolving neutrally since the two lineages diverged).  This yielded
**223 HTT candidate pairs** spanning **20 clade pairs**
(`HTT/scripts/06_htt_candidates.sh`).

#### 15.5 Flanking locus conservation — non-vertical transmission test

For each HTT candidate pair, 10,000 bp of flanking sequence on each side of the
TE insertion was extracted.  Flanks from the putative acceptor were aligned to the
putative donor genome and proteome with MUMmer v3.1 (`nucmer`, mode `--mum`,
minimum match 40 bp, minimum cluster 100 bp; `HTT/scripts/07_flankid.sh`).
Alignments shorter than 200 bp were discarded.  For non-HTT control pairs, an
equal number of same-family pairs from within the same clade were matched and
subjected to the same flanking alignment analysis.

A flanking conservation score for each pair was computed as the total aligned
length across both flanks divided by 20,000 bp.  Pairs where the flanking
conservation was exactly 0 % and TE identity was > 75 % — the strictest evidence
for non-vertical transfer — were retained as the **high-confidence HTT set**
(91 pairs; `htt_candidates_flank_0pct_teid_gt75.tsv`).

#### 15.6 Community detection — grouping into independent transfer events

Because multiple HTT pairs within the same biological event inflate the candidate
count, HTT pairs were clustered into independent transfer **communities** using
the approach of Romeijn *et al.* (2025).

For each pair of clades sharing at least one HTT candidate, an intra-clade
all-vs-all BLASTN was performed (`-dbsize 10,000,000,000 -evalue 10
-max_target_seqs 600,000`) to obtain within-clade identity scores.  HTT pairs
were represented as nodes in a graph; two nodes were connected by an edge if: (a)
they involved the same TE copy, OR (b) the between-clade identity of the TE pair
exceeded the within-clade identity of the same TE (i.e. the cross-species copy
is more similar than copies within the same species, ruling out independent
parallel evolution).  The Louvain community detection algorithm
(`igraph` v0.11.8, `cluster_fast_greedy`) was applied to identify communities.
Communities with inter-community edge density > 0.05 were merged to minimise
false positives.  This yielded **20 HTT communities** across 20 clade pairs
(`HTT/scripts/08_htt_clustering.sh`; `htt_pairs_with_community.tsv`).

#### 15.7 Minimal HTT events — Dollo parsimony on the species tree

To account for the possibility that a single ancestral transfer produced multiple
descendant lineages that each appear as independent transfers, a second-level
network analysis was applied.  Communities were treated as nodes; communities
sharing TE copies were connected.  Single-linkage clustering (networkX v2.5.1)
identified clusters of related communities.  For each cluster, the species tree
was pruned to include only the species and collapsed nodes present in that cluster
(ete3 v3.1.3).  Per node of the pruned tree, the minimum number of independent
transfers required to explain the observed pattern parsimoniously was counted:
if multiple transfers from different clade pairs involved the same TE community,
one event sufficed; if the same clade pair was involved in multiple independent
communities, that exact count was recorded.  Summing across all nodes yielded a
**minimum of 15 independent HTT events** across the EDF dataset
(`HTT/scripts/09_minimal_htt_events.sh`; Fig. 6A).

**HTT node distribution.**  Of the 15 minimal events, the largest node
(node_32) accounts for 6 events, followed by node_38 (3 events); the remaining
7 events are distributed singly across 7 additional tree nodes.

#### 15.8 TE type and clade enrichment

For each HTT community, the TE superfamily classification of the involved elements
was extracted from the `all_tes_info.tsv` table.  The 20 clade pairs involved span
clades CH, MU, MO, AS, BL, BS, AP, EN, GL, HO, KX, NU, and ZO, with MO–MU being
the most represented pair (involved in multiple communities).  Fisher's exact tests
(BH-FDR 5 %) were used to assess enrichment of specific TE superfamilies and clade
pairs relative to their expected frequencies under a null model of random pairing
(`fig_htt_te_enrichment.py`).

#### 15.9 HTT visualisation

Chord diagrams illustrating the clade-pair network of HTT events were generated
with the `circlize` R package v0.4.16 (`HTT/scripts/plot_htt_chord.R`,
`plot_htt_chord_validated.R`).  Fig. 6 shows the validated high-confidence events;
Supplementary Fig. S8 shows the full community set.  A locus-level synteny
comparison for selected HTT pairs was generated with MUMmer + custom R plotting
(`HTT/scripts/plot_htt_synteny.R`, `plot_htt_synteny_v2.R`).

---

### 16. Statistical framework

All statistical tests were two-tailed.  Multiple comparisons were corrected by the
Benjamini–Hochberg (BH) procedure at FDR ≤ 5 % whenever ≥ 5 hypotheses were tested
simultaneously (scipy.stats.false_discovery_control(method='bh')).  Effect sizes
are reported as: Spearman rho or Pearson r (correlations), odds ratio with
95 % Wilson confidence intervals (Fisher's exact), and rank-biserial correlation r
(Wilcoxon / Mann–Whitney tests).

Unless otherwise stated, all analyses were implemented in Python 3.11.3 with the
following library versions: numpy 1.26.2, pandas 2.1.4, scipy 1.11.4, scikit-learn
1.3.2, matplotlib 3.8.2, seaborn 0.13.2, pyBigWig 0.3.18, ete3 3.1.3,
networkx 3.2.1.  R 4.4.0 was used for visualisation with ggplot2 4.0.0, ggtree 3.10,
circlize 0.4.16, and ComplexHeatmap 2.18.0.  Network visualisation was performed with
Cytoscape 3.10.2 with the yFiles Organic layout.

---

### 17. Published data accession IDs

All publicly available datasets re-analysed in this study are listed below.
Genome assemblies and gene annotations are provided in the supplementary genome
catalogue (Source Data; §17 of Methods).

#### 17.1 Chromatin and epigenomic datasets

| Species | Strain/Assembly | Data type | Accession | Reference |
|---|---|---|---|---|
| *Rhizopus microsporus* | ATCC11559 | H3K9me3 ChIP-seq (IP + Input) | — (local; reads in `data/Extra/ChIPRMH3K9me3/`) | This study |
| *Rhizopus microsporus* | ATCC11559 | MNase-seq | — (legacy WIG; `data/nucleosomes/Rhizopus_microsporus.wig`) | This study |
| *Neurospora crassa* | OR74A | MNase-seq | SRR13067301, SRR13067302, SRR13067303 | Catania *et al.*, 2020 |
| *Botrytis cinerea* | B05.10 | MAINE-seq (MNase equivalent) | SRR13224644, SRR13224647, SRR13224648 | Catania *et al.*, 2020 |
| *Microbotryum eburneum* | — | ATAC-seq | SRR28069093, SRR28069094, SRR28069106 | Mondo *et al.*, 2024 |
| *Microbotryum eburneum* | — | Nuc-seq (in vitro) | SRR28069066, SRR28069068 | Mondo *et al.*, 2024 |
| *Rhodotorula toruloides* | — | Nuc-seq (in vitro) | SRR28069068 | Mondo *et al.*, 2024 |
| *Catenaria anguillulae* | PL171 | Nuc-seq (in vitro) | SRR7476903 | Mondo *et al.*, 2018 |
| *Ustilago hordei* | — | Nuc-seq (in vitro) | SRR7472154 | Mondo *et al.*, 2018 |
| *Morchella importuna* | — | Nuc-seq (in vitro) | SRR7472065 | Mondo *et al.*, 2018 |

#### 17.2 RNA-seq datasets for ChimeraTE analysis

| Species | Data type | SRA Accession | Reference |
|---|---|---|---|
| *Mucor lusitanicus* (MU_Mucci3) | mRNA-seq (interleaved, 3 reps) | — (MycoMobilome; local) | This study |
| *Rhizopus microsporus* (MU_Rhimi59) | mRNA-seq (interleaved, 3 reps) | — (MycoMobilome; local) | This study |
| *Phycomyces blakesleeanus* (MU_Phybl4) | mRNA-seq (interleaved, 3 reps) | — (MycoMobilome; local) | This study |
| *Amoeboaphelidium protococcarum* (AP_Amoprot) | mRNA-seq | SRR18053693 | Rybakov *et al.*, 2022 |
| *Blastobotrys brendanlejeunia* (BL_Blabri1) | mRNA-seq | SRR6057017 | Mondo *et al.*, 2017 |
| *Basidiobolus meristosporus* (BS_Basme92) | mRNA-seq | SRR10567091 | Corsaro *et al.*, 2019 |
| *Batrachochytrium dendrobatidis* (CH_Batde5) | mRNA-seq | SRR27683879 | JGI/NCBI |
| *Batrachochytrium salamandrivorans* (CH_Batsa1) | mRNA-seq | SRR3706726 | Stöhr *et al.*, 2016 |
| *Piromyces finnis* (CH_Pirfi3) | mRNA-seq | SRR1748130 | Mondo *et al.*, 2017 |
| *Enterocytozoon bieneusi* (EN_Entmai1) | mRNA-seq | SRR9001799 | NCBI/public |
| *Zoophthora radicans* (EN_Zoorad1) | mRNA-seq | SRR8189665 | NCBI/public |
| *Geosiphon pyriformis* (GL_Geopyr1) | mRNA-seq | SRR12018970 | NCBI/public |
| *Rhizophagus clarus* (GL_Rhicl1) | mRNA-seq | DRR121634 | NCBI/DDBJ |
| *Rhizophagus irregularis* (GL_RhiirA4) | mRNA-seq | SRR31872297 | NCBI/public |
| *Coelomomyces lativittatus* (KX_Coemoj1) | mRNA-seq | SRR7140835 | NCBI/public |
| *Kickxella alabastrina* (KX_Kicala1) | mRNA-seq | SRR6057250 | Mondo *et al.*, 2017 |
| *Encephalitozoon cuniculi* (MI_Enccu1) | mRNA-seq | SRR769608 | NCBI/public |
| *Spraguea lophii* (MI_Sprlo1) | mRNA-seq | SRR916958 | NCBI/public |
| *Vavraia culicis* (MI_Vavcu1) | mRNA-seq | SRR1047418 | NCBI/public |
| *Discodinia discoidea* (MO_Disdec1) | mRNA-seq | SRR13647699 | NCBI/public |
| *Mortierella elongata* (MO_Morel2) | mRNA-seq | DRR045521 | NCBI/DDBJ |
| *Fonticula alba* (NU_Fontalba) | mRNA-seq | SRR554356 | NCBI/public |
| *Pilobolus crystallinus* (ZO_Pipcy3) | mRNA-seq | SRR5506705 | NCBI/public |

#### 17.3 DAP-seq transcription factor binding data

| Species | Data type | Source |
|---|---|---|
| *Rhizopus microsporus* ATCC11559 | DAP-seq peaks for 58 TFs (MEME motif files) | This study; `data/TFsDAP/motifs/` |

#### 17.4 PacBio kinetics (6mA methylation)

| Species | Data type | Source |
|---|---|---|
| *Phycomyces blakesleeanus* | PacBio CLR + IPD kinetics | This study |
| *Rhizopus microsporus* | PacBio CLR + IPD kinetics | This study |

#### 17.5 Small RNA sequencing and RNAi transcriptomes

Public sRNA-seq and paired mRNA-seq data for *Mucor* species were retrieved from
the following BioProjects to characterise TE silencing by the RNAi pathway:

| BioProject | GEO | Species / Strain | Data type | SRA accessions | Condition | Reference |
|---|---|---|---|---|---|---|
| PRJNA1220135 | GSE288970 | *Mucor lusitanicus* MU636 (WT) + MU412 (r3b2Δ) | sRNA-seq (ncRNA-Seq, NovaSeq 6000) | SRR32258784–SRR32258789 (6 runs) | WT vs r3b2Δ, mycelium | Tahiri, Lax *et al.*, 2025 (PMID 40141380) |
| PRJNA1220135 | GSE288970 | *Mucor lusitanicus* MU636 (WT) + MU412 (r3b2Δ) | mRNA-seq (NovaSeq 6000) | SRR32258790–SRR32258795 (6 runs) | WT vs r3b2Δ, mycelium | Tahiri, Lax *et al.*, 2025 |
| PRJNA1263436 | GSE297279 | *Mucor lusitanicus* MU636 (WT) + MU412 (r3b2Δ) | sRNA-seq (ncRNA-Seq, NovaSeq 6000) | SRR33582828–SRR33582835 (8 runs) | WT vs r3b2Δ, spores ± macrophages | Tahiri, Lax *et al.*, 2026 (PMID 41481777) |
| PRJNA1263436 | GSE297279 | *Mucor lusitanicus* MU636 (WT) + MU412 (r3b2Δ) | mRNA-seq (NovaSeq 6000) | SRR33582836–SRR33582847 (12 runs) | WT vs r3b2Δ, spores ± macrophages | Tahiri, Lax *et al.*, 2026 |
| PRJNA903107 | — | *Mucor lusitanicus* CBS277.49 (MU402 + RNAi mutants) | sRNA-seq (miRNA-Seq, NextSeq 500) | SRR22390021–SRR22390049 (19 runs) | WT + RNAi mutant series | Navarro-Mendoza *et al.*, 2023 (PMID 37014858) |
| PRJNA903107 | — | *Mucor lusitanicus* CBS277.49 (MU402 + RNAi mutants) | mRNA-seq (paired-end, NextSeq 500) | SRR22390026–SRR22390078 (20 runs) | WT + RNAi mutant series | Navarro-Mendoza *et al.*, 2023 |
| PRJNA453739 | GSE113706 | *Mucor circinelloides* WT + pyrF/pyrG epimutants | sRNA-seq (HiSeq 2000) | SRR7068267–SRR7068276 (10 runs) | WT vs RNAi epimutants | Chang *et al.*, 2019 (PMID 30742617) |
| PRJNA243024 | GSE56353 | *Mucor lusitanicus* WT + FK506 epimutants | sRNA-seq (HiSeq 2000) | SRR1209506–SRR1209511 (6 runs) | WT vs FK506-resistant epimutants | Calo *et al.*, 2014 (PMID 25079329) |

GSE288970 (PRJNA1220135) was the primary dataset for TE silencing analyses: sRNA reads were mapped to the *M. lusitanicus* CBS277.49 genome (GCF_001624055.1) and the EarlGrey TE annotation to quantify siRNA targeting of TE families; paired mRNA-seq was used to measure TE derepression in r3b2Δ mutants (35 direct NCRIP targets identified by Tahiri *et al.*).



## Part II — Figure Legends

### Figure 1. Genomic diversity and synteny across 320 EDF genomes

**(A)** Circular maximum-likelihood phylogenetic tree of all 320 genomes (Section 2;
LG+G4 model, IQ-TREE 2; `Section2_TElandscape/figures/circular_tree_v4.4.pdf`).
Tips are coloured by clade following the 18-clade colour scheme; bootstrap support
≥ 95 is indicated at internal nodes.  The 18 clade abbreviations are labelled on the
outer ring.
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
individual points.  Statistical comparisons: Kruskal–Wallis H test across 18 clades;
post-hoc pairwise Dunn test, Benjamini–Hochberg (BH) FDR correction; significance
brackets shown for clade pairs with BH-corrected q < 0.05.
Source data: Fig2C_kimura_ri.

**(D)** Clustered Kimura divergence summary across 238 species, showing the mean
Kimura divergence per TE superfamily × species combination as a compact heat-map
(Section 6.2).  Species and superfamilies are ordered by hierarchical clustering.
Source data: Fig2D_kimura_summary.

**(E)** Two evidence streams: statistical proximity × structural pack-elements —
panel D of `expansion_analysis/results/figures/fig2_te_proximity_mechanism.pdf`
(Section 13.4).  Scatter plot of 41 candidate domain families: x-axis = Gene–TE
proximity significance (−log₁₀ BH-corrected Mann–Whitney U p, expanded-clade gene
copies vs. a size-matched random background in the same genome); y-axis = log₂(n+1)
gene copies fully embedded within a TE annotation.  Dashed reference lines indicate
the BH-FDR 5% significance threshold.  Top-right quadrant: domains statistically
over-represented near TEs AND with TE-embedded copies; labels show domain name and
TE class (colour).
Source data: Fig2E_proximity_scatter.

**(F)** TE subfamily in symmetric flanks (41 domains) — panel C of the same file
(Section 13.4).  Bubble chart: rows = TE superfamilies; columns = host-domain
families; bubble size ∝ % of flanking TE copies; colour = TE class.
Source data: Fig2F_proximity_flanks.

**(G)** Sex locus (MAT) analysis — panels A and C of
`MAT_locus/figures/sex_locus_figure.pdf`.
Panel A: genomic architecture of the MAT locus in representative species, showing
the gene content and orientation of the mating-type region.
Panel C: TE class composition (bar chart by superfamily) of the TE copies
found within or immediately flanking the MAT locus regions shown in panel A.
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
enrichment level (Section 10.3).  Format as in (C).  Kolmogorov–Smirnov test
comparing the distributions of the lowest H3K9me3 class (Q1) versus the highest
(Q4); D statistic and BH-corrected p shown in the panel.
Source data: Fig3E_gene_dist_h3k9.

**(F)** Adjacent gene expression by H3K9me3 peak presence/absence (Section 10.4).
Format as in (D).  BH-corrected Mann–Whitney U test per quartile comparison;
significance brackets shown above boxes.
Source data: Fig3F_adj_gene_expr_h3k9.

**(G)** TE-associated methylated adenine cluster (MAC) expression at transcription
start sites (Section 12.1).  Violin plot comparing expression of TEs with ≥ 1 MAC
in the TSS window (−150 to +400 bp) versus TEs without MACs, for *Phycomyces
blakesleeanus* and *Rhizopus microsporus* grown under mycelium dark conditions.
MAC-associated TEs show distinct expression profiles relative to non-MAC TEs.
Mann–Whitney U test (two-sided) comparing MAC-TSS versus non-MAC-TSS TE expression;
BH-corrected p shown above the violins.
Source data: Fig3G_mac_tss.

**(H)** Metaplot of 5mC levels (Q1 vs. Q4 methylation bins) across TcMar-Tc1 and
MULE-MuDR superfamily insertions (Section 9.5).  Each line represents one
methylation quartile bin; shaded region: TE body; bands: ± 1 s.e.m.
Source data: Fig3H_metaplot.

**(I)** TE families enriched at methylated transcription start sites (Section 10.4).
Bar chart of TE superfamilies most frequently overlapping methylated TSS windows,
per species.  Enrichment significance assessed by hypergeometric test against the
background frequency of each TE superfamily at all TSS; BH-FDR 5%; only superfamilies
with q < 0.05 and enrichment ratio > 2 are labelled.
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

**(A)** Co-option gain/loss tree (Section 13.1;
`Co-option/figures/figA_tree_cooption_v2.pdf`).  Phylogenetic tree at clade level
(tips labelled by clade name).  At each ancestral node, the number of co-option
orthogroup families conserved from the ancestor and the number of newly acquired
families are annotated.  At the tip of each branch, the total number of distinct
orthogroups with co-opted TE-derived genes within that clade is indicated.
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
BH correction; colour: median −log₁₀ FDR.  Fisher's exact test per TF × TE-class
combination (58 TFs × 11 TE classes); BH-FDR 5% across all 638 combinations;
only combinations with BH-corrected q < 0.05 and log₂ odds ratio > 0 are shown.
Source data: Fig4I_tf_enrichment.

**(J)** deepTools metaplot of HMGb3 (HMGb_1228) DAP-seq signal over LINE/I-Jockey
TE insertions (Section 13.3).  Solid line: mean signal; shaded band: ± 1 s.e.m.;
vertical dashed lines: TE body boundaries.  Enrichment significance: permutation-based
p-value comparing the observed mean signal over all LINE/I-Jockey loci against 1,000
shuffled-peak controls; one-sided p shown in the panel.
Source data: Fig4J_hmgb3_metaplot.

**(K)** MEME motif logo for the HMGb3 binding motif enriched within LINE/I-Jockey
TE sequences (Section 13.3).  Nucleotide height proportional to information content.

**(L)** Clade-level heat-map of HMGb3 occupancy across Jockey superfamily insertions
(Section 13.3).  Rows: EDF clades; columns: TE subfamilies within Jockey; colour:
mean normalised DAP-seq signal.
Source data: Fig4L_hmgb3_heatmap.

---

### Figure 5. Whole-genome duplication and TE dynamics in Mucoromycota

**(A)** Gene duplication rate tree for Mucoromycota (≤ 200-scaffold filter;
Section 14.9).  The species tree topology was inferred with IQ-TREE2 under the
LG+G4 model (1,000 UFBoot replicates; Section 2) and time-calibrated using a
separate reference tree provided by the authors; branch lengths are therefore
proportional to absolute time (Myr).  Normalised duplication rates per
phylogenetic node (Section 14.3) were computed as the number of duplication
events mapped to each node from OrthoFinder's `Duplications.tsv`, divided by the
number of eligible orthogroups (≥ 2 sequences, max. 5 copies per species to
exclude species-specific expansions).  The tree is drawn with ete3 v3.1.3 in
rectangular layout (`TreeStyle`); internal nodes and terminal branches are
coloured on a continuous heat-map scale from low (blue) to high (red) normalised
duplication rate, giving a visual proxy for duplication events per million years.
Nodes identified as duplication-rate outliers (Robust Z-score ≥ 1.5 of the
node-level duplication rate distribution binned by tree depth; Section 14.3) are
marked with filled circles whose diameter is proportional to the Robust Z-score.
WGD-candidate species were identified by integrating three independent lines of
evidence (Section 14.6): (i) elevated normalised duplication rate (R ≥ 0.30),
(ii) significant Ks peaks in the paranome (wgd v2 `dmd` + Gaussian mixture model;
Section 14.4), and (iii) collinearity-based multiplicon detection (wgd v2 `syn`;
Section 14.5).  Species were classified as **Strong WGD**, **Weak/ancient WGD**,
**Duplication-only**, or **No WGD** (Section 14.6, Table §14.6), with the
classification indicated by a coloured strip beside the tree (red = strong,
orange = weak, yellow = questionable, grey = none).  An adjacent bar chart shows
the number of co-linear paralog pairs (multiplicon anchor points, wgd v2 `syn`;
Section 14.5) per species.
Source data: Fig5A_paralogs.

**(B)** Syntenic block structure and TE content (Section 14.5).  Ribbon diagram
or dot-plot showing syntenic block coverage relative to total genome length, with
TE density per block overlaid as a colour track, per WGD species.
Source data: Fig5B_syntenic_blocks.

**(C)** TE enrichment at WGD synteny breakpoints — log₂ fold-change over randomised
windows, per TE superfamily, per species (Section 14.7).  Colour: BH-corrected
significance; dot size: absolute log₂FC.  Statistical test: Wilcoxon rank-sum test
comparing TE density (bp per window) at synteny breakpoint windows versus 10,000
randomly placed windows of equal size per species; BH-FDR 5% across all
superfamily × species combinations.
Source data: Fig5C_te_breakpoints.

**(D)** TE insertion age at WGD breakpoints vs. enrichment (n = 18 species; Section
14.7–14.8; panel C of `fig_te_dynamics_combined.pdf`).  Scatter plot: x-axis =
ΔKimura80 (mean Kimura divergence at breakpoints − genome-wide mean; Δ < 0 means
younger insertions at breakpoints); y-axis = TE enrichment at breakpoints (log₂FC
over randomised windows).  Spearman ρ = −0.58, p = 0.012*.  Red points: species
with significant enrichment AND younger TE age at breakpoints.
Source data: Fig5D_te_dynamics.

**(E)** TE–Ks age correlation — panel A: Spearman correlation between anchor-pair Ks
(WGD age proxy) and TE superfamily burst mode, for species with reliable Ks estimates
(Section 14.8).  Each point: one species × superfamily; colour: TE class.  Spearman
rank correlation computed per TE superfamily across all WGD species; BH-corrected
p per superfamily listed in source data Fig5E_te_ks_corr.
Source data: Fig5E_te_ks_corr.

---

### Figure 6. Horizontal transposon transfer (HTT)

**(A)** HTT chord diagrams (Section 15).  Large panel: all 223 HTT candidate links
between EDF clades (`HTT/figures/htt_chord_v2.pdf`).  Small inset: validated HTT
events only (`HTT/figures/figHTT_chord_validated.pdf`).  Arc width ∝ number of
candidates; clade colours follow the 18-clade scheme.
Source data: Fig6A_htt_chord.

**(B)** HTT summary — panel A of `HTT/figures/figHTT_summary.pdf`: Ks distribution
of 223 HTT candidates, stacked by TE class (Section 15.1).
Source data: Fig6B_htt_ks.

**(C)** HTT summary — panel B: heat-map of HTT candidate counts per clade pair
(triangular matrix; Section 15.1).  Cell colour: log₁₀(n + 1) candidates.
Source data: Fig6C_htt_heatmap.

**(D)** HTT LINE synteny — synteny plot showing shared locus architecture (flanking
gene conservation) for validated LINE-class HTT events (Section 15.3).
Source data: Fig6D_htt_line_synteny.


### Supplementary Figures

**Figure S1. Genomic feature variation across EDF clades and pairwise synteny.**

**(A)** Box plots of mean introns per gene in multi-exon genes per clade across 18 EDF
clades (HO+AM n = 29, NU n = 2, AP n = 2, CR n = 2, MI n = 22, BL n = 5, OL n = 1,
CH n = 30, NE n = 7, BS n = 2, EN n = 5, ZO n = 4, KX n = 13, GL n = 7, MO n = 26,
MU n = 63, BA n = 14, AS n = 13; Section 3.1; source:
`GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY_V4.txt`).  Species counts shown below each
box; box = IQR, whiskers = 1.5 × IQR.
Source data: FigS1A_introns_per_clade.

**(B)** Horizontal stacked bar chart of the distribution of corrected pairwise synteny
ratios (0–1.0) per species pair, sorted by total bar length to produce a triangular
pattern; bars coloured by comparison group (Section 4.2).
Source data: FigS1B_synteny_distribution.

**(C)** Box plots of mean pairwise corrected synteny ratio per clade (0.0–1.0).
NU shows a bimodal distribution with one near-zero outlier.  EN has a median close
to 0.
Source data: FigS1C_synteny_per_clade.

---

**Figure S2. TE landscape — per-clade and per-subdivision summaries.**

**(A)** Violin plots of total TE content (% genome) per clade (Section 7; source:
`Section2_TElandscape/figure_species_data.tsv`).  Outlier species labeled:
*Jaapia flammicorona*, *Jaapia lactiflua*, *Endogone* sp., *Scutellospora mucronatum*,
*Linderina selenospora*.  Entomophthoromycotina (EN) spans the widest range (~0–95%).
Source data: FigS2A_te_content_per_clade.

**(B)** 100% stacked bar charts of TE class composition (DNA, LTR, LINE, SINE, Rolling
Circle, Penelope, Unclassified) per species within 10 EDF-only subdivisions:
Blastocladiomycota (5 spp.), Chytridiomycetes (~15 spp.), Entomophthoromycotina
(5 spp.), Glomeromycotina (~7 spp.), Kickxellomycotina (~10 spp.), Microsporidia
(~12 spp.), Mortierellomycotina (~20 spp.), Neocallimastigomycetes (~6 spp.),
Zoopagomycotina (4 spp.), Mucoromycota (30+ spp.; Section 7.1; source:
`Section2_TElandscape/figure_species_data.tsv`).
Source data: FigS2B_te_class_edf_clades.

**(C)** Box plots of Shannon TE diversity index H′ per clade (~0.5–3.5; Section 7.2;
source: `Section2_TElandscape/B01_shannon_diversity.tsv`).  EN has the lowest median
(~1.0); MU shows the widest within-clade range.
Source data: FigS2C_shannon_per_clade.

---

**Figure S3. TE content correlates with genomic architecture; domain expansion
screening; MAT locus architecture.**

**(A)** Scatter plot of genome size (Mb, log scale, x-axis) versus total TE content
(% of genome, y-axis; Section 3.2; source:
`Section2_TElandscape/figure_species_data.tsv`).  Points coloured by clade; linear
regression shown.  **Pearson's r = 0.72, p = 6.75 × 10⁻²⁵.**
Source data: FigS3A_genome_size_te_scatter.

**(B)** Pearson r heatmap between eight TE category abundances (rows: Total TE, DNA,
LTR, LINE, SINE, Rolling Circle, Penelope, Unclassified) and three genomic properties
(columns: Genome Size, Intron Density, Gene Density; Section 3.2).  Selected values:
Total TE × Genome Size = +0.76\*, Total TE × Gene Density = −0.48\*, LTR × Genome
Size = +0.69\*, LINE × Genome Size = +0.55\*.  Asterisks = BH-corrected significance
(q < 0.05; source: `GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY_V4.txt`).
Source data: FigS3B_te_genomic_correlations.

**(C)** Bubble scatter plot of domain-expanded gene families.  x-axis: gene copies in
expanded clade (log scale); y-axis: clade groupings.  Bubbles coloured by clade and
labeled with domain names including Dynein heavy chain, FAD domain, Poly-adenylate
binding, Serpin, Trypsin, SPRY domain (Section 13.4; source:
`expansion_analysis/results/clade_expansion_table.tsv`).
Source data: FigS3C_domain_expansion_bubble.

**(D)** Horizontal bar chart of expanded domains per clade, split by species-specific
expansions (orange) and cladal expansions (blue).  Neocallimastigomycetes (NE) has the
longest cladal bar (~90 domains; Section 13.4).
Source data: FigS3D_domain_expansion_counts.

**(E)** Linear synteny diagrams of the MAT locus for eight species: *Actinomucor
caerulea*, *Chaetomium thermophila*, *Spizellomyces megalocarpus*, *Rhizopus
spectabilis*, *Zygorhynchus mexicana*, *Smittium vasiformis*, *Hesseltinella pulchrum*,
*Pilobolus articulosus* (Section 5; source: `MAT_locus/results/`).  Labeled genes:
*sagA*, *sexP*, *sexM*, *rnhA*, *tnlA*.  Scale bar = 1 kb; shaded regions indicate
MAT idiomorphs.
Source data: FigS3E_mat_locus_synteny.

---

**Figure S4. Nucleosome repeat length (NRL) analysis over TE insertions across
species and clades.**

**(A)** Top: schematic of NRL measurement and Δ NRL / Δ CV(linkers) computation
for TEs versus protein-coding genes (PCGs).  Bottom: Cleveland-style lollipop dot
plots of Δ NRL (bp) and Δ CV(linkers) with 95% confidence intervals for three
exemplar species (*Neurospora crassa*, *Botrytis cinerea*, *Rhizopus microsporus*),
broken down by TE type (LINE, Unknown, LTR, DNA, PCG; Section 11).
Source data: FigS4A_nrl_exemplar_species.

**(B)** Strip plots of Δ NRL (bp) and Δ CV(linkers) across 13 clades:
Zoopagomycotina, Ascomycota, Kickxellomycotina, Entomophthoromycotina,
Chytridiomycetes, Basidiobolomycota, Glomeromycotina, Mucoromycina,
Mortierellomycotina, Basidiomycota, Microsporidia, Blastocladiomycota,
Neocallimastigomycetes.  Large filled dots = clade medians; small dots = individual
species values.
Source data: FigS4B_nrl_per_clade.

**(C)** Paired jitter plots of Δ NRL (bp) and Δ CV(linkers) broken down by TE type
(LINE, Unknown, SINE, LTR, DNA) relative to PCG baseline.  Small dots coloured by
clade (individual species); large dots = clade medians.  95% CIs shown.
Source data: FigS4C_nrl_per_te_type.

---

**Figure S5. Epigenetic regulation of TEs — supplementary analyses.**

**(A)** Box plots comparing total TE genomic content (%) between Mucoromycota species
carrying both DNMT1 and DIM2 orthologs versus all other Mucoromycota (Section 8.3;
source: `DAM_search/results/methylase_cooption_summary.tsv`).
**Wilcoxon rank-sum test: p = 4.33 × 10⁻³.**  DNMT1+DIM2 species have a higher
median TE content (~35%).
Source data: FigS5A_dnmt_te_mucoromycota.

**(B)** Violin plots with overlaid jitter of coverage-weighted 5mC methylation by TE
age class: Young (n = 6,206 elements), Mid (n = 8,924), Old (n = 29,444; Section 9.3;
source: `TE_analysis/WGBS/te_methylation_per_copy.tsv`).
**Kruskal–Wallis test: p < 2.2 × 10⁻³⁰⁰.**  Young TEs are the least methylated; Old
TEs are the most heavily methylated.
Source data: FigS5B_5mc_te_age_class.

**(C)** Five horizontal log₂ fold-enrichment bar charts showing TE class enrichment in
each epigenetic mark × species combination: *Phycomyces blakesleeanus* (5mC), *Mucor
lusitanicus* (H3K9me3), *Rhizopus microsporus* (H3K9me3), *P. blakesleeanus* (6mA),
*R. microsporus* (6mA; Section 9.5).
Source data: FigS5C_enrichment_per_species.

**(D)** Per-TE-family box plots of epigenetic enrichment across five mark × species
contexts — H3K9me3 × *M. lusitanicus* (TcMar-Tc1, CMC-Chapaev-3, RTE-RTE, L1-CRE,
Crypton-F, TcMar-Pogo); H3K9me3 × *R. microsporus* (TcMar-Tc1, Merlin, Gypsy,
RTE-BovB, Alu, PiggyBac); 5mC × *P. blakesleeanus* (LINE, Crypton-F, Gypsy, RTE-BovB,
R1, CMC-EnSpm); 6mA × *R. microsporus* (MULE-MuDR, Crypton-F, TcMar-ISRm11,
TcMar-Tc1, Gypsy, CMC-EnSpm); 6mA × *P. blakesleeanus* (Crypton-V, Zseupton,
TcMar-Sagan, RTE-BovB, MULE-MuDR, CMC-EnSpm).
Source data: FigS5D_epigenetic_per_family.

**(E)** Two-circle Venn diagram: TE loci enriched for H3K9me3 (n = 180 exclusive;
3% of the H3K9me3-enriched set) versus TE loci targeted by sRNA (n = 7,587 exclusive;
82% of the sRNA-targeted set); **overlap = 784 TE loci** (Section 9.5; source:
`TE_analysis/sRNA/`).
Source data: FigS5E_h3k9me3_srna_venn.

**(F)** Pie charts per TE superfamily showing the proportion of TE copies targeted
exclusively by sRNA versus co-targeted by both sRNA and H3K9me3.  Superfamilies (n
copies): Gypsy (n = 1,833; 41%/54%), Copia (n = 489; 42%/57%), Other LTR (n = 592;
43%/57%), PIF-Harbinger (n = 158; 33%/66%), MULE (n = 192; 38%/62%), EnSpm
(n = 307; 41%/59%), PiggyBac (n = 71; 42%/55%), hAT (n = 311; 43%/57%), Polinton
(n = 149; 44%/55%).
Source data: FigS5F_srna_h3k9_pie.

---

**Figure S6. Correlations between epigenetic mark dynamics and TE transcription.**

**(A)** 5mC context (*P. blakesleeanus*): Forest plot of Pearson r between
condition-matched changes in TE methylation (Δ5mC) and TE expression (Δ RNA) for
11 TE families; filled red circles = p < 0.05 (uncorrected); dot size ∝ n.
Four condition comparisons shown (Light/Dark Myc, Spo/Myc Dark, Light/Dark Spo,
Spo/Myc Light).  Highlighted scatters:
**TcMar-Sagan: r = 0.41, p = 0.023, n = 24; RTE-BovB: r = −0.36, p = 0.04, n = 24**
(Section 12).
Source data: FigS6A_5mc_rna_corr.

**(B)** 6mA context (*R. microsporus*): Forest plot of Pearson r for 10 TE families;
filled blue circles = p < 0.05.  Highlighted scatters:
**CMC-EnSpm: r = 0.40, p = 0.001, n = 104; Helitron: r = 0.37, p = 0.018, n = 33**
(Section 12).
Source data: FigS6B_6ma_rna_corr_rmicro.

**(C)** Additional 6mA context (Light/Dark condition only): Forest plot for 5 TE
families; all comparisons non-significant.  Shown for reference:
Gypsy r = 0.19, p = 0.066, n = 63; TcMar-Tc1 r = −0.09, p = 0.227, n = 72.
Source data: FigS6C_6ma_rna_corr_alt.

---

**Figure S7. TE co-option — supplementary analyses.**

**(A)** Box plots of the percentage of orthogroups containing ≥1 co-opted gene per
clade across 13 EDF clades (HO+AM, MI, BL, CH, NE, EN, ZO, KX, GL, MO, MU, BA, AS;
Section 13.1).  Individual species data points overlaid; clades coloured distinctly.
Source data: FigS7A_cooption_pct_per_clade.

**(B)** Paired box plots of introns per gene for co-opted (dark grey) versus
non-co-opted (light grey) genes per clade (15 clades; n values per pair shown;
y-axis 0–40+ introns/gene; Section 13.2).
Significance: ** p < 0.01, *** p < 0.001 (Wilcoxon signed-rank test, BH correction).
Source data: FigS7B_cooption_introns.

**(C)** Domain interaction network diagrams (node-link graphs) for **DUF4208**,
**DUF3535**, and **DUF4217**.  Central white node = focal DUF domain; surrounding
nodes = co-occurring partner domains (node size ∝ number of genes; colour = fraction
of clade species).  Scale bars per network: DUF4208 (1, 5, 54 genes), DUF3535 (1, 5,
57 genes), DUF4217 (1, 18, 255 genes; Section 13.3).
Source data: FigS7C_duf_interaction_networks.

**(D)** Horizontal bar chart of the number of genes carrying each ORF-derived TE
domain fusion across EDF species (~40 domain families shown; Section 13.5).
Highlighted domains: **FLYWCH**, **CotH**, **NACHT**.
Source data: FigS7D_orf_te_fusion_bar.

**(E)** Radial fan phylogenetic tree (right half) annotated with convergent fusion
events: near *Dichotomocladium decumbens* — **NACHT + Rvt1 fusion**; near *Linderina
selenospora* — **TnpTc5 + NACHT fusion** (Section 13.6).
Source data: FigS7E_nacht_fusion_phylogeny.

---

**Figure S8. Adjacent gene expression by TE family and sRNA abundance.**

Scatter plot of TE superfamily mean sRNA abundance (RPM, y-axis) versus mean
adjacent-gene expression log₂ RPKM (x-axis) per TE superfamily across all
ChimeraTE-analysed species (Section 12.2; source:
`TE_dynamics/48_adjacent_genes/48_family_summary_scatter_v2.pdf`).  Each point is one
TE superfamily; colour = TE class.  Spearman ρ and BH-corrected p shown; labelled
points: superfamilies with |ρ| > 0.5 and padj < 0.05.
Source data: FigS8_srna_adjacent_expr.

---

**Figure S9. DAM methyltransferase–RT fusion phylogeny and HMGb3 motifs in Gypsy
elements.**

**(A)** Circular maximum-likelihood phylogeny of all DAM+RT fusion proteins identified
across EDF species (Section 8.4; source:
`DAM_search/results/phylo/DAM_full.nwk` + `DAM_full_annot.tsv`).  Tips coloured by
host clade (brown = Glomeromycota, blue = Holozoa, orange = Kickxellomycotina, pink =
Mortierellomycotina, rust = Mucoromycota, yellow-green = Zoopagomycota).  Filled
circles on branches: **dark red = DIRS-class origin; pink/rose = Gypsy-class origin.**
Species with multiple copies labeled (*R. oryzae*, *H. radiata*, *M. mucedo* ×2,
*R. microsporus* ×4, *D. discoideum*, *G. rosea*).
Source data: FigS9A_dam_rt_phylo.

**(B)** Bifurcating cladogram of the Gypsy RT chromovirus clade (Section 8.4).
Branches coloured by sub-group: non-chromo-outgroup (blue), Chlamyvir (green),
chromo-outgroup (red), Reina (purple), unclassified (grey).  Three columns of
presence/absence squares on the right indicate **HLH1/Myc-B**, **HLH1/Myc-S**, and
**HLH1/Myc-11** homologs per element.
Source data: FigS9B_gypsy_rt_cladogram.

**(C)** Horizontal bar chart of structural lengths (0–7,000 bp; green bars) of
individual Gypsy TE elements.  **Orange dots** mark the position of the HMGb3-binding
motif within each element; motifs occur at consistent internal positions across elements
(Section 8.4).
Source data: FigS9C_gypsy_hmgb3_structure.

---

**Figure S10. WGD synteny breakpoint TE enrichment — species-level details.**

**(A)** Scatter plot of total TE content (%, x-axis) versus TE enrichment at WGD
synteny breakpoints (log₂ fold-change, y-axis) per WGD species (Section 14.3).
Dashed regression line; **regression coefficient β = 0.20.**  Reference lines at
y = 0 (no enrichment) and y = 1 (2-fold enrichment).  Labeled data points:
*R. delemar* (Rd), *F. linderi* (Fr), *P. articulosus* (Pa), *C. cucurbitarum* (Cc),
*C. brefeldii* (Cb).
Source data: FigS10A_te_vs_bp_enrich.

**(B)** Volcano plots of TE superfamily enrichment at synteny breakpoints (log₂FC,
x-axis; −log₁₀ adjusted p-value, y-axis; red dashed significance threshold) for five
species individually (Section 14.3).  Significant superfamilies labeled:
*C. brefeldii*: Zeuplon (−log₁₀p ≈ 10), DIRS (~4), TcMar-ISRm11, hAT-Charlie;
*C. cucurbitarum*: Ginger-2 (~6), Penelope (~3); *F. linderi*: none;
*P. articulosus*: none; *R. delemar*: L1-Tx1 (~6), TcMar-ISRm11 (~4), MULE-MuDR
(~3), TcMar-Pogo (~3).
Source data: FigS10B_bp_enrichment_volcano.

**(C)** Horizontal bar chart of the percentage of syntenic multiplicons containing each
significant TE superfamily, per species: *C. brefeldii*: Zeuplon (~22%), DIRS (~5%),
hAT-Charlie (~45%), TcMar-ISRm11 (~45%); *C. cucurbitarum*: Ginger-2 (~17%), Penelope
(~15%); *R. delemar*: TcMar-Pogo (~8%), MULE-MuDR (~17%), L1-Tx1 (~47%),
TcMar-ISRm11 (~60%).
Source data: FigS10C_bp_multiplicon_pct.

---

**Figure S11. Horizontal transposon transfer — TE class proportions and biases.**

**(A)** 100% stacked bar chart of the proportion of HTT events by TE class (DNA, LTR,
LINE, SINE, Rolling Circle) per clade.  Species counts per clade: HO n = 16, NU n = 1,
AP n = 17, BL n = 19, CH n = 63, BS n = 60, ZD n = 14, EN n = 14, KX n = 16,
MD n = 83, GL n = 3, MU n = 116, AS n = 37, BA n = 11 (Section 15).
Source data: FigS11A_htt_te_class_per_clade.

**(B)** Scatter plot of genomic TE class abundance (% of all TE, x-axis, 0–70%)
versus fraction of HTT events (%, y-axis, 0–70%) per TE class.  Dashed diagonal =
null expectation.  **LINE** (teal, ~65% of HTTs vs. ~25% genomic) is strongly
overrepresented above the diagonal.  **LTR** (~12% HTTs, ~32% genomic) and **DNA**
(~8% HTTs, ~30% genomic) are underrepresented.  **RC/Rolling Circle** (~10% HTTs,
~3% genomic) is overrepresented.  **SINE** (~5% HTTs, ~5% genomic) is near the
diagonal (Section 15).
Source data: FigS11B_htt_te_class_scatter.

---

## Part III — Source Data Index

The Excel workbook `EDF-TE_source_data.xlsx` contains one sheet per figure panel.
Each sheet provides the numerical data underlying one or more figure panels as
described below.

### Main figures

| Sheet name | Figure panel(s) | Description | Key columns | Primary data source |
|---|---|---|---|---|
| Fig1A_tree | Fig. 1A | Species tree tip metadata | species, clade, clade_colour, bootstrap_support | `SpeciesTree/` newick + clade assignments |
| Fig1B_pangenome | Fig. 1B | Pan-genome accumulation curve data (mean ± s.d. per group at each species count) | group_code, group_name, n_species, pangenome, core, shell, unique | `PanGenome/results/summary_final_counts.tsv`; per-group curves in `PanGenome/results/per_group/{group}_pangenome.tsv` (columns: n, pangenome_mean, pangenome_sd, core_mean, core_sd, shell_mean, shell_sd, unique_mean, unique_sd) |
| Fig1C_genomic_features | Fig. 1C | Per-species genomic feature metrics | species, clade, Genome_Size_bp, Total_Genes, Total_Introns, Avg_Exons_Per_Gene, Avg_Introns_All, Avg_Isoforms_Per_Gene | `GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY_V4.txt` |
| Fig1D_synteny | Fig. 1D | Corrected synteny ratio per species pair | species_A, species_B, clade, observed_collinear_pairs, null_mean_collinear, total_adjacent_slots, corrected_score, normalised_score | `Synteny/results/synteny_summary.tsv` |
| Fig2A_te_content | Fig. 2A | TE class composition + Shannon H′ + Recency Index per species | taxon, clade, pct_realTE, pct_ClassI, pct_ClassII, pct_TE_total, pct_unknown | `Section2_TElandscape/figure_species_data.tsv`; Shannon H′ and RI from `B01_shannon_diversity.tsv` and `C01_kimura_species.tsv` |
| Fig2B_pca_composition | Fig. 2B | PCoA of TE high-level class composition (Bray–Curtis) | species, clade, PCo1, PCo2, PCo3 | `Section2_TElandscape/G01_pcoa_coords.tsv` |
| Fig2C_kimura_ri | Fig. 2C | Recency Index per species and Kruskal–Wallis + Dunn test results | species, clade, recency_index, burst_score, mean_div_weighted, n_superfamilies, total_te_bp | `Section2_TElandscape/C01_kimura_species.tsv` |
| Fig2D_kimura_summary | Fig. 2D | Mean Kimura divergence per TE superfamily × species | species, superfamily, total_bp, pct_of_TE, recency_index, mean_div_weighted, n_families | `Section2_TElandscape/C01_kimura_superfamily.tsv` |
| Fig2E_proximity_scatter | Fig. 2E | Gene–TE proximity scatter (41 candidate domains) | pfam, domain_name, te_class, mw_pval_copynumber, n_pack_element_candidates, n_fusion_proteins, n_te_consensus_with_domain, mechanism | `expansion_analysis/results/figures/labeled_domains_metrics.tsv` |
| Fig2F_proximity_flanks | Fig. 2F | TE subfamily in symmetric flanks of domain-expanded genes | species, clade, pfam, domain_name, n_copies, pct_nearest_is_te, bg_pct_nearest_is_te, mw_pval_vs_bg | `expansion_analysis/results/proximity/proximity_summary.tsv` |
| Fig2G_sex_locus | Fig. 2G | Sex locus MAT architecture and TE composition | species, clade, mat_locus_genes, te_superfamily, te_pct | `MAT_locus/figures/sex_locus_figure.pdf` (numerical data in `MAT_locus/results/`) |
| Fig2H_sex_locus_synteny | Fig. 2H | Sex locus microsynteny between *Fenestella* sp. and *Phycomyces blakesleeanus* | gene_A, gene_B, species_A, species_B, link_colour, scaffold_A, scaffold_B | `Synteny/results/sex_locus_synteny_links.tsv` |
| Fig3A_dnmt | Fig. 3A | Methyltransferase domain screening across TE ORFs + per-species co-option summary | domain, pfam, n_species_with_domain, n_species_TE_context, n_species_fused_TE, TE_class_involved | `DAM_search/results/methylase_cooption_summary.tsv`; TE-embedded DNMT screening: `Co-option/te_dnmt_search/te_with_dnmt.tsv` (columns: te_id, has_C5_DNMT, has_TE_core, pfam_domains) |
| Fig3B_5mc_family_age | Fig. 3B | Coverage-weighted mean 5mC per TE family × Kimura age bin | species, family, age_bin_kimura, wm_5mC, n_copies_covered | `TE_dynamics/12_te_age_5mc/te_wm_per_copy.tsv` or `20_chiapello_extension/data/F3_anova_wm.tsv` |
| Fig3C_gene_dist_5mc | Fig. 3C | Gene–TE distance distribution by 5mC methylation quartile; KS D and BH-p for Q1 vs Q4 | species, distance_bp, quartile, n_copies | `TE_dynamics/31_meth_bins_gene_expr/wm_per_TE.tsv` + `31_meth_bins_gene_expr/te_bin_adjgene_expr.tsv` |
| Fig3D_adj_expr_meth | Fig. 3D | Adjacent gene expression (log₂ RPKM) by methylation quartile; BH-corrected Mann–Whitney U | species, gene_id, te_id, quartile, log2_rpkm | `TE_dynamics/31_meth_bins_gene_expr/te_bin_adjgene_expr.tsv` |
| Fig3E_gene_dist_h3k9 | Fig. 3E | Gene–TE distance by H3K9me3 enrichment class; KS D and BH-p for Q1 vs Q4 | species, distance_bp, h3k9_class, n_copies | `TE_dynamics/33_h3k9_bins_gene_expr/te_h3k9_adjgene_expr.tsv` |
| Fig3F_adj_expr_h3k9 | Fig. 3F | Adjacent gene expression by H3K9me3 peak presence/absence; BH-corrected Mann–Whitney U | species, gene_id, te_id, h3k9_peak_present, log2_rpkm | `TE_dynamics/33_h3k9_bins_gene_expr/te_h3k9_adjgene_expr.tsv` |
| Fig3G_mac_tss | Fig. 3G | TE expression grouped by MAC at TSS; Mann–Whitney U BH-p | species, condition, te_id, mac_at_tss, log2_rpkm | `TE_dynamics/45_mac_tss_expression/45_wilcoxon_results.tsv` |
| Fig3H_metaplot | Fig. 3H | Per-bin mean 5mC across TE bodies (Q1 vs Q4) | superfamily, bin_index, quartile, mean_5mC, sem_5mC, n_copies | `TE_dynamics/F1_meth_q1q4_metaplot/F1_meth_q1q4_profile_data.tsv` or `F1_meth_q1q4_binned.tsv` |
| Fig3I_mac_tss_families | Fig. 3I | TE families enriched at methylated TSS; hypergeometric q and enrichment ratio | species, superfamily, n_mac_tss_overlaps, enrichment_ratio, BH_q | `TE_dynamics/44b_mac_tss_families/44b_summary.tsv` |
| Fig3J_rnai_pies | Fig. 3J | RNAi pathway targeting proportions per TE superfamily | species, te_superfamily, rnai_status, proportion | `TE_dynamics/46_rnai_snapshot/` RNAi overview data |
| Fig3K_srna | Fig. 3K | sRNA RPM vs. TE expression per superfamily | species, superfamily, srna_rpm, te_rpkm, srna_status | `TE_dynamics/35_srna_mucor/plots/35g_coregulated_superfamily_stats.tsv` |
| Fig4A_cooption_tree | Fig. 4A | Co-option gain/loss per node and total orthogroups per clade | node_id, clade, n_ancestral_retained, n_newly_acquired, n_total_og_tip | `Co-option/figures/14_cooption_traits.tsv` |
| Fig4B_domain_barchart | Fig. 4B | Host-function Pfam domain composition of co-opted genes per clade | clade, pfam_domain, domain_desc, n_coopted_genes, fraction | `expansion_analysis/results/clade_expansion_table.tsv` |
| Fig4C_te_class_prop | Fig. 4C | TE class origin of co-opted genes per clade | clade, te_class, n_coopted, fraction | `Co-option/per_species/classification/*.true_TE_coopted_genes.tsv` (aggregated) |
| Fig4D_fusion_network | Fig. 4D | Fusion protein network edge list (host domain × TE domain co-occurrence) | te_domain_ipr, host_domain_ipr, n_proteins, n_species, species_list | `expansion_analysis/results/figures/Supplementary_Table_S_fusion_combos.tsv` |
| Fig4E_chimte | Fig. 4E | ChimeraTE chimeric gene–TE junction counts per species | species, clade, te_class, n_chimeric_junctions | `Co-option/chimte_summary.tsv` |
| Fig4F_flywch_trees | Fig. 4F | FLYWCH × TE-RT fusion protein phylogeny metadata | protein_id, species, clade, te_superfamily, fusion_type, n_fusion_proteins | `expansion_analysis/results/figures/labeled_domains_metrics.tsv` (PF04500 rows) + `expansion_analysis/results/phylo/PF04500*.nwk` |
| Fig4G_coth_trees | Fig. 4G | CotH × TE fusion protein phylogeny metadata | protein_id, species, clade, te_superfamily, fusion_type, n_fusion_proteins | `expansion_analysis/results/figures/labeled_domains_metrics.tsv` (PF08757 rows) + `expansion_analysis/results/phylo/CotH*.nwk` |
| Fig4H_dam_rt_tree | Fig. 4H | DAM methyltransferase + RT protein metadata for circular phylogeny | species, clade, te_class, source, DAM_evalue, RT_evalue, rep_orf, orf_length_aa, domain_architecture | `DAM_search/results/DAM_TE_summary_table_annotated.tsv`; tree + annotation: `DAM_search/results/phylo/DAM_full.nwk` + `DAM_full_annot.tsv` |
| Fig4I_tf_enrichment | Fig. 4I | TF peak enrichment in TE classes (DAP-seq, 58 TFs); Fisher's exact test BH-q | te_class, motif_id, motif_family, log2_odds_ratio, fisher_pval, BH_q | `TF_TE_analysis/03_fimo_enrichment/fimo_enrichment.tsv` |
| Fig4J_hmgb3_metaplot | Fig. 4J | HMGb3 DAP-seq signal metaplot over LINE/I-Jockey; permutation p | bin_index, mean_signal, sem_signal, te_class | `TF_TE_analysis/02_bigwig_signal/bw_TE_signal_class_mean.tsv` (HMGb_1228 rows) |
| Fig4L_hmgb3_heatmap | Fig. 4L | HMGb3 occupancy across Jockey clade × subfamily | clade, te_subfamily, mean_norm_signal | `TF_TE_analysis/02_bigwig_signal/` |
| Fig5A_paralogs | Fig. 5A | Gene duplication rates and co-linear paralog counts per WGD-candidate species | species, node_id, dup_rate_per_eligible_tree, n_colinear_paralog_pairs, wgd_classification, ks_median | `WGD/Modern_wgd2_pipeline/results/paper_no_filter/duplication_rates_by_node_normalized.tsv` + `synteny_summary.tsv` |
| Fig5B_syntenic_blocks | Fig. 5B | Syntenic block structure and TE content per WGD species | species, block_id, block_length_bp, te_bp, te_fraction | `WGD/Modern_wgd2_pipeline/results/paper_no_filter/synteny_summary.tsv` |
| Fig5C_te_breakpoints | Fig. 5C | TE enrichment at WGD synteny breakpoints; Wilcoxon rank-sum BH-q | species, superfamily, log2_fc, wilcoxon_p, BH_q, obs_te_bp_per_window, bg_te_bp_per_window | `WGD/Modern_wgd2_pipeline/results/paper_no_filter/te_breakpoint_enrichment/` per-species TSVs |
| Fig5D_te_age_breakpoints | Fig. 5D | TE insertion age (ΔKimura80) vs. enrichment at WGD breakpoints; Spearman ρ = −0.58, p = 0.012 | species, delta_kimura80, log2_fc_breakpoint, spearman_rho, spearman_p | `WGD/Modern_wgd2_pipeline/results/paper_no_filter/te_breakpoint_enrichment/` + Kimura data |
| Fig5E_te_ks_corr | Fig. 5E | TE superfamily burst mode vs. anchor-Ks WGD age; Spearman ρ per superfamily | species, superfamily, te_burst_ks_mode, anchor_ks_median, spearman_rho, BH_padj | `WGD/Modern_wgd2_pipeline/results/paper_no_filter/ks_analysis/` anchor Ks files + `C01_kimura_superfamily.tsv` |
| Fig6A_htt_chord | Fig. 6A | HTT validated-event counts per clade pair (chord diagram) | clade_A, clade_B, n_candidates, n_validated, arc_width | `HTT/06_htt_candidates/htt_candidates/htt_by_clade_pair.tsv` |
| Fig6B_htt_ks | Fig. 6B | Ks distribution of HTT candidates by TE class | pair_id, qclade, sclade, te_class, te_ks, is_htt | `HTT/06_htt_candidates/htt_candidates/htt_candidates.tsv` |
| Fig6C_htt_heatmap | Fig. 6C | HTT candidate count heat-map per clade pair (triangular matrix) | qclade, sclade, htt_count | `HTT/06_htt_candidates/htt_candidates/htt_by_clade_pair.tsv` |
| Fig6D_htt_line_synteny | Fig. 6D | LINE HTT locus flanking gene conservation | event_id, species_A, species_B, flank_gene_A, flank_gene_B, conservation_status | `HTT/09_minimal_htt_events/events/` per-component TSVs |

### Supplementary figures

| Sheet name | Figure panel(s) | Description | Key columns | Primary data source |
|---|---|---|---|---|
| FigS1A_introns_per_clade | Fig. S1A | Mean introns per multi-exon gene per species, grouped by clade | species, clade, Avg_Introns_All, Total_Genes, Total_Introns | `GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY_V4.txt` |
| FigS1B_synteny_distribution | Fig. S1B | Distribution of corrected pairwise synteny ratios per species pair | species_A, species_B, comparison_group, corrected_score | `Synteny/results/synteny_summary.tsv` |
| FigS1C_synteny_per_clade | Fig. S1C | Mean pairwise corrected synteny ratio per clade | species_A, species_B, clade, corrected_score, mean_clade_synteny | `Synteny/results/synteny_summary.tsv` |
| FigS2A_te_content_per_clade | Fig. S2A | Total TE content (%) per species grouped by clade | taxon, clade, pct_realTE, pct_TE_total | `Section2_TElandscape/figure_species_data.tsv` |
| FigS2B_te_class_edf_clades | Fig. S2B | TE class composition stacked bars per species within each EDF-only subdivision | taxon, clade, subdivision, te_class, pct_of_te_content | `Section2_TElandscape/figure_species_data.tsv` (EDF-subdivision filter) |
| FigS2C_shannon_per_clade | Fig. S2C | Shannon H′ TE diversity per species grouped by clade | species, clade, H_prime, J_evenness, n_superfamilies, dominant_superfamily | `Section2_TElandscape/B01_shannon_diversity.tsv` |
| FigS3A_genome_size_te_scatter | Fig. S3A | Genome size (Mb) vs. total TE content (%); Pearson r = 0.72, p = 6.75 × 10⁻²⁵ | taxon, clade, genome_size_mb, pct_realTE | `Section2_TElandscape/figure_species_data.tsv` + `GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY_V4.txt` |
| FigS3B_te_genomic_correlations | Fig. S3B | Pearson r heatmap: 8 TE classes × 3 genomic features; BH-corrected q | te_class, genomic_feature, pearson_r, p_raw, BH_q | `GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY_V4.txt` (computed pairwise) |
| FigS3C_domain_expansion_bubble | Fig. S3C | Domain-expanded gene copies per clade bubble scatter | pfam, domain_name, clade, n_copies_expanded, max_copies_species | `expansion_analysis/results/clade_expansion_table.tsv` |
| FigS3D_domain_expansion_counts | Fig. S3D | Count of species-specific and cladal expansion domains per clade | clade, n_species_specific_domains, n_cladal_domains | `expansion_analysis/results/clade_expansion_table.tsv` |
| FigS3E_mat_locus_synteny | Fig. S3E | Linear synteny of MAT locus in 8 EDF species | species, gene_name, start_bp, end_bp, strand, mat_region | `MAT_locus/results/` per-species GFF + synteny links |
| FigS4A_nrl_exemplar_species | Fig. S4A | Δ NRL (bp) and Δ CV(linkers) for N. crassa, B. cinerea, R. microsporus by TE type; 95% CIs | species, te_type, delta_NRL, delta_CV, CI_lower, CI_upper | `TE_analysis/MNase/nrl_per_species_te_type.tsv` |
| FigS4B_nrl_per_clade | Fig. S4B | Δ NRL and Δ CV(linkers) per clade across 13 clades; species + clade medians | species, clade, te_type, delta_NRL, delta_CV | `TE_analysis/MNase/nrl_per_species_te_type.tsv` (aggregated by clade) |
| FigS4C_nrl_per_te_type | Fig. S4C | Δ NRL and Δ CV(linkers) by TE type (LINE, Unknown, SINE, LTR, DNA) vs. PCG baseline | te_type, species, clade, delta_NRL, delta_CV | `TE_analysis/MNase/nrl_per_species_te_type.tsv` |
| FigS5A_dnmt_te_mucoromycota | Fig. S5A | TE content comparison between DNMT1+DIM2 and other Mucoromycota; p = 4.33 × 10⁻³ | species, has_dnmt1_dim2, total_te_pct | `DAM_search/results/methylase_cooption_summary.tsv` (Mucoromycota filter) |
| FigS5B_5mc_te_age_class | Fig. S5B | Coverage-weighted 5mC by TE age class (Young/Mid/Old); Kruskal–Wallis p < 2.2 × 10⁻³⁰⁰ | te_copy_id, age_class, n_elements, weighted_5mC | `TE_analysis/WGBS/te_methylation_per_copy.tsv` (age-class annotation) |
| FigS5C_enrichment_per_species | Fig. S5C | Log₂ TE class enrichment in 5mC/H3K9me3/6mA per species | species, mark, te_class, log2_fc_enrichment | Per-species enrichment tables in `TE_analysis/` subdirectories |
| FigS5D_epigenetic_per_family | Fig. S5D | Epigenetic mark levels per TE family for 5 mark × species contexts | species, mark, te_family, mean_enrichment, n_copies | `TE_analysis/WGBS/`, `TE_analysis/ChIP/`, `TE_analysis/N6mA/` per-family summaries |
| FigS5E_h3k9me3_srna_venn | Fig. S5E | Venn diagram: H3K9me3-enriched vs. sRNA-targeted TE loci; overlap = 784 | locus_id, h3k9me3_enriched, srna_targeted, n_exclusive_h3k9, n_exclusive_srna, n_overlap | `TE_analysis/ChIP/` + `TE_analysis/sRNA/` overlap table |
| FigS5F_srna_h3k9_pie | Fig. S5F | Proportion of TE copies co-targeted by sRNA+H3K9me3 vs. sRNA-only, per superfamily | superfamily, n_copies, pct_srna_only, pct_srna_h3k9 | `TE_analysis/sRNA/srna_h3k9_overlap_by_superfamily.tsv` |
| FigS6A_5mc_rna_corr | Fig. S6A | Pearson r (Δ5mC vs. Δ RNA) per TE family, 5mC context; TcMar-Sagan r = 0.41 p = 0.023 n = 24; RTE-BovB r = −0.36 p = 0.04 n = 24 | te_family, pearson_r, p_val, n_comparisons, condition | `TE_analysis/WGBS/condition_meth_rna_corr.tsv` |
| FigS6B_6ma_rna_corr_rmicro | Fig. S6B | Pearson r (Δ6mA vs. Δ RNA) per TE family, R. microsporus; CMC-EnSpm r = 0.40 p = 0.001 n = 104; Helitron r = 0.37 p = 0.018 n = 33 | te_family, pearson_r, p_val, n_comparisons, condition | `TE_analysis/N6mA/condition_6ma_rna_corr_rmicro.tsv` |
| FigS6C_6ma_rna_corr_alt | Fig. S6C | Pearson r (Δ6mA vs. Δ RNA) per TE family, additional context; all ns | te_family, pearson_r, p_val, n_comparisons | `TE_analysis/N6mA/condition_6ma_rna_corr_alt.tsv` |
| FigS7A_cooption_pct_per_clade | Fig. S7A | Percentage orthogroups with ≥1 co-opted gene per clade | clade, species, pct_ogs_coopted, n_ogs_total, n_ogs_coopted | `Co-option/per_species/classification/*.true_TE_coopted_genes.tsv` (aggregated) |
| FigS7B_cooption_introns | Fig. S7B | Introns per gene: co-opted vs. non-co-opted genes; Wilcoxon ** p < 0.01, *** p < 0.001 | clade, species, gene_id, cooption_status, n_introns | `Co-option/per_species/classification/` + gene annotation |
| FigS7C_duf_interaction_networks | Fig. S7C | Domain interaction networks for DUF4208, DUF3535, DUF4217; n genes per node | focal_domain, partner_domain, n_genes, clade_fraction | `expansion_analysis/results/domain_networks/` |
| FigS7D_orf_te_fusion_bar | Fig. S7D | Gene count with ORF-derived TE domain fusions per domain family (~40 domains) | domain_family, n_genes_with_fusion | `expansion_analysis/results/figures/Supplementary_Table_S_fusion_combos.tsv` |
| FigS7E_nacht_fusion_phylogeny | Fig. S7E | Phylogenetic tree metadata for NACHT+Rvt1 and TnpTc5+NACHT fusion events | protein_id, species, fusion_type, te_class, branch_annotation | `expansion_analysis/results/phylo/NACHT*.nwk` + annotation TSV |
| FigS8_srna_adjacent_expr | Fig. S8 | TE superfamily sRNA RPM vs. adjacent-gene expression; Spearman ρ | species, superfamily, mean_srna_rpm, mean_adj_log2_rpkm, spearman_rho, BH_padj | `TE_analysis/sRNA/srna_vs_adjacent_gene_family_summary.tsv` |
| FigS9A_dam_rt_phylo | Fig. S9A | DAM+RT fusion protein circular phylogeny metadata | protein_id, species, clade, te_class, DAM_evalue, RT_evalue, orf_length_aa | `DAM_search/results/phylo/DAM_full.nwk` + `DAM_full_annot.tsv` |
| FigS9B_gypsy_rt_cladogram | Fig. S9B | Gypsy RT chromovirus cladogram with HLH/Myc-B/S/11 presence/absence | element_id, gypsy_subgroup, HLH_MycB_present, HLH_MycS_present, HLH_Myc11_present | `TF_TE_analysis/` Gypsy RT sub-classification tables |
| FigS9C_gypsy_hmgb3_structure | Fig. S9C | Gypsy TE element structural lengths and HMGb3 motif positions | element_id, element_length_bp, hmgb3_motif_pos_bp | `TF_TE_analysis/02_bigwig_signal/` + HMGb3 motif scan |
| FigS10A_te_vs_bp_enrich | Fig. S10A | TE content (%) vs. TE enrichment at WGD breakpoints (log₂FC); β = 0.20 | species, total_te_pct, log2fc_at_bp, species_code | `WGD/Modern_wgd2_pipeline/results/paper_no_filter/te_breakpoint_enrichment/` |
| FigS10B_bp_enrichment_volcano | Fig. S10B | Volcano plots per species: TE superfamily log₂FC vs. −log₁₀ adj p at WGD breakpoints | species, superfamily, log2fc, neg_log10_padj, significant | `WGD/Modern_wgd2_pipeline/results/paper_no_filter/te_breakpoint_enrichment/` per-species TSVs |
| FigS10C_bp_multiplicon_pct | Fig. S10C | Percentage of syntenic multiplicons containing each significant TE superfamily per species | species, superfamily, pct_multiplicons | `WGD/Modern_wgd2_pipeline/results/paper_no_filter/synteny_summary.tsv` |
| FigS11A_htt_te_class_per_clade | Fig. S11A | Proportion of HTT events by TE class per clade; n species per clade | clade, te_class, n_htt, proportion, n_species | `HTT/06_htt_candidates/htt_candidates/htt_by_clade_pair.tsv` |
| FigS11B_htt_te_class_scatter | Fig. S11B | HTT fraction (%) vs. genomic TE abundance (%) by TE class; LINE overrepresented | te_class, htt_fraction_pct, genomic_pct | `HTT/06_htt_candidates/htt_candidates/htt_candidates.tsv` + `Section2_TElandscape/figure_species_data.tsv` |

All numerical values are the direct inputs to the corresponding figure panels.
Processed intermediate tables are available in the per-section `results/`
directories of the code repository; paths and formats are defined in `config.sh`.

---

## Part IV — Published Data Accession IDs

The table below lists all publicly available datasets downloaded and reanalysed in
this study.  Accession numbers refer to the NCBI Sequence Read Archive (SRA) unless
otherwise stated.

### Nucleosome positioning and chromatin accessibility

| Species | Strain / condition | Assay | Accession(s) | Reference |
|---|---|---|---|---|
| *Rhizopus microsporus* | ATCC11559, vegetative | MNase-seq | Legacy WIG (see §10 / §11.4) | Lax et al. (this study) |
| *Neurospora crassa* | OR74A, vegetative | MNase-seq | SRR13067301, SRR13067302, SRR13067303 | Gurrea-Rubio et al. 2021 |
| *Botrytis cinerea* | B05.10, vegetative | MAINE-seq | SRR13224644, SRR13224647, SRR13224648 | Bellacosa et al. 2022 |
| *Microbotryum eburneum* | STE-U 5482, vegetative | ATAC-seq | SRR28069093, SRR28069094, SRR28069106 | Cascales et al. 2024 |
| *M. eburneum* | STE-U 5482, vegetative | NucSeq (in vitro) | SRR28069066, SRR28069068 | Cascales et al. 2024 |
| *Rhodotorula toruloides* | IFO 0880 | NucSeq (in vitro) | SRR28069068 | Cascales et al. 2024 |
| *Catenaria anguillulae* | PL171 | NucSeq (in vitro) | SRR7476903 | Mondo et al. 2017 (JGI) |
| *Ustilago scitaminea* / *Sf1*-related | — | NucSeq (in vitro) | SRR7472154 | Mondo et al. 2017 (JGI) |
| *Marasmius fiardii* / *Bf1* | — | NucSeq (in vitro) | SRR7472065 | Mondo et al. 2017 (JGI) |

### Histone modification ChIP-seq

| Species | Strain / condition | Assay | Accession(s) | Reference |
|---|---|---|---|---|
| *Rhizopus microsporus* | ATCC11559, vegetative | H3K9me3 ChIP + Input | In-house PE Illumina — deposited in this study | Lax et al. (this study) |
| *Mucor lusitanicus* | CBS 277.49, vegetative | H3K9me3 | (see §10; deposited in this study) | Lax et al. (this study) |

### 5-methylcytosine bisulfite sequencing (WGBS)

| Species | Strain / condition | Accession(s) | Reference |
|---|---|---|---|
| *Rhizopus microsporus* | ATCC11559 | Deposited in this study | Lax et al. (this study) |
| *Phycomyces blakesleeanus* | NRRL 1555 | Deposited in this study | Lax et al. (this study) |
| *Mucor lusitanicus* | CBS 277.49 | Deposited in this study | Lax et al. (this study) |

### N6-methyladenine (PacBio SMRT)

| Species | Strain | Accession(s) | Reference |
|---|---|---|---|
| *Phycomyces blakesleeanus* | NRRL 1555 | Deposited in this study | Lax et al. (this study) |
| *Rhizopus microsporus* | ATCC11559 | Deposited in this study | Lax et al. (this study) |

### Small RNA (sRNA / piRNA / siRNA)

| Species | Strain / condition | Accession(s) | Reference |
|---|---|---|---|
| *Mucor lusitanicus* | CBS 277.49, vegetative | Deposited in this study | Lax et al. (this study) |

### mRNA-seq for chimeric TE transcript detection (ChimeraTE)

| Species | Clade | Accession(s) | Reference |
|---|---|---|---|
| *Mucor lusitanicus* | MU | Deposited / reanalysed | Lax et al. (this study) |
| *Rhizopus microsporus* | MU | Deposited / reanalysed | Lax et al. (this study) |
| *Phycomyces blakesleeanus* | MU | Deposited / reanalysed | Lax et al. (this study) |
| *Aphanomyces* sp. (AP_Amoprot) | AP | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Blastocladiella* sp. (BL_Blabri1) | BL | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Basidiobolus* sp. (BS_Basme92) | BS | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Batrachochytrium dendrobatidis* JEL423 | CH | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Batrachochytrium salamandrivorans* | CH | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Piromyces finnis* | CH | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Enterocytozoon bieneusi* H348 | EN | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Zoophthora radicans* | EN | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Geosiphon pyriformis* | GL | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Rhizophagus clarus* | GL | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Rhizophagus irregularis* DAOM197198 | GL | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Coemansia mojavensis* | KX | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Kickxella alabastrina* | KX | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Encephalitozoon cuniculi* | MI | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Spraguea lophii* | MI | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Vavraia culicis* | MI | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Diaphanoeca grandis* (outgroup) | MO | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Monosiga brevicollis* | MO | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Nuclearia* sp. (NU_Fontalba) | NU | SRA (EDF-TE dataset) | Lax et al. (this study) |
| *Pigoraptor* sp. (ZO_Pipcy3) | ZO | SRA (EDF-TE dataset) | Lax et al. (this study) |

*Note: Exact SRA accession numbers for the EDF-TE-specific RNA-seq libraries will be added upon data deposition prior to publication.*

