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
Z-scores were computed per metric across all 320 species (mean subtraction, division
by standard deviation with ddof = 1) for the heatmap display (Figure 1B).

Pan-genome accumulation curves were generated by 1,000 random permutations of
species-addition order applied to the OrthoFinder gene presence/absence matrix
(`section1_genome_phylogenetics/03_pangenome_analysis.py`).  At each depth (number
of species sampled), orthogroups were classified as core (present in all sampled
species), shell (present in 2 to n−1 sampled species), or unique/cloud (present in
exactly one sampled species).  The mean and standard deviation across 1,000 random
orderings were calculated per depth for each category.

---

### 4. Macrosynteny analysis

#### 4.1 Pairwise synteny detection with SynChro

Pairwise macrosyntenic blocks were identified among all 38 Mucoromycota genomes
(sub-clades MU, MO and GL) using SynChro v2.0 (Drillon et al. 2014, Genome Research).
SynChro's key methodological distinction from graph-based tools (e.g. i-ADHoRe,
DAGchainer) is its phylogenetic distance correction: for every ordered species pair
(A, B) the tool extracts the tip-to-tip patristic distance d(A, B) (in expected
substitutions per site) directly from the input Newick species tree (the LG+C60+F+G4
ML tree from Section 2).  This distance is used to calibrate the minimum BLASTP
bitscore threshold for declaring a protein pair an anchor, so that syntenic regions
are detected with statistically equivalent sensitivity across genome pairs spanning
different divergence times.  Without this correction, a fixed bitscore cutoff would
systematically under-detect collinearity between divergent genomes (expected pairwise
identity < 40 %) and produce an excess of false anchor pairs between closely related
genomes (expected identity > 85 %).  In practice, the SynChro correction scales the
expected score from the BLOSUM62 substitution matrix by d(A, B) to produce a
pair-specific threshold T(d), which is then applied as the minimum bitscore for RBH
pair acceptance.

#### 4.2 Protein-sequence comparison

All-against-all pairwise protein comparisons were performed with DIAMOND BLASTP
v2.1.8.  For each ordered genome pair (A→B and B→A), DIAMOND was run as:

```
diamond blastp \
  -q <proteome_A.faa> -d <proteome_B.dmnd> \
  --more-sensitive --comp-based-stats 0 \
  -p 16 -k 5 -e 1e-5 --id 30 \
  --query-cover 50 --subject-cover 50 \
  -f 6 qseqid sseqid pident length mismatch gapopen \
       qstart qend sstart send evalue bitscore
```

`--comp-based-stats 0` disables composition-based statistics (which can reduce
sensitivity for divergent comparisons); `-k 5` retains at most five alignments per
query.  Reciprocal best hits (RBH) were defined as protein pairs (proteinA, proteinB)
where proteinB was the top DIAMOND hit of proteinA in genome B, and proteinA was the
top DIAMOND hit of proteinB in genome A (by bitscore).  RBH pairs were then filtered
by the SynChro per-pair threshold T(d) as described above.

#### 4.3 Collinear block assembly

SynChro collinear block detection was run with:
- `minimalNumberOfAnchors = 5` (minimum anchor pairs per syntenic segment)
- `maximalGap = 25` (maximum number of non-anchor genes permissible between two
  consecutive anchors within a segment)
- Strand-consistency required: forward and reverse collinear segments are tracked
  independently and are not merged.

Scaffold edges and internal runs of N ≥ 100 consecutive ambiguous bases were treated
as assembly boundaries; anchor pairs bridging such gaps were removed before block
assembly.  The output `.syn` files record the genomic coordinates (scaffold, start,
end, gene count, anchor count) of every collinear segment in both genomes, together
with the list of anchor gene pairs.

#### 4.4 Chord diagram visualisation

Pairwise anchor pair counts (summed across all collinear segments per species pair)
were visualised as a circular chord diagram.  Arc width on the chord is proportional
to log₁₀(total anchor pairs + 1); species arcs on the outer ring are ordered
phylogenetically and coloured by clade following the canonical 18-clade colour scheme
(`section1_genome_phylogenetics/05_circos_visualization.py`, matplotlib 3.8).

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

### 7. TE landscape ordination

All ordination analyses operated on the composition matrix **X** (236 species ×
234 superfamilies) with entries $x_{ij}$ = total base pairs of superfamily $j$ in
species $i$ (computed by summing the `bp_i` values of all copies belonging to each
superfamily in the `.filteredRepeats.bed` file).  The matrix was constructed after
excluding MI_Encro1 and MI_Ordco1.  Ordination code:
`section2_te_landscape/03_te_landscape_ordination.py`.

#### 7.1 Hellinger-transformed Bray-Curtis PCoA

Relative abundances: $p_{ij} = x_{ij} / \sum_j x_{ij}$.
Hellinger transformation: $h_{ij} = \sqrt{p_{ij}}$.

$$BC_{ab} = \frac{\sum_j |h_{aj} - h_{bj}|}{\sum_j (h_{aj} + h_{bj})}$$

PCoA (classical metric MDS) was applied to the 236 × 236 distance matrix using
`sklearn.manifold.MDS(n_components=10, metric=True, dissimilarity='precomputed',
random_state=42, n_init=4)`.  Variance explained by axis $k$:
$\lambda_k / \sum_{\lambda_k > 0} \lambda_k$ where $\lambda_k$ are the positive
eigenvalues of the double-centred distance matrix.  **PCo1 = 24.9 %, PCo2 = 13.5 %.**

#### 7.2 Aitchison CLR PCoA

Laplace pseudocount of 1 added to all entries.  Centre log-ratio (CLR):

$$c_{ij} = \ln\!\left(\frac{x_{ij}+1}{g_i}\right), \quad g_i = \exp\!\left(\frac{1}{p}\sum_j \ln(x_{ij}+1)\right)$$

where $p = 234$.  Euclidean distances on the CLR matrix equal Aitchison distances.
PCoA applied via the same `sklearn.manifold.MDS` call.

#### 7.3 Jaccard P/A PCoA

Binary matrix $b_{ij} = \mathbf{1}[x_{ij} > 0]$.

$$J_{ab} = 1 - \frac{\sum_j b_{aj} b_{bj}}{\sum_j \max(b_{aj}, b_{bj})}$$

#### 7.4 ANOSIM

For each ordination, analysis of similarities (ANOSIM) tested clade-level group
separation (999 permutations; permutation blocks = species):

$$R = \frac{\bar{r}_B - \bar{r}_W}{N(N-1)/4}$$

where $\bar{r}_B$ = mean rank dissimilarity between clades, $\bar{r}_W$ = within
clades, $N = 236$.  **Bray-Curtis: R = 0.23, p < 0.001; Jaccard: R = 0.18, p < 0.001.**

#### 7.5 Shannon diversity

$$H_i = -\sum_j p_{ij} \log_2 p_{ij}, \quad 0 \log_2 0 := 0$$

H is in bits.  Species with zero total TE content assigned H = 0.

---

### 8. DNMT conservation and 5-methylcytosine analysis

#### 8.1 DNMT1/DIM2 orthogroup identification

A gene product was classified as a DNMT1/DIM2-type DNA methyltransferase if its
InterProScan 5.67-99.0 record contained **both** the following domain accessions
simultaneously within the same protein sequence:
- IPR005105 (DNA cytosine methyltransferase, catalytic domain; Pfam PF00145)
- IPR022686 (PWWP domain; Pfam PF00855)

The PWWP domain is the hallmark of DNMT3-type methyltransferases and DIM2-class
enzymes in fungi; requiring co-occurrence of both domains minimises false positives
from non-maintenance methyltransferases.  Per-species copy numbers were retrieved
from the OrthoFinder orthogroup membership table and projected on the species tree
with ete3 v3.1.3 using the `NodeStyle` colouring module
(`section3_te_epigenetics/01_dnmt_conservation_tree.py`).

#### 8.2 WGBS library processing

Paired-end WGBS reads (Illumina NovaSeq 6000, 2 × 150 bp) were quality-trimmed
with Trim Galore v0.6.10 (Cutadapt v4.6 backend):

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

#### 8.3 Per-copy methylation quantification

The cytosine report was intersected with the TE copy BED file using:

```
bedtools intersect -wa -wb -sorted \
  -a <te_copies.bed> -b <cytosine_report.bed>
```

For each TE copy, the methylation level was computed from covered CpG positions only:

$$m_i = \frac{\sum_{k \in i,\; \text{cov}_k \geq 5} C_k}{\sum_{k \in i,\; \text{cov}_k \geq 5} (C_k + U_k)}$$

where $C_k$ and $U_k$ are the methylated and unmethylated call counts at CpG
position $k$ within copy $i$, and the sum runs over positions with coverage ≥ 5×.
Copies with fewer than 3 such covered CpGs were excluded from all methylation
analyses.  A copy was classified as methylated if $m_i \geq 0.50$.  Spearman's
rank correlation was used to test the monotonic association between $m_i$ and the
Kimura divergence $K_i$ across all retained copies from a given species
(`section3_te_epigenetics/02_5mc_te_age_analysis.py`).

---

### 9. H3K9me3 ChIP-seq and chromatin-category analysis

#### 9.1 ChIP-seq read processing

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

#### 9.2 Peak calling

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

#### 9.3 Per-TE-copy H3K9me3 enrichment and chromatin categories

H3K9me3 enrichment per TE copy was computed as the mean bedgraph signal over the
copy BED coordinates, normalised to RPKM:

$$\text{RPKM}_i = \frac{\text{reads in copy } i \times 10^9}{\text{copy length (bp)} \times \text{total uniquely mapped reads}}$$

TE copies were classified into four mutually exclusive chromatin categories by
crossing two binary marks with thresholds:
- 5mC: methylated if $m_i \geq 0.50$ (from Section 8.3)
- H3K9me3: enriched if RPKM ≥ 1.0

| Category | 5mC | H3K9me3 |
|----------|-----|---------|
| unmethyl_noH3K9 | − | − |
| methyl_noH3K9 | + | − |
| unmethyl_H3K9me3 | − | + |
| methyl_H3K9me3 | + | + |

For each TE copy, the RPKM of the nearest annotated gene (in the same or opposite
strand) was obtained from RNA-seq data.

#### 9.4 RNA-seq for nearest-gene expression

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

---

### 10. Nucleosome positioning — meta-TE profiles

#### 10.1 MNase-seq processing

Micrococcal nuclease (MNase)-digested chromatin was sequenced on an Illumina
platform (2 × 75 bp paired-end).  Sequencing reads were trimmed with Trim Galore
v0.6.10 (`--paired --quality 20 --fastqc`) and mapped to the reference genome with
Bowtie2 v2.5.1 (`--very-sensitive --no-mixed --no-discordant -p 16 -X 200`).  The
`-X 200` parameter limits insert size to 200 bp to select for mononucleosome-protected
fragments (canonical nucleosome footprint ~147 bp plus ~20 bp linker).  Fragment
size distribution was verified with Picard CollectInsertSizeMetrics v3.1.1; samples
were retained only if the modal insert size fell in the range 140–165 bp and
≥ 80 % of fragments were in the 100–220 bp range.  PCR duplicates were removed with
`samtools markdup -r`.

#### 10.2 Dyad probability track computation

Nucleosome positioning and dyad probability scores were computed with DANPOS3 v3.1.1:

```
dpos <sorted.bam> \
  --paired 1 \
  --out <out_dir> \
  --clonalcut 5 \
  --wgsnvf 0
```

`--paired 1` activates paired-end mode, which uses fragment mid-points (dyads) rather
than read 5′ ends.  `--clonalcut 5` removes putative PCR duplicates at single-base
resolution by capping any position's read count at 5 per base.  `--wgsnvf 0` disables
the SNV filter.  The output `.wig` file provides per-base dyad probability scores
(sum normalised to 1 across the genome).

#### 10.3 Meta-TE profile computation

For each TE class, meta-TE nucleosome positioning profiles were computed as follows
(`section3_te_epigenetics/05_nucleosome_te_analysis.py`):

1. For each annotated TE copy of the given class, the genomic window
   [copy_midpoint − 2,000 bp, copy_midpoint + 2,000 bp] was extracted from the
   DANPOS3 Wig track using pyBigWig v0.3.18.
2. The extracted 4,001-position profile was resampled to exactly 150 equally
   spaced bins by linear interpolation (`numpy.interp`).  Bins 38–113 (approximately)
   correspond to the TE body, depending on the relative width of the TE copy vs. the
   ±2 kb window.
3. The 150-bin profiles of all copies were averaged elementwise.  Bins with fewer
   than 10 contributing copies were excluded from the average.
4. The resulting mean profile was plotted; the TE body boundaries (bins
   [floor(150 × 2000 / (4000 + copy_length)) to
    floor(150 × (2000 + copy_length) / (4000 + copy_length))]) are marked with a
   shaded rectangle.  Error bands show ±1 SEM across contributing copies per bin.

---

### 11. N6-methyladenine detection and small RNA silencing

#### 11.1 6mA detection via PacBio kinetics

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

#### 11.2 Small RNA library processing

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

### 12. TE co-option analysis

#### 12.1 Domain-based co-option screen

All OrthoFinder orthogroups were screened for the simultaneous presence of a set of
13 TE-derived InterPro domain families in at least one member protein
(`section4_te_cooption/01_cooption_orthogroup_pipeline.py`).  The query domain
list comprised: IPR001584 (RVT_1, reverse transcriptase core), IPR002104 (DDE
transposase domain), IPR021967 (IS200/IS605 transposase), IPR001878 (CCHC-type
zinc finger), IPR002510 (integrase catalytic domain), IPR000477 (reverse transcriptase
RNA-dependent), IPR004875 (IS3/IS481 transposase), IPR006538 (IS5 transposase),
IPR002951 (IS4 transposase), IPR011010 (DNA breaking-rejoining enzyme),
IPR013762 (integrase-like N-terminal domain), IPR000380 (Mutator-like transposase),
and IPR003086 (SWIM zinc finger).

Orthogroups with ≥ 1 member carrying any of these domains were further validated by
DIAMOND BLASTX v2.1.8 of the complete orthogroup protein sequences against the
cross-species TE consensus library (`--more-sensitive -e 1e-5 --id 30 -c 0.50`).
A protein was classified as **co-opted** — i.e. of cellular rather than TE-encoded
origin — if all three of the following criteria were met:
(i) protein length < 800 aa (full-length TE ORFs in fungi typically exceed 900 aa);
(ii) the GFF3 record of the gene contained ≥ 1 exon ≥ 100 aa in length;
(iii) no full-length hit (query coverage ≥ 80 %) to any consensus in the TE library.
Criterion (iii) distinguishes actively transcribed TE copies (full-length hit) from
truncated or domesticated genes (no full-length hit).

#### 12.2 ChimeraTE chimeric transcript detection

Chimeric gene–TE transcripts were detected with ChimeraTE v1.0 in genome mode
(`section4_te_cooption/02_chimte_analysis.py`).  Each species was processed with:

```
ChimeraTE pipeline \
  --mode genome \
  --genome <reference.fa> \
  --te_annotation <EarlGrey_filtered.bed> \
  --gene_annotation <annotation.gff3> \
  --bam <STAR_sorted.bam> \
  --min_chimeric_reads 3 \
  --te_distance 10000 \
  --threads 8
```

ChimeraTE identifies read pairs in which one mate maps to an annotated TE copy
region and the other mate maps to an annotated gene exon.  The following thresholds
were enforced:
- Minimum 3 read pairs supporting the same TE copy × gene combination.
- The TE-mapping mate must overlap the EarlGrey BED annotation by ≥ 1 bp.
- The gene-mapping mate must overlap an annotated exon by ≥ 1 bp.
- Maximum distance from the TE copy boundary to the gene transcription start site
  (TSS): 10,000 bp.  This excludes insertions in intergenic regions far from any
  gene, focusing on potentially functional TE–gene fusions.

**Shared chimeras** (the same TE family × gene orthogroup junction detected in ≥ 2
species) were identified by inner join on the (te_family_label, orthogroup_id)
column pair across per-species ChimeraTE output tables.

#### 12.3 TF binding site enrichment in TE sequences

Transcription factor binding site (TFBS) motifs from the JASPAR 2024 Core Fungi
database (386 non-redundant motifs; downloaded February 2024) were scanned against
all TE consensus sequences using FIMO v5.5.5 (MEME suite):

```
fimo --oc <out_dir> --thresh 1e-4 --qv-thresh 0.05 \
  --max-stored-scores 500000 \
  <jaspar_fungi_2024.meme> <te_consensi.fa>
```

`--thresh 1e-4` sets the nominal p-value threshold for hit reporting; `--qv-thresh 0.05`
applies an additional Benjamini–Hochberg q-value filter within each sequence.
`--max-stored-scores 500000` increases the default score storage limit to handle
large databases.

For each TF × TE class pair, enrichment was tested with a two-sided Fisher's exact
test on the 2×2 contingency table:

|  | TE class X | Other TE classes |
|--|-----------|-----------------|
| ≥ 1 FIMO hit | a | b |
| No FIMO hit | c | d |

P-values were corrected for multiple testing across all TF × TE class combinations
(n = 386 × number_of_TE_classes) using the Benjamini–Hochberg procedure (FDR 5 %).
The result for HMGb3 × LINE/I-Jockey (p = 0.00275, OR = 4.8, n = 4 species with
both data types) is presented as suggestive owing to the small sample size.

#### 12.4 Multi-domain TE fusion proteins

Proteins harbouring at least one TE-associated domain AND at least one host-function
domain were identified from InterProScan output
(`section4_te_cooption/04_domain_fusion_analysis.py`).  Four canonical fusion classes
were defined by the co-occurrence of the following IPR pair in the same protein:

| Fusion class | TE domain (IPR) | Host domain (IPR) |
|---|---|---|
| FLYWCH–MULE | IPR007527 (FLYWCH Zn-finger) | IPR002992 (Mutator transposase) |
| CotH–DDE | IPR033853 (CotH-like kinase) | IPR002104 (DDE transposase) |
| NACHT–RT | IPR007111 (NACHT ATPase) | IPR000477 (Reverse transcriptase) |
| PARP–Chromo | IPR012317 (PARP catalytic) | IPR000953 (Chromodomain) |

A protein is assigned to a fusion class if both partner IPR accessions appear
within the same protein record in the InterProScan output (same protein identifier,
any database source).

#### 12.5 DAM methyltransferases embedded within TEs

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

### 13. Whole-genome duplication analysis

WGD detection combined three independent lines of evidence computed for each
species: (i) synonymous substitution rate (Ks) distributions of the whole-paranome
(wgd v2), (ii) collinearity-based detection of multiplicons (wgd v2 `syn` module),
and (iii) gene tree reconciliation–based duplication counts from OrthoFinder.
Evidence was integrated into a four-category classification, described in §13.6.
All WGD analyses were restricted to the 131 species with ≤ 50 assembly scaffolds
(to ensure reliable collinearity detection); the full 320-species dataset was used
only for the OrthoFinder reconciliation arm (§13.3).

#### 13.1 Genome pre-filtering for collinearity-based analyses

Before running the wgd v2 synteny and Ks pipelines, assemblies were filtered to
retain only those with ≤ 50 scaffolds/contigs (`scripts/select_genomes_le50_contigs.py`).
Scaffold count was obtained from the FASTA index (`.fai`) generated with `samtools faidx`.
This threshold was chosen because wgd v2 `syn` is computationally intractable and
statistically unreliable on highly fragmented assemblies where true collinear blocks
are interrupted by assembly gaps.  The 131 assemblies passing this threshold cover
all 18 clades and include all species previously described as WGD candidates in the
literature.

#### 13.2 Whole-paranome Ks distribution — wgd v2 `dmd` + Ks pipeline

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

#### 13.3 OrthoFinder gene tree reconciliation — duplication mapping

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

#### 13.4 GMM fitting and WGD peak detection in the paranome Ks

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

#### 13.5 Collinearity analysis — wgd v2 `syn` and multiplicons

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

#### 13.6 Multi-evidence WGD classification

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

#### 13.7 TE breakpoint enrichment at WGD synteny boundaries

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

#### 13.8 TE burst — WGD age correlation

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

#### 13.9 Normalised gene duplication rates per phylogenetic node

The normalised rate $R_n$ (§13.3) was plotted on the species tree as a heat-map
column, using the same colour scale across all nodes.  For visualisation, ancestral
node rates were projected onto the LG+C60+F+G4 ML species tree with ete3 v3.1.3
(`NodeStyle` colouring, `TreeStyle` rectangular layout).  Duplication-peak nodes
(Robust Z ≥ 1.5, §13.3) are highlighted with a filled circle; the size is
proportional to $Z_n$.

---

### 14. Horizontal transposon transfer (HTT)

#### 14.1 HTT candidate identification

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

#### 14.2 HTT dating

For each confirmed HTT event, the transfer date was estimated as the median
synonymous Ks among all orthologue pairs between the donor and recipient species
at the relevant OrthoFinder orthogroup, computed with yn00 (PAML v4.10.10) as in
Section 13.1.  Ks was converted to Mya using the per-lineage $r_{\text{syn}}$
from Section 13.4.

#### 14.3 Flanking sequence locus conservation

For each HTT candidate TE locus, 10,000 bp of flanking sequence on each side was
extracted.  Both flanks were queried against the proteome and genome of the putative
donor species via BLASTN (≥ 60 % identity, ≥ 50 bp alignment length, E ≤ 10⁻⁵,
`-task blastn-short`):

$$f_{\text{cons}} = \frac{\sum_{\text{hits}} \text{alignment\_length}_{\text{hit}}}{\text{total flank length (bp)}}$$

Loci were classified as:
- **Conserved**: $f_{\text{cons}} \geq 0.30$
- **Partial**: $0.10 \leq f_{\text{cons}} < 0.30$
- **Non-conserved**: $f_{\text{cons}} < 0.10$

#### 14.4 HTT event counting — Dollo parsimony

HTT events were counted under a Dollo parsimony model (an element can be gained only
once, losses are free) using the ete3 v3.1.3 `TreeNode.get_descendant_evol_events()`
framework.  The full candidate set comprised 224 BLASTN-flagged TE family × clade-pair
combinations.  After applying the anomaly-score filter ($a \geq 0.10$) and the
locus-conservation classification (≥ partial), the parsimony-minimum number of
independent transfer events across the species tree was **11**.

---

### 15. Statistical framework

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

### Figure 1. TE landscapes and genomic context across 320 fungal and protistan genomes

**(A)** Cladogram of the 320 species organised into 18 clades (Section 1, Section 2).
Tip colours indicate the kingdom-level classification (Fungi, Oomycota, Chytridiomycota,
Microsporidia, Neocallimastigomycota).  Node shading indicates the presence (orange)
or absence (grey) of DNMT1/DIM2 in the MRCA of each clade.  Scale bar: 100 My.

**(B)** Genome size (Mb, log₁₀ scale) shown as a bar plot on the right of the
cladogram.  Error bars: not applicable (one assembly per species).  Source data:
Fig1B_genome_size.

**(C)** TE content as a percentage of the genome assembly size.  Bars are stacked by
Wicker et al. (2007) class: Class I LTR (purple), Class I non-LTR (blue), Class II
DNA transposons (orange), Helitrons (yellow), Unknown (grey).  Source data:
Fig1C_te_content_class.

**(D)** Macrosynteny chord diagram.  Coloured ribbons connect orthologous syntenic
blocks between species (Methods Section 4).  Ribbon width is proportional to the
number of SynChro anchor gene pairs in the block.  Only blocks with ≥ 5 anchor
pairs are shown.  Species are arranged in phylogenetic order on the outer ring.
Source data: Fig1D_synteny.

**(E)** Spearman rank correlation between genome size and TE content across all 320
species ($\rho$, colour-coded by clade; two-tailed; BH-corrected).  Dotted line:
$\rho = 0$.

---

### Figure 2. TE landscape diversity and ordination

**(A)** PCoA biplot of Hellinger-transformed Bray-Curtis dissimilarity across 236
species (Sections 5–7).  Each dot is one species; colour indicates clade.  Ellipses:
95 % confidence intervals per clade (multivariate normal).  Axis labels include the
percentage of variance explained.  Source data: Fig2A_te_content.

**(B)** ANOSIM results for each ordination metric (Bray-Curtis R = 0.23, Jaccard
R = 0.18; 999 permutations; p < 0.001 for both).  Shown as a dot-range plot with
observed R (dot) and permuted distribution (violin).

**(C)** Shannon diversity index $H_i$ plotted on the phylogeny.  Tip colours
indicate clade; heat-map column shows H value per species (log₂ scale).
Source data: Fig2C_shannon.

**(D)** Kimura divergence landscape selected examples.  Stacked histograms (bin
width = 0.005) for four representative species spanning high-TE (genome > 50 % TE)
and low-TE (< 5 % TE) content.  Colours: as in (C).  Source data: Fig2D_kimura.

---

### Figure 3. Epigenetic regulation of TEs in fungi

**(A)** Cladogram coloured by DNMT1/DIM2 presence (orange) or absence (grey) in
extant taxa (Section 8.1).  Ancestral-state reconstruction using parsimony (ete3
v3.1.3) indicates the predicted number of independent DNMT losses (filled squares
on nodes = loss events).  Source data: Fig3A_dnmt.

**(B)** Per-copy methylation level ($m_i$, x-axis) versus Kimura K2P distance
($K_i$, y-axis) scatter plot for all TE copies with ≥ 3 covered CpGs in
DNMT-positive species (Section 8.3).  Spearman $\rho$ shown per species.
Source data: Fig3B_methylation_age.

**(C)** H3K9me3 RPKM distribution across the four chromatin categories (Section 9.3)
shown as a violin+jitter plot.  Horizontal line: RPKM = 1.0 threshold.
Source data: Fig3C_H3K9.

**(D)** Nearest-gene RNA-seq RPKM (log₁₀, y-axis) per chromatin category (x-axis).
Pairwise Mann–Whitney U-tests (BH-corrected) indicated above brackets.
Source data: Fig3D_gene_expression.

**(E)** Meta-TE nucleosome positioning profile (DANPOS3 dyad score, Section 10)
averaged across all Class I LTR copies in DNMT-positive species.  Shaded region:
TE body.  Error band: ±1 SEM across copies.  Source data: Fig3E_nucleosome.

**(F)** Small RNA mapping (RPM, y-axis) per TE superfamily versus mean H3K9me3 RPKM
(x-axis) in piRNA-pathway competent species (Section 11.2).  Spearman $\rho$
and p-value shown.  Source data: Fig3F_srna.

---

### Figure 4. TE co-option diversity

**(A)** Heatmap of co-opted domain families across 18 clades.  Rows: 13 TE-derived
IPR families (Section 12.1).  Columns: clades.  Cell colour: fraction of species
within the clade harbouring ≥ 1 co-opted gene with that domain.  Hierarchical
clustering (Ward linkage, Jaccard distance) applied to rows and columns.

**(B)** Heatmap of statistically supported domain expansions in TE-rich species
(Sections 12.1 and 7.1 ordination PCo1 > 0).  Rows: Pfam domains with fold-change
vs. TE-poor species ≥ 2.0 AND Mann–Whitney BH < 0.05 in ≥ 1 clade.  Columns:
clades.  Cell colour: fraction of clade species with significant expansion.
Source data: Fig4B_cooption_heatmap.

**(C)** Shared chimeric gene–TE junctions detected by ChimeraTE (Section 12.2).
Each row is a TE family × gene orthogroup pair; columns are clades.  Shared
(present in ≥ 2 species in ≥ 1 clade) junctions only.  Source data: Fig4C_chimte.

**(D)** Multi-domain TE fusion network.  Nodes: protein families (diamonds) or TE
superfamilies (circles).  Edges: ≥ 1 protein detected harbouring both domains.
Edge width: number of species.  Node colour: domain type (TE vs. host).  Cytoscape
yFiles Organic layout.  Source data: Fig4D_fusions.

**(E)** Frequency of the four canonical fusion classes (Section 12.4) across the
18 clades.  Stacked bar: proportion of species with each fusion class.
Source data: Fig4E_fusions.

**(F)** DAM methyltransferase co-occurrence with AMT1 (Section 12.5).  Grouped
stacked bar per clade (n species): DAM-TE positive + AMT1 (dark), DAM-TE positive
only (medium), AMT1 only (light), neither (white).  Fisher's exact test inset:
OR = 5.19, p = 0.028.  Source data: Fig4F_dam_amt1.

---

### Figure 5. Whole-genome duplications and TE dynamics

**(A)** Paralogous Ks distributions for 12 WGD-candidate species.  Each panel shows
the histogram of Ks values (grey bars) and the fitted GMM (coloured Gaussian
components, BIC-selected k; Section 13.2).  WGD peak components annotated with
$\mu_k$ and weight.  Source data: Fig5A_paralogs.

**(B)** TE density at synteny breakpoints vs. random windows.  Paired dot-line plot
per WGD species (x: random, y: breakpoint); colour: Wilcoxon signed-rank result
(BH < 0.05 = orange, non-significant = grey).  Source data: Fig5B_breakpoints.

**(C)** Per-superfamily TE enrichment at breakpoints.  Heatmap: −log₁₀ BH-corrected
p-value; rows: TE superfamilies; columns: WGD species.  Only superfamilies
significant in ≥ 1 species shown.  Dot size: fold-enrichment.
Source data: Fig5C_te_class_bp.

**(D)** WGD age vs. TE burst KDE-mode scatter plot (Section 13.4).  Each point:
one WGD species × TE superfamily pair with significant Spearman $\rho$.  Colour:
superfamily class.  Dashed line: y = x.  Source data: Fig5D_dup_rate.

---

### Figure 6. Horizontal transposon transfer (HTT)

**(A)** HTT network.  Nodes: species (colour = clade).  Edges: confirmed HTT event
(anomaly score $a \geq 0.10$, Section 14).  Edge width: $a$ value.  Only the
minimum spanning tree per connected component is shown for clarity.  Layout:
phylogenetic circular (ete3).

**(B)** Anomaly score ($a$) distribution across all 224 candidate pairs (grey)
versus 11 confirmed events (orange).  Dashed line: threshold $a = 0.10$.

**(C)** HTT age distribution (Mya) for 11 confirmed events, estimated from median
orthologue Ks (Section 14.2).

**(D)** Locus conservation (Section 14.3) for confirmed HTT events.  Stacked bar:
conserved / partial / non-conserved.

**(E)** Dollo parsimony reconstruction of HTT events on the species tree (Section
14.4).  Nodes with inferred transfer events are highlighted.

**(F)** TE superfamily composition of the 11 confirmed HTT events.  Pie chart coloured
by Wicker class.


---

## Part III — Source Data Index

The Excel workbook `EDF-TE_source_data.xlsx` contains 28 sheets.  Each sheet
provides the numerical data underlying one or more figure panels as described below.

| Sheet name | Figure panel(s) | Description | Key columns |
|---|---|---|---|
| Fig1B_genome_size | Fig. 1B | Genome assembly statistics per species | species, clade, assembly_size_mb, scaffold_n50_kb, busco_complete_pct |
| Fig1C_te_content_class | Fig. 1C | TE content (bp and %) per class per species | species, clade, te_class, te_bp, te_pct |
| Fig1D_synteny | Fig. 1D | SynChro macrosynteny block table | species_A, species_B, block_id, n_anchors, block_length_bp, chr_A, start_A, end_A, chr_B, start_B, end_B |
| Fig1E_te_genome_corr | Fig. 1E | Spearman correlation: genome size vs. TE content | clade, rho, pvalue, BH_qvalue |
| Fig2A_te_content | Fig. 2A | TE superfamily bp matrix for PCoA | species, superfamily (234 columns), PCo1, PCo2, clade |
| Fig2B_anosim | Fig. 2B | ANOSIM results per ordination metric | metric, R_observed, pvalue, n_permutations |
| Fig2C_shannon | Fig. 2C | Shannon H per species | species, clade, shannon_H, total_te_bp |
| Fig2D_kimura | Fig. 2D | Kimura K2P per copy (subset for 4 representative species) | species, te_family, te_class, copy_id, K2P, copy_bp, recency_index |
| Fig3A_dnmt | Fig. 3A | DNMT1/DIM2 presence per species + ancestral states | species, clade, dnmt_copy_number, dnmt_present, parsimony_state |
| Fig3B_methylation_age | Fig. 3B | Per-copy methylation vs. K2P (DNMT+ species) | species, copy_id, te_class, m_i, K_i, n_CpG_covered |
| Fig3C_H3K9 | Fig. 3C | H3K9me3 RPKM per copy and category | species, copy_id, te_class, H3K9me3_RPKM, category_5mC, category_H3K9, chromatin_category |
| Fig3D_gene_expression | Fig. 3D | Nearest-gene RPKM per chromatin category | species, gene_id, chromatin_category, gene_RPKM |
| Fig3E_nucleosome | Fig. 3E | Meta-TE nucleosome profile (150 bins) | bin_index, mean_dyad_score, sem_dyad_score, n_copies, te_class |
| Fig3F_srna | Fig. 3F | sRNA RPM per TE family vs. H3K9me3 RPKM | species, te_family, sRNA_RPM, mean_H3K9_RPKM, n_copies |
| Fig4A_cooption_heatmap_raw | Fig. 4A | Co-opted domain presence per species | species, clade, ipr_id, ipr_name, n_coopted_genes, present |
| Fig4B_cooption_heatmap | Fig. 4B | Domain expansion heatmap data | pfam, domain_desc, clade, expansion_class, n_positive, pct_positive, fold_vs_rest, BH_qvalue |
| Fig4C_chimte | Fig. 4C | Shared ChimeraTE junctions | te_family, orthogroup_id, species_list, n_species, clade |
| Fig4D_fusions_network | Fig. 4D | Fusion protein network edges | protein_id, species, te_domain_ipr, host_domain_ipr, n_species |
| Fig4E_fusions | Fig. 4E | Four canonical fusion class frequencies per clade | fusion_class, clade, n_species_with_fusion, total_species, fraction |
| Fig4F_dam_amt1 | Fig. 4F | DAM-TE and AMT1 co-occurrence | species, clade, dam_te_present, amt1_status, fisher_OR, fisher_p |
| Fig5A_paralogs | Fig. 5A | GMM-fitted Ks distribution per WGD species | species, Ks, gmm_component, mu_k, weight_k, BIC_k |
| Fig5B_breakpoints | Fig. 5B | Breakpoint vs. random TE density (Wilcoxon) | species, breakpoint_id, te_density_breakpoint, te_density_random, wilcoxon_stat, wilcoxon_p, BH_q |
| Fig5C_te_class_bp | Fig. 5C | Per-superfamily breakpoint enrichment | species, superfamily, te_density_breakpoint, te_density_random, fold_enrichment, wilcoxon_p, BH_q |
| Fig5D_dup_rate | Fig. 5D | WGD age vs. TE burst KDE mode | species, wgd_age_mya, superfamily, te_burst_Ks_mode, spearman_rho, spearman_p, BH_q |
| Fig6A_htt_network | Fig. 6A | HTT event edge list | species_donor, species_recipient, te_family, anomaly_score, confidence |
| Fig6B_anomaly_scores | Fig. 6B | Full candidate anomaly score distribution | species_A, species_B, te_family, a_score, confirmed |
| Fig6C_htt_ages | Fig. 6C | HTT event age estimates | event_id, species_donor, species_recipient, te_family, Ks_median, age_mya |
| Fig6D_locus_conservation | Fig. 6D | Locus conservation classification | event_id, f_cons, conservation_class |

All values in this workbook are the direct numerical inputs to the corresponding
figure panels.  Processed intermediate tables are available in the per-section
`results/` directories of the code repository; their paths and formats are defined
in `config.sh`.

