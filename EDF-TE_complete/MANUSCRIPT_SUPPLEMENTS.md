# EDF-TE — Materials and Methods, Figure Legends, and Source Data

---

## PART I — MATERIALS AND METHODS

### 1. Genome dataset and taxonomic sampling

We assembled a dataset of 320 eukaryotic genomes spanning ~1 billion years of evolutionary divergence, covering 18 major lineages: Holozoa (HO), Nucleariida (NU), Aphelida (AP), Cryptomycetes (CR), Microsporidia (MI), Blastocladiomycetes (BL), Olpidiales (OL), Chytridiomycetes (CH), Neocallimastigomycota (NE), Bacillariophyta (BS), Entomophthorales (EN), Zoopagomycota (ZO), Kickxellomycota (KX), Glomeromycota (GL), Mortierellomycota (MO), Mucoromycota (MU), Basidiomycota (BA), and Ascomycota (AS). Neocallimastigomycota species (CH_Anasp1, CH_Neocon1, CH_Neolan1, CH_Orpsp1, CH_Pecora1, CH_Caecom1, CH_Pirfi3) were reclassified from the CH prefix and handled separately in all clade-level analyses. Genome assemblies were downloaded from the Joint Genome Institute (JGI) MycoCosm portal, NCBI RefSeq, and the European Nucleotide Archive. All genome FASTA files were pre-processed with `seqtk seq -A` to strip non-ACGTN characters, and chromosomal scaffolds shorter than 200 bp were removed. Protein annotations (GFF/GFF3) and functional domain annotations (InterProScan 5.67 output in IPR.tab format) were obtained from JGI or re-computed with InterProScan 5.67-99.0 using Pfam 36.0, PRINTS, PANTHER 18.0, TIGRFAM, and CDD databases.

Three additional Mucoromycota species sequenced as part of this study (*Mucor circinelloides* strain MU_Mucci3, *Rhizopus microsporus* strain MU_Rhimi59, and *Phycomyces blakesleeanus* strain MU_Phybl4) were assembled with Flye 2.9 from PacBio HiFi reads. CCS reads were generated after separating rDNA and mtDNA reads. Assembly completeness was assessed with BUSCO 5.8.2 (fungi_odb10 lineage dataset).

---

### 2. Species tree inference

A maximum-likelihood species tree was inferred from 238 fungal and protistan genomes (Section 1). Single-copy orthogroups were identified with OrthoFinder 3.1 (parameters: `-M msa -A mafft -T fasttree -S diamond`). From the resulting orthogroups, markers were selected that were present in ≥40% of species (≥95 of 238 genomes) and with ≤2 copies per species. Each marker was aligned with MAFFT v7.520 (`--auto`) and trimmed with trimAl v1.4.1 (`-automated1`). Markers with fewer than 60 post-trimming alignment columns, fewer than 25 parsimony-informative sites, or relative composition variability (RCV) >0.5 (PhyKIT v2.0) were excluded. Per-gene trees were inferred under the LG+G4 model with IQ-TREE 2.3.6; genes containing taxa with branch lengths >3.5 standard deviations above the per-gene mean were removed as likely paralogs. The remaining markers were concatenated into a supermatrix, and the partition file was submitted to IQ-TREE 2.3.6 for final inference under a partitioned LG+C60+F+G4 model with 1,000 ultrafast bootstrap replicates (`-B 1000`). The tree was rooted on the Holozoa outgroup. All steps were executed through a 10-step bash pipeline (`section1_genome_phylogenetics/01_species_tree_pipeline.sh`) run under the `phylo` micromamba environment.

---

### 3. Genomic feature characterisation

Six genome-level metrics were computed for each species (`section1_genome_phylogenetics/02_genomic_features.py`): (i) genome size (Mb), (ii) GC content (%), (iii) gene count, (iv) mean gene length (bp), (v) intron density (introns per gene), and (vi) N50 scaffold length (Mb). All values were extracted from the genome FASTA and the corresponding GFF/GFF3 annotation. Pan-genome accumulation curves (core, shell, unique gene counts as a function of the number of sampled genomes) were computed by 1,000-times random permutation of species addition order (`section1_genome_phylogenetics/03_pangenome_analysis.py`).

---

### 4. Synteny analysis

Pairwise macrosynteny was assessed within and across major clades using SynChro (Section 1, `04_synteny_analysis.sh`). For each Mucoromycota species pair, protein sequences were compared with DIAMOND BLASTP (identity ≥30%, e-value ≤10⁻⁵, coverage ≥50%), and reciprocal best hits were used to anchor syntenic blocks (minimum 5 anchors, maximum gap 25 genes). Syntenic block coordinates were converted to BED format and visualised as circular chord diagrams with a custom Python script using matplotlib (`05_circos_visualization.py`). The chord width is proportional to the number of syntenic anchor pairs.

---

### 5. TE annotation

Transposable elements were annotated in all 320 genomes with EarlGrey v1.6 (`section2_te_landscape/01_te_annotation.sh`). EarlGrey integrates RepeatModeler2 (de novo consensus discovery) with phylogenetically-aware classification, followed by RepeatMasker 4.1.5 re-annotation using the curated per-species library. All jobs were submitted to a SLURM cluster (partition `normal.168h`, 16 CPUs, 64 GB RAM per species). After per-species annotation, all consensus sequences were merged, clustered at 80% identity using MMseqs2 (`easy-cluster`, coverage 0.80, `--cov-mode 1`), and sequences with DIAMOND hits to RefSeq non-redundant proteins (identity ≥40%, e-value ≤10⁻⁵) were removed as host-derived. The resulting cross-species TE library was used to re-classify each species-level consensus. TE superfamily nomenclature follows the Wicker et al. (2007) hierarchical classification scheme. Two Microsporidia species lacking detectable TEs (MI_Encro1, MI_Ordco1) were excluded from ordination analyses.

---

### 6. TE landscape ordination

Species-level TE landscapes were represented as a matrix of 236 species × 234 superfamilies, with values equal to base pairs of TE sequence per unit genome (relative abundance). Three complementary ordinations were computed (`section2_te_landscape/03_te_landscape_ordination.py`): (i) principal coordinates analysis (PCoA) on Bray-Curtis dissimilarities computed from Hellinger-transformed relative abundances (PCo1 = 24.9%, PCo2 = 13.5% of variance); (ii) PCoA on Aitchison distances from centre log-ratio (CLR) transformed counts (Laplace pseudocount = 1 for zero-imputation); and (iii) PCoA on Jaccard distances from binarised presence/absence superfamily profiles. For each ordination, analysis of similarities (ANOSIM, 999 permutations) tested whether broad clade membership explained TE composition differences better than expected by chance (Bray-Curtis ANOSIM R = 0.23; Jaccard ANOSIM R = 0.18). Shannon diversity (H) was computed over superfamily-level proportional abundances.

---

### 7. Kimura divergence and TE age

TE insertion age distributions were derived from the substitution frequency matrices in RepeatMasker `.align` files, following the Kimura two-parameter (K2P) model implemented in EarlGrey's divergence pipeline (`section2_te_landscape/04_kimura_divergence.py`). For each species, the recency index (RI) was defined as the fraction of TE base pairs with Kimura divergence <0.10 (corresponding to young, recently active copies). Kimura landscapes were visualised as stacked histograms binned at 0.005 intervals per TE class.

---

### 8. DNMT conservation and 5-methylcytosine analysis

Orthologs of fungal DNA methyltransferases were identified via InterPro domain signatures: DNMT1/DIM2 were defined by simultaneous presence of IPR005105 (catalytic domain) and IPR022686 (PWWP domain). Presence/absence was scored per species and projected onto the species tree using ete3 (`section3_te_epigenetics/01_dnmt_conservation_tree.py`).

Whole-genome bisulphite sequencing (WGBS) data were processed with Bismark 0.24.2 (`--genome` option, paired-end). CpG methylation calls were extracted with `bismark_methylation_extractor`. TE copy-level methylation was computed by intersecting CpG cytosine positions (minimum 5× coverage) with TE BED coordinates (minimum 3 CpGs per copy, minimum 50% methylated cytosines at covered positions for a copy to be considered methylated) (`section3_te_epigenetics/02_5mc_te_age_analysis.py`). Methylation level was correlated with Kimura divergence using Spearman's rank correlation.

---

### 9. Chromatin and nucleosome analysis

H3K9me3 ChIP-seq enrichment over TEs was computed by intersecting per-base signal (bedgraph format from MACS2 callpeak with `--nomodel --bdg`) with TE BED intervals. TE copies were classified into four chromatin categories based on joint methylation and H3K9me3 status: unmethylated/no-H3K9me3, methylated/no-H3K9me3, unmethylated/H3K9me3, and methylated/H3K9me3 (`section3_te_epigenetics/03_chromatin_proximity_analysis.py`). RNA-seq expression of the nearest gene to each TE copy was quantified with featureCounts (Subread 2.0.6) and compared across chromatin categories.

Nucleosome positioning maps were obtained from DANPOS3 MNase-seq data. Meta-TE profiles (±2 kb flanks, 150-bin resolution) were generated by averaging per-nucleotide signal across all TE copies of each class (`section3_te_epigenetics/05_nucleosome_te_analysis.py`).

---

### 10. 6mA and sRNA silencing

PacBio CLR kinetic data were used to detect N6-methyladenine (6mA) marks. Inter-pulse duration (IPD) ratios relative to an in-silico reference were extracted using `ipdSummary`. Positions with IPD ratio ≥5.0 were classified as 6mA-positive. TE copies containing ≥1 6mA-positive site were compared for RNA-seq coverage against unmarked copies (`section3_te_epigenetics/04_6mA_expression_analysis.py`).

Small RNA (sRNA) reads (20–25 nt) were trimmed with Cutadapt 4.6 and mapped to the cross-species TE consensus library with Bowtie 1.3.1 (`-v 2 -M 1 --best --strata`). Reads per million (RPM) was calculated per TE family. Spearman correlation was tested between sRNA RPM and H3K9me3 enrichment per TE copy (`section3_te_epigenetics/06_srna_silencing_pipeline.py`).

---

### 11. TE co-option analysis

#### 11a. Domain-based co-option screen

All OrthoFinder orthogroups were screened for simultaneous presence of TE-derived InterPro domains (`section4_te_cooption/01_cooption_orthogroup_pipeline.py`). Thirteen TE-associated domain families were used as queries (IPR001584, IPR002104, IPR021967, IPR001878, IPR002510, IPR000477, IPR004875, IPR006538, IPR002951, IPR011010, IPR013762, IPR000380, IPR003086). Orthogroups with TE domains were validated with DIAMOND BLASTX (identity ≥30%, e-value ≤10⁻⁵) against the cross-species TE library to confirm homology. Co-opted genes were defined as those with TE domain evidence but lacking TE-length ORF structure (protein <800 aa and ≥1 exon with ≥100 aa).

#### 11b. ChimeraTE analysis

Chimeric TE–gene transcripts were detected with ChimeraTE (section 4, `02_chimte_analysis.py`) applied to RNA-seq data. ChimeraTE output directories were parsed to aggregate chimeric event counts per TE class and per species. Shared chimera events (identical donor TE family + recipient gene orthogroup) across species were identified using inner merge on the chimera call tables.

#### 11c. TF–TE interaction analysis

FIMO (MEME suite 5.5.5) was used to scan TE consensus sequences for transcription factor binding site motifs from the Jaspar 2024 Core vertebrate/fungi database (`section4_te_cooption/03_tf_te_interaction_analysis.py`). Hits with p ≤ 10⁻⁴ and q ≤ 0.05 were retained. Fisher's exact test was applied to test enrichment of each TF motif in each TE class (HMGb3 × LINE/I-Jockey: OR = ?, p = 0.00275, n = 4 species). A TF × TE class enrichment heatmap was plotted with hierarchical clustering (Ward linkage, Euclidean distance on –log10 p-values).

#### 11d. Multi-domain fusion proteins

Proteins were screened for four canonical TE–host domain fusions: (i) FLYWCH–MULE (IPR007527 + IPR002992), (ii) CotH–DDE (IPR033853 + IPR002104), (iii) NACHT–RT (IPR007111 + IPR000477), and (iv) PARP–Chromo (IPR012317 + IPR000953) (`section4_te_cooption/04_domain_fusion_analysis.py`). A protein was classified as a fusion if both partner domains appeared in the InterProScan output for the same protein identifier.

#### 11e. DAM-like methyltransferases and TE insertion

Proteins with HMMER hits (e-value ≤10⁻⁵) to the DAM methyltransferase HMM profile (Pfam PF05869) were identified in all 320 proteomes (`section4_te_cooption/05_dam_te_methyltransferase.py`). Each DAM-containing gene was cross-referenced with TE BED coordinates (bedtools closest, maximum distance 2 kb) to determine whether the gene was TE-associated. Fisher's exact test assessed co-occurrence of DAM-TE genes with AMT1-like adenine methyltransferases (IPR014885): OR = 5.19, p = 0.028. Fifteen species harboured DAM-TE genes (14 associated with DIRS-class elements; 1 with Gypsy LTR).

---

### 12. Whole-genome duplication detection

Paralogous gene pairs were identified within each Mucoromycota species from OrthoFinder 3.1 self-comparison Ks output. Synonymous substitution rates (Ks) were computed with the yn00 module of PAML 4.10.10 after codon-level alignment with PAL2NAL v14 applied to MAFFT protein alignments. Gaussian mixture models (GMM; up to five components, sklearn 1.4 GaussianMixture, BIC-selected optimal k) were fitted to per-species Ks distributions filtered to 0 < Ks ≤ 3.0 (`section5_wgd_te_dynamics/01_wgd_detection.py`). WGD events were called when a GMM component contained ≥10% of the data and had a mean Ks in the range 0.05–3.0.

TE density at synteny breakpoints (±5 kb from block ends) was compared with 100 size-matched random windows per breakpoint from the same chromosome using Wilcoxon signed-rank test (`section5_wgd_te_dynamics/02_te_breakpoint_enrichment.py`). Ks-estimated WGD times were converted to absolute age using per-lineage synonymous substitution rates calibrated against the species tree. TE burst timing (mode of Kimura divergence distribution) was correlated with WGD time using Spearman's rank correlation (`section5_wgd_te_dynamics/03_te_age_wgd_correlation.py`). Per-superfamily copy number was compared between WGD and non-WGD species using Mann-Whitney U tests with FDR correction (Benjamini-Hochberg) (`section5_wgd_te_dynamics/05_te_copy_number_evolution.py`).

---

### 13. Horizontal TE transfer detection and dating

Per-clade TE libraries were constructed by concatenating EarlGrey consensus sequences for all species within each clade, followed by within-clade clustering at 95% nucleotide identity (MMseqs2 `easy-cluster`, coverage 0.80) to remove redundancy from closely related species (`section6_htt/01_htt_phylo_prep.sh`). Clustered libraries were re-classified with RepeatClassifier (RepeatModeler 2.0.5).

Cross-clade TE similarity was assessed by all-against-all BLASTN between every pair of clade libraries (identity ≥85%, alignment length ≥100 bp, word size 11). For each hit pair, an anomaly score was computed as the difference between the expected nucleotide divergence (derived from calibrated tip-to-tip distances in the species tree) and the observed identity (1 − BLASTN pident/100). Candidate horizontal transfer events were defined as pairs with anomaly score ≥0.10 (`section6_htt/02_htt_candidate_detection.py`). This threshold corresponds to TE copies at least 10% more similar across clades than expected from the species divergence time.

Transfer events were dated by Ks-based calibration: for each candidate donor–recipient clade pair, the median Ks of single-copy OrthoFinder orthologs (Section 12 pipeline) provided the divergence time lower bound (`section6_htt/03_htt_ks_dating.py`). Insertion locus conservation was assessed by BLASTN of 10 kb flanking sequences of each TE copy in the donor against the recipient genome (60% identity, 50 bp minimum, fraction of flanking length with hits reported as conservation score; `section6_htt/04_htt_locus_conservation.py`). A minimum of 11 HTT events (parsimony-reconciled on the species tree) explains the 224 transfer candidates detected across the dataset.

---

### 14. Statistical analyses and data visualisation

All statistical tests are two-tailed unless otherwise specified. Multiple-testing correction used the Benjamini-Hochberg procedure (false discovery rate q < 0.05) wherever >5 simultaneous hypotheses were tested. All visualisations were produced in Python 3.11 with matplotlib 3.8 and seaborn 0.13, or in R 4.4 with ggplot2 4.0. Network visualisations used networkx 3.2 (Python) and layouts were exported to Cytoscape 3.10.2 for publication-quality rendering. All scripts, configuration files, and the full reproducible pipeline are available at https://github.com/cLax8M/EDF-TE in the `EDF-TE_complete/` directory.

---

## PART II — FIGURE LEGENDS

### Figure 1 — Phylogenomic framework and genomic feature landscape

**(A)** Maximum-likelihood species tree of 238 eukaryotic genomes inferred from a concatenated supermatrix of single-copy orthogroups (LG+C60+F+G4, 1,000 ultrafast bootstrap replicates). Branches are coloured by clade (see key). Node circles indicate bootstrap support (filled: ≥95; open: <95). **(B)** Heatmap of six genomic features (genome size, GC content, gene count, mean gene length, intron density, N50) per species, ordered as in (A). Values are z-score normalised per column. **(C)** Pan-genome accumulation curves for Mucoromycota species: core (present in all), shell (2–99%), and unique (singletons) orthogroups as a function of the number of sampled genomes (mean ± SD of 1,000 random sampling orders). **(D)** Circular chord diagram of macrosyntenic blocks among Mucoromycota species. Chord width is proportional to the number of syntenic anchor pairs. Chords are coloured by the taxonomic clade of the anchor-donating species.

---

### Figure 2 — TE landscape diversity across the fungal kingdom

**(A)** Cladogram of all 320 species (topology from Figure 1A, pruned) with stacked horizontal bars showing the proportion of genome occupied by each TE class (LTR, DNA, LINE, SINE, Rolling Circle, Penelope, Unclassified). Bar width is proportional to total TE content (% genome). Shannon diversity (H) is shown as a continuous colour track adjacent to the tip labels. **(B)** Principal coordinates analysis (PCoA) of the 236-species × 234-superfamily Bray-Curtis dissimilarity matrix (Hellinger-transformed abundances). Points are coloured by clade. PCo1 = 24.9%, PCo2 = 13.5% of total variance. Convex hulls enclose each clade. ANOSIM R = 0.23, p < 0.001 (999 permutations). **(C)** Kimura divergence landscapes (stacked histograms, bin width 0.005) for representative clades. The recency index (RI; fraction of TE bp with Kimura < 0.10) is shown for each. **(D)** Shannon diversity boxplot per clade, ordered phylogenetically. Individual species shown as jittered points.

---

### Figure 3 — TE epigenetic regulation

**(A)** Species tree of 320 taxa with binary presence/absence of DNMT1/DIM2 (defined by co-occurrence of IPR005105 + IPR022686) annotated as tip symbols. Node bars show orthogroup copy number (OrthoFinder 3.1). **(B)** Violin plot of TE copy methylation level (% CpGs methylated, WGBS) stratified by Kimura divergence quartile. Each point is a TE copy. Spearman ρ and p-value are shown. Copies with <3 CpGs at ≥5× coverage are excluded. **(C)** Box plots of RNA-seq expression (log2 RPKM + 1) of the nearest gene to TE copies in each of four chromatin categories: unmethylated/H3K9me3-negative (open), methylated/H3K9me3-negative (blue), unmethylated/H3K9me3-positive (orange), methylated/H3K9me3-positive (red). **(D)** Meta-TE nucleosome profiles (DANPOS3 MNase-seq, ±2 kb, 150 bins) per TE class in *Rhizopus microsporus* MU_Rhimi59. The TE body is shaded grey. **(E)** Proportion of sRNA (20–25 nt) reads mapping to each TE superfamily (top 15), normalised to RPM. **(F)** Scatter: H3K9me3 enrichment per TE copy (x-axis) vs sRNA RPM from the same family (y-axis). Spearman correlation and regression line shown.

---

### Figure 4 — TE co-option across fungal genomes

**(A)** Bubble plot of co-opted orthogroups per species: bubble area = number of species with the orthogroup, colour = dominant TE domain class. **(B)** Heatmap of co-option evidence per clade (rows: TE domain type; columns: 18 clades). Cell colour = fraction of species with ≥1 co-opted gene in each domain class. **(C)** Proportion of chimeric TE–gene transcripts per TE class (ChimeraTE analysis). **(D)** TF binding site enrichment in TE sequences: hierarchically clustered heatmap of –log10(Fisher p) for 20 TF motifs × 10 TE classes. Red cells indicate significant enrichment (FDR < 0.05). **(E)** Multi-domain fusion protein occurrence per species for the four fusion classes (FLYWCH–MULE, CotH–DDE, NACHT–RT, PARP–Chromo). **(F)** Fisher's exact test for DAM-TE × AMT1 co-occurrence. Observed vs expected counts shown as mosaic plot; OR = 5.19, p = 0.028. The 15 DAM-TE species are highlighted (14 DIRS, 1 Gypsy).

---

### Figure 5 — Whole-genome duplication and TE dynamics in Mucoromycota

**(A)** Per-species paralog Ks distributions for Mucoromycota species. Gaussian mixture model (GMM) fits are overlaid in red; WGD-associated peaks (weight ≥10%, Ks 0.05–3.0) are marked with blue dashed lines. **(B)** TE density (bp TE / bp window) at synteny breakpoints (red) compared with 100 size-matched random windows per breakpoint from the same chromosome (blue). Wilcoxon p-value shown. **(C)** Scatter plot of inferred WGD time (x-axis, Mya) vs TE burst peak time (y-axis, Mya), coloured by TE class. Spearman ρ shown. The dashed diagonal represents simultaneous WGD and TE activation. **(D)** Mucoromycota species tree with WGD events marked (gold stars at inferred nodes), TE content bars per species, and synteny multiplicon counts.

---

### Figure 6 — Horizontal TE transfer across the eukaryotic tree of life

**(A)** Circular chord diagram of 224 HTT candidate events. Clades are arranged around the perimeter (coloured by clade as in Figure 1). Arcs connect donor and recipient clades; arc width is proportional to log10(candidate count) and arc colour indicates the dominant TE class. **(B)** Stacked bar plot of HTT candidates per donor–recipient clade pair, coloured by TE class. The 20 most frequently exchanging pairs are shown. **(C)** Strip plot of BLASTN percent identity for each candidate, stratified by TE class. The 85% identity threshold is shown as a dashed line. **(D)** Parsimony-reconciled minimum HTT events per TE class (total: 11 events). **(E)** Locus conservation at HTT sites: proportion of candidates with conserved (≥30% flanking homology), partial (10–30%), or non-conserved (<10%) insertion loci. **(F)** Inter-clade heatmap of HTT candidate counts (colour scale: YlOrRd; cells = number of transfer candidates between each clade pair).

---

## PART III — SOURCE DATA

Source data are provided as tab-separated value (`.tsv`) files generated by the analysis scripts. All files are deposited in the EDF-TE Zenodo repository (DOI: _to be assigned_) and organised by analysis section. Column descriptions follow below.

---

### Source Data 1 — Genomic features and species metadata

**File:** `source_data/SD1_genomic_features.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code (format: CLADE_SpeciesID) |
| `clade` | Broad clade tag (2-letter code, see Methods §1) |
| `genome_size_mb` | Total assembly size in megabases |
| `gc_pct` | GC content (%) computed over unmasked assembly |
| `gene_count` | Number of protein-coding gene models |
| `mean_gene_length_bp` | Mean length of gene models (bp, exons only) |
| `intron_density` | Mean number of introns per gene |
| `n50_mb` | N50 scaffold length (Mb) |

Derived from: `section1_genome_phylogenetics/02_genomic_features.py`.

---

### Source Data 2 — TE content per species and superfamily

**File:** `source_data/SD2_te_superfamily_content.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code |
| `clade` | Clade tag |
| `superfamily` | TE superfamily name (Wicker et al. 2007 nomenclature) |
| `te_class` | TE class (LTR, DNA, LINE, SINE, Rolling Circle, Penelope, Unclassified) |
| `n_copies` | Number of annotated TE copies |
| `bp_total` | Total base pairs covered by this superfamily |
| `pct_genome` | Percentage of genome covered |

Derived from: EarlGrey v1.6 `*_summaryFiles/*familyLevelCount.txt` outputs, processed by `section2_te_landscape/02_te_content_cladogram.py`.

---

### Source Data 3 — TE landscape ordination coordinates

**File:** `source_data/SD3_pcoa_coordinates.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code |
| `clade` | Clade tag |
| `pcoa_bc_1` | Bray-Curtis PCoA axis 1 (Hellinger-transformed; variance explained = 24.9%) |
| `pcoa_bc_2` | Bray-Curtis PCoA axis 2 (variance explained = 13.5%) |
| `pcoa_ait_1` | Aitchison (CLR) PCoA axis 1 |
| `pcoa_ait_2` | Aitchison (CLR) PCoA axis 2 |
| `pcoa_jac_1` | Jaccard P/A PCoA axis 1 |
| `pcoa_jac_2` | Jaccard P/A PCoA axis 2 |
| `shannon_h` | Shannon diversity H over superfamily-level proportions |

Derived from: `section2_te_landscape/03_te_landscape_ordination.py`.

---

### Source Data 4 — Kimura divergence and recency index

**File:** `source_data/SD4_kimura_divergence.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code |
| `clade` | Clade tag |
| `te_class` | TE class |
| `kimura_bin` | Kimura K2P bin midpoint (0.005 step) |
| `bp_in_bin` | Total base pairs of TE copies in this bin |
| `recency_index` | Fraction of species TE bp with Kimura < 0.10 (per-class) |

Derived from: RepeatMasker `.align` outputs processed by `section2_te_landscape/04_kimura_divergence.py`.

---

### Source Data 5 — TE methylation and epigenetics

**File:** `source_data/SD5_te_methylation.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code |
| `te_copy_id` | TE copy identifier (chrom:start-end:family) |
| `te_class` | TE class |
| `kimura_div` | Kimura divergence of this TE copy |
| `n_cpg_covered` | Number of CpG sites with ≥5× coverage |
| `pct_methylated` | Percentage of covered CpGs that are methylated |
| `is_methylated` | Boolean: ≥50% methylated CpGs |
| `h3k9me3_enrichment` | Mean H3K9me3 signal (RPKM) over TE copy body |
| `nearest_gene_rpkm` | RNA-seq RPKM of nearest gene |
| `chromatin_category` | One of: unmethyl_noH3K9, methyl_noH3K9, unmethyl_H3K9me3, methyl_H3K9me3 |

Derived from: `section3_te_epigenetics/02_5mc_te_age_analysis.py` and `03_chromatin_proximity_analysis.py`.

---

### Source Data 6 — TE co-option and domain fusions

**File:** `source_data/SD6_te_cooption.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code |
| `orthogroup` | OrthoFinder orthogroup ID |
| `gene_id` | Gene model identifier |
| `te_domain_ipr` | InterPro accession of TE-derived domain |
| `te_domain_name` | Domain name |
| `diamond_pident` | DIAMOND BLASTX identity to TE consensus (%) |
| `diamond_evalue` | DIAMOND e-value |
| `protein_length_aa` | Protein length (amino acids) |
| `co_option_class` | Classification: co-opted (no full-length TE ORF) or TE-like |

Derived from: `section4_te_cooption/01_cooption_orthogroup_pipeline.py`.

**File:** `source_data/SD6b_domain_fusions.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code |
| `gene_id` | Gene model identifier |
| `fusion_class` | Fusion type (FLYWCH-MULE, CotH-DDE, NACHT-RT, PARP-Chromo) |
| `te_ipr` | TE-derived domain IPR accession |
| `partner_ipr` | Host partner domain IPR accession |

Derived from: `section4_te_cooption/04_domain_fusion_analysis.py`.

---

### Source Data 7 — DAM methyltransferase and AMT1 co-occurrence

**File:** `source_data/SD7_dam_te_amt1.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code |
| `gene_id` | DAM-domain gene identifier |
| `hmmer_evalue` | HMMER e-value against DAM HMM (Pfam PF05869) |
| `te_associated` | Boolean: gene within 2 kb of a TE copy |
| `te_family` | Nearest TE family (if te_associated = True) |
| `te_class` | TE class of nearest element |
| `amt1_present` | Boolean: species has ≥1 AMT1-class methyltransferase (IPR014885) |

Fisher's exact test result: OR = 5.19, p = 0.028 (DAM-TE × AMT1 presence). Derived from: `section4_te_cooption/05_dam_te_methyltransferase.py`.

---

### Source Data 8 — WGD detection and TE copy number

**File:** `source_data/SD8_wgd_summary.tsv`

| Column | Description |
|--------|-------------|
| `species` | Species code |
| `clade` | Clade tag |
| `n_pairs` | Number of paralogous Ks pairs used for GMM fitting |
| `n_gmm_k` | Number of GMM components selected by BIC |
| `n_peaks` | Number of WGD-associated peaks identified |
| `wgd_call` | Boolean: species assigned WGD event |
| `peak_ks` | Semicolon-delimited list of GMM peak Ks values |

Derived from: `section5_wgd_te_dynamics/01_wgd_detection.py`.

**File:** `source_data/SD8b_te_copy_evolution.tsv`

| Column | Description |
|--------|-------------|
| `superfamily` | TE superfamily |
| `mean_wgd` | Mean copy number in WGD species |
| `mean_nowgd` | Mean copy number in non-WGD species |
| `lfc` | log2(WGD / non-WGD) copy number fold change |
| `pval` | Mann-Whitney U p-value (two-sided) |

Derived from: `section5_wgd_te_dynamics/05_te_copy_number_evolution.py`.

---

### Source Data 9 — Horizontal TE transfer candidates

**File:** `source_data/SD9_htt_candidates.tsv`

| Column | Description |
|--------|-------------|
| `qseqid` | Query TE consensus identifier (format: SPECIES__family) |
| `sseqid` | Subject TE consensus identifier |
| `q_clade` | Donor clade |
| `s_clade` | Recipient clade |
| `pident` | BLASTN percent identity |
| `length` | BLASTN alignment length (bp) |
| `anomaly_score` | Expected − observed divergence (higher = more anomalous) |
| `median_ks` | Median Ks of single-copy orthologs between donor and recipient clade (transfer time proxy) |
| `conservation_fraction` | Fraction of 10 kb flanking sequence with ≥60% identity in recipient genome |
| `conservation_category` | conserved (≥30%), partial (10–30%), non-conserved (<10%) |
| `te_class` | TE class of the transferred element |

Summary: 224 candidate HTT events; 11 minimum parsimony-reconciled events. Derived from: `section6_htt/02_htt_candidate_detection.py`, `03_htt_ks_dating.py`, `04_htt_locus_conservation.py`, `05_htt_visualization.py`.

---

*All source data files are provided in tab-separated format with a header line. Missing values are encoded as `NA`. For questions regarding data access or analysis reproducibility, consult the pipeline documentation in `EDF-TE_complete/` at https://github.com/cLax8M/EDF-TE.*
