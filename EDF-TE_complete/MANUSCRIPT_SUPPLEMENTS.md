# EDF-TE — Materials and Methods, Figure Legends, and Source Data

> **Notation.**  Species codes: CLADE_SpeciesID (e.g. MU_Mucci3 = *Mucor circinelloides* v3).
> `ln` = natural log; `log₂` = base-2 log.  Micromamba environments defined in `config.sh`.

---

## PART I — MATERIALS AND METHODS

### 1. Genome dataset and taxonomic sampling

We assembled a dataset of 320 eukaryotic genomes spanning ~1 Gyr of divergence
across 18 broad lineages: Holozoa (HO), Nucleariida (NU), Aphelida (AP),
Cryptomycetes (CR), Microsporidia (MI), Blastocladiomycetes (BL), Olpidiales (OL),
Chytridiomycetes (CH), Neocallimastigomycota (NE), Bacillariophyta (BS),
Entomophthorales (EN), Zoopagomycota (ZO), Kickxellomycota (KX), Glomeromycota (GL),
Mortierellomycota (MO), Mucoromycota (MU), Basidiomycota (BA), and Ascomycota (AS).
Seven Neocallimastigomycota species (CH_Anasp1, CH_Neocon1, CH_Neolan1, CH_Orpsp1,
CH_Pecora1, CH_Caecom1, CH_Pirfi3) share the JGI CH prefix but were assigned to
the NE clade for all biological analyses.

Assemblies were downloaded from JGI MycoCosm, NCBI RefSeq, and ENA.  Scaffolds
shorter than 200 bp were discarded.  Protein annotations (GFF/GFF3) and functional
domain annotations (InterProScan 5.67-99.0, Pfam 36.0, PRINTS, PANTHER 18.0,
TIGRFAM, CDD) were obtained from JGI or recomputed locally.  Three Mucoromycota
assemblies (MU_Mucci3, MU_Rhimi59, MU_Phybl4) were generated from PacBio HiFi CCS
reads with Flye 2.9 and assessed with BUSCO 5.8.2 (fungi_odb10 dataset).

---

### 2. Species tree inference

Single-copy orthogroups were identified with OrthoFinder 3.1
(`-M msa -A mafft -T fasttree -S diamond`) across 238 genomes.  Markers were
retained if they satisfied all three criteria simultaneously:

1. **Taxon occupancy** ≥ 40 % (≥ 95 of 238 genomes present).
2. **Copy number** ≤ 2 per species (to limit hidden paralogy).
3. **Quality** after alignment and trimming: post-trimming length ≥ 60 amino acid
   columns; parsimony-informative sites (PIS) ≥ 25; relative composition variability
   (RCV; PhyKIT v2.0) ≤ 0.5.  RCV is defined as
   $\text{RCV} = \frac{1}{4L}\sum_i |f_i - \bar{f}|$,
   where $f_i$ is the frequency of amino acid *i* and *L* is the alignment length.

Each marker was aligned with MAFFT v7.520 (`--auto`) and trimmed with trimAl v1.4.1
(`-automated1`).  Per-gene trees under LG+G4 (IQ-TREE 2.3.6) were used to detect
long-branch paralogs: taxa with branch length > 3.5 standard deviations above the
per-gene mean were removed.  The remaining markers were concatenated into a
supermatrix (one partition per gene) and the final tree was inferred under
LG+C60+F+G4 with 1,000 ultrafast bootstrap replicates (`-B 1000 -bnni`).  A
guide-tree pass under LG+F+G4 preceded the final inference.  The tree was rooted
on Holozoa.  All steps are in
`section1_genome_phylogenetics/01_species_tree_pipeline.sh`.

---

### 3. Genomic feature characterisation

Six per-genome metrics were extracted with
`section1_genome_phylogenetics/02_genomic_features.py`:

| Metric | Calculation |
|--------|-------------|
| Genome size (Mb) | Sum of scaffold lengths ÷ 10⁶; N bases excluded |
| GC content (%) | (G + C) / (A + T + G + C) × 100 |
| Gene count | Number of protein-coding loci (GFF type = gene or mRNA) |
| Mean gene length (bp) | Mean of (tx_end − tx_start + 1) for longest isoform per gene |
| Intron density | Total intron features / total gene features in GFF |
| N50 (Mb) | Scaffold length at which ≥ 50 % of assembled bases are in scaffolds ≥ this length |

Pan-genome accumulation curves (core, shell, unique gene counts as a function of
sampling depth) were estimated by 1,000 random permutations of species-addition
order applied to the OrthoFinder presence/absence matrix
(`section1_genome_phylogenetics/03_pangenome_analysis.py`).  Core genes are
present in all sampled species; unique (cloud) genes in exactly one; shell
otherwise.

---

### 4. Synteny analysis

Pairwise macrosynteny for Mucoromycota species pairs was computed with SynChro v2.0.
Protein sequences were compared by DIAMOND BLASTP v2.1.8 (identity ≥ 30 %,
e-value ≤ 10⁻⁵, coverage ≥ 50 %, `-k 5`).  Reciprocal best hits anchored
syntenic blocks (minimum 5 anchors; maximum permissible inter-anchor gene gap = 25).
Syntenic block BED coordinates were used downstream for breakpoint analysis.  The
circular chord diagram (`section1_genome_phylogenetics/05_circos_visualization.py`)
shows each pairwise anchor count as a chord; chord width is proportional to the
number of anchor pairs.

---

### 5. TE annotation with EarlGrey

TEs were annotated in all 320 genomes with EarlGrey v1.6
(`section2_te_landscape/01_te_annotation.sh`).  EarlGrey runs sequentially:

1. **RepeatModeler2 v2.0.5** — de novo consensus discovery (RECON + RepeatScout
   + LTRharvest + LTR_retriever).
2. **Structural filtering** — consensi < 80 bp or > 40 kb discarded; TIRs and
   TSDs verified for DNA transposons; LTR boundaries verified with LTR_retriever.
3. **Phylogenetic classification** — each consensus classified against Dfam 3.7
   with nhmmer (HMMER 3.3.2, E-value ≤ 10⁻⁵); unclassified consensi passed to
   RepeatClassifier.
4. **RepeatMasker 4.1.5 re-annotation** — curated per-species library used to
   produce `*.filteredRepeats.bed` and `*.align` outputs.

SLURM array: 16 CPUs, 64 GB RAM, partition `normal.168h` per species.

After per-species annotation, all consensi were merged and clustered at 80 % identity
with MMseqs2 v14.7564 (`easy-cluster --min-seq-id 0.80 -c 0.80 --cov-mode 1`).
Sequences with DIAMOND BLASTX hits to RefSeq proteins (identity ≥ 40 %, e-value
≤ 10⁻⁵) were removed as host-derived.  TE superfamily nomenclature follows
Wicker et al. (2007).  MI_Encro1 and MI_Ordco1 (zero TE content) were excluded
from ordination analyses.

---

### 6. Kimura divergence landscapes and Recency Index

RepeatMasker produces a pairwise alignment (`.align`) of each TE insertion against
its consensus.  For every copy we compute the Kimura two-parameter (K2P) distance:

$$K = -\frac{1}{2}\ln(1 - 2p - q) - \frac{1}{4}\ln(1 - 2q)$$

where *p* = (transitions) / (aligned sites) and *q* = (transversions) / (aligned
sites), read directly from the per-copy substitution matrix in the `.align` file.
Copies with ≤ 10 aligned sites after gap-stripping are discarded.

Kimura landscapes are stacked histograms with bin width 0.005 (200 bins, K ∈ [0, 1]).
Each bin y-value is the sum of base pairs from copies in that bin, stacked by TE
class (`section2_te_landscape/04_kimura_divergence.py`).

The **Recency Index (RI)** for species *s*:

$$RI_s = \frac{\displaystyle\sum_{\substack{i \in s \\ K_i < 0.10}} bp_i}{\displaystyle\sum_{i \in s} bp_i}$$

where *i* indexes all annotated TE copies in *s*, $bp_i$ is the copy footprint
(BED end − start), and $K_i$ is its K2P distance.  The threshold K < 0.10
captures copies that have accumulated fewer substitutions than expected under
neutral decay over ~50 Myr at typical fungal substitution rates, i.e. copies
inserted relatively recently.

---

### 7. TE landscape ordination

The composition matrix **X** (236 species × 234 superfamilies) has entries
$x_{ij}$ = base pairs of superfamily *j* in species *i* (from
`section2_te_landscape/03_te_landscape_ordination.py`).

#### 7.1 Bray-Curtis PCoA (Hellinger-transformed)

Relative abundances: $p_{ij} = x_{ij} / \sum_j x_{ij}$.
Hellinger transformation: $h_{ij} = \sqrt{p_{ij}}$.

$$BC_{ab} = \frac{\sum_j |h_{aj} - h_{bj}|}{\sum_j (h_{aj} + h_{bj})}$$

PCoA (classical MDS) applied via `sklearn.manifold.MDS(metric=True)`.
Variance explained by axis *k*: $\lambda_k / \sum_{\lambda>0} \lambda$.
**PCo1 = 24.9 %, PCo2 = 13.5 %.**

#### 7.2 Aitchison PCoA (CLR-transformed)

Laplace pseudocount of 1 added to all entries.  Centre log-ratio (CLR):

$$c_{ij} = \ln\!\left(\frac{x_{ij}+1}{g_i}\right), \quad g_i = \exp\!\!\left(\frac{1}{p}\sum_j \ln(x_{ij}+1)\right)$$

where *p* = 234 superfamilies and $g_i$ is the geometric mean of the pseudocount-
augmented composition.  Euclidean distances on CLR vectors equal Aitchison distances.
PCoA applied to the resulting distance matrix.

#### 7.3 Jaccard P/A PCoA

Binary matrix $b_{ij} = \mathbf{1}[x_{ij} > 0]$.

$$J_{ab} = 1 - \frac{\sum_j b_{aj} b_{bj}}{\sum_j \max(b_{aj}, b_{bj})}$$

#### 7.4 ANOSIM

For each ordination, ANOSIM tests group separation (999 permutations):

$$R = \frac{\bar{r}_B - \bar{r}_W}{N(N-1)/4}$$

where $\bar{r}_B$ = mean rank dissimilarity between clades, $\bar{r}_W$ = within
clades, *N* = total species.  **Bray-Curtis R = 0.23, p < 0.001; Jaccard R = 0.18,
p < 0.001.**

#### 7.5 Shannon diversity

$$H_i = -\sum_j p_{ij} \log_2 p_{ij}, \quad p_{ij} = \frac{x_{ij}}{\sum_j x_{ij}}$$

The convention $0 \log_2 0 = 0$ is applied.  H is in bits.  Species with zero
total TE content are assigned H = 0.

---

### 8. DNMT conservation and 5-methylcytosine analysis

#### 8.1 DNMT orthogroup identification

A gene was classified as DNMT1/DIM2 if its InterProScan record contained **both**
IPR005105 (DNA cytosine methyltransferase catalytic domain, Pfam PF00145) **and**
IPR022686 (PWWP domain, Pfam PF00855).  Per-species copy numbers from OrthoFinder
were projected on the species tree with ete3 v3.1.3
(`section3_te_epigenetics/01_dnmt_conservation_tree.py`).

#### 8.2 WGBS processing

Paired-end WGBS reads were quality-trimmed with Trim Galore v0.6.10 (Cutadapt
backend, `--quality 20 --phred33 --paired`) and mapped with Bismark v0.24.2
(`--genome --non_directional --score_min L,0,-0.6`).  Duplicates removed with
`deduplicate_bismark`.  CpG methylation calls extracted with
`bismark_methylation_extractor --CpG --paired-end --comprehensive --cytosine_report`.

#### 8.3 Per-copy methylation

CpG cytosine calls were intersected with TE BED coordinates using
`bedtools intersect -wa -wb -sorted`.  Per-copy methylation level:

$$m_i = \frac{\text{methylated CpG calls in copy } i \text{ at coverage} \geq 5\times}{\text{total CpG calls in copy } i \text{ at coverage} \geq 5\times}$$

Copies require ≥ 3 covered CpGs; copies with $m_i \geq 0.50$ are classified
methylated.  Spearman ρ between $m_i$ and $K_i$ computed across all covered copies.
(`section3_te_epigenetics/02_5mc_te_age_analysis.py`).

---

### 9. H3K9me3 ChIP-seq and chromatin categories

ChIP-seq reads trimmed with Trim Galore v0.6.10, mapped with Bowtie2 v2.5.1
(`--very-sensitive -X 2000`).  Peak calling and per-base signal with MACS2 v2.2.9.1
(`callpeak --nomodel --extsize 200 --bdg -B -q 0.05`).  H3K9me3 enrichment per TE
copy = mean bedgraph value over copy coordinates, depth-normalised to RPKM.

TE copies classified into four chromatin categories (5mC threshold $m_i \geq 0.50$;
H3K9me3 threshold RPKM ≥ 1.0):

| Category | 5mC | H3K9me3 |
|----------|-----|---------|
| unmethyl_noH3K9 | − | − |
| methyl_noH3K9 | + | − |
| unmethyl_H3K9me3 | − | + |
| methyl_H3K9me3 | + | + |

RNA-seq reads mapped with STAR v2.7.11a.  Nearest-gene expression quantified with
featureCounts (Subread v2.0.6, `-s 0 -p -T 8`), reported as RPKM.  Distance to
nearest gene: `bedtools closest -d`
(`section3_te_epigenetics/03_chromatin_proximity_analysis.py`).

---

### 10. Nucleosome positioning — meta-TE profiles

MNase-seq mapped with Bowtie2 v2.5.1 (`--very-sensitive`), processed with
DANPOS3 v3.1.1 (`dpos --paired`), output per-base dyad probability (Wig).

For each TE class, meta-TE profiles were generated
(`section3_te_epigenetics/05_nucleosome_te_analysis.py`): the ±2,000 bp region
around each copy mid-point was extracted from the Wig signal, scaled to 150 equal
bins by linear interpolation, then averaged across all copies.  Bins overlapping
fewer than 10 copies excluded.  TE body boundaries within the scaled profile are
proportional to the ratio (copy length) / (window length).

---

### 11. 6mA detection and small RNA silencing

#### 11.1 N6-methyladenine via PacBio kinetics

PacBio CLR reads mapped with `pbmm2 align --preset SUBREAD --sort`.  IPD ratios
computed with `ipdSummary` (SMRTLink v12.0).  Positions with IPD ratio ≥ 5.0 and
coverage ≥ 25× called as 6mA.  TE copies with ≥ 1 called 6mA site are 6mA-positive.
Nearest-gene RNA-seq RPKM compared (6mA-positive vs. negative) by Wilcoxon rank-sum
test (`section3_te_epigenetics/04_6mA_expression_analysis.py`).

#### 11.2 Small RNA mapping

sRNA reads trimmed with Cutadapt v4.6 (`--minimum-length 18 --maximum-length 32`),
size-filtered to 20–25 nt.  Mapped to the TE consensus library with Bowtie v1.3.1:

```
bowtie -v 2 -M 1 --best --strata -p 8 -x <te_library_index> -q <reads.fq>
```

`-M 1` reports one random alignment for multimappers.  RPM per TE family:

$$\text{RPM}_j = \frac{\text{mapped reads to family } j}{\text{total mapped reads}} \times 10^6$$

Spearman ρ between per-family RPM and mean H3K9me3 RPKM tested per species
(`section3_te_epigenetics/06_srna_silencing_pipeline.py`).

---

### 12. TE co-option analysis

#### 12.1 Domain-based co-option screen

All OrthoFinder orthogroups were screened for the presence of TE-associated
InterPro domains in ≥ 1 member protein.  Thirteen TE-derived domain families were
used as queries (IPR001584 RVT_1, IPR002104 DDE transposase, IPR021967 IS200/IS605
transposase, IPR001878 CCHC Zn-finger, IPR002510 integrase catalytic domain,
IPR000477 reverse transcriptase, IPR004875 IS3/IS481 transposase, IPR006538 IS5
transposase, IPR002951 IS4 transposase, IPR011010 DNA breaking-rejoining enzyme,
IPR013762 integrase-like N-term, IPR000380 MULE transposase, IPR003086 SWIM Zn-finger).
Orthogroups with ≥ 1 hit were validated by DIAMOND BLASTX v2.1.8 (identity ≥ 30 %,
e-value ≤ 10⁻⁵, coverage ≥ 50 %) against the cross-species TE library.
A protein is **co-opted** (not a TE-encoded ORF) if: (i) protein length < 800 aa,
(ii) ≥ 1 exon ≥ 100 aa in the GFF, AND (iii) no full-length hit (query coverage
≥ 80 %) to any TE consensus (`section4_te_cooption/01_cooption_orthogroup_pipeline.py`).

#### 12.2 ChimeraTE analysis

Chimeric gene–TE transcripts were detected with ChimeraTE v1.0
(`section4_te_cooption/02_chimte_analysis.py`) applied to polyA-selected RNA-seq.
ChimeraTE was run in `--mode genome` using the EarlGrey-filtered TE BED, the species
GFF3, and a coordinate-sorted BAM from STAR v2.7.11a.

A read pair is chimeric if:
- One mate maps entirely within a TE copy (≥ 1 bp overlap with the TE BED annotation).
- The other mate maps to an annotated exon (≥ 1 bp overlap with the gene GFF).
- ≥ 3 read pairs support the same TE copy × gene junction.
- The TE copy and gene are within 10 kb.

Chimera calls from all species were aggregated by TE class and by OrthoFinder
orthogroup.  **Shared chimeras** (same TE family × same orthogroup in ≥ 2 species)
were identified by inner join on (te_family, orthogroup_id) columns.

#### 12.3 TF binding sites in TE sequences

FIMO v5.5.5 scanned all TE consensus sequences for TFBS motifs from JASPAR 2024
Core Fungi:

```
fimo --oc <out_dir> --thresh 1e-4 --qv-thresh 0.05 <motifs.meme> <te_consensi.fa>
```

For each TF × TE class pair, enrichment was tested with Fisher's exact test on
the table: (TE copies with ≥ 1 FIMO hit vs. no hit) × (TE class of interest vs.
all others).  BH correction applied.  HMGb3 × LINE/I-Jockey: p = 0.00275,
interpreted cautiously given n = 4 species
(`section4_te_cooption/03_tf_te_interaction_analysis.py`).

#### 12.4 Multi-domain fusion proteins

Proteins with both a TE-derived domain and a host functional domain were identified
from InterProScan output
(`section4_te_cooption/04_domain_fusion_analysis.py`):

| Fusion class | TE domain | Host domain |
|---|---|---|
| FLYWCH–MULE | IPR007527 (FLYWCH Zn-finger) | IPR002992 (MULE transposase) |
| CotH–DDE | IPR033853 (CotH-like) | IPR002104 (DDE transposase) |
| NACHT–RT | IPR007111 (NACHT ATPase) | IPR000477 (Reverse transcriptase) |
| PARP–Chromo | IPR012317 (PARP catalytic) | IPR000953 (Chromodomain) |

A protein is a fusion if **both** IPR accessions appear in its InterProScan record.

#### 12.5 DAM methyltransferases embedded in TEs

DAM-domain proteins identified by HMMER v3.3.2 hmmsearch against Pfam PF05869
(E-value ≤ 10⁻⁵, `--noali`).  A gene is **TE-embedded** if its genomic coordinates
fully overlap a TE BED annotation (`bedtools closest -d` distance = 0).
AMT1-class methyltransferases defined by IPR014885.  Fisher's exact test on
(DAM-TE present/absent) × (AMT1 present/absent) across 320 species:
**OR = 5.19, p = 0.028** (two-tailed).
Fifteen species harbour DAM-TE genes: 14 DIRS, 1 Gypsy
(`section4_te_cooption/05_dam_te_methyltransferase.py`).

---

### 13. Whole-genome duplication

#### 13.1 Paralog Ks distributions

Intraspecific paralogous gene pairs from OrthoFinder 3.1.  Protein alignments:
MAFFT v7.520 (`--auto`); codon-level alignment: PAL2NAL v14.  Ks computed with
PAML v4.10.10 yn00 module.  Filter: 0 < Ks ≤ 3.0.

#### 13.2 GMM fitting

Fitted with `sklearn.mixture.GaussianMixture`
(`covariance_type='full'`, `n_init=5`, `max_iter=200`, `random_state=42`), *k* ∈ {1,…,5}.
Model selected by lowest BIC:

$$\text{BIC}(k) = k \ln n - 2 \ln \hat{L}$$

A component is a **WGD peak** if weight ≥ 0.10 and mean $\mu \in [0.05, 3.0]$
(`section5_wgd_te_dynamics/01_wgd_detection.py`).

#### 13.3 TE density at breakpoints

Breakpoints = both termini of each SynChro syntenic block ± 5 kb.  TE density
per window:

$$d_W = \frac{\displaystyle\sum_i \text{overlap}(\text{TE}_i, W)}{\text{len}(W)}$$

Null: 100 randomly sampled non-syntenic windows of the same length on the same
chromosome.  Wilcoxon signed-rank test per species
(`section5_wgd_te_dynamics/02_te_breakpoint_enrichment.py`).

#### 13.4 TE age vs WGD timing

WGD Ks peaks translated to Mya using per-lineage synonymous substitution rates
from calibrated tree nodes.  TE burst time = mode of Kimura distribution (KDE,
Gaussian kernel, Scott's bandwidth rule).  Spearman ρ tested between WGD and TE
burst times across species with WGD calls
(`section5_wgd_te_dynamics/03_te_age_wgd_correlation.py`).  Per-superfamily copy
number compared (WGD vs non-WGD) by Mann-Whitney U, BH-corrected
(`section5_wgd_te_dynamics/05_te_copy_number_evolution.py`).

---

### 14. Horizontal TE transfer

#### 14.1 Per-clade TE libraries

EarlGrey consensi concatenated per clade (headers prefixed `SPECIES__family`),
then clustered at 95 % identity with MMseqs2 `easy-cluster`
(`--min-seq-id 0.95 -c 0.80 --cov-mode 1`).  Re-classified with RepeatClassifier
(RepeatModeler 2.0.5) (`section6_htt/01_htt_phylo_prep.sh`).

#### 14.2 Anomaly-score detection

All-against-all BLASTN between every clade-library pair
(identity ≥ 85 %, length ≥ 100 bp, word_size 11, E-value ≤ 10⁻⁵).  For each hit:

$$\text{anomaly} = d_{\text{expected}}(X,Y) - d_{\text{observed}}(A,B)$$

$d_{\text{expected}}$ = tip-to-tip branch length in calibrated species tree;
$d_{\text{observed}} = 1 - \text{pident}/100$.
Threshold: anomaly ≥ 0.10 (TE copies ≥ 10 % more similar than expected from
vertical inheritance) (`section6_htt/02_htt_candidate_detection.py`).

#### 14.3 Ks-based dating

Median Ks of single-copy OrthoFinder orthologs between clade representatives
(same yn00 pipeline as §13.1) provides a divergence time lower bound for each
candidate donor–recipient clade pair
(`section6_htt/03_htt_ks_dating.py`).

#### 14.4 Insertion locus conservation

For each candidate, 10 kb flanks around the TE copy (donor genome) were extracted
with `bedtools getfasta` and compared to the recipient genome by BLASTN
(identity ≥ 60 %, length ≥ 50 bp).  Conservation fraction:

$$f = \frac{\displaystyle\sum_{\text{hits}} \text{aln\_len}}{\text{total flank length}}, \quad f \leq 1$$

Classified: conserved (f ≥ 0.30), partial (0.10–0.30), non-conserved (f < 0.10)
(`section6_htt/04_htt_locus_conservation.py`).

#### 14.5 Parsimony reconciliation

224 candidate events mapped onto the species tree and reduced by Dollo parsimony,
yielding **11 minimum independent HTT events**
(`section6_htt/05_htt_visualization.py`).

---

### 15. Statistical analyses and software

All tests two-tailed unless stated.  Multiple-testing corrected with
Benjamini-Hochberg (FDR < 0.05) where ≥ 5 simultaneous hypotheses are tested.
Visualisation: Python 3.11, matplotlib 3.8, seaborn 0.13; R 4.4, ggplot2 4.0.
Phylogenetic trees: ete3 v3.1.3.  Networks: networkx 3.2 → Cytoscape 3.10.2.

Full pipeline: **https://github.com/cLax8M/EDF-TE** (`EDF-TE_complete/`).

---

## PART II — FIGURE LEGENDS

### Figure 1 — Phylogenomic framework and genomic features

**(A)** ML species tree of 238 eukaryotic genomes (LG+C60+F+G4, 1,000 ultrafast
bootstrap).  Branches coloured by clade.  Node symbols: filled = bootstrap ≥ 95,
open < 95.  **(B)** Heatmap of six per-species genomic metrics, z-score normalised
per column, ordered as in (A).  **(C)** Pan-genome accumulation curves for
Mucoromycota (mean ± 1 SD, 1,000 random addition orders).  **(D)** Circular chord
diagram of macrosyntenic blocks among Mucoromycota species.  Chord width ∝ number
of syntenic anchor pairs; chord colour = clade of anchor-donating genome.

---

### Figure 2 — TE landscape diversity

**(A)** Cladogram of 320 species (pruned from Figure 1A) with stacked horizontal
bars showing TE class composition (% genome).  Adjacent colour track: Shannon
diversity H (bits) per species.  **(B)** Bray-Curtis PCoA (Hellinger-transformed,
236 species, 234 superfamilies).  Points coloured by clade; convex hulls per clade.
PCo1 = 24.9 %, PCo2 = 13.5 % of variance.  ANOSIM R = 0.23, p < 0.001 (999
permutations).  **(C)** Kimura divergence landscapes (stacked histograms, bin width
0.005) for representative species per clade.  RI value annotated.  **(D)** Shannon
H per clade (phylogenetically ordered boxplots; individual species as jittered
points; whiskers = 1.5 × IQR).

---

### Figure 3 — Epigenetic regulation of TEs

**(A)** Species tree with DNMT1/DIM2 presence/absence (IPR005105 + IPR022686)
as tip symbols; node bars = orthogroup copy number.  **(B)** Violin plot of 5mC
per TE copy (WGBS; ≥ 3 CpGs at ≥ 5×) per Kimura divergence quartile; Spearman ρ.
**(C)** Nearest-gene RPKM per chromatin category (boxplots; Kruskal-Wallis +
Dunn p-values).  **(D)** Meta-TE nucleosome dyad-probability profiles (DANPOS3;
±2 kb, 150 bins; TE body shaded grey).  **(E)** sRNA (20–25 nt) RPM mapping to
top-15 TE superfamilies.  **(F)** Per-family sRNA RPM (y) vs mean H3K9me3 RPKM
(x); Spearman ρ; linear regression line.

---

### Figure 4 — TE co-option

**(A)** Bubble plot of co-opted orthogroups per species; bubble area = number of
species sharing the OG; colour = dominant TE domain class.  **(B)** Heatmap:
fraction of species per clade with ≥ 1 co-opted gene per TE domain family
(18 clades × 13 domain families).  **(C)** Stacked bars: proportion of chimeric
TE–gene transcripts (ChimeraTE) per TE class per species.  **(D)** TF × TE class
enrichment heatmap (−log₁₀ BH-FDR; Ward-linkage clustering; black outline = FDR
< 0.05).  **(E)** Fraction of species per clade with ≥ 1 gene per multi-domain
fusion class.  **(F)** Mosaic plot: DAM-TE × AMT1 contingency (Fisher OR = 5.19,
p = 0.028); DAM-TE species annotated by TE class.

---

### Figure 5 — WGD and TE dynamics in Mucoromycota

**(A)** Colinear paralog pair count (≤ 200 scaffolds) per Mucoromycota species,
ranked; coloured by WGD classification.  **(B)** TE density at synteny breakpoints
(red) vs 100 random windows per breakpoint (blue); Wilcoxon p-value shown.
**(C)** TE class composition (bp) at breakpoints vs genome background for four
significantly enriched species.  **(D)** Mucoromycota species tree with gene
duplication rate per branch (moving-average smoothed) and WGD events (gold stars).

---

### Figure 6 — Horizontal TE transfer

**(A)** Circular chord diagram of 224 HTT candidates; arc width ∝ log₁₀(count);
arc colour = dominant TE class.  **(B)** Stacked bar of HTT candidates per clade
pair (top-20), by TE class.  **(C)** BLASTN pident per candidate by TE class;
dashed line at 85 %.  **(D)** Minimum parsimony HTT events per TE class (total 11).
**(E)** HTT locus conservation category proportions.  **(F)** Inter-clade HTT
candidate count heatmap (YlOrRd scale).

---

## PART III — SOURCE DATA

Source data are provided as a single Excel workbook **`EDF-TE_source_data.xlsx`**
generated by `EDF-TE_complete/create_source_data.py`.

**Each worksheet = one quantitative figure panel.**  Each sheet contains the
exact numerical values used to draw that panel — enough to reproduce the figure
from scratch.  Row 1 = column headers; Row 2 = unit/description; data from Row 3.
Missing values = blank cell.  Run `python create_source_data.py --config config.sh`
to regenerate from pipeline outputs.

### Worksheet index

| Sheet | Panel | Content |
|-------|-------|---------|
| Fig1B_genomic_features | Fig 1B | species, clade + 6 raw metrics + 6 z-score normalised values |
| Fig1C_pangenome | Fig 1C | n_species_sampled, core_mean, core_sd, shell_mean, shell_sd, unique_mean, unique_sd |
| Fig1D_synteny | Fig 1D | sp_a, sp_b, n_anchor_pairs, n_multiplicons, n_segments |
| Fig2A_te_content | Fig 2A | species, clade, te_class, bp_total, pct_genome, shannon_h |
| Fig2B_pcoa | Fig 2B | species, clade, pcoa_bc_1, pcoa_bc_2, pcoa_ait_1, pcoa_ait_2, pcoa_jac_1, pcoa_jac_2, shannon_h |
| Fig2C_kimura | Fig 2C | species, clade, te_class, kimura_bin_midpoint, bp_in_bin, recency_index |
| Fig2D_shannon | Fig 2D | species, clade, shannon_h |
| Fig3A_dnmt | Fig 3A | species, clade, dnmt_present (0/1), orthogroup_copynum |
| Fig3B_methylation | Fig 3B | species, te_copy_id, te_class, kimura_div, kimura_quartile, pct_methylated, n_cpg_covered |
| Fig3C_chromatin | Fig 3C | species, te_copy_id, te_class, chromatin_category, nearest_gene_rpkm |
| Fig3D_nucleosome | Fig 3D | te_class, bin_index (1–150), mean_dyad_prob, n_copies_contributing |
| Fig3E_srna_rpm | Fig 3E | species, te_superfamily, rpm, total_mapped_reads |
| Fig3F_srna_h3k9 | Fig 3F | species, te_superfamily, rpm, mean_h3k9me3_rpkm |
| Fig4B_cooption_heatmap | Fig 4B | clade, te_domain_ipr, te_domain_name, frac_species_with_cooption |
| Fig4C_chimte | Fig 4C | species, clade, te_class, n_chimeric_events |
| Fig4D_tf_te | Fig 4D | tf_motif, te_class, neg_log10_fdr, fisher_or |
| Fig4E_fusions | Fig 4E | clade, fusion_class, frac_species_with_fusion |
| Fig4F_dam_amt1 | Fig 4F | species, clade, dam_te_present (0/1), amt1_present (0/1), te_class_dam, dam_evalue |
| Fig5A_paralogs | Fig 5A | species, n_colinear_paralog_pairs, n_multiplicons, wgd_classification |
| Fig5B_breakpoints | Fig 5B | species, n_breakpoints, obs_te_bp_per_window, bg_te_bp_per_window, bg_std, log2_fold_change, p_value_permutation |
| Fig5C_te_class_bp | Fig 5C | species, te_class, bp_at_breakpoints, bp_genome_background, log2_bp_vs_genome |
| Fig5D_dup_rate | Fig 5D | node_id, node_type, n_duplications, dup_rate_per_desc_leaf, descendant_leaves |
| Fig6A_htt_chord | Fig 6A | q_clade, s_clade, n_candidates, dominant_te_class |
| Fig6B_htt_class | Fig 6B | q_clade, s_clade, te_class, n_candidates |
| Fig6C_htt_identity | Fig 6C | q_clade, s_clade, te_class, pident (one row per candidate) |
| Fig6D_min_events | Fig 6D | te_class, n_minimum_events |
| Fig6E_locus_cons | Fig 6E | conservation_category, n_candidates, frac_total |
| Fig6F_htt_heatmap | Fig 6F | q_clade (rows), s_clade (columns), cell = n_candidates (wide format) |
