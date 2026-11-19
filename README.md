# EDF-TE
Code for the study of TE evolution in the Holomycota (EDF project).

## Repository structure

**Section 1** Comparative Genomics and Genome Evolution
    Phylogenomic reconstruction, genome features, pangenome dynamics and whole-genome synteny analysis across fungal species.
   
**Section 2** Transposable Element Landscape and Diversity
    TE annotation, compositional diversity, evolutionary dynamics, divergence, protein domain expansions and TE enrichment at mating-type loci.
    
**Section 3** Analysis of epigenetic mechanisms associated with transposable element (TE) silencing across fungal species, including:

* Conservation of DNA methyltransferases (DNMT1 and DIM2).
* Association between TE evolutionary age and DNA methylation (5mC).
* Small RNA (sRNA)-mediated TE silencing and its relationship with H3K9me3.
* TE association with heterochromatin marks (H3K9me3 and 5mC) and gene expression.
* Association between N6-methyladenine (6mA) and TE expression.
* Nucleosome positioning and chromatin organization at TE loci.
  
**Section 4** TE Functional Co-option and Regulatory Networks
    Identification of TE-derived genes, chimeric transcripts, TF binding, domain fusions, TE-associated methyltransferases and TF–TE–gene interaction networks.

**Section 5** Whole-Genome Duplication and TE-associated Genome Evolution
    WGD detection, synteny breakpoint analysis, TE insertion timing, WGD–TE associations and TE copy-number evolution.
    
**Section 6** Horizontal Transposable Element Transfer (HTT)
    Identification and evolutionary dating of candidate horizontal TE transfer events, genomic locus conservation and visualization.

## Data

Large data files (genomes, BAMs, TSVs, WIGs) are **not tracked** in this repository.
See `/legserv/Temp/Carlos/EDF-TE/data/` and `/legserv/Temp/Carlos/EDF-TE/results/` for all data.


## Dependencies

All R scripts require the `treewasEDFTE` conda environment (`TreeWAS/treewasEDFTE.yml`).
Nucleosome scripts require additional conda environments:
- `nucseq_pipeline`: sra-tools, bowtie, samtools, picard, danpos2, bbtools
- `nuclpred`: see `TE_analysis/NuclPred_v1/nucleosome_prediction_env.yml` + DNAshapeR (Bioconductor)

