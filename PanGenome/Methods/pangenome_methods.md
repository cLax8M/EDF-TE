# Pangenome Analysis – Methods

## Pangenome and gene family classification

To characterise the gene repertoire diversity across fungal lineages, we performed a pangenome analysis using the orthogroup assignments produced by OrthoFinder v2 (Emms & Kelly 2019). Gene families were defined as flat orthogroups from the OrthoFinder output file `Orthogroups.GeneCount.tsv`, which reports the number of genes per species per orthogroup across the full dataset of 238 fungal and non-fungal outgroup proteomes. A binary presence/absence matrix was derived from this count matrix, recording a gene family as present in a given genome if at least one gene was assigned to that orthogroup.

The analysis was restricted to 11 target fungal lineages: Blastocladiomycota (n = 5), Chytridiomycota (n = 37), Entomophthoromycota (n = 5), Glomeromycota (n = 7), Kickxellomycota (n = 13), Microsporidia (n = 22), Mortierellomycota (n = 26), Mucoromycota (n = 63), Zoopagomycota (n = 4), Cryptomycota (n = 2) and Basidiobolomycota (n = 2), comprising a total of 186 genomes. Species were assigned to lineages based on a two-letter prefix in their taxon identifier (e.g. MU for Mucoromycota). Genomes from Ascomycota, Basidiomycota, Holozoa, Aphelida, Nucleariida and *Olpidium* were excluded from the focal analysis.

## Pangenome accumulation curves

For each lineage and for the combined dataset (all 186 genomes pooled), we computed pangenome accumulation curves following the random subsampling approach described in Steindorff et al. (2024). At each genome count *n* (ranging from 1 to the total number of genomes *N* in the group), we performed 100 random subsamples without replacement and calculated four statistics per subsample:

- **Pangenome**: number of orthogroups present in at least one of the *n* sampled genomes.
- **Core genome**: number of orthogroups present in all *n* sampled genomes.
- **Shell genome**: number of orthogroups present in ≥ 2 but < *n* sampled genomes (undefined for *n* = 1 and *n* = 2, where it equals zero by definition).
- **Unique (accessory) genome**: number of orthogroups present in exactly one of the *n* sampled genomes.

The mean and standard deviation of each statistic across the 100 subsamples were recorded at each value of *n*. A fixed random seed (42) was used to ensure reproducibility. All computations were performed in Python 3.12 using NumPy 2.4 and pandas 3.0.

## Classification at full sampling depth

At the maximum sampling depth (*n* = *N*), gene families were classified into three mutually exclusive compartments based on their distribution across all genomes in the lineage:

- **Core**: present in all *N* genomes of the lineage.
- **Shell**: present in ≥ 2 but < *N* genomes.
- **Unique**: present in exactly one genome.

These three categories sum to the total pangenome size for each lineage.

## Assessment of pangenome openness

The openness of the pangenome was assessed qualitatively from the shape of the accumulation curve. A **closed pangenome** is characterised by a curve that reaches an asymptote before the maximum number of genomes is sampled, indicating that adding further genomes contributes few new gene families. An **open pangenome** shows no clear plateau, with each additional genome still contributing a substantial number of novel gene families. Heaps' law parameters were not formally fitted given the uneven sampling depths across lineages; however, the rate of new gene family discovery per additional genome was used as a qualitative indicator.

## Visualisation

Pangenome and core genome accumulation curves were plotted as line graphs with shaded bands representing ± 1 standard deviation around the mean across 100 subsamples, following the style of Steindorff et al. (2024, Fig. 2). Separate figures were generated for each of the 11 lineages as well as for the combined dataset. Each figure contains two panels: panel **a** shows the pangenome, shell and unique accumulation curves; panel **b** shows the core genome decay curve. All figures were produced in Python 3.12 using matplotlib 3.10 and exported as vector PDF files (TrueType font embedding, `pdf.fonttype = 42`) to allow editing in Adobe Illustrator, as well as PNG raster previews at 300 dpi.

All analysis scripts are available at https://github.com/cLax8M/EDF-TE under `PanGenome/scripts/`.

## References

- Emms DM, Kelly S. OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biology* 2019; **20**: 238. https://doi.org/10.1186/s13059-019-1832-y
- Steindorff AS, *et al.* Comparative genomic analysis of thermophilic fungi reveals convergent evolutionary adaptations and gene family expansions. *[journal]* 2024.
