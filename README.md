# EDF-TE
Code for the study of TE evolution in the Holomycota (EDF project).

## Repository structure

| Directory | Description |
|-----------|-------------|
| `SpeciesTree/` | Phylogenomic marker extraction, alignment, trimming, gene trees, species tree (IQ-TREE) |
| `PanGenome/` | Pangenome computation and visualization across 238 species |
| `Synteny/` | Collinearity and synteny analysis |
| `TreeWAS/` | Phylogenetic GWAS of genomic features |
| `DuplicationCG/` | Gene duplication statistics from OrthoFinder |
| `GenomicFeatures/` | Genome structure metrics (intron density, gene length, etc.) |
| `TE_analysis/` | EarlGrey library curation and annotation pipeline |
| `Cooption/` | TE co-option analysis: domain annotation, fusion detection, gene classification |
| `Section2_TElandscape/` | TE landscape analyses: diversity, Kimura divergence, PGLS, ancestral reconstruction, nucleosome insertion bias |

## Data

Large data files (genomes, BAMs, TSVs, WIGs) are **not tracked** in this repository.
See `/legserv/Temp/Carlos/EDF-TE/data/` and `/legserv/Temp/Carlos/EDF-TE/results/` for all data.

## Analysis modules

### TE Co-option (`Cooption/`)
Detects genes with TE-derived domains that have been co-opted for host functions.
Pipeline: InterProScan normalization (00) → domain annotation (01) → protein classification (02) →
TE/host domain splitting (03) → host domain extraction (04) → gene classification (05) → true co-option filter (06).

**Validation control:** SETMAR (human gene with Mariner transposase + SET domain) must appear in
`HO_Homosap.true_TE_coopted_genes.tsv` after running Phase 2.

Run Phase 2 for all species:
```bash
sbatch Cooption/scripts/run_cooption_phase2.sh
# Dry-run to preview:
bash Cooption/scripts/run_cooption_phase2.sh --dry-run
```

### TE Landscape (`Section2_TElandscape/`)

**Shannon diversity (B):** TE superfamily diversity per species (H', Pielou's J) + Pagel's λ
```bash
python Section2_TElandscape/scripts/B01_shannon_diversity.py <annotated_dir> <output.tsv>
Rscript Section2_TElandscape/scripts/B02_pagel_lambda.R <shannon.tsv> <tree.nwk> <output.tsv>
```

**Kimura divergence (C):** Recency index, burst score, weighted mean divergence per species/superfamily
```bash
python Section2_TElandscape/scripts/C01_kimura_metrics.py <annotated_dir> <species.tsv> <superfamily.tsv>
```

**TE ~ Genome features (D):** PGLS regression of TE% against intron density, genome size, gene count, etc.
```bash
python Section2_TElandscape/scripts/D01_merge_te_genomics.py <annotated_dir> <genomicfeatures.tsv> <merged.tsv>
Rscript Section2_TElandscape/scripts/D02_pgls_regression.R <merged.tsv> <tree.nwk> <output_dir/>
```

**Ancestral reconstruction (A):** Brownian motion ancestral state reconstruction of TE% on the species tree
```bash
Rscript Section2_TElandscape/scripts/A01_ancestral_reconstruction.R <annotated_dir> <tree.nwk> <output_dir/>
```

**Nucleosome insertion bias (F):** NuclPred predictions on 238 genomes + TE insertion occupancy analysis
```bash
# Setup: clone NuclPred and install conda env
git clone https://github.com/sjmondo/NuclPred_v1 TE_analysis/NuclPred_v1
conda env create -f TE_analysis/NuclPred_v1/nucleosome_prediction_env.yml
# Run NuclPred for all 238 genomes (SLURM array):
sbatch --array=1-238 Section2_TElandscape/scripts/run_nuclpred_array.sh
# Run Nuc-seq pipeline for 5 experimental species (SLURM array):
sbatch Section2_TElandscape/scripts/run_nucseq_array.sh
# Compute TE insertion bias per species:
python Section2_TElandscape/scripts/F3_te_nucleosome_bias.py <wig> <bed> <sp> <output.tsv>
# Validate NuclPred vs experimental data:
python Section2_TElandscape/scripts/F4_validate_nuclpred.py <predicted.wig> <observed.wig> <sp> <output_dir/>
```

## Dependencies

All R scripts require the `treewasEDFTE` conda environment (`TreeWAS/treewasEDFTE.yml`).
Nucleosome scripts require additional conda environments:
- `nucseq_pipeline`: sra-tools, bowtie, samtools, picard, danpos2, bbtools
- `nuclpred`: see `TE_analysis/NuclPred_v1/nucleosome_prediction_env.yml` + DNAshapeR (Bioconductor)

