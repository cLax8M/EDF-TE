# TreeWAS Analysis — Methods

## Overview

Phylogeny-aware genome-wide association study (pGWAS) to identify gene families and protein
domains associated with genome architecture variables across 183 fungal species. We follow
the methodology of Badet & Croll (2025, *Nature Communications*) using the `treeWAS` R
package (Collins & Didelot 2018, *PLOS Computational Biology*).

---

## Species Selection

Starting from 238 species with complete genomic annotation, we excluded:

| Exclusion criterion | Clades | N species removed |
|---|---|---|
| Outgroups (non-fungi) | HO (Holozoa) | 20 |
| Basal lineages (phylogenetically distant, inflating signal) | AS (Ascomycota), BA (Basidiomycota) | 15 + 14 = 29 |
| Outgroup | OL (*Olpidium*) | 1 |
| Assembly-fragmented (synteny = 0, unmeasurable) | BA_Rhisol1, CH_Batsa1, GL_Rhicl1, MU_Bifad1, MU_Mucam1, MU_Muccirc1 | 6 |

**Final dataset: 183 species**, spanning clades:
AP (2), BS (2), BL (5), CH (35), CR (2), EN (5), GL (6), KX (13), MI (22), MO (26), MU (60), NU (2), ZO (4).

---

## Phenotype Variables (Y)

Eight continuous traits capturing genome architecture:

| Variable | Transformation | Source |
|---|---|---|
| `log10_Genome_Size` | log₁₀(Genome_Size_bp) | GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY.txt |
| `Pct_Intergenic` | none | GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY.txt |
| `Avg_Introns_Multiexon` | none | GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY.txt |
| `Avg_Exons_Per_Gene` | none | GenomicFeatures_V3/ALL_SPECIES_ALL_SUMMARY.txt |
| `Pct_paralog_genes` | none | DuplicationCG/orthogroup_stats.tsv |
| `Pct_OGs_multicopy` | none | DuplicationCG/orthogroup_stats.tsv |
| `Percent_unknown` | none | unknown_function/unknown_function_summary.tsv |
| `mean_r_intra` | none | Synteny/perspecies_synteny.tsv |

`log10_Genome_Size` is included as a control variable, since many genome features correlate
with genome size. Badet & Croll (2025) log₁₀-transformed all continuous Y variables; we
apply this transformation only to genome size since the other variables are already on bounded
or near-normal scales.

---

## Genotype Matrices (X)

### X1: Orthogroup presence/absence

Source: OrthoFinder `Orthogroups.GeneCount.tsv` (148,965 OGs × 238 species).

1. Restrict to 183 included species.
2. Binarise: ≥1 gene copy → 1, else 0.
3. Filter: retain OGs present in ≥5% of species (≥10/183).

Estimated retained OGs: ~15,000–25,000.

### X2: InterPro domain presence/absence

Source: `data/functional_interpro/*_IPR.tab` (238 species, InterProScan output).

1. For each species, extract unique `iprId` values (one row per protein × domain hit).
2. Build species × IPR binary matrix.
3. Restrict to 183 included species.
4. Filter: retain IPR terms present in ≥5% of species (≥10/183).

Estimated retained IPR terms: ~5,000–8,000.

---

## Phylogenetic Tree

**Provisional tree** (used here): OrthoFinder STAG/ASTRAL consensus tree
(`Orthofinder_run/Species_Tree/SpeciesTree_rooted.txt`). Tip labels carry `.aa` suffixes
(e.g. `MU_Mucmuc1.aa`) stripped during preprocessing. Exception: `CH_Batde5` has no `.aa`
suffix — handled explicitly.

**Final tree**: IQ-TREE 2 maximum-likelihood tree (`iqtree_final`, SLURM job 2537000) will
replace the provisional tree once complete. Re-run `00_prep_data.py` with `--tree` pointing
to the final tree.

The provisional tree is pruned to the 183 included species using `ete3`.

---

## treeWAS Analysis

We run **16 analyses** (8 phenotypes × 2 X matrices):

```R
treeWAS(
  snps      = X_matrix,          # species × features binary matrix
  phen      = Y_vector,          # named numeric vector
  tree      = pruned_tree,       # ape phylo object
  n.snps.sim = 10000,            # null simulation size
  test      = c("terminal", "simultaneous", "subsequent")
)
```

### The three tests

| Test | Interpretation |
|---|---|
| **Terminal** | Association at tips; detects recent correlations |
| **Simultaneous** | Changes in X and Y occur on same branch; detects co-evolution |
| **Subsequent** | X change precedes or follows Y change; detects lagged evolution |

A feature is called a **robust hit** if significant in ≥2/3 tests — following Badet & Croll (2025).

---

## TE-Related Pfam Flags

Post-hoc annotation flags hits derived from transposable elements:

| Pfam ID | Description |
|---|---|
| PF00078 | Reverse transcriptase (LTR retrotransposons) |
| PF03184 | DDE transposase (class II TEs) |
| PF13456 | Reverse transcriptase-like |
| PF14214 | Helitron helicase |
| PF13683 | Gypsy integrase |
| PF01497 | Mariner transposase |
| PF13952 | AAA+ATPase (TE mobilization) |

---

## Output Structure

```
TreeWAS/results/
├── phenotypes.tsv              # 183 × 8 phenotype table
├── X_OG.tsv                   # 183 × N_OG binary matrix
├── X_IPR.tsv                  # 183 × N_IPR binary matrix
├── tree_pruned.nwk             # pruned Newick tree
└── {phenotype}_{matrix}/
    ├── treeWAS_results.rds     # full treeWAS output object
    ├── terminal_hits.tsv
    ├── simultaneous_hits.tsv
    └── subsequent_hits.tsv
```

Final annotated table: `TreeWAS/results/treewas_annotated_hits.tsv`

---

## References

- Collins C, Didelot X (2018). A phylogenetic method to perform genome-wide association
  studies in microbes that accounts for population structure and recombination. *PLOS
  Computational Biology*, 14(2): e1005958.
- Badet T, Croll D (2025). Genome architecture differences across the fungal tree of life are
  linked to transposable element families. *Nature Communications*, 16: XXXX.
