# Synteny Analysis — Plan & Ideas

**Project**: EDF-TE (Transposable Element evolution across the Eukaryote-Dikarya-Fungi tree)  
**Analysis**: Within-clade gene microsynteny (collinear ortholog pair frequency)  
**Last updated**: 2026-02-20

---

## 1. Motivation

Synteny conservation is a proxy for genome architectural stability. In this project we want to know:
1. **How stable is gene order within each of the 17 fungal/eukaryotic clades?**  
2. **Is there a correlation between TE load and synteny disruption?** High TE content → frequent ectopic recombination and chromosome rearrangement → lower synteny ratios. This would be a key result linking the TE analysis to genome evolution.
3. **Which lineages have experienced the most genomic rearrangements?** Candidates: Microsporidia (heavily reduced genomes), Mucoromycota (high TE load in some species), vs. Ascomycota (generally stable genomes).
4. **How does synteny scale with phylogenetic distance?** Within tightly clustered clades (e.g. Lucilia spp. in HO) vs. broader ones (e.g. full HO with Metazoa + Amoebozoa).

The approach closely follows **Sebé-Pedrós et al. 2013** (Figure 3C), which applied the same method to 10 unicellular holozoans. We extend it to all 17 clades of the EDF dataset.

---

## 2. Scope Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Pairwise scope | **Within-clade only** | Full all-vs-all (28,920 pairs) is computationally prohibitive; inter-clade synteny is near-zero for most combinations anyway |
| Ortholog source | **OrthoFinder pairwise 1:1 tables** (`Orthologues/`) | Global single-copy OGs set is empty for 241 divergent species; pairwise 1:1 tables are the correct and available input |
| Synteny definition | **Adjacent ortholog pairs on same scaffold** | Strand-agnostic (orientation-independent), per Putnam et al. 2007 |
| Normalization | **Per-clade maximum** | Clades differ vastly in divergence scale; a global normalizer would compress most clades to near-zero |
| Null distribution | **100× bedtools shuffle** per species, `-chrom` flag (stay on same scaffold), fixed seed per replicate | Exact replication of paper's method |
| Minimum scaffold genes | **≥ 3 genes** | Single- or two-gene scaffolds cannot contribute adjacency information |

---

## 3. Clade Overview

| Clade | Code | N spp | N pairs | Notes |
|---|---|---|---|---|
| Aphelida | AP | 2 | 1 | Single pair only |
| Ascomycota | AS | ~13 | ~78 | Well-studied; expected high synteny within Saccharomycotina |
| Basidiomycota | BA | ~10 | ~45 | |
| Blastocladiomycota | BL | 4 | 6 | |
| Basidiobolomycota | BS | 2 | 1 | Single pair only |
| Chytridiomycota | CH | ~10 | ~45 | |
| Cryptomycota/Rozellida | CR | 2 | 1 | Single pair only |
| Entomophthorales | EN | 5 | 10 | |
| Glomeromycota | GL | 6 | 15 | |
| **Holozoa** | **HO** | **20** | **190** | Most complex; very heterogeneous phylogenetic distances |
| Kickxellomycotina | KX | ~6 | ~15 | |
| Microsporidia | MI | ~6 | ~15 | Expected lowest synteny (reduced genomes) |
| Mortierellomycota | MO | ~6 | ~15 | |
| Mucoromycota | MU | ~7 | ~21 | Includes 3 in-house PacBio assemblies |
| Nucleariida | NU | 2 | 1 | Single pair only |
| Olpidium | OL | 1 | 0 | **Skip** — single species |
| Zoopagomycota | ZO | 4 | 6 | |

> **Species counts to be verified**: run `python 02_extract_orthologs.py --all --dry_run` to check exact pair counts and OrthoFinder file availability.

---

## 4. The HO Group

### 4.1 Composition

The HO group (Holozoa outgroups) contains 20 species from very different phylogenetic backgrounds:

| Subgroup | Species |
|---|---|
| **Insecta — Diptera** | *Drosophila melanogaster*, *Musca domestica*, *Calliphora vicina*, *Lucilia cuprina*, *Lucilia sericata* |
| **Insecta — Lepidoptera** | *Plutella xylostella* (focal species) |
| **Insecta — Hemiptera** | *Magicicada septendecula* |
| **Other Metazoa** | *Hydra vulgaris* (Cnidaria), *Trichoplax adhaerens* (Placozoa), *Ciona intestinalis* (Tunicata), *Mnemiopsis leidyi* (Ctenophora) |
| **Unicellular Holozoa** | *Monosiga brevicollis*, *Salpingoeca rosetta* (Choanoflagellata); *Capsaspora owczarzaki* (Filasterea); *Sphaeroforma arctica* (Ichthyosporea) |
| **Other Eukaryotes** | *Naegleria gruberi* (Excavata), *Acanthamoeba castellanii*, *Dictyostelium discoideum* (Amoebozoa), *Thecamonas trahens* (Apusozoa) |

### 4.2 Synteny expectations

- **Lucilia cuprina × Lucilia sericata**: highest r_ij expected (congeneric flies); use as **positive control** for step sanity check
- **Diptera block**: Drosophila, Musca, Calliphora, Lucilia spp. should show moderate-to-high synteny
- **Trans-metazoan comparisons** (e.g. Hydra × Drosophila): low but detectable synteny (macrosynteny blocks known from literature)
- **Metazoa × Unicellular Holozoa**: near-zero (no conserved microsynteny expected)
- **Holozoa × Amoebozoa/Excavata**: zero (completely unrelated genome architectures)

The HO heatmap will therefore be **bimodal**: a high-synteny insect block + a low/zero block for non-insect comparisons. This is biologically expected and should be noted in the paper.

### 4.3 Article definition (draft)

*"The HO group comprises 20 non-Holomycota eukaryotes included as phylogenetically distant outgroups. The group intentionally over-represents insects (7 species: 5 Diptera, 1 Lepidoptera, 1 Hemiptera) to provide a detailed comparative context for the focal species* Plutella xylostella. *The remaining 13 species span the full diversity of Holozoa (4 unicellular holozoans from Choanoflagellata, Filasterea and Ichthyosporea) and other eukaryotic supergroups (Excavata, Amoebozoa, Apusozoa). This design maximizes phylogenetic coverage while anchoring comparisons to the lepidopteran outgroup."*

---

## 5. Pipeline Overview

```
GFF files                    OrthoFinder Orthologues/
    |                                |
    v                                v
[01] gff_to_bed.py          [02] extract_orthologs.py
    (→ .bed, .genome)              (→ 1:1 pairs per clade)
         \                          /
          \________________________/
                     |
                     v
            [03] collinear_pairs.py
              (→ c_ij, N_ij per pair)
                     |
                     v
            [04] shuffle_null.sh
              (100× bedtools shuffle → s_ij)
                     |
                     v
            [05] synteny_ratio.py
              (→ r_ij normalized within clade)
                     |
                     v
            [06] plot_synteny.py
              (→ heatmap PDF + PNG per clade)
                     |
                     v (optional)
            [07] circos_prep.py
              (→ Circos input for selected pairs)
```

---

## 6. Expected Figures

| Figure | Description |
|---|---|
| **Fig X** | 17 heatmaps (one per clade), arranged as a grid or supplementary multipanel |
| **Fig X+1** | Scatter plot: synteny ratio vs. TE content per species (linking SyntenyTE analysis) |
| **Fig X+2** | Circos plot for the most synteny-conserved pair in HO (*Lucilia* spp.) and in a key fungal clade (TBD) |
| **Supp** | Assembly quality metrics (N50, scaffold count) vs. synteny ratio per species |

---

## 7. TODOs

- [ ] **Run dry_run for step 02** to verify OrthoFinder file paths and species counts per clade
- [ ] **Test pipeline end-to-end with HO_Naegrub + HO_Sphaearti** (smallest HO GFF files) before full run
- [ ] **Extract assembly quality stats** from `Samples.xlsx` (N50, scaffold count) to flag fragmented assemblies in heatmaps
- [ ] **Verify HO species list** against actual GFF files in `genes_annot_gff/` (some species may have `_c` suffix in filenames)
- [ ] **Investigate MU clade size** (up to 64 spp = 2016 pairs): consider submitting step 04 as a SLURM array job
- [ ] **Correlate r_ij with TE load** after both analyses are complete
- [ ] **Add HO sub-clade labels** to heatmap (insect vs. metazoa vs. unicellular) using seaborn row/col colors
- [ ] **Validate interpretation** of near-zero r_ij pairs (biology vs. assembly fragmentation)

---

## 8. Software & Versions

| Tool | Version | Purpose |
|---|---|---|
| Python | 3.12 | Pipeline scripts |
| pandas | ≥ 2.0 | Data manipulation |
| numpy | ≥ 1.25 | Numerical ops |
| seaborn | ≥ 0.13 | Heatmaps |
| matplotlib | ≥ 3.8 | Plotting |
| bedtools | ≥ 2.24 | Shuffle step |
| GNU parallel | latest | Parallelization of shuffle step |
| Circos | 0.67 | Optional syntenic circle plots |
| OrthoFinder | 2.x | Ortholog source (pre-computed) |

---

## 9. References

- Sebé-Pedrós A, Irimia M, del Campo J, et al. (2013). *Regulated aggregative multicellularity in a close unicellular relative of metazoa.* eLife 2:e01287. [Paper containing Fig 3C and the method used here]
- Putnam NH, et al. (2007). *Sea anemone genome reveals ancestral eumetazoan gene repertoire and genomic organization.* Science 317:86–94.
- Quinlan AR, Hall IM. (2010). *BEDTools: a flexible suite of utilities for comparing genomic features.* Bioinformatics 26:841–842.
- Krzywinski M, et al. (2009). *Circos: an information aesthetic for comparative genomics.* Genome Research 19:1639–1645.
- Emms DM, Kelly S. (2019). *OrthoFinder: phylogenetic orthology inference for comparative genomics.* Genome Biology 20:238.
