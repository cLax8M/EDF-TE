# Nucleosome analysis (Section 3 – epigenetics)

Complete, reproducible pipeline for the EDF-TE nucleosome analyses described
in Methods §11 of the manuscript.  These scripts mirror the analyses that
were carried out during the study (see `Section2_TElandscape/scripts/` in
the project tree) but are packaged here as self-contained, general-purpose
tools that only need the standard scientific Python stack
(numpy, scipy, pandas, matplotlib).

## Workflow at a glance

| Script | Analysis | Key inputs | Key outputs |
|---|---|---|---|
| `01_nuclpred_analysis.py` | Genome-wide NuclPred occupancy → 4 chromatin states + TE-vs-background bias | NuclPred WIG, FAI, TE BED, (GFF) | `*_window_metrics.tsv.gz`, `*_chromatin_states.tsv`, `*_te_bias.tsv` |
| `02_danpos_nucleosome_metrics.py` | Per-nucleosome DANPOS metrics by context (TE/gene/intergenic) | `positions.xls`, TE BED, GFF | `*_nucleosomes.tsv.gz`, `*_summary.tsv` |
| `03_nucleosome_spacing.py` | NRL / spacing from autocorrelation, TE-rich vs TE-poor | WIG, FAI, TE BED | `*_window_spacing.tsv`, `*_spacing_summary.tsv` |
| `04_te_vs_nonte.py` | TE vs gene vs background occupancy (multi-species) | WIG(s), FAI, TE BED, GFF | `te_vs_nonte_summary.tsv` |
| `05_mac_te_profiles.py` | Meta-TE (MAC) occupancy profiles by class | WIG, FAI, TE BED | `*_meta_te_profiles.tsv`, `*_meta_te_summary.tsv` |
| `06_te_fuzziness.py` | Nucleosome fuzziness at TE vs non-TE, and vs H3K9me3 | `positions.xls`, TE BED, GFF, (H3K9 BED) | `*_fuzziness_nucleosomes.tsv.gz`, `*_fuzziness_summary.tsv` |

`helpers/nucleosome_common.py` contains the shared WIG/DANPOS/TE loading
code, window metrics (mean occupancy, CV, autocorrelation NRL/phasing) and
the statistics helpers (BH–FDR, Cliff's delta).

## Statistics included

* Two-sided Mann–Whitney U (TE vs gene / TE vs background) with
  Benjamini–Hochberg FDR correction.
* Kruskal–Wallis across contexts, followed by pairwise Mann–Whitney U
  (Script 06).
* Wilcoxon signed-rank for flank-vs-body meta-TE profiles (Script 05) and
  for cross-species TE+/TE− deltas (Script 04).
* Cliff's delta effect size for every pairwise comparison.
* Spearman correlation and BH-corrected tests in the sRNA–H3K9me3 analyses
  (see `../06_srna_silencing_pipeline.py` and `../07_srna_h3k9me3_statistics.py`).

## Example invocations

```bash
# 01 — one species, NuclPred prediction
python 01_nuclpred_analysis.py \
    --wig    /path/nucleosomes/nuclpred/MU_Rhimi59_nuclpred.wig \
    --fai    /path/genomes/MU_Rhimi59.fasta.fai \
    --te_bed /path/annotated/MU_Rhimi59_EarlGrey/MU_Rhimi59_summaryFiles/MU_Rhimi59.filteredRepeats.bed \
    --gff    /path/genes_annot_gff/MU_Rhimi59.gff \
    --species MU_Rhimi59 --out_dir out/01

# 02 — DANPOS per-nucleosome metrics (N. crassa MNase in vivo)
python 02_danpos_nucleosome_metrics.py \
    --positions /path/AS_Ncrass2_dedup.smooth.positions.xls \
    --te_bed    /path/AS_Ncrass2.filteredRepeats.bed \
    --gff       /path/AS_Ncrass2.gff \
    --label AS_Ncrass2_MNase --out_dir out/02

# 03 — NRL / spacing
python 03_nucleosome_spacing.py \
    --wig /path/Rhizopus_microsporus.wig --fai /path/MU_Rhimi59.fasta.fai \
    --te_bed /path/MU_Rhimi59.filteredRepeats.bed \
    --species MU_Rhimi59 --out_dir out/03

# 04 — TE vs non-TE (multiple species)
python 04_te_vs_nonte.py \
    --wig /path/MU_Rhimi59.wig  --fai /path/MU_Rhimi59.fai  --te_bed /path/MU_Rhimi59.bed  --gff /path/MU_Rhimi59.gff  --species MU_Rhimi59 \
    --wig /path/AS_Ncrass2.wig   --fai /path/AS_Ncrass2.fai   --te_bed /path/AS_Ncrass2.bed   --gff /path/AS_Ncrass2.gff   --species AS_Ncrass2 \
    --out_dir out/04

# 05 — meta-TE profiles
python 05_mac_te_profiles.py \
    --wig /path/MU_Rhimi59.wig --fai /path/MU_Rhimi59.fai \
    --te_bed /path/MU_Rhimi59.filteredRepeats.bed \
    --species MU_Rhimi59 --window 2000 --out_dir out/05

# 06 — fuzziness at TE (+ H3K9me3)
python 06_te_fuzziness.py \
    --positions /path/AS_Ncrass2_dedup.smooth.positions.xls \
    --te_bed /path/AS_Ncrass2.filteredRepeats.bed \
    --gff /path/AS_Ncrass2.gff \
    --h3k9_bed /path/H3K9_WT_MACS2_broad_peaks.bed \
    --label AS_Ncrass2 --out_dir out/06
```

## SLURM array driver

`run_nucleosome_analysis.sh` runs 01/03/05 per species in an array job using
the species list, matching the original `run_nucl_te_interplay.sh` workflow.
Edit the path variables at the top and submit with:

```bash
sbatch --array=1-217 run_nucleosome_analysis.sh
```

## Notes

* Chromosome names in the WIG may carry a `SPECIES.` prefix; this is
  stripped automatically to match the FAI / BED naming.
* Scripts 02 and 06 need DANPOS `positions.xls` (MNase/MAINE in-vivo
  datasets); script 05 also accepts in-vivo DANPOS smooth WIGs.
* ATAC-seq signal is accessibility (opposite direction to occupancy) — flag
  `--normalize minmax` and interpret accordingly when using ATAC tracks.
