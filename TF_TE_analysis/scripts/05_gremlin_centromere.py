#!/usr/bin/env python3
"""
05_gremlin_centromere.py
─────────────────────────────────────────────────────────────────────────────
Grem-LINE1 putative centromere analysis in Mucoromycota.

Reference strain (Navarro-Mendoza 2019, Curr Biol):
  MU402  = Mucor circinelloides f. lusitanicus CBS277.49 — JGI assembly.
           Genome: /legserv/Temp/Carlos/EDF-TE/data/Extra/mu402GenomeJGI.fa
           GremLINE1_table_S3.bed = known centromere positions in MU402 coords.
  MU_Mucci3  = CBSv3 assembly of the SAME organism (newer/better assembly).
               Used as calibration species (9 expected chromosomes → 9 peaks).
               Gremlin families: rnd-1_family-198/180/206/106 (in CBSv3).
               NOT the same as MU_Muccirc1.
  MU_Muccirc1 = a different Mucor circinelloides strain (fragmented KE_ scaffolds,
                ~150 gremlins spread across hundreds of tiny contigs → 0 centromeres
                predicted, expected for a fragmented assembly).

Gremlin identification strategy (paper-validated):
  Primary  : tBLASTn ORF2 proteins → whole genome, pident ≥ 30%, qcov ≥ 10%
             (thresholds from Navarro-Mendoza 2019: ≥30% id, ≥50% coverage).
             Hit coordinates loaded directly via load_genome_hits_direct().
             Used for Analyses 1 (density peaks) + flanking/inverted-repeat check.
  Fallback : EarlGrey de-novo family names (kimura80 available, used for
             Analyses 2/3/4 that require divergence values or small-scaffold
             coverage).

Analyses
  1.  Centromere prediction — density-based (1 peak per scaffold); predicts
      chromosome / centromere count across 63 MU species.
      Flanking check: each peak is tested for inverted-repeat configuration
      (gremlins on BOTH sides, oriented AWAY from core — Navarro-Mendoza 2019).
      Calibration: MU_Mucci3 (CBS277.49 CBSv3) → 9 putative centromeres = 9 chr.

  2.  H3K9me3 peak-centric validation in Rhizopus microsporus (MU_Rhimi59).
      H3K9me3-first strategy: start from MACS2 broad peaks (heterochromatin
      domains) and ask which ones contain Grem-LINE1 elements.
      Statistics: binomial enrichment test, Mann-Whitney (peak score gremlin+
      vs gremlin-), Spearman ρ (peak score vs n_gremlin per peak).
      Visualization: scatter (score × n_gremlin), violin (score gremlin±),
      IGV mini-tracks for top gremlin-containing peaks.

  3.  Cross-species conservation in Mucorales.
      For Mucor & Rhizopus spp.: centromere count, orientation bias, K80 at
      centromere peaks vs genome-wide (dynamic vs. stable?).

  4.  K80 divergence comparison — Grem-LINE1 vs other TE classes.
      Violin plots across all 63 MU species.
      Secondary age axis with Szantho et al. 2025 calibration lines.

  5.  Integrated TF evidence on LINE/I (general analysis, not gremlin-specific).
      Re-uses: 01_peak_overlap matrix, 02_bigwig_signal mean, 03_fimo_enrichment.
      Shows top TFs from all three evidence sources for LINE/I elements.

Calibrations: Szantho et al. 2025, Nat. Ecol. Evol. 9:1989

Output
  <OUTDIR>/gremlin_clusters/       per-species BED + summary TSV
  <OUTDIR>/h3k9me3_igv/            per-scaffold IGV data
  <OUTDIR>/k80_comparison/         per-class K80 npy + TSV
  <OUTDIR>/gremlin_centromere_analysis.pdf

Usage
  python 05_gremlin_centromere.py [--skip-clusters] [--skip-k80]
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from scipy.stats import binomtest, mannwhitneyu, spearmanr
import pyBigWig

# ── CONFIG ────────────────────────────────────────────────────────────────────
ANNOT_ROOT = Path("/data/carlos/EDF-TE/TE_analysis/annotated")
OUTDIR     = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/05_gremlin_centromere")

PEAK_MATRIX_TSV = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/"
                       "01_peak_overlap/peak_TE_overlap_matrix.tsv")
BW_SIGNAL_TSV   = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/"
                       "02_bigwig_signal/bw_TE_signal_class_mean.tsv")
FIMO_TSV        = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/"
                       "03_fimo_enrichment/fimo_enrichment.tsv")

H3K9_BW  = Path("/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/"
                "Epigenetics/RmH3K9me3/bbC_WT_H3K9_CPM10.bigwig")
H3K9_BED = Path("/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/"
                "Epigenetics/RmH3K9me3/H3K9_WT_MACS2_broad_peaks.bed")
# Genome index for MU_Rhimi59 (used to compute genome coverage fraction)
RHIMI59_FAI = Path("/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/"
                   "MU_Rhimi59.fasta.fai")

GREMLIN_FAMILIES_TSV    = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/tblastn_results/"
                               "gremlin_families_tblastn_denovo.tsv")
GREMLIN_FAMILIES_GENOME = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/genome_tblastn_results/"
                               "gremlin_families_genome.tsv")
GREMLIN_FAMILIES_BLASTN = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/blast_results/"
                               "gremlin_families_per_species.tsv")

# Direct genome hit BED files from 05d tBLASTn (ORF2 proteins vs whole-genome)
# These bypass EarlGrey classification — positions are identified directly from
# homology to Grem-LINE1 ORF2 proteins (Navarro-Mendoza 2019 sequences).
GENOME_HITS_DIR  = Path("/legserv/Temp/Carlos/EDF-TE/data/gremlin/genome_tblastn_results")
PIDENT_STRICT    = 30.0   # paper threshold: ≥30% identity (Navarro-Mendoza 2019)
GREMLIN_MERGE_GAP = 500   # bp: merge overlapping HSPs → single element position

# Centromere prediction parameters
DENSITY_WINDOW    = 100_000   # 100 kb sliding window for density peak
DENSITY_STEP      = 10_000    # 10 kb step
MIN_GREM_SCAFFOLD = 5         # min gremlin elements per scaffold -> putative centromere
                              # Calibration: MU_Mucci3 (CBSv3 = M. circ. f. lusitanicus
                              # CBS277.49, i.e. paper strain) → 9 scaffolds with >=5 = 9 chr.
                              # Note: MU_Muccirc1 = DIFFERENT strain with fragmented assembly
                              # (~150 gremlins on KE_ scaffolds, max 4/scaffold → 0 peaks).

# Flanking / inverted-repeat filter parameters
MIN_FLANK_GREM   = 2    # min merged gremlin elements on EACH side of peak center
MIN_ORIENT_FRAC  = 0.60  # fraction threshold to call strand bias (inverted repeat)

# IGV parameters
IGV_FLANK    = 100_000   # +/-100 kb around centromere peak
IGV_BINS     = 300
MAX_IGV_PANELS = 6

# Szantho et al. 2025 node ages (midpoint Ma; 95% HPD low-high)
CLADE_AGES = {
    "LCA Mucor.+Dikary.": (1107, 1080, 1140),
    "Crown\nMucoromycota":  ( 946,  678, 1213),
    "Glomeromycotina":      ( 494,  408,  580),
    "Endogonales":          ( 208,   76,  340),
}

MU_RATES = {
    "3e-9/yr (slow)": 3e-9,
    "1e-8/yr (fast)": 1e-8,
}

MU_SPECIES = [
    "MU_Abscae1","MU_Absrep1","MU_Amyrou1","MU_Apoossi","MU_Bacci1",
    "MU_Baclam1","MU_Benpoi1","MU_Bifad1","MU_Blatri","MU_Calth1",
    "MU_Chabre1","MU_Chlapade","MU_Chocucu1","MU_Cirumb1","MU_Cokrec1",
    "MU_Cunech1","MU_Dicele1","MU_Dicrob1","MU_Ellano2","MU_Endsp1",
    "MU_Fenlin","MU_Gilper1","MU_Gonbut1","MU_Halrad1","MU_Helpul1",
    "MU_Hesve2","MU_Jimfl","MU_Jimlac1","MU_Kircor1","MU_Lenpar1",
    "MU_Lenzyc1","MU_Liccor1","MU_Licra1","MU_Mucam1","MU_Mucci3",
    "MU_Muccirc1","MU_Mucend1","MU_Muclan1","MU_Mucmuc1","MU_Mucpir1",
    "MU_Mucrac1","MU_Mycafr1","MU_Mycind1","MU_Parpar1","MU_Phaart1",
    "MU_Phybl4","MU_Pilano1","MU_Pilumb1","MU_PnitS607","MU_Radspe1",
    "MU_Rhimi59","MU_Rhior3","MU_Rhipu1","MU_Rhisto1","MU_Sakvas1",
    "MU_Spifus1","MU_Spoumb1","MU_Synrac1","MU_Syzme1","MU_Thaele1",
    "MU_Umban1","MU_Umbra1","MU_Zycmex1",
]

MUCORALES_PANEL = [
    "MU_Mucci3","MU_Muccirc1","MU_Mucmuc1","MU_Mucpir1",
    "MU_Muclan1","MU_Mucend1","MU_Mucam1","MU_Mucrac1",
    "MU_Rhimi59","MU_Rhior3","MU_Rhipu1","MU_Rhisto1",
]

CLADE_MAP = {
    "MU_Endsp1": "Endogonales",  "MU_Jimfl": "Endogonales",
    "MU_Jimlac1": "Endogonales",
    "MU_Lenpar1": "Glomeromycotina", "MU_Lenzyc1": "Glomeromycotina",
    "MU_Calth1":  "Glomeromycotina",
}
CLADE_COLORS = {
    "Endogonales":     "#7B68EE",
    "Glomeromycotina": "#20B2AA",
    "Mucorales":       "#E07B39",
}
SKIP_TE = {"Simple_repeat","Low_complexity","Satellite","rRNA","snRNA",
           "Unknown","Unspecified","Other"}


# ── HELPERS ───────────────────────────────────────────────────────────────────

def gff_path(species):
    return ANNOT_ROOT / f"{species}_EarlGrey" / f"{species}_summaryFiles" / \
           f"{species}.filteredRepeats.gff"


def load_gff_as_bed(path):
    pat_name   = re.compile(r"Name=([^;]+)")
    pat_kimura = re.compile(r"kimura80=([0-9.]+)")
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            m_n = pat_name.search(p[8])
            m_k = pat_kimura.search(p[8])
            rows.append((
                p[0], int(p[3]) - 1, int(p[4]), p[2], p[6],
                m_n.group(1) if m_n else p[2],
                float(m_k.group(1)) if m_k else float("nan"),
            ))
    return pd.DataFrame(rows, columns=["chrom","start","end","te_class",
                                       "strand","te_name","kimura80"])


def _merge_genome_hits(df, gap=GREMLIN_MERGE_GAP):
    """
    Merge overlapping / nearby tBLASTn HSPs on the same (chrom, strand) into
    single element positions.  Returns a DataFrame with columns compatible
    with the EarlGrey GFF format used everywhere else in this script.
    """
    rows = []
    for (chrom, strand), sub in df.groupby(["chrom", "strand"]):
        sub = sub.sort_values("start").reset_index(drop=True)
        cur_s = int(sub.loc[0, "start"])
        cur_e = int(sub.loc[0, "end"])
        cur_p = float(sub.loc[0, "pident"])
        for i in range(1, len(sub)):
            s = int(sub.loc[i, "start"])
            e = int(sub.loc[i, "end"])
            p = float(sub.loc[i, "pident"])
            if s <= cur_e + gap:
                cur_e = max(cur_e, e)
                cur_p = max(cur_p, p)
            else:
                rows.append({"chrom": chrom, "start": cur_s, "end": cur_e,
                             "strand": strand, "pident": cur_p,
                             "te_name": "GremLINE1", "te_class": "LINE/I",
                             "kimura80": float("nan")})
                cur_s, cur_e, cur_p = s, e, p
        rows.append({"chrom": chrom, "start": cur_s, "end": cur_e,
                     "strand": strand, "pident": cur_p,
                     "te_name": "GremLINE1", "te_class": "LINE/I",
                     "kimura80": float("nan")})
    if not rows:
        return pd.DataFrame(columns=["chrom","start","end","te_class",
                                     "strand","te_name","kimura80"])
    result = pd.DataFrame(rows)
    return result[["chrom","start","end","te_class","strand","te_name","kimura80"]]


def load_genome_hits_direct(species, pident_min=None):
    """
    Load per-species genome tBLASTn BED (from script 05d) and return
    deduplicated Grem-LINE1 element positions filtered by strict identity.
    Returns None if no BED file exists for this species (fall back to EarlGrey).
    """
    if pident_min is None:
        pident_min = PIDENT_STRICT
    bed = GENOME_HITS_DIR / f"{species}_genome_hits.bed"
    if not bed.exists():
        return None
    df = pd.read_csv(bed, sep="\t", header=None,
                     names=["chrom","start","end","query","pident","strand"])
    df = df[df["pident"] >= pident_min].copy()
    if df.empty:
        return pd.DataFrame(columns=["chrom","start","end","te_class",
                                     "strand","te_name","kimura80"])
    return _merge_genome_hits(df)


def load_gremlin_families():
    families = {}

    def _load(path, label):
        if not path.exists():
            return
        print(f"  Loading gremlin families from {label}: {path.name}")
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return
        for sp, grp in df.groupby("species"):
            families.setdefault(sp, set()).update(grp["family_name"].tolist())

    _load(GREMLIN_FAMILIES_TSV,    "tBLASTn de-novo (05c)")
    _load(GREMLIN_FAMILIES_GENOME, "genome tBLASTn (05d)")
    _load(GREMLIN_FAMILIES_BLASTN, "blastn fallback (05b)")
    if families:
        print(f"  Loaded gremlin families for {len(families)} species")
    else:
        print("  INFO: no gremlin families -- will use LINE/I class as fallback")
    return families


def filter_gremlins(gff_df, species, gfams, use_genome_hits=True):
    # 1. Prefer direct genome-hit positions (pident >= PIDENT_STRICT, paper threshold)
    #    These bypass EarlGrey classification entirely.
    #    Pass use_genome_hits=False for analyses requiring kimura80 or full scaffold
    #    coverage (e.g. Analysis 2 H3K9me3 on small centromeric scaffolds, Analysis 4 K80).
    if use_genome_hits:
        ghits = load_genome_hits_direct(species)
        if ghits is not None and not ghits.empty:
            return ghits
    # 2. EarlGrey-based fallback: use tBLASTn-assigned family names
    fams = gfams.get(species)
    if fams:
        mask = gff_df["te_name"].isin(fams)
        result = gff_df[mask].copy()
        if result.empty:
            fams_low = {f.lower() for f in fams}
            mask = gff_df["te_name"].str.lower().isin(fams_low)
            result = gff_df[mask].copy()
        return result
    return gff_df[gff_df["te_class"].str.startswith("LINE/I")].copy()


# ── ANALYSIS 1: Centromere-peak prediction ────────────────────────────────────

def _flanking_check(gremlins_sub, center):
    """
    Test whether gremlin elements form an inverted-repeat flanking arrangement
    around a putative centromere core, as described by Navarro-Mendoza 2019.

    Criteria:
      has_both_flanks : ≥ MIN_FLANK_GREM elements on EACH side of center
      inverted_repeat : elements on left are predominantly + strand while
                        elements on right are predominantly – strand (or vice
                        versa), consistent with orientation AWAY from the core.

    Returns (has_both_flanks, inverted_repeat, n_left, n_right).
    """
    left  = gremlins_sub[gremlins_sub["start"] <  center]
    right = gremlins_sub[gremlins_sub["start"] >= center]
    n_left, n_right = len(left), len(right)
    has_both = n_left >= MIN_FLANK_GREM and n_right >= MIN_FLANK_GREM
    inverted = False
    if has_both:
        left_plus  = (left["strand"]  == "+").mean()
        right_plus = (right["strand"] == "+").mean()
        # Inverted repeat: one side + dominated, opposite side – dominated
        inverted = (
            (left_plus >= MIN_ORIENT_FRAC and right_plus <= 1.0 - MIN_ORIENT_FRAC) or
            (left_plus <= 1.0 - MIN_ORIENT_FRAC and right_plus >= MIN_ORIENT_FRAC)
        )
    return has_both, inverted, int(n_left), int(n_right)


def _scaffold_gremlin_peak(positions, chrom_len,
                           window=DENSITY_WINDOW, step=DENSITY_STEP):
    """
    Find position of max Grem-LINE1 density in a sliding window.
    Uses O(N log N) two-pointer approach for speed (safe for Jimfl ~ 22k elements).
    """
    if len(positions) == 0:
        return None, 0, None, None
    pos = np.sort(positions)
    best_count = 0
    best_lo_idx = 0
    left = 0
    for right in range(len(pos)):
        # advance left until the window fits
        while pos[right] - pos[left] > window:
            left += 1
        count = right - left + 1
        if count > best_count:
            best_count = count
            best_lo_idx = left
    # center = midpoint of the best window
    best_hi_idx = best_lo_idx + best_count - 1
    best_center = int((pos[best_lo_idx] + pos[best_hi_idx]) // 2)
    best_start  = max(0, best_center - window // 2)
    best_end    = min(chrom_len, best_center + window // 2)
    return best_center, best_count, best_start, best_end


def analyse_centromere_peaks(species_list, gfams):
    """
    Per species: find 1 density-peak per scaffold (putative centromere).
    Scaffolds with >= MIN_GREM_SCAFFOLD gremlin elements qualify.
    n_putative_centromeres = n_qualifying_scaffolds.
    """
    cluster_dir = OUTDIR / "gremlin_clusters"
    cluster_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for sp in species_list:
        gp = gff_path(sp)
        if not gp.exists():
            print(f"  WARN: {gp} not found, skipping {sp}")
            continue
        gff_df   = load_gff_as_bed(gp)
        gremlins = filter_gremlins(gff_df, sp, gfams)
        n_total  = len(gremlins)

        centro_rows = []
        for scaf, sub in gremlins.groupby("chrom"):
            # Qualify on TOTAL scaffold count, not peak-window count.
            # This correctly identifies chromosome-bearing scaffolds in
            # chromosome-level assemblies (e.g. Mucci3 -> 9 scaffolds with >=5).
            if len(sub) < MIN_GREM_SCAFFOLD:
                continue
            chrom_len = int(gff_df[gff_df["chrom"] == scaf]["end"].max()) + 1
            positions = sub["start"].values
            center, peak_count, win_s, win_e = _scaffold_gremlin_peak(positions, chrom_len)
            in_win  = sub[(sub["start"] >= win_s) & (sub["end"] <= win_e)]
            n_plus  = int((in_win["strand"] == "+").sum())
            n_minus = int((in_win["strand"] == "-").sum())
            n_str   = n_plus + n_minus
            orient_p = (binomtest(max(n_plus, n_minus), n_str, 0.5,
                                  alternative="greater").pvalue
                        if n_str >= 2 else 1.0)
            k80_win = in_win["kimura80"].dropna().values
            has_both, inverted, n_left, n_right = _flanking_check(sub, center)
            centro_rows.append({
                "species": sp, "chrom": scaf,
                "center": center, "win_start": win_s, "win_end": win_e,
                "n_gremlin_window": peak_count, "n_gremlin_scaf": len(sub),
                "n_plus": n_plus, "n_minus": n_minus,
                "orient_pval": orient_p, "orient_biased": orient_p < 0.05,
                "k80_mean_window": float(np.mean(k80_win)) if len(k80_win) else np.nan,
                "has_both_flanks": has_both,
                "inverted_repeat":  inverted,
                "n_left_flank":    n_left,
                "n_right_flank":   n_right,
            })

        centro_df = pd.DataFrame(centro_rows)
        n_centro  = len(centro_df)
        k80_all   = gremlins["kimura80"].dropna().values

        if not centro_df.empty:
            out = cluster_dir / f"{sp}_centromes.bed"
            centro_df[["chrom","win_start","win_end","n_gremlin_window"]].to_csv(
                out, sep="\t", index=False, header=False)

        rows.append({
            "species":               sp,
            "n_gremlin_total":       n_total,
            "n_putative_centromeres": n_centro,
            "n_orient_biased":       int(centro_df["orient_biased"].sum()) if not centro_df.empty else 0,
            "n_both_flanks":         int(centro_df["has_both_flanks"].sum()) if not centro_df.empty else 0,
            "n_inverted_repeat":     int(centro_df["inverted_repeat"].sum())  if not centro_df.empty else 0,
            "k80_mean_all":          float(np.mean(k80_all)) if len(k80_all) else np.nan,
            "k80_mean_centromere":   float(centro_df["k80_mean_window"].mean()) if not centro_df.empty else np.nan,
        })
        print(f"  {sp}: {n_total} gremlin -> {n_centro} putative centromeres "
              f"({int(centro_df['has_both_flanks'].sum()) if not centro_df.empty else 0} both-flanks, "
              f"{int(centro_df['inverted_repeat'].sum()) if not centro_df.empty else 0} inverted-repeat)")

    summary = pd.DataFrame(rows)
    summary.to_csv(cluster_dir / "centromere_prediction_summary.tsv", sep="\t", index=False)

    # Combined centromere table
    all_c = []
    for sp in species_list:
        p = cluster_dir / f"{sp}_centromes.bed"
        if p.exists():
            d = pd.read_csv(p, sep="\t", header=None,
                            names=["chrom","win_start","win_end","n_gremlin_window"])
            d.insert(0, "species", sp)
            all_c.append(d)
    if all_c:
        pd.concat(all_c).to_csv(cluster_dir / "all_centromere_peaks.tsv",
                                sep="\t", index=False)
    return summary


# ── ANALYSIS 2: H3K9me3 peak-centric validation ──────────────────────────────

def analyse_h3k9me3_igv(gfams):
    """
    H3K9me3 scaffold-classification analysis in MU_Rhimi59.

    Key discovery from spatial analysis:
    - Large scaffolds (1-3 Mb) = chromosome arm scaffolds.  Gremlins
      accumulate in the scaffold INTERIOR (centromeric position) while H3K9me3
      is preferentially enriched at scaffold ENDS (subtelomeric).
    - Small scaffolds (<200 kb) can be centromeric sequences — centromere DNA
      is so repetitive that it assembles in short, discontinuous contigs.
      These show uniformly HIGH H3K9me3 throughout AND contain gremlin elements.

    Architecture model:
      [subtelo H3K9]──[GREMLINS pericent.]──[centromeric scaffold H3K9+Grem.]
                                         ──[GREMLINS pericent.]──[subtelo H3K9]

    Classification criteria:
      "centromeric scaffold" : size < CENTRO_MAX_KB kb, H3K9me3 mean > CENTRO_H3K9
      "arm scaffold"         : size >= CENTRO_MAX_KB kb, gremlin+ (n >= MIN_GREM_SCAFFOLD)

    Statistics:
      1. Mann-Whitney: H3K9me3 mean in centromeric scaffolds vs arm scaffolds.
      2. Fraction of genome in each scaffold class.
      3. Number of centromeric scaffold candidates.

    Outputs:
      Dict with scaffold classification table, per-class H3K9me3 arrays,
      IGV track data for the best centromeric scaffold, and a full-scaffold
      profile for the best arm scaffold.
    """
    igv_dir = OUTDIR / "h3k9me3_igv"
    igv_dir.mkdir(parents=True, exist_ok=True)

    sp = "MU_Rhimi59"
    gp = gff_path(sp)
    if not (gp.exists() and H3K9_BW.exists()):
        print("  WARN: Rhimi59 GFF or H3K9 bigwig not found -- skipping")
        return {}

    CENTRO_MAX_KB = 500     # scaffolds < 500 kb may be assembled centromere fragments
    CENTRO_H3K9   = 3.0    # CPM threshold for "uniformly H3K9me3-positive" scaffold

    print(f"  Loading {sp} GFF…")
    gff_df   = load_gff_as_bed(gp)
    # Analysis 2 requires EarlGrey elements (small centromeric scaffolds + kimura80)
    gremlins = filter_gremlins(gff_df, sp, gfams, use_genome_hits=False)
    print(f"  {len(gremlins)} Grem-LINE1 elements in {sp}")

    # ── Genome and scaffold sizes ─────────────────────────────────────────────
    bw = pyBigWig.open(str(H3K9_BW))
    bw_chroms = bw.chroms()
    if RHIMI59_FAI.exists():
        fai = pd.read_csv(RHIMI59_FAI, sep="\t", header=None,
                          names=["chrom","size","offset","bpline","byteline"])
        chrom_sizes = dict(zip(fai["chrom"], fai["size"]))
    else:
        chrom_sizes = dict(bw_chroms)

    grem_per_scaf = gremlins.groupby("chrom").size()

    # ── Per-scaffold classification ───────────────────────────────────────────
    rows = []
    for sc, n_grem in grem_per_scaf.items():
        if sc not in bw_chroms:
            continue
        size = chrom_sizes.get(sc, bw_chroms[sc])
        # H3K9me3 mean over whole scaffold (10 bins for CV)
        n_bins = max(1, min(10, size // 5_000))
        v = bw.stats(sc, 0, size, type="mean", nBins=n_bins)
        v = np.array([x if x is not None else 0.0 for x in v])
        h3k9_mean = float(np.mean(v))
        h3k9_cv   = (float(np.std(v) / max(np.mean(v), 0.001))
                     if n_bins > 1 else 0.0)

        sub_g = gremlins[gremlins["chrom"] == sc]["start"]
        grem_mid_rel = float(sub_g.median() / size) if len(sub_g) else 0.5

        size_kb = size / 1e3
        if size_kb < CENTRO_MAX_KB and h3k9_mean >= CENTRO_H3K9:
            scaf_type = "centromeric"
        elif n_grem >= MIN_GREM_SCAFFOLD:
            scaf_type = "arm"
        else:
            scaf_type = "arm_low"

        rows.append({
            "chrom": sc, "size_kb": size_kb, "n_gremlin": n_grem,
            "h3k9_mean": h3k9_mean, "h3k9_cv": h3k9_cv,
            "grem_mid_rel": grem_mid_rel, "scaffold_type": scaf_type,
        })

    scaf_df = pd.DataFrame(rows).sort_values("size_kb")
    scaf_df.to_csv(igv_dir / "rhimi59_scaffold_classification.tsv", sep="\t", index=False)

    arm_h3k9   = scaf_df.loc[scaf_df["scaffold_type"].isin(["arm","arm_low"]),
                              "h3k9_mean"].values
    centro_h3k9 = scaf_df.loc[scaf_df["scaffold_type"] == "centromeric",
                               "h3k9_mean"].values
    n_centro   = int((scaf_df["scaffold_type"] == "centromeric").sum())
    n_arm      = int(scaf_df["scaffold_type"].isin(["arm","arm_low"]).sum())

    mw_p = float("nan")
    if len(centro_h3k9) and len(arm_h3k9):
        _, mw_p = mannwhitneyu(centro_h3k9, arm_h3k9, alternative="greater")
    print(f"  Scaffold classification: {n_centro} centromeric, {n_arm} arm")
    print(f"  H3K9me3 mean — centromeric: {float(np.mean(centro_h3k9)):.2f} CPM  "
          f"arm: {float(np.mean(arm_h3k9)):.2f} CPM  MW p={mw_p:.4e}")

    # ── IGV track for best centromeric scaffold ───────────────────────────────
    centro_scafs = (scaf_df[scaf_df["scaffold_type"] == "centromeric"]
                   .nlargest(MAX_IGV_PANELS, "n_gremlin"))
    centro_tracks = []
    for _, row in centro_scafs.iterrows():
        sc = row["chrom"]
        if sc not in bw_chroms:
            continue
        size = int(chrom_sizes.get(sc, bw_chroms[sc]))
        bins = max(50, min(IGV_BINS, size // 200))
        v = bw.stats(sc, 0, size, type="mean", nBins=bins)
        h3k9v = np.array([x if x is not None else 0.0 for x in v])
        grem_in = gremlins[gremlins["chrom"] == sc].copy()
        centro_tracks.append({
            "chrom": sc, "win_start": 0, "win_end": size,
            "n_gremlin": int(row["n_gremlin"]),
            "h3k9_mean": float(row["h3k9_mean"]),
            "h3k9_vals": h3k9v, "gremlin_df": grem_in,
            "scaffold_type": "centromeric",
        })
        np.save(igv_dir / f"{sp}_{sc}_centro_h3k9.npy", h3k9v)

    # ── Full-scaffold profile for best arm scaffold ───────────────────────────
    arm_tracks = []
    best_arm = (scaf_df[scaf_df["scaffold_type"] == "arm"]
                .nlargest(min(2, MAX_IGV_PANELS - 1), "n_gremlin"))
    for _, row in best_arm.iterrows():
        sc = row["chrom"]
        if sc not in bw_chroms:
            continue
        size = int(chrom_sizes.get(sc, bw_chroms[sc]))
        v = bw.stats(sc, 0, size, type="mean", nBins=IGV_BINS)
        h3k9v = np.array([x if x is not None else 0.0 for x in v])
        grem_in = gremlins[gremlins["chrom"] == sc].copy()
        arm_tracks.append({
            "chrom": sc, "win_start": 0, "win_end": size,
            "n_gremlin": int(row["n_gremlin"]),
            "h3k9_mean": float(row["h3k9_mean"]),
            "h3k9_vals": h3k9v, "gremlin_df": grem_in,
            "scaffold_type": "arm",
        })
        np.save(igv_dir / f"{sp}_{sc}_arm_h3k9.npy", h3k9v)

    bw.close()
    print(f"  IGV data: {len(centro_tracks)} centromeric scaffold(s), "
          f"{len(arm_tracks)} arm scaffold(s)")

    return {
        "scaf_df":       scaf_df,
        "centro_h3k9":   centro_h3k9,
        "arm_h3k9":      arm_h3k9,
        "n_centro":      n_centro,
        "n_arm":         n_arm,
        "mw_p":          mw_p,
        "centro_tracks": centro_tracks,
        "arm_tracks":    arm_tracks,
    }


# ── ANALYSIS 3: Cross-species conservation ───────────────────────────────────

def analyse_cross_species(gfams, summary_df):
    """
    For MUCORALES_PANEL species: compute centromere count, orientation bias,
    K80 at centromere windows vs genome-wide.
    """
    rows = []
    for sp in MUCORALES_PANEL:
        row = summary_df[summary_df["species"] == sp]
        if row.empty:
            continue
        r = row.iloc[0]
        gp = gff_path(sp)
        if not gp.exists():
            continue

        gff_df   = load_gff_as_bed(gp)
        # Analysis 3 needs kimura80 values → use EarlGrey
        gremlins = filter_gremlins(gff_df, sp, gfams, use_genome_hits=False)
        k80_all  = gremlins["kimura80"].dropna().values

        k80_centro, k80_outside = [], []
        centro_bed = OUTDIR / "gremlin_clusters" / f"{sp}_centromes.bed"
        if centro_bed.exists() and len(gremlins):
            cb = pd.read_csv(centro_bed, sep="\t", header=None,
                             names=["chrom","win_start","win_end","n_grem"])
            for _, cr in cb.iterrows():
                mask_in  = ((gremlins["chrom"] == cr["chrom"]) &
                            (gremlins["start"] >= cr["win_start"]) &
                            (gremlins["end"]   <= cr["win_end"]))
                k80_centro.extend(gremlins[mask_in]["kimura80"].dropna().tolist())
                k80_outside.extend(gremlins[~mask_in]["kimura80"].dropna().tolist())

        genus = ("Mucor" if "Muc" in sp else
                 "Rhizopus" if "Rhi" in sp else "Other")
        rows.append({
            "species": sp, "genus": genus,
            "n_putative_centromeres": int(r.get("n_putative_centromeres", 0)),
            "n_gremlin_total":        int(r.get("n_gremlin_total", 0)),
            "pct_orient_biased":      (100.0 * r.get("n_orient_biased", 0) /
                                       max(1, r.get("n_putative_centromeres", 1))),
            "k80_mean_all":       float(np.mean(k80_all))    if len(k80_all)    else np.nan,
            "k80_mean_centromere":float(np.mean(k80_centro)) if k80_centro      else np.nan,
            "k80_mean_outside":   float(np.mean(k80_outside)) if k80_outside    else np.nan,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUTDIR / "cross_species_conservation.tsv", sep="\t", index=False)
    return df


# ── ANALYSIS 4: K80 divergence comparison ────────────────────────────────────

def analyse_k80_comparison(species_list, gfams):
    """Pool K80 values by TE group across all 63 MU species."""
    k80_dir = OUTDIR / "k80_comparison"
    k80_dir.mkdir(parents=True, exist_ok=True)

    groups = {k: [] for k in ["Grem-LINE1","Other LINE/I","LINE/L1",
                               "LTR/Gypsy","LTR/Copia","LTR/DIRS","DNA"]}

    for sp in species_list:
        gp = gff_path(sp)
        if not gp.exists():
            continue
        gff_df   = load_gff_as_bed(gp)
        # Analysis 4 needs kimura80 values from EarlGrey annotations
        gremlins = filter_gremlins(gff_df, sp, gfams, use_genome_hits=False)
        fams     = gfams.get(sp)

        gk80 = gremlins["kimura80"].dropna().values
        if len(gk80):
            groups["Grem-LINE1"].extend(gk80.tolist())

        li = gff_df[gff_df["te_class"].str.startswith("LINE/I")].copy()
        if fams:
            fl = {f.lower() for f in fams}
            li = li[~li["te_name"].str.lower().isin(fl)]
        groups["Other LINE/I"].extend(li["kimura80"].dropna().tolist())

        ll1 = gff_df[gff_df["te_class"].str.startswith("LINE/L1")]
        groups["LINE/L1"].extend(ll1["kimura80"].dropna().tolist())

        for cls, key in [("LTR/Gypsy","LTR/Gypsy"),("LTR/Copia","LTR/Copia"),
                         ("LTR/DIRS","LTR/DIRS")]:
            sub = gff_df[gff_df["te_class"] == cls]
            groups[key].extend(sub["kimura80"].dropna().tolist())

        dna = gff_df[gff_df["te_class"].str.startswith("DNA")]
        groups["DNA"].extend(dna["kimura80"].dropna().tolist())

    result = {}
    for k, v in groups.items():
        arr = np.array(v, dtype=float)
        arr = arr[(arr >= 0) & (arr <= 1.5)]
        result[k] = arr
        n = len(arr)
        m = float(np.mean(arr)) if n else 0
        print(f"  {k}: n={n:,}  mean={m:.3f}")
        safe = k.replace("/","_").replace(" ","_")
        np.save(k80_dir / f"k80_{safe}.npy", arr)

    ages_by_rate = {}
    for label, mu in MU_RATES.items():
        gk = result.get("Grem-LINE1", np.array([]))
        ages_by_rate[label] = (gk / (2.0 * mu)) / 1e6

    return result, ages_by_rate


# ── ANALYSIS 5: Integrated TF evidence on LINE/I ─────────────────────────────

def analyse_tf_integrated():
    """
    Load scripts 01-03 outputs and extract LINE/I TF evidence.
    Returns (peaks_df, signal_df, fimo_df).
    """
    TOP_N = 15
    peaks_df = signal_df = fimo_df = pd.DataFrame()

    if PEAK_MATRIX_TSV.exists():
        mat = pd.read_csv(PEAK_MATRIX_TSV, sep="\t", index_col=0)
        li_cols = [c for c in mat.columns
                   if str(c).startswith("LINE/I") and "Jockey" not in str(c)]
        if li_cols:
            mat["n_peaks_LINE_I"] = mat[li_cols].sum(axis=1)
            peaks_df = (mat["n_peaks_LINE_I"].reset_index()
                        .rename(columns={"index": "tf_label"})
                        .sort_values("n_peaks_LINE_I", ascending=False)
                        .head(TOP_N))

    if BW_SIGNAL_TSV.exists():
        bmat = pd.read_csv(BW_SIGNAL_TSV, sep="\t", index_col=0)
        li_rows = [r for r in bmat.index
                   if str(r).startswith("LINE/I") and "Jockey" not in str(r)]
        if li_rows:
            signal = bmat.loc[li_rows].mean(axis=0)
            signal_df = (signal.reset_index()
                                .rename(columns={"index": "tf_label", 0: "signal_LINE_I"})
                                .sort_values("signal_LINE_I", ascending=False)
                                .head(TOP_N))

    if FIMO_TSV.exists():
        f3 = pd.read_csv(FIMO_TSV, sep="\t")
        li_f = f3[f3["te_class"].str.startswith("LINE/I") &
                  ~f3["te_class"].str.contains("Jockey")].copy()
        sig = li_f[(li_f["fisher_pval"] < 0.05) &
                   (li_f["log2_odds_ratio"] > 0)].copy()
        if not sig.empty:
            fimo_df = (sig.groupby("motif_family")
                          .agg(log2OR=("log2_odds_ratio","max"),
                               pval=("fisher_pval","min"))
                          .reset_index()
                          .sort_values("log2OR", ascending=False)
                          .head(TOP_N))

    return peaks_df, signal_df, fimo_df


# ── FIGURES ───────────────────────────────────────────────────────────────────

def make_figure(summary_df, h3k9_result, cross_df, k80_groups, ages_by_rate,
                peaks_df, signal_df, fimo_df):
    fig = plt.figure(figsize=(22, 26))
    gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.45,
                            height_ratios=[0.85, 1.1, 1.0])
    ax_A = fig.add_subplot(gs[0, :])

    # Panel B layout: 2 stat axes + up to MAX_IGV_PANELS-2 IGV tracks
    n_centro_tracks = len(h3k9_result.get("centro_tracks", [])) if h3k9_result else 0
    n_arm_tracks    = len(h3k9_result.get("arm_tracks",    [])) if h3k9_result else 0
    n_b_slots = 2 + min(n_centro_tracks + n_arm_tracks, MAX_IGV_PANELS - 2)
    n_cols_b = 2
    n_rows_b = max(1, (n_b_slots + n_cols_b - 1) // n_cols_b)
    gs_B = gridspec.GridSpecFromSubplotSpec(n_rows_b, n_cols_b, subplot_spec=gs[1, 0],
                                            hspace=0.55, wspace=0.35)
    ax_C = fig.add_subplot(gs[1, 1])
    ax_D = fig.add_subplot(gs[2, 0])

    gs_E  = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2, 1], wspace=0.55)
    ax_E1 = fig.add_subplot(gs_E[0])
    ax_E2 = fig.add_subplot(gs_E[1])
    ax_E3 = fig.add_subplot(gs_E[2])

    for lbl, ax in zip("ACDE", [ax_A, ax_C, ax_D, ax_E1]):
        ax.text(-0.03, 1.04, lbl, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="top")

    # B label on first B cell
    ax_b0 = fig.add_subplot(gs_B[0, 0])
    ax_b0.text(-0.12, 1.06, "B", transform=ax_b0.transAxes,
               fontsize=13, fontweight="bold", va="top")
    ax_b0.set_visible(False)

    _panel_A(ax_A, summary_df)

    b_axes = []
    for i in range(n_b_slots):
        r, c = divmod(i, n_cols_b)
        b_axes.append(fig.add_subplot(gs_B[r, c]))
    _panel_B(b_axes, h3k9_result)


    _panel_C(ax_C, cross_df)
    _panel_D(ax_D, k80_groups, ages_by_rate)
    _panel_E(ax_E1, ax_E2, ax_E3, peaks_df, signal_df, fimo_df)

    fig.suptitle(
        "Grem-LINE1 as putative centromere markers in Mucoromycota — "
        "H3K9me3 validation, cross-species conservation & TF enrichment",
        fontsize=11, y=1.002)
    out = OUTDIR / "gremlin_centromere_analysis.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {out}")


def _panel_A(ax, df):
    ax.set_title("A  Grem-LINE1 elements and putative centromere count across 63 Mucoromycota",
                 fontsize=9, loc="left")
    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    df = df.sort_values("n_putative_centromeres", ascending=False).reset_index(drop=True)
    x  = np.arange(len(df))
    w  = 0.55

    colors = [CLADE_COLORS.get(CLADE_MAP.get(sp, "Mucorales"), CLADE_COLORS["Mucorales"])
              for sp in df["species"]]
    ax.bar(x, df["n_gremlin_total"], width=w, color=colors, alpha=0.75,
           label="Grem-LINE1 elements")

    ax2 = ax.twinx()
    ax2.plot(x, df["n_putative_centromeres"], "o-", color="black",
             ms=3, lw=1.2, label="Putative centromeres (any)")
    if "n_both_flanks" in df.columns:
        ax2.plot(x, df["n_both_flanks"], "s--", color="#D62728",
                 ms=3.5, lw=1.0, label="Both-flanks (high-conf.)")

    ax.set_xticks(x)
    ax.set_xticklabels(df["species"].str.replace("MU_",""), rotation=90, fontsize=5)
    ax.set_ylabel("Grem-LINE1 elements (bars)", fontsize=8, color=CLADE_COLORS["Mucorales"])
    ax2.set_ylabel(f"Putative centromeres (line)\n"
                   f"[scaffolds with >={MIN_GREM_SCAFFOLD} gremlin copies]", fontsize=8)
    ax.tick_params(axis="y", labelsize=7)
    ax2.tick_params(axis="y", labelsize=7)

    idx = df.index[df["species"] == "MU_Mucci3"].tolist()
    if idx:
        i = idx[0]
        ax.annotate("Mucor lusitn.\n~9 chr",
                    xy=(i, df.loc[i, "n_gremlin_total"]),
                    xytext=(max(0, i-8), df["n_gremlin_total"].max() * 0.80),
                    fontsize=6, ha="left", color="darkred",
                    arrowprops=dict(arrowstyle="-", color="darkred", lw=0.7))

    handles = [mpatches.Patch(color=v, label=k) for k, v in CLADE_COLORS.items()]
    ax2_handles = [
        plt.Line2D([0],[0], marker="o", color="black", ms=4, lw=1.2, label="Putative centromeres"),
    ]
    if "n_both_flanks" in df.columns:
        ax2_handles.append(
            plt.Line2D([0],[0], marker="s", color="#D62728", ms=4, lw=1.0,
                       ls="--", label="Both-flanks (high-conf.)")
        )
    ax.legend(handles=handles + ax2_handles, fontsize=6, loc="upper right")
    ax.set_xlabel("Species (sorted by putative centromere count)", fontsize=8)


def _panel_B(axes, h3k9_result):
    """
    Panel B: H3K9me3 scaffold-classification validation (MU_Rhimi59).
    axes[0]  — violin: H3K9me3 mean in centromeric scaffolds vs arm scaffolds
    axes[1]  — scatter: scaffold size vs H3K9me3 mean, colored by scaffold_type, sized by n_gremlin
    axes[2+] — IGV full-scaffold profiles (centromeric scaffold + arm scaffold)
    """
    if not axes or not h3k9_result:
        return

    scaf_df     = h3k9_result.get("scaf_df")
    if scaf_df is None or scaf_df.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=7)
        return

    centro_h3k9 = h3k9_result.get("centro_h3k9", np.array([]))
    arm_h3k9    = h3k9_result.get("arm_h3k9", np.array([]))
    mw_p        = h3k9_result.get("mw_p", float("nan"))
    n_centro    = h3k9_result.get("n_centro", 0)
    n_arm       = h3k9_result.get("n_arm", 0)
    centro_tracks = h3k9_result.get("centro_tracks", [])
    arm_tracks    = h3k9_result.get("arm_tracks", [])

    TYPE_COLORS = {"centromeric": "#E07B39", "arm": "#888888", "arm_low": "#BBBBBB"}
    p_fmt  = lambda p: f"{p:.2e}" if p < 0.001 else f"{p:.3f}"
    sig_txt = ("***" if mw_p < 0.001 else "**" if mw_p < 0.01
               else "*" if mw_p < 0.05 else "ns")

    # ── B1: violin H3K9me3 mean per scaffold class ────────────────────────────
    ax1 = axes[0]
    ax1.set_title("B  H3K9me3 scaffold classification (MU_Rhimi59)",
                  fontsize=8, loc="left")
    if len(arm_h3k9) and len(centro_h3k9):
        parts = ax1.violinplot([arm_h3k9, centro_h3k9], positions=[0, 1],
                               showmedians=True, showextrema=False)
        for pc, col in zip(parts["bodies"], ["#888888", "#E07B39"]):
            pc.set_facecolor(col); pc.set_alpha(0.65)
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.5)
        y_max = max(float(centro_h3k9.max()), float(arm_h3k9.max()), 1.0)
        ax1.annotate("", xy=(1, y_max * 1.04), xytext=(0, y_max * 1.04),
                     arrowprops=dict(arrowstyle="-", lw=1))
        ax1.text(0.5, y_max * 1.07, f"{sig_txt} p={p_fmt(mw_p)}",
                 ha="center", fontsize=6)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels([f"Arm scaffolds\n(n={n_arm})",
                          f"Centromeric\ncandidates\n(n={n_centro})"], fontsize=6.5)
    ax1.set_ylabel("H3K9me3 mean CPM (whole scaffold)", fontsize=7)
    ax1.tick_params(labelsize=6)
    ax1.text(0.97, 0.03,
             "Centromeric: < 500 kb, H3K9me3 ≥ 3 CPM\n"
             "Arm: ≥ 500 kb, ≥5 gremlins",
             transform=ax1.transAxes, ha="right", va="bottom", fontsize=5,
             bbox=dict(facecolor="lightyellow", alpha=0.8, edgecolor="none"))

    # ── B2: scatter scaffold size vs H3K9me3 mean ─────────────────────────────
    if len(axes) >= 2:
        ax2 = axes[1]
        colors_sc = [TYPE_COLORS.get(t, "#888888") for t in scaf_df["scaffold_type"]]
        sizes_sc  = np.clip(scaf_df["n_gremlin"] * 8, 12, 120)
        ax2.scatter(scaf_df["size_kb"], scaf_df["h3k9_mean"],
                    c=colors_sc, s=sizes_sc, alpha=0.75, linewidths=0)
        ax2.set_xscale("log")
        ax2.set_xlabel("Scaffold size (kb, log)", fontsize=7)
        ax2.set_ylabel("H3K9me3 mean CPM", fontsize=7)
        ax2.set_title("Size vs H3K9me3 per scaffold", fontsize=7)
        ax2.tick_params(labelsize=6)
        ax2.legend(handles=[
            mpatches.Patch(color="#E07B39", label="Centromeric candidate"),
            mpatches.Patch(color="#888888", label="Arm scaffold"),
        ], fontsize=5.5, loc="upper right", framealpha=0.7)
        ax2.axvline(500, color="grey", lw=0.8, ls="--", alpha=0.6)
        ax2.axhline(3.0, color="grey", lw=0.8, ls="--", alpha=0.6)

    # ── B3+: full-scaffold IGV profiles ───────────────────────────────────────
    all_tracks = centro_tracks[:1] + arm_tracks[:max(0, len(axes) - 3)]
    for ax, track in zip(axes[2:], all_tracks):
        win_s   = track["win_start"]; win_e = track["win_end"]
        h3k9v   = track["h3k9_vals"]; grem_df = track["gremlin_df"]
        sc_type = track["scaffold_type"]
        chrom   = track["chrom"]
        size_kb = (win_e - win_s) / 1e3

        x_pos = np.linspace(win_s, win_e, len(h3k9v))
        y_max = max(float(np.max(h3k9v)), 1.0)
        col_fill = "#E07B39" if sc_type == "centromeric" else "#2171B5"

        ax.fill_between(x_pos / 1e6, h3k9v, alpha=0.60, color=col_fill)
        ax.plot(x_pos / 1e6, h3k9v, color=col_fill, lw=0.6, alpha=0.9)

        if not grem_df.empty:
            y_base   = -y_max * 0.14; y_height = y_max * 0.10
            for _, el in grem_df.iterrows():
                col  = "#8B0000" if el["strand"] == "+" else "#4B0082"
                rect = mpatches.Rectangle(
                    (el["start"] / 1e6, y_base),
                    (el["end"] - el["start"]) / 1e6, y_height,
                    color=col, alpha=0.85, linewidth=0)
                ax.add_patch(rect)

        lbl = "centromeric scaffold" if sc_type == "centromeric" else "arm scaffold"
        ax.set_title(f"{chrom}  ({lbl})\n"
                     f"{size_kb:.0f} kb, {track['n_gremlin']} gremlins, "
                     f"H3K9me3={track['h3k9_mean']:.1f} CPM",
                     fontsize=6.5, pad=2)
        ax.set_xlabel("Position (Mb)", fontsize=6.5)
        ax.set_ylabel("H3K9me3 CPM", fontsize=6.5)
        ax.tick_params(labelsize=6)
        ax.set_xlim(win_s / 1e6, win_e / 1e6)
        ax.set_ylim(bottom=-y_max * 0.18)


def _panel_C(ax, df):
    ax.set_title("C  Centromere conservation — Mucor & Rhizopus species",
                 fontsize=9, loc="left")
    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return

    sp_labels = df["species"].str.replace("MU_","")
    x  = np.arange(len(df))
    w  = 0.5

    color_map = {"Mucor":"#E07B39", "Rhizopus":"#4878CF", "Other":"#888888"}
    colors = [color_map.get(g, "#888888") for g in df["genus"]]

    ax.bar(x, df["n_putative_centromeres"], width=w, color=colors, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(sp_labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Putative centromeres", fontsize=8)
    ax.tick_params(axis="y", labelsize=7)

    ax2 = ax.twinx()
    has_k80 = df[["k80_mean_centromere","k80_mean_all"]].notna().any(axis=1)
    if has_k80.any():
        ax2.plot(x[has_k80], df["k80_mean_centromere"][has_k80], "D",
                 color="#D62728", ms=5, zorder=5)
        ax2.plot(x[has_k80], df["k80_mean_all"][has_k80], "s",
                 color="#2CA02C", ms=5, zorder=5)
        ax2.set_ylabel("Mean K80 divergence", fontsize=8)
        ax2.tick_params(axis="y", labelsize=7)

    for i, row in df.reset_index(drop=True).iterrows():
        pct = row.get("pct_orient_biased", 0)
        if pct > 0:
            ax.text(i, row["n_putative_centromeres"] + 0.1, f"{pct:.0f}%",
                    ha="center", fontsize=5.5, color="darkred")

    genus_h = [mpatches.Patch(color=v, label=k)
               for k, v in color_map.items() if k in df["genus"].values]
    k80_h = ([plt.Line2D([0],[0],marker="D",color="w",mfc="#D62728",ms=5,
                          label="K80 centromere window"),
               plt.Line2D([0],[0],marker="s",color="w",mfc="#2CA02C",ms=5,
                          label="K80 genome-wide")]
              if has_k80.any() else [])
    ax.legend(handles=genus_h + k80_h, fontsize=6, loc="upper right", framealpha=0.7)
    ax.text(0.02, 0.02, "% = orientation-biased (binomial p<0.05)",
            transform=ax.transAxes, fontsize=5.5, color="darkred", va="bottom")


def _panel_D(ax, k80_groups, ages_by_rate):
    ax.set_title("D  K80 divergence: Grem-LINE1 vs other TE classes (63 Mucoromycota)",
                 fontsize=9, loc="left")
    if not k80_groups:
        ax.text(0.5, 0.5, "No K80 data", ha="center", va="center", transform=ax.transAxes)
        return

    labels  = list(k80_groups.keys())
    data    = [k80_groups[k] for k in labels]
    medians = [np.median(d) if len(d) else 0 for d in data]
    order   = np.argsort(medians)[::-1]
    labels  = [labels[i] for i in order]
    data    = [data[i] for i in order]

    vc = {"Grem-LINE1":"#E07B39","Other LINE/I":"#FF7F0E","LINE/L1":"#AEC7E8",
          "LTR/Gypsy":"#1F77B4","LTR/Copia":"#98DF8A","LTR/DIRS":"#2CA02C","DNA":"#BCBD22"}
    colors = [vc.get(k, "gray") for k in labels]

    parts = ax.violinplot(data, positions=range(len(labels)),
                          showmedians=True, showextrema=False)
    for pc, col in zip(parts["bodies"], colors):
        pc.set_facecolor(col)
        pc.set_alpha(0.7)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.5)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel("Kimura80 divergence vs consensus", fontsize=8)
    ax.set_ylim(0, 1.2)

    g1 = k80_groups.get("Grem-LINE1", np.array([]))
    g2 = k80_groups.get("Other LINE/I", np.array([]))
    if len(g1) > 2 and len(g2) > 2:
        _, pval = mannwhitneyu(g1, g2, alternative="two-sided")
        ax.text(0.02, 0.97, f"Grem vs other LINE/I: MW p={pval:.2e}",
                transform=ax.transAxes, ha="left", va="top", fontsize=7,
                bbox=dict(facecolor="lightyellow", alpha=0.8, edgecolor="none"))

    # Estimated age secondary axis (right)
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    k80_ticks = [0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0]
    fast_mu = list(MU_RATES.values())[1]
    slow_mu = list(MU_RATES.values())[0]
    af = [(k / (2 * fast_mu)) / 1e6 for k in k80_ticks]
    as_ = [(k / (2 * slow_mu)) / 1e6 for k in k80_ticks]
    ax2.set_yticks(k80_ticks)
    ax2.set_yticklabels([f"{a:.0f}-{b:.0f}" for a, b in zip(af, as_)], fontsize=5.5)
    ax2.set_ylabel("Estimated age Ma (fast-slow mu)", fontsize=7)


def _panel_E(ax1, ax2, ax3, peaks_df, signal_df, fimo_df):
    for ax in [ax1, ax2, ax3]:
        ax.tick_params(labelsize=6)

    # E1: peaks
    ax1.set_title("E  DAP peaks\non LINE/I", fontsize=8, loc="left")
    if not peaks_df.empty:
        top = peaks_df.head(15).sort_values("n_peaks_LINE_I")
        ax1.barh(range(len(top)), top["n_peaks_LINE_I"], color="#4878CF", alpha=0.85)
        ax1.set_yticks(range(len(top)))
        labels = top["tf_label"].str.replace(r"_\d+$","", regex=True)
        ax1.set_yticklabels(labels, fontsize=5.5)
        ax1.set_xlabel("# peaks on LINE/I", fontsize=7)
    else:
        ax1.text(0.5, 0.5, "No data", ha="center", va="center",
                 transform=ax1.transAxes, fontsize=8)

    # E2: bigwig signal
    ax2.set_title("BigWig\nsignal", fontsize=8, loc="left")
    if not signal_df.empty:
        ycol = [c for c in signal_df.columns if c != "tf_label"][0]
        top  = signal_df.head(15).sort_values(ycol)
        ax2.barh(range(len(top)), top[ycol], color="#E07B39", alpha=0.85)
        ax2.set_yticks(range(len(top)))
        labels = top["tf_label"].str.replace(r"_\d+$","", regex=True)
        ax2.set_yticklabels(labels, fontsize=5.5)
        ax2.set_xlabel("Mean CPM on LINE/I", fontsize=7)
    else:
        ax2.text(0.5, 0.5, "No data", ha="center", va="center",
                 transform=ax2.transAxes, fontsize=8)

    # E3: FIMO
    ax3.set_title("FIMO motifs\nlog2OR", fontsize=8, loc="left")
    if not fimo_df.empty:
        top    = fimo_df.head(15).sort_values("log2OR")
        colors = ["#D62728" if p < 0.01 else "#1F77B4" for p in top["pval"]]
        ax3.barh(range(len(top)), top["log2OR"], color=colors, alpha=0.85)
        ax3.set_yticks(range(len(top)))
        ax3.set_yticklabels(top["motif_family"], fontsize=5.5)
        ax3.axvline(0, color="black", lw=0.8)
        ax3.set_xlabel("log2 odds ratio (FIMO)", fontsize=7)
        ax3.legend(handles=[mpatches.Patch(color="#D62728", label="p<0.01"),
                             mpatches.Patch(color="#1F77B4", label="p<0.05")],
                   fontsize=5.5, loc="lower right")
    else:
        ax3.text(0.5, 0.5, "No significant\nenrichments", ha="center", va="center",
                 transform=ax3.transAxes, fontsize=8)


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Grem-LINE1 centromere analysis")
    parser.add_argument("--skip-clusters", action="store_true",
                        help="Re-use existing centromere prediction TSV")
    parser.add_argument("--skip-k80", action="store_true",
                        help="Re-use existing K80 .npy files")
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "gremlin_clusters").mkdir(exist_ok=True)
    (OUTDIR / "h3k9me3_igv").mkdir(exist_ok=True)
    (OUTDIR / "k80_comparison").mkdir(exist_ok=True)

    print("=== Loading Grem-LINE1 BLAST family assignments ===")
    gfams = load_gremlin_families()

    # Analysis 1
    summary_tsv = OUTDIR / "gremlin_clusters" / "centromere_prediction_summary.tsv"
    if not args.skip_clusters or not summary_tsv.exists():
        print("\n=== Analysis 1: Centromere peak prediction (63 MU species) ===")
        summary_df = analyse_centromere_peaks(MU_SPECIES, gfams)
    else:
        print(f"Loading existing summary: {summary_tsv}")
        summary_df = pd.read_csv(summary_tsv, sep="\t")

    mucci = summary_df[summary_df["species"] == "MU_Mucci3"]
    if not mucci.empty:
        n_c = int(mucci.iloc[0]["n_putative_centromeres"])
        print(f"\n  CALIBRATION -- MU_Mucci3 (CBSv3 = paper strain): "
              f"{n_c} putative centromeres (expected ~9 chromosomes)")

    # Analysis 2
    print("\n=== Analysis 2: H3K9me3 peak-centric validation (MU_Rhimi59) ===")
    h3k9_result = analyse_h3k9me3_igv(gfams)

    # Analysis 3
    print("\n=== Analysis 3: Cross-species conservation (Mucorales panel) ===")
    cross_df = analyse_cross_species(gfams, summary_df)

    # Analysis 4
    k80_npy = OUTDIR / "k80_comparison" / "k80_Grem-LINE1.npy"
    if not args.skip_k80 or not k80_npy.exists():
        print("\n=== Analysis 4: K80 divergence comparison (63 MU species) ===")
        k80_groups, ages_by_rate = analyse_k80_comparison(MU_SPECIES, gfams)
    else:
        print(f"Loading existing K80 arrays...")
        k80_groups = {}
        for k in ["Grem-LINE1","Other LINE/I","LINE/L1",
                  "LTR/Gypsy","LTR/Copia","LTR/DIRS","DNA"]:
            safe = k.replace("/","_").replace(" ","_")
            p = OUTDIR / "k80_comparison" / f"k80_{safe}.npy"
            if p.exists():
                k80_groups[k] = np.load(p)
        ages_by_rate = {lbl: (k80_groups.get("Grem-LINE1", np.array([]))
                               / (2*mu)) / 1e6
                        for lbl, mu in MU_RATES.items()}

    # Analysis 5
    print("\n=== Analysis 5: Integrated TF evidence on LINE/I ===")
    peaks_df, signal_df, fimo_df = analyse_tf_integrated()

    # Figure
    print("\n=== Generating figure ===")
    make_figure(summary_df, h3k9_result, cross_df, k80_groups, ages_by_rate,
                peaks_df, signal_df, fimo_df)
    print("Done.")


if __name__ == "__main__":
    main()
