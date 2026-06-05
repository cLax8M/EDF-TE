#!/usr/bin/env python3
"""
create_source_data.py
---------------------
Reads pipeline result TSVs and generates the Source Data Excel workbook
EDF-TE_source_data.xlsx.  Each worksheet corresponds to one figure panel
as described in MANUSCRIPT_SUPPLEMENTS.md Part III.

Usage:
    python create_source_data.py [--outdir <dir>]

Paths are resolved from config.sh variables, which can be overridden by
environment variables of the same name.

Required packages:
    pandas >= 2.0, openpyxl >= 3.1, scipy, numpy
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")


# ---------------------------------------------------------------------------
# Path configuration  (mirrors config.sh; override with env vars)
# ---------------------------------------------------------------------------

BASE_DIR = Path(os.environ.get("BASE_DIR", "/data/carlos/EDF-TE"))
RESULTS_DIR = Path(os.environ.get("RESULTS_DIR", BASE_DIR / "TE_analysis" / "DAM_search" / "results"))
EXP_RESULTS = Path(os.environ.get("EXP_RESULTS", BASE_DIR / "TE_analysis" / "expansion_analysis" / "results"))
WGD_DIR = Path(os.environ.get("WGD_DIR", "/legserv/Temp/Carlos/EDF-TE/WGD/Modern_wgd2_pipeline/results/paper_no_filter"))
WGD_BP_DIR = WGD_DIR / "te_breakpoint_enrichment"
TE_ANNOT_DIR = Path(os.environ.get("TE_ANNOT_DIR", BASE_DIR / "TE_analysis" / "annotated"))
TF_TE_RESULTS = Path(os.environ.get("TF_TE_RESULTS", "/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis"))
S2_DIR = Path(os.environ.get("S2_DIR", BASE_DIR / "EDF-TE" / "EDF-TE_complete" / "section2_te_landscape"))
S3_DIR = Path(os.environ.get("S3_DIR", BASE_DIR / "EDF-TE" / "EDF-TE_complete" / "section3_te_epigenetics"))
S4_DIR = Path(os.environ.get("S4_DIR", BASE_DIR / "EDF-TE" / "EDF-TE_complete" / "section4_te_cooption"))
S6_DIR = Path(os.environ.get("S6_DIR", BASE_DIR / "EDF-TE" / "EDF-TE_complete" / "section6_htt"))
HTT_DIR = Path(os.environ.get("HTT_DIR", "/legserv/Temp/Carlos/EDF-TE/results/HTT"))
TE_LANDSCAPE_DIR = Path(os.environ.get("TE_LANDSCAPE_DIR", "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape"))
PANGENOME_DIR = Path(os.environ.get("PANGENOME_DIR", "/legserv/Temp/Carlos/EDF-TE/results/PanGenome/results"))
GENOMIC_FEAT_DIR = Path(os.environ.get("GENOMIC_FEAT_DIR", "/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures_V3"))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

UNIT_ROW_MARKER = "__UNITS__"


def _add_units(df: pd.DataFrame, units: dict) -> pd.DataFrame:
    """Insert a units/description row immediately below the header."""
    unit_row = {c: units.get(c, "") for c in df.columns}
    unit_df = pd.DataFrame([unit_row], columns=df.columns)
    return pd.concat([unit_df, df], ignore_index=True)


def safe_read(path: Path, **kwargs) -> pd.DataFrame | None:
    """Read TSV, return None and warn if missing."""
    if not path.exists():
        print(f"  [WARN] Missing input: {path}", file=sys.stderr)
        return None
    return pd.read_csv(path, sep="\t", **kwargs)


def write_sheet(writer: pd.ExcelWriter, df: pd.DataFrame, sheet_name: str,
                units: dict | None = None) -> None:
    """Write df to sheet, optionally inserting a units row after the header."""
    if df is None or df.empty:
        placeholder = pd.DataFrame([["Data not available"]])
        placeholder.to_excel(writer, sheet_name=sheet_name, index=False, header=False)
        return
    if units:
        df = _add_units(df, units)
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    # Auto-fit column widths
    ws = writer.sheets[sheet_name]
    for col_cells in ws.columns:
        max_len = max(len(str(cell.value)) if cell.value is not None else 0
                      for cell in col_cells)
        ws.column_dimensions[col_cells[0].column_letter].width = min(max_len + 2, 50)


# ---------------------------------------------------------------------------
# Sheet builders
# ---------------------------------------------------------------------------

def build_fig1b(writer):
    """Fig1B_pangenome: pan-genome accumulation curve data per taxonomic group."""
    # Primary: summary_final_counts.tsv (one row per clade: n_species, pangenome, core, shell, unique)
    summary_path = PANGENOME_DIR / "summary_final_counts.tsv"
    df = safe_read(summary_path)
    # Also gather per-group accumulation curves
    group_dir = PANGENOME_DIR / "per_group"
    frames = []
    if group_dir.is_dir():
        for tsv in sorted(group_dir.glob("*_pangenome.tsv")):
            gdf = safe_read(tsv)
            if gdf is not None:
                gdf["group_code"] = tsv.stem.replace("_pangenome", "")
                frames.append(gdf)
    if frames:
        acc = pd.concat(frames, ignore_index=True)
        if df is not None:
            # Merge group metadata onto accumulation curves
            meta_cols = [c for c in df.columns if c in ("group_code", "group_name", "n_species")]
            if meta_cols:
                acc = acc.merge(df[[*meta_cols]], on="group_code", how="left")
        df = acc
    units = {
        "group_code": "clade code",
        "group_name": "clade full name",
        "n": "number of species sampled (accumulation step)",
        "pangenome_mean": "mean HOG count (pan-genome)",
        "pangenome_sd": "s.d. HOG count (pan-genome)",
        "core_mean": "mean HOG count (core)",
        "core_sd": "s.d. HOG count (core)",
        "shell_mean": "mean HOG count (shell)",
        "shell_sd": "s.d. HOG count (shell)",
        "unique_mean": "mean HOG count (unique)",
        "unique_sd": "s.d. HOG count (unique)",
    }
    write_sheet(writer, df, "Fig1B_pangenome", units)


def build_fig1c(writer):
    """Fig1C_genomic_features: per-species genomic feature metrics (exon/intron content)."""
    path = GENOMIC_FEAT_DIR / "ALL_SPECIES_ALL_SUMMARY_V4.txt"
    df = safe_read(path)
    if df is None:
        alt = GENOMIC_FEAT_DIR / "ALL_SPECIES_ALL_SUMMARY.txt"
        df = safe_read(alt)
    if df is None:
        write_sheet(writer, None, "Fig1C_genomic_features")
        return
    # Rename Species column to species for consistency
    if "Species" in df.columns:
        df = df.rename(columns={"Species": "species"})
    # Compute percentages from raw counts
    if "Total_Exons" in df.columns and "Genome_Size_bp" in df.columns:
        df["genome_size_mb"] = df["Genome_Size_bp"] / 1e6
    units = {
        "species": "species code",
        "Genome_Size_bp": "bp",
        "genome_size_mb": "Mb",
        "Total_Genes": "count",
        "Total_Introns": "count",
        "Total_Exons": "count",
        "Avg_Exons_Per_Gene": "mean exons per gene",
        "Avg_Isoforms_Per_Gene": "mean isoforms per gene",
        "Avg_Introns_All": "mean introns per multi-exon gene",
    }
    write_sheet(writer, df, "Fig1C_genomic_features", units)


def build_fig1d(writer):
    """Fig1D_synteny: pairwise anchor counts."""
    path = WGD_DIR / "synteny_summary.tsv"
    df = safe_read(path)
    units = {
        "species_id": "JGI species code",
        "syn_completed": "bool",
        "n_multiplicons": "count",
        "n_anchorpoints": "count",
        "n_segments": "count",
        "n_files_syn": "count",
    }
    write_sheet(writer, df, "Fig1D_synteny", units)


def build_fig2a(writer):
    """Fig2A_te_content: per-species TE bp per subclass (aggregated EarlGrey RepeatLandscape)."""
    # Primary: aggregate EarlGrey per-species summary_table.tsv files
    # Path pattern: TE_ANNOT_DIR/{species}_EarlGrey/{species}_RepeatLandscape/{species}_summary_table.tsv
    rows = []
    if TE_ANNOT_DIR.is_dir():
        for sp_dir in sorted(TE_ANNOT_DIR.iterdir()):
            if not sp_dir.is_dir():
                continue
            sp = sp_dir.name  # e.g. MU_Pilumb1_EarlGrey
            sp_code = sp.replace("_EarlGrey", "")
            # Try RepeatLandscape summary_table.tsv first
            rl_dir = sp_dir / f"{sp_code}_RepeatLandscape"
            stab = rl_dir / f"{sp_code}_summary_table.tsv"
            if not stab.exists():
                # fallback: divergence_summary_table in summaryFiles
                sf_dir = sp_dir / f"{sp_code}_summaryFiles"
                stab = sf_dir / f"{sp_code}_divergence_summary_table.tsv"
            if not stab.exists():
                continue
            sdf = safe_read(stab)
            if sdf is None:
                continue
            # Aggregate total_bp per subclass
            if "subclass" in sdf.columns and "total_bp" in sdf.columns:
                agg = sdf.groupby("subclass")["total_bp"].sum().reset_index()
                agg["species"] = sp_code
                rows.append(agg)
    if rows:
        long_df = pd.concat(rows, ignore_index=True)
        # Pivot to wide: species × subclass
        wide = long_df.pivot_table(index="species", columns="subclass",
                                    values="total_bp", aggfunc="sum", fill_value=0)
        wide = wide.reset_index()
        wide.columns.name = None
        df = wide
    else:
        df = None
    units = {
        "species": "species code",
        "[TE subclass]": "total bp masked by EarlGrey RepeatMasker (per TE subclass)",
    }
    write_sheet(writer, df, "Fig2A_te_content", units)


def build_fig2b(writer):
    """Fig2B_pca_composition: PCoA of TE class composition (Bray-Curtis) per species."""
    path = TE_LANDSCAPE_DIR / "G01_pcoa_coords.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "clade": "18-clade classification",
        "PCo1": "principal coordinate 1 (Bray-Curtis; arbitrary units)",
        "PCo2": "principal coordinate 2 (arbitrary units)",
        "PCo3": "principal coordinate 3 (arbitrary units)",
    }
    write_sheet(writer, df, "Fig2B_pca_composition", units)


def build_fig2c(writer):
    """Fig2C_kimura_ri: Recency Index per species (input to clade boxplot)."""
    path = TE_LANDSCAPE_DIR / "C01_kimura_species.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "recency_index": "fraction [0,1] — proportion of TE bp in youngest Kimura quartile",
        "burst_score": "burst score (0=ancient, 1=very recent)",
        "dominant_superfamily": "most abundant TE superfamily by bp",
        "mean_div_weighted": "coverage-weighted mean Kimura divergence",
        "n_superfamilies": "count of distinct annotated superfamilies",
        "total_te_bp": "total TE bp (EarlGrey)",
    }
    write_sheet(writer, df, "Fig2C_kimura_ri", units)


def build_fig2d(writer):
    """Fig2D_kimura_summary: mean Kimura divergence per TE superfamily x species."""
    path = TE_LANDSCAPE_DIR / "C01_kimura_superfamily.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "superfamily": "TE superfamily (EarlGrey annotation)",
        "total_bp": "total bp of this superfamily",
        "pct_of_TE": "% of total TE bp",
        "recency_index": "fraction [0,1]",
        "mean_div_weighted": "coverage-weighted mean Kimura divergence",
        "n_families": "number of families within this superfamily",
    }
    write_sheet(writer, df, "Fig2D_kimura_summary", units)


def build_fig3a(writer):
    """Fig3A_dnmt: DNMT presence/absence per species + TE-embedded DNMT screening."""
    # Primary: per-species methyltransferase co-option summary
    path = RESULTS_DIR / "methylase_cooption_summary.tsv"
    df = safe_read(path)
    # Supplement with TE-embedded DNMT screening results
    te_dnmt_path = Path("/legserv/Temp/Carlos/EDF-TE/results/Co-option/te_dnmt_search/te_with_dnmt.tsv")
    te_dnmt = safe_read(te_dnmt_path)
    if te_dnmt is not None and df is None:
        df = te_dnmt
    units = {
        "domain": "methyltransferase domain name (DAM, DNMT, MTA70)",
        "pfam": "Pfam accession(s)",
        "function": "biological function",
        "n_species_with_domain": "count of species with this domain",
        "n_species_TE_context": "count of species where domain found in TE ORF",
        "n_species_fused_TE": "count of species with host gene fused to TE structural domain",
        "TE_class_involved": "TE classes harbouring the domain",
        "co-option_type": "co-option mechanism description",
        "fused_species": "species codes with fused proteins",
        "clades": "clades represented",
        "note": "additional notes",
        # te_with_dnmt.tsv columns (if using alternate source)
        "te_id": "TE consensus identifier (EarlGrey)",
        "has_C5_DNMT": "True if TE ORF encodes C5-DNMT domain (PF00145/PF17980)",
        "has_TE_core": "True if TE consensus encodes a TE structural domain",
        "pfam_domains": "semicolon-separated Pfam domain hits in this TE ORF",
    }
    write_sheet(writer, df, "Fig3A_dnmt", units)


def build_fig3b(writer):
    """Fig3B_methylation: per-copy 5mC and Kimura divergence."""
    path = BASE_DIR / "TE_analysis" / "WGBS" / "te_methylation_per_copy.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "te_copy_id": "EarlGrey copy identifier",
        "te_class": "TE class",
        "kimura_div": "K2P distance",
        "kimura_quartile": "Q1-Q4",
        "pct_methylated": "% [0,100]",
        "n_cpg_covered": "count of CpGs with coverage >= 5x",
    }
    write_sheet(writer, df, "Fig3B_methylation", units)


def build_fig3c(writer):
    """Fig3C_chromatin: nearest-gene RPKM per chromatin category."""
    path = BASE_DIR / "TE_analysis" / "ChIP" / "te_chromatin_gene_rpkm.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "te_copy_id": "copy identifier",
        "te_class": "TE class",
        "chromatin_category": "unmethyl_noH3K9 | methyl_noH3K9 | unmethyl_H3K9me3 | methyl_H3K9me3",
        "nearest_gene_rpkm": "RPKM",
    }
    write_sheet(writer, df, "Fig3C_chromatin", units)


def build_fig3d(writer):
    """Fig3D_nucleosome: meta-TE dyad probability (150 bins)."""
    path = BASE_DIR / "TE_analysis" / "MNase" / "meta_te_dyad_profile.tsv"
    df = safe_read(path)
    units = {
        "te_class": "TE class",
        "bin_index": "1–150 (position relative to TE midpoint; ±2 kb window)",
        "mean_dyad_prob": "mean dyad probability (DANPOS3 output)",
        "n_copies_contributing": "count",
    }
    write_sheet(writer, df, "Fig3D_nucleosome", units)


def build_fig3e(writer):
    """Fig3E_srna_rpm: sRNA RPM per superfamily."""
    path = BASE_DIR / "TE_analysis" / "sRNA" / "srna_te_rpm.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "te_superfamily": "TE superfamily",
        "rpm": "reads per million",
        "total_mapped_reads": "count",
    }
    write_sheet(writer, df, "Fig3E_srna_rpm", units)


def build_fig3f(writer):
    """Fig3F_srna_h3k9: per-family sRNA RPM vs H3K9me3 RPKM."""
    path = BASE_DIR / "TE_analysis" / "sRNA" / "srna_vs_h3k9_per_family.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "te_superfamily": "TE superfamily",
        "rpm": "sRNA reads per million",
        "mean_h3k9me3_rpkm": "mean H3K9me3 RPKM over copies of this family",
    }
    write_sheet(writer, df, "Fig3F_srna_h3k9", units)


def build_fig2e(writer):
    """Fig2E_proximity_scatter: gene-TE proximity scatter for 41 candidate domain families."""
    path = EXP_RESULTS / "figures" / "labeled_domains_metrics.tsv"
    df = safe_read(path)
    if df is not None:
        keep = [c for c in [
            "pfam", "name", "target_clade", "mechanism", "expansion_type",
            "mw_pval_copynumber", "n_pack_element_candidates", "n_fusion_proteins",
            "n_te_consensus_with_domain", "n_copies_symmetric_TE_flanking",
            "prox_mw_pval", "pct5kb_clade", "fold_change",
        ] if c in df.columns]
        df = df[keep]
    units = {
        "pfam": "Pfam accession",
        "name": "domain name",
        "target_clade": "clade with significant expansion",
        "mechanism": "expansion mechanism classification",
        "mw_pval_copynumber": "Mann-Whitney U p (copy number; expanded vs rest)",
        "n_pack_element_candidates": "gene copies fully embedded in TE (pack elements)",
        "n_fusion_proteins": "proteins with both host and TE domain in same ORF",
        "prox_mw_pval": "Mann-Whitney U p (proximity to nearest TE vs background)",
        "pct5kb_clade": "% domain-containing genes within 5 kb of a TE (clade)",
        "fold_change": "fold change in copy number vs non-expanded clades",
    }
    write_sheet(writer, df, "Fig2E_proximity_scatter", units)


def build_fig2f(writer):
    """Fig2F_proximity_flanks: TE subfamily in symmetric flanks of expanded genes."""
    path = EXP_RESULTS / "proximity" / "proximity_summary.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "clade": "18-clade",
        "pfam": "Pfam accession",
        "domain_name": "domain name",
        "n_copies": "total gene copies of this domain in this species",
        "median_dist_target": "median bp to nearest TE (domain-containing genes)",
        "median_dist_bg": "median bp to nearest TE (background genes)",
        "pct_nearest_is_te": "% domain genes with a TE as nearest neighbour",
        "bg_pct_nearest_is_te": "% background genes with a TE as nearest neighbour",
        "mw_pval_vs_bg": "Mann-Whitney U p (distance to TE; domain vs background)",
    }
    write_sheet(writer, df, "Fig2F_proximity_flanks", units)


def build_fig4b(writer):
    """Fig4B_cooption_heatmap: fraction co-opted per clade × domain."""
    path = EXP_RESULTS / "clade_expansion_table.tsv"
    df = safe_read(path)
    if df is not None:
        # Select the columns relevant for the figure
        keep = [c for c in ["pfam", "domain_desc", "clade", "clade_name",
                             "expansion_class", "n_clade_species", "n_positive",
                             "pct_positive", "median_copies_positive",
                             "fold_vs_rest", "in_master"] if c in df.columns]
        df = df[keep]
        if "n_clade_species" in df.columns and "n_positive" in df.columns:
            df["frac_species_with_expansion"] = df["n_positive"] / df["n_clade_species"]
    units = {
        "pfam": "Pfam accession",
        "domain_desc": "domain description",
        "clade": "18-clade code",
        "clade_name": "clade full name",
        "expansion_class": "cladal | species-specific",
        "n_clade_species": "total species in clade",
        "n_positive": "species with expansion",
        "pct_positive": "% [0,100]",
        "median_copies_positive": "median copy number in expanded species",
        "fold_vs_rest": "fold change vs non-expanded species",
        "frac_species_with_expansion": "fraction [0,1]",
    }
    write_sheet(writer, df, "Fig4B_cooption_heatmap", units)


def build_fig4c(writer):
    """Fig4C_chimte: chimeric event counts per TE class per species."""
    path = BASE_DIR / "EDF-TE" / "Cooption" / "chimte_summary.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "clade": "18-clade",
        "te_class": "TE class",
        "n_chimeric_events": "count of chimeric read-pair junctions",
    }
    write_sheet(writer, df, "Fig4C_chimte", units)


def build_fig4d(writer):
    """Fig4D_tf_te: TF × TE class enrichment (DAP-seq bigwig signal, Rhizopus microsporus)."""
    # Real data: bw_TE_signal_class_mean.tsv — te_class × 58 TF columns (mean bigwig log2 enrichment)
    path = TF_TE_RESULTS / "02_bigwig_signal" / "bw_TE_signal_class_mean.tsv"
    df = safe_read(path)
    if df is None:
        # Fallback: fimo enrichment TSV
        path2 = TF_TE_RESULTS / "03_fimo_enrichment" / "fimo_enrichment.tsv"
        df = safe_read(path2)
    units = {
        "te_class": "TE class / superfamily (EarlGrey annotation)",
        "[TF_name]": "mean bigwig signal enrichment over TE copies (log2 FC vs input; 58 DAP-seq TFs, R. microsporus ATCC11559)",
    }
    write_sheet(writer, df, "Fig4D_tf_te", units)


def build_fig4e(writer):
    """Fig4E_fusions: fusion gene counts per clade × fusion class."""
    path = EXP_RESULTS / "fusion_summary.tsv"
    df = safe_read(path)
    if df is not None:
        # fusion_summary has: pfam, domain_name, n_proteins_with_domain,
        # n_fusions_with_TE, n_species_with_fusion, fusion_rate
        # Also try protein_fusion_screen for more detail
        screen = safe_read(EXP_RESULTS / "protein_fusion_screen.tsv")
        if screen is not None:
            df = screen
    units = {
        "pfam": "Pfam accession of host domain",
        "domain_name": "domain name",
        "n_proteins_with_domain": "count",
        "n_fusions_with_TE": "proteins with TE domain co-occurrence",
        "n_species_with_fusion": "count",
        "fusion_rate": "n_fusions / n_proteins_with_domain",
    }
    write_sheet(writer, df, "Fig4E_fusions", units)


def build_fig4f(writer):
    """Fig4F_flywch_trees: FLYWCH x TE-RT fusion protein metadata (Fig4F)."""
    path = EXP_RESULTS / "figures" / "labeled_domains_metrics.tsv"
    df = safe_read(path)
    if df is not None:
        # Filter to FLYWCH domain (PF04500)
        if "pfam" in df.columns:
            df = df[df["pfam"] == "PF04500"].copy()
    units = {
        "pfam": "Pfam accession (PF04500 = FLYWCH)",
        "name": "domain name",
        "target_clade": "NE = Neocallimastigomycota",
        "mechanism": "expansion mechanism",
        "n_fusion_proteins": "proteins with FLYWCH + TE structural domain",
        "n_pack_element_candidates": "FLYWCH gene copies embedded in TE",
        "prox_mw_pval": "Mann-Whitney U p (proximity to TE)",
    }
    write_sheet(writer, df, "Fig4F_flywch_trees", units)


def build_fig4g(writer):
    """Fig4G_coth_trees: CotH x TE fusion protein metadata (Fig4G)."""
    path = EXP_RESULTS / "figures" / "labeled_domains_metrics.tsv"
    df = safe_read(path)
    if df is not None:
        if "pfam" in df.columns and "name" in df.columns:
            # CotH domain: PF08757 (spore coat protein, Mucorales-expanded)
            mask = df["pfam"].str.startswith("PF08757", na=False) | \
                   df["name"].str.contains("CotH", case=False, na=False)
            df = df[mask].copy()
    units = {
        "pfam": "Pfam accession (PF08757 = spore coat protein CotH)",
        "target_clade": "clade with CotH expansion",
        "mechanism": "expansion mechanism",
        "n_fusion_proteins": "proteins with CotH + TE structural domain",
        "n_pack_element_candidates": "CotH gene copies embedded in TE",
        "prox_mw_pval": "Mann-Whitney U p (proximity to TE)",
    }
    write_sheet(writer, df, "Fig4G_coth_trees", units)


def build_fig4h(writer):
    """Fig4H_dam_rt_tree: DAM+RT protein metadata for circular phylogeny."""
    dam_path = RESULTS_DIR / "DAM_TE_summary_table_annotated.tsv"
    df = safe_read(dam_path)
    if df is None:
        # Fallback: phase2 master table
        df = safe_read(RESULTS_DIR / "phase2_master_table.tsv")
    if df is not None:
        cols_keep = ["species", "te_class", "source", "DAM_evalue", "RT_evalue",
                     "rep_orf", "orf_length_aa", "domain_architecture"]
        df = df[[c for c in cols_keep if c in df.columns]]
    # Also try to attach tree annotation if available
    annot_path = RESULTS_DIR / "phylo" / "DAM_full_annot.tsv"
    annot = safe_read(annot_path)
    if annot is not None and df is not None and "species" in annot.columns:
        df = df.merge(annot, on="species", how="left", suffixes=("", "_annot"))
    units = {
        "species": "species code",
        "te_class": "TE class of DAM-embedding element (LTR/DIRS or LTR/Gypsy)",
        "source": "Phase1 | Phase2 discovery phase",
        "DAM_evalue": "HMMER E-value vs PF05869 (DAM methyltransferase)",
        "RT_evalue": "HMMER E-value for reverse transcriptase domain (RVT_1)",
        "rep_orf": "RepBase consensus ORF label",
        "orf_length_aa": "TE-embedded ORF length (aa)",
        "domain_architecture": "InterProScan domain chain (e.g. RVT_1 -> Dam)",
    }
    write_sheet(writer, df, "Fig4H_dam_rt_tree", units)


def build_fig5a(writer):
    """Fig5A_paralogs: colinear paralog count + anchor-Ks per WGD species."""
    # --- Part 1: synteny / classification summary ---
    path = WGD_DIR / "wgd_candidates_integrated.tsv"
    df = safe_read(path)
    if df is None:
        syn = safe_read(WGD_DIR / "synteny_summary.tsv")
        wgd = safe_read(WGD_DIR / "wgd_classification.tsv")
        if syn is not None and wgd is not None:
            df = syn.merge(wgd, on="species_id", how="outer")
    if df is not None:
        if "species_id" in df.columns:
            df = df.rename(columns={"species_id": "species"})
        cols = ["species", "n_anchorpoints", "n_multiplicons", "n_segments",
                "classification", "confidence", "dup_rate_per_eligible_tree",
                "has_dup_signal", "has_synteny", "event_timing",
                "z_multiplicons", "z_segments", "max_peak_z"]
        df = df[[c for c in cols if c in df.columns]]
        df = df.rename(columns={"n_anchorpoints": "n_colinear_paralog_pairs",
                                 "classification": "wgd_classification"})

    # --- Part 2: aggregate per-species anchor-Ks files ---
    ks_dir = WGD_DIR / "ks_analysis"
    ks_frames = []
    if ks_dir.is_dir():
        for sp_dir in sorted(ks_dir.iterdir()):
            if not sp_dir.is_dir():
                continue
            sp = sp_dir.name
            anchor_f = sp_dir / f"{sp}.anchor_ks.tsv"
            if not anchor_f.exists():
                # try any tsv in dir
                tsv_files = list(sp_dir.glob("*.tsv"))
                anchor_f = tsv_files[0] if tsv_files else None
            if anchor_f and anchor_f.exists():
                kdf = safe_read(anchor_f)
                if kdf is not None:
                    kdf["species"] = sp
                    # rename dS/Ks columns
                    kdf = kdf.rename(columns={"dS": "Ks", "dN": "Ka"})
                    ks_frames.append(kdf)
    if ks_frames:
        ks_all = pd.concat(ks_frames, ignore_index=True)
        ks_all = ks_all[pd.to_numeric(ks_all.get("Ks", pd.Series(dtype=float)),
                                      errors="coerce").between(0.001, 3.0)]
        # Compute per-species Ks stats
        if "Ks" in ks_all.columns and "species" in ks_all.columns:
            ks_stats = (ks_all.groupby("species")["Ks"]
                        .agg(n_ks_pairs="count",
                             ks_median="median",
                             ks_mean="mean",
                             ks_mode_approx=lambda x: float(
                                 x.value_counts(bins=60).idxmax().mid))
                        .reset_index())
            if df is not None and "species" in df.columns:
                df = df.merge(ks_stats, on="species", how="left")
            else:
                df = ks_stats

    units = {
        "species": "WGD-candidate species code",
        "n_colinear_paralog_pairs": "anchor gene pairs in wgd syn multiplicons",
        "n_multiplicons": "collinear multiplicons (paralogons)",
        "n_segments": "collinear segments",
        "wgd_classification": "strong | uncertain | weak | no_WGD",
        "dup_rate_per_eligible_tree": "OrthoFinder normalised duplication rate",
        "has_dup_signal": "bool — OrthoFinder duplication peak detected",
        "has_synteny": "bool — wgd syn multiplicons > 0",
        "event_timing": "ancient_deep | recent_terminal | mixed_deep_and_recent",
        "z_multiplicons": "robust Z-score of n_multiplicons",
        "z_segments": "robust Z-score of n_segments",
        "max_peak_z": "max robust Z-score of duplication rate across nodes",
        "n_ks_pairs": "anchor-Ks pairs used for age estimation",
        "ks_median": "median anchor-Ks (paralogous divergence)",
        "ks_mean": "mean anchor-Ks",
        "ks_mode_approx": "approximate Ks mode (histogram peak, 60 bins)",
    }
    write_sheet(writer, df, "Fig5A_paralogs", units)


def build_fig5b(writer):
    """Fig5B_syntenic_blocks: syntenic block structure and TE content."""
    path = WGD_DIR / "synteny_summary.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "n_breakpoints": "count",
        "obs_te_bp_per_window": "mean bp of TE per ±5 kb breakpoint window",
        "bg_te_bp_per_window": "mean bp of TE per random window (n=100)",
        "bg_std": "std of random window TE bp",
        "z_score": "(obs − bg_mean) / bg_std",
        "p_value_permutation": "permutation p-value (100 random sets)",
        "log2_fold_change": "log2(obs / bg_mean)",
        "window_bp": "bp",
        "n_permutations": "count",
    }
    write_sheet(writer, df, "Fig5B_syntenic_blocks", units)


def build_fig5c(writer):
    """Fig5C_te_breakpoints: TE enrichment at WGD breakpoints — log2FC + Wilcoxon BH-q."""
    # Aggregate per-species per_superfamily results from te_breakpoint_enrichment/
    frames = []
    if WGD_BP_DIR.is_dir():
        for tsv in sorted(WGD_BP_DIR.glob("*.tsv")):
            sdf = safe_read(tsv)
            if sdf is None:
                continue
            sp = tsv.stem
            sdf["species"] = sp
            frames.append(sdf)
    if frames:
        df = pd.concat(frames, ignore_index=True)
    else:
        # Fallback: try a single summary file
        df = safe_read(WGD_BP_DIR / "enrichment_summary.tsv")
    if df is None:
        write_sheet(writer, None, "Fig5C_te_breakpoints")
        return
    # Ensure species column exists and select informative columns
    cols = [c for c in [
        "species", "superfamily", "te_class", "obs_te_bp_per_window",
        "bg_te_bp_per_window", "log2_fold_change", "p_value_permutation",
        "wilcoxon_p", "BH_q", "n_breakpoints", "n_permutations",
    ] if c in df.columns]
    if cols:
        df = df[cols]
    units = {
        "species": "species code",
        "superfamily": "TE superfamily (EarlGrey)",
        "te_class": "TE class",
        "obs_te_bp_per_window": "mean TE bp per breakpoint window (±5 kb)",
        "bg_te_bp_per_window": "mean TE bp per random window (100 permutations)",
        "log2_fold_change": "log2(obs / bg_mean)",
        "p_value_permutation": "permutation p-value (one-sided)",
        "wilcoxon_p": "Wilcoxon rank-sum p (observed vs 10,000 random windows)",
        "BH_q": "Benjamini-Hochberg FDR-corrected q",
        "n_breakpoints": "count of synteny breakpoints",
        "n_permutations": "number of random window permutations",
    }
    write_sheet(writer, df, "Fig5C_te_breakpoints", units)


def build_fig5d(writer):
    """Fig5D_te_age_breakpoints: TE insertion age vs enrichment at WGD breakpoints."""
    # Primary: paper-like normalised rates
    path = WGD_DIR / "duplication_rates_by_node_normalized.tsv"
    df = safe_read(path)
    if df is None:
        # Fallback: raw per-node rates
        path = WGD_DIR / "duplication_node_rates.tsv"
        df = safe_read(path)
    if df is None:
        path = WGD_DIR / "species_duplication_peak_signals.tsv"
        df = safe_read(path)
    units = {
        "node_id": "OrthoFinder / IQ-TREE node label",
        "node_type": "terminal | internal",
        "n_desc_species": "number of descendant species",
        "n_duplications": "count of duplication events mapped to this node",
        "n_eligible_trees": "orthogroups eligible for duplication at this node",
        "dup_rate_per_eligible_tree": "n_duplications / n_eligible_trees (normalised)",
        "max_species_copies_filter": "max allowed copies per species (expansion filter)",
    }
    write_sheet(writer, df, "Fig5D_te_age_breakpoints", units)


def build_fig6a(writer):
    """Fig6A_htt_chord: validated HTT event counts per clade pair (chord diagram)."""
    path = HTT_DIR / "06_htt_candidates" / "htt_candidates" / "htt_by_clade_pair.tsv"
    df = safe_read(path)
    if df is None:
        path2 = S6_DIR / "results" / "htt_clade_pairs.tsv"
        df = safe_read(path2)
    units = {
        "qclade": "query/donor clade",
        "sclade": "subject/recipient clade",
        "htt_count": "count of HTT candidate pairs (is_htt=True)",
    }
    write_sheet(writer, df, "Fig6A_htt_chord", units)


def build_fig6b(writer):
    """Fig6B_htt_ks: Ks distribution of HTT candidates stacked by TE class."""
    path = HTT_DIR / "06_htt_candidates" / "htt_candidates" / "htt_candidates.tsv"
    df = safe_read(path)
    if df is not None:
        # Only true HTT pairs; keep informative columns
        if "is_htt" in df.columns:
            df = df[df["is_htt"] == True].copy()
        cols = [c for c in ["pair_id", "qclade", "sclade", "te_ks", "ka",
                             "mrca_node", "busco_ks_threshold"] if c in df.columns]
        df = df[cols]
    units = {
        "pair_id": "unique HTT candidate pair identifier",
        "qclade": "donor clade",
        "sclade": "recipient clade",
        "te_ks": "Ks of TE copies between species pair",
        "ka": "Ka (non-synonymous rate)",
        "mrca_node": "MRCA node of the two species",
        "busco_ks_threshold": "BUSCO-based Ks threshold for this MRCA node",
    }
    write_sheet(writer, df, "Fig6B_htt_ks", units)


def build_fig6c(writer):
    """Fig6C_htt_heatmap: HTT candidate count matrix per clade pair."""
    path = HTT_DIR / "06_htt_candidates" / "htt_candidates" / "htt_by_clade_pair.tsv"
    df = safe_read(path)
    if df is not None and "qclade" in df.columns and "sclade" in df.columns:
        pivot = (df.groupby(["qclade", "sclade"])["htt_count"]
                 .sum().unstack(fill_value=0).reset_index())
        pivot = pivot.rename(columns={"qclade": "donor_clade"})
        df = pivot
    units = None  # wide format; column headers are clade codes
    write_sheet(writer, df, "Fig6C_htt_heatmap", units)


def build_fig6d(writer):
    """Fig6D_htt_line_synteny: LINE HTT locus flanking gene conservation."""
    events_dir = HTT_DIR / "09_minimal_htt_events" / "events"
    frames = []
    if events_dir.is_dir():
        for tsv in sorted(events_dir.glob("*_per_node.tsv")):
            edf = safe_read(tsv)
            if edf is not None:
                edf["component"] = tsv.stem.replace("_per_node", "")
                frames.append(edf)
    df = pd.concat(frames, ignore_index=True) if frames else None
    if df is None:
        # Fallback: minimal events summary
        df = safe_read(HTT_DIR / "09_minimal_htt_events" / "summary" / "minimal_events_by_component.tsv")
    units = {
        "component": "HTT community component ID",
        "component_id": "HTT community component ID",
        "n_communities": "count of HTT communities in this component",
        "n_clades": "number of distinct clades sharing this TE",
        "n_distinct_clade_pairs": "unique clade pairs connected",
        "min_events_component": "minimum HTT transfer events (Dollo parsimony)",
    }
    write_sheet(writer, df, "Fig6D_htt_line_synteny", units)


def build_fig6e(writer):
    """Fig6E_htt_te_enrichment: TE class enrichment among HTT candidates (hypergeometric)."""
    path = HTT_DIR / "figures" / "figHTT_te_enrichment.tsv"
    df = safe_read(path)
    units = {
        "te_class": "TE class",
        "obs": "observed count of HTT candidates in this TE class",
        "bg": "background genome bp in this TE class",
        "expected": "expected count based on background frequency",
        "obs_exp": "observed / expected ratio",
        "log2_obs_exp": "log2(observed / expected)",
        "bg_frac": "fraction of background in this TE class",
        "obs_frac": "fraction of HTT candidates in this TE class",
        "residual": "(obs - expected) / sqrt(expected) — Pearson residual",
    }
    write_sheet(writer, df, "Fig6E_htt_te_enrichment", units)


def build_fig6f(writer):
    """Fig6F_htt_full: full HTT community set (all 223 candidates) for Supp Fig S7."""
    path = HTT_DIR / "06_htt_candidates" / "htt_candidates" / "htt_candidates.tsv"
    df = safe_read(path)
    if df is not None and "is_htt" in df.columns:
        df = df[df["is_htt"] == True].copy()
    if df is not None:
        cols = [c for c in ["pair_id", "qclade", "sclade", "te_ks", "mrca_node",
                             "qgenome", "sgenome", "qseqid", "sseqid"] if c in df.columns]
        df = df[cols]
    units = {
        "pair_id": "unique HTT candidate pair ID",
        "qclade": "donor clade",
        "sclade": "recipient clade",
        "te_ks": "Ks of the TE copies between the two species",
        "mrca_node": "MRCA node of the species pair",
        "qgenome": "donor genome species code",
        "sgenome": "recipient genome species code",
    }
    write_sheet(writer, df, "Fig6F_htt_full", units)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

SHEETS = [
    ("Fig1B_pangenome", build_fig1b),
    ("Fig1C_genomic_features", build_fig1c),
    ("Fig1D_synteny", build_fig1d),
    ("Fig2A_te_content", build_fig2a),
    ("Fig2B_pca_composition", build_fig2b),
    ("Fig2C_kimura_ri", build_fig2c),
    ("Fig2D_kimura_summary", build_fig2d),
    ("Fig2E_proximity_scatter", build_fig2e),
    ("Fig2F_proximity_flanks", build_fig2f),
    ("Fig3A_dnmt", build_fig3a),
    ("Fig3B_5mc_family_age", build_fig3b),
    ("Fig3C_gene_dist_5mc", build_fig3c),
    ("Fig3D_adj_expr_meth", build_fig3d),
    ("Fig3E_gene_dist_h3k9", build_fig3e),
    ("Fig3F_adj_expr_h3k9", build_fig3f),
    ("Fig4B_domain_barchart", build_fig4b),
    ("Fig4C_te_class_prop", build_fig4c),
    ("Fig4D_fusion_network", build_fig4d),
    ("Fig4E_chimte", build_fig4e),
    ("Fig4F_flywch_trees", build_fig4f),
    ("Fig4G_coth_trees", build_fig4g),
    ("Fig4H_dam_rt_tree", build_fig4h),
    ("Fig5A_paralogs", build_fig5a),
    ("Fig5B_syntenic_blocks", build_fig5b),
    ("Fig5C_te_breakpoints", build_fig5c),
    ("Fig5D_te_age_breakpoints", build_fig5d),
    ("Fig6A_htt_chord", build_fig6a),
    ("Fig6B_htt_ks", build_fig6b),
    ("Fig6C_htt_heatmap", build_fig6c),
    ("Fig6D_htt_line_synteny", build_fig6d),
    ("Fig6E_htt_te_enrichment", build_fig6e),
    ("Fig6F_htt_full", build_fig6f),
]


def main():
    parser = argparse.ArgumentParser(
        description="Generate EDF-TE source data Excel workbook."
    )
    parser.add_argument(
        "--outdir", default=".",
        help="Directory to write EDF-TE_source_data.xlsx (default: current dir)"
    )
    args = parser.parse_args()

    outpath = Path(args.outdir) / "EDF-TE_source_data.xlsx"
    print(f"Writing source data to: {outpath}")

    with pd.ExcelWriter(outpath, engine="openpyxl") as writer:
        for sheet_name, build_fn in SHEETS:
            print(f"  Building {sheet_name} …")
            try:
                build_fn(writer)
            except Exception as exc:
                print(f"  [ERROR] {sheet_name}: {exc}", file=sys.stderr)
                write_sheet(writer, None, sheet_name)

    print("Done.")


if __name__ == "__main__":
    main()
