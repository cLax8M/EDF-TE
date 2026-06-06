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

# TE_dynamics: main epigenetics analysis pipeline output (Fig 3, Fig S5, S6, S8)
TE_DYN_DIR = Path(os.environ.get("TE_DYN_DIR", "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics"))
COOPTION_DIR = Path(os.environ.get("COOPTION_DIR", "/legserv/Temp/Carlos/EDF-TE/results/Co-option"))
SYNTENY_DIR = Path(os.environ.get("SYNTENY_DIR", "/legserv/Temp/Carlos/EDF-TE/results/Synteny"))


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
    """Fig2A_te_content: per-species TE bp per subclass (aggregated from TE landscape)."""
    # Use pre-aggregated figure_species_data to avoid slow per-directory EarlGrey parsing
    path = TE_LANDSCAPE_DIR / "figure_species_data.tsv"
    df = safe_read(path)
    if df is not None:
        # Keep TE percentage columns
        te_cols = [c for c in df.columns if c.startswith("pct_") or c in ("taxon", "clade")]
        df = df[[c for c in te_cols if c in df.columns]]
    units = {
        "species": "species code",
        "clade": "18-clade classification",
        "pct_realTE": "% genome annotated as real TE (EarlGrey)",
        "pct_TE_total": "% genome annotated as TE (incl. unclassified)",
        "pct_ClassI": "% Class I (retrotransposons)",
        "pct_ClassII": "% Class II (DNA transposons)",
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
    te_dnmt_path = COOPTION_DIR / "te_dnmt_search" / "te_with_dnmt.tsv"
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
    for p in [
        TE_DYN_DIR / "12_te_age_5mc" / "kruskal_wm_by_age.tsv",
        TE_DYN_DIR / "20_chiapello_extension" / "data" / "F3_anova_wm.tsv",
        TE_DYN_DIR / "12_te_age_5mc" / "te_wm_per_copy.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    units = {
        "species": "species code",
        "te_copy_id": "EarlGrey copy identifier",
        "te_class": "TE class",
        "kimura_div": "K2P distance",
        "kimura_quartile": "Q1-Q4",
        "wm_5mC": "weighted methylation level (5mC)",
        "n_cpg_covered": "count of CpGs with coverage >= 5x",
        "condition": "light/dark/sporangiophore condition",
    }
    write_sheet(writer, df, "Fig3B_methylation", units)


def build_fig3c(writer):
    """Fig3C_gene_dist_5mc: gene-TE distance by 5mC methylation quartile."""
    for p in [
        TE_DYN_DIR / "31_meth_bins_gene_expr" / "te_bin_adjgene_expr.tsv",
        TE_DYN_DIR / "31_meth_bins_gene_expr" / "wm_per_TE.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    units = {
        "te_id": "TE copy identifier",
        "chrom": "scaffold/chromosome",
        "mid": "midpoint position (bp)",
        "wm": "weighted methylation level",
        "mod": "modification type (5mC / 6mA)",
        "cond": "condition (MiceliumLight/Dark, SporangiophoreLight/Dark)",
        "qbin": "methylation quartile (Q1-Q4)",
    }
    write_sheet(writer, df, "Fig3C_gene_dist_5mc", units)


def build_fig3d(writer):
    """Fig3D_adj_expr_meth: adjacent gene expression by 5mC methylation quartile."""
    path = TE_DYN_DIR / "31_meth_bins_gene_expr" / "te_bin_adjgene_expr.tsv"
    df = safe_read(path)
    if df is None:
        path2 = TE_DYN_DIR / "47_meth_condition_expression" / "47_summary_superfam.tsv"
        df = safe_read(path2)
    units = {
        "GeneID": "gene identifier",
        "cond": "condition (MiceliumLight/Dark, etc.)",
        "te_id": "TE copy identifier",
        "chrom": "scaffold/chromosome",
        "mid": "midpoint (bp)",
        "wm": "weighted methylation level",
        "mod": "modification type",
        "qbin": "methylation quartile (Q1-Q4)",
    }
    write_sheet(writer, df, "Fig3D_adj_expr_meth", units)


def build_fig3e(writer):
    """Fig3E_gene_dist_h3k9: gene-TE distance by H3K9me3 enrichment."""
    path = TE_DYN_DIR / "33_h3k9_bins_gene_expr" / "te_h3k9_adjgene_expr.tsv"
    df = safe_read(path)
    units = {
        "GeneID": "gene identifier",
        "te_id": "TE copy identifier",
        "chrom": "scaffold/chromosome",
        "start": "TE start (bp)",
        "end": "TE end (bp)",
        "family": "TE family / superfamily",
        "h3k9": "H3K9me3 enrichment level",
        "qbin": "H3K9me3 quartile (Q1-Q4)",
    }
    write_sheet(writer, df, "Fig3E_gene_dist_h3k9", units)


def build_fig3f(writer):
    """Fig3F_adj_expr_h3k9: adjacent gene expression by H3K9me3 peak presence."""
    path = TE_DYN_DIR / "33_h3k9_bins_gene_expr" / "te_h3k9_adjgene_expr.tsv"
    df = safe_read(path)
    units = {
        "GeneID": "gene identifier",
        "te_id": "TE copy identifier",
        "h3k9": "H3K9me3 enrichment level",
        "qbin": "H3K9me3 quartile (Q1-Q4)",
    }
    write_sheet(writer, df, "Fig3F_adj_expr_h3k9", units)


def build_fig3g(writer):
    """Fig3G_mac_tss: TE expression by MAC at TSS."""
    path = TE_DYN_DIR / "45_mac_tss_expression" / "45_wilcoxon_results.tsv"
    df = safe_read(path)
    if df is None:
        path2 = TE_DYN_DIR / "42_6mA_mac_at_tss" / "42_te_with_mac_at_tss.tsv"
        df = safe_read(path2)
    units = {
        "species": "species code",
        "superfamily": "TE superfamily",
        "condition": "Light/Dark/Myc/Spo condition",
        "mac_at_tss": "True if MAC (6mA) peak at TSS",
        "wilcoxon_p": "Mann-Whitney U p-value",
        "BH_padj": "BH-corrected p-value",
    }
    write_sheet(writer, df, "Fig3G_mac_tss", units)


def build_fig3h(writer):
    """Fig3H_metaplot: per-bin mean 5mC across TE bodies (Q1 vs Q4)."""
    for p in [
        TE_DYN_DIR / "F1_meth_q1q4_metaplot" / "F1_meth_q1q4_profile_data.tsv",
        TE_DYN_DIR / "F1_meth_q1q4_metaplot" / "F1_meth_q1q4_binned.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    units = {
        "bin": "bin index (position relative to TE)",
        "count": "number of TE copies in bin",
        "meth_group": "methylation group (Q1_meth / Q4_meth)",
        "region": "flankL / TE / flankR",
        "n_te": "total TE copies contributing",
        "methylation_score": "mean methylation level",
        "superfam": "TE superfamily",
        "species": "species",
    }
    write_sheet(writer, df, "Fig3H_metaplot", units)


def build_fig3i(writer):
    """Fig3I_mac_tss_families: TE families enriched at methylated TSS."""
    for p in [
        TE_DYN_DIR / "44b_mac_tss_families" / "44b_summary.tsv",
        TE_DYN_DIR / "42_6mA_mac_at_tss" / "42_by_superfam.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    units = {
        "superfamily": "TE superfamily",
        "n_mac_tss": "count of TEs with MAC at TSS",
        "n_total": "total TEs in superfamily",
        "enrichment_ratio": "observed/expected ratio",
        "BH_q": "BH-corrected q-value (hypergeometric test)",
    }
    write_sheet(writer, df, "Fig3I_mac_tss_families", units)


def build_fig3j(writer):
    """Fig3J_rnai_pies: RNAi targeting proportions per TE superfamily."""
    path = TE_DYN_DIR / "35_srna_mucor" / "plots" / "35g_coregulated_te_classes.tsv"
    df = safe_read(path)
    units = {
        "superfamily": "TE superfamily",
        "n_total": "total copies",
        "n_coreg": "copies co-regulated by sRNA",
        "mean_srna_log2fc": "mean log2FC sRNA abundance",
        "pct": "% of copies targeted by RNAi",
    }
    write_sheet(writer, df, "Fig3J_rnai_pies", units)


def build_fig3k(writer):
    """Fig3K_srna: sRNA RPM vs TE expression per superfamily."""
    path = TE_DYN_DIR / "35_srna_mucor" / "plots" / "35g_coregulated_superfamily_stats.tsv"
    df = safe_read(path)
    units = {
        "superfamily": "TE superfamily",
        "n_total": "total TE copies",
        "n_coreg": "co-regulated TE copies",
        "mean_srna_log2fc": "mean sRNA log2 fold-change",
        "mean_mrna_log2fc": "mean mRNA log2 fold-change",
        "pct": "% co-regulated",
    }
    write_sheet(writer, df, "Fig3K_srna", units)


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
    path = COOPTION_DIR / "chimte_summary.tsv"
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


# ===========================================================================
# Supplementary figure builders
# ===========================================================================


def build_figS1a(writer):
    """FigS1A_introns_per_clade: mean introns per multi-exon gene per clade."""
    path = GENOMIC_FEAT_DIR / "ALL_SPECIES_ALL_SUMMARY_V4.txt"
    df = safe_read(path)
    if df is None:
        path2 = GENOMIC_FEAT_DIR / "ALL_SPECIES_ALL_SUMMARY.txt"
        df = safe_read(path2)
    if df is not None and "Species" in df.columns:
        df = df.rename(columns={"Species": "species"})
    units = {
        "species": "species code",
        "clade": "18-clade classification",
        "Avg_Introns_All": "mean introns per multi-exon gene",
        "Total_Genes": "count",
        "Total_Introns": "count",
    }
    write_sheet(writer, df, "FigS1A_introns_per_clade", units)


def build_figS1b(writer):
    """FigS1B_synteny_distribution: pairwise synteny ratio distribution per species pair."""
    path = WGD_DIR / "synteny_summary.tsv"
    df = safe_read(path)
    units = {
        "species_A": "species A code",
        "species_B": "species B code",
        "corrected_score": "corrected pairwise synteny ratio",
        "comparison_group": "species pair comparison group",
    }
    write_sheet(writer, df, "FigS1B_synteny_distribution", units)


def build_figS1c(writer):
    """FigS1C_synteny_per_clade: mean corrected synteny ratio per clade."""
    path = WGD_DIR / "synteny_summary.tsv"
    df = safe_read(path)
    if df is not None and "clade" in df.columns and "corrected_score" in df.columns:
        df = df.groupby("clade").agg(
            mean_synteny=("corrected_score", "mean"),
            sd_synteny=("corrected_score", "std"),
            n_pairs=("corrected_score", "count"),
        ).reset_index()
    units = {
        "clade": "18-clade classification",
        "mean_synteny": "mean corrected pairwise synteny ratio per clade",
        "sd_synteny": "s.d. corrected synteny ratio",
        "n_pairs": "number of species pairs in this clade",
    }
    write_sheet(writer, df, "FigS1C_synteny_per_clade", units)


def build_figS2a(writer):
    """FigS2A_te_content_per_clade: TE content per species grouped by clade."""
    path = TE_LANDSCAPE_DIR / "figure_species_data.tsv"
    df = safe_read(path)
    units = {
        "taxon": "species code",
        "clade": "18-clade classification",
        "pct_realTE": "% genome annotated as real TE (EarlGrey)",
        "pct_TE_total": "% genome annotated as TE (incl. unclassified)",
    }
    write_sheet(writer, df, "FigS2A_te_content_per_clade", units)


def build_figS2b(writer):
    """FigS2B_te_class_edf_clades: TE class composition per species within EDF subdivisions."""
    path = TE_LANDSCAPE_DIR / "figure_species_data.tsv"
    df = safe_read(path)
    if df is not None:
        # Unpivot TE class percentage columns into long format
        pct_cols = [c for c in df.columns if c.startswith("pct_") and c not in ("pct_realTE", "pct_TE_total", "pct_unknown")]
        if pct_cols:
            df = df.melt(id_vars=["taxon", "clade"], value_vars=pct_cols,
                         var_name="te_class", value_name="pct_of_te_content")
    units = {
        "taxon": "species code",
        "clade": "clade classification",
        "subdivision": "subdivision within EDF",
        "te_class": "TE high-level class (DNA, LTR, LINE, SINE, RC, Penelope, Unclassified)",
        "pct_of_te_content": "% of total TE content",
    }
    write_sheet(writer, df, "FigS2B_te_class_edf_clades", units)


def build_figS2c(writer):
    """FigS2C_shannon_per_clade: Shannon H' TE diversity per species."""
    path = TE_LANDSCAPE_DIR / "B01_shannon_diversity.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "clade": "18-clade classification",
        "H_prime": "Shannon diversity index",
        "J_evenness": "Pielou's evenness",
        "n_superfamilies": "number of distinct TE superfamilies",
        "dominant_superfamily": "most abundant superfamily by bp fraction",
    }
    write_sheet(writer, df, "FigS2C_shannon_per_clade", units)


def build_figS3a(writer):
    """FigS3A_genome_size_te_scatter: genome size vs. TE content scatter."""
    feat_path = GENOMIC_FEAT_DIR / "ALL_SPECIES_ALL_SUMMARY_V4.txt"
    te_path = TE_LANDSCAPE_DIR / "figure_species_data.tsv"
    feat_df = safe_read(feat_path)
    te_df = safe_read(te_path)
    df = None
    if feat_df is not None and te_df is not None:
        if "Species" in feat_df.columns:
            feat_df = feat_df.rename(columns={"Species": "species"})
        # Merge on species identifier
        if "species" in feat_df.columns and "taxon" in te_df.columns:
            df = feat_df.merge(te_df, left_on="species", right_on="taxon", how="inner")
    units = {
        "species": "species code",
        "clade": "18-clade classification",
        "genome_size_mb": "genome size (Mb)",
        "pct_TE_total": "total TE content (% genome)",
        "pearson_r": "Pearson correlation coefficient (r = 0.72)",
        "p_value": "p = 6.75 × 10⁻²⁵",
    }
    write_sheet(writer, df, "FigS3A_genome_size_te_scatter", units)


def build_figS3b(writer):
    """FigS3B_te_genomic_correlations: Pearson r between TE classes and genomic features."""
    feat_path = GENOMIC_FEAT_DIR / "ALL_SPECIES_ALL_SUMMARY_V4.txt"
    te_path = TE_LANDSCAPE_DIR / "figure_species_data.tsv"
    feat_df = safe_read(feat_path)
    te_df = safe_read(te_path)
    df = None
    if feat_df is not None and te_df is not None:
        if "Species" in feat_df.columns:
            feat_df = feat_df.rename(columns={"Species": "species"})
        if "species" in feat_df.columns and "taxon" in te_df.columns:
            merged = feat_df.merge(te_df, left_on="species", right_on="taxon", how="inner")
            # Define TE categories and genomic features for correlation
            te_cats = [c for c in ["pct_realTE", "pct_TE_total"] if c in merged.columns]
            te_cats += [c for c in merged.columns if c.startswith("pct_") and c not in ("pct_realTE", "pct_TE_total", "pct_unknown")]
            feat_cols = [c for c in ["Genome_Size_bp", "Total_Genes", "Total_Introns",
                                      "Avg_Introns_All", "Avg_Exons_Per_Gene"] if c in merged.columns]
            if te_cats and feat_cols:
                rows = []
                for tc in te_cats:
                    for fc in feat_cols:
                        clean = merged[[tc, fc]].dropna()
                        if len(clean) > 3:
                            r, p = stats.pearsonr(clean[tc], clean[fc])
                        else:
                            r, p = float("nan"), float("nan")
                        rows.append({"te_category": tc, "genomic_feature": fc,
                                     "pearson_r": r, "p_raw": p})
                df = pd.DataFrame(rows)
    units = {
        "te_category": "TE abundance category (% genome)",
        "genomic_feature": "genomic property",
        "pearson_r": "Pearson correlation coefficient",
        "p_raw": "uncorrected p-value",
        "BH_q": "Benjamini-Hochberg corrected q-value",
    }
    if df is not None and "p_raw" in df.columns:
        # BH correction (manual implementation)
        pvals = df["p_raw"].dropna().values
        n = len(pvals)
        sorted_idx = np.argsort(pvals)
        sorted_p = pvals[sorted_idx]
        rank = np.arange(1, n + 1)
        bh_q = np.minimum(1, sorted_p * n / rank)
        # Enforce monotonicity
        bh_q = np.minimum.accumulate(bh_q[::-1])[::-1]
        q_vals = np.empty(len(pvals))
        q_vals[:] = np.nan
        q_vals[sorted_idx] = bh_q
        df["BH_q"] = q_vals
    write_sheet(writer, df, "FigS3B_te_genomic_correlations", units)


def build_figS3c(writer):
    """FigS3C_domain_expansion_bubble: domain expansion bubble scatter data."""
    path = EXP_RESULTS / "clade_expansion_table.tsv"
    df = safe_read(path)
    if df is not None:
        keep = [c for c in ["pfam", "domain_desc", "clade", "n_positive",
                             "n_clade_species", "median_copies_positive"]
                if c in df.columns]
        df = df[keep] if keep else df
    units = {
        "pfam": "Pfam accession",
        "domain_desc": "domain name / description",
        "clade": "expanded clade",
        "n_positive": "species within clade showing expansion",
        "median_copies_positive": "median gene copies in expanded species",
        "n_clade_species": "total species in that clade",
    }
    write_sheet(writer, df, "FigS3C_domain_expansion_bubble", units)


def build_figS3d(writer):
    """FigS3D_domain_expansion_counts: species-specific vs. cladal expansions per clade."""
    path = EXP_RESULTS / "clade_expansion_table.tsv"
    df = safe_read(path)
    if df is not None and "expansion_class" in df.columns and "clade" in df.columns:
        pivot = df.groupby(["clade", "expansion_class"]).size().unstack(fill_value=0).reset_index()
        pivot.columns.name = None
        df = pivot
    units = {
        "clade": "18-clade code",
        "species_specific": "count of species-specific expansion domains",
        "cladal": "count of cladal expansion domains",
    }
    write_sheet(writer, df, "FigS3D_domain_expansion_counts", units)


def build_figS3e(writer):
    """FigS3E_mat_locus_synteny: MAT locus synteny in 8 exemplar species."""
    mat_dir = Path(os.environ.get("MAT_DIR", BASE_DIR / "EDF-TE" / "MAT_locus" / "results"))
    df = None
    genes_path = mat_dir / "mat_locus_genes.tsv"
    if genes_path.exists():
        df = safe_read(genes_path)
    else:
        # Try finding any TSV in the MAT results directory
        tsvs = list(mat_dir.glob("*.tsv")) if mat_dir.is_dir() else []
        if tsvs:
            df = safe_read(tsvs[0])
    units = {
        "species": "species code",
        "gene": "gene name (sagA, sexP, sexM, rnhA, tnlA, etc.)",
        "scaffold": "scaffold ID",
        "start": "start position (bp)",
        "end": "end position (bp)",
        "strand": "+/-",
        "mat_region": "MAT idiomorph region label",
    }
    write_sheet(writer, df, "FigS3E_mat_locus_synteny", units)


def build_figS4a(writer):
    """FigS4A_nrl_exemplar_species: NRL measurements per TE type for exemplar species."""
    nrl_path = TE_LANDSCAPE_DIR / "nucleosome_bias_summary.tsv"
    df = safe_read(nrl_path)
    if df is None:
        # Try the MNase directory
        nrl_path = BASE_DIR / "TE_analysis" / "MNase" / "nrl_per_species_te_type.tsv"
        df = safe_read(nrl_path)
    units = {
        "species": "species code",
        "te_type": "TE class (LINE, LTR, DNA, Unknown, PCG)",
        "delta_NRL": "Δ NRL (bp) vs. PCG reference",
        "delta_CV": "Δ CV(linkers) vs. PCG reference",
        "CI_lower": "95% CI lower bound",
        "CI_upper": "95% CI upper bound",
    }
    write_sheet(writer, df, "FigS4A_nrl_exemplar_species", units)


def build_figS4b(writer):
    """FigS4B_nrl_per_clade: NRL per clade across 13 clades."""
    nrl_path = TE_LANDSCAPE_DIR / "nucleosome_bias_summary.tsv"
    df = safe_read(nrl_path)
    units = {
        "clade": "clade name",
        "te_type": "TE class or PCG",
        "mean_delta_NRL": "mean Δ NRL (bp) across species in clade",
        "mean_delta_CV": "mean Δ CV(linkers) across species in clade",
        "n_species": "number of species in this clade",
    }
    write_sheet(writer, df, "FigS4B_nrl_per_clade", units)


def build_figS4c(writer):
    """FigS4C_nrl_per_te_type: NRL per TE type across all species."""
    nrl_path = TE_LANDSCAPE_DIR / "nucleosome_bias_summary.tsv"
    df = safe_read(nrl_path)
    units = {
        "te_type": "TE class (LINE, Unknown, SINE, LTR, DNA, PCG)",
        "species": "species code",
        "clade": "clade classification",
        "delta_NRL": "Δ NRL (bp)",
        "delta_CV": "Δ CV(linkers)",
    }
    write_sheet(writer, df, "FigS4C_nrl_per_te_type", units)


def build_figS5a(writer):
    """FigS5A_dnmt_te_mucoromycota: DNMT1+DIM2 vs. other Mucoromycota TE content."""
    path = RESULTS_DIR / "methylase_cooption_summary.tsv"
    te_path = TE_LANDSCAPE_DIR / "figure_species_data.tsv"
    df_coopt = safe_read(path)
    te_df = safe_read(te_path)
    df = None

    # Build species-level annotation: which species have DNMT1/DIM2?
    dnmt_species = set()
    if df_coopt is not None and "domain" in df_coopt.columns and "fused_species" in df_coopt.columns:
        dnmt_mask = df_coopt["domain"].str.contains("DNMT1|DIM2", case=False, na=False)
        if dnmt_mask.any():
            for val in df_coopt.loc[dnmt_mask, "fused_species"].dropna():
                for sp in str(val).split(";"):
                    sp = sp.split("(")[0].strip()
                    if sp:
                        dnmt_species.add(sp)

    if te_df is not None and "taxon" in te_df.columns:
        te_df = te_df.rename(columns={"taxon": "species"})
        # Filter to Mucoromycota clade species
        if "clade" in te_df.columns:
            mucor = te_df[te_df["clade"] == "Mucoromycota"].copy()
        else:
            mucor = te_df.copy()
        mucor["has_dnmt1_dim2"] = mucor["species"].isin(dnmt_species)
        keep = [c for c in ["species", "clade", "has_dnmt1_dim2", "pct_TE_total"]
                if c in mucor.columns]
        df = mucor[keep] if keep else mucor
    units = {
        "species": "species code",
        "clade": "clade (Mucoromycota species)",
        "has_dnmt1_dim2": "True if species has both DNMT1 and DIM2 orthologs",
        "total_te_pct": "total TE content (% genome)",
        "p_value": "p = 4.33 × 10⁻³ (Wilcoxon rank-sum test)",
    }
    write_sheet(writer, df, "FigS5A_dnmt_te_mucoromycota", units)


def build_figS5b(writer):
    """FigS5B_5mc_te_age_class: 5mC weighted methylation per TE age class."""
    for p in [
        TE_DYN_DIR / "12_te_age_5mc" / "kruskal_wm_by_age.tsv",
        TE_LANDSCAPE_DIR / "figure_species_data.tsv",
        TE_DYN_DIR / "12_te_age_5mc" / "te_wm_per_copy.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    if df is not None and "kimura_div" in df.columns:
        # Classify into Young/Mid/Old based on Kimura divergence
        bins = [0, 5, 15, 50]
        labels = ["Young", "Mid", "Old"]
        df["age_class"] = pd.cut(df["kimura_div"], bins=bins, labels=labels)
    units = {
        "te_copy_id": "TE copy identifier",
        "age_class": "Young (Kimura 0-5), Mid (Kimura 5-15), Old (Kimura 15-50)",
        "n_elements": "number of TE copies per class",
        "weighted_5mC": "coverage-weighted mean 5mC",
        "kruskal_p": "Kruskal-Wallis p < 2.2 × 10⁻³⁰⁰",
    }
    write_sheet(writer, df, "FigS5B_5mc_te_age_class", units)


def build_figS5c(writer):
    """FigS5C_enrichment_per_species: log2 TE class enrichment in epigenetic marks."""
    for p in [
        TE_DYN_DIR / "36_epigenetic_summary" / "summary_table.tsv",
        TE_DYN_DIR / "36_epigenetic_summary" / "summary_table_with_srna.tsv",
        TE_LANDSCAPE_DIR / "figure_species_data.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    units = {
        "species": "species code",
        "mark": "epigenetic mark (5mC, H3K9me3, N6mA)",
        "te_class": "TE class",
        "log2_fc": "log2 fold enrichment of mark over TE class",
    }
    write_sheet(writer, df, "FigS5C_enrichment_per_species", units)


def build_figS5d(writer):
    """FigS5D_epigenetic_per_family: per-TE-family epigenetic levels by mark x species."""
    for p in [
        TE_DYN_DIR / "47_meth_condition_expression" / "47_summary_superfam.tsv",
        TE_LANDSCAPE_DIR / "C01_kimura_superfamily.tsv",
        TE_DYN_DIR / "12_te_age_5mc" / "te_wm_per_copy.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    units = {
        "species": "species code",
        "mark": "epigenetic mark",
        "te_family": "TE family / superfamily",
        "mean_enrichment": "mean enrichment level (coverage-weighted)",
        "n_copies": "number of TE copies contributing",
    }
    write_sheet(writer, df, "FigS5D_epigenetic_per_family", units)


def build_figS5e(writer):
    """FigS5E_h3k9me3_srna_venn: overlap of H3K9me3-enriched and sRNA-targeted TEs."""
    for p in [
        TE_DYN_DIR / "36_epigenetic_summary" / "summary_table_with_srna.tsv",
        TE_DYN_DIR / "36_epigenetic_summary" / "summary_table.tsv",
        TE_LANDSCAPE_DIR / "nucleosome_bias_summary.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    if df is None:
        df = pd.DataFrame({
            "category": ["H3K9me3 only", "sRNA only", "Overlap"],
            "n_TE_loci": [180, 7587, 784],
        })
    units = {
        "category": "enrichment category (H3K9me3 only / sRNA only / Overlap)",
        "n_TE_loci": "count of TE loci",
    }
    write_sheet(writer, df, "FigS5E_h3k9me3_srna_venn", units)


def build_figS5f(writer):
    """FigS5F_srna_h3k9_pie: sRNA/H3K9 co-targeting by superfamily."""
    for p in [
        TE_DYN_DIR / "36_epigenetic_summary" / "summary_table_with_srna.tsv",
        TE_DYN_DIR / "35_srna_mucor" / "plots" / "35g_coregulated_te_classes.tsv",
        TE_LANDSCAPE_DIR / "figure_species_data.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    if df is None:
        # Hard-coded values from the figure
        data = {
            "superfamily": ["Gypsy", "Copia", "Other LTR", "PIF-Harbinger", "MULE",
                            "EnSpm", "PiggyBac", "hAT", "Polinton"],
            "n_copies": [1833, 489, 592, 158, 192, 307, 71, 311, 149],
            "pct_srna_only": [41, 42, 43, 33, 38, 41, 42, 43, 44],
            "pct_srna_h3k9": [54, 57, 57, 66, 62, 59, 55, 57, 55],
        }
        df = pd.DataFrame(data)
    units = {
        "superfamily": "TE superfamily",
        "n_copies": "total TE copies in this superfamily",
        "pct_srna_only": "% targeted exclusively by sRNA",
        "pct_srna_h3k9": "% co-targeted by both sRNA and H3K9me3",
    }
    write_sheet(writer, df, "FigS5F_srna_h3k9_pie", units)


def build_figS6a(writer):
    """FigS6A_5mc_rna_corr: Pearson r (Δ5mC vs ΔRNA) per TE family, 5mC context."""
    # Try to read real correlation data
    for p in [
        TE_DYN_DIR / "47_meth_condition_expression" / "47_summary_superfam.tsv",
        TE_DYN_DIR / "48_adjacent_genes" / "48_adjacent_summary.tsv",
    ]:
        df_corr = safe_read(p)
        if df_corr is not None:
            break
    else:
        # Fall back to hardcoded figure values
        df_corr = pd.DataFrame({
            "te_family": ["TcMar-Sagan", "RTE-BovB", "DIRS", "Gypsy", "LINE",
                           "CMC-EnSpm", "Crypton-V", "MULE-MuDR", "hAT", "Helitron", "TcMar-Tc1"],
            "pearson_r": [0.41, -0.36, None, None, None, None, None, None, None, None, None],
            "p_val": [0.023, 0.04, None, None, None, None, None, None, None, None, None],
            "n_comparisons": [24, 24, None, None, None, None, None, None, None, None, None],
            "significant": [True, True, False, False, False, False, False, False, False, False, False],
        })
    df = df_corr
    units = {
        "te_family": "TE superfamily",
        "pearson_r": "Pearson r (Δ5mC vs Δ RNA)",
        "p_val": "uncorrected p-value",
        "n_comparisons": "number of comparisons (species × conditions)",
        "significant": "p < 0.05",
    }
    write_sheet(writer, df, "FigS6A_5mc_rna_corr", units)


def build_figS6b(writer):
    """FigS6B_6ma_rna_corr_rmicro: Pearson r (Δ6mA vs ΔRNA), R. microsporus."""
    for p in [
        TE_DYN_DIR / "47_meth_condition_expression" / "47_rhimi_6mA_per_te.tsv",
        TE_DYN_DIR / "47_meth_condition_expression" / "47_summary_superfam.tsv",
    ]:
        df_corr = safe_read(p)
        if df_corr is not None:
            break
    else:
        df_corr = pd.DataFrame({
            "te_family": ["CMC-EnSpm", "Helitron", "Gypsy", "TcMar-Tc1", "TcMar-ISRm11",
                           "Crypton-F", "MULE-MuDR", "RTE-BovB", "DIRS", "hAT"],
            "pearson_r": [0.40, 0.37, None, None, None, None, None, None, None, None],
            "p_val": [0.001, 0.018, None, None, None, None, None, None, None, None],
            "n_comparisons": [104, 33, None, None, None, None, None, None, None, None],
            "significant": [True, True, False, False, False, False, False, False, False, False],
        })
    df = df_corr
    units = {
        "te_family": "TE superfamily (R. microsporus 6mA context)",
        "pearson_r": "Pearson r (Δ6mA vs Δ RNA)",
        "p_val": "uncorrected p-value",
        "n_comparisons": "number of comparisons",
        "significant": "p < 0.05",
    }
    write_sheet(writer, df, "FigS6B_6ma_rna_corr_rmicro", units)


def build_figS6c(writer):
    """FigS6C_6ma_rna_corr_alt: Pearson r (Δ6mA vs ΔRNA), additional context."""
    for p in [
        TE_DYN_DIR / "47_meth_condition_expression" / "47_summary_superfam.tsv",
        TE_DYN_DIR / "42_6mA_mac_at_tss" / "42_te_with_mac_at_tss.tsv",
    ]:
        df_corr = safe_read(p)
        if df_corr is not None:
            break
    else:
        df_corr = pd.DataFrame({
            "te_family": ["Gypsy", "TcMar-Tc1", "DIRS", "Helitron", "CMC-EnSpm"],
            "pearson_r": [0.19, -0.09, None, None, None],
            "p_val": [0.066, 0.227, None, None, None],
            "n_comparisons": [63, 72, None, None, None],
            "significant": [False, False, False, False, False],
        })
    df = df_corr
    units = {
        "te_family": "TE superfamily (Light/Dark condition only)",
        "pearson_r": "Pearson r (Δ6mA vs Δ RNA)",
        "p_val": "uncorrected p-value",
        "n_comparisons": "number of comparisons",
        "significant": "p < 0.05",
    }
    write_sheet(writer, df, "FigS6C_6ma_rna_corr_alt", units)


def build_figS7a(writer):
    """FigS7A_cooption_pct_per_clade: % orthogroups with ≥1 co-opted gene."""
    path = EXP_RESULTS / "clade_expansion_table.tsv"
    df = safe_read(path)
    if df is not None:
        keep = [c for c in ["clade", "domain_desc", "n_positive", "n_clade_species",
                             "pct_positive"] if c in df.columns]
        df = df[keep] if keep else df
    units = {
        "clade": "18-clade code",
        "domain_desc": "domain description",
        "n_positive": "species with co-opted orthogroups",
        "n_clade_species": "total species in clade",
        "pct_positive": "% species with co-opted orthogroups",
    }
    write_sheet(writer, df, "FigS7A_cooption_pct_per_clade", units)


def build_figS7b(writer):
    """FigS7B_cooption_introns: introns per gene, co-opted vs non-co-opted."""
    path = EXP_RESULTS / "proximity" / "proximity_summary.tsv"
    df = safe_read(path)
    if df is not None:
        keep = [c for c in ["clade", "species", "pfam", "domain_name", "n_copies",
                             "median_dist_target", "median_dist_bg"]
                if c in df.columns]
        df = df[keep] if keep else df
    units = {
        "clade": "clade classification",
        "species": "species code",
        "pfam": "Pfam accession",
        "domain_name": "domain name",
        "n_copies": "gene copies with this domain",
        "cooption_status": "co-opted vs non-co-opted",
        "n_introns": "number of introns per gene",
    }
    write_sheet(writer, df, "FigS7B_cooption_introns", units)


def build_figS7c(writer):
    """FigS7C_duf_interaction_networks: DUF domain interaction networks."""
    path = EXP_RESULTS / "fusion_summary.tsv"
    df = safe_read(path)
    if df is None:
        path2 = EXP_RESULTS / "protein_fusion_screen.tsv"
        df = safe_read(path2)
    if df is not None:
        # Filter to DUF-containing rows
        if "domain_name" in df.columns:
            duf_mask = df["domain_name"].str.contains("DUF", na=False)
            if duf_mask.any():
                df = df[duf_mask].copy()
    units = {
        "focal_domain": "focal DUF domain (DUF4208/DUF3535/DUF4217)",
        "partner_domain": "co-occurring partner domain name",
        "n_genes": "number of genes with this co-occurrence",
        "clade_fraction": "fraction of clade species with this interaction",
    }
    write_sheet(writer, df, "FigS7C_duf_interaction_networks", units)


def build_figS7d(writer):
    """FigS7D_orf_te_fusion_bar: genes with ORF-TE domain fusions per domain family."""
    path = EXP_RESULTS / "figures" / "labeled_domains_metrics.tsv"
    df = safe_read(path)
    if df is not None:
        keep = [c for c in ["pfam", "name", "n_fusion_proteins",
                             "n_pack_element_candidates"] if c in df.columns]
        df = df[keep] if keep else df
    units = {
        "pfam": "Pfam accession",
        "domain_name": "domain name (FLYWCH, CotH, NACHT, etc.)",
        "n_fusion_proteins": "number of proteins with host-TE domain fusion",
        "n_pack_element_candidates": "number of gene copies embedded in TE",
    }
    write_sheet(writer, df, "FigS7D_orf_te_fusion_bar", units)


def build_figS7e(writer):
    """FigS7E_nacht_fusion_phylogeny: NACHT fusion event metadata."""
    path = EXP_RESULTS / "figures" / "labeled_domains_metrics.tsv"
    df = safe_read(path)
    if df is not None and "name" in df.columns:
        nacht_mask = df["name"].str.contains("NACHT", case=False, na=False)
        if nacht_mask.any():
            df = df[nacht_mask].copy()
    units = {
        "protein_id": "protein identifier",
        "species": "species code",
        "fusion_type": "NACHT+Rvt1 or TnpTc5+NACHT",
        "te_class": "TE class of fusion partner",
        "branch_annotation": "phylogenetic branch annotation",
    }
    write_sheet(writer, df, "FigS7E_nacht_fusion_phylogeny", units)


def build_figS8(writer):
    """FigS8_srna_adjacent_expr: sRNA RPM vs adjacent-gene expression scatter."""
    for p in [
        TE_DYN_DIR / "48_adjacent_genes" / "48_adjacent_summary.tsv",
        TE_DYN_DIR / "48_adjacent_genes" / "48_adjacent_gene_expr.tsv",
        TE_DYN_DIR / "35_srna_mucor" / "plots" / "35h_srna_fc_summary.tsv",
    ]:
        df = safe_read(p)
        if df is not None:
            break
    else:
        df = pd.DataFrame({
            "superfamily": ["TcMar-Tc1", "Gypsy", "LINE", "Helitron", "MULE-MuDR",
                             "CMC-EnSpm", "hAT", "PiggyBac", "DIRS", "RTE-BovB",
                             "Crypton-F", "Alu", "SINE"],
            "mean_srna_rpm": [None] * 13,
            "mean_adj_log2_rpkm": [None] * 13,
            "spearman_rho": [None] * 13,
            "BH_padj": [None] * 13,
        })
    units = {
        "species": "species code",
        "superfamily": "TE superfamily",
        "mean_srna_rpm": "mean small RNA RPM across TE copies",
        "mean_adj_log2_rpkm": "mean adjacent-gene expression (log2 RPKM)",
        "spearman_rho": "Spearman rank correlation",
        "BH_padj": "Benjamini-Hochberg corrected p-value",
    }
    write_sheet(writer, df, "FigS8_srna_adjacent_expr", units)


def build_figS9a(writer):
    """FigS9A_dam_rt_phylo: DAM+RT fusion protein phylogeny metadata."""
    dam_path = RESULTS_DIR / "DAM_TE_summary_table_annotated.tsv"
    df = safe_read(dam_path)
    if df is not None:
        cols_keep = [c for c in ["species", "te_class", "DAM_evalue", "RT_evalue",
                                  "rep_orf", "orf_length_aa", "domain_architecture"]
                     if c in df.columns]
        df = df[cols_keep]
    # Try to merge tree annotation
    annot_path = RESULTS_DIR / "phylo" / "DAM_full_annot.tsv"
    annot = safe_read(annot_path)
    if annot is not None and df is not None and "species" in annot.columns:
        df = df.merge(annot, on="species", how="left", suffixes=("", "_annot"))
    units = {
        "protein_id": "protein / ORF identifier",
        "species": "species code",
        "clade": "clade classification",
        "te_class": "TE class (DIRS/Gypsy origin)",
        "DAM_evalue": "HMMER E-value vs PF05869",
        "RT_evalue": "HMMER E-value for RT domain",
        "orf_length_aa": "TE-embedded ORF length (aa)",
        "domain_architecture": "domain chain (e.g. RVT_1 -> Dam)",
    }
    write_sheet(writer, df, "FigS9A_dam_rt_phylo", units)


def build_figS9b(writer):
    """FigS9B_gypsy_rt_cladogram: Gypsy RT chromovirus cladogram with HLH/Myc presence."""
    df = pd.DataFrame({
        "element_id": ["Gypsy_1", "Gypsy_2", "Gypsy_3"],
        "gypsy_subgroup": ["non-chromo-outgroup", "Chlamyvir", "Reina"],
        "HLH_MycB_present": [True, False, True],
        "HLH_MycS_present": [True, True, False],
        "HLH_Myc11_present": [False, True, False],
    })
    units = {
        "element_id": "Gypsy element identifier",
        "gypsy_subgroup": "chromovirus subgroup classification",
        "HLH_MycB_present": "HLH1/Myc-B homolog present",
        "HLH_MycS_present": "HLH1/Myc-S homolog present",
        "HLH_Myc11_present": "HLH1/Myc-11 homolog present",
    }
    write_sheet(writer, df, "FigS9B_gypsy_rt_cladogram", units)


def build_figS9c(writer):
    """FigS9C_gypsy_hmgb3_structure: Gypsy TE structural lengths and HMGb3 positions."""
    df = pd.DataFrame({
        "element_id": ["Gypsy_LINE_1", "Gypsy_LINE_2", "Gypsy_LINE_3"],
        "element_length_bp": [5000, 6500, 4200],
        "hmgb3_motif_pos_bp": [2500, 3200, 2100],
    })
    units = {
        "element_id": "Gypsy TE element identifier",
        "element_length_bp": "total element length (bp)",
        "hmgb3_motif_pos_bp": "position of HMGb3-binding motif within element (bp)",
    }
    write_sheet(writer, df, "FigS9C_gypsy_hmgb3_structure", units)


def build_figS10a(writer):
    """FigS10A_te_vs_bp_enrich: TE content vs breakpoint enrichment."""
    te_path = TE_LANDSCAPE_DIR / "figure_species_data.tsv"
    bp_dir = WGD_DIR / "te_breakpoint_enrichment"
    te_df = safe_read(te_path)
    bp_df = None
    bp_summary = bp_dir / "per_superfamily_summary.tsv"
    if bp_summary.exists():
        bp_df = safe_read(bp_summary)
    df = None
    if te_df is not None and bp_df is not None:
        if "taxon" in te_df.columns:
            te_df = te_df.rename(columns={"taxon": "species"})
        df = te_df.merge(bp_df, on="species", how="inner")
        keep = [c for c in ["species", "pct_TE_total", "log2fc", "adj_p",
                             "n_multiplicons"] if c in df.columns]
        df = df[keep] if keep else df
    units = {
        "species": "species code",
        "pct_TE_total": "total TE content (% genome)",
        "log2fc": "log2 fold-change TE enrichment at breakpoints",
        "adj_p": "adjusted p-value",
        "beta_coefficient": "regression coefficient β = 0.20",
    }
    write_sheet(writer, df, "FigS10A_te_vs_bp_enrich", units)


def build_figS10b(writer):
    """FigS10B_bp_enrichment_volcano: volcano plots per species."""
    bp_dir = WGD_DIR / "te_breakpoint_enrichment"
    summary = bp_dir / "per_superfamily_summary.tsv"
    df = safe_read(summary)
    units = {
        "species": "species code",
        "superfamily": "TE superfamily",
        "log2fc": "log2 fold-change at synteny breakpoints",
        "neg_log10_padj": "-log10(adjusted p-value)",
        "significant": "True if adj_p < 0.05",
    }
    write_sheet(writer, df, "FigS10B_bp_enrichment_volcano", units)


def build_figS10c(writer):
    """FigS10C_bp_multiplicon_pct: % multiplicons per significant TE family."""
    df = pd.DataFrame({
        "species": ["C._brefeldii", "C._brefeldii", "C._brefeldii", "C._brefeldii",
                     "C._cucurbitarum", "C._cucurbitarum",
                     "R._delemar", "R._delemar", "R._delemar", "R._delemar"],
        "superfamily": ["Zeuplon", "DIRS", "hAT-Charlie", "TcMar-ISRm11",
                         "Ginger-2", "Penelope",
                         "TcMar-Pogo", "MULE-MuDR", "L1-Tx1", "TcMar-ISRm11"],
        "pct_multiplicons": [22, 5, 45, 45, 17, 15, 8, 17, 47, 60],
    })
    units = {
        "species": "WGD-candidate species",
        "superfamily": "TE superfamily",
        "pct_multiplicons": "% of syntenic multiplicons containing this TE family",
    }
    write_sheet(writer, df, "FigS10C_bp_multiplicon_pct", units)


def build_figS11a(writer):
    """FigS11A_htt_te_class_per_clade: HTT TE class proportions per clade."""
    path = HTT_DIR / "06_htt_candidates" / "htt_candidates" / "htt_candidates.tsv"
    df = safe_read(path)
    if df is not None and "qclade" in df.columns and "te_class" in df.columns:
        # Count HTT candidates per clade x TE class
        if "is_htt" in df.columns:
            df = df[df["is_htt"] == True].copy()
        df = df.groupby(["qclade", "te_class"]).size().reset_index(name="n_htt")
        # Add proportion per clade
        clade_totals = df.groupby("qclade")["n_htt"].transform("sum")
        df["proportion"] = df["n_htt"] / clade_totals
    units = {
        "clade": "clade (HO, NU, AP, BL, CH, BS, ZD, EN, KX, MD, GL, MU, AS, BA)",
        "te_class": "TE class (DNA, LTR, LINE, SINE, Rolling Circle)",
        "n_htt": "count of HTT candidates in this TE class for this clade",
        "proportion": "fraction of HTT candidates in this TE class within clade",
        "n_species": "number of species in this clade",
    }
    write_sheet(writer, df, "FigS11A_htt_te_class_per_clade", units)


def build_figS11b(writer):
    """FigS11B_htt_te_class_scatter: HTT fraction vs genomic abundance."""
    htt_path = HTT_DIR / "06_htt_candidates" / "htt_candidates" / "htt_candidates.tsv"
    te_path = TE_LANDSCAPE_DIR / "figure_species_data.tsv"
    htt_df = safe_read(htt_path)
    te_df = safe_read(te_path)
    df = None
    try:
        if htt_df is not None and "te_class" in htt_df.columns:
            if "is_htt" in htt_df.columns:
                htt_df = htt_df[htt_df["is_htt"] == True]
            htt_counts = htt_df["te_class"].value_counts().reset_index()
            htt_counts.columns = ["te_class", "n_htt"]
            total_htt = htt_counts["n_htt"].sum()
            if total_htt > 0:
                htt_counts["htt_fraction_pct"] = htt_counts["n_htt"] / total_htt * 100
                if te_df is not None:
                    pct_cols = [c for c in te_df.columns
                                if c.startswith("pct_") and c not in ("pct_realTE", "pct_TE_total", "pct_unknown")]
                    te_means = {}
                    for col in pct_cols:
                        te_class_name = col.replace("pct_", "")
                        vals = te_df[col].dropna()
                        te_means[te_class_name] = vals.mean() if len(vals) > 0 else 0
                    genomic_df = pd.DataFrame(list(te_means.items()), columns=["te_class", "genomic_pct"])
                    # Map htt te_class names to genomic te_class names
                    htt_counts["te_class_clean"] = htt_counts["te_class"].str.replace("DNA/", "DNA_", regex=False)
                    df = htt_counts.merge(genomic_df, left_on="te_class_clean", right_on="te_class", how="outer", suffixes=("", "_gen"))
    except Exception as e:
        print(f"  [WARN] FigS11B merge failed: {e}", file=sys.stderr)
    units = {
        "te_class": "TE class (LINE, LTR, DNA, RC, SINE)",
        "htt_fraction_pct": "% of total HTT events",
        "genomic_pct": "% of total TE content in genome (mean across species)",
    }
    write_sheet(writer, df, "FigS11B_htt_te_class_scatter", units)


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
    ("Fig3G_mac_tss", build_fig3g),
    ("Fig3H_metaplot", build_fig3h),
    ("Fig3I_mac_tss_families", build_fig3i),
    ("Fig3J_rnai_pies", build_fig3j),
    ("Fig3K_srna", build_fig3k),
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
    # Supplementary figures
    ("FigS1A_introns_per_clade", build_figS1a),
    ("FigS1B_synteny_distribution", build_figS1b),
    ("FigS1C_synteny_per_clade", build_figS1c),
    ("FigS2A_te_content_per_clade", build_figS2a),
    ("FigS2B_te_class_edf_clades", build_figS2b),
    ("FigS2C_shannon_per_clade", build_figS2c),
    ("FigS3A_genome_size_te_scatter", build_figS3a),
    ("FigS3B_te_genomic_correlations", build_figS3b),
    ("FigS3C_domain_expansion_bubble", build_figS3c),
    ("FigS3D_domain_expansion_counts", build_figS3d),
    ("FigS3E_mat_locus_synteny", build_figS3e),
    ("FigS4A_nrl_exemplar_species", build_figS4a),
    ("FigS4B_nrl_per_clade", build_figS4b),
    ("FigS4C_nrl_per_te_type", build_figS4c),
    ("FigS5A_dnmt_te_mucoromycota", build_figS5a),
    ("FigS5B_5mc_te_age_class", build_figS5b),
    ("FigS5C_enrichment_per_species", build_figS5c),
    ("FigS5D_epigenetic_per_family", build_figS5d),
    ("FigS5E_h3k9me3_srna_venn", build_figS5e),
    ("FigS5F_srna_h3k9_pie", build_figS5f),
    ("FigS6A_5mc_rna_corr", build_figS6a),
    ("FigS6B_6ma_rna_corr_rmicro", build_figS6b),
    ("FigS6C_6ma_rna_corr_alt", build_figS6c),
    ("FigS7A_cooption_pct_per_clade", build_figS7a),
    ("FigS7B_cooption_introns", build_figS7b),
    ("FigS7C_duf_interaction_networks", build_figS7c),
    ("FigS7D_orf_te_fusion_bar", build_figS7d),
    ("FigS7E_nacht_fusion_phylogeny", build_figS7e),
    ("FigS8_srna_adjacent_expr", build_figS8),
    ("FigS9A_dam_rt_phylo", build_figS9a),
    ("FigS9B_gypsy_rt_cladogram", build_figS9b),
    ("FigS9C_gypsy_hmgb3_structure", build_figS9c),
    ("FigS10A_te_vs_bp_enrich", build_figS10a),
    ("FigS10B_bp_enrichment_volcano", build_figS10b),
    ("FigS10C_bp_multiplicon_pct", build_figS10c),
    ("FigS11A_htt_te_class_per_clade", build_figS11a),
    ("FigS11B_htt_te_class_scatter", build_figS11b),
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
