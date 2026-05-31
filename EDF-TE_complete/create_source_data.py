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
TE_ANNOT_DIR = Path(os.environ.get("TE_ANNOT_DIR", BASE_DIR / "TE_analysis" / "EarlGrey"))
S2_DIR = Path(os.environ.get("S2_DIR", BASE_DIR / "EDF-TE" / "EDF-TE_complete" / "section2_te_landscape"))
S3_DIR = Path(os.environ.get("S3_DIR", BASE_DIR / "EDF-TE" / "EDF-TE_complete" / "section3_te_epigenetics"))
S4_DIR = Path(os.environ.get("S4_DIR", BASE_DIR / "EDF-TE" / "EDF-TE_complete" / "section4_te_cooption"))
S6_DIR = Path(os.environ.get("S6_DIR", BASE_DIR / "EDF-TE" / "EDF-TE_complete" / "section6_htt"))


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
    """Fig1B_genomic_features: per-species 6 metrics + z-scores."""
    path = BASE_DIR / "EDF-TE" / "Section2_TElandscape" / "genomic_features.tsv"
    df = safe_read(path)
    if df is None:
        # Try alternate location
        alt = BASE_DIR / "EDF-TE" / "GenomicFeatures" / "genomic_features.tsv"
        df = safe_read(alt)
    if df is None:
        write_sheet(writer, None, "Fig1B_genomic_features")
        return
    metrics = [c for c in df.columns if c not in ("species", "clade")]
    for m in metrics:
        if pd.api.types.is_numeric_dtype(df[m]):
            df[f"{m}_zscore"] = stats.zscore(df[m].fillna(df[m].median()), ddof=1)
    units = {
        "species": "species code",
        "clade": "18-clade classification",
        "genome_size_mb": "Mb",
        "gc_pct": "%",
        "gene_count": "count",
        "mean_gene_length_bp": "bp",
        "intron_density": "introns/gene",
        "n50_mb": "Mb",
    }
    write_sheet(writer, df, "Fig1B_genomic_features", units)


def build_fig1c(writer):
    """Fig1C_pangenome: accumulation curve mean ± SD."""
    path = BASE_DIR / "EDF-TE" / "PanGenome" / "pangenome_accumulation.tsv"
    df = safe_read(path)
    if df is None:
        alt = BASE_DIR / "PanGenome" / "pangenome_accumulation.tsv"
        df = safe_read(alt)
    units = {
        "n_species_sampled": "count",
        "core_mean": "gene count",
        "core_sd": "gene count",
        "shell_mean": "gene count",
        "shell_sd": "gene count",
        "unique_mean": "gene count",
        "unique_sd": "gene count",
    }
    write_sheet(writer, df, "Fig1C_pangenome", units)


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
    """Fig2A_te_content: per-species TE bp and % per class."""
    path = BASE_DIR / "TE_analysis" / "EarlGrey" / "te_content_summary.tsv"
    df = safe_read(path)
    if df is None:
        # Try to aggregate from expansion results
        alt = EXP_RESULTS / "expansion_stats.tsv"
        df = safe_read(alt)
    units = {
        "species": "species code",
        "clade": "18-clade",
        "te_class": "TE class (Wicker 2007)",
        "bp_total": "bp",
        "pct_genome": "%",
        "shannon_h": "bits",
    }
    write_sheet(writer, df, "Fig2A_te_content", units)


def build_fig2b(writer):
    """Fig2B_pcoa: PCoA coordinates per species."""
    path = BASE_DIR / "TE_analysis" / "EarlGrey" / "pcoa_coordinates.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "clade": "18-clade",
        "pcoa_bc_1": "arbitrary units",
        "pcoa_bc_2": "arbitrary units",
        "pcoa_ait_1": "arbitrary units",
        "pcoa_ait_2": "arbitrary units",
        "pcoa_jac_1": "arbitrary units",
        "pcoa_jac_2": "arbitrary units",
        "shannon_h": "bits",
    }
    write_sheet(writer, df, "Fig2B_pcoa", units)


def build_fig2c(writer):
    """Fig2C_kimura: Kimura histogram bin values per species."""
    path = BASE_DIR / "TE_analysis" / "EarlGrey" / "kimura_bins.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "clade": "18-clade",
        "te_class": "TE class",
        "kimura_bin_midpoint": "K2P distance (dimensionless)",
        "bp_in_bin": "bp",
        "recency_index": "fraction [0,1]",
    }
    write_sheet(writer, df, "Fig2C_kimura", units)


def build_fig2d(writer):
    """Fig2D_shannon: Shannon H per species."""
    path = BASE_DIR / "TE_analysis" / "EarlGrey" / "shannon_diversity.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "clade": "18-clade",
        "shannon_h": "bits",
        "clade_order": "integer (phylogenetic order of clade)",
    }
    write_sheet(writer, df, "Fig2D_shannon", units)


def build_fig3a(writer):
    """Fig3A_dnmt: DNMT presence/absence per species."""
    path = RESULTS_DIR / "methylase_cooption_summary.tsv"
    df = safe_read(path)
    units = {
        "species": "species code",
        "clade": "18-clade",
        "dnmt_present": "0/1 (IPR005105+IPR022686 both present)",
        "orthogroup_copynum": "count",
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


def build_fig4b(writer):
    """Fig4B_cooption_heatmap: fraction co-opted per clade × domain."""
    path = EXP_RESULTS / "clade_expansion_table.tsv"
    df = safe_read(path)
    units = {
        "clade": "18-clade code",
        "te_domain_ipr": "InterPro accession",
        "te_domain_name": "domain name",
        "frac_species_with_cooption": "fraction [0,1]",
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
    """Fig4D_tf_te: TF × TE class enrichment matrix."""
    path = BASE_DIR / "TF_TE_analysis" / "tf_te_enrichment_matrix.tsv"
    df = safe_read(path)
    units = {
        "tf_motif": "JASPAR motif ID / TF name",
        "te_class": "TE class",
        "neg_log10_fdr": "−log10(BH-FDR p-value)",
        "fisher_or": "Fisher exact odds ratio",
    }
    write_sheet(writer, df, "Fig4D_tf_te", units)


def build_fig4e(writer):
    """Fig4E_fusions: fusion gene counts per clade × fusion class."""
    path = EXP_RESULTS / "fusion_summary.tsv"
    df = safe_read(path)
    units = {
        "clade": "18-clade",
        "fusion_class": "FLYWCH-MULE | CotH-DDE | NACHT-RT | PARP-Chromo",
        "frac_species_with_fusion": "fraction [0,1]",
    }
    write_sheet(writer, df, "Fig4E_fusions", units)


def build_fig4f(writer):
    """Fig4F_dam_amt1: DAM-TE × AMT1 contingency + per-species annotation."""
    dam_path = RESULTS_DIR / "DAM_TE_summary_table_annotated.tsv"
    amt_path = RESULTS_DIR / "ipr_AMT1_presence.tsv"
    dam = safe_read(dam_path)
    amt = safe_read(amt_path)
    if dam is None or amt is None:
        write_sheet(writer, None, "Fig4F_dam_amt1")
        return
    # Merge on species
    df = dam.merge(amt, on="species", how="outer")
    cols_keep = ["species", "clade", "AMT1_status", "te_class",
                 "DAM_evalue", "domain_architecture"]
    # Rename for clarity
    df = df.rename(columns={"AMT1_status": "amt1_present"})
    # Keep only available columns
    df = df[[c for c in cols_keep if c in df.columns]]
    # Create numeric dam_te_present column if not present
    if "dam_te_present" not in df.columns:
        # Any row in the DAM table implies dam_te_present = 1
        all_species_path = RESULTS_DIR.parent / "species_list.txt"
        if all_species_path.exists():
            all_sp = pd.read_csv(all_species_path, header=None, names=["species"])
            df_full = all_sp.merge(df, on="species", how="left")
            df_full["dam_te_present"] = df_full["te_class"].notna().astype(int)
            df = df_full
    units = {
        "species": "species code",
        "clade": "18-clade",
        "dam_te_present": "0/1",
        "amt1_present": "present | absent",
        "te_class": "TE class of DAM-embedding element (DIRS or Gypsy)",
        "DAM_evalue": "HMMER E-value vs PF05869",
        "domain_architecture": "InterProScan domain string",
    }
    write_sheet(writer, df, "Fig4F_dam_amt1", units)


def build_fig5a(writer):
    """Fig5A_paralogs: colinear paralog count per Mucoromycota species."""
    path = WGD_DIR / "wgd_candidates_integrated.tsv"
    df = safe_read(path)
    if df is None:
        syn = safe_read(WGD_DIR / "synteny_summary.tsv")
        wgd = safe_read(WGD_DIR / "wgd_classification.tsv")
        if syn is not None and wgd is not None:
            df = syn.merge(wgd, on="species_id", how="outer")
    if df is not None:
        # Select Mucoromycota (MU_ prefix) and ≤200 scaffolds (scaffold_n column if present)
        if "species_id" in df.columns:
            df = df.rename(columns={"species_id": "species"})
        if "species" in df.columns:
            df = df[df["species"].str.startswith("MU_", na=False)].copy()
        cols = ["species", "n_anchorpoints", "n_multiplicons", "n_segments",
                "classification", "confidence", "has_dup_signal", "has_synteny"]
        df = df[[c for c in cols if c in df.columns]]
        df = df.rename(columns={"n_anchorpoints": "n_colinear_paralog_pairs",
                                 "classification": "wgd_classification"})
    units = {
        "species": "Mucoromycota species code (MU_ prefix)",
        "n_colinear_paralog_pairs": "count of syntenic anchor pairs (colinear paralogs)",
        "n_multiplicons": "count of syntenic multiplicons",
        "n_segments": "count of collinear segments",
        "wgd_classification": "strong | moderate | no_WGD (from GMM + synteny)",
        "has_dup_signal": "bool (GMM Ks peak detected)",
        "has_synteny": "bool (synteny evidence present)",
    }
    write_sheet(writer, df, "Fig5A_paralogs", units)


def build_fig5b(writer):
    """Fig5B_breakpoints: TE density at breakpoints vs random per species."""
    path = WGD_DIR / "enrichment_summary.tsv"
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
    write_sheet(writer, df, "Fig5B_breakpoints", units)


def build_fig5c(writer):
    """Fig5C_te_class_bp: TE class at breakpoints vs background per species."""
    path = WGD_DIR / "per_superfamily_summary.tsv"
    df = safe_read(path)
    if df is None:
        write_sheet(writer, None, "Fig5C_te_class_bp")
        return
    # per_superfamily_summary is wide: species × TE superfamily
    # Melt to long format
    id_vars = [c for c in df.columns if c in ("species", "species_id")]
    if not id_vars:
        write_sheet(writer, df, "Fig5C_te_class_bp")
        return
    id_col = id_vars[0]
    melted = df.melt(id_vars=[id_col], var_name="te_class",
                     value_name="bp_at_breakpoints")
    melted = melted.rename(columns={id_col: "species"})
    units = {
        "species": "species code",
        "te_class": "TE class / superfamily",
        "bp_at_breakpoints": "bp of TE class at breakpoint windows",
    }
    write_sheet(writer, melted, "Fig5C_te_class_bp", units)


def build_fig5d(writer):
    """Fig5D_dup_rate: gene duplication rate per branch."""
    path = WGD_DIR / "duplication_node_rates.tsv"
    df = safe_read(path)
    units = {
        "node_id": "IQ-TREE node label",
        "node_type": "internal | leaf",
        "n_duplications": "count of duplication events on branch",
        "dup_rate_per_desc_leaf": "duplications per descendant leaf (normalised)",
        "descendant_leaves": "count",
        "tree_depth_scaled": "scaled branch depth [0,1]",
    }
    write_sheet(writer, df, "Fig5D_dup_rate", units)


def build_fig6a(writer):
    """Fig6A_htt_chord: HTT candidate counts per clade pair."""
    path = S6_DIR / "results" / "htt_clade_pairs.tsv"
    df = safe_read(path)
    if df is None:
        path2 = BASE_DIR / "EDF-TE" / "HTT" / "htt_clade_pairs.tsv"
        df = safe_read(path2)
    units = {
        "q_clade": "query clade (donor)",
        "s_clade": "subject clade (recipient)",
        "n_candidates": "count of HTT candidates",
        "dominant_te_class": "most frequent TE class in this clade pair",
    }
    write_sheet(writer, df, "Fig6A_htt_chord", units)


def build_fig6b(writer):
    """Fig6B_htt_class: HTT candidates per clade pair × TE class."""
    path = S6_DIR / "results" / "htt_candidates_annotated.tsv"
    df = safe_read(path)
    if df is None:
        path2 = BASE_DIR / "EDF-TE" / "HTT" / "htt_candidates_annotated.tsv"
        df = safe_read(path2)
    if df is not None and "q_clade" in df.columns and "te_class" in df.columns:
        agg = (df.groupby(["q_clade", "s_clade", "te_class"])
               .size().reset_index(name="n_candidates"))
        df = agg
    units = {
        "q_clade": "donor clade",
        "s_clade": "recipient clade",
        "te_class": "TE class",
        "n_candidates": "count",
    }
    write_sheet(writer, df, "Fig6B_htt_class", units)


def build_fig6c(writer):
    """Fig6C_htt_identity: BLASTN pident per HTT candidate × TE class."""
    path = S6_DIR / "results" / "htt_candidates_annotated.tsv"
    df = safe_read(path)
    if df is None:
        path2 = BASE_DIR / "EDF-TE" / "HTT" / "htt_candidates_annotated.tsv"
        df = safe_read(path2)
    if df is not None:
        cols = ["q_clade", "s_clade", "te_class", "pident", "anomaly_score", "alignment_length"]
        df = df[[c for c in cols if c in df.columns]]
    units = {
        "q_clade": "donor clade",
        "s_clade": "recipient clade",
        "te_class": "TE class",
        "pident": "BLASTN percent identity [85,100]",
        "anomaly_score": "d_expected − d_observed",
        "alignment_length": "bp",
    }
    write_sheet(writer, df, "Fig6C_htt_identity", units)


def build_fig6d(writer):
    """Fig6D_min_events: minimum parsimony HTT events per TE class."""
    path = S6_DIR / "results" / "htt_min_events.tsv"
    df = safe_read(path)
    if df is None:
        path2 = BASE_DIR / "EDF-TE" / "HTT" / "htt_min_events.tsv"
        df = safe_read(path2)
    units = {
        "te_class": "TE class",
        "n_minimum_events": "Dollo parsimony minimum independent transfer events",
    }
    write_sheet(writer, df, "Fig6D_min_events", units)


def build_fig6e(writer):
    """Fig6E_locus_cons: HTT locus conservation fraction counts."""
    path = S6_DIR / "results" / "htt_locus_conservation.tsv"
    df = safe_read(path)
    if df is None:
        path2 = BASE_DIR / "EDF-TE" / "HTT" / "htt_locus_conservation.tsv"
        df = safe_read(path2)
    if df is not None and "f_cons" in df.columns:
        # Summarise counts per category
        bins = pd.cut(df["f_cons"], bins=[-np.inf, 0.10, 0.30, np.inf],
                      labels=["non_conserved", "partial", "conserved"])
        counts = bins.value_counts().reset_index()
        counts.columns = ["conservation_category", "n_candidates"]
        counts["frac_total"] = counts["n_candidates"] / counts["n_candidates"].sum()
        df = counts
    units = {
        "conservation_category": "conserved (>=0.30) | partial (0.10-0.30) | non_conserved (<0.10)",
        "n_candidates": "count",
        "frac_total": "fraction [0,1]",
    }
    write_sheet(writer, df, "Fig6E_locus_cons", units)


def build_fig6f(writer):
    """Fig6F_htt_heatmap: inter-clade HTT candidate count matrix (wide)."""
    path = S6_DIR / "results" / "htt_candidates_annotated.tsv"
    df = safe_read(path)
    if df is None:
        path2 = BASE_DIR / "EDF-TE" / "HTT" / "htt_candidates_annotated.tsv"
        df = safe_read(path2)
    if df is not None and "q_clade" in df.columns and "s_clade" in df.columns:
        pivot = (df.groupby(["q_clade", "s_clade"]).size()
                 .unstack(fill_value=0).reset_index())
        pivot = pivot.rename(columns={"q_clade": "donor_clade"})
        df = pivot
    units = None  # wide format; headers are self-explanatory clade names
    write_sheet(writer, df, "Fig6F_htt_heatmap", units)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

SHEETS = [
    ("Fig1B_genomic_features", build_fig1b),
    ("Fig1C_pangenome", build_fig1c),
    ("Fig1D_synteny", build_fig1d),
    ("Fig2A_te_content", build_fig2a),
    ("Fig2B_pcoa", build_fig2b),
    ("Fig2C_kimura", build_fig2c),
    ("Fig2D_shannon", build_fig2d),
    ("Fig3A_dnmt", build_fig3a),
    ("Fig3B_methylation", build_fig3b),
    ("Fig3C_chromatin", build_fig3c),
    ("Fig3D_nucleosome", build_fig3d),
    ("Fig3E_srna_rpm", build_fig3e),
    ("Fig3F_srna_h3k9", build_fig3f),
    ("Fig4B_cooption_heatmap", build_fig4b),
    ("Fig4C_chimte", build_fig4c),
    ("Fig4D_tf_te", build_fig4d),
    ("Fig4E_fusions", build_fig4e),
    ("Fig4F_dam_amt1", build_fig4f),
    ("Fig5A_paralogs", build_fig5a),
    ("Fig5B_breakpoints", build_fig5b),
    ("Fig5C_te_class_bp", build_fig5c),
    ("Fig5D_dup_rate", build_fig5d),
    ("Fig6A_htt_chord", build_fig6a),
    ("Fig6B_htt_class", build_fig6b),
    ("Fig6C_htt_identity", build_fig6c),
    ("Fig6D_min_events", build_fig6d),
    ("Fig6E_locus_cons", build_fig6e),
    ("Fig6F_htt_heatmap", build_fig6f),
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
