#!/usr/bin/env python3
"""
Generate supplementary Excel tables for HTT (Horizontal Transposon Transfer) analysis.

Output: HTT_supplementary_tables.xlsx with 5 sheets:
  S1 - All 223 validated HTT pairs (one row per species pair)
  S2 - HTT communities summary (15 communities, one row each)
  S3 - HTT pair count matrix (clade × clade)
  S4 - Per-species HTT participation
  S5 - Minimal HTT events summary
"""

import os
import re
import sys
import glob
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

# ── Paths ────────────────────────────────────────────────────────────────────
RESULTS = "/legserv/Temp/Carlos/EDF-TE/results/HTT"
PAIRS_FILE    = f"{RESULTS}/08_htt_clustering/summary/htt_pairs_with_community.tsv"
SUMMARY_FILE  = f"{RESULTS}/08_htt_clustering/summary/htt_cluster_summary.tsv"
CANDIDATES    = f"{RESULTS}/06_htt_candidates/htt_candidates/htt_candidates.tsv"
CLADE_FILE    = f"{RESULTS}/clade_assignment.tsv"
BED_DIR       = f"{RESULTS}/04_te_annotation"
MIN_EVENTS    = f"{RESULTS}/09_minimal_htt_events/summary/minimal_events_by_component.tsv"
MIN_BY_NODE   = f"{RESULTS}/09_minimal_htt_events/summary/minimal_events_by_node.tsv"
COMM_TO_COMP  = f"{RESULTS}/09_minimal_htt_events/network/community_to_component.tsv"
OUT_XLSX      = f"{RESULTS}/HTT_supplementary_tables.xlsx"

# ── Clade full names (from plot_htt_chord.R — note BA/BS convention used in analysis) ─
CLADE_NAMES = {
    "HO": "Holozoa",
    "CH": "Chytridiomycota",
    "BL": "Blastocladiomycota",
    "BA": "Basidiomycota",
    "BS": "Basidiobolales",
    "MO": "Mortierellomycotina",
    "MU": "Mucoromycotina",
    "GL": "Glomeromycotina",
    "KX": "Kickxellomycotina",
    "EN": "Entomophthoromycotina",
    "ZO": "Zoopagomycota",
    "NU": "Nuclearida",
    "AP": "Aphelida",
    "MI": "Microsporidia",
    "CR": "Rozellomycota",
    "OL": "Olpidiomycota",
    "AS": "Ascomycota",
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_locus(locus_str):
    """Parse 'SPECIES|chr:start-end:strand' → (species, chr, start, end, strand)."""
    m = re.match(r"(.+)\|(.+):(\d+)-(\d+):([+-])", locus_str)
    if not m:
        return (locus_str, "", 0, 0, "")
    return m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)


def load_bed_index():
    """Load all filteredRepeats.bed files into a dict: species → list of (chr,s,e,te_type)."""
    idx = {}
    for path in glob.glob(f"{BED_DIR}/*.filteredRepeats.bed"):
        sp = os.path.basename(path).replace(".filteredRepeats.bed", "")
        records = []
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                try:
                    records.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
                except ValueError:
                    continue
        idx[sp] = records
    return idx


def lookup_te_type(bed_idx, species, chrom, start, end):
    """Return TE type for an exact coordinate match in the bed index."""
    for (c, s, e, te) in bed_idx.get(species, []):
        if c == chrom and s == start and e == end:
            return te
    return "unknown"


# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading data...")
comm_pairs_df = pd.read_csv(PAIRS_FILE, sep="\t")   # 49 pairs with community IDs
summary_df    = pd.read_csv(SUMMARY_FILE, sep="\t")
cands_df      = pd.read_csv(CANDIDATES,  sep="\t")  # 223 validated HTT pairs
clade_df      = pd.read_csv(CLADE_FILE,  sep="\t")
min_comp_df   = pd.read_csv(MIN_EVENTS,  sep="\t")
min_node_df   = pd.read_csv(MIN_BY_NODE, sep="\t")
comm_comp_df  = pd.read_csv(COMM_TO_COMP, sep="\t")

clade_map = dict(zip(clade_df["sample"], clade_df["clade"]))

# Build community lookup: pair_id → community_id
pair_to_comm = dict(zip(comm_pairs_df["pair_id"], comm_pairs_df["community_id"]))

print("Building TE annotation index (this may take a moment)...")
bed_idx = load_bed_index()
print(f"  Loaded bed files for {len(bed_idx)} species.")

# ── Resolve TE type for each pair ─────────────────────────────────────────────
def get_species_from_locus(loc):
    return loc.split("|")[0]

def get_te_type_for_pair(qseqid, sseqid):
    sp_q, chr_q, s_q, e_q, _ = parse_locus(qseqid)
    te_q = lookup_te_type(bed_idx, sp_q, chr_q, s_q, e_q)
    if te_q != "unknown":
        return te_q
    sp_s, chr_s, s_s, e_s, _ = parse_locus(sseqid)
    return lookup_te_type(bed_idx, sp_s, chr_s, s_s, e_s)

cands_df["te_type"] = cands_df.apply(
    lambda r: get_te_type_for_pair(r["qseqid"], r["sseqid"]), axis=1
)
cands_df["community_id"] = cands_df["pair_id"].map(pair_to_comm).fillna("—")

# ── Sheet S1: ALL 223 validated pairs ─────────────────────────────────────────
print("Building Sheet S1: All validated HTT pairs (223)...")

def expand_locus(loc):
    sp, chrom, s, e, strand = parse_locus(loc)
    return pd.Series({"species": sp, "chromosome": chrom,
                      "start": s, "end": e, "strand": strand})

q_exp = cands_df["qseqid"].apply(expand_locus).add_prefix("A_")
s_exp = cands_df["sseqid"].apply(expand_locus).add_prefix("B_")

s1_df = pd.concat([
    cands_df[["pair_id", "community_id"]].reset_index(drop=True),
    q_exp.reset_index(drop=True),
    pd.Series(cands_df["qclade"].values,
              name="Clade_A_code"),
    pd.Series([CLADE_NAMES.get(c, c) for c in cands_df["qclade"]],
              name="Clade_A_name"),
    s_exp.reset_index(drop=True),
    pd.Series(cands_df["sclade"].values,
              name="Clade_B_code"),
    pd.Series([CLADE_NAMES.get(c, c) for c in cands_df["sclade"]],
              name="Clade_B_name"),
    pd.Series(cands_df["te_type"].values,  name="TE_type"),
    pd.Series(cands_df["te_ks"].values,    name="TE_Ks"),
    pd.Series(cands_df["busco_ks_threshold"].values, name="BUSCO_Ks_threshold"),
    pd.Series(cands_df["mrca_node"].values, name="MRCA_node"),
], axis=1).sort_values(["community_id", "pair_id"])

# ── Sheet S2: Communities summary (15) ────────────────────────────────────────
print("Building Sheet S2: Communities summary (15)...")

comm_pairs_df["te_type"] = comm_pairs_df.apply(
    lambda r: get_te_type_for_pair(r["qseqid"], r["sseqid"]), axis=1
)

rows_s2 = []
for comm_id, grp in comm_pairs_df.groupby("community_id"):
    clade_a = grp["cladeA"].iloc[0]
    clade_b = grp["cladeB"].iloc[0]
    sp_a = sorted(set(grp["qseqid"].apply(get_species_from_locus)))
    sp_b = sorted(set(grp["sseqid"].apply(get_species_from_locus)))
    te_counts = grp["te_type"].value_counts()
    te_type = te_counts.index[0] if len(te_counts) else "unknown"
    rows_s2.append({
        "Community_ID":       comm_id,
        "Clade_A_code":       clade_a,
        "Clade_A_name":       CLADE_NAMES.get(clade_a, clade_a),
        "Clade_B_code":       clade_b,
        "Clade_B_name":       CLADE_NAMES.get(clade_b, clade_b),
        "n_species_pairs":    len(grp),
        "n_species_cladeA":   len(sp_a),
        "n_species_cladeB":   len(sp_b),
        "Species_cladeA":     "; ".join(sp_a),
        "Species_cladeB":     "; ".join(sp_b),
        "TE_type":            te_type,
    })

s2_df = pd.DataFrame(rows_s2).sort_values("Community_ID")

# ── Sheet S3: Clade pair count matrix (from all 223 pairs) ────────────────────
print("Building Sheet S3: Clade pair matrix...")

all_clades = sorted(set(cands_df["qclade"]) | set(cands_df["sclade"]))
matrix = pd.DataFrame(0, index=all_clades, columns=all_clades)
for _, row in cands_df.iterrows():
    a, b = row["qclade"], row["sclade"]
    matrix.loc[a, b] += 1
    if a != b:
        matrix.loc[b, a] += 1

matrix.index   = [f"{c} ({CLADE_NAMES.get(c, c)})" for c in matrix.index]
matrix.columns = [f"{c} ({CLADE_NAMES.get(c, c)})" for c in matrix.columns]

# ── Sheet S4: Per-species participation (from all 223 pairs) ──────────────────
print("Building Sheet S4: Species participation...")

sp_counts = {}
for _, row in cands_df.iterrows():
    sp_q = get_species_from_locus(row["qseqid"])
    sp_s = get_species_from_locus(row["sseqid"])
    sp_counts.setdefault(sp_q, {"as_A": 0, "as_B": 0})
    sp_counts.setdefault(sp_s, {"as_A": 0, "as_B": 0})
    sp_counts[sp_q]["as_A"] += 1
    sp_counts[sp_s]["as_B"] += 1

rows_s4 = []
for sp, cnts in sorted(sp_counts.items()):
    clade_code = clade_map.get(sp, "?")
    rows_s4.append({
        "Species":          sp,
        "Clade_code":       clade_code,
        "Clade_name":       CLADE_NAMES.get(clade_code, clade_code),
        "n_pairs_as_A":     cnts["as_A"],
        "n_pairs_as_B":     cnts["as_B"],
        "n_pairs_total":    cnts["as_A"] + cnts["as_B"],
    })
s4_df = pd.DataFrame(rows_s4).sort_values("n_pairs_total", ascending=False)

# ── Sheet S5: Minimal HTT events ──────────────────────────────────────────────
print("Building Sheet S5: Minimal HTT events...")

# Build per-component clade and community info
# community summary dict: comm_id → (cladeA, cladeB, cladeA_name, cladeB_name)
comm_info = {}
for _, r in comm_pairs_df.drop_duplicates("community_id").iterrows():
    comm_info[r["community_id"]] = (r["cladeA"], r["cladeB"])

comp_rows = []
for comp_id, grp in comm_comp_df.groupby("component_id"):
    comms = sorted(grp["community_id"].tolist())
    clades_seen = set()
    pair_labels = []
    for c in comms:
        if c in comm_info:
            a, b = comm_info[c]
            clades_seen.update([a, b])
            pair_labels.append(f"{a}–{b}")
    clade_names_list = sorted(
        f"{c} ({CLADE_NAMES.get(c, c)})" for c in clades_seen
    )
    # Get min_events from min_comp_df
    row_mc = min_comp_df[min_comp_df["component_id"] == comp_id]
    min_ev = int(row_mc["min_events_component"].iloc[0]) if len(row_mc) else "?"
    n_comms = int(row_mc["n_communities"].iloc[0]) if len(row_mc) else len(comms)
    comp_rows.append({
        "Component_ID":          comp_id,
        "min_HTT_events":        min_ev,
        "n_communities":         n_comms,
        "Communities":           "; ".join(comms),
        "n_clades_involved":     len(clades_seen),
        "Clades_involved":       "; ".join(clade_names_list),
        "n_distinct_clade_pairs": len(pair_labels),
        "Clade_pairs":           "; ".join(pair_labels),
    })

s5a_df = pd.DataFrame(comp_rows).sort_values("Component_ID")

s5b_df = min_node_df.rename(columns={
    "node":              "MRCA_node",
    "min_events_total":  "min_HTT_events_at_node",
})

# ── Write Excel ───────────────────────────────────────────────────────────────
print(f"Writing {OUT_XLSX} ...")

with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
    s1_df.to_excel(writer, sheet_name="S1_All_HTT_pairs",        index=False)
    s2_df.to_excel(writer, sheet_name="S2_Communities_summary",  index=False)
    matrix.to_excel(writer, sheet_name="S3_Clade_pair_matrix")
    s4_df.to_excel(writer, sheet_name="S4_Species_participation", index=False)
    s5a_df.to_excel(writer, sheet_name="S5_Minimal_events",       index=False)
    s5b_df.to_excel(writer, sheet_name="S5b_Events_by_node",      index=False)

# ── Apply formatting ───────────────────────────────────────────────────────────
print("Applying formatting...")
wb = load_workbook(OUT_XLSX)

HEADER_FILL  = PatternFill("solid", fgColor="1F4E79")
HEADER_FONT  = Font(color="FFFFFF", bold=True)
EVEN_FILL    = PatternFill("solid", fgColor="D6E4F0")
ODD_FILL     = PatternFill("solid", fgColor="FFFFFF")
BORDER       = Border(
    bottom=Side(style="thin", color="BBBBBB"),
    right=Side(style="thin", color="BBBBBB"),
)
MATRIX_FILL  = PatternFill("solid", fgColor="1F4E79")
NONZERO_FILL = PatternFill("solid", fgColor="AED6F1")
ZERO_FILL    = PatternFill("solid", fgColor="F0F4F8")

def format_sheet(ws, freeze_row=2, col_widths=None):
    for col_idx, cell in enumerate(ws[1], 1):
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
    for row_idx, row in enumerate(ws.iter_rows(min_row=2), 1):
        fill = EVEN_FILL if row_idx % 2 == 0 else ODD_FILL
        for cell in row:
            cell.fill = fill
            cell.border = BORDER
            cell.alignment = Alignment(vertical="top", wrap_text=False)
    if col_widths:
        for col_letter, width in col_widths.items():
            ws.column_dimensions[col_letter].width = width
    if freeze_row:
        ws.freeze_panes = f"A{freeze_row}"


# S1
ws1 = wb["S1_All_HTT_pairs"]
format_sheet(ws1, col_widths={
    "A": 18, "B": 18, "C": 20, "D": 12, "E": 8, "F": 10, "G": 8,
    "H": 12, "I": 24, "J": 20, "K": 12, "L": 8, "M": 10, "O": 8,
    "P": 12, "Q": 24, "R": 10, "S": 12, "T": 14, "U": 18,
})

# S2
ws2 = wb["S2_Communities_summary"]
format_sheet(ws2, col_widths={
    "A": 18, "B": 12, "C": 24, "D": 12, "E": 24,
    "F": 14, "G": 14, "H": 14, "I": 40, "J": 40, "K": 12,
})

# S3 — matrix with colour coding
ws3 = wb["S3_Clade_pair_matrix"]
n_rows = ws3.max_row
n_cols = ws3.max_column
for row in ws3.iter_rows(min_row=1, max_row=n_rows, min_col=1, max_col=n_cols):
    for cell in row:
        if cell.row == 1 or cell.column == 1:
            cell.fill = MATRIX_FILL
            cell.font = Font(color="FFFFFF", bold=True)
            cell.alignment = Alignment(horizontal="center")
        else:
            try:
                val = int(cell.value) if cell.value is not None else 0
                cell.fill = NONZERO_FILL if val > 0 else ZERO_FILL
                cell.alignment = Alignment(horizontal="center")
                if val == 0:
                    cell.value = ""
            except (TypeError, ValueError):
                pass
for col_idx in range(1, n_cols + 1):
    ws3.column_dimensions[get_column_letter(col_idx)].width = 22
ws3.freeze_panes = "B2"

# S4
ws4 = wb["S4_Species_participation"]
format_sheet(ws4, col_widths={
    "A": 22, "B": 12, "C": 24, "D": 14, "E": 14, "F": 14,
})

# S5
ws5 = wb["S5_Minimal_events"]
format_sheet(ws5, col_widths={
    "A": 20, "B": 16, "C": 16, "D": 60, "E": 16, "F": 70, "G": 22, "H": 60,
})
ws5b = wb["S5b_Events_by_node"]
format_sheet(ws5b, col_widths={"A": 14, "B": 22})

# Sheet tab colours
wb["S1_All_HTT_pairs"].sheet_properties.tabColor        = "1F4E79"
wb["S2_Communities_summary"].sheet_properties.tabColor  = "2E75B6"
wb["S3_Clade_pair_matrix"].sheet_properties.tabColor    = "2E86C1"
wb["S4_Species_participation"].sheet_properties.tabColor= "5DADE2"
wb["S5_Minimal_events"].sheet_properties.tabColor       = "76B041"
wb["S5b_Events_by_node"].sheet_properties.tabColor      = "A9CCE3"

wb.save(OUT_XLSX)
print(f"\nDone! Saved: {OUT_XLSX}")
print(f"  S1 All HTT pairs:        {len(s1_df)} validated pairs")
print(f"  S2 Communities summary:  {len(s2_df)} HTT communities")
print(f"  S3 Clade pair matrix:    {len(all_clades)} clades")
print(f"  S4 Species:              {len(s4_df)} species with HTT events")
print(f"  S5 Minimal events:       {min_comp_df['min_events_component'].sum()} events across {len(min_comp_df)} components")
