#!/usr/bin/env python3
"""
Figure: Sex locus summary + synteny diagrams with TE content
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow, FancyArrowPatch
import matplotlib.gridspec as gridspec
from pathlib import Path

WORKDIR = Path("/data/carlos/EDF-TE/TE_analysis/sex_locus")

# ─── Load data ────────────────────────────────────────────────────────────────
# 2026-04-24: Use blastp-derived master v3 (sex_locus_master_v3.tsv) — built
# from genome-wide blastp of user references vs every Mucorales proteome.
#   ▸ Sex call    : top hit of {sexP_user,sexM_user} with E<1e-10, alen≥80,
#                   Δbitscore≥20  (else _tentative)
#   ▸ Locus       : sex_gene_start − 5 kb  →  sex_gene_end + 5 kb
#   ▸ TE filter   : drops Simple_repeat, Low_complexity, Unknown, Tandem,
#                   Satellite, rRNA/tRNA/snRNA  (only true TEs are summed)
#   ▸ Orthologue  : E<1e-15  ∧  alen≥50%×qlen  ∧  pident≥25%  ∧  same scaffold
#                   ∧  ≤50 kb from sex-gene midpoint  (true sex-locus block,
#                   not a paralog elsewhere in the genome).
_V3_MASTER = WORKDIR / "blastp_proteomes" / "sex_locus_master_v3.tsv"
_V3_IDIO   = WORKDIR / "blastp_proteomes" / "idiomorphic_regions_v3.bed"
_V3_TE     = WORKDIR / "blastp_proteomes" / "TE_in_idio_v3.tsv"
_V3_VIEW   = WORKDIR / "blastp_proteomes" / "TE_in_view_v3.tsv"
_V3_FL     = WORKDIR / "blastp_proteomes" / "full_length_TEs_v3.tsv"
_V3_ORTHO  = WORKDIR / "blastp_proteomes" / "locus_orthologs_v3.tsv"
if _V3_MASTER.exists():
    master = pd.read_csv(_V3_MASTER, sep="\t")
    print(f"[INFO] Master v3: {_V3_MASTER.name}  ({len(master)} rows)")
else:
    master  = pd.read_csv(WORKDIR / "sex_locus_master_table_clean.tsv", sep="\t")
tblastn = pd.read_csv(WORKDIR / "tblastn_sex_markers_per_species.tsv", sep="\t")
if _V3_TE.exists():
    te_det = pd.read_csv(_V3_TE, sep="\t")
    print(f"[INFO] TE detail v3: {_V3_TE.name}  ({len(te_det)} rows)")
else:
    te_det = pd.read_csv(WORKDIR / "TE_in_locus_detailed.tsv", sep="\t")
# Wider TE detail (±15 kb of sex midpoint) used for figure rendering
if _V3_VIEW.exists():
    te_view = pd.read_csv(_V3_VIEW, sep="\t")
    print(f"[INFO] TE view v3: {_V3_VIEW.name}  ({len(te_view)} rows, ±15 kb)")
else:
    te_view = te_det.copy()
# Full-length TE annotations (≥ 3.5 kb single annotation in ±15 kb)
if _V3_FL.exists():
    fl_te = pd.read_csv(_V3_FL, sep="\t")
    print(f"[INFO] Full-length TEs v3: {_V3_FL.name}  ({len(fl_te)} rows)")
else:
    fl_te = pd.DataFrame(columns=["species","te_class","te_start","te_end","te_len","dist_to_sex"])
if _V3_ORTHO.exists():
    locus_orthos = pd.read_csv(_V3_ORTHO, sep="\t")
    print(f"[INFO] Locus orthologs v3: {_V3_ORTHO.name}  ({len(locus_orthos)} rows)")
else:
    locus_orthos = pd.DataFrame(columns=["species","marker","scaffold","start",
                                          "end","strand","pident","alen_frac",
                                          "evalue","bitscore"])

# 2026-04-24: also load user-provided helicase + tptA tBLASTn hits
_USER_MARKERS = WORKDIR / "user_refs" / "markers_user_helicase_tptA.tsv"
if _USER_MARKERS.exists():
    _um = pd.read_csv(_USER_MARKERS, sep="\t")
    # Normalize: keep one strand-corrected start/end and merge into `tblastn`
    _um["start"]  = _um[["sstart","send"]].min(axis=1)
    _um["end"]    = _um[["sstart","send"]].max(axis=1)
    _um["strand"] = np.where(_um["sstart"] < _um["send"], "+", "-")
    _um = _um.rename(columns={"length":"alen"})
    _um = _um[["species","marker","scaffold","start","end","strand","pident","evalue","bitscore"]]
    print(f"[INFO] Loaded {len(_um)} user-marker HSPs from {_USER_MARKERS.name}")
    user_markers = _um
else:
    user_markers = pd.DataFrame(columns=["species","marker","scaffold","start","end","strand","pident","evalue","bitscore"])

# 2026-04-24: Append MU402 helicase + tptA HSPs (relabel species code = Mlus402)
_MU402_MARKERS = WORKDIR / "user_refs" / "mu402_user_markers.tsv"
if _MU402_MARKERS.exists():
    _mm = pd.read_csv(_MU402_MARKERS, sep="\t")
    _mm["species"] = "Mlus402"
    _mm["start"]   = _mm[["sstart","send"]].min(axis=1)
    _mm["end"]     = _mm[["sstart","send"]].max(axis=1)
    _mm["strand"]  = np.where(_mm["sstart"] < _mm["send"], "+", "-")
    _mm = _mm[["species","marker","scaffold","start","end","strand","pident","evalue","bitscore"]]
    user_markers = pd.concat([user_markers, _mm], ignore_index=True)
    print(f"[INFO] Added {len(_mm)} MU402 marker HSPs")

# 2026-04-24: load gene annotations (JGI-style GFF) for background grey arrows
GFF_DIR = Path("/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff")
def _load_gff_genes(species, scaffold, vs, ve):
    """Return list of (start, end, strand, proteinId, product) for CDSs in window.
    Tries JGI-style GFF first (CDS rows + proteinId attribute), then GFF3 (gene rows)."""
    if scaffold is None or pd.isna(vs) or pd.isna(ve):
        return []
    candidates = [GFF_DIR / f"MU_{species}.gff", GFF_DIR / f"MU_{species}.gff3"]
    fp = next((p for p in candidates if p.exists()), None)
    if fp is None:
        return []
    genes = {}  # pid -> [s, e, strand, product]
    gff3_rows = []
    try:
        with open(fp) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[0] != scaffold:
                    continue
                try:
                    s, e = int(f[3]), int(f[4])
                except ValueError:
                    continue
                if e < vs or s > ve:
                    continue
                strand = f[6]
                attrs  = f[8]
                if f[2] == "CDS" and "proteinId" in attrs:
                    pid = None; prod = ""
                    for tok in attrs.split(";"):
                        tok = tok.strip()
                        if tok.startswith("proteinId "):
                            pid = tok.split()[1].strip().strip('"')
                        elif tok.startswith("product_name "):
                            prod = tok.split(" ", 1)[1].strip().strip('"')
                    if pid is None: continue
                    if pid in genes:
                        genes[pid][0] = min(genes[pid][0], s)
                        genes[pid][1] = max(genes[pid][1], e)
                    else:
                        genes[pid] = [s, e, strand, prod]
                elif f[2] == "gene":
                    pid = None
                    for tok in attrs.split(";"):
                        if tok.startswith("ID="):
                            pid = tok.split("=", 1)[1].strip()
                            break
                    gff3_rows.append((s, e, strand, pid or "", ""))
    except Exception as ex:
        print(f"[WARN] GFF parse failed for {species}: {ex}")
        return []
    if genes:
        return sorted([(s, e, st, pid, pr) for pid, (s, e, st, pr) in genes.items()])
    return sorted(gff3_rows)

idio    = pd.read_csv(WORKDIR / "idiomorphic_regions.bed", sep="\t", header=None,
                      names=["scaffold","start","end","label"])
idio["species"] = idio["label"].str.split("|").str[0]
# 2026-04-24: append v3 (blastp-derived) idiomorphic windows
if _V3_IDIO.exists():
    _idio_v2 = pd.read_csv(_V3_IDIO, sep="\t", header=None,
                            usecols=[0,1,2,3],
                            names=["scaffold","start","end","label"])
    _idio_v2["species"] = _idio_v2["label"].str.split("|").str[0]
    # Replace any species present in v3 (newer + canonical)
    idio = idio[~idio["species"].isin(_idio_v2["species"])].reset_index(drop=True)
    idio = pd.concat([idio, _idio_v2], ignore_index=True)
    idio["start"] = pd.to_numeric(idio["start"], errors="coerce").astype("Int64")
    idio["end"]   = pd.to_numeric(idio["end"],   errors="coerce").astype("Int64")
    print(f"[INFO] Merged {len(_idio_v2)} v3 idiomorphic regions")

# 2026-04-24: Drop Mucend1 entirely (assembly artifact — sex locus split across
# scaffolds 6/10/4; user requested removal from all downstream analyses/figure).
_DROP_SP = {"Mucend1"}
master       = master[~master["species"].isin(_DROP_SP)].reset_index(drop=True)
tblastn      = tblastn[~tblastn["species"].isin(_DROP_SP)].reset_index(drop=True)
te_det       = te_det[~te_det["species"].isin(_DROP_SP)].reset_index(drop=True)
idio         = idio[~idio["species"].isin(_DROP_SP)].reset_index(drop=True)
user_markers = user_markers[~user_markers["species"].isin(_DROP_SP)].reset_index(drop=True)
print(f"[INFO] Dropped species: {sorted(_DROP_SP)} (assembly-artifact split locus)")

# 2026-04-24: Inject Rhizopus microsporus (Rhimi59) — sexP confirmed by user
# blastp of sexP/sexM refs vs JGI proteinId 1981299 (annotated as 'expressed
# protein' / helicase-like, but actually a fused gene model containing
# helicase + rnhA + sexP HMG domain). The HMG portion (AA 42-323) maps to
# scaffold_1:367788-368324 (-) with sexP_user, e=3.7e-32, bs=124.
def _inject_rhimi59(master, idio, te_det, tblastn, user_markers):
    sp = "Rhimi59"
    if sp in set(master["species"]):
        # Already present (probably as not_identified) → upgrade in place
        master = master[master["species"] != sp].reset_index(drop=True)
    sex_scaf, sex_s, sex_e, sex_str = "scaffold_1", 367788, 368324, "-"
    # Idiomorphic region: rnhA-helicase-sexP cassette + small flanks.
    idio_scaf, idio_s, idio_e = "scaffold_1", 358000, 373000
    region_len = idio_e - idio_s
    # TE intervals from per-window RepeatMasker table
    te_rows = []
    rm = WORKDIR / "synteny_figure" / "rmask" / f"{sp}_TE.tsv"
    if rm.exists():
        try:
            df = pd.read_csv(rm, sep="\t", header=None,
                             names=["region","s_local","e_local","pid","alen","te","fid","ev","bs"])
            df["win_start"] = df["region"].str.extract(r":(\d+)-").astype(int)
            df["abs_s"] = df["win_start"] + df[["s_local","e_local"]].min(axis=1)
            df["abs_e"] = df["win_start"] + df[["s_local","e_local"]].max(axis=1)
            df = df[(df["abs_e"] >= idio_s) & (df["abs_s"] <= idio_e)]
            for _, r in df.iterrows():
                cls = str(r["te"]).split("#")[-1] if "#" in str(r["te"]) else str(r["te"])
                te_rows.append({
                    "species": sp, "region": "idiomorphic",
                    "te_class": cls, "overlap_bp": int(r["abs_e"] - r["abs_s"]),
                    "te_start": int(r["abs_s"]), "te_end": int(r["abs_e"]),
                })
        except Exception as ex:
            print(f"[WARN] Rhimi59 TE parse failed: {ex}")
    te_bp = sum(r["overlap_bp"] for r in te_rows)
    # Master row (P confirmed)
    new_row = {c: np.nan for c in master.columns}
    new_row.update({
        "species": sp, "mating_type": "P",
        "sex_scaffold": sex_scaf, "sex_start": float(sex_s), "sex_end": float(sex_e),
        "sex_strand": sex_str, "sex_size_aa": 178.0,
        "TE_bp_idiomorphic": float(te_bp),
        "TE_pct_idiomorphic": float(100.0 * te_bp / region_len),
    })
    master = pd.concat([master, pd.DataFrame([new_row])], ignore_index=True)
    # Idiomorphic region row
    idio = pd.concat([idio, pd.DataFrame([{
        "scaffold": idio_scaf, "start": idio_s, "end": idio_e,
        "label": f"{sp}|idiomorphic", "species": sp,
    }])], ignore_index=True)
    # te_det rows
    if te_rows:
        te_det = pd.concat([te_det, pd.DataFrame(te_rows)], ignore_index=True)
    # Helicase user-marker (upstream half of fused gene model 1981299)
    new_um = pd.DataFrame([{
        "species": sp, "marker": "helicase", "scaffold": sex_scaf,
        "start": 363146, "end": 365805, "strand": "-",
        "pident": np.nan, "evalue": np.nan, "bitscore": np.nan,
    }])
    user_markers = pd.concat([user_markers, new_um], ignore_index=True)
    print(f"[INFO] Injected Rhimi59 sexP locus (scaffold_1:{sex_s}-{sex_e}, "
          f"TE in idio = {te_bp} bp / {region_len} bp = {100.0*te_bp/region_len:.2f}%)")
    return master, idio, te_det, tblastn, user_markers

master, idio, te_det, tblastn, user_markers = _inject_rhimi59(
    master, idio, te_det, tblastn, user_markers)

# 2026-04-24: Inject Mucor lusitanicus CBS v3 (MU402) — sexM strain.
# tBLASTn coordinates (user-supplied sexP/sexM/helicase/tptA queries vs mu402GenomeJGI.fa).
def _inject_mu402(master, idio, te_det, tblastn):
    sp = "Mlus402"
    # Best sexM HSP: scaffold_1 : 2225624-2225884 (39.1% pident, 1.78e-5)
    sex_scaf, sex_s, sex_e, sex_str = "scaffold_1", 2225624, 2225884, "+"
    # Idiomorphic region defined as tptA-end → helicase-end on the syntenic block.
    idio_scaf, idio_s, idio_e = "scaffold_1", 2223400, 2230700
    # TE table from the windowed RM run (window starts at 2,175,952 on scaffold_1)
    te_rows = []
    rm = WORKDIR / "synteny_figure" / "rmask" / "mu402_TE.tsv"
    if rm.exists():
        try:
            df = pd.read_csv(rm, sep="\t", header=None,
                             names=["region","s_local","e_local","pid","alen","te","fid","ev","bs"])
            # parse window absolute start from region string "...scaffold_1:2175952-2275952"
            df["win_start"] = df["region"].str.extract(r":(\d+)-").astype(int)
            df["abs_s"] = df["win_start"] + df[["s_local","e_local"]].min(axis=1)
            df["abs_e"] = df["win_start"] + df[["s_local","e_local"]].max(axis=1)
            df = df[(df["abs_e"] >= idio_s) & (df["abs_s"] <= idio_e)]
            for _, r in df.iterrows():
                # Class formatting: rnd-1_family-12#RC/Helitron -> RC/Helitron
                cls = str(r["te"]).split("#")[-1] if "#" in str(r["te"]) else str(r["te"])
                te_rows.append({
                    "species": sp, "region": "idiomorphic",
                    "te_class": cls, "overlap_bp": int(r["abs_e"] - r["abs_s"]),
                    "te_start": int(r["abs_s"]), "te_end": int(r["abs_e"]),
                })
        except Exception as ex:
            print(f"[WARN] MU402 TE parse failed: {ex}")
    te_bp = sum(r["overlap_bp"] for r in te_rows)
    region_len = max(1, idio_e - idio_s)
    # master row
    new_row = {c: np.nan for c in master.columns}
    new_row.update({
        "species": sp, "mating_type": "M",
        "sex_scaffold": sex_scaf, "sex_start": float(sex_s), "sex_end": float(sex_e),
        "sex_strand": sex_str, "sex_size_aa": 205.0,
        "TE_bp_idiomorphic": float(te_bp),
        "TE_pct_idiomorphic": float(100.0 * te_bp / region_len),
    })
    master = pd.concat([master, pd.DataFrame([new_row])], ignore_index=True)
    # idio row
    idio = pd.concat([idio, pd.DataFrame([{
        "scaffold": idio_scaf, "start": idio_s, "end": idio_e,
        "label": f"{sp}|idiomorphic", "species": sp,
    }])], ignore_index=True)
    # te_det rows
    if te_rows:
        te_det = pd.concat([te_det, pd.DataFrame(te_rows)], ignore_index=True)
    # tBLASTn marker rows (rnhA/tptA on scaffold_1 syntenic block)
    new_tb = pd.DataFrame([
        {"species": sp, "marker": "tptA", "scaffold": sex_scaf,
         "start": 2223422, "end": 2224642, "strand": "-", "pident": 62.7, "evalue": 1.5e-156, "bitscore": 493.0},
    ])
    tblastn = pd.concat([tblastn, new_tb], ignore_index=True)
    return master, idio, te_det, tblastn

master, idio, te_det, tblastn = _inject_mu402(master, idio, te_det, tblastn)

# ── 2026-05-01: replace the "idiomorphic region" TE% with a fixed ±10 kb
# flanking window centred on the sex-gene midpoint (FLANK_BP = 10 000 bp,
# total window = 20 000 bp).  This makes panel A directly comparable across
# species (the idiomorphic region had heterogeneous lengths from blastp
# orthology calls).  TE intervals are taken from te_view (collected at
# ±15 kb in build_master_v3.py, so it fully contains the new ±10 kb window).
FLANK_BP = 10000
def _compute_flank_te_pct(master_df, te_view_df, flank=FLANK_BP):
    rows = []
    for _, r in master_df.iterrows():
        sp = r["species"]
        if pd.isna(r.get("sex_start")) or pd.isna(r.get("sex_end")):
            rows.append(None); continue
        sex_mid = (float(r["sex_start"]) + float(r["sex_end"])) / 2.0
        ws, we = sex_mid - flank, sex_mid + flank
        sub = te_view_df[te_view_df["species"] == sp]
        if len(sub) == 0:
            rows.append(0.0); continue
        # Clip each TE interval to the ±flank window and sum bp
        s = sub["te_start"].astype(float).clip(lower=ws)
        e = sub["te_end"].astype(float).clip(upper=we)
        ovl = (e - s + 1).clip(lower=0)
        rows.append(float(100.0 * ovl.sum() / (2 * flank)))
    return rows

master["TE_pct_flank10"] = _compute_flank_te_pct(master, te_view)
master["TE_bp_flank10"]  = (master["TE_pct_flank10"].fillna(0) * (2 * FLANK_BP) / 100.0).round().astype("Int64")
print(f"[INFO] Computed TE_pct_flank10 (±{FLANK_BP} bp window around sex-gene midpoint) "
      f"for {master['TE_pct_flank10'].notna().sum()} species")

identified_full = master[master["mating_type"].isin(["P", "M"])].copy()  # only CONFIRMED sexM / sexP
# 2026-04-28: Fenlin1 (second Fenlin strain) is excluded from the GENERAL figure;
# it is shown separately in the dedicated Fenlin synteny figure.  Data structures
# above (master/idio/te_det/tblastn) keep Fenlin1 so the strain-comparison figure
# can still access it.
_DROP_FROM_MAIN = {"Fenlin1"}
identified_full = identified_full[~identified_full["species"].isin(_DROP_FROM_MAIN)].copy()
identified_full = identified_full.sort_values("TE_pct_flank10", ascending=False)

# Keep ONLY the top-10 species by TE_pct_flank10 (user request 2026-05-01)
TOP_N = 10
identified = identified_full.head(TOP_N).copy()
# 2026-04-24: also force-include Mlus402 (M. lusitanicus CBS v3) regardless of rank
_force = ["Mlus402"]
for _sp in _force:
    if _sp in identified["species"].values:
        continue
    extra = identified_full[identified_full["species"] == _sp]
    if len(extra):
        identified = pd.concat([identified, extra], ignore_index=True)
# Re-sort ascending so panel A reads bottom→top from low→high TE%
identified = identified.sort_values("TE_pct_flank10", ascending=True).reset_index(drop=True)
print(f"[INFO] Keeping top-{TOP_N} species (out of {len(identified_full)} confirmed sexM/sexP) by TE% (±{FLANK_BP} bp):")
for _, r in identified.iterrows():
    print(f"        {r['species']:>10s}  {r['mating_type']}  {r['TE_pct_flank10']:.1f}% TE")

# ── Emit a TSV with the genomic coordinates for the selected 10 ────────────
# Save in BOTH WORKDIR (data dir) and /legserv MAT_locus/tables (visible tree).
TABLES_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/MAT_locus/tables")
TABLES_DIR.mkdir(parents=True, exist_ok=True)
TOP10_OUT_LOCAL = WORKDIR / "sex_locus_top10_coords.tsv"
TOP10_OUT       = TABLES_DIR / "sex_locus_top10_coords.tsv"
_top10_cols = ["species", "mating_type", "sex_scaffold", "sex_start", "sex_end",
               "sex_strand", "sex_size_aa",
               "TE_bp_idiomorphic", "TE_pct_idiomorphic",
               "TE_bp_flank10", "TE_pct_flank10"]
_idio_lookup = (idio.drop_duplicates("species", keep="last")
                    .set_index("species")[["scaffold", "start", "end"]]
                    .to_dict("index"))
_top10_rows = []
for _, r in identified.iterrows():
    sp = r["species"]
    idio_r = _idio_lookup.get(sp, {})
    _top10_rows.append({
        **{c: r.get(c) for c in _top10_cols},
        "idiomorphic_scaffold": idio_r.get("scaffold"),
        "idiomorphic_start":    idio_r.get("start"),
        "idiomorphic_end":      idio_r.get("end"),
    })
_df_top10 = pd.DataFrame(_top10_rows)
_df_top10.to_csv(TOP10_OUT_LOCAL, sep="\t", index=False)
_df_top10.to_csv(TOP10_OUT,       sep="\t", index=False)
print(f"[OK] Top-10 coordinates → {TOP10_OUT}")
print(f"                          {TOP10_OUT_LOCAL}")

# ── Emit an AA fasta restricted to the 10 selected species ─────────────────
# Source: per-species sex-candidate proteins curated by the HMG/tBLASTN pipeline.
# Headers carry the mating type and TE% so the file is self-describing.
import re as _re
_AA_SRC = Path("/legserv/Temp/Carlos/EDF-TE/results/MAT_locus/data/reference_loci"
               "/all_hmg_sex_candidates_clean.faa")
TOP10_FAA_LOCAL = WORKDIR / "sex_locus_top10_proteins.faa"
TOP10_FAA       = TABLES_DIR / "sex_locus_top10_proteins.faa"
_keep_species = set(identified["species"])
_mt_map = dict(zip(identified["species"], identified["mating_type"]))
_te_map = dict(zip(identified["species"], identified["TE_pct_flank10"]))
if _AA_SRC.exists():
    _written = 0
    _seen_species = set()
    with open(_AA_SRC) as _fh, \
         open(TOP10_FAA_LOCAL, "w") as _out_l, \
         open(TOP10_FAA,       "w") as _out_g:
        _take = False
        for _line in _fh:
            if _line.startswith(">"):
                _hdr_full = _line[1:].rstrip()
                _hdr_sp   = _hdr_full.split("|")[0]
                _take = _hdr_sp in _keep_species
                if _take:
                    _mt  = _mt_map.get(_hdr_sp, "?")
                    _tep = _te_map.get(_hdr_sp, float("nan"))
                    _new = (f">{_hdr_full}  mating_type={_mt}  "
                            f"TE_pct_idio={_tep:.1f}\n")
                    _out_l.write(_new); _out_g.write(_new)
                    _seen_species.add(_hdr_sp)
                    _written += 1
            elif _take:
                _out_l.write(_line); _out_g.write(_line)
    print(f"[OK] AA fasta ({_written} sequences from {len(_seen_species)} species) →")
    print(f"        {TOP10_FAA}")
    print(f"        {TOP10_FAA_LOCAL}")
    _missing = _keep_species - _seen_species
    if _missing:
        print(f"[WARN] No AA sequence found for: {sorted(_missing)}")
else:
    print(f"[WARN] AA source not found, skipped: {_AA_SRC}")

# Canonical TE-class palette (matches cladogram figure)
import sys as _sys, os as _os
_sys.path.insert(0, "/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts")
from clade_palette import TE_CLASS_COLORS as _TE_PALETTE  # noqa: E402

def _canonical_class(name: str) -> str:
    """Map any granular TE class string to canonical cladogram class."""
    if not isinstance(name, str):
        return "Unclassified"
    head = name.split("/")[0].strip()
    if head in ("LTR",):                       return "LTR"
    if head in ("LINE", "Non-LTR"):            return "LINE"
    if head in ("DNA",):                       return "DNA"
    if head in ("SINE",):                      return "SINE"
    if head in ("RC", "Helitron"):             return "Rolling Circle"
    if head in ("PLE", "Penelope"):            return "Penelope"
    if head in ("Simple_repeat", "Low_complexity", "Satellite", "Unknown"):
        return "Unclassified"
    return "Unclassified"

# ─── Color scheme ─────────────────────────────────────────────────────────────
MT_COLOR = {
    "P":           "#E8671A",   # orange
    "P_tentative": "#F4A86A",   # light orange
    "M":           "#C0392B",   # red
    "M_tentative": "#E57373",   # light red
}

TE_CLASS_COLOR = {
    "LTR/Gypsy":       "#8E44AD",
    "LTR/Copia":       "#9B59B6",
    "LTR/ERV1":        "#A569BD",
    "LTR/DIRS":        "#7D3C98",
    "LTR/Pao":         "#BB8FCE",
    "LTR/ERVK":        "#884EA0",
    "LINE/L1":         "#2471A3",
    "LINE/L1-Tx1":     "#2E86C1",
    "LINE/L2":         "#85C1E9",
    "LINE/CR1":        "#5DADE2",
    "LINE/I":          "#1ABC9C",
    "LINE/Penelope":   "#48C9B0",
    "LINE/RTE-BovB":   "#76D7C4",
    "DNA/TcMar-Tc1":   "#D4AC0D",
    "DNA/TcMar-ISRm11":"#F1C40F",
    "DNA/hAT-Tip100":  "#F39C12",
    "DNA/Maverick":    "#E67E22",
    "DNA/Zator":       "#CA6F1E",
    "DNA/Kolobok-T2":  "#784212",
    "DNA/hAT-Ac":      "#D35400",
    "DNA/Zisupton":    "#BA4A00",
    "DNA/CMC-EnSpm":   "#A04000",
    "DNA/MULE-MuDR":   "#7E5109",
    "DNA/Merlin":      "#6E2F1A",
    "DNA/PIF-Harbinger":"#935116",
    "DNA/P":           "#7B7D7D",
    "DNA/Crypton-V":   "#626567",
    "RC/Helitron":     "#17A589",
    "PLE":             "#0E6655",
    "PLE/Naiad":       "#0B5345",
    "SINE/tRNA":       "#F1948A",
    "SINE/tRNA-Core":  "#EC7063",
    "Simple_repeat":   "#D5D8DC",
    "Low_complexity":  "#BFC9CA",
    "Unknown":         "#ABB2B9",
}

def te_color(cls):
    return _TE_PALETTE.get(_canonical_class(cls), "#bdbdbd")

# Collapse TE class to main family for legend
def te_main_class(cls):
    if cls in ("Simple_repeat","Low_complexity"):
        return "Simple/Low-complexity"
    p = cls.split("/")[0]
    return p

# ─── PANEL A: Bar chart — TE% per species ─────────────────────────────────────
def panel_a(ax):
    sp_list = identified["species"].tolist()
    te_pct  = identified["TE_pct_flank10"].fillna(0).tolist()
    mt_list = identified["mating_type"].tolist()
    colors  = [MT_COLOR[m] for m in mt_list]

    y = np.arange(len(sp_list))
    ax.barh(y, te_pct, color=colors, height=0.7, linewidth=0.3, edgecolor="white")
    ax.set_yticks(y)
    ax.set_yticklabels([f"$\\it{{{s}}}$" for s in sp_list], fontsize=6.5)
    ax.set_xlabel(f"TE content within \u00b1{FLANK_BP//1000} kb of sex gene (%)", fontsize=8)
    ax.set_xlim(0, max(75, max(te_pct) * 1.10 if te_pct else 75))
    ax.axvline(0, color="k", lw=0.5)
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)

    # Annotate bars with TE%
    for i, (v, sp, mt) in enumerate(zip(te_pct, sp_list, mt_list)):
        if v > 1:
            ax.text(v + 0.5, i, f"{v:.1f}%", va="center", fontsize=5.5, color="black")

    # Legend
    handles = [
        mpatches.Patch(color=MT_COLOR["P"], label="sexP confirmed"),
        mpatches.Patch(color=MT_COLOR["M"], label="sexM confirmed"),
    ]
    ax.legend(handles=handles, fontsize=6.5, loc="lower right",
              frameon=True, framealpha=0.9, edgecolor="lightgray")
    ax.set_title(f"A   TE density (\u00b1{FLANK_BP//1000} kb)",
                 fontsize=9, loc="left", pad=4, fontweight="bold")


# ─── Helper: draw gene arrow ──────────────────────────────────────────────────
def draw_gene(ax, x0, x1, y, color, label=None, lw=1.5, arrow_frac=0.12, fontsize=6):
    """Draw a thick arrow representing a gene. x0/x1 in data coords."""
    width = 0.5
    hw    = 0.7
    hl    = (x1 - x0) * arrow_frac
    hl    = min(hl, abs(x1 - x0) * 0.4)
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                 arrowprops=dict(arrowstyle=f"-|>", color=color,
                                 lw=lw,
                                 mutation_scale=10,
                                 shrinkA=0, shrinkB=0))
    ax.plot([x0, x1], [y, y], color=color, lw=3.5, solid_capstyle="butt", zorder=3)
    if label:
        xm = (x0 + x1) / 2
        ax.text(xm, y + 0.45, label, ha="center", va="bottom", fontsize=fontsize,
                fontstyle="italic", color=color)


def draw_arrow_gene(ax, start, end, y, color, label=None, fontsize=6.5, height=0.55):
    """Draw a filled arrow gene model."""
    from matplotlib.patches import FancyArrow
    dx = end - start
    ax.add_patch(FancyArrow(start, y, dx, 0,
                             width=height, head_width=height*1.5,
                             head_length=abs(dx)*0.15,
                             length_includes_head=True,
                             color=color, zorder=4, linewidth=0))
    if label:
        xm = (start + end) / 2
        ax.text(xm, y + height * 1.3, label, ha="center", va="bottom",
                fontsize=fontsize, fontstyle="italic", color=color, clip_on=True)


# ─── PANEL B: Synteny locus diagrams ─────────────────────────────────────────
# User request 2026-04-22: panel B must show the SAME 10 species as panel A,
# i.e. the top-10 confirmed sexP/sexM species ranked by TE % in the
# idiomorphic region. Order = highest TE% on top → lowest at the bottom.
REPRESENTATIVES = [(sp, "top10") for sp in identified.sort_values(
    "TE_pct_flank10", ascending=False)["species"].tolist()]

def panel_b(ax_list, te_detail):
    """Draw one locus diagram per axis."""
    for idx, (sp, stype) in enumerate(REPRESENTATIVES):
        ax = ax_list[idx]
        rows_sp = master[master["species"] == sp]
        if len(rows_sp) == 0 or rows_sp.iloc[0]["mating_type"] not in ("P", "M"):
            ax.set_visible(False)
            continue
        row = rows_sp.iloc[0]

        scaf = row["sex_scaffold"]
        if pd.isna(scaf):
            ax.set_visible(False)
            continue

        # Idio region (optional — only used for the light grey shading;
        # since 2026-05-01 the figure no longer depends on it for the metric)
        idio_row = idio[idio["species"]==sp]
        if len(idio_row) == 0:
            reg_start = reg_end = None
        else:
            # take last (avoid duplicate for Phybl4)
            _ir = idio_row.iloc[-1]
            reg_start = _ir["start"]
            reg_end   = _ir["end"]

        # Get rnhA + tptA from tblastn on locus scaffold
        hits = tblastn[(tblastn["species"]==sp) & (tblastn["scaffold"]==scaf)]

        # 2026-04-24: FIXED view window — same width for all species, centred
        # on the sexP/M gene midpoint, sized to the most spacious well-annotated
        # locus (PnitS607: tptA…sexP…helicase block ≈ ±15 kb). This guarantees
        # that all panels are drawn on the SAME scale and TE/gene density can
        # be compared visually 1-to-1.
        FIXED_VIEW_FLANK = 15000
        sex_mid = int((row["sex_start"] + row["sex_end"]) / 2)
        # 2026-05-01: do NOT clamp view_start to 0 — clamping breaks the
        # alignment of species whose sex locus lies near the start of a
        # scaffold (e.g. Lenpar1, sex_mid ≈ 4 kb).  The fixed-width window
        # must always be 2 × FIXED_VIEW_FLANK so all panel-B tracks share
        # the same scale; off-scaffold positions simply have no features
        # to draw, which is the correct behaviour.
        view_start = sex_mid - FIXED_VIEW_FLANK
        view_end   = sex_mid + FIXED_VIEW_FLANK

        # We still load helicase / tptA HSPs (only used to draw arrows below;
        # NO longer used to size the view).
        rnhA_hits_all = tblastn[(tblastn["species"]==sp) & (tblastn["scaffold"]==scaf) & (tblastn["marker"]=="rnhA")]
        tptA_hits_all = tblastn[(tblastn["species"]==sp) & (tblastn["scaffold"]==scaf) & (tblastn["marker"]=="tptA")]
        _anchor_lo = sex_mid - 50000
        _anchor_hi = sex_mid + 50000
        heli_all = user_markers[(user_markers["species"]==sp) &
                                 (user_markers["scaffold"]==scaf) &
                                 (user_markers["marker"]=="helicase") &
                                 (user_markers["end"]   >= _anchor_lo) &
                                 (user_markers["start"] <= _anchor_hi)]
        tptA_user_all = user_markers[(user_markers["species"]==sp) &
                                       (user_markers["scaffold"]==scaf) &
                                       (user_markers["marker"]=="tptA") &
                                       (user_markers["end"]   >= _anchor_lo) &
                                       (user_markers["start"] <= _anchor_hi)]

        # Normalize so view_start=0
        def norm(x):
            return x - view_start
        vlen = view_end - view_start

        ax.set_xlim(0, vlen)
        ax.set_ylim(-1.2, 2.2)
        ax.axis("off")

        # Backbone line
        ax.plot([0, vlen], [0, 0], color="#555555", lw=1.5, zorder=1)

        # Idiomorphic region shading (only if known)
        if reg_start is not None and reg_end is not None:
            ax.axvspan(norm(reg_start), norm(reg_end),
                       ymin=0.25, ymax=0.75, alpha=0.12, color="#888888", zorder=0)

        # TE blocks within the WIDER view window (±15 kb of sex midpoint).
        # We draw all TEs intersecting the view (te_view), and outline in red
        # any single annotation ≥ 3.5 kb (= "full-length-like" = potentially
        # autonomous copy).
        sp_te = te_view[(te_view["species"]==sp)]
        sp_fl = fl_te[fl_te["species"]==sp]
        fl_set = set(zip(sp_fl["te_start"], sp_fl["te_end"]))
        for _, tr in sp_te.iterrows():
            ts = norm(tr["te_start"])
            te = norm(tr["te_end"])
            w  = te - ts
            tc = te_color(tr["te_class"])
            is_fl = (tr["te_start"], tr["te_end"]) in fl_set
            ax.add_patch(mpatches.Rectangle((ts, -0.22), w, 0.44,
                                             color=tc, zorder=2, linewidth=0))
            if is_fl:
                # Outline + star marker on top of full-length TE
                ax.add_patch(mpatches.Rectangle((ts, -0.22), w, 0.44,
                                                 fill=False, edgecolor="black",
                                                 linewidth=0.8, zorder=3))
                ax.plot([(ts+te)/2], [0.35], marker="*", color="black",
                         markersize=7, zorder=6, markeredgewidth=0.3)

        # TE blocks in sex gene CDS
        sp_cds = te_detail[(te_detail["species"]==sp) & (te_detail["region"]=="sex_gene_CDS")]
        for _, tr in sp_cds.iterrows():
            ts = norm(tr["te_start"])
            te = norm(tr["te_end"])
            w  = te - ts
            tc = te_color(tr["te_class"])
            ax.add_patch(mpatches.Rectangle((ts, -0.22), w, 0.44,
                                             color=tc, alpha=0.7, zorder=5, linewidth=0.5,
                                             edgecolor="red"))

        # ── 2026-04-24: Background grey gene arrows from JGI/GFF annotations ──
        bg_genes = _load_gff_genes(sp, scaf, view_start, view_end)
        sex_s = row.get("sex_start"); sex_e = row.get("sex_end")
        for tup in bg_genes:
            gs, ge, gst = tup[0], tup[1], tup[2]
            pid = tup[3] if len(tup) > 3 else ""
            # Skip if it overlaps the sexP/M gene by >50% (it'll be drawn in colour)
            if pd.notna(sex_s) and pd.notna(sex_e):
                ov = max(0, min(ge, sex_e) - max(gs, sex_s))
                if ov > 0.5 * (ge - gs):
                    continue
            xs, xe = norm(gs), norm(ge)
            if gst == "-":
                xs, xe = xe, xs
            draw_arrow_gene(ax, xs, xe, 1.55, "#9aa0a6", label=None, height=0.32)
            # Tiny proteinId label below the arrow
            if pid:
                ax.text((xs + xe) / 2, 1.12, str(pid), ha="center", va="top",
                        fontsize=4.0, color="#666", clip_on=True)

        # ── 2026-04-24: Helicase HSP clusters (user-provided reference) ──
        if len(heli_all):
            # Cluster HSPs whose gap < 1.5 kb into one gene call
            hh = heli_all.sort_values("start").reset_index(drop=True)
            clusters = []
            cs, ce, css = hh.iloc[0]["start"], hh.iloc[0]["end"], hh.iloc[0]["strand"]
            for _, h in hh.iloc[1:].iterrows():
                if h["start"] - ce <= 1500:
                    ce = max(ce, h["end"])
                else:
                    clusters.append((cs, ce, css))
                    cs, ce, css = h["start"], h["end"], h["strand"]
            clusters.append((cs, ce, css))
            for (cs, ce, css) in clusters:
                xs, xe = norm(cs), norm(ce)
                if css == "-":
                    xs, xe = xe, xs
                draw_arrow_gene(ax, xs, xe, 0.8, "#8E44AD", label="helicase",
                                fontsize=5.5, height=0.45)

        # ── 2026-04-24: Locus orthologs (v3) — sagA / glrA / noname ──
        # Drawn from locus_orthologs_v3.tsv (4-criterion BBH-like filter:
        # E<1e-15, alen≥50% qlen, pident≥25%, ≤50 kb of sex midpoint).
        if len(locus_orthos):
            ORTH_COLORS = {"sagA":   "#16A085",
                           "glrA":   "#E74C3C",
                           "noname": "#9B59B6"}
            extra = locus_orthos[(locus_orthos["species"]==sp) &
                                  (locus_orthos["scaffold"]==scaf) &
                                  (locus_orthos["marker"].isin(ORTH_COLORS)) &
                                  (locus_orthos["end"]   >= view_start) &
                                  (locus_orthos["start"] <= view_end)]
            for _, eg in extra.iterrows():
                xs, xe = norm(eg["start"]), norm(eg["end"])
                if eg["strand"] == "-":
                    xs, xe = xe, xs
                draw_arrow_gene(ax, xs, xe, 0.8,
                                ORTH_COLORS[eg["marker"]],
                                label=eg["marker"], fontsize=5.5, height=0.45)

        # tptA gene (green)
        tptA_h = hits[hits["marker"]=="tptA"].sort_values("bitscore", ascending=False)
        if len(tptA_h) > 0:
            t = tptA_h.iloc[0]
            ts, te = norm(t["start"]), norm(t["end"])
            if t["strand"] == "-":
                ts, te = te, ts
            draw_arrow_gene(ax, ts, te, 0.8, "#27AE60", label="tptA", fontsize=5.5, height=0.45)

        # rnhA gene (blue)
        rnhA_h = hits[hits["marker"]=="rnhA"].sort_values("bitscore", ascending=False)
        if len(rnhA_h) > 0:
            r = rnhA_h.iloc[0]
            rs, re = norm(r["start"]), norm(r["end"])
            if r["strand"] == "-":
                rs, re = re, rs
            draw_arrow_gene(ax, rs, re, 0.8, "#2980B9", label="rnhA", fontsize=5.5, height=0.45)

        # sex gene (orange=P, red=M)
        mt = row["mating_type"]
        sc = "#C0392B" if "M" in mt else "#E8671A"
        sx_start = norm(row["sex_start"])
        sx_end   = norm(row["sex_end"])
        if row["sex_strand"] == "-":
            sx_start, sx_end = sx_end, sx_start
        draw_arrow_gene(ax, sx_start, sx_end, 0.8, sc, label="sex", fontsize=5.5, height=0.45)

        # Scale bar  = 1kb
        sb_x = vlen - 1200
        ax.plot([sb_x, sb_x + 1000], [-0.85, -0.85], color="black", lw=1.5)
        ax.text(sb_x + 500, -1.05, "1 kb", ha="center", va="top", fontsize=5.5, clip_on=False)

        # TE% label (±10 kb flanking window — same metric as panel A)
        te_pct = row.get("TE_pct_flank10")
        if pd.notna(te_pct):
            # place label below the sex-gene midpoint
            ax.text(norm(sex_mid), -0.55,
                    f"{te_pct:.1f}% TE", ha="center", va="top", fontsize=5.5,
                    color="#555555", clip_on=True)

        # Species label — in axes fraction coords so it never overlaps
        mt_label = {"P":"P","M":"M","P_tentative":"P?","M_tentative":"M?",
                    "not_identified":"—"}.get(mt, mt)
        ax.text(-0.01, 0.5, f"$\\it{{{sp}}}$ ({mt_label})",
                ha="right", va="center", fontsize=6.5, transform=ax.transAxes)

    # Gene legend on first axis
    ax0 = ax_list[0]
    leg_patches = [
        mpatches.Patch(color="#E8671A", label="sexP gene"),
        mpatches.Patch(color="#C0392B", label="sexM gene"),
        mpatches.Patch(color="#27AE60", label="tptA"),
        mpatches.Patch(color="#2980B9", label="rnhA"),
        mpatches.Patch(color="#8E44AD", label="helicase"),
        mpatches.Patch(color="#16A085", label="sagA"),
        mpatches.Patch(color="#E74C3C", label="glrA"),
        mpatches.Patch(color="#9B59B6", label="noname"),
        mpatches.Patch(color="#9aa0a6", label="adjacent gene (JGI)"),
        mpatches.Patch(color="#ABB2B9", label="TE (idiomorphic)"),
    ]
    ax0.legend(handles=leg_patches, fontsize=5.5, loc="upper right",
               frameon=True, framealpha=0.9, edgecolor="lightgray", ncol=10,
               bbox_to_anchor=(1.0, 2.5))

    # Title
    ax_list[0].set_title("B   Locus synteny (representative species)", fontsize=9,
                          loc="left", pad=28, fontweight="bold")


# ─── PANEL C: TE class composition — stacked bar ──────────────────────────────
def panel_c(ax):
    """Stacked bar chart of TE bp per class within \u00b1FLANK_BP of sex-gene midpoint, per species."""
    sp_order = identified["species"].tolist()  # same order as panel A

    # Aggregate TE bp by species and main class, restricted to the same
    # \u00b110 kb flanking window used for panel A (consistency requested
    # by user 2026-05-01: drop the heterogeneous "idiomorphic region").
    sp_te = te_view.copy()
    # Clip TE intervals to per-species \u00b1FLANK_BP window around sex_mid
    sex_mid_map = {}
    for _, _r in master.iterrows():
        if pd.notna(_r.get("sex_start")) and pd.notna(_r.get("sex_end")):
            sex_mid_map[_r["species"]] = (float(_r["sex_start"]) + float(_r["sex_end"])) / 2.0
    sp_te["_sex_mid"] = sp_te["species"].map(sex_mid_map)
    sp_te = sp_te.dropna(subset=["_sex_mid"]).copy()
    sp_te["_ws"] = sp_te["_sex_mid"] - FLANK_BP
    sp_te["_we"] = sp_te["_sex_mid"] + FLANK_BP
    s = sp_te[["te_start","_ws"]].max(axis=1)
    e = sp_te[["te_end",  "_we"]].min(axis=1)
    sp_te["_ovl_flank"] = (e - s + 1).clip(lower=0)
    sp_te = sp_te[sp_te["_ovl_flank"] > 0]
    sp_te["main_class"] = sp_te["te_class"].apply(te_main_class)
    agg = sp_te.groupby(["species","main_class"])["_ovl_flank"].sum().unstack(fill_value=0)

    # Add missing species as zeros
    for sp in sp_order:
        if sp not in agg.index:
            agg.loc[sp] = 0
    agg = agg.loc[sp_order]

    # Collapse small classes into "Other"
    class_totals = agg.sum()
    major = class_totals[class_totals > 200].index.tolist()
    other = [c for c in agg.columns if c not in major]
    if other:
        agg["Other"] = agg[other].sum(axis=1)
        agg = agg.drop(columns=other)
    classes = agg.columns.tolist()

    # Color map for main classes (cladogram TE palette)
    MAIN_COLOR = {
        "LTR":                 _TE_PALETTE["LTR"],
        "LINE":                _TE_PALETTE["LINE"],
        "DNA":                 _TE_PALETTE["DNA"],
        "RC":                  _TE_PALETTE["Rolling Circle"],
        "PLE":                 _TE_PALETTE["Penelope"],
        "Penelope":            _TE_PALETTE["Penelope"],
        "SINE":                _TE_PALETTE["SINE"],
        "Simple/Low-complexity":_TE_PALETTE["Unclassified"],
        "Unknown":             _TE_PALETTE["Unclassified"],
        "Unclassified":        _TE_PALETTE["Unclassified"],
        "Other":               "#cccccc",
    }
    colors = [MAIN_COLOR.get(c, "#CCCCCC") for c in classes]

    y   = np.arange(len(sp_order))
    bot = np.zeros(len(sp_order))
    for cls, col in zip(classes, colors):
        vals = agg[cls].values
        ax.barh(y, vals, left=bot, color=col, height=0.7, label=cls, linewidth=0)
        bot += vals

    ax.set_yticks(y)
    ax.set_yticklabels([])  # already labeled in panel A
    ax.set_xlabel(f"TE bp within \u00b1{FLANK_BP//1000} kb of sex gene", fontsize=8)
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)

    handles = [mpatches.Patch(color=MAIN_COLOR.get(c,"#CCC"), label=c) for c in classes]
    ax.legend(handles=handles, fontsize=6, loc="lower right",
              frameon=True, framealpha=0.9, edgecolor="lightgray", ncol=1)
    ax.set_title("C   TE class composition", fontsize=9, loc="left", pad=4, fontweight="bold")


# ─── Assemble figure ──────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
fig.patch.set_facecolor("white")

n_sp = len(identified)  # 37

# Layout: left col = panels A+C (stacked), right col = panel B (locus diagrams)
gs_outer = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1.1, 1.6], wspace=0.35)

# Left: A on top, C on bottom
gs_left = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[0], wspace=0.05)
ax_a = fig.add_subplot(gs_left[0])
ax_c = fig.add_subplot(gs_left[1])

# Right: one subplot per representative species (12 rows)
n_rep = len(REPRESENTATIVES)
gs_right = gridspec.GridSpecFromSubplotSpec(n_rep, 1, subplot_spec=gs_outer[1], hspace=0.05)
ax_b_list = [fig.add_subplot(gs_right[i]) for i in range(n_rep)]

# Draw panels
panel_a(ax_a)
panel_c(ax_c)
panel_b(ax_b_list, te_det)

# Add figure title
fig.suptitle("Sex determination locus structure and TE content across Mucorales",
             fontsize=11, fontweight="bold", y=0.99)

outfile = WORKDIR / "sex_locus_figure.pdf"
fig.savefig(outfile, dpi=200, bbox_inches="tight")
print(f"Saved: {outfile}")

outfile_png = WORKDIR / "sex_locus_figure.png"
fig.savefig(outfile_png, dpi=150, bbox_inches="tight")
print(f"Saved: {outfile_png}")

# ── Also publish to the visible MAT_locus/figures tree (= /Volumes/Temp/...)
FIGURES_DIR = Path("/legserv/Temp/Carlos/EDF-TE/results/MAT_locus/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
for _ext in ("pdf", "png"):
    _dst = FIGURES_DIR / f"sex_locus_figure.{_ext}"
    fig.savefig(_dst, dpi=(200 if _ext == "pdf" else 150), bbox_inches="tight")
    print(f"Saved: {_dst}")
plt.close()
