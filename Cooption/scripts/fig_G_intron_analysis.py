#!/usr/bin/env python3
"""
fig_G_intron_coopted_vs_noncoopted.py
Intron count comparison: co-opted vs non-co-opted genes, by clade.

Two approaches:
  A) Per-gene: join intron_counts.txt gene IDs with classification gene IDs.
     Works where gene IDs are numeric (JGI) or transcript IDs match.
  B) Genome-level fallback: use Avg_Introns_All from D01 table per species,
     combined with the fraction of co-opted genes.

Output:
  figG_introns_pergene_by_clade.pdf/.png   (if per-gene data available for ≥3 clades)
  figG_introns_species_level.pdf/.png      (genome-average fallback / context plot)
"""
import os, glob, collections
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu

# ─── Paths ──────────────────────────────────────────────────────────────────
INTRON_DIR    = "/legserv/Temp/Carlos/EDF-TE/results/GenomicFeatures_V2"
CLASS_GLOB    = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/per_species/classification/*.true_TE_coopted_genes.tsv"
ALL_CLASS_GLOB = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/per_species/classification/*.gene_classes.tsv"
TRAITS_FILE   = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures/14_cooption_traits.tsv"
D01_FILE      = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/D01_te_genomics_merged_provisional.tsv"
CLADE_FILE    = "/legserv/Temp/Carlos/EDF-TE/results/HTT/clade_assignment.tsv"
BED_DIR       = "/legserv/Temp/Carlos/EDF-TE/results/Synteny/bed"
OUT_DIR       = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"

CLADE_LABELS = {
    "AP":"Aphelida","AS":"Ascomycota","BA":"Basidiomycota",
    "BL":"Blastocladiomycota","BS":"Basidiobolomycota","CH":"Chytridiomycota",
    "CR":"Cryptomycota","EN":"Entomophthoromycota","GL":"Glomeromycota",
    "HO":"Holozoa","KX":"Kickxellomycotina","MI":"Microsporidia",
    "MO":"Mortierellomycota","MU":"Mucoromycotina","NU":"Nuclearida",
    "OL":"Olpidiomycota","ZO":"Zoopagomycota",
}

CLADE_COLORS = {
    "AS":"#e8887a","BA":"#4472c4","MU":"#e07840","MO":"#f4a8b0",
    "GL":"#8b5e48","KX":"#f0a820","ZO":"#c8dc70","EN":"#e01010",
    "BS":"#d0b4e0","OL":"#3668b0","CH":"#e0d040","NE":"#9060b0",
    "BL":"#7030a0","MI":"#3aab98","CR":"#6aafe6","AP":"#5c9e3a",
    "NU":"#1a2850","HO":"#b4d0e8",
}

# ─── Load clade info ──────────────────────────────────────────────────────────
clades = pd.read_csv(CLADE_FILE, sep="\t", header=None, names=["species","clade"])
sp2clade = dict(zip(clades["species"], clades["clade"]))

# ─── Attempt per-gene intron matching ────────────────────────────────────────
# Strategy:
# 1. Direct join on gene_id
# 2. "gene_" prefix additions
# 3. id_map bridge: class_id → id_map(protein_id/JGI-numeric → gene_lookup_key)
#                   → match gene_lookup_key in intron_df (with/without "gene-" prefix)

per_gene_results = []  # list of dicts: {species, clade, group, median_introns, n}
species_done = set()
species_method = {}  # for diagnostics

print("[INFO] Attempting per-gene intron matching...")

def _build_idmap_lookup(sp):
    """Load id_map and return (direct_map, jgi_map).
    Filters to rows where matched=='yes' or lookup_key looks like a real locus ID."""
    idmap_file = os.path.join(BED_DIR, f"{sp}.id_map.tsv")
    if not os.path.exists(idmap_file):
        return {}, {}
    idmap = pd.read_csv(idmap_file, sep="\t", dtype=str)
    # Keep only matched=yes rows if that column exists
    if "matched" in idmap.columns:
        idmap = idmap[idmap["matched"] == "yes"]
    elif len(idmap.columns) >= 2:
        idmap = idmap.iloc[:, [0, 1]]
        idmap.columns = ["protein_id", "gene_lookup_key"]
    if "protein_id" not in idmap.columns or "gene_lookup_key" not in idmap.columns:
        if len(idmap.columns) >= 2:
            idmap = idmap.iloc[:, [0, 1]]
            idmap.columns = ["protein_id", "gene_lookup_key"]
        else:
            return {}, {}
    idmap = idmap[idmap["gene_lookup_key"].notna()]
    idmap = idmap[~idmap["gene_lookup_key"].isin(["protein", "gene_lookup_key", ""])]
    direct = dict(zip(idmap["protein_id"].str.strip(), idmap["gene_lookup_key"].str.strip()))
    # JGI format: extract pipe-separated numeric field (position [2])
    jgi = {}
    for pid, lk in direct.items():
        parts = pid.split("|")
        if len(parts) >= 3:
            jgi[parts[2]] = lk
    return direct, jgi

# Load all class files
class_files = sorted(glob.glob(ALL_CLASS_GLOB))
for cf in class_files:
    sp = os.path.basename(cf).split(".")[0]
    clade = sp2clade.get(sp, sp.split("_")[0])

    # Find intron file
    intron_file = os.path.join(INTRON_DIR, f"{sp}_ALL_intron_counts.txt")
    # Handle _c suffix mismatch (e.g. NU_Parvatl_c → NU_Parvatl)
    if not os.path.exists(intron_file):
        sp_base = sp.replace("_c", "")
        intron_file = os.path.join(INTRON_DIR, f"{sp_base}_ALL_intron_counts.txt")
    if not os.path.exists(intron_file):
        continue

    try:
        intron_raw = pd.read_csv(intron_file, sep="\t", dtype=str)
        intron_raw.columns = ["gene_id", "n_introns"]
        intron_raw["n_introns"] = pd.to_numeric(intron_raw["n_introns"], errors="coerce")
        intron_raw = intron_raw.dropna(subset=["n_introns"])
        intron_raw["n_introns"] = intron_raw["n_introns"].astype(int)
        # Normalised key: strip leading "gene-" or "gene:" for lookup
        intron_raw["key_norm"] = (intron_raw["gene_id"]
                                  .str.replace(r"^gene[-:]", "", regex=True))
        # Build fast lookup dicts
        intron_by_id   = dict(zip(intron_raw["gene_id"],   intron_raw["n_introns"]))
        intron_by_norm = dict(zip(intron_raw["key_norm"],  intron_raw["n_introns"]))

        class_df = pd.read_csv(cf, sep="\t", dtype=str)
        class_df = class_df.astype(object)  # avoid pandas StringDtype NaN issues
        class_df["gene_id"] = class_df["gene_id"].astype(str).str.strip()

        # ── Attempt 1: direct join ────────────────────────────────────────
        merged = class_df.merge(
            intron_raw[["gene_id","n_introns"]], on="gene_id", how="inner")
        method = "direct"

        # ── Attempt 2: "gene_" prefix variations ─────────────────────────
        if len(merged) < 10:
            intron_g = intron_raw.copy()
            intron_g["gene_id"] = "gene_" + intron_g["key_norm"]
            class_g  = class_df.copy()
            class_g["gene_id"] = "gene_" + class_g["gene_id"].str.replace(r"^gene_", "", regex=True)
            merged = class_g.merge(intron_g[["gene_id","n_introns"]], on="gene_id", how="inner")
            if len(merged) >= 10:
                method = "gene_prefix"

        # ── Attempt 3: id_map bridge ──────────────────────────────────────
        if len(merged) < 10:
            # Try both sp and sp with _c stripped (e.g. NU_Parvatl_c → NU_Parvatl)
            sp_for_idmap = sp if os.path.exists(os.path.join(BED_DIR, f"{sp}.id_map.tsv")) \
                           else sp.replace("_c", "")
            direct_map, jgi_map = _build_idmap_lookup(sp_for_idmap)
            if direct_map or jgi_map:
                class_df2 = class_df.copy()
                # Convert to object dtype to avoid StringDtype NaN assignment issues
                class_df2 = class_df2.astype(object)
                # Map class gene_id → gene_lookup_key
                class_df2["lookup_key"] = class_df2["gene_id"].map(direct_map)
                # Fill with JGI numeric mapping where direct fails
                mapped_jgi = class_df2["gene_id"].map(jgi_map)
                class_df2["lookup_key"] = class_df2["lookup_key"].combine_first(mapped_jgi)
                class_df2 = class_df2.dropna(subset=["lookup_key"])
                if len(class_df2) >= 10:
                    # Intron lookup: try lookup_key directly, then with "gene-" prefix
                    def _lookup_n(k):
                        k = str(k)
                        if k in intron_by_id:
                            return intron_by_id[k]
                        if k in intron_by_norm:
                            return intron_by_norm[k]
                        if "gene-" + k in intron_by_id:
                            return intron_by_id["gene-" + k]
                        return np.nan
                    class_df2["n_introns"] = class_df2["lookup_key"].apply(_lookup_n)
                    merged2 = class_df2.dropna(subset=["n_introns"]).copy()
                    merged2["n_introns"] = merged2["n_introns"].astype(int)
                    if len(merged2) >= 10:
                        merged = merged2
                        method = "idmap"

        if len(merged) < 10:
            continue

        n_coopted    = merged[merged["class"].isin(
            ["TE_coopted_candidate", "true_coopted"])]["n_introns"].values
        n_noncoopted = merged[merged["class"] == "non_TE"]["n_introns"].values
        if len(n_coopted) == 0:
            n_coopted = merged[~merged["class"].isin(["non_TE", "host_only"])]["n_introns"].values

        if len(n_coopted) < 5 or len(n_noncoopted) < 5:
            continue

        per_gene_results.append({"species": sp, "clade": clade, "group": "co-opted",
                                  "n": len(n_coopted),
                                  "median_introns": float(np.median(n_coopted)),
                                  "values": n_coopted})
        per_gene_results.append({"species": sp, "clade": clade, "group": "non-co-opted",
                                  "n": len(n_noncoopted),
                                  "median_introns": float(np.median(n_noncoopted)),
                                  "values": n_noncoopted})
        species_done.add(sp)
        species_method[sp] = method

    except Exception as e:
        pass

print(f"[INFO] Per-gene match succeeded for {len(species_done)} species")
from collections import Counter
method_counts = Counter(species_method.values())
print(f"[INFO] Methods: {dict(method_counts)}")

# ─── Per-gene figure (if enough data) ────────────────────────────────────────
if len(species_done) >= 5:
    # Build clade-level pooled arrays from raw per-gene values
    # per_gene_results entries have keys: species, clade, group, n, median_introns, values
    from collections import defaultdict
    clade_vals = defaultdict(lambda: {"co-opted": [], "non-co-opted": []})
    clade_species_count = defaultdict(set)

    for r in per_gene_results:
        cl = r["clade"]
        if cl == "MI":
            continue
        clade_vals[cl][r["group"]].extend(r["values"].tolist())
        clade_species_count[cl].add(r["species"])

    # Keep clades that have both groups with data
    valid_clades = [cl for cl in clade_vals
                    if len(clade_vals[cl]["co-opted"]) >= 5
                    and len(clade_vals[cl]["non-co-opted"]) >= 5]

    print(f"[INFO] Valid clades for per-gene plot: {valid_clades}")

    # df_pg still needed for stats TSV
    rows = []
    for r in per_gene_results:
        rows.append({"species": r["species"], "clade": r["clade"],
                     "group": r["group"], "median_introns": r["median_introns"]})
    df_pg = pd.DataFrame(rows)

    # Order clades by median co-opted introns (descending)
    clade_diffs = {}
    for cl in valid_clades:
        m_co  = float(np.median(clade_vals[cl]["co-opted"]))
        m_non = float(np.median(clade_vals[cl]["non-co-opted"]))
        clade_diffs[cl] = m_co - m_non
    clade_order = sorted(valid_clades, key=lambda c: clade_diffs.get(c, 0), reverse=True)

    fig, ax = plt.subplots(figsize=(max(10, len(clade_order)*1.4), 6))
    width = 0.35

    for i, cl in enumerate(clade_order):
        color = CLADE_COLORS.get(cl, "#888888")
        vals_co  = clade_vals[cl]["co-opted"]
        vals_non = clade_vals[cl]["non-co-opted"]
        n_sp = len(clade_species_count[cl])

        bp1 = ax.boxplot(vals_co,  positions=[i - width/2], widths=0.28,
                         patch_artist=True, boxprops=dict(facecolor=color, alpha=0.85),
                         medianprops=dict(color="black", linewidth=1.8),
                         whiskerprops=dict(linestyle="--"),
                         showfliers=False)
        bp2 = ax.boxplot(vals_non, positions=[i + width/2], widths=0.28,
                         patch_artist=True, boxprops=dict(facecolor=color, alpha=0.30),
                         medianprops=dict(color="black", linewidth=1.8),
                         whiskerprops=dict(linestyle="--"),
                         showfliers=False)

        # Significance annotation
        try:
            U, p = mannwhitneyu(vals_co, vals_non, alternative="two-sided")
            sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        except Exception:
            sig = ""
        y_top = max(np.percentile(vals_co, 75), np.percentile(vals_non, 75)) + 0.5
        ax.text(i, y_top, sig, ha="center", va="bottom", fontsize=8)
        ax.text(i, -0.8, f"n={n_sp}", ha="center", va="top", fontsize=7, color="grey")

    xlabels = [CLADE_LABELS.get(cl, cl) for cl in clade_order]
    ax.set_xticks(range(len(clade_order)))
    ax.set_xticklabels(xlabels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Intron count per gene (all genes pooled)", fontsize=10)
    ax.set_title("Intron number: co-opted vs non-co-opted genes per clade\n(pooled across species; n = number of species)", fontsize=10)
    ax.set_ylim(bottom=-1.2)

    co_patch  = mpatches.Patch(facecolor="grey", alpha=0.85, label="Co-opted genes")
    non_patch = mpatches.Patch(facecolor="grey", alpha=0.30, label="Non-co-opted genes")
    ax.legend(handles=[co_patch, non_patch], fontsize=9)
    plt.tight_layout()
    fig.savefig(f"{OUT_DIR}/figG_introns_pergene_by_clade.pdf", dpi=150, bbox_inches="tight")
    fig.savefig(f"{OUT_DIR}/figG_introns_pergene_by_clade.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] figG_introns_pergene_by_clade.pdf/.png")

    # Save stats
    stat_rows = []
    for sp in species_done:
        co_vals  = next((r["values"] for r in per_gene_results if r["species"]==sp and r["group"]=="co-opted"), None)
        non_vals = next((r["values"] for r in per_gene_results if r["species"]==sp and r["group"]=="non-co-opted"), None)
        if co_vals is not None and non_vals is not None and len(co_vals) >= 3 and len(non_vals) >= 3:
            try:
                U, p = mannwhitneyu(co_vals, non_vals, alternative="two-sided")
            except:
                U, p = np.nan, np.nan
            clade = sp2clade.get(sp, sp.split("_")[0])
            stat_rows.append({
                "species": sp, "clade": clade,
                "n_coopted": len(co_vals),
                "n_noncoopted": len(non_vals),
                "median_coopted": float(np.median(co_vals)),
                "median_noncoopted": float(np.median(non_vals)),
                "mw_U": U, "mw_p": p,
                "signif": p < 0.05 if not np.isnan(p) else False
            })
    pd.DataFrame(stat_rows).sort_values("mw_p").to_csv(
        f"{OUT_DIR}/figG_intron_stats.tsv", sep="\t", index=False)
    print(f"[OK] figG_intron_stats.tsv")
else:
    print(f"[WARN] Only {len(species_done)} species matched — skipping per-gene figure")

# ─── Species-level figure using D01 genome averages + traits ─────────────────
print("[INFO] Building species-level intron context figure...")
d01 = pd.read_csv(D01_FILE, sep="\t")
traits = pd.read_csv(TRAITS_FILE, sep="\t")

df_sp = traits.merge(d01[["species","Avg_Introns_All"]], on="species", how="left")
df_sp["clade"] = df_sp["species"].str.split("_").str[0]
df_sp["clade_label"] = df_sp["clade"].map(CLADE_LABELS).fillna(df_sp["clade"])
df_sp = df_sp.merge(clades, on="species", how="left")
df_sp["clade"] = df_sp["clade_y"].fillna(df_sp["clade_x"])
df_sp["clade_label"] = df_sp["clade"].map(CLADE_LABELS).fillna(df_sp["clade"])

clade_n2 = df_sp["clade"].value_counts()
df_sp2 = df_sp[df_sp["clade"].isin(clade_n2[clade_n2 >= 3].index)].dropna(
    subset=["Avg_Introns_All","pct_coopted"])

fig2, ax2 = plt.subplots(figsize=(8, 6))
for code, grp in df_sp2.groupby("clade"):
    color = CLADE_COLORS.get(code, "#888888")
    label = CLADE_LABELS.get(code, code)
    ax2.scatter(grp["Avg_Introns_All"], grp["pct_coopted"],
                color=color, s=25, alpha=0.75, label=label)

ax2.set_xlabel("Genome-wide avg introns per gene (Avg_Introns_All)", fontsize=10)
ax2.set_ylabel("% co-opted genes per species", fontsize=10)
ax2.set_title("TE co-option level vs genome-wide intron count\n(one dot = one species)", fontsize=11)
ax2.legend(fontsize=7, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
plt.tight_layout()
fig2.savefig(f"{OUT_DIR}/figG_introns_species_level.pdf", dpi=150, bbox_inches="tight")
fig2.savefig(f"{OUT_DIR}/figG_introns_species_level.png", dpi=150, bbox_inches="tight")
plt.close(fig2)
print(f"[OK] figG_introns_species_level.pdf/.png")

print("[DONE] Fig G complete.")
