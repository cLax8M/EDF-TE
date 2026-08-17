#!/usr/bin/env python3
"""
fig_C_domain_barchart.py
Bar chart of most frequent TE-host domain pairs (panel d style).
For the top N TE superfamilies: how many genes + how many species carry each.
Colors:
  - DNA transposons (Transposase, CMC-EnSpm, DNA/hAT, Helitron, Helicase): blues/greens
  - Retrotransposons (Reverse_transcriptase, LTR*, Gypsy*, LINE*, Integrase, Gag, Copia): reds/oranges
  - Other: grey
"""
import sys
import os
import glob
import collections
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ─── Paths ──────────────────────────────────────────────────────────────────
COOPTION_DIR   = "/legserv/Temp/Carlos/EDF-TE/results/Co-option"
TE_MASTER      = f"{COOPTION_DIR}/config/TE_domains_master.tsv"
OUT_DIR        = f"{COOPTION_DIR}/figures"
FUSIONS_GLOB   = f"{COOPTION_DIR}/per_species/fusions/*.TE_host_fusions.tsv"
CLASS_GLOB     = f"{COOPTION_DIR}/per_species/classification/*.true_TE_coopted_genes.tsv"
ENRICH_FILE    = f"{OUT_DIR}/step6_domain_enrichment.tsv"

TOP_N          = 40   # number of domain pairs to show
SPECIES_PREFIX = 2    # code length (e.g. "AP_")

# ─── TE type classification ───────────────────────────────────────────────────
DNA_SF   = {"Transposase", "CMC-EnSpm", "DNA/hAT", "Helitron", "Helicase", "CMC-EnSpm"}
RT_SF    = {"Reverse_transcriptase", "LTR", "LTR/Copia", "LTR/Gypsy",
            "LINE", "LINE/LTR", "Gypsy/Chromovirus", "Integrase", "Gag", "Protease"}

def te_category(sf):
    if sf in DNA_SF:
        return "DNA transposon"
    elif sf in RT_SF:
        return "Retrotransposon"
    else:
        return "Other TE"

CATEGORY_COLORS = {
    "DNA transposon":  "#2166ac",
    "Retrotransposon": "#d6604d",
    "Other TE":        "#969696",
}

# ─── Load TE master ───────────────────────────────────────────────────────────
te_master = pd.read_csv(TE_MASTER, sep="\t")
te_core_ids = set(te_master.loc[te_master["TE_status"] == "TE_core", "domain_id"].astype(str))
te_id2sf    = dict(zip(te_master["domain_id"].astype(str),
                        te_master["TE_superfamily"].astype(str)))
print(f"[INFO] TE core domains: {len(te_core_ids)}")

# ─── Load true co-opted gene keys ──────────────────────────────────────────
true_genes = collections.defaultdict(set)  # species -> set of gene_ids
for f in sorted(glob.glob(CLASS_GLOB)):
    sp = os.path.basename(f).split(".")[0]
    df = pd.read_csv(f, sep="\t", usecols=["gene_id"])
    true_genes[sp].update(df["gene_id"].astype(str))

all_true_genes = {g for s in true_genes.values() for g in s}
print(f"[INFO] True co-opted genes loaded: {len(all_true_genes)} across {len(true_genes)} species")

# ─── Load step6 enrichment for non-TE domain names ───────────────────────────
enrich = pd.read_csv(ENRICH_FILE, sep="\t")
# Keep only enriched, non-TE domains
enrich_host = enrich[enrich["effect"] == "enriched"].copy()
# We'll use n_true_with_domain and domain name directly

# ─── Domain name canonicalization (merge related names) ──────────────────────
# Groups similar domains so "Helicase C-terminal", "helicase superfamily c-terminal", etc.
# map to a single canonical name.
DOMAIN_MERGE = [
    # pattern (case-insensitive substring) → canonical name
    # ─ Helicases ─────────────────────────────────────────────────────────────
    (["helicase superfamily c-terminal", "helicase conserved c-terminal",
      "helicase_c", "helicase c-term", "helicase associated domain",
      "ha2", "snf2", "hsa domain", "hsa"], "Helicase C-terminal domain"),
    (["helicase superfamily n-terminal", "helicase, superfamily i/ii",
      "helicase_n", "helicase n-term", "aaa+ atpase",
      "sfii helicase", "dead-like helicases superfamily",
      "dead_atp_helicase"], "Helicase N-terminal / ATPase"),
    (["dead/deah box helicase", "dead-box helicase", "dead box",
      "dead/deah", "dead_cter", "dead_nter", "dead box helicase",
      "deah-box", "dead-box rna helicase", "dead-box subfamily",
      "q_motif", "q motif"], "DEAD/DEAH box helicase"),
    # ─ Zinc fingers ──────────────────────────────────────────────────────────
    (["zinc finger, c2h2", "zinc finger c2h2", "c2hcznfinger",
      "zinc finger, c2h2-like"], "Zinc finger C2H2"),
    (["zinc finger, cchc", "zinc finger cchc", "znf_cchc"], "Zinc finger CCHC"),
    (["zinc knuckle", "znf_c2hc"], "Zinc knuckle"),
    (["zinc finger ring", "ring-type", "zinc ring",
      "ring finger", "ring domain", "znf_rbz",
      "zinc finger, rbz type"], "RING finger"),
    (["flywch"], "FLYWCH zinc finger"),
    # ─ RNA-binding / recognition ──────────────────────────────────────────────
    (["rrm", "rna recognition motif", "rna-recognition", "rna_bind",
      "rna-binding"], "RNA recognition motif"),
    # ─ TE-core (skip list) ───────────────────────────────────────────────────
    (["reverse transcriptase", "rvt_1", "rvt_2", "rvt_thumb"], "Reverse transcriptase"),
    (["integrase", "integrase core"], "Integrase (TE)"),
    (["transposase", "transposase_1", "transposase_2", "dde endonuclease"], "Transposase"),
    # ─ DNA-binding / chromatin ───────────────────────────────────────────────
    (["myb/sant", "myb domain", "myb-like"], "Myb/SANT DNA-binding"),
    (["chromo", "chromodomain", "chromatin organization modifier"], "Chromodomain"),
    (["at hook"], "AT-hook"),
    (["bromo", "bromodomain"], "Bromodomain"),
    (["phd finger", "plant homeodomain", "phd-finger"], "PHD finger"),
    (["pwwp"], "PWWP domain"),
    # ─ Repeats ───────────────────────────────────────────────────────────────
    (["wd40", "wd repeat", "wd_repeat"], "WD40 repeat"),
    (["leucine-rich repeat", "leucine rich repeat", "lrr",
      "c-terminal to leucine-rich"], "Leucine-rich repeat"),
    (["ankyrin", "ank repeat", "ank_repeat", "ank_rep"], "Ankyrin repeat"),
    (["fn3", "fibronectin type iii", "fibronectin type-iii"], "Fibronectin type III"),
    (["kelch", "kelch repeat"], "Kelch repeat"),
    (["armadillo", "arm repeat"], "Armadillo repeat"),
    (["tpr repeat", "tetratricopeptide"], "TPR repeat"),
    (["pentatricopeptide", "ppr repeat"], "PPR repeat"),
    # ─ Proteases ─────────────────────────────────────────────────────────────
    (["aspartyl protease", "acid protease", "asp_protease",
      "asp_prot_retrov", "retroviral-type"], "Aspartyl protease"),
    # ─ Ubiquitin system ──────────────────────────────────────────────────────
    (["ubiquitin", "ubq", "ubiquitin-like"], "Ubiquitin"),
    (["f-box", "f box"], "F-box domain"),
    # ─ Signalling ────────────────────────────────────────────────────────────
    (["kinase", "pkinase", "protein kinase"], "Protein kinase"),
    # Specific phosphatases BEFORE the generic pattern
    (["inositol polyphosphate", "inositol phosphatase"], "Inositol polyphosphate phosphatase (APE)"),
    (["calcineurin", "pp2c", "phosphoesterase"], "Calcineurin / PP2C phosphatase"),
    (["phosphatase"], "Phosphatase"),
    (["amp-binding", "amp binding", "adenylate", "adenyltransferase"], "AMP-binding enzyme"),
    # ─ Other functional ──────────────────────────────────────────────────────
    (["oligonucleotide/oligosaccharide", "ob-fold", "ob fold",
      "oligonucleotide-binding", "obfold", "duf1605"], "OB-fold"),
    (["dshct", "nuc185"], "DSHCT / NUC185 domain"),
    # inositol moved above to Signalling section (before generic phosphatase)
    (["cenp-b", "centromere protein b"], "CENP-B DNA-binding"),
    (["sec63", "brl domain"], "Sec63/Brl domain"),
]

# Suffixes to strip before pattern matching (InterPro/Pfam descriptors)
_STRIP_SUFFIXES = [
    " profile.", " signature.", " motif profile.", " motif signature.",
    " domain profile.", " domain signature.", " family profile.",
    " family signature.", " add an annotation", " domain superfamily.",
    " superfamily.", " subfamily.", ",  add an annotation",
    "; add an annotation",
]

def canonicalize_domain(name):
    nl = name.lower().strip()
    # Strip common InterPro/Pfam suffixes
    for sfx in _STRIP_SUFFIXES:
        if nl.endswith(sfx):
            nl = nl[: -len(sfx)].strip()
            break
    for patterns, canonical in DOMAIN_MERGE:
        for pat in patterns:
            if pat in nl:
                return canonical
    # Return the cleaned name (suffix stripped) if no merge found
    return name.strip().rstrip(".")

# ─── Scan fusion files for host-domain × TE-superfamily pairs ────────────────
# For each true co-opted gene, record what TE domains and host domains coexist
pair_gene_counter   = collections.Counter()   # (host_domain_name, te_sf) -> n_genes
pair_species_counter = collections.Counter()  # (host_domain_name, te_sf) -> n_species
pair_te_cat  = {}   # (host_domain, te_sf) -> category

print("[INFO] Scanning fusion files...")
for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    df = pd.read_csv(fpath, sep="\t", dtype=str).fillna("")

    # Keep only true co-opted genes
    true_sp = true_genes.get(sp, set())
    df_true = df[df["gene_id"].isin(true_sp)]

    # Separate TE-core domains from host domains per gene
    gene_te_sfs     = collections.defaultdict(set)
    gene_host_names = collections.defaultdict(set)

    for _, row in df_true.iterrows():
        gid    = str(row["gene_id"])
        did    = str(row["domain_id"])
        dname  = str(row["domain_name"]).strip() if str(row["domain_name"]) not in ("", "nan") else did
        dclass = str(row.get("domain_class", ""))

        if did in te_core_ids or dclass in ("TE_core", "TE_associated"):
            sf = te_id2sf.get(did, "Unknown")
            gene_te_sfs[gid].add(sf)
        else:
            # host domain — skip TE-derived names that sneak in as false host domains
            if dname and dname not in ("-", "nan", "", r"\N", "no description"):
                canon = canonicalize_domain(dname)
                if canon not in {"Reverse transcriptase", "Transposase", "Integrase (TE)"}:
                    gene_host_names[gid].add(canon)

    # Count pairs
    gene_pairs_counted = collections.defaultdict(set)
    for gid in set(gene_te_sfs) & set(gene_host_names):
        for te_sf in gene_te_sfs[gid]:
            for hd in gene_host_names[gid]:
                key = (hd, te_sf)
                pair_gene_counter[key] += 1
                gene_pairs_counted[key].add(sp)

    for key, sps in gene_pairs_counted.items():
        pair_species_counter[key] += len(sps)
        pair_te_cat[key] = te_category(key[1])

print(f"[INFO] Total unique host-domain × TE-superfam pairs: {len(pair_gene_counter)}")

# ─── Clade tracking per pair ──────────────────────────────────────────────────
# Load clade assignments
clade_df = pd.read_csv(f"{COOPTION_DIR}/../HTT/clade_assignment.tsv",
                       sep="\t", header=0, names=["species","clade"])
sp2clade = dict(zip(clade_df["species"], clade_df["clade"]))
NE_SP = {"CH_Anasp1","CH_Neocon1","CH_Nealan1","CH_Orpsp1","CH_Pecora1","CH_Caecom1","CH_Pirfi3"}
for s in NE_SP:
    sp2clade[s] = "NE"

CLADE_COLORS_PY = {
    "AS":"#e8887a","BA":"#4472c4","MU":"#e07840","MO":"#f4a8b0","GL":"#8b5e48",
    "KX":"#f0a820","ZO":"#c8dc70","EN":"#e01010","BS":"#d0b4e0","OL":"#3668b0",
    "CH":"#e0d040","NE":"#9060b0","BL":"#7030a0","MI":"#3aab98","CR":"#a8c8e8",
    "AP":"#5c9e3a","NU":"#1a2850","HO":"#b4d0e8",
}
CLADE_ORDER_PY = ["HO","NU","AP","CR","MI","BL","OL","CH","NE","BS",
                  "EN","ZO","KX","GL","MO","MU","BA","AS"]
FUNGI_CLADES = {"AS","BA","MU","MO","GL","KX","ZO","EN","BS","NE",
                "CH","OL","BL","MI","CR","AP"}
BASAL_FUNGI_CLADES = {"MU","MO","GL","KX","ZO","EN","BS","NE","CH","OL","BL","MI","CR","AP"}

# ─── Re-aggregate fusions with canonicalized names AND per-clade counts ───────
# We need to re-scan; accumulate clade-level species counts per pair
print("[INFO] Re-aggregating with canonical domain names and clade info...")

pair_gene_counter2   = collections.Counter()
pair_clade_sp_counter = collections.defaultdict(lambda: collections.Counter())  # pair -> clade -> n_species
pair_clade_gene_counter = collections.defaultdict(lambda: collections.Counter()) # pair -> clade -> n_genes

for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    clade = sp2clade.get(sp, sp[:2])
    df = pd.read_csv(fpath, sep="\t", dtype=str).fillna("")
    true_sp = true_genes.get(sp, set())
    df_true = df[df["gene_id"].isin(true_sp)]

    gene_te_sfs     = collections.defaultdict(set)
    gene_host_names = collections.defaultdict(set)
    for _, row in df_true.iterrows():
        gid   = str(row["gene_id"])
        did   = str(row["domain_id"])
        dname = str(row["domain_name"]).strip() if str(row["domain_name"]) not in ("", "nan") else did
        dclass = str(row.get("domain_class", ""))
        if did in te_core_ids or dclass in ("TE_core", "TE_associated"):
            gene_te_sfs[gid].add(te_id2sf.get(did, "Unknown"))
        else:
            _bad = {"", "-", "nan", r"\N", "\\N", "no description"}
            if dname.strip() not in _bad and len(dname.strip()) > 0:
                canon = canonicalize_domain(dname)
                # Skip TE-derived domain names appearing as false "host" domains
                # (happens when a TE domain ID is not in te_core_ids but has a TE name)
                if canon not in {"Reverse transcriptase", "Transposase", "Integrase (TE)"}:
                    gene_host_names[gid].add(canon)

    seen_pairs_sp = set()
    for gid in set(gene_te_sfs) & set(gene_host_names):
        for te_sf in gene_te_sfs[gid]:
            for hd in gene_host_names[gid]:
                key = (hd, te_sf)
                pair_gene_counter2[key] += 1
                pair_clade_gene_counter[key][clade] += 1
                if key not in seen_pairs_sp:
                    pair_clade_sp_counter[key][clade] += 1
                    seen_pairs_sp.add(key)

# ─── Build summary table ──────────────────────────────────────────────────────
rows = []
for (hd, te_sf), n_genes in pair_gene_counter2.items():
    n_sp = sum(pair_clade_sp_counter[(hd, te_sf)].values())
    rows.append({
        "host_domain": hd, "te_superfam": te_sf,
        "te_category": te_category(te_sf),
        "n_genes": n_genes, "n_species": n_sp,
    })

df_pairs = pd.DataFrame(rows).sort_values("n_genes", ascending=False)
_bad = {"", "-", "nan", r"\N", "\\N"}
df_pairs = df_pairs[df_pairs["host_domain"].str.strip().apply(lambda x: x not in _bad and len(x) > 0)]
df_pairs.to_csv(f"{OUT_DIR}/figC_domain_tepair_summary.tsv", sep="\t", index=False)
print(f"[INFO] Saved figC_domain_tepair_summary.tsv")

# ─── Helper: stacked clade bar chart ─────────────────────────────────────────
def plot_stacked_clade_panel(ax, df_top, value_col, xlabel, title, use_clades):
    """Draw horizontal stacked bars colored by clade fraction."""
    n = len(df_top)
    lefts = np.zeros(n)
    for clade in CLADE_ORDER_PY:
        if clade not in use_clades:
            continue
        col_c = f"clade_{clade}"
        if col_c not in df_top.columns:
            continue
        vals = df_top[col_c].values[::-1]
        ax.barh(range(n), vals, left=lefts, color=CLADE_COLORS_PY[clade],
                edgecolor="none", height=0.75, label=clade)
        lefts += vals
    ax.set_yticks(range(n))
    ax.set_yticklabels(df_top["label"].values[::-1], fontsize=7)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=9, loc="left")

def make_df_top(df_pairs, subset_clades=None, top_n=TOP_N):
    """Build top-N table with per-clade columns for species counts."""
    if subset_clades is not None:
        # Filter pairs: only genes/species from those clades
        tmp_rows = []
        for _, row in df_pairs.iterrows():
            key = (row["host_domain"], row["te_superfam"])
            n_g = sum(v for c, v in pair_clade_gene_counter[key].items() if c in subset_clades)
            n_s = sum(v for c, v in pair_clade_sp_counter[key].items() if c in subset_clades)
            if n_g > 0:
                tmp_rows.append({**row.to_dict(), "n_genes": n_g, "n_species": n_s})
        dft = pd.DataFrame(tmp_rows).sort_values("n_genes", ascending=False).head(top_n).copy()
    else:
        dft = df_pairs.head(top_n).copy()
    dft["label"] = dft["host_domain"].str[:38] + "  |  " + dft["te_superfam"]
    # Add per-clade species columns
    for clade in CLADE_ORDER_PY:
        dft[f"clade_{clade}"] = dft.apply(
            lambda r: pair_clade_sp_counter[(r["host_domain"], r["te_superfam"])].get(clade, 0), axis=1)
    return dft

# ─── Figure 1: All species ────────────────────────────────────────────────────
df_top_all = make_df_top(df_pairs, subset_clades=None)
fig, axes = plt.subplots(1, 2, figsize=(16, 10))

# Panel A: n genes (TE-class colored)
ax1 = axes[0]
df_top_all["color"] = df_top_all["te_category"].map(CATEGORY_COLORS)
ax1.barh(range(TOP_N), df_top_all["n_genes"].values[::-1],
         color=df_top_all["color"].values[::-1], edgecolor="white", linewidth=0.4)
ax1.set_yticks(range(TOP_N))
ax1.set_yticklabels(df_top_all["label"].values[::-1], fontsize=7)
ax1.set_xlabel("Number of co-opted genes", fontsize=9)
ax1.set_title("A   All species — by gene count", fontsize=9, loc="left")
legend_handles = [mpatches.Patch(facecolor=v, label=k) for k, v in CATEGORY_COLORS.items()]
ax1.legend(handles=legend_handles, fontsize=7, title="TE class", loc="lower right")

# Panel B: n species stacked by clade
plot_stacked_clade_panel(axes[1], df_top_all, "n_species",
    "Number of species", "B   All species — by clade", set(CLADE_ORDER_PY))
# Clade legend (compact)
handles_cl = [mpatches.Patch(facecolor=CLADE_COLORS_PY[c], label=c)
              for c in CLADE_ORDER_PY if c in set(CLADE_ORDER_PY)]
axes[1].legend(handles=handles_cl, fontsize=6, ncol=3, title="Clade",
               loc="lower right", framealpha=0.8)

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/figC_domain_barchart.pdf", dpi=150, bbox_inches="tight")
plt.savefig(f"{OUT_DIR}/figC_domain_barchart.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] figC_domain_barchart.pdf/.png saved (all species)")

# ─── Figure 2: Fungi only ─────────────────────────────────────────────────────
df_top_fungi = make_df_top(df_pairs, subset_clades=FUNGI_CLADES)
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 10))
df_top_fungi["color"] = df_top_fungi["te_category"].map(CATEGORY_COLORS)
axes2[0].barh(range(len(df_top_fungi)), df_top_fungi["n_genes"].values[::-1],
              color=df_top_fungi["color"].values[::-1], edgecolor="white", linewidth=0.4)
axes2[0].set_yticks(range(len(df_top_fungi)))
axes2[0].set_yticklabels(df_top_fungi["label"].values[::-1], fontsize=7)
axes2[0].set_xlabel("Number of co-opted genes", fontsize=9)
axes2[0].set_title("A   Fungi only — by gene count", fontsize=9, loc="left")
axes2[0].legend(handles=legend_handles, fontsize=7, title="TE class", loc="lower right")
plot_stacked_clade_panel(axes2[1], df_top_fungi, "n_species",
    "Number of species", "B   Fungi only — by clade", FUNGI_CLADES)
axes2[1].legend(handles=[mpatches.Patch(facecolor=CLADE_COLORS_PY[c], label=c)
                          for c in CLADE_ORDER_PY if c in FUNGI_CLADES],
                fontsize=6, ncol=3, title="Clade", loc="lower right", framealpha=0.8)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/figC_domain_barchart_fungi.pdf", dpi=150, bbox_inches="tight")
plt.savefig(f"{OUT_DIR}/figC_domain_barchart_fungi.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] figC_domain_barchart_fungi.pdf/.png saved (Fungi only)")

# ─── Figure 3: Basal Fungi only ───────────────────────────────────────────────
df_top_basal = make_df_top(df_pairs, subset_clades=BASAL_FUNGI_CLADES)
fig3, axes3 = plt.subplots(1, 2, figsize=(16, 10))
df_top_basal["color"] = df_top_basal["te_category"].map(CATEGORY_COLORS)
axes3[0].barh(range(len(df_top_basal)), df_top_basal["n_genes"].values[::-1],
              color=df_top_basal["color"].values[::-1], edgecolor="white", linewidth=0.4)
axes3[0].set_yticks(range(len(df_top_basal)))
axes3[0].set_yticklabels(df_top_basal["label"].values[::-1], fontsize=7)
axes3[0].set_xlabel("Number of co-opted genes", fontsize=9)
axes3[0].set_title("A   Basal Fungi only — by gene count", fontsize=9, loc="left")
axes3[0].legend(handles=legend_handles, fontsize=7, title="TE class", loc="lower right")
plot_stacked_clade_panel(axes3[1], df_top_basal, "n_species",
    "Number of species", "B   Basal Fungi only — by clade", BASAL_FUNGI_CLADES)
axes3[1].legend(handles=[mpatches.Patch(facecolor=CLADE_COLORS_PY[c], label=c)
                          for c in CLADE_ORDER_PY if c in BASAL_FUNGI_CLADES],
                fontsize=6, ncol=3, title="Clade", loc="lower right", framealpha=0.8)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/figC_domain_barchart_basalfungi.pdf", dpi=150, bbox_inches="tight")
plt.savefig(f"{OUT_DIR}/figC_domain_barchart_basalfungi.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"[OK] figC_domain_barchart_basalfungi.pdf/.png saved (Basal Fungi only)")

print("[DONE] figC complete.")

# ─── Helper: is this domain helicase-related? ────────────────────────────────
HELICASE_KEYWORDS = [
    "helicase", "dead/deah", "dead box", "deah", "dead_", "sfii",
    "aaa+ atpase", "aaa+", "atp-dependent",
    "snf2", "q_motif", "snf2_n", "swi/snf",  # SNF2 family helicase subdomains
]

def is_helicase_domain(name):
    nl = name.lower()
    # Catch standalone "DEAD" (DEAD-box helicase abbreviation)
    if nl.strip() == "dead":
        return True
    return any(kw in nl for kw in HELICASE_KEYWORDS)

df_pairs_nh = df_pairs[~df_pairs["host_domain"].apply(is_helicase_domain)].copy()
# Also rebuild pair counters filtered
pair_clade_gene_nh = {k: v for k, v in pair_clade_gene_counter.items()
                      if not is_helicase_domain(k[0])}
pair_clade_sp_nh   = {k: v for k, v in pair_clade_sp_counter.items()
                      if not is_helicase_domain(k[0])}

def make_df_top_nh(dft_base, subset_clades=None, top_n=TOP_N):
    """Like make_df_top but using non-helicase counters."""
    if subset_clades is not None:
        tmp_rows = []
        for _, row in dft_base.iterrows():
            key = (row["host_domain"], row["te_superfam"])
            n_g = sum(v for c, v in pair_clade_gene_nh.get(key, {}).items()
                      if c in subset_clades)
            n_s = sum(v for c, v in pair_clade_sp_nh.get(key, {}).items()
                      if c in subset_clades)
            if n_g > 0:
                tmp_rows.append({**row.to_dict(), "n_genes": n_g, "n_species": n_s})
        dft = pd.DataFrame(tmp_rows).sort_values("n_genes", ascending=False).head(top_n).copy()
    else:
        dft = dft_base.head(top_n).copy()
    dft["label"] = dft["host_domain"].str[:38] + "  |  " + dft["te_superfam"]
    for clade in CLADE_ORDER_PY:
        dft[f"clade_{clade}"] = dft.apply(
            lambda r: pair_clade_sp_nh.get((r["host_domain"], r["te_superfam"]), {}).get(clade, 0),
            axis=1)
    return dft

def _save_nh_figure(df_top, subset_clades, title_a, title_b, outbase):
    fig, axes = plt.subplots(1, 2, figsize=(16, 10))
    df_top["color"] = df_top["te_category"].map(CATEGORY_COLORS)
    axes[0].barh(range(len(df_top)), df_top["n_genes"].values[::-1],
                 color=df_top["color"].values[::-1], edgecolor="white", linewidth=0.4)
    axes[0].set_yticks(range(len(df_top)))
    axes[0].set_yticklabels(df_top["label"].values[::-1], fontsize=7)
    axes[0].set_xlabel("Number of co-opted genes", fontsize=9)
    axes[0].set_title(title_a, fontsize=9, loc="left")
    legend_handles = [mpatches.Patch(facecolor=v, label=k) for k, v in CATEGORY_COLORS.items()]
    axes[0].legend(handles=legend_handles, fontsize=7, title="TE class", loc="lower right")
    plot_stacked_clade_panel(axes[1], df_top, "n_species",
        "Number of species", title_b, subset_clades)
    axes[1].legend(
        handles=[mpatches.Patch(facecolor=CLADE_COLORS_PY[c], label=c)
                 for c in CLADE_ORDER_PY if c in subset_clades],
        fontsize=6, ncol=3, title="Clade", loc="lower right", framealpha=0.8)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/{outbase}.pdf", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT_DIR}/{outbase}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[OK] {outbase}.pdf/.png saved")

# All species (no-helicase)
_save_nh_figure(make_df_top_nh(df_pairs_nh, None),
                set(CLADE_ORDER_PY),
                "A   All species — by gene count (no helicase)",
                "B   All species — by clade (no helicase)",
                "figC_domain_barchart_nohelicase")

# Fungi only (no-helicase)
_save_nh_figure(make_df_top_nh(df_pairs_nh, FUNGI_CLADES),
                FUNGI_CLADES,
                "A   Fungi only — by gene count (no helicase)",
                "B   Fungi only — by clade (no helicase)",
                "figC_domain_barchart_fungi_nohelicase")

# Basal Fungi only (no-helicase)
_save_nh_figure(make_df_top_nh(df_pairs_nh, BASAL_FUNGI_CLADES),
                BASAL_FUNGI_CLADES,
                "A   Basal Fungi only — by gene count (no helicase)",
                "B   Basal Fungi only — by clade (no helicase)",
                "figC_domain_barchart_basalfungi_nohelicase")

print("[DONE] figC no-helicase versions complete.")

# ─── Versions excluding Helitron TE entirely ──────────────────────────────────
def is_helitron(te_sf):
    return "helitron" in te_sf.lower()

df_pairs_nte = df_pairs[~df_pairs["te_superfam"].apply(is_helitron)].copy()
pair_clade_gene_nte = {k: v for k, v in pair_clade_gene_counter.items()
                       if not is_helitron(k[1])}
pair_clade_sp_nte   = {k: v for k, v in pair_clade_sp_counter.items()
                       if not is_helitron(k[1])}

def make_df_top_nte(dft_base, subset_clades=None, top_n=TOP_N):
    if subset_clades is not None:
        tmp_rows = []
        for _, row in dft_base.iterrows():
            key = (row["host_domain"], row["te_superfam"])
            n_g = sum(v for c, v in pair_clade_gene_nte.get(key, {}).items()
                      if c in subset_clades)
            n_s = sum(v for c, v in pair_clade_sp_nte.get(key, {}).items()
                      if c in subset_clades)
            if n_g > 0:
                tmp_rows.append({**row.to_dict(), "n_genes": n_g, "n_species": n_s})
        dft = pd.DataFrame(tmp_rows).sort_values("n_genes", ascending=False).head(top_n).copy()
    else:
        dft = dft_base.head(top_n).copy()
    dft["label"] = dft["host_domain"].str[:38] + "  |  " + dft["te_superfam"]
    for clade in CLADE_ORDER_PY:
        dft[f"clade_{clade}"] = dft.apply(
            lambda r: pair_clade_sp_nte.get((r["host_domain"], r["te_superfam"]), {}).get(clade, 0),
            axis=1)
    return dft

def _save_nte_figure(df_top, subset_clades, title_a, title_b, outbase):
    fig, axes = plt.subplots(1, 2, figsize=(16, 10))
    df_top["color"] = df_top["te_category"].map(CATEGORY_COLORS)
    axes[0].barh(range(len(df_top)), df_top["n_genes"].values[::-1],
                 color=df_top["color"].values[::-1], edgecolor="white", linewidth=0.4)
    axes[0].set_yticks(range(len(df_top)))
    axes[0].set_yticklabels(df_top["label"].values[::-1], fontsize=7)
    axes[0].set_xlabel("Number of co-opted genes", fontsize=9)
    axes[0].set_title(title_a, fontsize=9, loc="left")
    legend_handles = [mpatches.Patch(facecolor=v, label=k) for k, v in CATEGORY_COLORS.items()]
    axes[0].legend(handles=legend_handles, fontsize=7, title="TE class", loc="lower right")
    plot_stacked_clade_panel(axes[1], df_top, "n_species",
        "Number of species", title_b, subset_clades)
    axes[1].legend(
        handles=[mpatches.Patch(facecolor=CLADE_COLORS_PY[c], label=c)
                 for c in CLADE_ORDER_PY if c in subset_clades],
        fontsize=6, ncol=3, title="Clade", loc="lower right", framealpha=0.8)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/{outbase}.pdf", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT_DIR}/{outbase}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[OK] {outbase}.pdf/.png saved")

_save_nte_figure(make_df_top_nte(df_pairs_nte, None),
                 set(CLADE_ORDER_PY),
                 "A   All species — by gene count (no Helitron)",
                 "B   All species — by clade (no Helitron)",
                 "figC_domain_barchart_noHelitron")
_save_nte_figure(make_df_top_nte(df_pairs_nte, FUNGI_CLADES),
                 FUNGI_CLADES,
                 "A   Fungi only — by gene count (no Helitron)",
                 "B   Fungi only — by clade (no Helitron)",
                 "figC_domain_barchart_fungi_noHelitron")
_save_nte_figure(make_df_top_nte(df_pairs_nte, BASAL_FUNGI_CLADES),
                 BASAL_FUNGI_CLADES,
                 "A   Basal Fungi only — by gene count (no Helitron)",
                 "B   Basal Fungi only — by clade (no Helitron)",
                 "figC_domain_barchart_basalfungi_noHelitron")
print("[DONE] figC no-Helitron versions complete.")

# =============================================================================
# ─── Alternative Panel B: proportion of clade species with co-option ─────────
# For each (domain_pair, clade): value = n_species_with_cooption / total_in_clade
# Three visual alternatives for figC_domain_barchart_basalfungi_nohelicase
# =============================================================================

import matplotlib.colors as mcolors

# ── Count total analysed species per clade (basal fungi + Dikarya, from fusion files) ──
_PANELB_CLADES = BASAL_FUNGI_CLADES | {"AS", "BA"}
total_clade_sp = collections.Counter()
for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    clade = sp2clade.get(sp, sp[:2])
    if clade in _PANELB_CLADES:
        total_clade_sp[clade] += 1
print(f"[INFO] Total species per clade (basal fungi + Dikarya): {dict(total_clade_sp)}")

# ── Build proportion matrix for basalfungi + Dikarya nohelicase ──
_subset = _PANELB_CLADES
# Row order matches figC_domain_barchart_basalfungi_nohelicase (sorted by basal-fungi gene count)
_df_top = make_df_top_nh(df_pairs_nh, BASAL_FUNGI_CLADES, top_n=TOP_N)
# Custom row order: BA, AS, MU, GL, MO, KX, EN, ZO, BS, OL, CH, NE, BL, MI, CR, AP, NU, HO+AB
DESIRED_CLADE_ORDER = ["BA","AS","MU","GL","MO","KX","EN","ZO","BS","OL",
                       "CH","NE","BL","MI","CR","AP","NU"]
# Include HO only if present and merge with AB if AB exists
_subset_clades_ordered = [c for c in DESIRED_CLADE_ORDER if c in _subset and total_clade_sp.get(c, 0) > 0]
# Add "HO" at the end if present
if "HO" in _subset and total_clade_sp.get("HO", 0) > 0:
    _subset_clades_ordered.append("HO")
_ordered_clades = _subset_clades_ordered
_pairs_labels   = _df_top["label"].values[::-1].tolist()   # bottom→top (same as barh)
_pairs_keys     = list(zip(_df_top["host_domain"].values[::-1], _df_top["te_superfam"].values[::-1]))

def _prop_matrix():
    """Return (n_pairs × n_clades) proportion matrix."""
    mat = np.zeros((len(_pairs_keys), len(_ordered_clades)))
    for pi, key in enumerate(_pairs_keys):
        for ci, clade in enumerate(_ordered_clades):
            n_sp   = pair_clade_sp_nh.get(key, {}).get(clade, 0)
            tot    = total_clade_sp.get(clade, 0)
            mat[pi, ci] = n_sp / tot if tot > 0 else 0.0
    return mat

_mat = _prop_matrix()

# ──────────────────────────────────────────────────────────────────────────────
# Alt B-1: Bubble / dot matrix
# ──────────────────────────────────────────────────────────────────────────────
def _plot_bubble_matrix(ax_or_fig=None):
    fig_bub, ax_bub = plt.subplots(figsize=(14, 12))
    max_radius = 18   # pt² for 100%
    for pi in range(len(_pairs_keys)):
        for ci, clade in enumerate(_ordered_clades):
            prop = _mat[pi, ci]
            if prop <= 0:
                continue
            size = max(prop * max_radius * 100, 4)
            ax_bub.scatter(ci, pi, s=size, color=CLADE_COLORS_PY[clade],
                           edgecolors="grey", linewidths=0.3, alpha=0.85, zorder=3)
            if prop >= 0.15:
                ax_bub.text(ci, pi, f"{int(round(prop*100))}%",
                            ha="center", va="center", fontsize=5.5,
                            color="white" if prop > 0.55 else "black")
    # Reference bubble legend
    for ref_pct, ref_lbl in [(0.25, "25%"), (0.50, "50%"), (1.0, "100%")]:
        ax_bub.scatter([], [], s=ref_pct * max_radius * 100,
                       color="grey", alpha=0.6, label=ref_lbl)
    ax_bub.legend(title="% of clade", fontsize=7, title_fontsize=7,
                  loc="lower right", framealpha=0.8, scatterpoints=1)
    ax_bub.set_xticks(range(len(_ordered_clades)))
    ax_bub.set_xticklabels(_ordered_clades, fontsize=8, rotation=45, ha="right")
    ax_bub.set_yticks(range(len(_pairs_labels)))
    ax_bub.set_yticklabels(_pairs_labels, fontsize=7)
    ax_bub.set_xlim(-0.7, len(_ordered_clades) - 0.3)
    ax_bub.set_ylim(-0.7, len(_pairs_labels) - 0.3)
    ax_bub.grid(True, color="0.90", zorder=0, linewidth=0.5)
    ax_bub.set_title("B (alt-1)   Proportion of clade species with co-option\n"
                     "Basal Fungi + Dikarya — no helicase domains", fontsize=9, loc="left")
    plt.tight_layout()
    out = f"{OUT_DIR}/figC_panelB_bubble_basalfungi_nohelicase"
    fig_bub.savefig(out + ".pdf", dpi=150, bbox_inches="tight")
    fig_bub.savefig(out + ".png", dpi=150, bbox_inches="tight")
    plt.close(fig_bub)
    print(f"[OK] {out}.pdf/.png saved (Alt B-1: bubble matrix, basal fungi + Dikarya)")

_plot_bubble_matrix()

# ──────────────────────────────────────────────────────────────────────────────
# Alt B-2: Heatmap (% of clade) — alternative version with custom colors and row order
# ──────────────────────────────────────────────────────────────────────────────
def _plot_heatmap_alt():
    fig_hm, ax_hm = plt.subplots(figsize=(12, 12))
    # Custom colormap: white → #53B6A4 (teal)
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_teal", ["white", "#53B6A4"])
    vmax = min(100, np.percentile(_mat[_mat > 0] * 100, 95) * 1.1
               if np.any(_mat > 0) else 100)
    norm = plt.Normalize(vmin=0, vmax=vmax)

    n_rows = len(_pairs_keys)
    n_cols = len(_ordered_clades)

    # Draw each cell as an individual Rectangle (editable vector objects in PDF)
    for pi in range(n_rows):
        for ci in range(n_cols):
            v = _mat[pi, ci] * 100
            color = cmap(norm(v)) if v > 0 else "white"
            rect = mpatches.Rectangle((ci - 0.5, pi - 0.5), 1, 1,
                                       facecolor=color, edgecolor="none",
                                       linewidth=0)
            ax_hm.add_patch(rect)

    # Annotate cells
    thresh = (vmax - 0) * 0.6 + 0
    for pi in range(n_rows):
        for ci in range(n_cols):
            v = _mat[pi, ci] * 100
            if v > 0:
                txt_col = "white" if v > thresh else "black"
                ax_hm.text(ci, pi, f"{v:.0f}", ha="center", va="center",
                           fontsize=5.5, color=txt_col)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax_hm, label="% of clade species with co-option",
                 fraction=0.03, pad=0.02)

    ax_hm.set_xlim(-0.5, n_cols - 0.5)
    ax_hm.set_ylim(-0.5, n_rows - 0.5)
    ax_hm.set_aspect("auto")
    ax_hm.set_xticks(range(n_cols))
    ax_hm.set_xticklabels(_ordered_clades, fontsize=8, rotation=45, ha="right")
    ax_hm.set_yticks(range(n_rows))
    ax_hm.set_yticklabels(_pairs_labels, fontsize=7)
    ax_hm.set_title("B (alt-2)   % of clade species with co-option\n"
                    "Basal Fungi + Dikarya — no helicase domains\n"
                    "Row order: BA,AS,MU,GL,MO,KX,EN,ZO,BS,OL,CH,NE,BL,MI,CR,AP,NU,HO",
                    fontsize=9, loc="left")
    plt.tight_layout()
    out = f"{OUT_DIR}/figC_panelB_heatmap_basalfungi_nohelicase_alt"
    fig_hm.savefig(out + ".pdf", dpi=150, bbox_inches="tight")
    fig_hm.savefig(out + ".png", dpi=150, bbox_inches="tight")
    plt.close(fig_hm)
    print(f"[OK] {out}.pdf/.png saved (Alt B-2: heatmap alt — custom colors and row order, editable rectangles)")

_plot_heatmap_alt()

# ──────────────────────────────────────────────────────────────────────────────
# Alt B-3: Stacked bar of proportion (each clade segment = % of clade)
#          Combined with panel A in a 2-panel figure
# ──────────────────────────────────────────────────────────────────────────────
def _plot_prop_stacked():
    fig_ps, axes_ps = plt.subplots(1, 2, figsize=(16, 12))

    # Panel A — same gene count bars as original
    df_rev = _df_top.iloc[::-1].reset_index(drop=True)
    df_rev["color"] = df_rev["te_category"].map(CATEGORY_COLORS)
    axes_ps[0].barh(range(len(df_rev)), df_rev["n_genes"].values,
                    color=df_rev["color"].values, edgecolor="white", linewidth=0.4)
    axes_ps[0].set_yticks(range(len(df_rev)))
    axes_ps[0].set_yticklabels(df_rev["label"].values, fontsize=7)
    axes_ps[0].set_xlabel("Number of co-opted genes", fontsize=9)
    axes_ps[0].set_title("A   Basal Fungi + Dikarya — by gene count (no helicase)", fontsize=9, loc="left")
    axes_ps[0].legend(
        handles=[mpatches.Patch(facecolor=v, label=k) for k, v in CATEGORY_COLORS.items()],
        fontsize=7, title="TE class", loc="lower right")

    # Panel B — proportion stacked bars (already in bottom→top order via _mat)
    ax_b = axes_ps[1]
    lefts = np.zeros(len(_pairs_keys))
    for ci, clade in enumerate(_ordered_clades):
        vals = _mat[:, ci]   # already in bottom→top pair order
        ax_b.barh(range(len(_pairs_keys)), vals * 100, left=lefts,
                  color=CLADE_COLORS_PY[clade], edgecolor="none", height=0.75, label=clade)
        lefts += vals * 100
    ax_b.set_yticks(range(len(_pairs_labels)))
    ax_b.set_yticklabels(_pairs_labels, fontsize=7)
    ax_b.set_xlabel("% of clade species with co-option (summed across clades)", fontsize=9)
    ax_b.set_title("B (alt-3)   Proportion of each clade's species with co-option",
                   fontsize=9, loc="left")
    ax_b.legend(
        handles=[mpatches.Patch(facecolor=CLADE_COLORS_PY[c], label=c)
                 for c in _ordered_clades],
        fontsize=6, ncol=3, title="Clade", loc="lower right", framealpha=0.8)

    plt.tight_layout()
    out = f"{OUT_DIR}/figC_panelB_propstack_basalfungi_nohelicase"
    fig_ps.savefig(out + ".pdf", dpi=150, bbox_inches="tight")
    fig_ps.savefig(out + ".png", dpi=150, bbox_inches="tight")
    plt.close(fig_ps)
    print(f"[OK] {out}.pdf/.png saved (Alt B-3: proportion stacked bars)")

_plot_prop_stacked()

print("[DONE] Alternative panel B figures complete.")
