#!/usr/bin/env python3
"""
fig_EF_cytoscape_inputs.py
Generate Cytoscape-ready edge and node tables for figE and figF.

figE — TE-centric networks (3 top TE-core domains):
  For each TE domain: central TE node, host-domain nodes, species nodes.
  Each species node connects to the TE domain it co-opts with.
  Species are sized by n_genes and coloured by clade.

figF — DUF-centric networks (top 4 DUF domains):
  Each DUF domain connects to all co-occurring partner domains (incl. TE).
  Each (species × partner_domain) is a unique node with a shared domain_group
  attribute so Cytoscape can group/aggregate identical domain names.

Output files in {OUT_DIR}/cytoscape/:
  figE_{TE_ID}_edges.tsv   figE_{TE_ID}_nodes.tsv
  figF_{DUF_ID}_edges.tsv  figF_{DUF_ID}_nodes.tsv
"""

import os, sys, glob, collections
import pandas as pd

# ─── Paths ───────────────────────────────────────────────────────────────────
COOPTION_DIR  = "/legserv/Temp/Carlos/EDF-TE/results/Co-option"
TE_MASTER     = f"{COOPTION_DIR}/config/TE_domains_master.tsv"
FUSIONS_GLOB  = f"{COOPTION_DIR}/per_species/fusions/*.TE_host_fusions.tsv"
CLASS_GLOB    = f"{COOPTION_DIR}/per_species/classification/*.true_TE_coopted_genes.tsv"
DUF_FILE      = f"{COOPTION_DIR}/figures/13_duf_te_fusions.tsv"
CLADE_FILE    = "/legserv/Temp/Carlos/EDF-TE/results/HTT/clade_assignment.tsv"
OUT_DIR       = f"{COOPTION_DIR}/cytoscape"
os.makedirs(OUT_DIR, exist_ok=True)

TOP_TE_ANCHORS = 3
TOP_DUF        = 4
MIN_GENES_EDGE = 2      # minimum co-occurring genes to include edge (figE/figF)
MAX_SP_PER_TE  = 60     # cap species per TE anchor for readability

# ─── Canonical clade palette ─────────────────────────────────────────────────
CLADE_COLORS = {
    "AS":"#e8887a","BA":"#4472c4","MU":"#e07840","MO":"#f4a8b0","GL":"#8b5e48",
    "KX":"#f0a820","ZO":"#c8dc70","EN":"#e01010","BS":"#d0b4e0","OL":"#3668b0",
    "CH":"#e0d040","NE":"#9060b0","BL":"#7030a0","MI":"#3aab98","CR":"#a8c8e8",
    "AP":"#5c9e3a","NU":"#1a2850","HO":"#b4d0e8",
}
CLADE_NAMES = {
    "AS":"Ascomycota","BA":"Basidiomycota","MU":"Mucoromycotina",
    "MO":"Mortierellomycotina","GL":"Glomeromycotina","KX":"Kickxellomycotina",
    "ZO":"Zoopagomycotina","EN":"Entomophthoromycotina","BS":"Basidiobolomycota",
    "OL":"Olpidiomycota","CH":"Chytridiomycetes","NE":"Neocallimastigomycetes",
    "BL":"Blastocladiomycota","MI":"Microsporidia","CR":"Cryptomycota",
    "AP":"Aphelidiomycota","NU":"Nucleariida","HO":"Holozoa+Amoebozoa",
}
NE_SPECIES = {"CH_Anasp1","CH_Neocon1","CH_Neolan1","CH_Orpsp1",
              "CH_Pecora1","CH_Caecom1","CH_Pirfi3"}

def get_clade(sp):
    return "NE" if sp in NE_SPECIES else sp[:2]

# ─── Load TE master ──────────────────────────────────────────────────────────
te_master    = pd.read_csv(TE_MASTER, sep="\t")
te_core_ids  = set(te_master.loc[te_master["TE_status"]=="TE_core","domain_id"].astype(str))
te_assoc_ids = set(te_master.loc[te_master["TE_status"]=="TE_associated","domain_id"].astype(str))
te_all_ids   = te_core_ids | te_assoc_ids
id2name      = dict(zip(te_master["domain_id"].astype(str), te_master["domain_name"].astype(str)))
id2sf        = dict(zip(te_master["domain_id"].astype(str), te_master["TE_superfamily"].astype(str)))

# ─── Load true co-opted genes ────────────────────────────────────────────────
true_genes = collections.defaultdict(set)
for f in sorted(glob.glob(CLASS_GLOB)):
    sp = os.path.basename(f).split(".")[0]
    df = pd.read_csv(f, sep="\t", usecols=["gene_id"])
    true_genes[sp].update(df["gene_id"].astype(str))

# ─── Scan fusions: TE-core × host-domain × species counts ────────────────────
# te_host_sp[te_id][(host_domain, species)] = n_genes
te_host_sp    = collections.defaultdict(lambda: collections.Counter())
te_gene_count = collections.Counter()

print("[INFO] Scanning fusion files ...")
for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    true_sp = true_genes.get(sp, set())
    if not true_sp:
        continue
    df = pd.read_csv(fpath, sep="\t", dtype=str).fillna("")
    df_true = df[df["gene_id"].isin(true_sp)]

    for gid, gdf in df_true.groupby("gene_id"):
        te_in_gene   = set()
        host_in_gene = set()
        for _, row in gdf.iterrows():
            did   = str(row["domain_id"])
            dname = str(row["domain_name"]).strip()
            dclass = str(row.get("domain_class",""))
            if not dname or dname in ("nan","-","",r"\N"):
                dname = did
            if did in te_core_ids or dclass == "TE_core":
                te_in_gene.add(did)
            elif did not in te_assoc_ids and dclass not in ("TE_core","TE_associated"):
                if dname not in ("nan","-","",r"\N"):
                    host_in_gene.add(dname[:60])
        for te_id in te_in_gene:
            te_gene_count[te_id] += 1
            for hd in host_in_gene:
                te_host_sp[te_id][(hd, sp)] += 1

print(f"[INFO] {len(te_gene_count)} TE-core domains with gene associations")

# Top TE anchors
top_te_ids = [tid for tid, _ in te_gene_count.most_common(TOP_TE_ANCHORS)]
print(f"[INFO] Top {TOP_TE_ANCHORS} TE anchors: {top_te_ids}")

# ────────────────────────────────────────────────────────────────────────────
# figE Cytoscape files (one per TE anchor)
# ────────────────────────────────────────────────────────────────────────────
for te_id in top_te_ids:
    te_nm   = id2name.get(te_id, te_id)
    te_sf   = id2sf.get(te_id, "Unknown")
    counter = te_host_sp[te_id]  # (host_domain, species) -> n_genes

    # Aggregate: (host_domain, species) -> n_genes
    # Also compute totals per host_domain and per species
    host_genes = collections.Counter()   # host_domain -> total n_genes
    sp_genes   = collections.Counter()   # species -> total n_genes

    for (hd, sp), n in counter.items():
        if n < MIN_GENES_EDGE:
            continue
        host_genes[hd] += n
        sp_genes[sp] += n

    # Limit to top species for readability
    top_species = [sp for sp, _ in sp_genes.most_common(MAX_SP_PER_TE)]

    edges = []
    nodes = {}

    # Central TE domain node
    te_node_id = f"TE::{te_id}::{te_nm}"
    nodes[te_node_id] = {
        "id": te_node_id, "label": te_nm[:40], "type": "TE_domain",
        "domain_id": te_id, "superfamily": te_sf,
        "n_genes": te_gene_count[te_id],
        "clade": "", "clade_name": "", "clade_color": "#d94801",
        "species": "", "domain_group": te_nm[:40],
    }

    for (hd, sp), n in counter.items():
        if n < MIN_GENES_EDGE or sp not in top_species:
            continue
        clade = get_clade(sp)

        # Host-domain node (shared across species that have same domain)
        hd_node_id = f"HOST::{hd}"
        if hd_node_id not in nodes:
            nodes[hd_node_id] = {
                "id": hd_node_id, "label": hd[:40], "type": "host_domain",
                "domain_id": "", "superfamily": "",
                "n_genes": host_genes[hd],
                "clade": "", "clade_name": "", "clade_color": "#4393c3",
                "species": "", "domain_group": hd[:40],
            }

        # Species node (unique per species)
        sp_node_id = f"SP::{sp}"
        if sp_node_id not in nodes:
            nodes[sp_node_id] = {
                "id": sp_node_id, "label": sp, "type": "species",
                "domain_id": "", "superfamily": "",
                "n_genes": sp_genes[sp],
                "clade": clade, "clade_name": CLADE_NAMES.get(clade,""),
                "clade_color": CLADE_COLORS.get(clade,"#aaaaaa"),
                "species": sp, "domain_group": "",
            }
        else:
            # Update n_genes (cumulative across host domains)
            nodes[sp_node_id]["n_genes"] = max(nodes[sp_node_id]["n_genes"], sp_genes[sp])

        # Edges: TE → host_domain and host_domain → species
        edges.append({"source": te_node_id, "target": hd_node_id,
                      "interaction": "TE_cooccurs_with_host", "weight": host_genes[hd],
                      "clade": ""})
        edges.append({"source": hd_node_id, "target": sp_node_id,
                      "interaction": "domain_in_species", "weight": n,
                      "clade": clade})

    safe = te_id.replace("/","_").replace(":","_").replace(" ","_")
    edge_df = pd.DataFrame(edges).drop_duplicates()
    node_df = pd.DataFrame(nodes.values())
    edge_df.to_csv(f"{OUT_DIR}/figE_{safe}_edges.tsv", sep="\t", index=False)
    node_df.to_csv(f"{OUT_DIR}/figE_{safe}_nodes.tsv", sep="\t", index=False)
    print(f"[OK] figE_{safe}: {len(node_df)} nodes, {len(edge_df)} edges")

# ────────────────────────────────────────────────────────────────────────────
# figF — DUF networks: DUF → all partner domains (incl. TE) × species
# ────────────────────────────────────────────────────────────────────────────
duf_df = pd.read_csv(DUF_FILE, sep="\t")
duf_agg = duf_df.groupby(["duf_domain_id","duf_domain_name"]).agg(
    n_sp_total=("n_species","sum"),
    n_genes_total=("n_genes_copresent","sum")
).reset_index().sort_values("n_sp_total", ascending=False)
top_dufs = duf_agg.head(TOP_DUF)
top_duf_ids = set(top_dufs["duf_domain_id"].astype(str))
duf_id2name = dict(zip(top_dufs["duf_domain_id"].astype(str),
                       top_dufs["duf_domain_name"].astype(str)))

print(f"\n[INFO] Top {TOP_DUF} DUF anchors: {list(top_duf_ids)}")

# Scan fusions for DUF co-occurrences (include ALL domains, even TE)
# duf_partner_sp[duf_id][(partner_name, partner_type, species)] = n_genes
duf_partner_sp = collections.defaultdict(lambda: collections.Counter())
duf_gene_count = collections.Counter()

print("[INFO] Scanning fusion files for DUF ...")
for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    true_sp = true_genes.get(sp, set())
    if not true_sp:
        continue
    df = pd.read_csv(fpath, sep="\t", dtype=str).fillna("")
    df_true = df[df["gene_id"].isin(true_sp)]

    for gid, gdf in df_true.groupby("gene_id"):
        duf_in_gene   = {}   # duf_id -> duf_name
        other_in_gene = {}   # domain_name -> type (TE_domain|host_domain)
        for _, row in gdf.iterrows():
            did   = str(row["domain_id"])
            dname = str(row["domain_name"]).strip()
            dclass = str(row.get("domain_class",""))
            if not dname or dname in ("nan","-","",r"\N"):
                dname = did
            if did in top_duf_ids:
                duf_in_gene[did] = duf_id2name.get(did, did)
            elif did in te_core_ids or dclass == "TE_core":
                te_nm_local = id2name.get(did, did)
                other_in_gene[te_nm_local[:60]] = "TE_domain"
            else:
                if dname not in ("nan","-","",r"\N"):
                    other_in_gene[dname[:60]] = "host_domain"

        for duf_id in duf_in_gene:
            duf_gene_count[duf_id] += 1
            for pname, ptype in other_in_gene.items():
                duf_partner_sp[duf_id][(pname, ptype, sp)] += 1

# Build Cytoscape files for each DUF
for duf_id in duf_agg.head(TOP_DUF)["duf_domain_id"].astype(str).tolist():
    duf_nm = duf_id2name.get(duf_id, duf_id)
    counter = duf_partner_sp[duf_id]

    # Aggregate partner totals
    partner_total = collections.Counter()
    sp_total      = collections.Counter()
    for (pname, ptype, sp), n in counter.items():
        if n < MIN_GENES_EDGE:
            continue
        partner_total[pname] += n
        sp_total[sp] += n

    edges = []
    nodes = {}

    # Central DUF node
    duf_node_id = f"DUF::{duf_id}::{duf_nm}"
    nodes[duf_node_id] = {
        "id": duf_node_id, "label": duf_nm[:40], "type": "DUF_domain",
        "domain_id": duf_id, "domain_group": duf_nm[:40],
        "n_genes": duf_gene_count.get(duf_id, 0),
        "clade": "", "clade_name": "", "clade_color": "#7b3f00",
        "species": "",
    }

    for (pname, ptype, sp), n in counter.items():
        if n < MIN_GENES_EDGE:
            continue
        clade = get_clade(sp)

        # Partner × species node — unique per (domain_name, species)
        p_sp_node_id = f"{ptype}::{sp}::{pname}"
        if p_sp_node_id not in nodes:
            nodes[p_sp_node_id] = {
                "id": p_sp_node_id,
                "label": f"{sp} | {pname[:30]}",
                "type": ptype,
                "domain_id": "",
                "domain_group": pname[:60],  # same domain_group for same domain name
                "n_genes": n,
                "clade": clade,
                "clade_name": CLADE_NAMES.get(clade,""),
                "clade_color": CLADE_COLORS.get(clade,"#aaaaaa"),
                "species": sp,
            }
        else:
            nodes[p_sp_node_id]["n_genes"] += n

        edges.append({
            "source": duf_node_id,
            "target": p_sp_node_id,
            "interaction": f"DUF_cooccurs_{ptype}",
            "weight": n,
            "clade": clade,
            "partner_type": ptype,
            "domain_group": pname[:60],
        })

    safe = duf_id.replace("/","_").replace(":","_").replace(" ","_")
    edge_df = pd.DataFrame(edges).drop_duplicates()
    node_df = pd.DataFrame(nodes.values())
    edge_df.to_csv(f"{OUT_DIR}/figF_{safe}_edges.tsv", sep="\t", index=False)
    node_df.to_csv(f"{OUT_DIR}/figF_{safe}_nodes.tsv", sep="\t", index=False)
    print(f"[OK] figF_{safe}: {len(node_df)} nodes, {len(edge_df)} edges")

print("\n[DONE] Cytoscape input files written to:", OUT_DIR)
print("Import in Cytoscape via: File > Import > Network from File (edges.tsv)")
print("Then: File > Import > Table from File (nodes.tsv) to add node attributes")
print("\nFor figF: use 'domain_group' column as shared attribute to group identical domains.")
print("Style: Map clade_color -> Fill Color, n_genes -> Size, type -> Shape")
