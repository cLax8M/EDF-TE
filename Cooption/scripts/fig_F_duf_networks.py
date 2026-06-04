#!/usr/bin/env python3
"""
fig_F_duf_networks.py
Co-occurrence network for the top 4 DUF domains most frequently associated with TEs.

For each of the 4 DUF anchors:
  - Show which non-TE / non-DUF domains co-occur in true co-opted genes
  - Node size = n co-occurring genes; edge width = n genes
  - Color edges by clade of origin

Outputs:
  figF_duf_network_{DUF_ID}.pdf/.png
  figF_duf_networks_combined.pdf/.png
"""
import os, sys, glob, collections
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx

# ─── Paths ──────────────────────────────────────────────────────────────────
COOPTION_DIR = "/legserv/Temp/Carlos/EDF-TE/results/Co-option"
TE_MASTER    = f"{COOPTION_DIR}/config/TE_domains_master.tsv"
DUF_FILE     = f"{COOPTION_DIR}/figures/13_duf_te_fusions.tsv"
FUSIONS_GLOB = f"{COOPTION_DIR}/per_species/fusions/*.TE_host_fusions.tsv"
CLASS_GLOB   = f"{COOPTION_DIR}/per_species/classification/*.true_TE_coopted_genes.tsv"
CLADE_FILE   = "/legserv/Temp/Carlos/EDF-TE/results/HTT/clade_assignment.tsv"
OUT_DIR      = f"{COOPTION_DIR}/figures"

TOP_DUF      = 4
TOP_PARTNERS = 15  # top co-occurring domains per DUF anchor

CLADE_COLORS = {
    "HO":"#e07b54","MI":"#b0b0b0","NU":"#d4b8e0","AP":"#5ba3d9",
    "GL":"#2e7d32","MO":"#81c784","MU":"#c8e6c9","ZO":"#ff8a65",
    "EN":"#ffb74d","BS":"#f06292","KX":"#ba68c8","BL":"#64b5f6",
    "CH":"#4dd0e1","CR":"#fff176","OL":"#a5d6a7","AS":"#ff7043",
    "BA":"#8d6e63",
}

# ─── TE master ────────────────────────────────────────────────────────────────
te_master = pd.read_csv(TE_MASTER, sep="\t")
te_all_ids = set(te_master["domain_id"].astype(str))
id2name_te = dict(zip(te_master["domain_id"].astype(str),
                       te_master["domain_name"].astype(str)))

# ─── Clade info ──────────────────────────────────────────────────────────────
clades = pd.read_csv(CLADE_FILE, sep="\t", header=None, names=["species","clade"])
sp2clade = dict(zip(clades["species"], clades["clade"]))

# ─── Find top 4 DUF anchors ───────────────────────────────────────────────────
duf_df = pd.read_csv(DUF_FILE, sep="\t")
# Aggregate by DUF: sum n_species across all TE superfams
duf_agg = duf_df.groupby(["duf_domain_id","duf_domain_name"]).agg(
    n_species_total=("n_species","sum"),
    n_genes_total=("n_genes_copresent","sum")
).reset_index().sort_values("n_species_total", ascending=False)

top_dufs = duf_agg.head(TOP_DUF)
print("[INFO] Top DUF domains:")
print(top_dufs[["duf_domain_id","duf_domain_name","n_species_total","n_genes_total"]].to_string(index=False))
top_duf_ids = set(top_dufs["duf_domain_id"].astype(str))
duf_id2name = dict(zip(top_dufs["duf_domain_id"].astype(str),
                       top_dufs["duf_domain_name"].astype(str)))

# ─── Load true co-opted genes ────────────────────────────────────────────────
true_genes = collections.defaultdict(set)
for f in sorted(glob.glob(CLASS_GLOB)):
    sp = os.path.basename(f).split(".")[0]
    df_c = pd.read_csv(f, sep="\t", usecols=["gene_id"])
    true_genes[sp].update(df_c["gene_id"].astype(str))

# ─── Scan fusions for DUF co-occurrences ─────────────────────────────────────
# duf_partner_counts[duf_id][(partner_domain, clade)] = n_genes
duf_partner_simple = collections.defaultdict(lambda: collections.Counter())  # duf_id → {partner: n_genes}

print("[INFO] Scanning fusion files...")
for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    true_sp = true_genes.get(sp, set())
    clade   = sp2clade.get(sp, sp.split("_")[0])
    if not true_sp:
        continue
    df = pd.read_csv(fpath, sep="\t", dtype=str).fillna("")
    df_true = df[df["gene_id"].isin(true_sp)]

    # Per gene: collect DUF ids and non-DUF/non-TE domains
    for gid, gdf in df_true.groupby("gene_id"):
        duf_in_gene  = set()
        non_duf_names = set()
        for _, row in gdf.iterrows():
            did   = str(row["domain_id"])
            dname = str(row["domain_name"]).strip()
            if not dname or dname in ("nan","-",""):
                dname = did
            if did in top_duf_ids:
                duf_in_gene.add(did)
            elif did not in te_all_ids:
                # non-TE, non-DUF = potential partner — skip TE-derived names and uninformative entries
                _skip = {"", "-", "nan", r"\N", "no description",
                         "reverse transcriptase", "reverse_transcriptase",
                         "transposase", "integrase", "integrase, core domain"}
                if dname not in ("nan","-","",r"\N") and dname.lower() not in _skip:
                    non_duf_names.add(dname[:60])

        for duf_id in duf_in_gene:
            for pname in non_duf_names:
                duf_partner_simple[duf_id][pname] += 1

print(f"[INFO] DUF anchors with partners: {list(duf_partner_simple.keys())}")

COMM_PALETTE_F = ["#66bb6a", "#ab47bc", "#ff8a65", "#ef5350",
                  "#42a5f5", "#ffca28", "#26c6da", "#8d6e63", "#ec407a", "#78909c"]

# ─── Combined force-directed community network ────────────────────────────────
duf_list = top_dufs["duf_domain_id"].astype(str).tolist()
MIN_F = 2  # min n_genes for an edge

G_f = nx.Graph()
for duf_id in duf_list:
    duf_nm = duf_id2name.get(duf_id, duf_id)
    n_sp = top_dufs.loc[top_dufs["duf_domain_id"].astype(str) == duf_id,
                        "n_species_total"].values[0]
    G_f.add_node(duf_nm, ntype="DUF", size=int(n_sp))

for duf_id in duf_list:
    duf_nm = duf_id2name.get(duf_id, duf_id)
    for pname, cnt in duf_partner_simple[duf_id].most_common(TOP_PARTNERS):
        if cnt < MIN_F or not pname or pname in (r"\N", "", "nan", "-"):
            continue
        if pname not in G_f:
            G_f.add_node(pname, ntype="partner", size=cnt)
        else:
            G_f.nodes[pname]["size"] = G_f.nodes[pname].get("size", 0) + cnt
        if G_f.has_edge(duf_nm, pname):
            G_f[duf_nm][pname]["weight"] += cnt
        else:
            G_f.add_edge(duf_nm, pname, weight=cnt)

print(f"[INFO] DUF combined network: {G_f.number_of_nodes()} nodes, "
      f"{G_f.number_of_edges()} edges")

if G_f.number_of_nodes() >= 3:
    communities_f = list(nx.community.greedy_modularity_communities(G_f, weight="weight"))
    node_comm_f = {}
    for ci, comm in enumerate(communities_f):
        for nd in comm:
            node_comm_f[nd] = ci

    pos_f = nx.spring_layout(G_f, k=3.0, iterations=150, seed=42, weight="weight")

    fig_f, ax_f = plt.subplots(figsize=(15, 13))
    ax_f.set_facecolor("#f8f8f8")

    # Community ellipses
    for ci, comm in enumerate(communities_f):
        color = COMM_PALETTE_F[ci % len(COMM_PALETTE_F)]
        pts = np.array([pos_f[n] for n in comm if n in pos_f])
        if len(pts) < 2:
            continue
        cx, cy = pts.mean(axis=0)
        sx = max(pts[:, 0].std(), 0.08) * 3.2
        sy = max(pts[:, 1].std(), 0.08) * 3.2
        el = mpatches.Ellipse((cx, cy), width=sx * 2, height=sy * 2,
                               alpha=0.13, color=color, zorder=0)
        ax_f.add_patch(el)
        ax_f.text(cx, cy + sy + 0.05, f"Cluster {ci + 1}",
                  ha="center", va="bottom", fontsize=9.5, fontweight="bold",
                  color=color, zorder=5)

    # Edges
    max_wf = max((d for _, _, d in G_f.edges(data="weight")), default=1)
    for u, v, w in G_f.edges(data="weight"):
        x0, y0 = pos_f[u]; x1, y1 = pos_f[v]
        lw = 0.4 + 3.5 * (w / max_wf)
        ax_f.plot([x0, x1], [y0, y1], color="#aaaaaa", lw=lw, alpha=0.40, zorder=1)

    # Nodes
    max_sz_f = max(d.get("size", 1) for _, d in G_f.nodes(data=True))
    for nd, (x, y) in pos_f.items():
        ntype = G_f.nodes[nd].get("ntype", "partner")
        ci = node_comm_f.get(nd, 0)
        color = COMM_PALETTE_F[ci % len(COMM_PALETTE_F)]
        sz = G_f.nodes[nd].get("size", 1)
        r = (0.05 + 0.22 * (sz / max_sz_f) ** 0.5) * (1.6 if ntype == "DUF" else 1.0)
        circle = plt.Circle((x, y), r, color=color,
                             ec="black" if ntype == "DUF" else "white",
                             lw=1.8 if ntype == "DUF" else 0.5,
                             alpha=0.93, zorder=3 if ntype == "DUF" else 2)
        ax_f.add_patch(circle)
        label = nd[:28] if ntype == "DUF" else nd[:24]
        ax_f.text(x, y, label, ha="center", va="center",
                  fontsize=7.5 if ntype == "DUF" else 5.5,
                  fontweight="bold" if ntype == "DUF" else "normal",
                  color="white" if ntype == "DUF" else "#1a1a2e", zorder=4)

    ax_f.autoscale()
    ax_f.axis("off")
    ax_f.set_title(
        "DUF domain co-option network\n"
        "Node size \u221d co-opted species \u00b7 edge width \u221d co-occurrences \u00b7 color = community",
        fontsize=12)
    fig_f.tight_layout()
    fig_f.savefig(f"{OUT_DIR}/figF_duf_networks_combined.pdf", dpi=150, bbox_inches="tight")
    fig_f.savefig(f"{OUT_DIR}/figF_duf_networks_combined.png", dpi=150, bbox_inches="tight")
    plt.close(fig_f)
    print(f"[OK] figF_duf_networks_combined.pdf/.png (force-directed community network)")
else:
    print("[WARN] Not enough nodes for DUF combined network")

# ─── Individual star figures per DUF ─────────────────────────────────────────
def draw_duf_network(duf_id, partner_counter, ax, title=""):
    duf_nm = duf_id2name.get(duf_id, duf_id)
    top_partners = partner_counter.most_common(TOP_PARTNERS)
    if not top_partners:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
        return

    G = nx.Graph()
    G.add_node(duf_nm, ntype="DUF", size=1)
    max_w = max(c for _, c in top_partners)
    for pname, cnt in top_partners:
        G.add_node(pname, ntype="partner", size=cnt)
        G.add_edge(duf_nm, pname, weight=cnt)

    # Force-directed layout
    pos = nx.spring_layout(G, k=2.5, iterations=80, seed=0, weight="weight")

    # Edges
    for pname, cnt in top_partners:
        x0, y0 = pos[duf_nm]; x1, y1 = pos[pname]
        ax.plot([x0, x1], [y0, y1], color="#aaaaaa",
                lw=0.3 + 3.5 * cnt / max_w, alpha=0.45, zorder=1)

    # Nodes
    max_sz = max(d.get("size", 1) for _, d in G.nodes(data=True))
    for nd, (x, y) in pos.items():
        ntype = G.nodes[nd].get("ntype", "partner")
        sz = G.nodes[nd].get("size", 1)
        r = (0.06 + 0.20 * (sz / max_sz) ** 0.5) * (1.5 if ntype == "DUF" else 1.0)
        color = "#7b3f00" if ntype == "DUF" else "#2ca02c"
        circle = plt.Circle((x, y), r, color=color,
                             ec="black" if ntype == "DUF" else "white",
                             lw=1.5 if ntype == "DUF" else 0.4,
                             alpha=0.92, zorder=2)
        ax.add_patch(circle)
        label = nd[:25]
        ax.text(x, y, label, ha="center", va="center",
                fontsize=7.5 if ntype == "DUF" else 5.5,
                fontweight="bold" if ntype == "DUF" else "normal",
                color="white", zorder=3)
    ax.autoscale()
    ax.set_title(title or duf_nm[:30], fontsize=8)
    ax.axis("off")


for duf_id in duf_list:
    duf_nm = duf_id2name.get(duf_id, duf_id)
    n_sp   = top_dufs.loc[top_dufs["duf_domain_id"].astype(str) == duf_id,
                          "n_species_total"].values[0]
    fig_i, ax_i = plt.subplots(figsize=(9, 8))
    draw_duf_network(duf_id, duf_partner_simple[duf_id], ax_i,
                     title=f"{duf_nm}\n(n_species = {int(n_sp)}; top {TOP_PARTNERS} partners)")
    safe = duf_id.replace("/", "_").replace(":", "_")
    fig_i.tight_layout()
    fig_i.savefig(f"{OUT_DIR}/figF_duf_network_{safe}.pdf", dpi=150, bbox_inches="tight")
    fig_i.savefig(f"{OUT_DIR}/figF_duf_network_{safe}.png", dpi=150, bbox_inches="tight")
    plt.close(fig_i)
    print(f"  [OK] figF_duf_network_{safe}.pdf/.png")

    # ── Export node/edge TSVs for Cytoscape/Gephi ─────────────────────────────
    partners = duf_partner_simple[duf_id].most_common()
    if partners:
        node_rows = [{"domain": duf_nm, "count": sum(c for _, c in partners), "ntype": "DUF", "community": 0}]
        for pname, cnt in partners[:TOP_PARTNERS]:
            node_rows.append({"domain": pname, "count": cnt, "ntype": "partner", "community": 1})
        edge_rows = [{"domain_u": duf_nm, "domain_v": pname, "weight": cnt}
                     for pname, cnt in partners[:TOP_PARTNERS]]
        pd.DataFrame(node_rows).to_csv(f"{OUT_DIR}/step7_article_network_{safe}.nodes.tsv", sep="\t", index=False)
        pd.DataFrame(edge_rows).to_csv(f"{OUT_DIR}/step7_article_network_{safe}.edges.tsv", sep="\t", index=False)
        print(f"  [OK] TSVs: step7_article_network_{safe}.nodes/edges.tsv")

    # ── Combined DUF network TSVs ──────────────────────────────────────────────
if G_f.number_of_nodes() > 0:
    if G_f.number_of_nodes() >= 3:
        communities_gf = list(nx.community.greedy_modularity_communities(G_f, weight="weight"))
        nc_map = {nd: ci for ci, comm in enumerate(communities_gf) for nd in comm}
    else:
        nc_map = {}
    combined_nodes = [{"domain": nd, "count": int(a.get("size", 0)),
                       "ntype": a.get("ntype", "partner"), "community": nc_map.get(nd, 0)}
                      for nd, a in G_f.nodes(data=True)]
    combined_edges = [{"domain_u": u, "domain_v": v, "weight": int(a.get("weight", 0))}
                      for u, v, a in G_f.edges(data=True)]
    pd.DataFrame(combined_nodes).to_csv(f"{OUT_DIR}/step7_article_network_DUF_combined.nodes.tsv", sep="\t", index=False)
    pd.DataFrame(combined_edges).to_csv(f"{OUT_DIR}/step7_article_network_DUF_combined.edges.tsv", sep="\t", index=False)
    print("[OK] TSVs: step7_article_network_DUF_combined.nodes/edges.tsv")

print("[DONE] Fig F complete.")
