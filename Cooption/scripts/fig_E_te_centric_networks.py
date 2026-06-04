#!/usr/bin/env python3
"""
fig_E_te_centric_networks.py
TE-domain-centric co-occurrence networks.

Inverted from old approach:
  - Take top 3 TE-core domains (by gene count in true co-opted set)
  - For each, build a star/ego network showing co-occurring NON-TE domains
  - Node size = n co-occurrences; edge weight = n genes with this pair
  - Separate figure per TE anchor domain

Output:
  figE_te_network_{TE_DOMAIN}.pdf/.png  (one per anchor)
  figE_te_centric_summary.tsv
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
FUSIONS_GLOB = f"{COOPTION_DIR}/per_species/fusions/*.TE_host_fusions.tsv"
CLASS_GLOB   = f"{COOPTION_DIR}/per_species/classification/*.true_TE_coopted_genes.tsv"
OUT_DIR      = f"{COOPTION_DIR}/figures"

TOP_TE_ANCHORS  = 8    # top TE-core domains to use as anchors
TOP_HOST_SHOWN  = 20   # top co-occurring host domains per anchor

# ─── TE master ───────────────────────────────────────────────────────────────
te_master = pd.read_csv(TE_MASTER, sep="\t")
te_core_ids  = set(te_master.loc[te_master["TE_status"] == "TE_core",  "domain_id"].astype(str))
te_assoc_ids = set(te_master.loc[te_master["TE_status"] == "TE_associated","domain_id"].astype(str))
te_all_ids   = te_core_ids | te_assoc_ids

id2name = dict(zip(te_master["domain_id"].astype(str),
                   te_master["domain_name"].astype(str)))
id2sf   = dict(zip(te_master["domain_id"].astype(str),
                   te_master["TE_superfamily"].astype(str)))

# ─── Load true co-opted gene keys ────────────────────────────────────────────
true_genes = collections.defaultdict(set)
for f in sorted(glob.glob(CLASS_GLOB)):
    sp = os.path.basename(f).split(".")[0]
    df = pd.read_csv(f, sep="\t", usecols=["gene_id"])
    true_genes[sp].update(df["gene_id"].astype(str))
print(f"[INFO] Co-opted genes: {sum(len(v) for v in true_genes.values())}")

# ─── Collect per-gene TE-core and non-TE domains ─────────────────────────────
# te_core_genes[te_domain_id] = {gene_id: set of host domain names}
te_to_host = collections.defaultdict(lambda: collections.Counter())  # te_id → {host_name: n_genes}
te_gene_count = collections.Counter()  # te_id → n_genes carrying it

print("[INFO] Scanning fusion files for TE-centric co-occurrences...")
for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    true_sp = true_genes.get(sp, set())
    if not true_sp:
        continue
    df = pd.read_csv(fpath, sep="\t", dtype=str).fillna("")
    df_true = df[df["gene_id"].isin(true_sp)]

    # Group by gene
    for gid, gdf in df_true.groupby("gene_id"):
        te_core_in_gene  = set()
        host_names_in_gene = set()
        for _, row in gdf.iterrows():
            did   = str(row["domain_id"])
            dname = str(row["domain_name"]).strip()
            if not dname or dname in ("nan", "-", "", r"\N"):
                dname = did
            dclass = str(row.get("domain_class", ""))
            if did in te_core_ids or dclass == "TE_core":
                te_core_in_gene.add(did)
            elif did not in te_assoc_ids and dclass not in ("TE_core","TE_associated"):
                # host domain — skip TE-derived names and uninformative entries
                _skip = {"", "-", "nan", r"\N", "no description",
                         "reverse transcriptase", "reverse_transcriptase",
                         "transposase", "integrase", "integrase, core domain"}
                if dname not in ("-", "nan", "", r"\N") and dname.lower() not in _skip:
                    host_names_in_gene.add(dname[:60])  # truncate long names

        for te_id in te_core_in_gene:
            te_gene_count[te_id] += 1
            for hd in host_names_in_gene:
                te_to_host[te_id][hd] += 1

print(f"[INFO] TE-core domains with gene associations: {len(te_gene_count)}")

# ─── Pick top TE anchors ──────────────────────────────────────────────────────
top_te_ids = [tid for tid, _ in te_gene_count.most_common(TOP_TE_ANCHORS)]
print(f"[INFO] Top {TOP_TE_ANCHORS} TE-core anchor domains:")
for tid in top_te_ids:
    nm = id2name.get(tid, tid)
    sf = id2sf.get(tid, "?")
    print(f"  {tid}: {nm} ({sf}), n_genes={te_gene_count[tid]}")

# ─── Summary TSV ─────────────────────────────────────────────────────────────
rows = []
for te_id, host_cnt in te_to_host.items():
    te_nm = id2name.get(te_id, te_id)
    te_sf = id2sf.get(te_id, "Unknown")
    for hd, n in host_cnt.items():
        rows.append({"te_domain_id": te_id, "te_domain_name": te_nm,
                     "te_superfamily": te_sf, "host_domain": hd, "n_genes": n})
pd.DataFrame(rows).sort_values("n_genes", ascending=False).to_csv(
    f"{OUT_DIR}/figE_te_centric_summary.tsv", sep="\t", index=False)
print(f"[OK] figE_te_centric_summary.tsv")

# ─── Network figure per anchor ────────────────────────────────────────────────
def make_te_network(te_id, host_counter, te_gene_count, top_n=TOP_HOST_SHOWN, out_dir=OUT_DIR):
    te_nm = id2name.get(te_id, te_id)
    te_sf = id2sf.get(te_id, "Unknown")
    n_genes_anchor = te_gene_count[te_id]

    # top N host domains
    top_hosts = host_counter.most_common(top_n)
    if not top_hosts:
        print(f"  [SKIP] {te_id}: no host domains")
        return

    G = nx.Graph()
    # Center: TE node
    G.add_node(te_nm, node_type="TE", size=n_genes_anchor, sf=te_sf)
    for hd, cnt in top_hosts:
        G.add_node(hd, node_type="host", size=cnt)
        G.add_edge(te_nm, hd, weight=cnt)

    # Layout: star with TE center
    pos = {te_nm: np.array([0.0, 0.0])}
    n_hosts = len(top_hosts)
    for i, (hd, _) in enumerate(top_hosts):
        angle = 2 * np.pi * i / n_hosts
        pos[hd] = np.array([1.5 * np.cos(angle), 1.5 * np.sin(angle)])

    fig, ax = plt.subplots(figsize=(11, 9))

    # Edges — width proportional to weight
    max_w = max(cnt for _, cnt in top_hosts)
    for hd, cnt in top_hosts:
        nx.draw_networkx_edges(G, pos, edgelist=[(te_nm, hd)],
                               width=0.5 + 4.5 * cnt / max_w,
                               alpha=0.5, edge_color="#888888", ax=ax)

    # TE center node
    nx.draw_networkx_nodes(G, pos, nodelist=[te_nm],
                           node_size=1200, node_color="#d94801",
                           alpha=0.95, ax=ax)
    # Host nodes
    host_sizes = [max(40, 500 * cnt / max_w) for _, cnt in top_hosts]
    nx.draw_networkx_nodes(G, pos, nodelist=[hd for hd, _ in top_hosts],
                           node_size=host_sizes, node_color="#4393c3",
                           alpha=0.85, ax=ax)

    # Labels
    nx.draw_networkx_labels(G, pos, labels={te_nm: te_nm[:30]},
                            font_size=9, font_color="white",
                            font_weight="bold", ax=ax)
    host_labels = {hd: hd[:28] for hd, _ in top_hosts}
    nx.draw_networkx_labels(G, pos, labels=host_labels,
                            font_size=6.5, font_color="#1a1a2e", ax=ax)

    ax.set_title(f"TE domain: {te_nm}\n"
                 f"Superfamily: {te_sf}  |  genes with this TE domain: {n_genes_anchor:,}\n"
                 f"Top {top_n} co-occurring host domains (edge width = n co-opted genes)",
                 fontsize=10)
    ax.axis("off")

    te_patch  = mpatches.Patch(color="#d94801", label=f"TE domain: {te_nm[:20]}")
    h_patch   = mpatches.Patch(color="#4393c3", label="Host domain (size = n genes)")
    ax.legend(handles=[te_patch, h_patch], fontsize=9, loc="lower left")

    safe_name = te_id.replace("/","_").replace(":","_")
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_dir}/figE_te_network_{safe_name}.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] figE_te_network_{safe_name}.pdf/.png")


for te_id in top_te_ids:
    make_te_network(te_id, te_to_host[te_id], te_gene_count)

# ─── Combined force-directed community network ───────────────────────────────
MIN_COOCCUR = 3   # min n_genes to include host node in combined graph
TOP_PER_TE  = 20  # host domains per TE anchor

COMM_PALETTE = ["#66bb6a", "#ab47bc", "#ff8a65", "#ef5350",
                "#42a5f5", "#ffca28", "#26c6da", "#8d6e63", "#ec407a", "#78909c"]

G_all = nx.Graph()
for te_id in top_te_ids:
    te_nm = id2name.get(te_id, te_id)
    te_sf = id2sf.get(te_id, "Unknown")
    G_all.add_node(te_nm, node_type="TE", size=te_gene_count[te_id], sf=te_sf)

for te_id in top_te_ids:
    te_nm = id2name.get(te_id, te_id)
    for hd, cnt in te_to_host[te_id].most_common(TOP_PER_TE):
        if cnt < MIN_COOCCUR or not hd or hd in (r"\N", "", "nan", "-"):
            continue
        if hd not in G_all:
            G_all.add_node(hd, node_type="host", size=cnt)
        else:
            G_all.nodes[hd]["size"] = G_all.nodes[hd].get("size", 0) + cnt
        if G_all.has_edge(te_nm, hd):
            G_all[te_nm][hd]["weight"] += cnt
        else:
            G_all.add_edge(te_nm, hd, weight=cnt)

print(f"[INFO] Combined network: {G_all.number_of_nodes()} nodes, "
      f"{G_all.number_of_edges()} edges")

if G_all.number_of_nodes() >= 3:
    # Community detection
    communities = list(nx.community.greedy_modularity_communities(G_all, weight="weight"))
    node_comm = {}
    for ci, comm in enumerate(communities):
        for nd in comm:
            node_comm[nd] = ci

    # Force-directed layout
    pos = nx.spring_layout(G_all, k=3.0, iterations=150, seed=42, weight="weight")

    fig_c, ax_c = plt.subplots(figsize=(15, 13))
    ax_c.set_facecolor("#f8f8f8")

    # Community background ellipses
    for ci, comm in enumerate(communities):
        color = COMM_PALETTE[ci % len(COMM_PALETTE)]
        pts = np.array([pos[n] for n in comm if n in pos])
        if len(pts) < 2:
            continue
        cx, cy = pts.mean(axis=0)
        sx = max(pts[:, 0].std(), 0.08) * 3.2
        sy = max(pts[:, 1].std(), 0.08) * 3.2
        el = matplotlib.patches.Ellipse((cx, cy), width=sx * 2, height=sy * 2,
                                         alpha=0.13, color=color, zorder=0)
        ax_c.add_patch(el)
        ax_c.text(cx, cy + sy + 0.05, f"Cluster {ci + 1}",
                  ha="center", va="bottom", fontsize=9.5, fontweight="bold",
                  color=color, zorder=5)

    # Edges
    max_w = max((d for _, _, d in G_all.edges(data="weight")), default=1)
    for u, v, w in G_all.edges(data="weight"):
        x0, y0 = pos[u]; x1, y1 = pos[v]
        lw = 0.4 + 3.5 * (w / max_w)
        ax_c.plot([x0, x1], [y0, y1], color="#aaaaaa", lw=lw, alpha=0.40, zorder=1)

    # Nodes
    max_sz = max(d.get("size", 1) for _, d in G_all.nodes(data=True))
    for nd, (x, y) in pos.items():
        ntype = G_all.nodes[nd].get("node_type", "host")
        ci = node_comm.get(nd, 0)
        color = COMM_PALETTE[ci % len(COMM_PALETTE)]
        sz = G_all.nodes[nd].get("size", 1)
        r = (0.05 + 0.22 * (sz / max_sz) ** 0.5) * (1.6 if ntype == "TE" else 1.0)
        circle = plt.Circle((x, y), r, color=color,
                             ec="black" if ntype == "TE" else "white",
                             lw=1.8 if ntype == "TE" else 0.5,
                             alpha=0.93, zorder=3 if ntype == "TE" else 2)
        ax_c.add_patch(circle)
        label = nd[:28] if ntype == "TE" else nd[:24]
        ax_c.text(x, y, label, ha="center", va="center",
                  fontsize=7.5 if ntype == "TE" else 5.5,
                  fontweight="bold" if ntype == "TE" else "normal",
                  color="white" if ntype == "TE" else "#1a1a2e", zorder=4)

    ax_c.autoscale()
    ax_c.axis("off")
    ax_c.set_title(
        "TE-domain co-option network\n"
        "Node size \u221d co-opted genes \u00b7 edge width \u221d co-occurrences \u00b7 color = community",
        fontsize=12)
    fig_c.tight_layout()
    fig_c.savefig(f"{OUT_DIR}/figE_te_networks_combined.pdf", dpi=150, bbox_inches="tight")
    fig_c.savefig(f"{OUT_DIR}/figE_te_networks_combined.png", dpi=150, bbox_inches="tight")
    plt.close(fig_c)
    print(f"[OK] figE_te_networks_combined.pdf/.png (force-directed community network)")
else:
    print("[WARN] Not enough nodes for combined network")
print("[DONE] Fig E complete.")
