#!/usr/bin/env python3


# builds TF–TE–gene interaction networks integrating multiple analysis layers:

# TF binding (FIMO), chimera formation, and TE domain co-option

# imports
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# clade colors (kept consistent with upstream analyses)
CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

# node type colors for network preview
NODE_COLORS = {
    "TE":        "#e41a1c",
    "TF":        "#4daf4a",
    "gene":      "#377eb8",
    "cooption":  "#984ea3",
}

# parse arguments
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tf_te_tsv",     required=True)
    p.add_argument("--chimte_tsv",    required=True)
    p.add_argument("--cooption_tsv",  required=True)
    p.add_argument("--clade_map",     required=True)
    p.add_argument("--focus_species", default="",
                   help="Comma-separated species to include (default: all)")
    p.add_argument("--out_dir",       default="network_out")
    return p.parse_args()

# build interaction network from multiple evidence layers
def build_network(tf_te: pd.DataFrame, chimte: pd.DataFrame,
                   cooption: pd.DataFrame,
                   focus_species: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = {}  # name → {type, label, ...}
    edges = []  # list of {source, target, edge_type, species}

    def add_node(name, node_type):
        if name not in nodes:
            nodes[name] = {"id": name, "type": node_type}

    # TF → TE edges
    if not tf_te.empty and "te_class" in tf_te.columns:
        for _, row in tf_te.iterrows():
            if focus_species and row.get("species") not in focus_species:
                continue
            if pd.isna(row.get("te_class")):
                continue
            motif = str(row.get("motif_id", row.get("motif_alt_id", "unknown_TF")))
            te    = str(row["te_class"])
            add_node(motif, "TF")
            add_node(te, "TE")
            edges.append({"source": motif, "target": te,
                          "edge_type": "TF_binds_TE",
                          "species": row.get("species", "")})

    # TE → gene edges via chimeric transcripts
    if not chimte.empty and "te_class" in chimte.columns:
        gene_col = next((c for c in ["gene_id", "gene", "transcript"] if c in chimte.columns), None)
        for _, row in chimte.iterrows():
            if focus_species and row.get("species") not in focus_species:
                continue
            te   = str(row.get("te_class", "Unclassified"))
            gene = str(row[gene_col]) if gene_col else "unknown_gene"
            add_node(te, "TE")
            add_node(gene, "gene")
            edges.append({"source": te, "target": gene,
                          "edge_type": "chimeric_fusion",
                          "species": row.get("species", "")})

    # OG → TE edges (co-option)
    if not cooption.empty and "te_domains" in cooption.columns:
        for _, row in cooption[cooption.get("blast_confirmed", pd.Series(True, index=cooption.index))].iterrows():
            og = row.get("orthogroup", "unknown_OG")
            for dom in str(row["te_domains"]).split(","):
                dom = dom.strip()
                if dom:
                    add_node(dom, "cooption")
                    add_node(og, "cooption")
                    edges.append({"source": og, "target": dom,
                                  "edge_type": "TE_domain_cooption",
                                  "species": ""})

    nodes_df = pd.DataFrame(list(nodes.values()))
    edges_df = pd.DataFrame(edges)
    return nodes_df, edges_df

# network preview plot 
def plot_network_preview(nodes: pd.DataFrame, edges: pd.DataFrame,
                          out_dir: Path) -> None:
    """Simple spring-layout preview (not Cytoscape quality)."""
    if nodes.empty or edges.empty:
        return

    try:
        import networkx as nx
    except ImportError:
        print("networkx not available — skipping network preview plot.")
        return

    G = nx.DiGraph()
    for _, n in nodes.iterrows():
        G.add_node(n["id"], node_type=n["type"])
    for _, e in edges.iterrows():
        if G.has_node(e["source"]) and G.has_node(e["target"]):
            G.add_edge(e["source"], e["target"], edge_type=e["edge_type"])

    # Limit to largest connected component
    G_un = G.to_undirected()
    largest_cc = max(nx.connected_components(G_un), key=len)
    G = G.subgraph(largest_cc).copy()

    pos = nx.spring_layout(G, k=0.5, seed=42)
    node_colors = [NODE_COLORS.get(G.nodes[n].get("node_type", "gene"), "#aaaaaa")
                   for n in G.nodes]

    fig, ax = plt.subplots(figsize=(12, 12))
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=60, alpha=0.8)
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#cccccc",
                           arrows=True, arrowsize=10, width=0.5)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=4)

    patches = [mpatches.Patch(color=c, label=k) for k, c in NODE_COLORS.items()]
    ax.legend(handles=patches, fontsize=8, frameon=False)
    ax.axis("off")
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "te_interaction_network.pdf")
    fig.savefig(out_dir / "te_interaction_network.png", dpi=200)
    plt.close(fig)
    print("Saved: te_interaction_network.pdf")



# main


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    focus_species = [s.strip() for s in args.focus_species.split(",") if s.strip()]

    def safe_load(path):
        p = Path(path)
        return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()

    tf_te    = safe_load(args.tf_te_tsv)
    chimte   = safe_load(args.chimte_tsv)
    cooption = safe_load(args.cooption_tsv)

    print("Building interaction network ...")
    nodes, edges = build_network(tf_te, chimte, cooption, focus_species)
    print(f"  {len(nodes)} nodes, {len(edges)} edges")

    nodes.to_csv(out_dir / "network_nodes.tsv", sep="\t", index=False)
    edges.to_csv(out_dir / "network_edges.tsv", sep="\t", index=False)

    plot_network_preview(nodes, edges, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
