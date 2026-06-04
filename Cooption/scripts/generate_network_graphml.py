#!/usr/bin/env python3
"""
generate_network_graphml.py
Convert step7 node/edge TSVs to GraphML files importable by Cytoscape or Gephi.

For each anchor domain (Chromo, PHD) produces:
  step7_article_network_<ANCHOR>.graphml

Node attributes: id (domain name), count (co-occurrence count), community
Edge attributes: weight (co-occurrence weight)

Usage:
  python3 generate_network_graphml.py [--figures-dir PATH] [--out-dir PATH]
"""

from __future__ import annotations

import argparse
import os
import glob
import xml.etree.ElementTree as ET
from xml.dom import minidom

import pandas as pd

FIGURES_DIR = "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"


def build_graphml(nodes_df: pd.DataFrame, edges_df: pd.DataFrame) -> str:
    """Build a GraphML XML string from node and edge DataFrames."""
    # Root
    root = ET.Element("graphml", {
        "xmlns": "http://graphml.graphdrawing.org/graphml",
        "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
        "xsi:schemaLocation": (
            "http://graphml.graphdrawing.org/graphml "
            "http://graphml.graphdrawing.org/graphml/graphml.xsd"
        ),
    })

    # Key declarations
    ET.SubElement(root, "key", {
        "id": "d_count", "for": "node",
        "attr.name": "count", "attr.type": "integer"
    })
    ET.SubElement(root, "key", {
        "id": "d_community", "for": "node",
        "attr.name": "community", "attr.type": "integer"
    })
    ET.SubElement(root, "key", {
        "id": "d_weight", "for": "edge",
        "attr.name": "weight", "attr.type": "double"
    })

    graph = ET.SubElement(root, "graph", {
        "id": "G", "edgedefault": "undirected"
    })

    # Nodes
    node_ids = {}
    for _, row in nodes_df.iterrows():
        nid = str(row["domain"])
        node_ids[nid] = nid
        node_el = ET.SubElement(graph, "node", {"id": nid})
        d = ET.SubElement(node_el, "data", {"key": "d_count"})
        d.text = str(int(row.get("count", 0)))
        d2 = ET.SubElement(node_el, "data", {"key": "d_community"})
        d2.text = str(int(row.get("community", 0)))

    # Edges
    for i, row in edges_df.iterrows():
        u = str(row["domain_u"])
        v = str(row["domain_v"])
        # Add nodes if not already present (can happen in edge-only entries)
        for nid in (u, v):
            if nid not in node_ids:
                ET.SubElement(graph, "node", {"id": nid})
                node_ids[nid] = nid
        edge_el = ET.SubElement(graph, "edge", {
            "id": f"e{i}", "source": u, "target": v
        })
        dw = ET.SubElement(edge_el, "data", {"key": "d_weight"})
        dw.text = str(float(row.get("weight", 1)))

    # Pretty-print
    raw = ET.tostring(root, encoding="unicode")
    dom = minidom.parseString(raw)
    return dom.toprettyxml(indent="  ", encoding=None)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures-dir", default=FIGURES_DIR)
    ap.add_argument("--out-dir", default=None,
                    help="Output directory (default: same as figures-dir)")
    args = ap.parse_args()

    out_dir = args.out_dir or args.figures_dir
    os.makedirs(out_dir, exist_ok=True)

    # Find all anchor names from nodes files
    pattern = os.path.join(args.figures_dir, "step7_article_network_*.nodes.tsv")
    node_files = sorted(glob.glob(pattern))

    if not node_files:
        print(f"[WARN] No node TSV files found in {args.figures_dir}")
        return

    for nf in node_files:
        anchor = os.path.basename(nf).replace("step7_article_network_", "").replace(".nodes.tsv", "")
        ef = nf.replace(".nodes.tsv", ".edges.tsv")
        if not os.path.exists(ef):
            print(f"[SKIP] No edges file for {anchor}")
            continue

        nodes_df = pd.read_csv(nf, sep="\t") if os.path.getsize(nf) > 1 else pd.DataFrame(columns=["domain","count","community"])
        edges_df = pd.read_csv(ef, sep="\t") if os.path.getsize(ef) > 1 else pd.DataFrame(columns=["domain_u","domain_v","weight"])
        if nodes_df.empty:
            print(f"[SKIP] Empty nodes for {anchor}")
            continue

        print(f"[{anchor}] {len(nodes_df)} nodes, {len(edges_df)} edges")

        graphml_str = build_graphml(nodes_df, edges_df)

        out_path = os.path.join(out_dir, f"step7_article_network_{anchor}.graphml")
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(graphml_str)
        print(f"[OK] Saved {out_path}")


if __name__ == "__main__":
    main()
