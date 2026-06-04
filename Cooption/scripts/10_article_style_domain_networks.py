#!/usr/bin/env python3

"""Step 7 article-style domain co-occurrence networks (Fig.5c-like).

Builds accessory-domain networks centered on anchor core domains.
- Node: accessory domain co-occurring with anchor
- Node size: number of genes where accessory co-occurs with anchor
- Edge: co-occurrence between accessory domains in anchor-positive genes
- Communities: asynchronous label propagation (networkx)

Outputs:
- Combined PDF/PNG panel for all anchors
- Per-anchor node and edge TSVs with community labels
"""

from __future__ import annotations

import argparse
import glob
import itertools
import os
import re
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Set, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Polygon

try:
    from scipy.spatial import ConvexHull
    HAS_HULL = True
except Exception:
    HAS_HULL = False

BAD = {"", "-", "\\N", "unknown", "nan", "none", "no description", "n/a", "na"}


DOMAIN_CANON_RULES: List[Tuple[re.Pattern[str], str]] = [
    (re.compile(r"^(?:\\\\?n|n/?a|na)$", re.IGNORECASE), ""),
    (re.compile(r"integrase|\brve\b", re.IGNORECASE), "rve"),
    (re.compile(r"reverse transcriptase|\brvt[_ -]?1\b|\brt_pol\b", re.IGNORECASE), "RVT_1"),
    (re.compile(r"\brvp\b|aspartyl protease", re.IGNORECASE), "RVP"),
    (re.compile(r"dead[-/ ]?like helicase|dead/deah box helicase|\bdead\b", re.IGNORECASE), "DEAD_helicase"),
    (re.compile(r"helicase conserved c-terminal|\bhelicase_c\b", re.IGNORECASE), "Helicase_C"),
    (re.compile(r"snf2 family n-terminal", re.IGNORECASE), "SNF2_N"),
    (re.compile(r"chromo", re.IGNORECASE), "Chromo"),
    (re.compile(r"phd", re.IGNORECASE), "PHD"),
    (re.compile(r"\bset\b|set domain", re.IGNORECASE), "SET"),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Article-style domain co-occurrence networks")
    p.add_argument("--cooption-dir", default="/legserv/Temp/Carlos/EDF-TE/results/Co-option")
    p.add_argument("--out-dir", default="/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures")
    p.add_argument("--anchors", default="Chromo,PHD")
    p.add_argument("--true-only", action="store_true", help="Restrict to true_TE_coopted genes")
    p.add_argument("--min-node", type=int, default=4)
    p.add_argument("--min-edge", type=int, default=2)
    p.add_argument("--top-labels", type=int, default=20)
    return p.parse_args()


def clean_domain(x: str) -> str:
    s = re.sub(r"\s+", " ", str(x).strip()).strip(" .;")
    if s.lower() in BAD:
        return ""
    for patt, repl in DOMAIN_CANON_RULES:
        if patt.search(s):
            return repl
    return s


def load_true_gene_keys(class_glob: str) -> Set[str]:
    keys: Set[str] = set()
    for path in sorted(glob.glob(class_glob)):
        sp = os.path.basename(path).replace(".true_TE_coopted_genes.tsv", "")
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if "gene_id" not in df.columns:
            continue
        for gid in df["gene_id"].astype(str):
            g = gid.strip()
            if g:
                keys.add(f"{sp}::{g}")
    return keys


def load_gene_domains(fusion_glob: str, true_keys: Set[str], true_only: bool) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = defaultdict(set)
    for path in sorted(glob.glob(fusion_glob)):
        sp = os.path.basename(path).replace(".TE_host_fusions.tsv", "")
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if "gene_id" not in df.columns or "domain_name" not in df.columns:
            continue
        df["gene_id"] = df["gene_id"].astype(str).str.strip()
        df["domain_name"] = df["domain_name"].map(clean_domain)
        df = df[(df["gene_id"] != "") & (df["domain_name"] != "")]
        if df.empty:
            continue
        grp = df.groupby("gene_id", as_index=False)["domain_name"].unique()
        for _, r in grp.iterrows():
            key = f"{sp}::{r['gene_id']}"
            if true_only and key not in true_keys:
                continue
            out[key] = set(r["domain_name"])
    return out


def has_anchor(domains: Iterable[str], anchor: str) -> bool:
    a = anchor.lower()
    return any(a in d.lower() for d in domains)


def build_anchor_graph(gene_domains: Dict[str, Set[str]], anchor: str, min_node: int, min_edge: int):
    node_count = Counter()
    edge_count = Counter()
    n_anchor_genes = 0

    for _, doms in gene_domains.items():
        if not has_anchor(doms, anchor):
            continue
        n_anchor_genes += 1
        acc = sorted({d for d in doms if anchor.lower() not in d.lower()})
        if not acc:
            continue
        for d in acc:
            node_count[d] += 1
        for u, v in itertools.combinations(acc, 2):
            if u == v:
                continue
            if u > v:
                u, v = v, u
            edge_count[(u, v)] += 1

    keep_nodes = {d for d, n in node_count.items() if n >= min_node}
    G = nx.Graph(anchor=anchor, n_anchor_genes=n_anchor_genes)
    for d in keep_nodes:
        G.add_node(d, count=int(node_count[d]))
    for (u, v), w in edge_count.items():
        if w < min_edge:
            continue
        if u in keep_nodes and v in keep_nodes:
            G.add_edge(u, v, weight=int(w))

    if G.number_of_nodes() > 0:
        comms = list(nx.algorithms.community.asyn_lpa_communities(G, weight="weight", seed=42))
    else:
        comms = []
    node_to_comm = {}
    for i, c in enumerate(comms, start=1):
        for n in c:
            node_to_comm[n] = i
    nx.set_node_attributes(G, node_to_comm, "community")

    return G


def graph_layout(G: nx.Graph) -> Dict[str, Tuple[float, float]]:
    if G.number_of_nodes() == 0:
        return {}
    try:
        pos = nx.nx_pydot.pydot_layout(G, prog="neato")
        if pos:
            return pos
    except Exception:
        pass
    return nx.spring_layout(G, k=0.35, seed=42, iterations=300, weight="weight")


def export_graph_tables(G: nx.Graph, out_prefix: str) -> None:
    nodes = []
    for n, a in G.nodes(data=True):
        nodes.append({"domain": n, "count": int(a.get("count", 0)), "community": int(a.get("community", 0))})
    edges = []
    for u, v, a in G.edges(data=True):
        edges.append({"domain_u": u, "domain_v": v, "weight": int(a.get("weight", 0))})
    ndf = pd.DataFrame(nodes)
    if not ndf.empty:
        ndf = ndf.sort_values(["count", "domain"], ascending=[False, True])
    ndf.to_csv(out_prefix + ".nodes.tsv", sep="\t", index=False)
    edf = pd.DataFrame(edges)
    if not edf.empty:
        edf = edf.sort_values(["weight", "domain_u", "domain_v"], ascending=[False, True, True])
    edf.to_csv(out_prefix + ".edges.tsv", sep="\t", index=False)


def draw_community_halos(ax, G: nx.Graph, pos: Dict[str, Tuple[float, float]], comm_to_color: Dict[int, tuple]) -> None:
    for comm, color in comm_to_color.items():
        comm_nodes = [n for n, a in G.nodes(data=True) if a.get("community", 0) == comm]
        if len(comm_nodes) < 2:
            continue
        pts = np.array([pos[n] for n in comm_nodes], dtype=float)
        if len(comm_nodes) >= 3 and HAS_HULL:
            try:
                hull = ConvexHull(pts)
                poly = Polygon(pts[hull.vertices], closed=True, facecolor=color, edgecolor=color, alpha=0.10, linewidth=0.8, zorder=0)
                ax.add_patch(poly)
                continue
            except Exception:
                pass
        center = pts.mean(axis=0)
        radius = np.max(np.linalg.norm(pts - center, axis=1)) * 1.35
        circ = Circle(center, max(radius, 20.0), facecolor=color, edgecolor=color, alpha=0.08, linewidth=0.8, zorder=0)
        ax.add_patch(circ)


def add_size_legend(ax) -> None:
    vals = [5, 20, 50]
    sizes = [22 + 16 * np.sqrt(v) for v in vals]
    handles = [plt.scatter([], [], s=s, facecolor="#c9c9c9", edgecolor="black", linewidth=0.25) for s in sizes]
    ax.legend(handles, [str(v) for v in vals], title="Co-occurrences", loc="lower right", frameon=False, scatterpoints=1, fontsize=7, title_fontsize=8)


def plot_anchor(ax, G: nx.Graph, anchor: str, top_labels: int) -> None:
    if G.number_of_nodes() == 0:
        ax.set_title(f"{anchor}: no data")
        ax.axis("off")
        return

    pos = graph_layout(G)
    nodes = list(G.nodes())
    counts = np.array([G.nodes[n].get("count", 1) for n in nodes], dtype=float)
    comms = np.array([G.nodes[n].get("community", 0) for n in nodes], dtype=int)
    uniq_comms = sorted(set(comms))
    cmap = plt.get_cmap("tab20", max(3, len(uniq_comms) + 1))
    comm_to_color = {c: cmap(i % cmap.N) for i, c in enumerate(uniq_comms)}
    colors = [comm_to_color[c] for c in comms]

    draw_community_halos(ax, G, pos, comm_to_color)

    edge_w = np.array([G[u][v].get("weight", 1) for u, v in G.edges()], dtype=float)
    edge_lw = 0.25 + 1.8 * (edge_w / edge_w.max() if edge_w.size else edge_w)

    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_lw.tolist() if edge_w.size else 0.3, alpha=0.24, edge_color="#7a7a7a")
    node_size = 22 + 16 * np.sqrt(counts)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_size.tolist(), node_color=colors, linewidths=0.25, edgecolors="black", alpha=0.92)

    anchor_xy = np.mean(np.array(list(pos.values())), axis=0)
    ax.scatter([anchor_xy[0]], [anchor_xy[1]], s=320, c="#d9d9d9", edgecolors="black", linewidths=0.8, zorder=4)
    ax.text(anchor_xy[0], anchor_xy[1], anchor, ha="center", va="center", fontsize=8, fontweight="bold", zorder=5)

    # label top domains by anchor co-occurrence count
    top_nodes = sorted(nodes, key=lambda n: G.nodes[n].get("count", 0), reverse=True)[:top_labels]
    labels = {n: n for n in top_nodes}
    texts = nx.draw_networkx_labels(G, pos, labels=labels, font_size=7, ax=ax)
    for txt in texts.values():
        txt.set_bbox(dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))

    ax.set_title(f"Domain co-occurrence with {anchor}\n"
                 f"nodes={G.number_of_nodes()} edges={G.number_of_edges()} genes={G.graph.get('n_anchor_genes', 0)}")
    add_size_legend(ax)
    ax.axis("off")


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    anchors = [a.strip() for a in args.anchors.split(",") if a.strip()]
    true_glob = os.path.join(args.cooption_dir, "per_species", "classification", "*.true_TE_coopted_genes.tsv")
    fusion_glob = os.path.join(args.cooption_dir, "per_species", "fusions", "*.TE_host_fusions.tsv")

    true_keys = load_true_gene_keys(true_glob) if args.true_only else set()
    gene_domains = load_gene_domains(fusion_glob, true_keys=true_keys, true_only=args.true_only)

    graphs = []
    for a in anchors:
        G = build_anchor_graph(gene_domains, anchor=a, min_node=args.min_node, min_edge=args.min_edge)
        graphs.append((a, G))
        pref = os.path.join(args.out_dir, f"step7_article_network_{a}")
        export_graph_tables(G, pref)

    n = max(1, len(graphs))
    fig, axes = plt.subplots(1, n, figsize=(7.2 * n, 7.2), squeeze=False)
    axes = axes[0]
    for i, (a, G) in enumerate(graphs):
        plot_anchor(axes[i], G, a, top_labels=args.top_labels)

    mode = "true_only" if args.true_only else "all_fusions"
    fig.suptitle("Step 7: Accessory-domain co-occurrence networks by core anchor\n"
                 "Communities by label propagation; layout via NEATO (fallback: spring)", y=1.02)
    fig.tight_layout()

    out_pdf = os.path.join(args.out_dir, f"step7_article_style_networks_{mode}.pdf")
    out_png = os.path.join(args.out_dir, f"step7_article_style_networks_{mode}.png")
    out_pdf_v2 = os.path.join(args.out_dir, f"step7_article_style_networks_{mode}_v2.pdf")
    out_png_v2 = os.path.join(args.out_dir, f"step7_article_style_networks_{mode}_v2.png")
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf_v2)
    fig.savefig(out_png_v2, dpi=300)
    plt.close(fig)

    print(f"[OK] Genes used: {len(gene_domains)}")
    print(f"[OK] Wrote: {out_pdf}")
    print(f"[OK] Wrote: {out_png}")
    print(f"[OK] Wrote: {out_pdf_v2}")
    print(f"[OK] Wrote: {out_png_v2}")
    for a, G in graphs:
        print(f"[OK] {a}: nodes={G.number_of_nodes()} edges={G.number_of_edges()} genes={G.graph.get('n_anchor_genes', 0)}")


if __name__ == "__main__":
    main()
