#!/usr/bin/env python3

"""Generate Fig.5-style candidate panels for TE-host cooption results.

Outputs:
1) Co-occurrence bubble panels for anchor domains (default: Chromo, PHD)
2) Architecture-by-clade bubble matrix for top recurrent fusion architectures
3) TSV summaries backing both figures
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import pandas as pd

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False


BAD_DOMAIN_NAMES = {
    "",
    "-",
    "\\N",
    "no description",
    "unknown",
    "nan",
    "none",
}


# Collapse recurring InterPro/CDD/Pfam textual variants into stable labels.
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


@dataclass
class GeneDomains:
    species: str
    gene_id: str
    domains: Set[str]
    te_core_domains: Set[str]
    te_assoc_domains: Set[str]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build Fig.5-style cooption candidate figures")
    p.add_argument(
        "--cooption-dir",
        default="/legserv/Temp/Carlos/EDF-TE/results/Co-option",
        help="Root Co-option results directory",
    )
    p.add_argument(
        "--out-dir",
        default="/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures",
        help="Output directory for figures and TSVs",
    )
    p.add_argument(
        "--anchors",
        default="Chromo,PHD",
        help="Comma-separated anchor domain keywords",
    )
    p.add_argument(
        "--top-partners",
        type=int,
        default=18,
        help="Top partner domains per anchor to plot",
    )
    p.add_argument(
        "--top-architectures",
        type=int,
        default=20,
        help="Top architectures to show in architecture matrix",
    )
    return p.parse_args()


def normalize_domain_name(name: str) -> str:
    x = str(name).strip()
    if x.lower() in BAD_DOMAIN_NAMES:
        return ""
    x = re.sub(r"\s+", " ", x)
    x = x.strip(" .;")
    if not x:
        return ""
    for patt, repl in DOMAIN_CANON_RULES:
        if patt.search(x):
            return repl
    return x


def species_from_path(path: str) -> str:
    base = os.path.basename(path)
    return base.replace(".TE_host_fusions.tsv", "")


def load_te_domain_master(master_path: str) -> Tuple[Set[str], Set[str]]:
    te_core: Set[str] = set()
    te_assoc: Set[str] = set()
    if not os.path.exists(master_path):
        return te_core, te_assoc

    df = pd.read_csv(master_path, sep="\t", dtype=str).fillna("")
    for _, row in df.iterrows():
        status = row.get("TE_status", "").strip()
        did = row.get("domain_id", "").strip()
        dname = normalize_domain_name(row.get("domain_name", ""))
        if status == "TE_core":
            if did:
                te_core.add(did)
            if dname:
                te_core.add(dname)
        elif status == "TE_associated":
            if did:
                te_assoc.add(did)
            if dname:
                te_assoc.add(dname)
    return te_core, te_assoc


def load_gene_domain_map(fusion_glob: str, te_core_ref: Set[str], te_assoc_ref: Set[str]) -> List[GeneDomains]:
    items: List[GeneDomains] = []
    for path in sorted(glob.glob(fusion_glob)):
        sp = species_from_path(path)
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if "gene_id" not in df.columns:
            continue

        gene_domains: Dict[str, Set[str]] = defaultdict(set)
        gene_te_core: Dict[str, Set[str]] = defaultdict(set)
        gene_te_assoc: Dict[str, Set[str]] = defaultdict(set)

        for _, row in df.iterrows():
            gid = str(row.get("gene_id", "")).strip()
            if not gid:
                continue
            dname = normalize_domain_name(row.get("domain_name", ""))
            did = normalize_domain_name(row.get("domain_id", ""))
            if dname:
                gene_domains[gid].add(dname)

            # mark TE class from either ID or name
            if did in te_core_ref or dname in te_core_ref:
                gene_te_core[gid].add(dname or did)
            elif did in te_assoc_ref or dname in te_assoc_ref:
                gene_te_assoc[gid].add(dname or did)

        for gid, domains in gene_domains.items():
            if not domains:
                continue
            items.append(
                GeneDomains(
                    species=sp,
                    gene_id=gid,
                    domains=domains,
                    te_core_domains=gene_te_core.get(gid, set()),
                    te_assoc_domains=gene_te_assoc.get(gid, set()),
                )
            )
    return items


def match_anchor(domains: Iterable[str], anchor: str) -> bool:
    needle = anchor.lower()
    return any(needle in d.lower() for d in domains)


def partner_class(domain: str, te_core_ref: Set[str], te_assoc_ref: Set[str]) -> str:
    if domain in te_core_ref:
        return "TE_core"
    if domain in te_assoc_ref:
        return "TE_associated"
    return "host_like"


def summarize_anchor_partners(
    genes: Sequence[GeneDomains],
    anchors: Sequence[str],
    te_core_ref: Set[str],
    te_assoc_ref: Set[str],
) -> pd.DataFrame:
    rows: List[dict] = []
    for anchor in anchors:
        counter = Counter()
        for g in genes:
            if not match_anchor(g.domains, anchor):
                continue
            for d in g.domains:
                if anchor.lower() in d.lower():
                    continue
                counter[d] += 1
        for dom, n in counter.items():
            rows.append(
                {
                    "anchor": anchor,
                    "partner_domain": dom,
                    "count": n,
                    "partner_class": partner_class(dom, te_core_ref, te_assoc_ref),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["anchor", "count"], ascending=[True, False])


def architecture_signature(g: GeneDomains, max_terms: int = 4) -> str:
    ordered = sorted(g.domains)
    if len(ordered) > max_terms:
        ordered = ordered[:max_terms]
    return " + ".join(ordered)


def summarize_architecture_matrix(genes: Sequence[GeneDomains], top_n: int) -> pd.DataFrame:
    arch_counter = Counter(architecture_signature(g) for g in genes)
    top_arch = [a for a, _ in arch_counter.most_common(top_n)]

    rows: List[dict] = []
    for g in genes:
        arch = architecture_signature(g)
        if arch not in top_arch:
            continue
        clade = g.species.split("_", 1)[0]
        rows.append({"architecture": arch, "species": g.species, "clade": clade})

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    agg = (
        df.groupby(["architecture", "clade"], as_index=False)
        .agg(n_genes=("species", "size"), n_species=("species", "nunique"))
    )

    order_map = {a: i for i, a in enumerate(top_arch)}
    agg["arch_order"] = agg["architecture"].map(order_map)
    agg = agg.sort_values(["arch_order", "n_genes"], ascending=[True, False])
    return agg


def plot_anchor_bubbles(df: pd.DataFrame, anchors: Sequence[str], top_k: int, out_png: str, out_pdf: str) -> None:
    if not HAS_MPL:
        return
    if df.empty:
        return
    color_map = {"TE_core": "#f28e2b", "TE_associated": "#e15759", "host_like": "#4e79a7"}
    fig, axes = plt.subplots(1, len(anchors), figsize=(7 * len(anchors), 7), squeeze=False)
    axes = axes[0]

    for i, anchor in enumerate(anchors):
        ax = axes[i]
        sub = df[df["anchor"] == anchor].head(top_k).copy()
        if sub.empty:
            ax.set_title(f"{anchor}: no data")
            ax.axis("off")
            continue

        sub = sub.sort_values("count", ascending=True).reset_index(drop=True)
        y = range(len(sub))
        size = (sub["count"] ** 0.8) * 9
        colors = [color_map.get(c, "#999999") for c in sub["partner_class"]]

        ax.scatter(sub["count"], y, s=size, c=colors, alpha=0.82, edgecolor="black", linewidth=0.3)
        ax.set_yticks(list(y))
        ax.set_yticklabels(sub["partner_domain"], fontsize=8)
        ax.set_xlabel("Number of co-occurrences")
        ax.set_title(f"Domain co-occurrence with {anchor}")
        ax.grid(axis="x", alpha=0.25, linestyle="--")

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markeredgecolor="black", markersize=9, label=lab)
        for lab, c in color_map.items()
    ]
    fig.legend(handles=handles, labels=[h.get_label() for h in handles], loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_architecture_matrix(df: pd.DataFrame, out_png: str, out_pdf: str) -> None:
    if not HAS_MPL:
        return
    if df.empty:
        return
    clades = sorted(df["clade"].unique())
    arch = list(dict.fromkeys(df.sort_values("arch_order")["architecture"].tolist()))

    fig, ax = plt.subplots(figsize=(max(10, len(clades) * 1.0), max(7, len(arch) * 0.38)))

    x_pos = {c: i for i, c in enumerate(clades)}
    y_pos = {a: i for i, a in enumerate(arch)}

    for _, r in df.iterrows():
        x = x_pos[r["clade"]]
        y = y_pos[r["architecture"]]
        n_genes = int(r["n_genes"])
        n_species = int(r["n_species"])
        size = (n_genes ** 0.75) * 14
        color = "#2a9d8f" if n_species > 1 else "#bdbdbd"
        ax.scatter(x, y, s=size, c=color, alpha=0.85, edgecolor="black", linewidth=0.35)

    ax.set_xticks(list(x_pos.values()))
    ax.set_xticklabels(clades, rotation=45, ha="right")
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels(arch, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Clade")
    ax.set_ylabel("Fusion architecture (top recurrent)")
    ax.set_title("TE-host fusion architecture recurrence by clade")
    ax.grid(alpha=0.22, linestyle="--")

    legend_elements = [
        plt.Line2D([0], [0], marker="o", color="w", label="Single species", markerfacecolor="#bdbdbd", markeredgecolor="black", markersize=8),
        plt.Line2D([0], [0], marker="o", color="w", label=">1 species", markerfacecolor="#2a9d8f", markeredgecolor="black", markersize=8),
    ]
    ax.legend(handles=legend_elements, frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    master = os.path.join(args.cooption_dir, "config", "TE_domains_master.tsv")
    te_core_ref, te_assoc_ref = load_te_domain_master(master)

    fusion_glob = os.path.join(args.cooption_dir, "per_species", "fusions", "*.TE_host_fusions.tsv")
    genes = load_gene_domain_map(fusion_glob, te_core_ref, te_assoc_ref)
    if not genes:
        raise SystemExit("No fusion rows detected in per_species/fusions/*.TE_host_fusions.tsv")

    anchors = [a.strip() for a in args.anchors.split(",") if a.strip()]
    partner_df = summarize_anchor_partners(genes, anchors, te_core_ref, te_assoc_ref)
    arch_df = summarize_architecture_matrix(genes, top_n=args.top_architectures)

    partner_tsv = os.path.join(args.out_dir, "fig5_candidate_anchor_cooccurrence.tsv")
    arch_tsv = os.path.join(args.out_dir, "fig5_candidate_architecture_matrix.tsv")
    partner_df.to_csv(partner_tsv, sep="\t", index=False)
    arch_df.to_csv(arch_tsv, sep="\t", index=False)

    panel_c_png = os.path.join(args.out_dir, "fig5_candidate_panelC_anchor_cooccurrence.png")
    panel_e_png = os.path.join(args.out_dir, "fig5_candidate_panelE_architecture_matrix.png")
    panel_c_pdf = os.path.join(args.out_dir, "fig5_candidate_panelC_anchor_cooccurrence.pdf")
    panel_e_pdf = os.path.join(args.out_dir, "fig5_candidate_panelE_architecture_matrix.pdf")
    plot_anchor_bubbles(partner_df, anchors, top_k=args.top_partners, out_png=panel_c_png, out_pdf=panel_c_pdf)
    plot_architecture_matrix(arch_df, out_png=panel_e_png, out_pdf=panel_e_pdf)

    print(f"[OK] Genes processed: {len(genes)}")
    print(f"[OK] Wrote: {partner_tsv}")
    print(f"[OK] Wrote: {arch_tsv}")
    if HAS_MPL:
        print(f"[OK] Wrote: {panel_c_png}")
        print(f"[OK] Wrote: {panel_e_png}")
        print(f"[OK] Wrote: {panel_c_pdf}")
        print(f"[OK] Wrote: {panel_e_pdf}")
    else:
        print("[WARN] matplotlib not available; PNG figures were skipped")
        print(f"[INFO] Install matplotlib, then rerun to render: {panel_c_png}, {panel_e_png}, {panel_c_pdf} and {panel_e_pdf}")


if __name__ == "__main__":
    main()
