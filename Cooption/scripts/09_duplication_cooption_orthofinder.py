#!/usr/bin/env python3

"""Step 8-9: Co-option and duplication using OrthoFinder outputs.

Outputs TSV + PDF summaries of:
- co-opted genes by duplication class (single/duplicated/multicopy)
- expected background distribution from orthogroup copy counts
- mixed-fate orthogroups (co-opted and non-co-opted copies within species)
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Co-option vs duplication analysis")
    p.add_argument("--cooption-dir", default="/legserv/Temp/Carlos/EDF-TE/results/Co-option")
    p.add_argument("--orthogroups", default="/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups/Orthogroups.tsv")
    p.add_argument("--genecount", default="/legserv/Temp/Carlos/EDF-TE/results/Orthofinder_run/Orthogroups/Orthogroups.GeneCount.tsv")
    p.add_argument("--out-dir", default="/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures")
    return p.parse_args()


def dup_class(n: int) -> str:
    if n <= 1:
        return "single_copy"
    if n == 2:
        return "duplicated"
    return "multicopy"


def species_from_true(path: str) -> str:
    return os.path.basename(path).replace(".true_TE_coopted_genes.tsv", "")


def gene_aliases(token: str) -> Set[str]:
    t = str(token).strip()
    if not t:
        return set()
    out = {t}
    if "|" in t:
        parts = [p.strip() for p in t.split("|") if p.strip()]
        out.update(parts)
        out.add(parts[-1])
    out.add(t.split()[0])
    out.add(t.replace('"', ""))
    if "." in t:
        out.add(t.rsplit(".", 1)[0])
    for m in re.findall(r"\d+", t):
        out.add(m)
    return {x for x in out if x}


def load_coopted_genes(true_glob: str) -> Tuple[Dict[str, Set[str]], Dict[str, Dict[str, str]]]:
    coopted: Dict[str, Set[str]] = defaultdict(set)
    alias_map: Dict[str, Dict[str, str]] = defaultdict(dict)

    for path in sorted(glob_paths(true_glob)):
        sp = species_from_true(path)
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if "gene_id" not in df.columns:
            continue
        for gid in df["gene_id"].astype(str):
            g = gid.strip()
            if not g:
                continue
            coopted[sp].add(g)
            for a in gene_aliases(g):
                alias_map[sp][a] = g
    return coopted, alias_map


def glob_paths(pattern: str) -> List[str]:
    import glob

    return glob.glob(pattern)


def load_copy_classes(genecount_tsv: str) -> Tuple[pd.DataFrame, Dict[str, Counter]]:
    df = pd.read_csv(genecount_tsv, sep="\t", dtype=str)
    if "Orthogroup" not in df.columns:
        raise ValueError("Orthogroups.GeneCount.tsv missing Orthogroup column")

    species_cols = [c for c in df.columns if c not in {"Orthogroup", "Total"}]
    clean_species = {c: c.replace(".aa", "") for c in species_cols}

    # Vectorized melt — much faster than iterrows for large matrices
    num_df = df[["Orthogroup"] + species_cols].copy()
    for c in species_cols:
        num_df[c] = pd.to_numeric(num_df[c], errors="coerce").fillna(0).astype(int)

    long_df = num_df.melt(id_vars="Orthogroup", var_name="col", value_name="copy_n")
    long_df["species"] = long_df["col"].map(clean_species)
    long_df["copy_n"] = long_df["copy_n"].astype(int)
    long_df["dup_class"] = long_df["copy_n"].map(dup_class)
    long_df = long_df.rename(columns={"Orthogroup": "orthogroup"})[["orthogroup", "species", "copy_n", "dup_class"]]

    bg_counts: Dict[str, Counter] = defaultdict(Counter)
    for row in long_df[long_df["copy_n"] > 0][["species", "dup_class", "copy_n"]].itertuples(index=False, name=None):
        bg_counts[row[0]][row[1]] += row[2]

    return long_df, bg_counts


def map_coopted_to_og(
    orthogroups_tsv: str,
    copy_df: pd.DataFrame,
    coopted: Dict[str, Set[str]],
    alias_map: Dict[str, Dict[str, str]],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    species = sorted(coopted.keys())
    species_col_to_name = {}
    for sp in species:
        species_col_to_name[f"{sp}.aa"] = sp
        species_col_to_name[sp] = sp

    copy_lookup = {(row[0], row[1]): int(row[2]) for row in copy_df[["orthogroup", "species", "copy_n"]].itertuples(index=False, name=None)}

    match_rows = []
    mixed_rows = []

    with open(orthogroups_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            og = row["Orthogroup"]
            for col, sp in species_col_to_name.items():
                if col not in row:
                    continue
                cell = str(row[col]).strip()
                if not cell:
                    continue

                genes = [x.strip() for x in cell.split(",") if x.strip()]
                if not genes:
                    continue

                matched = []
                for g in genes:
                    aliases = gene_aliases(g)
                    found = None
                    amap = alias_map.get(sp, {})
                    for a in aliases:
                        if a in amap:
                            found = amap[a]
                            break
                    if found is not None:
                        matched.append(found)

                n_total = len(genes)
                n_match = len(set(matched))
                if n_match > 0:
                    cn = copy_lookup.get((og, sp), n_total)
                    dc = dup_class(cn)
                    for g in sorted(set(matched)):
                        match_rows.append(
                            {
                                "species": sp,
                                "orthogroup": og,
                                "gene_id": g,
                                "copy_n_species_og": int(cn),
                                "dup_class": dc,
                            }
                        )

                if n_match > 0 and n_match < n_total:
                    mixed_rows.append(
                        {
                            "species": sp,
                            "orthogroup": og,
                            "n_total_genes_in_species_og": n_total,
                            "n_coopted_genes_in_species_og": n_match,
                            "n_noncoopted_genes_in_species_og": n_total - n_match,
                        }
                    )

    return pd.DataFrame(match_rows), pd.DataFrame(mixed_rows)


def summarize(
    matched_df: pd.DataFrame,
    bg_counts: Dict[str, Counter],
    coopted: Dict[str, Set[str]],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if matched_df.empty:
        empty = pd.DataFrame()
        return empty, empty, empty

    co_counts_sp = (
        matched_df.groupby(["species", "dup_class"], as_index=False)
        .agg(n_coopted_genes=("gene_id", "count"))
        .sort_values(["species", "dup_class"])
    )

    bg_rows = []
    for sp, c in bg_counts.items():
        for dc in ["single_copy", "duplicated", "multicopy"]:
            bg_rows.append({"species": sp, "dup_class": dc, "n_bg_genes": int(c.get(dc, 0))})
    bg_df = pd.DataFrame(bg_rows)

    merged = pd.merge(co_counts_sp, bg_df, on=["species", "dup_class"], how="outer").fillna(0)

    co_global = matched_df["dup_class"].value_counts().reindex(["single_copy", "duplicated", "multicopy"]).fillna(0).astype(int)
    bg_global = Counter()
    for c in bg_counts.values():
        for k in ["single_copy", "duplicated", "multicopy"]:
            bg_global[k] += int(c.get(k, 0))

    rows = []
    for dc in ["single_copy", "duplicated", "multicopy"]:
        rows.append({"dup_class": dc, "n_coopted": int(co_global[dc]), "n_background": int(bg_global[dc])})
    global_df = pd.DataFrame(rows)

    return co_counts_sp, merged, global_df


def stats_table(global_df: pd.DataFrame) -> pd.DataFrame:
    if global_df.empty:
        return pd.DataFrame()
    order = ["single_copy", "duplicated", "multicopy"]
    a = [int(global_df.loc[global_df["dup_class"] == x, "n_coopted"].iloc[0]) for x in order]
    b = [int(global_df.loc[global_df["dup_class"] == x, "n_background"].iloc[0]) for x in order]

    chi2, p, _, _ = chi2_contingency([a, b])
    rows = [{"test": "chi2_2x3", "statistic": chi2, "p_value": p, "details": "coopted vs background across 3 duplication classes"}]

    # Pairwise OR multicopy vs single; duplicated vs single
    for x in ["duplicated", "multicopy"]:
        i = order.index(x)
        j = order.index("single_copy")
        table = [[a[i], a[j]], [b[i], b[j]]]
        orr, fp = fisher_exact(table)
        rows.append({"test": f"fisher_{x}_vs_single", "statistic": float(orr), "p_value": float(fp), "details": str(table)})

    return pd.DataFrame(rows)


def plot_distribution(global_df: pd.DataFrame, out_pdf: str) -> None:
    if global_df.empty:
        return
    sub = global_df.copy()
    sub["pct_coopted"] = sub["n_coopted"] / max(1, sub["n_coopted"].sum()) * 100
    sub["pct_bg"] = sub["n_background"] / max(1, sub["n_background"].sum()) * 100

    x = np.arange(len(sub))
    w = 0.38
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.bar(x - w / 2, sub["pct_coopted"], width=w, label="Co-opted", color="#e76f51")
    ax.bar(x + w / 2, sub["pct_bg"], width=w, label="Background", color="#457b9d")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["dup_class"], rotation=15)
    ax.set_ylabel("% genes")
    ax.set_title("Step 8: Duplication class distribution of co-opted genes")
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_mixed(mixed_df: pd.DataFrame, out_pdf: str) -> None:
    if mixed_df.empty:
        return
    tmp = mixed_df.copy()
    tmp["clade"] = tmp["species"].str.split("_", n=1).str[0]
    cl = tmp.groupby("clade", as_index=False).agg(n_mixed_ogs=("orthogroup", "nunique"), n_species=("species", "nunique"))
    cl["mixed_ogs_per_species"] = cl["n_mixed_ogs"] / cl["n_species"].replace(0, np.nan)
    cl = cl.sort_values("mixed_ogs_per_species", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.bar(cl["clade"], cl["mixed_ogs_per_species"], color="#2a9d8f", alpha=0.85)
    ax.set_ylabel("Mixed OGs per species")
    ax.set_xlabel("Clade")
    ax.set_title("Step 9: Mixed-fate orthogroups (TE-derived + TE-free copies)")
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    true_glob = os.path.join(args.cooption_dir, "per_species", "classification", "*.true_TE_coopted_genes.tsv")
    coopted, alias_map = load_coopted_genes(true_glob)
    copy_df, bg_counts = load_copy_classes(args.genecount)
    matched_df, mixed_df = map_coopted_to_og(args.orthogroups, copy_df, coopted, alias_map)

    co_counts_sp, merged_sp, global_df = summarize(matched_df, bg_counts, coopted)
    st = stats_table(global_df)

    out_mapped = os.path.join(args.out_dir, "step8_coopted_genes_with_dup_class.tsv")
    out_merged = os.path.join(args.out_dir, "step8_duplication_distribution_by_species.tsv")
    out_global = os.path.join(args.out_dir, "step8_duplication_distribution_global.tsv")
    out_stats = os.path.join(args.out_dir, "step8_duplication_stats.tsv")
    out_mixed = os.path.join(args.out_dir, "step9_mixed_fate_orthogroups.tsv")

    matched_df.to_csv(out_mapped, sep="\t", index=False)
    merged_sp.to_csv(out_merged, sep="\t", index=False)
    global_df.to_csv(out_global, sep="\t", index=False)
    st.to_csv(out_stats, sep="\t", index=False)
    mixed_df.to_csv(out_mixed, sep="\t", index=False)

    pdf_dist = os.path.join(args.out_dir, "step8_duplication_distribution_global.pdf")
    pdf_mixed = os.path.join(args.out_dir, "step9_mixed_fate_orthogroups_by_clade.pdf")
    plot_distribution(global_df, pdf_dist)
    plot_mixed(mixed_df, pdf_mixed)

    n_co = sum(len(v) for v in coopted.values())
    n_match = 0 if matched_df.empty else matched_df["gene_id"].nunique()
    print(f"[OK] Co-opted genes input: {n_co}")
    print(f"[OK] Co-opted genes mapped to orthogroups: {n_match}")
    print(f"[OK] Wrote: {out_mapped}")
    print(f"[OK] Wrote: {out_merged}")
    print(f"[OK] Wrote: {out_global}")
    print(f"[OK] Wrote: {out_stats}")
    print(f"[OK] Wrote: {out_mixed}")
    if os.path.exists(pdf_dist):
        print(f"[OK] Wrote: {pdf_dist}")
    if os.path.exists(pdf_mixed):
        print(f"[OK] Wrote: {pdf_mixed}")


if __name__ == "__main__":
    main()
