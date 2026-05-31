#!/usr/bin/env python3

"""Step 6: Functional characterization of TE co-option domains.

Builds domain enrichment statistics for non-TE domains in true co-opted genes
versus background genes, then writes paper-ready TSV + PDF figures.
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import re
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


BAD = {"", "-", "\\N", "unknown", "nan", "none", "no description", "n/a", "na"}


def parse_args() -> argparse.Namespace:
	p = argparse.ArgumentParser(description="Domain enrichment in true TE co-option genes")
	p.add_argument("--cooption-dir", default="/legserv/Temp/Carlos/EDF-TE/results/Co-option")
	p.add_argument("--out-dir", default="/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures")
	p.add_argument("--top", type=int, default=30, help="Top enriched domains for PDF")
	return p.parse_args()


def clean_domain(x: str) -> str:
	s = re.sub(r"\s+", " ", str(x).strip()).strip(" .;")
	if s.lower() in BAD:
		return ""
	return s


def bh_fdr(pvals: Iterable[float]) -> List[float]:
	vals = list(pvals)
	n = len(vals)
	order = np.argsort(vals)
	adj = np.empty(n, dtype=float)
	prev = 1.0
	for rank, idx in enumerate(order[::-1], start=1):
		p = vals[idx]
		q = p * n / (n - rank + 1)
		prev = min(prev, q)
		adj[idx] = min(prev, 1.0)
	return adj.tolist()


def load_true_domains(class_glob: str) -> Tuple[Set[str], Dict[str, Set[str]]]:
	true_genes: Set[str] = set()
	dom_by_gene: Dict[str, Set[str]] = defaultdict(set)

	for path in sorted(glob.glob(class_glob)):
		sp = os.path.basename(path).replace(".true_TE_coopted_genes.tsv", "")
		df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
		if df.empty or "gene_id" not in df.columns:
			continue
		for _, r in df.iterrows():
			gid = str(r["gene_id"]).strip()
			if not gid:
				continue
			key = f"{sp}::{gid}"
			true_genes.add(key)
			doms = [clean_domain(d) for d in str(r.get("domains", "")).split(",")]
			for d in doms:
				if d:
					dom_by_gene[key].add(d)
	return true_genes, dom_by_gene


def load_background_non_te(dom_glob: str) -> Dict[str, Set[str]]:
	out: Dict[str, Set[str]] = defaultdict(set)
	for path in sorted(glob.glob(dom_glob)):
		sp = os.path.basename(path).replace(".domains.tsv", "")
		df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
		if "gene_id" not in df.columns or "domain_name" not in df.columns:
			continue

		if "domain_class" in df.columns:
			mask = ~df["domain_class"].isin(["TE_core", "TE_associated"])
			df = df.loc[mask].copy()

		# Vectorized cleanup + groupby is much faster than iterrows on large tables.
		df["domain_name"] = df["domain_name"].map(clean_domain)
		df["gene_id"] = df["gene_id"].astype(str).str.strip()
		df = df[(df["gene_id"] != "") & (df["domain_name"] != "")]
		if df.empty:
			continue

		grp = df.groupby("gene_id", as_index=False)["domain_name"].unique()
		for _, r in grp.iterrows():
			out[f"{sp}::{r['gene_id']}"] = set(r["domain_name"])
	return out


def run_enrichment(true_genes: Set[str], true_dom: Dict[str, Set[str]], bg_dom: Dict[str, Set[str]]) -> pd.DataFrame:
	universe = set(bg_dom.keys())
	true_u = sorted(universe & true_genes)
	bg_u = sorted(universe - true_genes)

	c_true = Counter()
	c_bg = Counter()

	for g in true_u:
		for d in bg_dom.get(g, set()):
			c_true[d] += 1
	for g in bg_u:
		for d in bg_dom.get(g, set()):
			c_bg[d] += 1

	n_true = len(true_u)
	n_bg = len(bg_u)
	rows = []
	for d in sorted(set(c_true) | set(c_bg)):
		a = int(c_true.get(d, 0))
		c = int(c_bg.get(d, 0))
		b = n_true - a
		d2 = n_bg - c
		if a + c == 0:
			continue
		odds, p = fisher_exact([[a, b], [c, d2]], alternative="two-sided")
		rows.append(
			{
				"domain": d,
				"n_true_with_domain": a,
				"n_bg_with_domain": c,
				"pct_true": (a / n_true * 100.0) if n_true else 0.0,
				"pct_bg": (c / n_bg * 100.0) if n_bg else 0.0,
				"odds_ratio": float(odds) if np.isfinite(odds) else np.nan,
				"p_value": float(p),
				"n_true_genes": n_true,
				"n_bg_genes": n_bg,
			}
		)

	df = pd.DataFrame(rows)
	if df.empty:
		return df
	df["fdr_bh"] = bh_fdr(df["p_value"].tolist())
	df["effect"] = np.where(df["odds_ratio"] > 1, "enriched", "depleted")
	df = df.sort_values(["fdr_bh", "odds_ratio", "n_true_with_domain"], ascending=[True, False, False])
	return df


def plot_top_enriched(df: pd.DataFrame, out_pdf: str, top: int) -> None:
	sub = df[(df["odds_ratio"] > 1) & (df["n_true_with_domain"] >= 5)].copy().head(top)
	if sub.empty:
		return
	sub = sub.iloc[::-1].copy()
	sub["score"] = -np.log10(np.clip(sub["fdr_bh"].values, 1e-300, 1.0))

	fig, ax = plt.subplots(figsize=(11, max(8, len(sub) * 0.28)))
	bars = ax.barh(sub["domain"], sub["score"], color="#2a9d8f", alpha=0.85)
	ax.set_xlabel("-log10(FDR)")
	ax.set_ylabel("Domain")
	ax.set_title("Step 6: Enriched non-TE domains in true TE-coopted genes")
	ax.grid(axis="x", linestyle="--", alpha=0.3)
	for b, pct in zip(bars, sub["pct_true"]):
		ax.text(b.get_width() + 0.05, b.get_y() + b.get_height() / 2, f"{pct:.1f}%", va="center", fontsize=8)
	fig.tight_layout()
	fig.savefig(out_pdf)
	plt.close(fig)


def plot_true_vs_bg(df: pd.DataFrame, out_pdf: str, top: int) -> None:
	sub = df[(df["odds_ratio"] > 1) & (df["n_true_with_domain"] >= 5)].copy().head(top)
	if sub.empty:
		return
	sub = sub.iloc[::-1].copy()

	y = np.arange(len(sub))
	h = 0.38
	fig, ax = plt.subplots(figsize=(11, max(8, len(sub) * 0.28)))
	ax.barh(y - h / 2, sub["pct_true"], height=h, color="#e76f51", label="True co-opted")
	ax.barh(y + h / 2, sub["pct_bg"], height=h, color="#457b9d", label="Background")
	ax.set_yticks(y)
	ax.set_yticklabels(sub["domain"])
	ax.set_xlabel("% genes carrying domain")
	ax.set_title("Step 6: Domain prevalence in co-opted genes vs background")
	ax.grid(axis="x", linestyle="--", alpha=0.3)
	ax.legend(frameon=False)
	fig.tight_layout()
	fig.savefig(out_pdf)
	plt.close(fig)


def main() -> None:
	args = parse_args()
	os.makedirs(args.out_dir, exist_ok=True)

	true_glob = os.path.join(args.cooption_dir, "per_species", "classification", "*.true_TE_coopted_genes.tsv")
	dom_glob = os.path.join(args.cooption_dir, "per_species", "domains", "*.domains.tsv")

	true_genes, _true_dom = load_true_domains(true_glob)
	bg_dom = load_background_non_te(dom_glob)
	enr = run_enrichment(true_genes, _true_dom, bg_dom)

	out_tsv = os.path.join(args.out_dir, "step6_domain_enrichment.tsv")
	enr.to_csv(out_tsv, sep="\t", index=False)

	out_pdf1 = os.path.join(args.out_dir, "step6_domain_enrichment_top.pdf")
	out_pdf2 = os.path.join(args.out_dir, "step6_domain_prevalence_true_vs_bg.pdf")
	if not enr.empty:
		plot_top_enriched(enr, out_pdf1, top=args.top)
		plot_true_vs_bg(enr, out_pdf2, top=args.top)

	print(f"[OK] Wrote: {out_tsv}")
	if os.path.exists(out_pdf1):
		print(f"[OK] Wrote: {out_pdf1}")
	if os.path.exists(out_pdf2):
		print(f"[OK] Wrote: {out_pdf2}")


if __name__ == "__main__":
	main()
