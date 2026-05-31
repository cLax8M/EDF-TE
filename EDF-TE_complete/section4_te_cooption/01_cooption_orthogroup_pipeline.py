#!/usr/bin/env python3
"""
01_cooption_orthogroup_pipeline.py  —  TE-derived gene co-option analysis.

Identifies orthogroups with evidence of TE-domain origin and traces their
evolutionary history across the fungal/protistan phylogeny.

Pipeline:
  1. Screen protein orthogroups for TE-related Pfam/InterPro domains
     (transposase, RT, integrase, chromo, RNase H, etc.)
  2. Confirm TE origin by reciprocal BLAST against TE consensus library
  3. Classify co-opted genes by TE class and domain architecture
  4. Score per-clade presence/absence and infer minimum gain/loss events
  5. Produce domain-frequency bubble plot and orthogroup heatmap

Inputs:
  --orthogroup_seqs  Directory with one FASTA per orthogroup
  --ipr_dir          Directory with per-species IPR annotation TSVs
  --te_lib           TE consensus library FASTA
  --species_tree     Newick species tree
  --clade_map        Two-column TSV: species_code<TAB>clade_tag
  --out_dir          Output directory
  --blast_pid        Minimum % identity for TE-origin confirmation (default: 30)
  --blast_cov        Minimum query coverage for TE-origin confirmation (default: 60)
  --te_domains       Comma-separated IPR IDs to screen (or use default set)
  --threads          DIAMOND threads (default: 8)

Outputs:
  cooption_orthogroups.tsv
  cooption_domain_summary.tsv
  figures/domain_bubble_plot.pdf
  figures/cooption_heatmap.pdf
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Default TE-related IPR domain signatures
TE_DOMAINS_DEFAULT = {
    "IPR001584": "Integrase_catalytic",
    "IPR002104": "Chromo_domain",
    "IPR002710": "Pol_plant",
    "IPR003308": "Transposase_Mu",
    "IPR005093": "Gag_plant_retrovirus",
    "IPR005484": "RT_thumb",
    "IPR007721": "DDE_superfam",
    "IPR008906": "HAT_dimerisation",
    "IPR010994": "RNase_HI_RT",
    "IPR013149": "ADK_C",
    "IPR017716": "Tc1_transposase",
    "IPR021967": "Transposase_MULE",
    "IPR036397": "Gag_zinc_finger",
}

CLADE_ORDER = [
    "HO", "NU", "AP", "CR", "MI", "BL", "OL", "CH", "NE",
    "BS", "EN", "ZO", "KX", "GL", "MO", "MU", "BA", "AS",
]

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orthogroup_seqs", required=True)
    p.add_argument("--ipr_dir",         required=True)
    p.add_argument("--te_lib",          required=True)
    p.add_argument("--species_tree",    required=True)
    p.add_argument("--clade_map",       required=True)
    p.add_argument("--out_dir",         default="cooption_out")
    p.add_argument("--blast_pid",       type=float, default=30.0)
    p.add_argument("--blast_cov",       type=float, default=60.0)
    p.add_argument("--te_domains",      default=None,
                   help="Comma-separated IPR IDs to flag as TE-related")
    p.add_argument("--threads",         type=int, default=8)
    return p.parse_args()


def load_ipr_domain_presence(ipr_dir: Path, species_list: list[str],
                              te_ipr: set[str]) -> pd.DataFrame:
    """Build gene→domain table for TE-related domains."""
    rows = []
    for sp in species_list:
        files = list(Path(ipr_dir).glob(f"{sp}*.ipr.tab")) + \
                list(Path(ipr_dir).glob(f"{sp}*.tsv"))
        if not files:
            continue
        df = pd.read_csv(files[0], sep="\t", header=None, usecols=[0, 3, 4, 5],
                         names=["gene", "db", "ipr_id", "ipr_desc"]).dropna()
        sub = df[df["ipr_id"].isin(te_ipr)]
        if not sub.empty:
            sub = sub.copy()
            sub["species"] = sp
            rows.append(sub)
    if not rows:
        return pd.DataFrame(columns=["gene", "ipr_id", "ipr_desc", "species"])
    return pd.concat(rows, ignore_index=True)


def build_diamond_db(te_lib: str, db_path: str) -> None:
    subprocess.run(
        ["diamond", "makedb", "--in", te_lib, "--db", db_path],
        check=True, capture_output=True)


def blast_orthogroup_vs_te(og_fasta: Path, db: str, pid: float,
                             cov: float, threads: int) -> bool:
    """Return True if any protein in the orthogroup BLAST-hits the TE library."""
    res = subprocess.run(
        ["diamond", "blastp", "-q", str(og_fasta), "--db", db,
         "--outfmt", "6", "--id", str(pid), "--query-cover", str(cov),
         "--max-target-seqs", "1", "--threads", str(threads),
         "--quiet"],
        capture_output=True, text=True)
    return bool(res.stdout.strip())


def screen_orthogroups(og_dir: Path, te_ipr_genes: pd.DataFrame,
                        db: str, pid: float, cov: float, threads: int,
                        species_list: list[str]) -> pd.DataFrame:
    """Screen all orthogroups for TE domain evidence + BLAST confirmation."""
    # Map gene → species from IPR table
    gene2sp = dict(zip(te_ipr_genes["gene"], te_ipr_genes["species"]))
    gene2ipr = te_ipr_genes.groupby("gene")["ipr_id"].apply(list).to_dict()

    rows = []
    og_fastas = sorted(og_dir.glob("OG*.fa")) + sorted(og_dir.glob("OG*.fasta"))

    for og_fa in og_fastas:
        og_id = og_fa.stem

        # Check if any sequence in OG has a TE domain
        ipr_hits = []
        with open(og_fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    gid = line[1:].split("|")[-1].strip()
                    if gid in gene2ipr:
                        ipr_hits.extend(gene2ipr[gid])

        if not ipr_hits:
            continue

        # Confirm with BLAST
        blast_confirmed = blast_orthogroup_vs_te(og_fa, db, pid, cov, threads)

        rows.append({
            "orthogroup":       og_id,
            "te_domains":       ",".join(sorted(set(ipr_hits))),
            "blast_confirmed":  blast_confirmed,
            "n_sequences":      sum(1 for l in open(og_fa) if l.startswith(">")),
        })

    return pd.DataFrame(rows)


def add_clade_presence(og_table: pd.DataFrame, orthogroup_seqs_dir: Path,
                        clade_map: dict) -> pd.DataFrame:
    """Annotate which clades have each co-opted orthogroup."""
    clade_cols = {c: [] for c in CLADE_ORDER}

    for _, row in og_table.iterrows():
        og_id = row["orthogroup"]
        fa = next(iter(Path(orthogroup_seqs_dir).glob(f"{og_id}*.fa*")), None)
        present_clades = set()
        if fa:
            with open(fa) as fh:
                for line in fh:
                    if line.startswith(">"):
                        sp = line[1:].split("|")[0].strip()
                        c  = clade_map.get(sp)
                        if c:
                            present_clades.add(c)
        for c in CLADE_ORDER:
            clade_cols[c].append(c in present_clades)

    for c in CLADE_ORDER:
        og_table[c] = clade_cols[c]
    return og_table


# =============================================================================
# Figures
# =============================================================================

def plot_domain_bubble(og_table: pd.DataFrame, te_ipr_desc: dict,
                        out_dir: Path) -> None:
    """Bubble plot: domain × clade, size = number of orthogroups."""
    confirmed = og_table[og_table["blast_confirmed"]]
    # Expand: one row per domain per clade
    rows = []
    for _, r in confirmed.iterrows():
        for dom in r["te_domains"].split(","):
            for clade in CLADE_ORDER:
                if r.get(clade, False):
                    rows.append({"domain": dom.strip(), "clade": clade})
    if not rows:
        print("No confirmed co-options for bubble plot.")
        return

    agg = pd.DataFrame(rows).groupby(["domain", "clade"]).size().reset_index(name="n")
    domains  = sorted(agg["domain"].unique())
    clades   = [c for c in CLADE_ORDER if c in agg["clade"].values]
    d_idx    = {d: i for i, d in enumerate(domains)}
    c_idx    = {c: i for i, c in enumerate(clades)}

    fig, ax = plt.subplots(figsize=(max(6, len(clades) * 0.7),
                                    max(4, len(domains) * 0.5)))
    for _, r in agg.iterrows():
        x = c_idx[r["clade"]]
        y = d_idx[r["domain"]]
        color = CLADE_COLORS.get(r["clade"], "#999")
        ax.scatter(x, y, s=r["n"] * 40, color=color, alpha=0.7, edgecolors="white")

    ax.set_xticks(range(len(clades)))
    ax.set_xticklabels(clades, fontsize=9, rotation=45)
    ax.set_yticks(range(len(domains)))
    ax.set_yticklabels([te_ipr_desc.get(d, d) for d in domains], fontsize=7)
    ax.set_xlabel("Clade", fontsize=10)
    ax.set_ylabel("TE domain", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "domain_bubble_plot.pdf")
    fig.savefig(out_dir / "domain_bubble_plot.png", dpi=200)
    plt.close(fig)
    print("Saved: domain_bubble_plot.pdf")


def plot_cooption_heatmap(og_table: pd.DataFrame, out_dir: Path) -> None:
    confirmed = og_table[og_table["blast_confirmed"]]
    mat = confirmed[[c for c in CLADE_ORDER if c in confirmed.columns]].astype(int)
    if mat.empty:
        return

    fig, ax = plt.subplots(figsize=(max(6, len(mat.columns) * 0.7),
                                    max(4, len(mat) * 0.25 + 2)))
    im = ax.imshow(mat.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, fontsize=8, rotation=45)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(confirmed["orthogroup"].values, fontsize=5)
    fig.colorbar(im, ax=ax, shrink=0.4, label="Clade presence")
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "cooption_heatmap.pdf")
    plt.close(fig)
    print("Saved: cooption_heatmap.pdf")


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))
    species_list = list(clade_map.keys())

    te_ipr = (set(args.te_domains.split(",")) if args.te_domains
              else set(TE_DOMAINS_DEFAULT.keys()))
    te_ipr_desc = TE_DOMAINS_DEFAULT if not args.te_domains else {k: k for k in te_ipr}

    print("Loading IPR domain annotations ...")
    ipr_genes = load_ipr_domain_presence(Path(args.ipr_dir), species_list, te_ipr)
    print(f"  {len(ipr_genes)} gene-domain rows with TE signatures")

    print("Building DIAMOND database ...")
    db_path = str(out_dir / "te_lib_diamond")
    build_diamond_db(args.te_lib, db_path)

    print("Screening orthogroups ...")
    og_table = screen_orthogroups(
        Path(args.orthogroup_seqs), ipr_genes, db_path,
        args.blast_pid, args.blast_cov, args.threads, species_list)
    n_conf = og_table["blast_confirmed"].sum()
    print(f"  {len(og_table)} OGs with TE domains; {n_conf} BLAST-confirmed")

    og_table = add_clade_presence(og_table, Path(args.orthogroup_seqs), clade_map)
    og_table.to_csv(out_dir / "cooption_orthogroups.tsv", sep="\t", index=False)

    # Domain summary
    dom_rows = []
    for _, r in og_table[og_table["blast_confirmed"]].iterrows():
        for d in r["te_domains"].split(","):
            dom_rows.append({"domain": d.strip(),
                             "desc": te_ipr_desc.get(d.strip(), d.strip()),
                             "orthogroup": r["orthogroup"]})
    pd.DataFrame(dom_rows).to_csv(out_dir / "cooption_domain_summary.tsv",
                                   sep="\t", index=False)

    plot_domain_bubble(og_table, te_ipr_desc, fig_dir)
    plot_cooption_heatmap(og_table, fig_dir)
    print("Done.")


if __name__ == "__main__":
    main()
