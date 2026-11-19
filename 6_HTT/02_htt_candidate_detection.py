#!/usr/bin/env python3

# identifies candidate horizontal TE transfer events across clades

# imports
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLADE_COLORS = {
    "HO": "#e6194b", "NU": "#f58231", "AP": "#ffe119", "CR": "#bfef45",
    "MI": "#3cb44b", "BL": "#42d4f4", "OL": "#4363d8", "CH": "#911eb4",
    "NE": "#f032e6", "BS": "#a9a9a9", "EN": "#469990", "ZO": "#9a6324",
    "KX": "#800000", "GL": "#808000", "MO": "#000075", "MU": "#e6beff",
    "BA": "#aaffc3", "AS": "#ffd8b1",
}

# parse arguments
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lib_dir",      required=True)
    p.add_argument("--species_tree", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="htt_candidates_out")
    p.add_argument("--identity",     type=float, default=85.0)
    p.add_argument("--aln_len",      type=int,   default=100)
    p.add_argument("--anomaly_th",   type=float, default=0.10)
    p.add_argument("--n_threads",    type=int,   default=8)
    return p.parse_args()

# get pairwise species distances from phylogeny
def load_pairwise_distances(nwk_file: str) -> dict[tuple[str, str], float]:
  
    try:
        from ete3 import Tree
    except ImportError:
        print("ete3 not found; install: pip install ete3", file=sys.stderr)
        return {}

    t = Tree(nwk_file, format=1)
    leaves = t.get_leaves()
    dist = {}
    for i, a in enumerate(leaves):
        for b in leaves[i + 1:]:
            d = a.get_distance(b)
            dist[(a.name, b.name)] = d
            dist[(b.name, a.name)] = d
    return dist

# run blast between two TE libraries
def run_blast_pairwise(query_fa: Path, subject_fa: Path,
                        identity: float, aln_len: int,
                        n_threads: int) -> pd.DataFrame:
    """BLAST query against subject; return hits above thresholds."""
    cols = ["qseqid", "sseqid", "pident", "length", "qcovs", "evalue", "bitscore"]
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=True) as tmp:
        cmd = [
            "blastn", "-query", str(query_fa), "-subject", str(subject_fa),
            "-outfmt", "6 " + " ".join(cols),
            "-perc_identity", str(identity),
            "-word_size", "11",
            "-num_threads", str(n_threads),
            "-out", tmp.name,
        ]
        ret = subprocess.run(cmd, capture_output=True)
        if ret.returncode != 0 or Path(tmp.name).stat().st_size == 0:
            return pd.DataFrame(columns=cols)
        df = pd.read_csv(tmp.name, sep="\t", header=None, names=cols)
    df = df[df["length"] >= aln_len].copy()
    return df

# extract species code from TE header
def sp_from_header(header: str) -> str:
    """Extract species prefix from TE copy header (format: SPECIES__family)."""
    return header.split("__")[0] if "__" in header else header.split("_")[0]

# get clade from species
def clade_from_sp(sp: str, clade_map: dict) -> str:
    return clade_map.get(sp, "?")


# Figures

# plot anomaly score distribution
def plot_anomaly_dist(candidates: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(candidates["anomaly_score"], bins=60, color="#e31a1c", alpha=0.7)
    ax.set_xlabel("Anomaly score (expected − observed divergence)", fontsize=10)
    ax.set_ylabel("Candidate pairs", fontsize=10)
    ax.axvline(0.10, color="blue", linestyle="--", linewidth=1, label="Threshold 0.10")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "htt_anomaly_dist.pdf")
    plt.close(fig)
    print("Saved: htt_anomaly_dist.pdf")


# main


def main():
    # parse args
    args = parse_args()
    # define output dirs
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # load metadata
    clade_df  = pd.read_csv(args.clade_map, sep="\t", header=None,
                             names=["species", "clade"])
    clade_map = dict(zip(clade_df["species"], clade_df["clade"]))
    
    # load species distances
    pairwise_dist = load_pairwise_distances(args.species_tree)
    # load TE libraries
    lib_dir  = Path(args.lib_dir)
    lib_files = sorted(lib_dir.glob("*.fa"))

    if not lib_files:
        print(f"No .fa files found in {lib_dir}", file=sys.stderr)
        sys.exit(1)
    
    # compare TE libraries across clades
    all_hits = []
    for i, q_fa in enumerate(lib_files):
        q_clade = q_fa.stem
        for j, s_fa in enumerate(lib_files):
            if j <= i:
                continue
            s_clade = s_fa.stem
            if q_clade == s_clade:
                continue
            print(f"  BLAST {q_clade} vs {s_clade} ...", end=" ")
            hits = run_blast_pairwise(q_fa, s_fa, args.identity, args.aln_len,
                                      args.n_threads)
            if hits.empty:
                print("0 hits")
                continue
            hits["q_clade"] = q_clade
            hits["s_clade"] = s_clade
            hits["q_sp"]    = hits["qseqid"].apply(sp_from_header)
            hits["s_sp"]    = hits["sseqid"].apply(sp_from_header)
            # calculate anomaly score
            hits["expected_div"] = hits.apply(
                lambda r: pairwise_dist.get((r["q_sp"], r["s_sp"]), np.nan), axis=1)
            hits["observed_div"] = 1.0 - hits["pident"] / 100.0
            hits["anomaly_score"] = hits["expected_div"] - hits["observed_div"]
            all_hits.append(hits)
            print(f"{len(hits)} hits")

    if not all_hits:
        print("No BLAST hits found.", file=sys.stderr)
        sys.exit(0)

    # identify candidate HTT events
    df = pd.concat(all_hits, ignore_index=True)
    candidates = df[df["anomaly_score"] >= args.anomaly_th].copy()
    print(f"\n{len(candidates)} candidate pairs with anomaly_score ≥ {args.anomaly_th}")
    candidates.to_csv(out_dir / "htt_candidates.tsv", sep="\t", index=False)

    # summarize candidate events
    summary = candidates.groupby(["q_clade", "s_clade"]).agg(
        n_pairs=("anomaly_score", "count"),
        max_identity=("pident", "max"),
        mean_anomaly=("anomaly_score", "mean"),
    ).reset_index()
    summary.to_csv(out_dir / "htt_candidate_summary.tsv", sep="\t", index=False)

    # generate figures
    plot_anomaly_dist(candidates, fig_dir)
    print(f"Summary: {len(summary)} clade-pair groups with HTT candidates")
    print("Done.")


if __name__ == "__main__":
    main()
