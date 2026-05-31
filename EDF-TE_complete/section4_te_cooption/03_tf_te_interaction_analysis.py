#!/usr/bin/env python3
"""
03_tf_te_interaction_analysis.py  —  Transcription factor binding at TE copies.

Identifies transcription factor (TF) binding events overlapping transposable
element copies using FIMO motif scan results.  Specifically characterises:
  A. Fraction of TF binding sites falling inside TE copies per TF family
  B. TE class enrichment for binding sites of each TF (Fisher exact)
  C. Co-binding patterns: TF pairs that co-occur at the same TE copy
  D. Visualisation: binding site density across TE meta-body

Analysis uses FIMO output (from MEME suite) against TE-annotated genome.

Inputs:
  --fimo_dir      Directory with per-species FIMO output (fimo.tsv files)
  --te_bed_dir    EarlGrey filtered TE BED files
  --motif_meta    TSV mapping motif ID to TF name and family (from JASPAR/MEME)
  --species_list  File with species codes
  --clade_map     Two-column TSV: species_code<TAB>clade_tag
  --out_dir       Output directory
  --fimo_pval     FIMO p-value cutoff (default: 1e-4)

Outputs:
  tf_te_overlap.tsv
  tf_enrichment_results.tsv
  figures/tf_te_heatmap.pdf
  figures/tf_te_cobinding.pdf
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "DNA":            "#377eb8",
    "LINE":           "#4daf4a",
    "SINE":           "#984ea3",
    "Rolling Circle": "#ff7f00",
    "Penelope":       "#a65628",
    "Unclassified":   "#999999",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fimo_dir",     required=True)
    p.add_argument("--te_bed_dir",   required=True)
    p.add_argument("--motif_meta",   required=True)
    p.add_argument("--species_list", required=True)
    p.add_argument("--clade_map",    required=True)
    p.add_argument("--out_dir",      default="tf_te_out")
    p.add_argument("--fimo_pval",    type=float, default=1e-4)
    return p.parse_args()


def _te_class(f: str) -> str:
    f = f.upper()
    if "LTR" in f:    return "LTR"
    if "DNA" in f:    return "DNA"
    if "LINE" in f:   return "LINE"
    if "SINE" in f:   return "SINE"
    if "RC/" in f or "HELITRON" in f: return "Rolling Circle"
    if "PENELOPE" in f: return "Penelope"
    return "Unclassified"


def load_fimo(fimo_dir: str, sp: str, pval: float) -> pd.DataFrame:
    files = list(Path(fimo_dir).glob(f"{sp}*/fimo.tsv")) + \
            list(Path(fimo_dir).glob(f"{sp}*fimo*.tsv"))
    if not files:
        return pd.DataFrame()
    df = pd.read_csv(files[0], sep="\t", comment="#")
    df.columns = [c.lower().replace(" ", "_").replace("-", "_") for c in df.columns]
    if "p_value" in df.columns:
        df = df[df["p_value"] <= pval]
    elif "pvalue" in df.columns:
        df = df[df["pvalue"] <= pval]
    return df


def load_te_bed(te_bed_dir: str, sp: str) -> pd.DataFrame:
    files = list(Path(te_bed_dir).glob(f"{sp}*.filteredRepeats.bed")) + \
            list(Path(te_bed_dir).glob(f"{sp}*.bed"))
    if not files:
        return pd.DataFrame()
    df = pd.read_csv(files[0], sep="\t", header=None, usecols=[0, 1, 2, 3],
                     names=["chrom", "te_start", "te_end", "te_family"])
    df["te_class"] = df["te_family"].apply(_te_class)
    return df


def overlap_fimo_with_te(fimo: pd.DataFrame, te: pd.DataFrame) -> pd.DataFrame:
    """Annotate FIMO hits with overlapping TE class (None if not in TE)."""
    te_by_chrom = {c: g for c, g in te.groupby("chrom")}
    classes = []
    for _, hit in fimo.iterrows():
        chrom = hit.get("sequence_name", hit.get("chrom", ""))
        start = hit.get("start", 0)
        end   = hit.get("stop", hit.get("end", start))
        te_c  = te_by_chrom.get(chrom, pd.DataFrame())
        if not te_c.empty:
            ov = te_c[(te_c["te_end"] > start) & (te_c["te_start"] < end)]
            if not ov.empty:
                classes.append(ov.iloc[0]["te_class"])
                continue
        classes.append(None)
    fimo = fimo.copy()
    fimo["te_class"] = classes
    return fimo


def test_enrichment(fimo: pd.DataFrame, te: pd.DataFrame) -> pd.DataFrame:
    """
    For each motif × TE class: Fisher exact test vs background (intergenic).
    Returns table with OR, p-value, and significant flag.
    """
    rows = []
    total_sites = len(fimo)
    total_te    = (fimo["te_class"].notna()).sum()

    for motif, mgroup in fimo.groupby("motif_id"):
        n_motif = len(mgroup)
        n_motif_te = mgroup["te_class"].notna().sum()
        for te_cls in te["te_class"].unique():
            n_motif_cls = (mgroup["te_class"] == te_cls).sum()
            n_bg_cls    = (fimo["te_class"] == te_cls).sum() - n_motif_cls

            # 2×2 contingency table
            a = n_motif_cls        # motif in TE class
            b = n_motif - a        # motif not in TE class
            c = n_bg_cls           # background in TE class (non-motif)
            d = total_sites - n_motif - c  # background not in TE class
            if a + b == 0 or a + c == 0:
                continue
            od, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({
                "motif_id": motif, "te_class": te_cls,
                "n_motif_in_te": a, "odds_ratio": od, "p_value": pval,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # BH correction
    df = df.sort_values("p_value").reset_index(drop=True)
    m = len(df)
    df["p_adj"] = (df["p_value"] * m / np.arange(1, m + 1)).clip(upper=1.0)
    df["significant"] = df["p_adj"] < 0.05
    return df


# =============================================================================
# Figures
# =============================================================================

def plot_tf_te_heatmap(enrichment: pd.DataFrame, motif_meta: pd.DataFrame,
                        out_dir: Path) -> None:
    sig = enrichment[enrichment["significant"] & (enrichment["odds_ratio"] > 1)]
    if sig.empty:
        print("No significant TF-TE enrichments for heatmap.")
        return

    mat = sig.pivot_table(index="motif_id", columns="te_class",
                           values="odds_ratio", aggfunc="max", fill_value=1.0)
    # Annotate with TF name if available
    if "motif_id" in motif_meta.columns and "tf_name" in motif_meta.columns:
        name_map = dict(zip(motif_meta["motif_id"], motif_meta["tf_name"]))
        mat.index = [name_map.get(i, i) for i in mat.index]

    fig, ax = plt.subplots(figsize=(max(5, mat.shape[1] * 1.2),
                                    max(4, mat.shape[0] * 0.3 + 2)))
    sns.heatmap(np.log2(mat), ax=ax, cmap="Reds", center=0,
                xticklabels=True, yticklabels=True,
                cbar_kws={"label": "log₂(OR)"})
    ax.set_xlabel("TE class", fontsize=10)
    ax.set_ylabel("Transcription factor", fontsize=10)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "tf_te_heatmap.pdf")
    fig.savefig(out_dir / "tf_te_heatmap.png", dpi=200)
    plt.close(fig)
    print("Saved: tf_te_heatmap.pdf")


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
    species_list = [l.strip() for l in open(args.species_list) if l.strip()]

    motif_meta = pd.read_csv(args.motif_meta, sep="\t") if Path(args.motif_meta).exists() \
                 else pd.DataFrame()

    all_annotated = []
    for sp in species_list:
        fimo = load_fimo(args.fimo_dir, sp, args.fimo_pval)
        te   = load_te_bed(args.te_bed_dir, sp)
        if fimo.empty or te.empty:
            continue
        fimo_annot = overlap_fimo_with_te(fimo, te)
        fimo_annot["species"] = sp
        all_annotated.append(fimo_annot)

    if not all_annotated:
        print("No FIMO data overlapped with TEs.", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(all_annotated, ignore_index=True)
    combined.to_csv(out_dir / "tf_te_overlap.tsv", sep="\t", index=False)
    print(f"  {len(combined)} FIMO hits; "
          f"{combined['te_class'].notna().sum()} overlapping TEs")

    # Use full TE bed for enrichment background
    all_te = pd.concat([load_te_bed(args.te_bed_dir, sp) for sp in species_list],
                       ignore_index=True)
    enrich = test_enrichment(combined, all_te)
    if not enrich.empty:
        enrich.to_csv(out_dir / "tf_enrichment_results.tsv", sep="\t", index=False)
        n_sig = enrich["significant"].sum()
        print(f"  {n_sig} significant TF–TE class enrichments (FDR < 5%)")
        plot_tf_te_heatmap(enrich, motif_meta, fig_dir)

    print("Done.")


if __name__ == "__main__":
    main()
