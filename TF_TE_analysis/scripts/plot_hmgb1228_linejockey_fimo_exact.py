#!/usr/bin/env python3
"""
Exact FIMO-based motif mapping + HMGb_1228 metaplot on LINE/I-Jockey copies.

Produces:
1) Custom metaplot with per-copy tracks and exact FIMO motif positions.
2) TSV with exact motif coordinates in genome and TE-oriented coordinates.
3) BED + FASTA inputs reused by deepTools commands.
"""

from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from Bio import SeqIO
import pyBigWig

BW_PATH = Path("/legserv/Temp/Carlos/EDF-TE/data/TFsDAP/bigwigs/") / \
    "Rhimi59_2_1891228_R_microsporus_ATCC11559_DAPi709D03-ORGi5043-96P422_dedup_norm.bw"
GENOME_FASTA = Path("/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Rhimi59.fasta")
TE_BED = Path("/data/carlos/EDF-TE/TE_analysis/annotated/MU_Rhimi59_EarlGrey/") / \
    "MU_Rhimi59_mergedRepeats/looseMerge/MU_Rhimi59.filteredRepeats.bed"
MOTIF_MEME = Path("/legserv/Temp/Carlos/EDF-TE/data/TFsDAP/motifs/HMGb3.meme")

OUTDIR = Path("/legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

N_BINS = 200
MIN_LEN = 200
FIMO_THRESH = "1e-4"
TOP_Q_LABELS_PER_TE = 2


def reverse_complement(seq: str) -> str:
    tab = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(tab)[::-1]


def load_line_jockey_copies(te_bed: Path, min_len: int) -> pd.DataFrame:
    rows = []
    with te_bed.open() as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            chrom, start, end, te_class, score, strand = p[:6]
            if te_class != "LINE/I-Jockey":
                continue
            start, end = int(start), int(end)
            length = end - start
            if length >= min_len:
                rows.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": length,
                    "name": f"{chrom}:{start}-{end}|{strand}",
                    "label": f"{chrom}:{start}-{end} ({strand})",
                })
    return pd.DataFrame(rows).sort_values("length", ascending=False).reset_index(drop=True)


def write_bed_and_fasta(te_df: pd.DataFrame, genome_fa: Path, out_bed: Path, out_fa: Path):
    fasta = SeqIO.to_dict(SeqIO.parse(str(genome_fa), "fasta"))

    with out_bed.open("w") as b, out_fa.open("w") as f:
        for _, r in te_df.iterrows():
            b.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['name']}\t0\t{r['strand']}\n")

            seq = str(fasta[r["chrom"]].seq[r["start"]:r["end"]]).upper()
            if r["strand"] == "-":
                seq = reverse_complement(seq)
            f.write(f">{r['name']}\n{seq}\n")


def run_fimo_exact(motif_meme: Path, fasta_in: Path, outdir: Path, thresh: str):
    if outdir.exists():
        # Keep previous run cleanly reproducible.
        subprocess.run(["rm", "-rf", str(outdir)], check=True)

    cmd = [
        "fimo",
        "--oc", str(outdir),
        "--thresh", str(thresh),
        str(motif_meme),
        str(fasta_in),
    ]
    subprocess.run(cmd, check=True)


def parse_fimo_hits(fimo_tsv: Path, te_df: pd.DataFrame) -> pd.DataFrame:
    # FIMO tsv has comment lines starting with '#'.
    fimo = pd.read_csv(fimo_tsv, sep="\t", comment="#")
    if fimo.empty:
        return pd.DataFrame(columns=[
            "te_name", "te_label", "te_chrom", "te_start", "te_end", "te_strand", "te_length",
            "motif_start_te_oriented", "motif_end_te_oriented",
            "motif_start_genome", "motif_end_genome", "motif_strand_genome",
            "p_value", "q_value", "score", "matched_sequence",
            "motif_start_pct", "motif_end_pct",
        ])

    # Ensure numeric columns are numeric.
    for c in ["start", "stop", "score", "p-value", "q-value"]:
        if c in fimo.columns:
            fimo[c] = pd.to_numeric(fimo[c], errors="coerce")

    te_map = {r["name"]: r for _, r in te_df.iterrows()}

    rows = []
    for _, h in fimo.iterrows():
        name = str(h["sequence_name"])

        # Two possible formats:
        # 1) sequence_name keeps FASTA header name used by this script.
        # 2) FIMO decodes genomic-style IDs and returns chrom + genomic start/stop.
        if name in te_map:
            r = te_map[name]
            s0 = int(h["start"]) - 1
            e0 = int(h["stop"])
        else:
            chrom = name
            g_start = int(h["start"]) - 1
            g_end = int(h["stop"])
            sub = te_df[(te_df["chrom"] == chrom) &
                        (te_df["start"] <= g_start) &
                        (te_df["end"] >= g_end)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            name = r["name"]

            if r["strand"] == "+":
                s0 = g_start - int(r["start"])
                e0 = g_end - int(r["start"])
            else:
                s0 = int(r["end"]) - g_end
                e0 = int(r["end"]) - g_start

        if r["strand"] == "+":
            g_start = r["start"] + s0
            g_end = r["start"] + e0
            mstrand = h["strand"]
        else:
            # Sequence was reverse-complemented to TE orientation.
            g_start = r["end"] - e0
            g_end = r["end"] - s0
            mstrand = "+" if h["strand"] == "-" else "-"

        rows.append({
            "te_name": name,
            "te_label": r["label"],
            "te_chrom": r["chrom"],
            "te_start": int(r["start"]),
            "te_end": int(r["end"]),
            "te_strand": r["strand"],
            "te_length": int(r["length"]),
            "motif_start_te_oriented": int(s0),
            "motif_end_te_oriented": int(e0),
            "motif_start_genome": int(g_start),
            "motif_end_genome": int(g_end),
            "motif_strand_genome": mstrand,
            "p_value": float(h.get("p-value", np.nan)),
            "q_value": float(h.get("q-value", np.nan)),
            "score": float(h.get("score", np.nan)),
            "matched_sequence": str(h.get("matched_sequence", "")),
            "motif_start_pct": 100.0 * s0 / r["length"],
            "motif_end_pct": 100.0 * e0 / r["length"],
        })

    out = pd.DataFrame(rows).sort_values(["te_chrom", "te_start", "motif_start_genome"]).reset_index(drop=True)
    return out


def collapse_overlapping_hits(hit_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse overlapping motif hits per TE, keeping best q-value/score for each cluster."""
    if hit_df.empty:
        return hit_df.copy()

    clusters = []
    for te_name, sub in hit_df.groupby("te_name"):
        s = sub.sort_values(["motif_start_te_oriented", "motif_end_te_oriented"]).reset_index(drop=True)
        cur_start = int(s.loc[0, "motif_start_te_oriented"])
        cur_end = int(s.loc[0, "motif_end_te_oriented"])
        cur_rows = [s.loc[0]]

        for i in range(1, len(s)):
            r = s.loc[i]
            rs = int(r["motif_start_te_oriented"])
            re = int(r["motif_end_te_oriented"])
            if rs <= cur_end:
                cur_end = max(cur_end, re)
                cur_rows.append(r)
            else:
                cdf = pd.DataFrame(cur_rows)
                best = cdf.sort_values(["q_value", "p_value", "score"], ascending=[True, True, False]).iloc[0]
                clusters.append({
                    "te_name": te_name,
                    "te_label": best["te_label"],
                    "te_chrom": best["te_chrom"],
                    "te_start": int(best["te_start"]),
                    "te_end": int(best["te_end"]),
                    "te_strand": best["te_strand"],
                    "te_length": int(best["te_length"]),
                    "cluster_start_te_oriented": int(cur_start),
                    "cluster_end_te_oriented": int(cur_end),
                    "cluster_start_pct": 100.0 * cur_start / int(best["te_length"]),
                    "cluster_end_pct": 100.0 * cur_end / int(best["te_length"]),
                    "best_motif_start_genome": int(best["motif_start_genome"]),
                    "best_motif_end_genome": int(best["motif_end_genome"]),
                    "best_motif_strand_genome": best["motif_strand_genome"],
                    "best_p_value": float(best["p_value"]),
                    "best_q_value": float(best["q_value"]),
                    "best_score": float(best["score"]),
                    "best_matched_sequence": best["matched_sequence"],
                    "n_raw_hits_in_cluster": int(len(cdf)),
                })
                cur_start, cur_end = rs, re
                cur_rows = [r]

        cdf = pd.DataFrame(cur_rows)
        best = cdf.sort_values(["q_value", "p_value", "score"], ascending=[True, True, False]).iloc[0]
        clusters.append({
            "te_name": te_name,
            "te_label": best["te_label"],
            "te_chrom": best["te_chrom"],
            "te_start": int(best["te_start"]),
            "te_end": int(best["te_end"]),
            "te_strand": best["te_strand"],
            "te_length": int(best["te_length"]),
            "cluster_start_te_oriented": int(cur_start),
            "cluster_end_te_oriented": int(cur_end),
            "cluster_start_pct": 100.0 * cur_start / int(best["te_length"]),
            "cluster_end_pct": 100.0 * cur_end / int(best["te_length"]),
            "best_motif_start_genome": int(best["motif_start_genome"]),
            "best_motif_end_genome": int(best["motif_end_genome"]),
            "best_motif_strand_genome": best["motif_strand_genome"],
            "best_p_value": float(best["p_value"]),
            "best_q_value": float(best["q_value"]),
            "best_score": float(best["score"]),
            "best_matched_sequence": best["matched_sequence"],
            "n_raw_hits_in_cluster": int(len(cdf)),
        })

    return pd.DataFrame(clusters).sort_values(["te_chrom", "te_start", "cluster_start_te_oriented"]).reset_index(drop=True)


def build_signal_profiles(te_df: pd.DataFrame, bw_path: Path):
    bw = pyBigWig.open(str(bw_path))
    profiles = []
    per_te = []

    for _, r in te_df.iterrows():
        chrom, start, end, strand = r["chrom"], int(r["start"]), int(r["end"]), r["strand"]
        vals = np.array(bw.values(chrom, start, end, numpy=True), dtype=float)
        vals = np.nan_to_num(vals, nan=0.0)
        if strand == "-":
            vals = vals[::-1]

        x_old = np.linspace(0, 1, num=len(vals), endpoint=False)
        x_new = np.linspace(0, 1, num=N_BINS, endpoint=False)
        y_new = np.interp(x_new, x_old, vals)
        profiles.append(y_new)

        per_te.append({
            "name": r["name"],
            "label": r["label"],
            "length": int(r["length"]),
            "signal": vals,
        })

    bw.close()
    return np.vstack(profiles), per_te


def plot_custom_metaplot(te_df: pd.DataFrame, per_te: list, prof: np.ndarray, hit_df: pd.DataFrame,
                         cluster_df: pd.DataFrame, out_prefix: Path):
    mean_prof = prof.mean(axis=0)
    sd_prof = prof.std(axis=0)

    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 1, height_ratios=[1.1, 1.6], hspace=0.25)

    ax1 = fig.add_subplot(gs[0, 0])
    x = np.linspace(0, 100, N_BINS, endpoint=False)
    ax1.plot(x, mean_prof, color="#2C7FB8", lw=2.2, label="Mean HMGb_1228 coverage")
    ax1.fill_between(x, mean_prof - sd_prof, mean_prof + sd_prof,
                     color="#2C7FB8", alpha=0.22, lw=0, label="±1 SD")
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("Relative position on LINE/I-Jockey copy (%)")
    ax1.set_ylabel("BigWig signal")
    ax1.set_title("HMGb_1228 coverage metaplot on LINE/I-Jockey (4 copies >=200 bp)", loc="left", fontsize=11)
    ax1.legend(frameon=False, fontsize=9)
    ax1.spines[["top", "right"]].set_visible(False)

    ax2 = fig.add_subplot(gs[1, 0])
    y_gap = 1.4
    y_base = np.arange(len(per_te))[::-1] * y_gap

    cl_map = {}
    if not cluster_df.empty:
        for te_name, sub in cluster_df.groupby("te_name"):
            cl_map[te_name] = sub

    for i, rec in enumerate(per_te):
        y0 = y_base[i]
        sig = rec["signal"]
        L = rec["length"]
        xs = np.linspace(0, 100, num=L, endpoint=False)

        if sig.max() > 0:
            sig_plot = sig / sig.max() * 0.8
        else:
            sig_plot = sig
        ax2.plot(xs, y0 + sig_plot, color="#4A4A4A", lw=1.0)
        ax2.hlines(y0, 0, 100, color="#BEBEBE", lw=0.7)

        sub = cl_map.get(rec["name"], None)
        n_hits = 0
        if sub is not None:
            n_hits = len(sub)
            for _, h in sub.iterrows():
                x1 = 100.0 * h["cluster_start_te_oriented"] / L
                x2 = 100.0 * h["cluster_end_te_oriented"] / L
                ax2.add_patch(plt.Rectangle((x1, y0 - 0.18), x2 - x1, 0.22,
                                            color="#D7301F", alpha=0.9))

            # Annotate top hits by q-value for readability
            top = sub.sort_values("best_q_value").head(TOP_Q_LABELS_PER_TE)
            for _, h in top.iterrows():
                xm = 50.0 * (h["cluster_start_te_oriented"] + h["cluster_end_te_oriented"]) / L
                ax2.text(xm, y0 + 0.84, f"q={h['best_q_value']:.1e}",
                         fontsize=6.5, rotation=25, ha="center", va="bottom", color="#222222")

        ax2.text(101.2, y0 + 0.15, f"{rec['label']}  L={L}bp  clusters={n_hits}",
                 fontsize=8, va="center", ha="left")

    ax2.set_xlim(0, 130)
    ax2.set_ylim(-0.8, y_base[0] + 1.2)
    ax2.set_yticks([])
    ax2.set_xlabel("Relative position on each TE copy (%)")
    ax2.set_title("Per-copy coverage and exact FIMO HMGb3 motif locations", loc="left", fontsize=10)
    ax2.spines[["top", "right", "left"]].set_visible(False)

    ax2.add_patch(plt.Rectangle((0.5, -0.63), 2.8, 0.2, color="#D7301F", alpha=0.9))
    ax2.text(3.8, -0.52, "Exact FIMO HMGb3 cluster", fontsize=8, va="center")

    fig.suptitle("MU_Rhimi59 LINE/I-Jockey copies and HMGb_1228 occupancy", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    for ext in ("pdf", "svg"):
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=250, bbox_inches="tight")
        print(f"Saved: {out}")


def main():
    te_df = load_line_jockey_copies(TE_BED, MIN_LEN)
    if te_df.empty:
        raise RuntimeError("No LINE/I-Jockey copies >=200bp found")

    out_bed = OUTDIR / "LINE_I-Jockey_4copies_ge200.bed"
    out_fa = OUTDIR / "LINE_I-Jockey_4copies_ge200.oriented.fa"
    write_bed_and_fasta(te_df, GENOME_FASTA, out_bed, out_fa)
    print(f"Saved: {out_bed}")
    print(f"Saved: {out_fa}")

    fimo_out = OUTDIR / "fimo_HMGb3_LINE-I-Jockey_exact"
    run_fimo_exact(MOTIF_MEME, out_fa, fimo_out, FIMO_THRESH)

    fimo_tsv = fimo_out / "fimo.tsv"
    hit_df = parse_fimo_hits(fimo_tsv, te_df)
    hit_tsv = OUTDIR / "hmgb1228_LINE-I-Jockey_motif_hits_exact_fimo.tsv"
    hit_df.to_csv(hit_tsv, sep="\t", index=False)
    print(f"Saved: {hit_tsv}")

    cluster_df = collapse_overlapping_hits(hit_df)
    cluster_tsv = OUTDIR / "hmgb1228_LINE-I-Jockey_motif_clusters_exact_fimo.tsv"
    cluster_df.to_csv(cluster_tsv, sep="\t", index=False)
    print(f"Saved: {cluster_tsv}")

    prof, per_te = build_signal_profiles(te_df, BW_PATH)
    out_prefix = OUTDIR / "hmgb1228_LINE-I-Jockey_metaplot_exactFIMO"
    plot_custom_metaplot(te_df, per_te, prof, hit_df, cluster_df, out_prefix)


if __name__ == "__main__":
    main()
