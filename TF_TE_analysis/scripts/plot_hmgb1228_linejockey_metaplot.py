#!/usr/bin/env python3
"""
Metaplot of HMGb_1228 coverage on LINE/I-Jockey copies (MU_Rhimi59),
with motif positions (HMGb3 PWM scan) per TE copy.

Outputs:
  /legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/figures/hmgb1228_LINE-I-Jockey_metaplot.pdf
  /legserv/Temp/Carlos/EDF-TE/results/TF_TE_analysis/figures/hmgb1228_LINE-I-Jockey_metaplot.svg
"""

from pathlib import Path
import re
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
REL_SCORE_THRESH = 0.82


def parse_meme_pwm(meme_file: Path, motif_name: str = "HMGb3"):
    lines = meme_file.read_text().splitlines()
    bg = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    # Parse background frequencies if present.
    for i, ln in enumerate(lines):
        if ln.startswith("Background letter frequencies"):
            if i + 1 < len(lines):
                parts = lines[i + 1].strip().split()
                if len(parts) >= 8:
                    bg = {parts[0]: float(parts[1]), parts[2]: float(parts[3]),
                          parts[4]: float(parts[5]), parts[6]: float(parts[7])}
            break

    # Parse motif matrix.
    start = None
    w = None
    for i, ln in enumerate(lines):
        if ln.strip() == f"MOTIF {motif_name}":
            for j in range(i + 1, min(i + 8, len(lines))):
                if "letter-probability matrix" in lines[j]:
                    start = j + 1
                    m = re.search(r"w=\s*(\d+)", lines[j])
                    if m:
                        w = int(m.group(1))
                    break
            break

    if start is None or w is None:
        raise RuntimeError(f"Could not parse motif {motif_name} from {meme_file}")

    pwm = []
    for k in range(start, start + w):
        vals = lines[k].strip().split()
        if len(vals) < 4:
            raise RuntimeError("Malformed PWM row in MEME file")
        pwm.append([float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3])])

    pwm = np.array(pwm, dtype=float)  # shape: (w, 4) A,C,G,T
    return pwm, bg


def reverse_complement(seq: str) -> str:
    tab = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tab)[::-1]


def score_windows(seq: str, pwm: np.ndarray, bg: dict):
    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    bg_vec = np.array([bg.get("A", 0.25), bg.get("C", 0.25), bg.get("G", 0.25), bg.get("T", 0.25)], dtype=float)
    bg_vec = np.clip(bg_vec, 1e-6, 1.0)

    w = pwm.shape[0]
    seq = seq.upper()

    # Log-likelihood ratio per position.
    llr = np.log2(np.clip(pwm, 1e-6, 1.0) / bg_vec[None, :])
    max_score = llr.max(axis=1).sum()
    min_score = llr.min(axis=1).sum()

    hits = []
    for i in range(0, len(seq) - w + 1):
        kmer = seq[i:i + w]
        if any(ch not in base_to_idx for ch in kmer):
            continue

        idx = np.array([base_to_idx[ch] for ch in kmer], dtype=int)
        s_f = llr[np.arange(w), idx].sum()

        rc = reverse_complement(kmer)
        idx_rc = np.array([base_to_idx[ch] for ch in rc], dtype=int)
        s_r = llr[np.arange(w), idx_rc].sum()

        if s_f >= s_r:
            s = s_f
            strand = "+"
        else:
            s = s_r
            strand = "-"

        rel = (s - min_score) / (max_score - min_score + 1e-9)
        if rel >= REL_SCORE_THRESH:
            hits.append((i, i + w, strand, s, rel))

    # Merge overlapping hits, keep best per overlap cluster.
    if not hits:
        return []
    hits.sort(key=lambda x: (x[0], -x[4]))
    merged = []
    cur = list(hits[0])
    for h in hits[1:]:
        if h[0] < cur[1]:
            if h[4] > cur[4]:
                cur = list(h)
            else:
                cur[1] = max(cur[1], h[1])
        else:
            merged.append(tuple(cur))
            cur = list(h)
    merged.append(tuple(cur))
    return merged


def load_line_jockey_copies(te_bed: Path, min_len: int):
    rows = []
    with te_bed.open() as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            chrom, start, end, te_class, score, strand = p[:6]
            if te_class != "LINE/I-Jockey":
                continue
            start = int(start)
            end = int(end)
            length = end - start
            if length >= min_len:
                rows.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": length,
                    "label": f"{chrom}:{start}-{end} ({strand})",
                })
    df = pd.DataFrame(rows).sort_values("length", ascending=False).reset_index(drop=True)
    return df


def main():
    te_df = load_line_jockey_copies(TE_BED, MIN_LEN)
    if te_df.empty:
        raise RuntimeError("No LINE/I-Jockey copies >=200bp found")

    pwm, bg = parse_meme_pwm(MOTIF_MEME, motif_name="HMGb3")

    fasta = SeqIO.to_dict(SeqIO.parse(str(GENOME_FASTA), "fasta"))
    bw = pyBigWig.open(str(BW_PATH))

    profiles = []
    per_te = []
    hit_rows = []

    for _, r in te_df.iterrows():
        chrom, start, end, strand = r["chrom"], int(r["start"]), int(r["end"]), r["strand"]
        seq = str(fasta[chrom].seq[start:end]).upper()
        vals = np.array(bw.values(chrom, start, end, numpy=True), dtype=float)
        vals = np.nan_to_num(vals, nan=0.0)

        # Orient to TE strand.
        if strand == "-":
            seq = reverse_complement(seq)
            vals = vals[::-1]

        # Motif hits in oriented TE coordinates.
        hits = score_windows(seq, pwm, bg)

        # Save motif hit coordinates in both TE-oriented and genome coordinates.
        for s, e, mstrand, score, rel in hits:
            if strand == "+":
                g_start = start + s
                g_end = start + e
                motif_genome_strand = mstrand
            else:
                # Convert from TE-oriented (reverse-complemented) back to genome coords.
                g_start = end - e
                g_end = end - s
                motif_genome_strand = "+" if mstrand == "-" else "-"

            hit_rows.append({
                "te_label": r["label"],
                "te_chrom": chrom,
                "te_start": start,
                "te_end": end,
                "te_strand": strand,
                "te_length": int(r["length"]),
                "motif_start_te_oriented": int(s),
                "motif_end_te_oriented": int(e),
                "motif_start_genome": int(g_start),
                "motif_end_genome": int(g_end),
                "motif_strand_genome": motif_genome_strand,
                "motif_score_llr": float(score),
                "motif_rel_score": float(rel),
                "motif_start_pct": float(100.0 * s / r["length"]),
                "motif_end_pct": float(100.0 * e / r["length"]),
            })

        # Resample to a common metaplot x-axis.
        x_old = np.linspace(0, 1, num=len(vals), endpoint=False)
        x_new = np.linspace(0, 1, num=N_BINS, endpoint=False)
        y_new = np.interp(x_new, x_old, vals)
        profiles.append(y_new)

        per_te.append({
            "label": r["label"],
            "length": int(r["length"]),
            "signal": vals,
            "hits": hits,
        })

    bw.close()

    prof = np.vstack(profiles)
    mean_prof = prof.mean(axis=0)
    sd_prof = prof.std(axis=0)

    # ---- Plot ----
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 1, height_ratios=[1.1, 1.6], hspace=0.25)

    # Top: metaplot mean +/- sd
    ax1 = fig.add_subplot(gs[0, 0])
    x = np.linspace(0, 100, N_BINS, endpoint=False)
    ax1.plot(x, mean_prof, color="#2C7FB8", lw=2.2, label="Mean HMGb_1228 coverage")
    ax1.fill_between(x, mean_prof - sd_prof, mean_prof + sd_prof,
                     color="#2C7FB8", alpha=0.22, lw=0, label="±1 SD")
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("Relative position on LINE/I-Jockey copy (%)")
    ax1.set_ylabel("BigWig signal")
    ax1.set_title("HMGb_1228 coverage metaplot on LINE/I-Jockey (4 copies ≥200 bp)", loc="left", fontsize=11)
    ax1.legend(frameon=False, fontsize=9)
    ax1.spines[["top", "right"]].set_visible(False)

    # Bottom: per-copy tracks + motif hits
    ax2 = fig.add_subplot(gs[1, 0])
    y_gap = 1.4
    y_base = np.arange(len(per_te))[::-1] * y_gap

    for i, rec in enumerate(per_te):
        y0 = y_base[i]
        sig = rec["signal"]
        L = rec["length"]
        xs = np.linspace(0, 100, num=L, endpoint=False)

        # Signal track (normalized per TE for readability).
        if sig.max() > 0:
            sig_plot = sig / sig.max() * 0.8
        else:
            sig_plot = sig
        ax2.plot(xs, y0 + sig_plot, color="#4A4A4A", lw=1.0)
        ax2.hlines(y0, 0, 100, color="#BEBEBE", lw=0.7)

        # Motif hits as red blocks in relative coordinates.
        for s, e, strand, score, rel in rec["hits"]:
            x1 = 100.0 * s / L
            x2 = 100.0 * e / L
            ax2.add_patch(plt.Rectangle((x1, y0 - 0.18), x2 - x1, 0.22,
                                        color="#D7301F", alpha=0.9))

        ax2.text(101.2, y0 + 0.15, f"{rec['label']}  L={L}bp  motifs={len(rec['hits'])}",
                 fontsize=8, va="center", ha="left")

    ax2.set_xlim(0, 130)
    ax2.set_ylim(-0.8, y_base[0] + 1.2)
    ax2.set_yticks([])
    ax2.set_xlabel("Relative position on each TE copy (%)")
    ax2.set_title("Per-copy coverage and HMGb3 motif locations", loc="left", fontsize=10)
    ax2.spines[["top", "right", "left"]].set_visible(False)

    # Mini legend for motif blocks.
    ax2.add_patch(plt.Rectangle((0.5, -0.63), 2.8, 0.2, color="#D7301F", alpha=0.9))
    ax2.text(3.8, -0.52, f"HMGb3 motif hit (PWM relative score >= {REL_SCORE_THRESH:.2f})",
             fontsize=8, va="center")

    fig.suptitle("MU_Rhimi59 LINE/I-Jockey copies (EarlGrey) and HMGb_1228 occupancy", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save motif coordinates table.
    if hit_rows:
        hit_df = pd.DataFrame(hit_rows).sort_values([
            "te_chrom", "te_start", "motif_start_genome"
        ])
    else:
        hit_df = pd.DataFrame(columns=[
            "te_label", "te_chrom", "te_start", "te_end", "te_strand", "te_length",
            "motif_start_te_oriented", "motif_end_te_oriented",
            "motif_start_genome", "motif_end_genome", "motif_strand_genome",
            "motif_score_llr", "motif_rel_score", "motif_start_pct", "motif_end_pct"
        ])
    hit_tsv = OUTDIR / "hmgb1228_LINE-I-Jockey_motif_hits.tsv"
    hit_df.to_csv(hit_tsv, sep="\t", index=False)
    print(f"Saved: {hit_tsv}")

    for ext in ("pdf", "svg"):
        out = OUTDIR / f"hmgb1228_LINE-I-Jockey_metaplot.{ext}"
        fig.savefig(out, dpi=250, bbox_inches="tight")
        print(f"Saved: {out}")


if __name__ == "__main__":
    main()
