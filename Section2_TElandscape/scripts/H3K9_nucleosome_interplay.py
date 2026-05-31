#!/usr/bin/env python3
"""
H3K9_nucleosome_interplay.py
================================================================
Integrated analysis for Rhizopus microsporus (MU_Rhimi59) joining:
  - DANPOS nucleosome occupancy WIG
  - H3K9me3 ChIP-seq broad peaks (MACS2)
  - EarlGrey TE annotation (filteredRepeats.bed)
  - rDNA / NOR (nucleolar) regions detected with barrnap

Four questions answered in one run:

  (A) Is nucleosome occupancy higher inside H3K9me3 broad peaks
      than in the rest of the genome?

  (B) Is nucleosome ARRANGEMENT (NRL, fuzziness, phasing) different
      in H3K9me3 regions vs the rest?

  (C) How does the TE-vs-NOR(occupied) bias change when stratifying
      by H3K9me3? (i.e. is the "TEs prefer NOR" signal really a
      "TEs prefer H3K9me3" signal?)

  (D) Where do TEs accumulate when we classify the genome into:
        - Phased-occupied  (high occupancy, well positioned, fixed)
        - Fuzzy-occupied   (high occupancy, disorganised)
        - Phased-depleted  (low occupancy, regular spacing — gene-like)
        - Open-fuzzy       (low occupancy, disorganised — NDRs)

Outputs written to:
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes_H3K9/

Required: numpy, scipy, pandas, matplotlib, pybedtools (optional),
barrnap on PATH (optional, for rDNA masking).
================================================================
"""

import argparse, os, re, subprocess, sys, math
from pathlib import Path
from collections import defaultdict

import numpy  as np
import pandas as pd
from scipy.signal import find_peaks, correlate
from scipy.stats  import mannwhitneyu

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── default I/O ──────────────────────────────────────────────────────────────
DEF_WIG   = "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"
DEF_PEAKS = "/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/Epigenetics/RmH3K9me3/H3K9_WT_MACS2_broad_peaks.bed"
DEF_TEBED = "/data/carlos/EDF-TE/TE_analysis/annotated/MU_Rhimi59_EarlGrey/MU_Rhimi59_summaryFiles/MU_Rhimi59.filteredRepeats.bed"
DEF_FAI   = "/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Rhimi59.fasta.fai"
DEF_FA    = "/legserv/Temp/Carlos/EDF-TE/data/genomes/MU_Rhimi59.fasta"
DEF_OUT   = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes_H3K9"

# ─────────────────────────────────────────────────────────────────────────────
def log(msg):
    print(f"[H3K9xNUC] {msg}", flush=True)


# ── 1. Genome / WIG loading ──────────────────────────────────────────────────
def load_chrom_sizes(fai):
    out = {}
    with open(fai) as fh:
        for ln in fh:
            f = ln.split("\t")
            out[f[0]] = int(f[1])
    return out


def load_wig_arrays(wig_path, chrom_sizes):
    """Return dict {chrom: np.float32 array of length size}, value 0 outside coverage."""
    arrs = {c: np.zeros(s, dtype=np.float32) for c, s in chrom_sizes.items()}
    chrom = None; mode = None; pos = 1; step = 1; span = 1
    n_pos = 0
    with open(wig_path) as fh:
        for ln in fh:
            ln = ln.rstrip()
            if not ln or ln[0] == "#" or ln.startswith("track"):
                continue
            if ln.startswith("fixedStep"):
                mode = "fixed"
                m = re.search(r"chrom=(\S+)", ln);  chrom = m.group(1)
                m = re.search(r"start=(\d+)", ln);  pos   = int(m.group(1)) if m else 1
                m = re.search(r"step=(\d+)",  ln);  step  = int(m.group(1)) if m else 1
                m = re.search(r"span=(\d+)",  ln);  span  = int(m.group(1)) if m else 1
                continue
            if ln.startswith("variableStep"):
                mode = "variable"
                m = re.search(r"chrom=(\S+)", ln);  chrom = m.group(1)
                m = re.search(r"span=(\d+)",  ln);  span  = int(m.group(1)) if m else 1
                continue
            if chrom not in arrs:
                continue
            if mode == "fixed":
                v = float(ln)
                if v != 0.0:
                    end = min(pos - 1 + span, arrs[chrom].size)
                    arrs[chrom][pos - 1:end] = v
                    n_pos += end - (pos - 1)
                pos += step
            elif mode == "variable":
                p, v = ln.split()
                p = int(p); v = float(v)
                if v != 0.0:
                    end = min(p - 1 + span, arrs[chrom].size)
                    arrs[chrom][p - 1:end] = v
                    n_pos += end - (p - 1)
    log(f"WIG loaded: non-zero positions ≈ {n_pos:,}")
    return arrs


# ── 2. BED helpers ───────────────────────────────────────────────────────────
def load_bed_intervals(path):
    """Return dict {chrom: list[(start, end, *rest)]} (0-based half-open)."""
    iv = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            iv[f[0]].append((int(f[1]), int(f[2]), *f[3:]))
    return iv


def mask_array(chrom_sizes, intervals):
    """Boolean mask per chrom marking True inside intervals."""
    masks = {c: np.zeros(s, dtype=bool) for c, s in chrom_sizes.items()}
    for c, ivs in intervals.items():
        if c not in masks:
            continue
        a = masks[c]
        for s, e, *_ in ivs:
            s = max(0, s); e = min(a.size, e)
            if e > s:
                a[s:e] = True
    return masks


# ── 3. rDNA detection ────────────────────────────────────────────────────────
def detect_rdna(fa, out_bed):
    if Path(out_bed).exists():
        log(f"rDNA BED already present: {out_bed}")
        return out_bed
    if subprocess.run("command -v barrnap", shell=True, capture_output=True).returncode != 0:
        log("WARN: barrnap not found; skipping rDNA masking.")
        Path(out_bed).write_text("")
        return out_bed
    log("Running barrnap (kingdom=euk) for rDNA detection...")
    cmd = f"barrnap --kingdom euk --threads 4 {fa} 2>/dev/null"
    gff = subprocess.check_output(cmd, shell=True, text=True)
    with open(out_bed, "w") as fh:
        for ln in gff.splitlines():
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.split("\t")
            fh.write(f"{f[0]}\t{int(f[3])-1}\t{f[4]}\trDNA_{f[8].split(';')[0]}\n")
    log(f"rDNA BED written: {out_bed}")
    return out_bed


# ── 4. (A) H3K9me3 vs background occupancy ───────────────────────────────────
def analysis_A(occ_arrays, peaks_mask, valid_mask, outdir):
    log("(A) Occupancy in H3K9me3 vs rest …")
    in_vals  = []
    out_vals = []
    for c, a in occ_arrays.items():
        v = valid_mask[c]
        in_vals .append(a[ peaks_mask[c] & v])
        out_vals.append(a[~peaks_mask[c] & v])
    inv  = np.concatenate(in_vals)
    outv = np.concatenate(out_vals)

    # Sub-sample for the test (Mann-Whitney scales O(n log n) but bp = 27M)
    rng = np.random.default_rng(0)
    n_sub = 200_000
    inv_s  = rng.choice(inv,  min(n_sub, inv.size),  replace=False)
    outv_s = rng.choice(outv, min(n_sub, outv.size), replace=False)
    U, p = mannwhitneyu(inv_s, outv_s, alternative="greater")

    summary = pd.DataFrame({
        "set":      ["H3K9me3", "non-H3K9me3"],
        "bp":       [inv.size, outv.size],
        "mean_occ": [float(inv.mean()), float(outv.mean())],
        "median_occ":[float(np.median(inv)), float(np.median(outv))],
        "p25":      [float(np.percentile(inv,25)), float(np.percentile(outv,25))],
        "p75":      [float(np.percentile(inv,75)), float(np.percentile(outv,75))],
    })
    summary.to_csv(Path(outdir, "A_occupancy_H3K9_vs_bg.tsv"), sep="\t", index=False)
    with open(Path(outdir, "A_mannwhitney.txt"), "w") as fh:
        fh.write(f"Mann-Whitney U one-sided (H3K9me3 > non-H3K9me3)\n"
                 f"U = {U:.3e}\np = {p:.3e}\n"
                 f"n_subsample = {n_sub}\n")

    # Boxplot on per-1kb means (cleaner visualisation than per-bp)
    means_in, means_out = [], []
    for c, a in occ_arrays.items():
        n = a.size
        nbins = n // 1000
        if nbins == 0: continue
        a_b   = a[:nbins*1000].reshape(nbins, 1000).mean(axis=1)
        m_b   = peaks_mask[c][:nbins*1000].reshape(nbins, 1000).mean(axis=1) > 0.5
        v_b   = valid_mask[c][:nbins*1000].reshape(nbins, 1000).mean(axis=1) > 0.5
        means_in .append(a_b[ m_b & v_b])
        means_out.append(a_b[~m_b & v_b])
    mi = np.concatenate(means_in); mo = np.concatenate(means_out)
    fig, ax = plt.subplots(figsize=(4,4))
    ax.boxplot([mo, mi], labels=["non-H3K9me3", "H3K9me3"], showfliers=False)
    ax.set_ylabel("Mean nucleosome occupancy (1 kb bins)")
    ax.set_title(f"A) Occupancy bias\np={p:.1e} (MWU one-sided)")
    fig.tight_layout()
    fig.savefig(Path(outdir, "A_occupancy_boxplot.pdf"))
    plt.close(fig)
    log(f"(A) done: H3K9 mean={inv.mean():.2f}  non-H3K9 mean={outv.mean():.2f}  p={p:.2e}")


# ── 5. (B) Nucleosome arrangement metrics per region ────────────────────────
def autocorr_peak(sig, lag_min=140, lag_max=240):
    """Return (lag_at_max, max_corr) of normalised autocorrelation in [lag_min, lag_max]."""
    sig = sig - sig.mean()
    if sig.std() == 0:
        return np.nan, np.nan
    n = sig.size
    f = np.fft.rfft(sig, n=2*n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac /= ac[0]
    if lag_max >= n:
        return np.nan, np.nan
    seg = ac[lag_min:lag_max+1]
    i = int(np.argmax(seg))
    return lag_min + i, float(seg[i])


def window_metrics(occ, win=1000):
    """Per-window (mean, cv, n_peaks, ac_lag, ac_value) along one chromosome array."""
    n = occ.size
    nb = n // win
    if nb == 0:
        return np.empty((0,5))
    out = np.zeros((nb, 5), dtype=np.float32)
    for i in range(nb):
        s = occ[i*win:(i+1)*win]
        m = s.mean()
        sd = s.std()
        cv = sd / m if m > 1e-6 else 0.0
        # peaks ≥ mean+0.5sd, ≥147 bp apart, prominence 0.3*sd
        peaks, _ = find_peaks(s, height=m + 0.5*sd, distance=147, prominence=max(1.0, 0.3*sd))
        npk = len(peaks)
        lag, val = autocorr_peak(s, 140, 240)
        out[i] = (m, cv, npk, lag if lag==lag else 0, val if val==val else 0)
    return out


def analysis_B(occ_arrays, peaks_mask, valid_mask, outdir, win=1000):
    log("(B) Nucleosome arrangement (NRL, fuzziness, phasing) …")
    rows_in, rows_out = [], []
    for c, a in occ_arrays.items():
        n = a.size; nb = n // win
        if nb == 0: continue
        m  = window_metrics(a, win=win)              # (nb, 5)
        pm = peaks_mask[c][:nb*win].reshape(nb, win).mean(axis=1) > 0.5
        vm = valid_mask[c][:nb*win].reshape(nb, win).mean(axis=1) > 0.5
        for i in range(nb):
            if not vm[i]: continue
            row = (c, i*win, (i+1)*win, *m[i].tolist(), int(pm[i]))
            (rows_in if pm[i] else rows_out).append(row)
    cols = ["chrom","start","end","mean_occ","cv","n_peaks","ac_lag","ac_val","H3K9me3"]
    df = pd.DataFrame(rows_in + rows_out, columns=cols)
    df.to_csv(Path(outdir, "B_window_metrics_1kb.tsv.gz"), sep="\t", index=False,
              compression="gzip")

    summary = df.groupby("H3K9me3").agg(
        n=("mean_occ","size"),
        mean_occ=("mean_occ","mean"),
        mean_cv=("cv","mean"),
        median_npeaks=("n_peaks","median"),
        median_aclag=("ac_lag","median"),
        mean_acval=("ac_val","mean"),
    )
    summary.to_csv(Path(outdir, "B_arrangement_summary.tsv"), sep="\t")
    log("(B) summary:\n" + summary.to_string())

    # 4-panel comparison
    fig, axs = plt.subplots(1, 4, figsize=(14,3.4))
    metrics = [("mean_occ","Mean occupancy"),
               ("cv","Fuzziness (CV)"),
               ("n_peaks","# called peaks / kb"),
               ("ac_val","Phasing (autocorr 140-240 bp)")]
    for ax,(k,lbl) in zip(axs, metrics):
        d_o = df.loc[df.H3K9me3==0, k].values
        d_i = df.loc[df.H3K9me3==1, k].values
        if d_i.size > 5 and d_o.size > 5:
            _, p = mannwhitneyu(d_i, d_o, alternative="two-sided")
        else: p = np.nan
        ax.boxplot([d_o, d_i], labels=["bg","H3K9"], showfliers=False)
        ax.set_title(f"{lbl}\np={p:.1e}")
    fig.tight_layout()
    fig.savefig(Path(outdir, "B_arrangement_boxplots.pdf"))
    plt.close(fig)

    # NRL distribution histograms
    fig, ax = plt.subplots(figsize=(4.5,3.5))
    bins = np.arange(140, 241, 4)
    for label, sub, col in [("bg",     df[df.H3K9me3==0],"steelblue"),
                            ("H3K9",   df[df.H3K9me3==1],"crimson")]:
        v = sub.ac_lag.values; v = v[v>0]
        ax.hist(v, bins=bins, density=True, alpha=0.55, label=f"{label} (n={v.size})",
                color=col)
    ax.set_xlabel("Autocorrelation peak lag (≈ NRL, bp)")
    ax.set_ylabel("Density")
    ax.legend()
    fig.tight_layout()
    fig.savefig(Path(outdir, "B_NRL_histogram.pdf"))
    plt.close(fig)
    log("(B) done.")
    return df


# ── 6. (C) TE enrichment stratified by H3K9me3 ──────────────────────────────
def analysis_C(occ_arrays, te_intervals, peaks_mask, valid_mask,
               outdir, ndr_pct=25, nor_pct=75):
    log("(C) TE enrichment in NDR/INT/NOR, stratified by H3K9me3 …")
    # Global thresholds from valid genome
    vals = [a[valid_mask[c]] for c, a in occ_arrays.items()]
    vals = np.concatenate(vals)
    ndr_t = np.percentile(vals, ndr_pct)
    nor_t = np.percentile(vals, nor_pct)
    log(f"   thresholds (valid genome): NDR < {ndr_t:.2f}   NOR > {nor_t:.2f}")

    rows = []
    for stratum_name, peak_sel in [("non_H3K9me3", False), ("H3K9me3", True)]:
        # zone masks restricted to stratum
        zone_bp = {"NDR":0, "INT":0, "NOR":0}
        zone_mask = {}
        for c, a in occ_arrays.items():
            v = valid_mask[c]
            sel = peaks_mask[c] if peak_sel else ~peaks_mask[c]
            sel = sel & v
            ndr = (a < ndr_t) & sel
            nor_ = (a > nor_t) & sel
            inter = sel & ~ndr & ~nor_
            zone_mask[c] = {"NDR":ndr, "NOR":nor_, "INT":inter}
            zone_bp["NDR"] += int(ndr.sum())
            zone_bp["NOR"] += int(nor_.sum())
            zone_bp["INT"] += int(inter.sum())
        total_bp = sum(zone_bp.values())
        if total_bp == 0:
            continue
        # TE bp per superfamily per zone
        sf_zone_bp = defaultdict(lambda: {"NDR":0, "NOR":0, "INT":0})
        sf_total   = defaultdict(int)
        for c, ivs in te_intervals.items():
            if c not in zone_mask: continue
            for s, e, *rest in ivs:
                if not rest: continue
                sf = rest[0]
                s = max(0, s); e = min(occ_arrays[c].size, e)
                if e <= s: continue
                # restrict to stratum first
                seg_sel = (peaks_mask[c][s:e] if peak_sel else ~peaks_mask[c][s:e]) & valid_mask[c][s:e]
                ndr = zone_mask[c]["NDR"][s:e] & seg_sel
                nor_ = zone_mask[c]["NOR"][s:e] & seg_sel
                inter = zone_mask[c]["INT"][s:e] & seg_sel
                sf_zone_bp[sf]["NDR"]  += int(ndr.sum())
                sf_zone_bp[sf]["NOR"]  += int(nor_.sum())
                sf_zone_bp[sf]["INT"]  += int(inter.sum())
                sf_total[sf]           += int(seg_sel.sum())
        for sf, zd in sf_zone_bp.items():
            if sf_total[sf] < 200: continue
            for z, bp in zd.items():
                if zone_bp[z] == 0: continue
                obs = bp / sf_total[sf]
                exp = zone_bp[z] / total_bp
                enr = obs / exp if exp > 0 else np.nan
                log2 = math.log2(enr) if enr and enr > 0 else np.nan
                rows.append({"stratum":stratum_name, "superfamily":sf, "zone":z,
                             "te_bp":bp, "stratum_bp":zone_bp[z],
                             "te_total_in_stratum":sf_total[sf],
                             "enrichment":enr, "log2_enrichment":log2})
    df = pd.DataFrame(rows)
    df.to_csv(Path(outdir, "C_TE_enrichment_by_H3K9.tsv"), sep="\t", index=False)

    # Plot: NOR log2 enrichment per superfamily, two strata
    sub = df[df.zone=="NOR"].pivot(index="superfamily", columns="stratum",
                                   values="log2_enrichment").dropna()
    if not sub.empty:
        fig, ax = plt.subplots(figsize=(4.5,4.5))
        ax.scatter(sub["non_H3K9me3"], sub["H3K9me3"], s=18, c="darkslateblue")
        lim = np.nanmax(np.abs(sub.values)) + 0.2
        ax.plot([-lim,lim],[-lim,lim],"--",c="grey",lw=0.7)
        ax.axhline(0,c="grey",lw=0.4); ax.axvline(0,c="grey",lw=0.4)
        ax.set_xlabel("log2(NOR enrichment) | non-H3K9me3")
        ax.set_ylabel("log2(NOR enrichment) | H3K9me3")
        ax.set_title("C) TE→NOR bias inside vs outside H3K9me3")
        fig.tight_layout()
        fig.savefig(Path(outdir,"C_NOR_enrichment_scatter.pdf"))
        plt.close(fig)
    log("(C) done.")
    return df


# ── 7. (D) Phased / fuzzy classification + TE enrichment ────────────────────
def analysis_D(occ_arrays, te_intervals, peaks_mask, valid_mask, win_metrics_df,
               outdir, win=1000):
    log("(D) Classifying windows into Phased/Fuzzy x Occupied/Depleted …")
    df = win_metrics_df.copy()
    # Median splits on valid genome
    occ_med = df.mean_occ.median()
    pos_med = df.ac_val.median()    # higher ac_val = more phased
    log(f"   medians:  occ={occ_med:.2f}  ac_val={pos_med:.3f}")

    def classify(r):
        hi_occ = r.mean_occ >  occ_med
        hi_pos = r.ac_val   >  pos_med
        if hi_occ and  hi_pos: return "Phased-occupied"
        if hi_occ and not hi_pos: return "Fuzzy-occupied"
        if (not hi_occ) and hi_pos: return "Phased-depleted"
        return "Open-fuzzy"
    df["state"] = df.apply(classify, axis=1)
    df.to_csv(Path(outdir,"D_window_states_1kb.tsv.gz"), sep="\t", index=False,
              compression="gzip")

    # Build per-chrom state arrays at bp resolution
    state_levels = ["Phased-occupied","Fuzzy-occupied","Phased-depleted","Open-fuzzy"]
    code = {s:i+1 for i,s in enumerate(state_levels)}      # 0 = invalid
    state_arr = {c: np.zeros(a.size, dtype=np.uint8) for c, a in occ_arrays.items()}
    for r in df.itertuples():
        a = state_arr[r.chrom]
        s = int(r.start); e = min(int(r.end), a.size)
        a[s:e] = code[r.state]

    # Enrichment per superfamily per state (only on valid genome)
    state_bp = {s:0 for s in state_levels}
    for c, sa in state_arr.items():
        v = valid_mask[c]
        for s, k in code.items():
            state_bp[s] += int(((sa==k) & v).sum())
    total_bp = sum(state_bp.values())

    sf_state_bp = defaultdict(lambda: {s:0 for s in state_levels})
    sf_total    = defaultdict(int)
    for c, ivs in te_intervals.items():
        if c not in state_arr: continue
        for s_, e_, *rest in ivs:
            if not rest: continue
            sf = rest[0]
            s_ = max(0, s_); e_ = min(state_arr[c].size, e_)
            if e_ <= s_: continue
            seg = state_arr[c][s_:e_]
            v   = valid_mask[c][s_:e_]
            for st, k in code.items():
                sf_state_bp[sf][st] += int(((seg==k) & v).sum())
            sf_total[sf] += int(v.sum())

    rows = []
    for sf, zd in sf_state_bp.items():
        if sf_total[sf] < 200: continue
        for st, bp in zd.items():
            if state_bp[st] == 0: continue
            obs = bp / sf_total[sf]
            exp = state_bp[st] / total_bp
            enr = obs/exp if exp > 0 else np.nan
            log2 = math.log2(enr) if enr and enr > 0 else np.nan
            rows.append({"superfamily":sf,"state":st,"te_bp":bp,
                         "state_bp":state_bp[st],"te_total_bp":sf_total[sf],
                         "enrichment":enr,"log2_enrichment":log2})
    enr_df = pd.DataFrame(rows)
    enr_df.to_csv(Path(outdir,"D_TE_enrichment_per_state.tsv"), sep="\t", index=False)

    # Cross-tab state × H3K9me3 (genome bp)
    cross = []
    for st, k in code.items():
        for hk_label, hk_sel in [("H3K9me3",True),("non_H3K9me3",False)]:
            bp = 0
            for c, sa in state_arr.items():
                pk = peaks_mask[c] if hk_sel else ~peaks_mask[c]
                bp += int(((sa==k) & pk & valid_mask[c]).sum())
            cross.append({"state":st,"stratum":hk_label,"bp":bp})
    pd.DataFrame(cross).to_csv(Path(outdir,"D_state_vs_H3K9_bp.tsv"), sep="\t", index=False)

    # Heatmap of log2 enrichment, top-25 superfamilies by total TE bp
    top_sf = (enr_df.groupby("superfamily").te_total_bp.first()
                    .sort_values(ascending=False).head(25).index)
    pivot = (enr_df[enr_df.superfamily.isin(top_sf)]
             .pivot(index="superfamily", columns="state", values="log2_enrichment")
             .reindex(columns=state_levels))
    fig, ax = plt.subplots(figsize=(6, max(4, 0.3*len(top_sf))))
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(state_levels))); ax.set_xticklabels(state_levels, rotation=30, ha="right")
    ax.set_yticks(range(pivot.shape[0])); ax.set_yticklabels(pivot.index, fontsize=7)
    fig.colorbar(im, ax=ax, label="log2 enrichment")
    ax.set_title("D) TE enrichment per nucleosome state\n(top-25 superfamilies)")
    fig.tight_layout()
    fig.savefig(Path(outdir,"D_TE_enrichment_heatmap.pdf"))
    plt.close(fig)

    # State composition pie + state×H3K9 stacked bar
    fig, axes = plt.subplots(1, 2, figsize=(9,4))
    sizes = [state_bp[s] for s in state_levels]
    axes[0].pie(sizes, labels=[f"{s}\n{b/1e6:.1f} Mb" for s,b in zip(state_levels,sizes)],
                autopct="%1.0f%%", startangle=90,
                colors=["#c0392b","#e67e22","#27ae60","#3498db"])
    axes[0].set_title("Genome bp per state")
    cdf = pd.DataFrame(cross).pivot(index="state",columns="stratum",values="bp")\
                              .reindex(state_levels)
    cdf_pct = cdf.div(cdf.sum(axis=1), axis=0) * 100
    cdf_pct.plot(kind="barh", stacked=True, ax=axes[1],
                 color=["#7d3c98","#bdc3c7"])
    axes[1].set_xlabel("% bp"); axes[1].set_title("H3K9me3 fraction per state")
    fig.tight_layout()
    fig.savefig(Path(outdir,"D_state_overview.pdf"))
    plt.close(fig)
    log("(D) done.")
    return enr_df


# ── 8. main ──────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wig",   default=DEF_WIG)
    ap.add_argument("--peaks", default=DEF_PEAKS)
    ap.add_argument("--te",    default=DEF_TEBED)
    ap.add_argument("--fai",   default=DEF_FAI)
    ap.add_argument("--fa",    default=DEF_FA)
    ap.add_argument("--outdir",default=DEF_OUT)
    ap.add_argument("--win",   type=int, default=1000)
    ap.add_argument("--skip_rdna", action="store_true")
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    log("Loading chrom sizes...")
    sizes = load_chrom_sizes(args.fai)
    log(f"  {len(sizes)} scaffolds, total {sum(sizes.values())/1e6:.2f} Mb")

    log("Loading WIG...")
    occ = load_wig_arrays(args.wig, sizes)

    log("Loading H3K9me3 broad peaks...")
    peaks = load_bed_intervals(args.peaks)
    pmask = mask_array(sizes, peaks)
    n_peak_bp = sum(int(m.sum()) for m in pmask.values())
    log(f"  {sum(len(v) for v in peaks.values())} peaks, {n_peak_bp/1e6:.2f} Mb")

    log("Loading TE BED...")
    te = load_bed_intervals(args.te)
    n_te = sum(len(v) for v in te.values())
    log(f"  {n_te} TE features")

    # Validity mask = positions with WIG coverage AND not in rDNA
    valid = {c: (a > 0) for c, a in occ.items()}
    n_valid_pre = sum(int(m.sum()) for m in valid.values())

    if not args.skip_rdna:
        rdna_bed = Path(args.outdir, "rdna_barrnap.bed")
        detect_rdna(args.fa, str(rdna_bed))
        if rdna_bed.stat().st_size > 0:
            rdna_iv  = load_bed_intervals(str(rdna_bed))
            rdna_msk = mask_array(sizes, rdna_iv)
            for c in valid:
                valid[c] &= ~rdna_msk[c]
    n_valid = sum(int(m.sum()) for m in valid.values())
    log(f"Valid bp: {n_valid_pre/1e6:.2f} Mb → {n_valid/1e6:.2f} Mb after rDNA mask")

    # ─── A ───
    analysis_A(occ, pmask, valid, args.outdir)

    # ─── B ───
    win_df = analysis_B(occ, pmask, valid, args.outdir, win=args.win)

    # ─── C ───
    analysis_C(occ, te, pmask, valid, args.outdir)

    # ─── D ───
    analysis_D(occ, te, pmask, valid, win_df, args.outdir, win=args.win)

    log("ALL DONE.")


if __name__ == "__main__":
    main()
