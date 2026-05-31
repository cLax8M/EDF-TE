#!/usr/bin/env python3
"""
chromatin_states_from_nuclpred.py
=================================

Definir estados de cromatina en CUALQUIER hongo a partir SÓLO del track
de ocupación nucleosomal predicha (nuclpred WIG) y luego analizar la
distribución/preferencias de TEs.

Calibración (Rhizopus microsporus): demostramos previamente que las
regiones de H3K9me3 son densas pero pobremente fasadas (alto mean_occ +
bajo ac_val). Por lo tanto el estado nucleosómico
``Fuzzy_high`` se interpreta como **heterochromatin-like (proxy)**.

Estados (4):
    Open          -> mean_occ < q33                    (cromatina abierta)
    Phased        -> mean_occ >= q33 and ac_val >= med  (eucromatina activa)
    Fuzzy_mid     -> q33 <= mean_occ < q67 and ac_val < med
    Fuzzy_high    -> mean_occ >= q67 and ac_val < med   (HETEROCHROMATIN proxy)

Entradas por especie:
    --species   código (e.g. MU_Rhimi59); se usa para localizar archivos
    --nuclpred  WIG fixedStep step=1 (chrom names PUEDEN llevar prefijo)
    --fai       genome.fasta.fai
    --te        BED de TEs (chrom start end superfam ...)
    --outdir    directorio salida

Salidas:
    chromatin_states.tsv.gz, chromatin_states_summary.tsv,
    TE_per_chromatin_state.tsv, TE_enrichment_per_state.tsv,
    overview.pdf, TE_heatmap.pdf
"""
import argparse, os, sys, gzip, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

# -------------------- WIG loading --------------------
def load_nuclpred_wig(path, chrom_sizes):
    """Lee fixedStep o variableStep WIG. Acepta nombres de cromosoma con prefijo
    'SPECIES.scaffold_X' o 'scaffold_X'. Devuelve dict chrom->np.float32."""
    arrs = {c: np.zeros(L, dtype=np.float32) for c, L in chrom_sizes.items()}
    cur = None; pos = 0; step = 1; span = 1; mode = None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for ln in fh:
            if not ln: continue
            if ln.startswith("fixedStep"):
                kv = dict(tok.split("=") for tok in ln.split()[1:] if "=" in tok)
                chrom = kv.get("chrom","")
                if "." in chrom and chrom not in chrom_sizes:
                    chrom = chrom.split(".",1)[1]
                cur = chrom if chrom in chrom_sizes else None
                pos = int(kv.get("start",1)) - 1
                step = int(kv.get("step",1))
                span = int(kv.get("span",1))
                mode = "fixed"
                continue
            if ln.startswith("variableStep"):
                kv = dict(tok.split("=") for tok in ln.split()[1:] if "=" in tok)
                chrom = kv.get("chrom","")
                if "." in chrom and chrom not in chrom_sizes:
                    chrom = chrom.split(".",1)[1]
                cur = chrom if chrom in chrom_sizes else None
                span = int(kv.get("span",1))
                mode = "var"
                continue
            if cur is None: continue
            if mode == "fixed":
                try: v = float(ln)
                except ValueError: continue
                if pos < len(arrs[cur]):
                    end = min(pos+span, len(arrs[cur]))
                    arrs[cur][pos:end] = v
                pos += step
            elif mode == "var":
                parts = ln.split()
                if len(parts) < 2: continue
                try:
                    p = int(parts[0]) - 1
                    v = float(parts[1])
                except ValueError:
                    continue
                if 0 <= p < len(arrs[cur]):
                    end = min(p+span, len(arrs[cur]))
                    arrs[cur][p:end] = v
    return arrs

# -------------------- per-window metrics --------------------
def autocorr_peak(sig, lag_min=140, lag_max=240):
    """Devuelve (lag_pico, valor_normalizado) en el rango lag_min..lag_max."""
    sig = sig - sig.mean()
    if sig.std() < 1e-6 or len(sig) < lag_max+10:
        return np.nan, 0.0
    n = len(sig)
    # FFT autocorr
    f = np.fft.rfft(sig, n=2*n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac = ac / ac[0]
    seg = ac[lag_min:lag_max+1]
    if len(seg) == 0: return np.nan, 0.0
    i = int(np.argmax(seg))
    return lag_min + i, float(seg[i])

def window_metrics(occ, win):
    L = len(occ)
    n = L // win
    means = np.zeros(n, dtype=np.float32)
    cvs   = np.zeros(n, dtype=np.float32)
    aclags= np.zeros(n, dtype=np.float32)
    acvals= np.zeros(n, dtype=np.float32)
    valids= np.zeros(n, dtype=bool)
    for i in range(n):
        s = i*win; e = s+win
        seg = occ[s:e]
        # validez: la ventana entera está dentro de un cromosoma (siempre cierto aquí)
        # Para nuclpred (track predicho), aceptamos ventanas con cualquier señal
        if seg.size < win: continue
        m = seg.mean()
        sd = seg.std()
        means[i] = m
        cvs[i] = sd / (m + 1e-6)
        lag, v = autocorr_peak(seg)
        aclags[i] = lag if lag is not None and not np.isnan(lag) else 0
        acvals[i] = v
        valids[i] = True
    return means, cvs, aclags, acvals, valids

# -------------------- state assignment --------------------
def assign_states_nuc(df, low_q=0.33, high_q=0.67, min_occ=1e-3):
    occ = df["mean_occ"].values
    ac  = df["ac_val"].values
    # ventanas informativas: con al menos un poco de señal predicha
    informative = occ > min_occ
    if informative.sum() < 100:
        # fallback: usar todas
        informative = np.ones_like(occ, dtype=bool)
    occ_lo, occ_hi = np.quantile(occ[informative], [low_q, high_q])
    ac_med         = np.median(ac[informative])
    state = np.full(len(df), "Fuzzy_mid", dtype=object)
    state[~informative]                                      = "No_signal"
    state[informative & (occ < occ_lo)]                      = "Open"
    state[informative & (occ >= occ_lo) & (ac >= ac_med)]    = "Phased"
    state[informative & (occ >= occ_hi) & (ac < ac_med)]     = "Fuzzy_high"
    df = df.copy()
    df["state"] = state
    return df, dict(occ_lo=float(occ_lo), occ_hi=float(occ_hi),
                    ac_med=float(ac_med),
                    n_informative=int(informative.sum()))

STATE_ORDER  = ["No_signal", "Open", "Phased", "Fuzzy_mid", "Fuzzy_high"]
STATE_COLORS = {"No_signal":"#CCCCCC", "Open":"#6BAED6", "Phased":"#74C476",
                "Fuzzy_mid":"#FDAE6B", "Fuzzy_high":"#E6550D"}

# -------------------- TE intersection --------------------
def intersect_te(te_df, win_df, win):
    bp_mat = {}; cnt_mat = {}
    chrom_starts = {c: win_df.loc[win_df["chrom"]==c,"start"].values for c in win_df["chrom"].unique()}
    chrom_states = {c: win_df.loc[win_df["chrom"]==c,"state"].values for c in win_df["chrom"].unique()}
    for row in te_df.itertuples(index=False):
        c, s, e, fam = row.chrom, int(row.start), int(row.end), row.superfam
        if c not in chrom_starts: continue
        starts = chrom_starts[c]; states = chrom_states[c]
        i0 = max(0, np.searchsorted(starts, s, side="right") - 1)
        i1 = np.searchsorted(starts, e, side="left")
        for i in range(i0, i1):
            ws = starts[i]; we = ws + win
            ov = max(0, min(e, we) - max(s, ws))
            if ov <= 0: continue
            st = states[i]
            bp_mat[(st, fam)] = bp_mat.get((st, fam), 0) + ov
        mid = (s+e)//2
        i = np.searchsorted(starts, mid, side="right") - 1
        if 0 <= i < len(starts) and starts[i] <= mid < starts[i]+win:
            st = states[i]
            cnt_mat[(st, fam)] = cnt_mat.get((st, fam), 0) + 1
    return bp_mat, cnt_mat

def enrichment_table(bp_mat, state_bp, fam_bp_total, total_valid_bp):
    rows = []
    for st in sorted(state_bp):
        for fa in sorted(fam_bp_total):
            obs = bp_mat.get((st, fa), 0)
            sbp = state_bp[st]
            fbp = fam_bp_total[fa]
            log2or = np.log2(((obs+1)/(sbp+1)) / ((fbp+1)/(total_valid_bp+1)))
            a = obs; b = sbp - obs
            c = fbp - obs; d = (total_valid_bp - sbp) - (fbp - obs)
            try:
                if min(a,b,c,d) >= 0 and (a+b)>0 and (c+d)>0:
                    _, p, *_ = chi2_contingency([[a,b],[c,d]])
                else: p = np.nan
            except Exception: p = np.nan
            rows.append(dict(state=st, superfam=fa, obs_bp=obs,
                             exp_bp=sbp*fbp/total_valid_bp if total_valid_bp else np.nan,
                             log2OR=log2or, p=p,
                             state_bp=sbp, fam_bp_total=fbp))
    return pd.DataFrame(rows)

# -------------------- main --------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--nuclpred", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--te", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--win", type=int, default=500)
    ap.add_argument("--top_fams", type=int, default=15)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[{args.species}] chrom sizes ...")
    fai = pd.read_csv(args.fai, sep="\t", header=None,
                      names=["chrom","length","off","lineb","lineby"])
    chrom_sizes = dict(zip(fai["chrom"], fai["length"]))
    total_len = sum(chrom_sizes.values())
    print(f"[{args.species}]   {len(chrom_sizes)} scaffolds, {total_len/1e6:.2f} Mb")

    print(f"[{args.species}] loading nuclpred WIG ...")
    occ = load_nuclpred_wig(args.nuclpred, chrom_sizes)
    nz = sum(int((a>0).sum()) for a in occ.values())
    print(f"[{args.species}]   non-zero positions: {nz/1e6:.2f} M")

    print(f"[{args.species}] computing window metrics ({args.win} bp) ...")
    rows = []
    for c, L in chrom_sizes.items():
        m, cv, lg, av, v = window_metrics(occ[c], args.win)
        n = len(m)
        for i in range(n):
            if not v[i]: continue
            rows.append((c, i*args.win, i*args.win+args.win,
                         m[i], cv[i], lg[i], av[i]))
    win_df = pd.DataFrame(rows, columns=["chrom","start","end","mean_occ","cv","ac_lag","ac_val"])
    print(f"[{args.species}]   valid windows: {len(win_df):,}")

    win_df, thr = assign_states_nuc(win_df)
    print(f"[{args.species}]   thresholds: {thr}")
    win_df.to_csv(os.path.join(args.outdir, "chromatin_states.tsv.gz"),
                  sep="\t", index=False, compression="gzip")
    with open(os.path.join(args.outdir, "thresholds.json"), "w") as fh:
        json.dump({**thr, "win": args.win, "species": args.species}, fh, indent=2)

    summary = (win_df.groupby("state").size().rename("n_windows").to_frame())
    summary["bp"] = summary["n_windows"] * args.win
    summary["pct_bp"] = 100*summary["bp"]/summary["bp"].sum()
    summary = summary.reindex([s for s in STATE_ORDER if s in summary.index])
    summary.to_csv(os.path.join(args.outdir, "chromatin_states_summary.tsv"), sep="\t")
    print(summary)

    # TE
    print(f"[{args.species}] loading TE BED ...")
    te = pd.read_csv(args.te, sep="\t", header=None,
                     names=["chrom","start","end","superfam","score","strand"])
    te = te[te["chrom"].isin(chrom_sizes)]
    print(f"[{args.species}]   {len(te):,} TE features")
    bp_mat, cnt_mat = intersect_te(te, win_df, args.win)

    state_bp = summary["bp"].to_dict()
    fam_bp_total = {}
    for (st, fa), bp in bp_mat.items():
        fam_bp_total[fa] = fam_bp_total.get(fa, 0) + bp
    total_valid_bp = sum(state_bp.values())

    rows = [dict(state=st, superfam=fa, bp=bp, n_TE=cnt_mat.get((st,fa),0))
            for (st,fa), bp in bp_mat.items()]
    te_long = pd.DataFrame(rows).sort_values(["state","bp"], ascending=[True,False])
    te_long.to_csv(os.path.join(args.outdir, "TE_per_chromatin_state.tsv"), sep="\t", index=False)

    enr = enrichment_table(bp_mat, state_bp, fam_bp_total, total_valid_bp)
    enr.to_csv(os.path.join(args.outdir, "TE_enrichment_per_state.tsv"), sep="\t", index=False)

    # ---- overview pie ----
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    vals = summary["bp"].values
    labs = [f"{s}\n{int(b/1e3)} kb ({p:.1f}%)" for s, b, p in
            zip(summary.index, vals, summary["pct_bp"])]
    ax.pie(vals, labels=labs,
           colors=[STATE_COLORS.get(s,"#888") for s in summary.index],
           startangle=90, wedgeprops=dict(edgecolor="w", linewidth=1.2),
           textprops=dict(fontsize=8))
    ax.set_title(f"{args.species} chromatin states ({args.win} bp)\n"
                 f"Σ classified = {vals.sum()/1e6:.2f} Mb of {total_len/1e6:.2f} Mb")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "chromatin_states_overview.pdf"))
    plt.close(fig)

    # ---- TE heatmap (top families × state) ----
    top_fams = (te_long.groupby("superfam")["bp"].sum()
                .sort_values(ascending=False).head(args.top_fams).index.tolist())
    states_present = [s for s in STATE_ORDER if s in enr["state"].unique()]
    pivot = (enr[enr["superfam"].isin(top_fams)]
             .pivot(index="superfam", columns="state", values="log2OR")
             .reindex(index=top_fams).reindex(columns=states_present))
    if pivot.size:
        fig, ax = plt.subplots(figsize=(5.5, 0.4*len(top_fams)+1.5))
        vmax = np.nanmax(np.abs(pivot.values))
        im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(pivot.shape[1])); ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(pivot.shape[0])); ax.set_yticklabels(pivot.index, fontsize=8)
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = pivot.values[i,j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6,
                            color="white" if abs(v)>0.6*vmax else "black")
        cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
        cb.set_label("log2(obs/exp) TE bp", fontsize=8)
        ax.set_title(f"{args.species}: TE preference per state")
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, "TE_chromatin_heatmap.pdf"))
        plt.close(fig)

    print(f"[{args.species}] DONE -> {args.outdir}")

if __name__ == "__main__":
    main()
