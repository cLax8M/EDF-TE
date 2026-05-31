#!/usr/bin/env python3
"""
chromatin_states_TE.py
======================

Definir tipos de cromatina en Rhizopus microsporus combinando:
  - H3K9me3 (broad peak: 0/1)
  - Ocupación nucleosomal (mean_occ → tertiles low/mid/high)
  - Fasing (ac_val ≥ mediana → phased ; resto → fuzzy)
y, después, analizar la distribución/preferencia de TEs por estado.

Entrada: B_window_metrics_1kb.tsv.gz (en realidad 500 bp) producido por
H3K9_nucleosome_interplay.py, y la BED de TEs filtrados.

Salidas (en --outdir):
  chromatin_states.tsv.gz         : ventanas con estado y subestado
  chromatin_states.bed            : BED itemRgb (cargable en IGV/JBrowse)
  chromatin_states_summary.tsv    : nº ventanas y bp por estado
  TE_per_chromatin_state.tsv      : conteos y bp de TE por (estado × superfamilia)
  TE_enrichment_per_state.tsv     : log2(obs/esp) y chi² por celda
  chromatin_states_overview.pdf   : composición global + por scaffold
  TE_chromatin_heatmap.pdf        : heatmap log2OR superfamilia × estado
  TE_chromatin_stacked.pdf        : composición de TE dentro de cada estado
"""
import argparse, os, sys, gzip
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

# ---------- chromatin state definition ----------
def assign_states(df, occ_low_q=0.33, occ_high_q=0.67):
    occ = df["mean_occ"].values
    ac  = df["ac_val"].values
    h3  = df["H3K9me3"].values.astype(int)

    occ_lo, occ_hi = np.quantile(occ, [occ_low_q, occ_high_q])
    ac_med         = np.median(ac)

    occ_cat = np.where(occ < occ_lo, "low",
              np.where(occ > occ_hi, "high", "mid"))
    phasing = np.where(ac >= ac_med, "phased", "fuzzy")

    # Final chromatin state (5 categorías biológicamente interpretables)
    state = np.full(len(df), "Intermediate", dtype=object)
    # 1. Heterochromatin H3K9me3 (compacto, denso, pobremente posicionado)
    state[(h3 == 1) & (occ_cat == "high")]                   = "Heterochromatin_H3K9me3"
    state[(h3 == 1) & (occ_cat != "high")]                   = "H3K9me3_lowOcc"
    # 2. Heterochromatin-like sin H3K9me3 (denso pero fuzzy → posiblemente H3K27me/otra marca)
    state[(h3 == 0) & (occ_cat == "high") & (phasing == "fuzzy")]  = "Dense_Fuzzy_noH3K9"
    # 3. Eucromatina activa: nucleosomas bien fasados
    state[(h3 == 0) & (occ_cat == "mid")  & (phasing == "phased")] = "Euchromatin_Phased"
    state[(h3 == 0) & (occ_cat == "high") & (phasing == "phased")] = "Euchromatin_Phased"
    # 4. Open / NDR-like (poca densidad)
    state[(h3 == 0) & (occ_cat == "low")]                          = "Open_NDR"

    df = df.copy()
    df["occ_cat"] = occ_cat
    df["phasing"] = phasing
    df["state"]   = state
    return df, dict(occ_lo=occ_lo, occ_hi=occ_hi, ac_med=ac_med)

STATE_ORDER = ["Open_NDR",
               "Euchromatin_Phased",
               "Intermediate",
               "Dense_Fuzzy_noH3K9",
               "H3K9me3_lowOcc",
               "Heterochromatin_H3K9me3"]

STATE_COLORS = {
    "Open_NDR":                "#6BAED6",   # azul
    "Euchromatin_Phased":      "#74C476",   # verde
    "Intermediate":            "#BDBDBD",   # gris
    "Dense_Fuzzy_noH3K9":      "#FDAE6B",   # naranja
    "H3K9me3_lowOcc":          "#9E9AC8",   # lila
    "Heterochromatin_H3K9me3": "#E6550D",   # rojo
}

def rgb_str(hex_color):
    h = hex_color.lstrip("#")
    return f"{int(h[0:2],16)},{int(h[2:4],16)},{int(h[4:6],16)}"

# ---------- TE intersection (pure pandas, no bedtools) ----------
def intersect_te_with_windows(te_df, win_df, win):
    """Asigna a cada ventana la suma de bp y el conteo de TE por superfamilia."""
    # te_df: chrom,start,end,superfam,score,strand
    # win_df: chrom,start,end,..,state
    # Para cada superfamilia -> bp por ventana usando aproximación: asignar TE a ventana(s) cuyo solapamiento >0
    out_rows = []
    win_df = win_df.reset_index(drop=True)
    win_df["wid"] = np.arange(len(win_df))
    win_idx = {(c,s): wid for c,s,wid in zip(win_df["chrom"], win_df["start"], win_df["wid"])}

    # vector form
    chrom_arr = win_df["chrom"].values
    start_arr = win_df["start"].values
    state_arr = win_df["state"].values

    bp_mat   = {}   # (state, fam) -> bp
    cnt_mat  = {}   # (state, fam) -> n_TEs (asigna al window que contiene su midpoint)

    # build per-chrom window lookup as sorted arrays
    chroms = win_df["chrom"].unique()
    chrom_starts = {c: win_df.loc[win_df["chrom"]==c, "start"].values for c in chroms}
    chrom_states = {c: win_df.loc[win_df["chrom"]==c, "state"].values for c in chroms}

    for row in te_df.itertuples(index=False):
        c, s, e, fam = row.chrom, int(row.start), int(row.end), row.superfam
        if c not in chrom_starts: continue
        starts = chrom_starts[c]; states = chrom_states[c]
        # find windows overlapping [s,e)
        i0 = np.searchsorted(starts, s, side="right") - 1
        i1 = np.searchsorted(starts, e, side="left")
        if i0 < 0: i0 = 0
        for i in range(i0, i1):
            ws = starts[i]; we = ws + win
            ov = max(0, min(e, we) - max(s, ws))
            if ov <= 0: continue
            st = states[i]
            bp_mat[(st, fam)] = bp_mat.get((st, fam), 0) + ov
        # midpoint count
        mid = (s + e) // 2
        i = np.searchsorted(starts, mid, side="right") - 1
        if 0 <= i < len(starts) and starts[i] <= mid < starts[i] + win:
            st = states[i]
            cnt_mat[(st, fam)] = cnt_mat.get((st, fam), 0) + 1

    return bp_mat, cnt_mat

# ---------- enrichment ----------
def enrichment_table(bp_mat, state_bp, fam_bp_total, total_valid_bp):
    """log2((te_bp_in_state / state_bp) / (te_bp_total / total_valid_bp))
    + chi² 2x2 (state vs rest, te vs rest_te)."""
    states = sorted(state_bp)
    fams   = sorted(fam_bp_total)
    rows = []
    for st in states:
        for fa in fams:
            obs   = bp_mat.get((st, fa), 0)
            sbp   = state_bp[st]
            fbp   = fam_bp_total[fa]
            exp   = sbp * fbp / total_valid_bp if total_valid_bp > 0 else np.nan
            log2or = np.log2(((obs + 1) / (sbp + 1)) /
                             ((fbp + 1) / (total_valid_bp + 1)))
            # chi² 2x2: te-in-state, non-te-in-state, te-in-rest, non-te-in-rest
            a = obs
            b = sbp - obs
            c2 = fbp - obs
            d  = (total_valid_bp - sbp) - (fbp - obs)
            try:
                if min(a, b, c2, d) >= 0 and (a+b)>0 and (c2+d)>0:
                    chi2, p, *_ = chi2_contingency([[a, b], [c2, d]])
                else:
                    p = np.nan
            except Exception:
                p = np.nan
            rows.append(dict(state=st, superfam=fa,
                             obs_bp=obs, exp_bp=exp,
                             log2OR=log2or, p=p,
                             state_bp=sbp, fam_bp_total=fbp))
    return pd.DataFrame(rows)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", default="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes_H3K9_500/B_window_metrics_1kb.tsv.gz")
    ap.add_argument("--te",      default="/data/carlos/EDF-TE/TE_analysis/annotated/MU_Rhimi59_EarlGrey/MU_Rhimi59_summaryFiles/MU_Rhimi59.filteredRepeats.bed")
    ap.add_argument("--win",     type=int, default=500)
    ap.add_argument("--outdir",  default="/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/chromatin_states_500")
    ap.add_argument("--top_fams", type=int, default=15)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[CHRSTATE] Loading window metrics: {args.metrics}")
    win_df = pd.read_csv(args.metrics, sep="\t")
    print(f"[CHRSTATE]   {len(win_df):,} windows")
    win_df, thr = assign_states(win_df)
    print(f"[CHRSTATE] thresholds: occ_lo={thr['occ_lo']:.1f}  occ_hi={thr['occ_hi']:.1f}  ac_med={thr['ac_med']:.4f}")

    # save annotated windows
    win_out = os.path.join(args.outdir, "chromatin_states.tsv.gz")
    win_df.to_csv(win_out, sep="\t", index=False, compression="gzip")
    # BED itemRgb
    bed_path = os.path.join(args.outdir, "chromatin_states.bed")
    with open(bed_path, "w") as fh:
        for r in win_df.itertuples(index=False):
            col = STATE_COLORS.get(r.state, "#888888")
            fh.write(f"{r.chrom}\t{r.start}\t{r.end}\t{r.state}\t0\t.\t{r.start}\t{r.end}\t{rgb_str(col)}\n")

    # state summary
    summary = (win_df.groupby("state").size().rename("n_windows").to_frame())
    summary["bp"] = summary["n_windows"] * args.win
    summary["pct_bp"] = 100 * summary["bp"] / summary["bp"].sum()
    summary = summary.reindex([s for s in STATE_ORDER if s in summary.index])
    summary.to_csv(os.path.join(args.outdir, "chromatin_states_summary.tsv"), sep="\t")
    print("[CHRSTATE] state composition (bp):")
    print(summary)

    # load TEs
    print("[CHRSTATE] Loading TE BED ...")
    te = pd.read_csv(args.te, sep="\t", header=None,
                     names=["chrom","start","end","superfam","score","strand"])
    # group rare families: keep all; collapse '/' subdivision? mantenemos superfamilia tal cual
    print(f"[CHRSTATE]   {len(te):,} TE features, {te['superfam'].nunique()} superfamilias")

    # intersect
    print("[CHRSTATE] Intersecting TEs with chromatin windows ...")
    bp_mat, cnt_mat = intersect_te_with_windows(te, win_df, args.win)

    # totals
    state_bp = summary["bp"].to_dict()
    fam_bp_total = {}
    for (st, fa), bp in bp_mat.items():
        fam_bp_total[fa] = fam_bp_total.get(fa, 0) + bp
    total_valid_bp = sum(state_bp.values())

    # TE per state long table
    rows = []
    for (st, fa), bp in bp_mat.items():
        rows.append(dict(state=st, superfam=fa, bp=bp,
                         n_TE=cnt_mat.get((st, fa), 0)))
    te_long = pd.DataFrame(rows).sort_values(["state","bp"], ascending=[True, False])
    te_long.to_csv(os.path.join(args.outdir, "TE_per_chromatin_state.tsv"), sep="\t", index=False)

    # enrichment
    enr = enrichment_table(bp_mat, state_bp, fam_bp_total, total_valid_bp)
    enr.to_csv(os.path.join(args.outdir, "TE_enrichment_per_state.tsv"), sep="\t", index=False)

    # ---- plots ----
    # 1. chromatin overview pie
    fig, ax = plt.subplots(figsize=(6.5, 6))
    vals = summary["bp"].values
    labs = [f"{s}\n{int(b/1e3)} kb ({p:.1f}%)" for s, b, p in
            zip(summary.index, vals, summary["pct_bp"])]
    cols = [STATE_COLORS.get(s, "#888") for s in summary.index]
    ax.pie(vals, labels=labs, colors=cols, startangle=90,
           wedgeprops=dict(edgecolor="w", linewidth=1.3),
           textprops=dict(fontsize=8))
    ax.set_title(f"Chromatin states (Rhizopus microsporus, {args.win} bp windows)\n"
                 f"Σ classified = {int(vals.sum()/1e6):.2f} Mb")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "chromatin_states_overview.pdf"))
    plt.close(fig)

    # 2. per-scaffold composition (top 25 scaffolds by length)
    sc_len = win_df.groupby("chrom").size().sort_values(ascending=False).head(25).index
    sub = win_df[win_df["chrom"].isin(sc_len)]
    comp = (sub.groupby(["chrom","state"]).size().unstack(fill_value=0))
    comp = comp.div(comp.sum(axis=1), axis=0)
    comp = comp.reindex(sc_len)
    comp = comp[[s for s in STATE_ORDER if s in comp.columns]]
    fig, ax = plt.subplots(figsize=(11, 5))
    bottom = np.zeros(len(comp))
    for s in comp.columns:
        ax.bar(comp.index, comp[s].values, bottom=bottom,
               color=STATE_COLORS.get(s,"#888"), label=s, width=0.85)
        bottom += comp[s].values
    ax.set_ylabel("Fraction of windows")
    ax.set_title("Chromatin-state composition per scaffold (top 25 by window count)")
    ax.set_xticklabels(comp.index, rotation=60, ha="right", fontsize=7)
    ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1.0), loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "chromatin_states_per_scaffold.pdf"))
    plt.close(fig)

    # 3. Heatmap superfamily x state (log2OR), top families
    top_fams = (te_long.groupby("superfam")["bp"].sum()
                .sort_values(ascending=False).head(args.top_fams).index.tolist())
    pivot = (enr[enr["superfam"].isin(top_fams)]
             .pivot(index="superfam", columns="state", values="log2OR")
             .reindex(index=top_fams)
             .reindex(columns=[s for s in STATE_ORDER if s in enr["state"].unique()]))
    fig, ax = plt.subplots(figsize=(7, 0.45*len(top_fams)+1.5))
    vmax = np.nanmax(np.abs(pivot.values))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(pivot.shape[1])); ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0])); ax.set_yticklabels(pivot.index, fontsize=8)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6,
                        color="white" if abs(v) > 0.6*vmax else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("log2 (obs / exp) of TE bp", fontsize=8)
    ax.set_title(f"TE preference per chromatin state (top {len(top_fams)} superfamilies)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "TE_chromatin_heatmap.pdf"))
    plt.close(fig)

    # 4. Stacked composition of TE bp inside each state
    pivot2 = (te_long[te_long["superfam"].isin(top_fams)]
              .pivot(index="state", columns="superfam", values="bp")
              .reindex([s for s in STATE_ORDER if s in te_long["state"].unique()])
              .fillna(0))
    pivot2 = pivot2.div(pivot2.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    bottom = np.zeros(len(pivot2))
    cmap = plt.get_cmap("tab20")
    for i, fa in enumerate(pivot2.columns):
        ax.bar(pivot2.index, pivot2[fa].values, bottom=bottom,
               color=cmap(i % 20), label=fa, width=0.7)
        bottom += pivot2[fa].values
    ax.set_ylabel("Fraction of TE bp")
    ax.set_title("TE composition within each chromatin state")
    ax.set_xticklabels(pivot2.index, rotation=30, ha="right", fontsize=8)
    ax.legend(fontsize=6, bbox_to_anchor=(1.01,1.0), loc="upper left", ncol=1)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "TE_chromatin_stacked.pdf"))
    plt.close(fig)

    print(f"[CHRSTATE] DONE. Outputs at {args.outdir}")

if __name__ == "__main__":
    main()
