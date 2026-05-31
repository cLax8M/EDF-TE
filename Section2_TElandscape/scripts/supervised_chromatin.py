#!/usr/bin/env python3
"""
supervised_chromatin.py
=======================
Aprende un PERFIL DE HETEROCROMATINA a partir de Rhimi59 (DANPOS MNase-seq +
H3K9me3 broad peaks) y lo aplica a TODAS las especies para reasignar estados
de cromatina con UMBRALES APRENDIDOS en lugar de cuantiles per-species.

Entrenamiento (Rhimi59):
  - features por ventana (500 bp): mean_occ, log_mean_occ, cv, ac_val, ac_lag
                                   + log mean_occ vecindario (smoothing local)
  - label: fracción de la ventana que solapa H3K9me3 broad peaks > 0.5
  - modelo: regresión logística (clase positiva = heterocromatina)

Salidas:
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/
    model_coefs.tsv            coeficientes y AUC
    roc.pdf                    curva ROC sobre Rhimi59
    feature_importance.pdf
    learned_thresholds.json    {p_open, p_het} y umbrales feature-wise
    Rhimi59_supervised_states.tsv.gz    (chrom,start,end,p_het,state_sup)
    confusion_supervised_vs_unsupervised.pdf
    cross_species/                       por especie:
       <SP>_supervised.tsv.gz
       <SP>_TE_enrichment.tsv
    composition_supervised.tsv + .pdf
    Het_TE_log2OR_heatmap.pdf            cross-species, estado Heterochromatin
"""
import argparse, os, sys, gzip, json, glob
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import (
    load_nuclpred_wig, window_metrics
)

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc

# ----------------------------- paths -----------------------------
GENOMES   = "/legserv/Temp/Carlos/EDF-TE/data/genomes"
TE_DIR    = "/data/carlos/EDF-TE/TE_analysis/annotated"
NUCL_DIR  = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"

# Rhimi59 ground truth + DANPOS WIG
RHIMI59             = "MU_Rhimi59"
RHIMI59_DANPOS_WIG  = "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"
H3K9_BED            = "/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/Epigenetics/RmH3K9me3/H3K9_WT_MACS2_broad_peaks.bed"

# Real DANPOS MNase-seq for the 5 additional species
EXTRA_MNASE = {
    "AS_Usnflo1":  "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/AS_Usnflo1/AS_Usnflo1_danpos/pooled/AS_Usnflo1_dedup.smooth.wig",
    "BA_Mereb1":   "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Mereb1/BA_Mereb1_danpos/pooled/BA_Mereb1_dedup.smooth.wig",
    "BA_Mrafri1":  "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Mrafri1/BA_Mrafri1_danpos/pooled/BA_Mrafri1_dedup.smooth.wig",
    "BA_Rhoto":    "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Rhoto/BA_Rhoto_danpos/pooled/BA_Rhoto_dedup.smooth.wig",
    "BL_Catan2":   "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BL_Catan2/BL_Catan2_danpos/pooled/BL_Catan2_dedup.smooth.wig",
}

OUT_ROOT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin"

WIN = 500

# ----------------------------- helpers -----------------------------
def chrom_sizes_from_fai(fai_path):
    fai = pd.read_csv(fai_path, sep="\t", header=None,
                      names=["chrom","length","off","lineb","lineby"])
    return dict(zip(fai["chrom"], fai["length"].astype(int)))


def windows_dataframe(arrs, win=WIN):
    rows = []
    for c, occ in arrs.items():
        L = len(occ)
        n = L // win
        if n == 0: continue
        m, cv, ac_lag, ac_v, ok = window_metrics(occ, win)
        starts = np.arange(n, dtype=np.int64) * win
        rows.append(pd.DataFrame({
            "chrom": c, "start": starts, "end": starts + win,
            "mean_occ": m, "cv": cv, "ac_lag": ac_lag, "ac_val": ac_v,
        }))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def add_neighborhood(df, k=5):
    """Mean-occ smoothed within each scaffold + per-species percentile rank."""
    out = []
    for c, sub in df.groupby("chrom", sort=False):
        m = sub["mean_occ"].values.astype(np.float32)
        n = len(m)
        if n == 0:
            sm = m
        else:
            kk = min(k, max(0, (n - 1) // 2))
            if kk == 0:
                sm = m.copy()
            else:
                kernel = np.ones(2*kk+1, dtype=np.float32) / (2*kk+1)
                sm = np.convolve(m, kernel, mode="same")
        sub = sub.copy(); sub["mean_occ_local"] = sm
        out.append(sub)
    df = pd.concat(out, ignore_index=True)
    # Per-SAMPLE percentile-rank normalization (scale-invariant: makes
    # DANPOS, real MNase, and nuclpred z-score comparable)
    def pct_rank(x):
        x = np.asarray(x, dtype=np.float64)
        # rank only over informative bp; non-informative -> 0
        info = x > 1e-3
        out = np.zeros_like(x)
        if info.sum() > 0:
            ranks = pd.Series(x[info]).rank(method="average").values
            out[info] = ranks / info.sum()
        return out.astype(np.float32)
    df["occ_pct"]       = pct_rank(df["mean_occ"].values)
    df["occ_local_pct"] = pct_rank(df["mean_occ_local"].values)
    return df


def label_h3k9(df, bed, win=WIN, frac_thr=0.5):
    peaks = pd.read_csv(bed, sep="\t", header=None, usecols=[0,1,2],
                        names=["chrom","start","end"])
    overlap = np.zeros(len(df), dtype=np.float32)
    by_chr = {c: peaks.loc[peaks["chrom"]==c, ["start","end"]].values
              for c in peaks["chrom"].unique()}
    chrom_arr = df["chrom"].values
    s_arr = df["start"].values; e_arr = df["end"].values
    for c, P in by_chr.items():
        mask = chrom_arr == c
        if not mask.any() or len(P)==0: continue
        ws = s_arr[mask]; we = e_arr[mask]
        ov = np.zeros(len(ws), dtype=np.int64)
        # for each peak, accumulate overlap into windows
        for ps, pe in P:
            o = np.minimum(we, pe) - np.maximum(ws, ps)
            o[o<0] = 0
            ov += o
        overlap[mask] = ov / win
    df = df.copy()
    df["frac_h3k9"] = overlap
    df["label"] = (overlap >= frac_thr).astype(int)
    return df

# ----------------------------- training -----------------------------
FEATS = ["occ_pct", "occ_local_pct", "cv", "ac_val"]


def featurize(df):
    df = df.copy()
    for c in ["occ_pct","occ_local_pct","cv","ac_val","ac_lag","mean_occ","mean_occ_local"]:
        if c in df.columns:
            df[c] = np.nan_to_num(df[c].values, nan=0.0, posinf=0.0, neginf=0.0)
    return df


def train_model(df_train):
    X = df_train[FEATS].values.astype(np.float32)
    y = df_train["label"].values.astype(int)
    scaler = StandardScaler().fit(X)
    Xs = scaler.transform(X)
    clf = LogisticRegression(max_iter=2000, class_weight="balanced",
                             solver="lbfgs", C=1.0)
    clf.fit(Xs, y)
    p = clf.predict_proba(Xs)[:,1]
    fpr, tpr, thr = roc_curve(y, p)
    auc_v = auc(fpr, tpr)
    # Youden cut for "Heterochromatin"
    j = tpr - fpr
    idx_het  = int(np.argmax(j))
    p_het    = float(thr[idx_het])
    # Define "Open" as low p < specificity-driven cutoff (1st quartile of negatives)
    p_open   = float(np.quantile(p[y==0], 0.25))
    return clf, scaler, dict(auc=float(auc_v), fpr=fpr.tolist(),
                             tpr=tpr.tolist(), thr=thr.tolist(),
                             p_het=p_het, p_open=p_open)


def predict_states(df, clf, scaler, p_open, p_het):
    df = df.copy()
    X = df[FEATS].values.astype(np.float32)
    Xs = scaler.transform(X)
    p = clf.predict_proba(Xs)[:,1]
    df["p_het"] = p
    state = np.full(len(df), "Intermediate", dtype=object)
    state[(df["mean_occ"].values <= 1e-3)] = "No_signal"
    state[(p < p_open) & (df["mean_occ"].values > 1e-3)] = "Open"
    state[p >= p_het] = "Heterochromatin"
    df["state_sup"] = state
    return df

STATE_ORDER  = ["No_signal","Open","Intermediate","Heterochromatin"]
STATE_COLORS = {"No_signal":"#CCCCCC","Open":"#6BAED6",
                "Intermediate":"#FDAE6B","Heterochromatin":"#E6550D"}

# ----------------------------- TE enrichment -----------------------------
def te_enrichment(df_states, te_bed, win=WIN):
    te = pd.read_csv(te_bed, sep="\t", header=None,
                     names=["chrom","start","end","superfam","score","strand"],
                     usecols=[0,1,2,3])
    chrom_starts = {c: df_states.loc[df_states["chrom"]==c,"start"].values
                    for c in df_states["chrom"].unique()}
    chrom_states = {c: df_states.loc[df_states["chrom"]==c,"state_sup"].values
                    for c in df_states["chrom"].unique()}
    bp_mat = {}
    for row in te.itertuples(index=False):
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
    state_bp = {st: int((df_states["state_sup"]==st).sum())*win for st in STATE_ORDER}
    fam_bp = {fa: int(np.minimum(te["end"]-te["start"],
                                 te["end"]-te["start"]).where(te["superfam"]==fa,0).sum())
              for fa in te["superfam"].unique()}
    total_bp = sum(state_bp.values())
    rows = []
    for (st, fa), obs in bp_mat.items():
        sbp = state_bp.get(st,0); fbp = fam_bp.get(fa,0)
        if sbp == 0 or fbp == 0 or total_bp == 0:
            log2or = np.nan
        else:
            log2or = np.log2(((obs+1)/(sbp+1)) / ((fbp+1)/(total_bp+1)))
        rows.append(dict(state=st, superfam=fa, obs_bp=obs,
                         state_bp=sbp, fam_bp=fbp, log2OR=log2or))
    return pd.DataFrame(rows)

# ----------------------------- per species worker -----------------------------
def process_species(args_tuple):
    sp, wig, fai_path, te_path, clf, scaler, p_open, p_het, outdir = args_tuple
    try:
        cs = chrom_sizes_from_fai(fai_path)
        arrs = load_nuclpred_wig(wig, cs)
        df = windows_dataframe(arrs, WIN)
        if len(df) == 0:
            return sp, False, "empty"
        df = add_neighborhood(df, k=5)
        df = featurize(df)
        df = predict_states(df, clf, scaler, p_open, p_het)
        out_tsv = os.path.join(outdir, f"{sp}_supervised.tsv.gz")
        df[["chrom","start","end","mean_occ","ac_val","p_het","state_sup"]].to_csv(
            out_tsv, sep="\t", index=False, compression="gzip")
        # composition
        comp = df["state_sup"].value_counts(normalize=True).reindex(STATE_ORDER, fill_value=0)
        comp.to_csv(os.path.join(outdir, f"{sp}_composition.tsv"), sep="\t")
        # TE enrichment
        if os.path.exists(te_path):
            enr = te_enrichment(df, te_path, WIN)
            enr.to_csv(os.path.join(outdir, f"{sp}_TE_enrichment.tsv"),
                       sep="\t", index=False)
        return sp, True, "OK"
    except Exception as e:
        return sp, False, str(e)

# ----------------------------- main pipeline -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=12)
    ap.add_argument("--skip_pan", action="store_true",
                    help="solo entrenamiento + 5 MNase, sin las 215 nuclpred")
    args = ap.parse_args()

    os.makedirs(OUT_ROOT, exist_ok=True)
    cs_dir = os.path.join(OUT_ROOT, "cross_species")
    os.makedirs(cs_dir, exist_ok=True)

    # ============ TRAIN on Rhimi59 ============
    print("[TRAIN] loading Rhimi59 DANPOS WIG ...", flush=True)
    cs_rhi = chrom_sizes_from_fai(os.path.join(GENOMES, f"{RHIMI59}.fasta.fai"))
    arrs   = load_nuclpred_wig(RHIMI59_DANPOS_WIG, cs_rhi)
    df_rhi = windows_dataframe(arrs, WIN)
    df_rhi = add_neighborhood(df_rhi, k=5)
    df_rhi = featurize(df_rhi)
    df_rhi = label_h3k9(df_rhi, H3K9_BED, WIN, frac_thr=0.5)
    print(f"[TRAIN] {len(df_rhi)} windows; {df_rhi['label'].sum()} positives ({df_rhi['label'].mean()*100:.1f}%)",
          flush=True)
    # use only windows with at least minimal signal for fitting
    train_mask = df_rhi["mean_occ"] > 1e-3
    df_train = df_rhi[train_mask].reset_index(drop=True)
    clf, scaler, info = train_model(df_train)
    p_open = info["p_open"]; p_het = info["p_het"]
    print(f"[TRAIN] AUC={info['auc']:.3f}  p_open={p_open:.3f}  p_het={p_het:.3f}",
          flush=True)
    # save model coefficients
    coefs = pd.DataFrame({"feature": FEATS,
                          "coef": clf.coef_[0].tolist(),
                          "scaler_mean": scaler.mean_.tolist(),
                          "scaler_scale": scaler.scale_.tolist()})
    coefs.to_csv(os.path.join(OUT_ROOT, "model_coefs.tsv"), sep="\t", index=False)
    coefs.loc[len(coefs)] = ["intercept", float(clf.intercept_[0]), 0.0, 1.0]
    coefs.to_csv(os.path.join(OUT_ROOT, "model_coefs.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT_ROOT, "learned_thresholds.json"), "w") as fh:
        json.dump({"auc": info["auc"], "p_open": p_open, "p_het": p_het,
                   "training_windows": int(len(df_train)),
                   "training_positives": int(df_train['label'].sum()),
                   "feature_means": dict(zip(FEATS, scaler.mean_.tolist())),
                   "feature_scales": dict(zip(FEATS, scaler.scale_.tolist())),
                   "feature_coefs":  dict(zip(FEATS, clf.coef_[0].tolist())),
                   "intercept":      float(clf.intercept_[0])},
                  fh, indent=2)

    # ROC plot
    fpr = np.array(info["fpr"]); tpr = np.array(info["tpr"])
    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(fpr, tpr, color="#E6550D", lw=2,
            label=f"AUC={info['auc']:.3f}")
    ax.plot([0,1],[0,1], "k--", lw=0.8)
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title(f"H3K9me3 classifier (Rhimi59)\nN={len(df_train):,} windows")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(OUT_ROOT,"roc.pdf")); plt.close(fig)

    # Feature importance plot
    fig, ax = plt.subplots(figsize=(5,3))
    ax.barh(FEATS, clf.coef_[0], color=["#6BAED6" if c<0 else "#E6550D" for c in clf.coef_[0]])
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("Logistic coefficient (standardized)")
    ax.set_title("Heterochromatin classifier weights")
    fig.tight_layout(); fig.savefig(os.path.join(OUT_ROOT,"feature_importance.pdf")); plt.close(fig)

    # Apply to Rhimi59 → save + confusion vs unsupervised
    df_rhi_pred = predict_states(df_rhi, clf, scaler, p_open, p_het)
    df_rhi_pred[["chrom","start","end","mean_occ","ac_val",
                 "p_het","state_sup","frac_h3k9","label"]].to_csv(
        os.path.join(OUT_ROOT, "Rhimi59_supervised_states.tsv.gz"),
        sep="\t", index=False, compression="gzip")

    # H3K9me3 fraction per supervised state (validation)
    val = df_rhi_pred.groupby("state_sup").agg(
        n=("label","size"),
        n_pos=("label","sum"),
        mean_frac_h3k9=("frac_h3k9","mean")).reindex(STATE_ORDER)
    val["pct_h3k9"] = val["n_pos"]/val["n"]*100
    val.to_csv(os.path.join(OUT_ROOT,"Rhimi59_state_validation.tsv"), sep="\t")
    print(val, flush=True)

    fig, ax = plt.subplots(figsize=(5,3))
    bars = ax.bar(val.index, val["pct_h3k9"],
                  color=[STATE_COLORS[s] for s in val.index])
    ax.set_ylabel("% windows H3K9me3+ (frac>0.5)")
    ax.set_title("Supervised states vs H3K9me3 (Rhimi59)")
    for b,v in zip(bars, val["pct_h3k9"]):
        ax.text(b.get_x()+b.get_width()/2, v+0.3, f"{v:.1f}",
                ha="center", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT_ROOT,"Rhimi59_state_validation.pdf"))
    plt.close(fig)

    # ============ APPLY to 5 real-MNase fungi ============
    print("[MNASE5] applying to 5 real MNase species ...", flush=True)
    mnase_jobs = []
    for sp, w in EXTRA_MNASE.items():
        fai = os.path.join(GENOMES, f"{sp}.fasta.fai")
        te  = os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles", f"{sp}.filteredRepeats.bed")
        if not (os.path.exists(fai) and os.path.exists(w)):
            print(f"  [MNASE5] skip {sp}: missing files"); continue
        mnase_jobs.append((sp, w, fai, te, clf, scaler, p_open, p_het, cs_dir))
    # also add Rhimi59 itself for cross-species composition
    mnase_jobs.append((RHIMI59, RHIMI59_DANPOS_WIG,
                       os.path.join(GENOMES, f"{RHIMI59}.fasta.fai"),
                       os.path.join(TE_DIR, f"{RHIMI59}_EarlGrey",
                                    f"{RHIMI59}_summaryFiles",
                                    f"{RHIMI59}.filteredRepeats.bed"),
                       clf, scaler, p_open, p_het, cs_dir))
    done = []
    with ProcessPoolExecutor(max_workers=min(args.nproc, len(mnase_jobs))) as ex:
        futs = [ex.submit(process_species, j) for j in mnase_jobs]
        for fu in as_completed(futs):
            sp, ok, msg = fu.result()
            print(f"  [MNASE5] {sp} -> {msg}", flush=True)
            if ok: done.append(sp)

    # ============ APPLY to all 215 nuclpred species ============
    if not args.skip_pan:
        print("[PAN] applying classifier to all nuclpred species ...", flush=True)
        wigs = sorted(glob.glob(os.path.join(NUCL_DIR, "*_nuclpred.wig")))
        pan_jobs = []
        for w in wigs:
            sp = os.path.basename(w).replace("_nuclpred.wig","")
            if sp == RHIMI59 or sp in EXTRA_MNASE: continue  # MNase has priority
            fai = os.path.join(GENOMES, f"{sp}.fasta.fai")
            te  = os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles", f"{sp}.filteredRepeats.bed")
            if not (os.path.exists(fai) and os.path.exists(te)): continue
            pan_jobs.append((sp, w, fai, te, clf, scaler, p_open, p_het, cs_dir))
        print(f"[PAN] {len(pan_jobs)} species queued", flush=True)
        with ProcessPoolExecutor(max_workers=args.nproc) as ex:
            futs = [ex.submit(process_species, j) for j in pan_jobs]
            for i, fu in enumerate(as_completed(futs), 1):
                sp, ok, msg = fu.result()
                if i % 10 == 0 or not ok:
                    print(f"  [PAN {i}/{len(pan_jobs)}] {sp} -> {msg}", flush=True)
                if ok: done.append(sp)

    # ============ AGGREGATE composition + Het TE log2OR ============
    comp_files = sorted(glob.glob(os.path.join(cs_dir, "*_composition.tsv")))
    rows = []
    for f in comp_files:
        sp = os.path.basename(f).replace("_composition.tsv","")
        s = pd.read_csv(f, sep="\t", index_col=0).iloc[:,0]
        rows.append([sp] + [float(s.get(st,0.0)) for st in STATE_ORDER])
    comp = pd.DataFrame(rows, columns=["species"]+STATE_ORDER)
    comp["group"] = comp["species"].str.split("_").str[0]
    comp.sort_values(["group","species"], inplace=True)
    comp.to_csv(os.path.join(OUT_ROOT,"composition_supervised.tsv"),
                sep="\t", index=False)

    # stacked composition plot
    fig, ax = plt.subplots(figsize=(max(8, 0.06*len(comp)), 4))
    bottom = np.zeros(len(comp))
    for st in STATE_ORDER:
        ax.bar(comp["species"], comp[st], bottom=bottom,
               color=STATE_COLORS[st], width=0.95, label=st)
        bottom += comp[st].values
    ax.set_xticks([]); ax.set_ylabel("Genome fraction")
    ax.set_title(f"Supervised chromatin composition (n={len(comp)} species)")
    ax.legend(loc="upper right", fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT_ROOT,"composition_supervised.pdf"))
    plt.close(fig)

    # cross-species Heterochromatin TE log2OR heatmap
    enr_files = sorted(glob.glob(os.path.join(cs_dir, "*_TE_enrichment.tsv")))
    hets = []
    for f in enr_files:
        sp = os.path.basename(f).replace("_TE_enrichment.tsv","")
        e = pd.read_csv(f, sep="\t")
        e = e[e["state"]=="Heterochromatin"][["superfam","log2OR","obs_bp"]]
        e = e.rename(columns={"log2OR": sp})
        hets.append(e[["superfam", sp]])
    if hets:
        from functools import reduce
        het_mat = reduce(lambda a,b: a.merge(b, on="superfam", how="outer"), hets)
        het_mat = het_mat.set_index("superfam")
        # focus on TE superfamilies seen in many species with reasonable signal
        good = het_mat.notna().sum(axis=1) >= max(10, len(het_mat.columns)//4)
        het_mat = het_mat.loc[good]
        # rank superfams by mean
        het_mat["__mean"] = het_mat.mean(axis=1, skipna=True)
        het_mat = het_mat.sort_values("__mean", ascending=False).drop(columns="__mean")
        # top 25 superfams
        het_mat = het_mat.head(25)
        het_mat.to_csv(os.path.join(OUT_ROOT,"Het_TE_log2OR_matrix.tsv"), sep="\t")

        # heatmap
        cols_order = comp["species"].tolist()
        het_plot = het_mat[[c for c in cols_order if c in het_mat.columns]]
        fig, ax = plt.subplots(figsize=(max(10, 0.06*het_plot.shape[1]),
                                        0.35*het_plot.shape[0]+1))
        m = het_plot.values
        vmax = np.nanpercentile(np.abs(m), 95) or 1.0
        im = ax.imshow(m, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       interpolation="nearest")
        ax.set_yticks(range(het_plot.shape[0])); ax.set_yticklabels(het_plot.index, fontsize=7)
        ax.set_xticks([]); ax.set_xlabel(f"Species (n={het_plot.shape[1]}, ordered by group)")
        ax.set_title("TE log2OR in Heterochromatin (supervised)")
        cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
        cb.set_label("log2OR")
        fig.tight_layout(); fig.savefig(os.path.join(OUT_ROOT,"Het_TE_log2OR_heatmap.pdf"))
        plt.close(fig)

    print(f"[DONE] outputs at {OUT_ROOT}")

if __name__ == "__main__":
    main()
