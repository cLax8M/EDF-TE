#!/usr/bin/env python3
"""
simple_2feature_compaction.py
=============================
Modelo simple e interpretable: compactación = OCUPANCY ALTA + CV BAJO.

Aprendido en Rhimi59 (DANPOS MNase + H3K9me3 broad peaks):
  - per-window features: mean_occ (rank percentil per-species)
                         cv (rank percentil per-species)
    Rank percentil hace la regla transferible entre escalas (DANPOS, MNase
    real, nuclpred z-score).
  - score = occ_pct - cv_pct       (alto -> compactado)
  - umbral aprendido: el que maximiza el Youden's J (TPR-FPR) frente a
    H3K9me3 en Rhimi59. Ventanas con occ ~ 0 -> NoSig (excluidas).
  - 2 estados: Open / Compacted

Outputs en /.../supervised_chromatin/simple2f/:
  learned_threshold.json
  rhimi59_score_vs_h3k9.pdf      (ROC + densidad del score por TE×H3K9)
  composition_2state.tsv / .pdf
  composition_2state_per_group.tsv / .pdf
  cross_species_2state/<SP>_states.tsv.gz                (chrom,start,end,score,state)
  cross_species_2state/<SP>_family_enrichment.tsv
  cross_species_2state/<SP>_superfam_enrichment.tsv
  Compacted_superfam_log2OR_heatmap.pdf
  Compacted_top_families_per_group.tsv / .pdf
"""
import os, sys, glob, gzip, json, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.metrics import roc_curve, auc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig, window_metrics

# ---------------------------- paths ----------------------------
GENOMES = "/legserv/Temp/Carlos/EDF-TE/data/genomes"
TE_DIR  = "/data/carlos/EDF-TE/TE_analysis/annotated"
NUCL_DIR= "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"

RHIMI59             = "MU_Rhimi59"
RHIMI59_DANPOS_WIG  = "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"
H3K9_BED            = "/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/Epigenetics/RmH3K9me3/H3K9_WT_MACS2_broad_peaks.bed"

EXTRA_MNASE = {
    "AS_Usnflo1":  "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/AS_Usnflo1/AS_Usnflo1_danpos/pooled/AS_Usnflo1_dedup.smooth.wig",
    "BA_Mereb1":   "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Mereb1/BA_Mereb1_danpos/pooled/BA_Mereb1_dedup.smooth.wig",
    "BA_Mrafri1":  "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Mrafri1/BA_Mrafri1_danpos/pooled/BA_Mrafri1_dedup.smooth.wig",
    "BA_Rhoto":    "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BA_Rhoto/BA_Rhoto_danpos/pooled/BA_Rhoto_dedup.smooth.wig",
    "BL_Catan2":   "/data/carlos/EDF-TE/TE_analysis/nucleosomes/nuc-seq/BL_Catan2/BL_Catan2_danpos/pooled/BL_Catan2_dedup.smooth.wig",
}

OUT_ROOT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/simple2f"
CS_DIR   = os.path.join(OUT_ROOT, "cross_species_2state")
os.makedirs(CS_DIR, exist_ok=True)

WIN = 500

NAME_RE = re.compile(r"Name=([^;]+)")

# ------------------------- helpers ----------------------------
def chrom_sizes_from_fai(fai):
    fai = pd.read_csv(fai, sep="\t", header=None,
                      names=["chrom","L","off","lb","lby"])
    return dict(zip(fai["chrom"], fai["L"].astype(int)))


def windows_metrics_df(arrs, win=WIN):
    rows = []
    for c, occ in arrs.items():
        n = len(occ) // win
        if n == 0: continue
        means, cvs, _aclag, _acv, _v = window_metrics(occ, win)
        rows.append(pd.DataFrame({
            "chrom": c,
            "start": np.arange(n) * win,
            "end":   np.arange(n) * win + win,
            "mean_occ": means,
            "cv": cvs,
        }))
    return pd.concat(rows, ignore_index=True)


def add_score(df):
    """Add per-species rank-percentile of occupancy and cv, plus score.

    score = occ_pct - cv_pct     (high => high occ AND low cv)
    Only informative windows (mean_occ > 1e-3) get a rank; others -> NaN."""
    df = df.copy()
    occ = df["mean_occ"].values
    cv  = df["cv"].values
    info = occ > 1e-3
    occ_pct = np.full(len(df), np.nan, dtype=np.float32)
    cv_pct  = np.full(len(df), np.nan, dtype=np.float32)
    if info.sum() > 0:
        occ_pct[info] = pd.Series(occ[info]).rank(method="average").values / info.sum()
        cv_pct[info]  = pd.Series(cv[info]).rank(method="average").values / info.sum()
    df["occ_pct"] = occ_pct
    df["cv_pct"]  = cv_pct
    df["score"]   = occ_pct - cv_pct
    df["informative"] = info
    return df


def overlap_fraction_per_window(bed_path, sizes, win=WIN):
    """Per-window fraction overlap with a BED."""
    bp = {c: np.zeros(L, dtype=np.uint8) for c, L in sizes.items()}
    if not os.path.exists(bed_path):
        return None
    with open(bed_path) as fh:
        for ln in fh:
            f = ln.split("\t")
            if len(f) < 3: continue
            c = f[0]
            if c not in bp: continue
            try:
                s = int(f[1]); e = int(f[2])
            except ValueError: continue
            e = min(e, len(bp[c])); s = max(s, 0)
            if s < e: bp[c][s:e] = 1
    out = {}
    for c, L in sizes.items():
        n = L // win
        if n == 0: out[c] = np.zeros(0, dtype=np.float32); continue
        a = bp[c][:n*win].reshape(n, win)
        out[c] = a.mean(axis=1).astype(np.float32)
    return out


def parse_gff_families(gff_path):
    rows = []
    if not os.path.exists(gff_path):
        return pd.DataFrame(columns=["chrom","start","end","superfam","family"])
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip(): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9: continue
            try:
                s = int(f[3]) - 1; e = int(f[4])
            except ValueError: continue
            m = NAME_RE.search(f[8])
            fam = m.group(1) if m else "NA"
            rows.append((f[0], s, e, f[2], fam))
    return pd.DataFrame(rows, columns=["chrom","start","end","superfam","family"])


# ------------------------- training ----------------------------
def train_threshold():
    """Learn score threshold on Rhimi59 (DANPOS + H3K9me3)."""
    print("[TRAIN] Rhimi59 windows", flush=True)
    sizes = chrom_sizes_from_fai(os.path.join(GENOMES, f"{RHIMI59}.fasta.fai"))
    arrs  = load_nuclpred_wig(RHIMI59_DANPOS_WIG, sizes)
    df = windows_metrics_df(arrs, WIN)
    df = add_score(df)
    h3 = overlap_fraction_per_window(H3K9_BED, sizes, WIN)
    flag = []
    for c, sub in df.groupby("chrom", sort=False):
        flag.append(h3[c][:len(sub)])
    df["frac_h3k9"] = np.concatenate(flag)
    df["label"]     = (df["frac_h3k9"] >= 0.5).astype(int)

    info = df[df["informative"]].copy()
    fpr, tpr, thr = roc_curve(info["label"].values, info["score"].values)
    roc_auc = auc(fpr, tpr)
    j = tpr - fpr
    idx = int(np.argmax(j))
    s_thr = float(thr[idx])
    print(f"[TRAIN] AUC={roc_auc:.3f}  score_thr={s_thr:.3f}", flush=True)

    out = {
        "rule": "score = occ_pct - cv_pct ; Compacted iff informative AND score >= score_thr",
        "score_threshold": s_thr,
        "auc": roc_auc,
        "n_train_windows": int(len(info)),
        "n_train_positives": int(info["label"].sum()),
        "training_youden_J": float(j[idx]),
    }
    with open(os.path.join(OUT_ROOT, "learned_threshold.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    # ROC + score density plot
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    ax.plot(fpr, tpr, color="#D62728", lw=2, label=f"AUC = {roc_auc:.3f}")
    ax.plot([0,1],[0,1], color="grey", lw=0.8, ls="--")
    ax.scatter([fpr[idx]], [tpr[idx]], color="black", zorder=5,
               label=f"thr = {s_thr:.3f} (Youden J = {j[idx]:.3f})")
    ax.set_xlabel("FPR (H3K9me3-)"); ax.set_ylabel("TPR (H3K9me3+)")
    ax.set_title("ROC – score = occ_pct − cv_pct  (Rhimi59)"); ax.legend(loc="lower right", fontsize=9)

    ax = axes[1]
    for lbl, color, name in [(0, "#7F7F7F", "H3K9me3-"),
                              (1, "#D62728", "H3K9me3+")]:
        d = info.loc[info["label"]==lbl, "score"].values
        ax.hist(d, bins=50, alpha=0.55, color=color, density=True,
                label=f"{name} (n={len(d)})")
    ax.axvline(s_thr, color="black", lw=1.2, ls="--", label=f"thr = {s_thr:.3f}")
    ax.set_xlabel("score = occ_pct − cv_pct"); ax.set_ylabel("density")
    ax.legend(fontsize=9); ax.set_title("Score distribution by H3K9me3 label")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "rhimi59_score_vs_h3k9.pdf"))
    plt.close(fig)

    # validation table
    info["state"] = np.where(info["score"] >= s_thr, "Compacted", "Open")
    val = info.groupby("state").agg(
        n=("label","size"), n_pos=("label","sum"),
        mean_frac_h3k9=("frac_h3k9","mean"))
    val["pct_h3k9"] = 100 * val["n_pos"] / val["n"]
    val.to_csv(os.path.join(OUT_ROOT, "rhimi59_validation.tsv"), sep="\t")
    print(val, flush=True)
    return s_thr


# ------------------------- per-species worker ----------------------------
def process_species(args):
    sp, wig, fai, gff, bed, score_thr = args
    try:
        sizes = chrom_sizes_from_fai(fai)
        if not os.path.exists(wig):
            return sp, None, None, None, "wig missing"
        arrs = load_nuclpred_wig(wig, sizes)
        df = windows_metrics_df(arrs, WIN)
        df = add_score(df)
        df["state"] = np.where(df["informative"] & (df["score"] >= score_thr),
                               "Compacted",
                               np.where(df["informative"], "Open", "NoSig"))
        # save
        df_out = df[["chrom","start","end","mean_occ","cv","score","state"]].copy()
        df_out.to_csv(os.path.join(CS_DIR, f"{sp}_states.tsv.gz"),
                      sep="\t", index=False, compression="gzip")
        info = df[df["informative"]]
        if len(info) == 0:
            return sp, None, None, None, "no informative"
        comp = info["state"].value_counts(normalize=True).reindex(
            ["Open","Compacted"], fill_value=0.0)

        # state per window per chrom
        sw = {c: sub["state"].values for c, sub in df.groupby("chrom", sort=False)}

        # per-family + per-superfam enrichment
        te = parse_gff_families(gff)
        if te.empty:
            return sp, comp, None, None, None

        bp = {}
        n_copies = {}
        for r in te.itertuples(index=False):
            n_copies[(r.family, r.superfam)] = n_copies.get((r.family, r.superfam), 0) + 1
            arr = sw.get(r.chrom)
            if arr is None: continue
            i0 = r.start // WIN
            i1 = (r.end - 1) // WIN
            i1 = min(i1, len(arr) - 1)
            for i in range(i0, i1 + 1):
                ws = i * WIN; we = ws + WIN
                ov = min(r.end, we) - max(r.start, ws)
                if ov <= 0: continue
                st = arr[i]
                if st == "NoSig": continue
                bp[(r.family, r.superfam, st)] = bp.get((r.family, r.superfam, st), 0) + ov

        win_state_bp = {st: int((info["state"]==st).sum()) * WIN
                        for st in ["Open","Compacted"]}
        tot_bp = sum(win_state_bp.values())
        if tot_bp == 0:
            return sp, comp, None, None, None
        p_comp_bg = win_state_bp["Compacted"] / tot_bp

        # family
        fam_rows = []
        fams = sorted({(f, s) for (f, s, _) in bp.keys()})
        for fam, sf in fams:
            t = sum(bp.get((fam, sf, st), 0) for st in ["Open","Compacted"])
            if t < 200: continue
            obs = bp.get((fam, sf, "Compacted"), 0)
            exp = t * p_comp_bg
            log2OR = float(np.log2((obs + 1) / (exp + 1)))
            fam_rows.append(dict(family=fam, superfam=sf,
                                 n_copies=n_copies.get((fam, sf), 0),
                                 bp_total=t, bp_compacted=obs,
                                 frac_compacted=obs/t,
                                 log2OR_compacted=log2OR))
        fam_df = pd.DataFrame(fam_rows).sort_values("log2OR_compacted", ascending=False)

        sf_bp = {}
        for (fam, sf, st), v in bp.items():
            sf_bp.setdefault(sf, {"Open":0, "Compacted":0})
            if st in sf_bp[sf]:
                sf_bp[sf][st] += v
        sf_rows = []
        for sf, dd in sf_bp.items():
            t = sum(dd.values())
            if t < 500: continue
            obs = dd.get("Compacted", 0)
            exp = t * p_comp_bg
            sf_rows.append(dict(superfam=sf, bp_total=t, bp_compacted=obs,
                                frac_compacted=obs/t,
                                log2OR_compacted=float(np.log2((obs+1)/(exp+1)))))
        sf_df = pd.DataFrame(sf_rows).sort_values("log2OR_compacted", ascending=False)

        if len(fam_df):
            fam_df.to_csv(os.path.join(CS_DIR, f"{sp}_family_enrichment.tsv"),
                          sep="\t", index=False)
        if len(sf_df):
            sf_df.to_csv(os.path.join(CS_DIR, f"{sp}_superfam_enrichment.tsv"),
                         sep="\t", index=False)
        return sp, comp, fam_df, sf_df, None
    except Exception as ex:
        return sp, None, None, None, str(ex)


# --------------------------------- main ---------------------------------
def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    score_thr = train_threshold()

    # Build job list
    jobs = []
    # Rhimi59 + 5 MNase-real
    jobs.append((RHIMI59,
                 RHIMI59_DANPOS_WIG,
                 os.path.join(GENOMES, f"{RHIMI59}.fasta.fai"),
                 os.path.join(TE_DIR, f"{RHIMI59}_EarlGrey", f"{RHIMI59}_summaryFiles",
                              f"{RHIMI59}.filteredRepeats.gff"),
                 os.path.join(TE_DIR, f"{RHIMI59}_EarlGrey", f"{RHIMI59}_summaryFiles",
                              f"{RHIMI59}.filteredRepeats.bed"),
                 score_thr))
    for sp, wig in EXTRA_MNASE.items():
        jobs.append((sp, wig,
                     os.path.join(GENOMES, f"{sp}.fasta.fai"),
                     os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                                  f"{sp}.filteredRepeats.gff"),
                     os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                                  f"{sp}.filteredRepeats.bed"),
                     score_thr))
    mnase_sp = {RHIMI59} | set(EXTRA_MNASE)
    for fai in sorted(glob.glob(os.path.join(GENOMES, "*.fasta.fai"))):
        sp = os.path.basename(fai).replace(".fasta.fai","")
        if sp in mnase_sp: continue
        wig = os.path.join(NUCL_DIR, f"{sp}_nuclpred.wig")
        if not os.path.exists(wig): continue
        jobs.append((sp, wig, fai,
                     os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                                  f"{sp}.filteredRepeats.gff"),
                     os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                                  f"{sp}.filteredRepeats.bed"),
                     score_thr))
    print(f"[RUN] {len(jobs)} species", flush=True)

    comps = {}; sf_all = {}
    n_done = 0
    with ProcessPoolExecutor(max_workers=16) as ex:
        futs = [ex.submit(process_species, j) for j in jobs]
        for fut in as_completed(futs):
            sp, comp, fam_df, sf_df, err = fut.result()
            n_done += 1
            if err:
                print(f"  [{sp}] ERROR: {err}", flush=True)
            if comp is not None:
                comps[sp] = comp
            if sf_df is not None and len(sf_df):
                for r in sf_df.itertuples(index=False):
                    sf_all.setdefault(r.superfam, {})[sp] = r.log2OR_compacted
            if n_done % 25 == 0:
                print(f"  [{n_done}/{len(jobs)}]", flush=True)

    # ---------------- composition ----------------
    comp_df = pd.DataFrame(comps).T.fillna(0.0)
    comp_df.index.name = "species"
    comp_df["group"] = [s.split("_",1)[0] for s in comp_df.index]
    comp_df = comp_df.sort_values(["group","Compacted"], ascending=[True, False])
    comp_df.to_csv(os.path.join(OUT_ROOT, "composition_2state.tsv"), sep="\t")
    print(f"[COMP] {len(comp_df)} species", flush=True)

    fig, ax = plt.subplots(figsize=(max(8, len(comp_df)*0.04), 4))
    bottom = np.zeros(len(comp_df))
    cmap = {"Open":"#6BAED6", "Compacted":"#E6550D"}
    for st in ["Open","Compacted"]:
        ax.bar(np.arange(len(comp_df)), comp_df[st].values,
               bottom=bottom, color=cmap[st], width=1.0,
               edgecolor="none", label=st)
        bottom += comp_df[st].values
    g = comp_df["group"].values; last = g[0]
    for i, gi in enumerate(g):
        if gi != last:
            ax.axvline(i-0.5, color="black", lw=0.4); last = gi
    ax.set_xlim(-0.5, len(comp_df)-0.5); ax.set_ylim(0,1)
    ax.set_xticks([]); ax.set_ylabel("Fraction (Open + Compacted = 1)")
    ax.set_title(f"Two-feature 2-state composition  (score_thr = {score_thr:.3f})")
    ax.legend(loc="upper right", framealpha=0.9, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "composition_2state.pdf"))
    plt.close(fig)

    grp = comp_df.groupby("group")[["Open","Compacted"]].mean().sort_values(
        "Compacted", ascending=False)
    grp.to_csv(os.path.join(OUT_ROOT, "composition_2state_per_group.tsv"), sep="\t")
    fig, ax = plt.subplots(figsize=(7, 0.32*len(grp)+1))
    bottom = np.zeros(len(grp)); yy = np.arange(len(grp))
    for st in ["Open","Compacted"]:
        ax.barh(yy, grp[st].values, left=bottom, color=cmap[st], label=st)
        bottom += grp[st].values
    ax.set_yticks(yy); ax.set_yticklabels(grp.index); ax.invert_yaxis()
    ax.set_xlim(0,1); ax.set_xlabel("Mean fraction (informative bp)")
    ax.set_title("Mean composition per taxonomic group")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "composition_2state_per_group.pdf"))
    plt.close(fig)

    # ---------------- superfam heatmap ----------------
    sf_mat = pd.DataFrame(sf_all).T.reindex(columns=comp_df.index)
    keep = sf_mat.notna().sum(axis=1) >= max(20, sf_mat.shape[1] // 5)
    sf_mat = sf_mat.loc[keep]
    sf_mat["__mean"] = sf_mat.mean(axis=1, skipna=True)
    sf_mat = sf_mat.sort_values("__mean", ascending=False).drop(columns="__mean")
    sf_mat.to_csv(os.path.join(OUT_ROOT, "Compacted_superfam_log2OR_matrix.tsv"), sep="\t")

    fig, ax = plt.subplots(figsize=(max(10, sf_mat.shape[1]*0.05),
                                    max(4, sf_mat.shape[0]*0.22)))
    vmax = float(np.nanpercentile(np.abs(sf_mat.values), 95)) or 1.0
    im = ax.imshow(sf_mat.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.set_yticks(np.arange(len(sf_mat))); ax.set_yticklabels(sf_mat.index, fontsize=7)
    ax.set_xticks([])
    grp_arr = np.array([s.split("_",1)[0] for s in sf_mat.columns])
    last = grp_arr[0]; xticks_pos=[]; xticks_lab=[]; start=0
    for i, gi in enumerate(np.append(grp_arr, ["__END__"])):
        if gi != last:
            ax.axvline(i-0.5, color="black", lw=0.4)
            xticks_pos.append((start+i-1)/2); xticks_lab.append(last)
            last = gi; start = i
    ax.set_xticks(xticks_pos); ax.set_xticklabels(xticks_lab, fontsize=8)
    ax.set_xlabel("Species (grouped)")
    ax.set_title(f"log2 OR of TE superfamily in Compacted state (score_thr = {score_thr:.3f})")
    plt.colorbar(im, ax=ax, shrink=0.6, label="log2 OR")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "Compacted_superfam_log2OR_heatmap.pdf"))
    plt.close(fig)

    # ---------------- top families per group ----------------
    fam_files = sorted(glob.glob(os.path.join(CS_DIR, "*_family_enrichment.tsv")))
    rows=[]
    for f in fam_files:
        sp = os.path.basename(f).replace("_family_enrichment.tsv","")
        gi = sp.split("_",1)[0]
        try: d = pd.read_csv(f, sep="\t")
        except Exception: continue
        if d.empty: continue
        d = d[(d["bp_total"]>=1000) & (d["log2OR_compacted"]>0.5)]
        d["species"]=sp; d["group"]=gi
        rows.append(d)
    if rows:
        fam_all = pd.concat(rows, ignore_index=True)
        fam_all.to_csv(os.path.join(OUT_ROOT, "Compacted_families_all.tsv"),
                       sep="\t", index=False)
        summary = (fam_all.groupby(["group","superfam"])
                          .agg(n_fam=("family","nunique"),
                               n_sp =("species","nunique"),
                               mean_log2OR=("log2OR_compacted","mean"),
                               median_frac_compacted=("frac_compacted","median"),
                               total_bp=("bp_total","sum"))
                          .reset_index()
                          .sort_values(["group","mean_log2OR"], ascending=[True, False]))
        summary.to_csv(os.path.join(OUT_ROOT, "Compacted_top_families_per_group.tsv"),
                       sep="\t", index=False)
        groups = sorted(summary["group"].unique())
        nc = 4; nr = (len(groups)+nc-1)//nc
        fig, axes = plt.subplots(nrows=nr, ncols=nc, figsize=(16, 2.5*nr), squeeze=False)
        for ax, gi in zip(axes.ravel(), groups):
            d = summary[(summary["group"]==gi) & (summary["n_fam"]>=2)].head(8)
            if d.empty:
                ax.set_axis_off(); ax.set_title(gi, fontsize=9); continue
            ax.barh(np.arange(len(d)), d["mean_log2OR"].values,
                    color="#E6550D", edgecolor="black", linewidth=0.3)
            ax.set_yticks(np.arange(len(d)))
            ax.set_yticklabels([f"{s} (n={n})" for s, n in
                                zip(d["superfam"], d["n_fam"])], fontsize=7)
            ax.invert_yaxis()
            ax.set_title(f"{gi} ({d['n_sp'].iloc[0]} sp)", fontsize=9)
            ax.set_xlabel("mean log2 OR (Compacted)", fontsize=7)
            ax.tick_params(axis="x", labelsize=7)
        for ax in axes.ravel()[len(groups):]:
            ax.set_axis_off()
        fig.suptitle("TE superfamilies most enriched in Compacted (per group, simple 2-feature model)",
                     fontsize=11)
        fig.tight_layout(rect=[0,0,1,0.97])
        fig.savefig(os.path.join(OUT_ROOT, "Compacted_top_families_per_group.pdf"))
        plt.close(fig)
    print("[DONE]", flush=True)


if __name__ == "__main__":
    main()
