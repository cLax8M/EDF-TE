#!/usr/bin/env python3
"""
calibrate_and_aggregate.py
==========================

(1) Calibra el modelo nucleosome-only en Rhizopus microsporus comparando
    el estado 'Fuzzy_high' con los picos broad de H3K9me3.
(2) Recorre todas las especies con (nuclpred WIG) ∩ (genome.fai) ∩ (TE BED)
    y ejecuta chromatin_states_from_nuclpred.py.
(3) Agrega:
       - composición de estados por especie
       - log2OR(superfamilia × estado) por especie
       - heatmap cross-species de log2OR(LTR/Gypsy etc. en Fuzzy_high)
"""
import os, sys, glob, subprocess, json, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PYTHON   = "/home/carlos/micromamba/envs/syntenyEDFTE/bin/python"
SCRIPT   = "/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts/chromatin_states_from_nuclpred.py"
NUCL_DIR = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosomes/nuclpred"
GENOMES  = "/legserv/Temp/Carlos/EDF-TE/data/genomes"
TE_DIR   = "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_ROOT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/nucleosome_states_pan"
WIN      = 500

H3K9_BED = "/legserv/Temp/Carlos/EDF-TE/data/Genome_Assemblies_Conditions/Epigenetics/RmH3K9me3/H3K9_WT_MACS2_broad_peaks.bed"
RHIMI59  = "MU_Rhimi59"
# Real MNase-seq DANPOS WIG (79 scaffolds), used ONLY for Rhimi59 calibration.
RHIMI59_DANPOS_WIG = "/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"

os.makedirs(OUT_ROOT, exist_ok=True)

# ---------- enumerate species ----------
def find_species():
    species = []
    for w in sorted(glob.glob(os.path.join(NUCL_DIR, "*_nuclpred.wig"))):
        sp = os.path.basename(w).replace("_nuclpred.wig","")
        fai = os.path.join(GENOMES, f"{sp}.fasta.fai")
        # buscar en TE_DIR la carpeta <sp>_EarlGrey
        te_glob = glob.glob(os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles", f"{sp}.filteredRepeats.bed"))
        if not os.path.exists(fai) or not te_glob: continue
        species.append((sp, w, fai, te_glob[0]))
    return species

# ---------- run one species ----------
def run_one(sp, wig, fai, te, outdir):
    if os.path.exists(os.path.join(outdir, "TE_enrichment_per_state.tsv")):
        return True
    os.makedirs(outdir, exist_ok=True)
    log = os.path.join(outdir, "run.log")
    cmd = [PYTHON, SCRIPT,
           "--species", sp, "--nuclpred", wig,
           "--fai", fai, "--te", te,
           "--outdir", outdir, "--win", str(WIN)]
    with open(log, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    return r.returncode == 0

# ---------- calibration on Rhizopus ----------
def calibrate():
    sp = RHIMI59
    outdir = os.path.join(OUT_ROOT, sp + "_DANPOS")
    if not os.path.exists(os.path.join(outdir, "chromatin_states.tsv.gz")):
        print(f"[CAL] running {sp} (DANPOS MNase-seq, 79 scaffolds) first ...")
        wig = RHIMI59_DANPOS_WIG
        fai = os.path.join(GENOMES, f"{sp}.fasta.fai")
        te  = os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles", f"{sp}.filteredRepeats.bed")
        run_one(sp, wig, fai, te, outdir)
    print("[CAL] computing overlap with H3K9me3 broad peaks ...")
    win = pd.read_csv(os.path.join(outdir, "chromatin_states.tsv.gz"), sep="\t")
    h3 = pd.read_csv(H3K9_BED, sep="\t", header=None,
                     names=["chrom","start","end"], usecols=[0,1,2])
    # build per-chrom interval list
    h3_by = {c: h3.loc[h3["chrom"]==c, ["start","end"]].values for c in h3["chrom"].unique()}
    # for each window, fraction overlap with H3K9me3
    frac = np.zeros(len(win))
    for idx, row in enumerate(win.itertuples(index=False)):
        if row.chrom not in h3_by: continue
        ws, we = row.start, row.end
        peaks = h3_by[row.chrom]
        if len(peaks) == 0: continue
        ov = np.maximum(0, np.minimum(peaks[:,1], we) - np.maximum(peaks[:,0], ws)).sum()
        frac[idx] = ov / (we - ws)
    win = win.copy()
    win["h3k9_frac"] = frac
    win["h3k9_pos"]  = (frac >= 0.5).astype(int)
    # confusion
    tab = pd.crosstab(win["state"], win["h3k9_pos"]).reindex(index=["No_signal","Open","Phased","Fuzzy_mid","Fuzzy_high"]).fillna(0).astype(int)
    if 1 not in tab.columns: tab[1] = 0
    if 0 not in tab.columns: tab[0] = 0
    tab["total"] = tab.sum(axis=1)
    tab["frac_H3K9"] = tab[1] / tab["total"]
    tab.to_csv(os.path.join(OUT_ROOT, "calibration_Rhimi59.tsv"), sep="\t")
    # log2OR fold of H3K9 in each state vs background
    bg_pos = win["h3k9_pos"].mean()
    tab["log2OR_H3K9"] = np.log2((tab["frac_H3K9"]+1e-9) / (bg_pos+1e-9))
    print(tab)
    # plot
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    bars = ax.bar(tab.index.astype(str), tab["frac_H3K9"]*100,
                  color=["#CCCCCC","#6BAED6","#74C476","#FDAE6B","#E6550D"])
    ax.axhline(bg_pos*100, color="k", ls="--", lw=1, label=f"genome bg ({bg_pos*100:.1f}%)")
    for b, v, lo in zip(bars, tab["frac_H3K9"]*100, tab["log2OR_H3K9"]):
        ax.text(b.get_x()+b.get_width()/2, v+0.5, f"{v:.1f}%\nlog2OR={lo:+.2f}",
                ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("% windows overlapping H3K9me3 (≥50%)")
    ax.set_title(f"Calibration in {sp}: nucleosome states vs H3K9me3 broad peaks")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "calibration_Rhimi59.pdf"))
    plt.close(fig)
    return tab

# ---------- aggregation ----------
def aggregate(species_done):
    # composition per species
    comp_rows = []
    for sp, _, _, _ in species_done:
        f = os.path.join(OUT_ROOT, sp, "chromatin_states_summary.tsv")
        if not os.path.exists(f): continue
        s = pd.read_csv(f, sep="\t").set_index("state")["bp"]
        s = s / s.sum()
        for st, v in s.items():
            comp_rows.append(dict(species=sp, state=st, frac=v))
    comp = pd.DataFrame(comp_rows).pivot(index="species", columns="state", values="frac").fillna(0)
    comp = comp[[c for c in ["No_signal","Open","Phased","Fuzzy_mid","Fuzzy_high"] if c in comp.columns]]
    comp.to_csv(os.path.join(OUT_ROOT, "composition_per_species.tsv"), sep="\t")

    # plot composition
    fig, ax = plt.subplots(figsize=(11, max(4, 0.18*len(comp))))
    bottom = np.zeros(len(comp))
    cols = {"No_signal":"#CCCCCC","Open":"#6BAED6","Phased":"#74C476",
            "Fuzzy_mid":"#FDAE6B","Fuzzy_high":"#E6550D"}
    comp_sorted = comp.sort_values("Fuzzy_high")
    for st in comp_sorted.columns:
        ax.barh(comp_sorted.index, comp_sorted[st].values, left=bottom,
                color=cols.get(st,"#888"), label=st, height=0.85)
        bottom += comp_sorted[st].values
    ax.set_xlim(0,1); ax.set_xlabel("Fraction of valid genome (bp)")
    ax.tick_params(axis="y", labelsize=5)
    ax.set_title("Chromatin-state composition (nucleosome-only) across fungi")
    ax.legend(fontsize=8, ncol=4, loc="lower center", bbox_to_anchor=(0.5,1.02))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "composition_per_species.pdf"))
    plt.close(fig)

    # cross-species enrichment (log2OR) for major superfamilies in Fuzzy_high
    majors = ["LTR/Gypsy","LTR/Copia","LINE","LINE/L1","LINE/RTE-BovB","DNA",
              "DNA/CMC-EnSpm","DNA/TcMar","DNA/TcMar-Tc1","DNA/hAT-Ac",
              "RC/Helitron","Unknown","Satellite","Simple_repeat",
              "DNA/Merlin","LTR/ERV1","LTR/ERVK","LINE/Penelope","DNA/PIF-Spy"]
    rows = []
    for sp, _, _, _ in species_done:
        f = os.path.join(OUT_ROOT, sp, "TE_enrichment_per_state.tsv")
        if not os.path.exists(f): continue
        d = pd.read_csv(f, sep="\t")
        d_fh = d[d["state"]=="Fuzzy_high"]
        for m in majors:
            sub = d_fh[d_fh["superfam"]==m]
            if len(sub) == 0: continue
            rows.append(dict(species=sp, superfam=m,
                             log2OR=float(sub["log2OR"].iloc[0]),
                             obs_bp=int(sub["obs_bp"].iloc[0])))
    cross = pd.DataFrame(rows)
    cross.to_csv(os.path.join(OUT_ROOT, "Fuzzy_high_TE_enrichment_cross_species.tsv"), sep="\t", index=False)

    if len(cross) > 0:
        pv = cross.pivot(index="species", columns="superfam", values="log2OR")
        pv = pv.reindex(columns=[m for m in majors if m in pv.columns])
        # only species with at least 5 valid families
        pv = pv.dropna(thresh=5)
        # sort species by phylo prefix (AS_/BA_/MU_…)
        pv = pv.sort_index()
        fig, ax = plt.subplots(figsize=(max(8, 0.45*pv.shape[1]+2), max(5, 0.18*pv.shape[0])))
        vmax = np.nanmax(np.abs(pv.values))
        im = ax.imshow(pv.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(pv.shape[1])); ax.set_xticklabels(pv.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(pv.shape[0])); ax.set_yticklabels(pv.index, fontsize=5)
        cb = fig.colorbar(im, ax=ax, fraction=0.012, pad=0.01)
        cb.set_label("log2(obs/exp) TE bp in Fuzzy_high (heterochromatin proxy)", fontsize=8)
        ax.set_title("TE preference for nucleosome-dense fuzzy regions across fungi")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_ROOT, "Fuzzy_high_TE_heatmap.pdf"))
        plt.close(fig)

# ---------- main ----------
if __name__ == "__main__":
    sps = find_species()
    print(f"[PAN] {len(sps)} species with all inputs available")
    with open(os.path.join(OUT_ROOT, "species_list.tsv"), "w") as fh:
        for sp, w, f, t in sps:
            fh.write(f"{sp}\t{w}\t{f}\t{t}\n")

    # calibration first
    print("=== CALIBRATION ===")
    calibrate()

    # process all (parallel)
    print("=== RUN ALL SPECIES (parallel) ===")
    from concurrent.futures import ProcessPoolExecutor, as_completed
    def _wrap(args):
        sp, w, f, t = args
        outdir = os.path.join(OUT_ROOT, sp)
        ok = run_one(sp, w, f, t, outdir)
        return (sp, w, f, t, ok)
    done = []
    with ProcessPoolExecutor(max_workers=int(os.environ.get("PAN_NPROC", 12))) as ex:
        futs = [ex.submit(_wrap, s) for s in sps]
        for i, fut in enumerate(as_completed(futs), 1):
            sp, w, f, t, ok = fut.result()
            print(f"[PAN {i}/{len(sps)}] {sp} -> {'OK' if ok else 'FAIL'}", flush=True)
            if ok: done.append((sp, w, f, t))

    print("=== AGGREGATE ===")
    aggregate(done)
    print(f"[PAN] DONE. Outputs at {OUT_ROOT}")
