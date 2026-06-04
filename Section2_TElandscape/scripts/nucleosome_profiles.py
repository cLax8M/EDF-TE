#!/usr/bin/env python3
"""
nucleosome_profiles.py
=======================
Compute average nucleosome-occupancy profiles over TE insertions and genes
(meta-TE / metagene style) for the 4 in vivo species.

Layout: one row per in vivo species (Rhimi59 MNase, Ncrass2 MNase,
        Bcine1 MAINE, Mereb1 ATAC), one column per feature class
        (LTR, LINE, DNA, PCG — plus Unknown if enough data).

Each profile:  [FLANK_BP upstream] [body scaled to N_BODY bins] [FLANK_BP downstream]
               shown as mean ± shaded 95 % CI (bootstrap SE × 1.96)

Signal is normalised per feature to its mean flank signal so curves are
comparable across species with very different absolute signal levels.

Outputs (invivo_landscape/):
  figH1_nucleosome_profiles.{pdf,png}
  figH1_profiles_summary.tsv   (mean profiles per species × class)
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chromatin_states_from_nuclpred import load_nuclpred_wig
from clade_palette import TE_CLASS_COLORS

# ----------------------------------------------------------------- paths ---
GFF_DIR  = "/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff"
FAI_DIR  = "/legserv/Temp/Carlos/EDF-TE/data/genomes"
TE_DATA  = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
TE_MV    = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
MV_NUC   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes/nuc-seq"
MV_ATC   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes/atac-seq"
OUT      = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/invivo_landscape"
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- params ---
FLANK_BP    = 2000    # flanking bp each side
N_BODY_BINS = 100     # bins for scaled body (TEs only)
N_FLANK_BINS= 50      # bins per flank
MIN_TE_LEN  = 300     # skip very short TEs (bp)
MAX_TE_LEN  = 80_000  # skip very long TEs (bp)
MIN_GENE_LEN= 300
MAX_GENE_LEN= 80_000
MAX_FEATURES= 5000    # max features per class to sample (for speed)

# For body-scaled (TEs): N_FLANK + N_BODY + N_FLANK bins
N_TOTAL_TE    = N_FLANK_BINS + N_BODY_BINS + N_FLANK_BINS   # 200

# For TSS-centered (PCG): fixed window [-2000 .. +4000] bp = 300 bp/bin resolution
TSS_UP_BP    = 2000   # bp upstream of TSS
TSS_DN_BP    = 4000   # bp downstream of TSS
BIN_SIZE_TSS = 20     # bp per bin
N_TOTAL_TSS  = (TSS_UP_BP + TSS_DN_BP) // BIN_SIZE_TSS  # 300

CLASSES_SHOW = ["LTR", "LINE", "DNA", "Unknown", "PCG"]

SF_COLORS = {
    "LTR":      TE_CLASS_COLORS["LTR"],
    "LINE":     TE_CLASS_COLORS["LINE"],
    "DNA":      TE_CLASS_COLORS["DNA"],
    "Helitron": TE_CLASS_COLORS.get("Rolling Circle","#ff7f00"),
    "PLE":      TE_CLASS_COLORS["PLE"],
    "SINE":     TE_CLASS_COLORS["SINE"],
    "Unknown":  TE_CLASS_COLORS["Unknown"],
    "PCG":      "#333333",
}

IN_VIVO_PANEL = [
    dict(label="MU_Rhimi59",     sp="MU_Rhimi59",  assay="MNase",
         wig="/legserv/Temp/Carlos/EDF-TE/data/nucleosomes/Rhizopus_microsporus.wig"),
    dict(label="AS_Ncrass2",     sp="AS_Ncrass2",  assay="MNase",
         wig=f"{MV_NUC}/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.wig"),
    dict(label="AS_Bcine1",      sp="AS_Bcine1",   assay="MAINE",
         wig=f"{MV_NUC}/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.wig"),
    dict(label="BA_Mereb1",      sp="BA_Mereb1",   assay="ATAC",
         wig=f"{MV_ATC}/BA_Mereb1/BA_Mereb1_atac.wig"),
]

# ---------------------------------------------------------------- helpers ---
def classify_te_sf(col4):
    raw = col4.strip(); c = raw.split("/")[0].rstrip("?").strip()
    if c.startswith("LTR"):                                      return "LTR"
    if c.startswith("LINE"):                                     return "LINE"
    if c.startswith("DNA"):                                      return "DNA"
    if c.startswith("RC") or c.lower().startswith("helitron"):  return "Helitron"
    if c.startswith("PLE") or c.startswith("Penelope"):         return "PLE"
    if c.startswith("SINE"):                                     return "SINE"
    if c in ("Unknown","Retroposon","Unspecified","Other"):      return "Unknown"
    return None

def te_bed_path(sp):
    a = TE_DATA.format(sp=sp); b = TE_MV.format(sp=sp)
    return a if os.path.exists(a) else (b if os.path.exists(b) else None)

def gff_path(sp):
    for ext in (".gff",".gff3",".gff2"):
        p = os.path.join(GFF_DIR, sp + ext)
        if os.path.exists(p): return p
    return None

def chrom_sizes_fai(sp):
    fai = os.path.join(FAI_DIR, f"{sp}.fasta.fai")
    if not os.path.exists(fai): return None
    df = pd.read_csv(fai, sep="\t", header=None, names=["c","L","o","lb","lby"])
    return dict(zip(df["c"], df["L"].astype(int)))


def extract_profile(occ_arr, start, end, chrom_len, resolution=1):
    """
    Body-scaled profile for TEs: [FLANK | body → N_BODY_BINS | FLANK]
    Returns array of length N_TOTAL_TE or None.
    """
    feat_len = end - start
    if feat_len < MIN_TE_LEN: return None

    L = len(occ_arr)
    left  = start - FLANK_BP
    right = end   + FLANK_BP
    if left < 0 or right > L: return None

    lf = occ_arr[left:start]
    lf_bins = np.interp(np.linspace(0, len(lf)-1, N_FLANK_BINS),
                        np.arange(len(lf)), lf)
    bd = occ_arr[start:end]
    bd_bins = np.interp(np.linspace(0, len(bd)-1, N_BODY_BINS),
                        np.arange(len(bd)), bd)
    rf = occ_arr[end:right]
    rf_bins = np.interp(np.linspace(0, len(rf)-1, N_FLANK_BINS),
                        np.arange(len(rf)), rf)

    profile = np.concatenate([lf_bins, bd_bins, rf_bins])
    flank_mean = np.mean(np.concatenate([lf_bins, rf_bins]))
    if flank_mean < 1e-6: return None
    return profile / flank_mean


def extract_tss_profile(occ_arr, tss, strand_plus):
    """
    TSS-centered profile for PCG genes (strand-aware).
    Extracts a fixed window: [TSS - TSS_UP_BP .. TSS + TSS_DN_BP] for + strand
                             [TES - TSS_DN_BP .. TES + TSS_UP_BP] (reversed) for - strand
    Returns array of length N_TOTAL_TSS or None.
    """
    L = len(occ_arr)
    if strand_plus:
        left  = tss - TSS_UP_BP
        right = tss + TSS_DN_BP
    else:
        # tss is the 3' end of feature; reverse so 5'→3' direction
        left  = tss - TSS_DN_BP
        right = tss + TSS_UP_BP

    if left < 0 or right > L: return None

    seg = occ_arr[left:right]
    if not strand_plus:
        seg = seg[::-1]   # flip to 5'→3'

    bins = np.interp(np.linspace(0, len(seg)-1, N_TOTAL_TSS),
                     np.arange(len(seg)), seg)

    # normalise to upstream flank
    up_bins = int(N_TOTAL_TSS * TSS_UP_BP / (TSS_UP_BP + TSS_DN_BP))
    flank_mean = np.mean(bins[:up_bins])
    if flank_mean < 1e-6: return None
    return bins / flank_mean


def collect_features_from_bed(bed_path, classes):
    """Returns dict class → list of (chrom, start, end)."""
    feats = {c: [] for c in classes}
    if not bed_path or not os.path.exists(bed_path): return feats
    with open(bed_path) as fh:
        for ln in fh:
            f = ln.split("\t")
            if len(f) < 4: continue
            sf = classify_te_sf(f[3])
            if sf not in feats: continue
            try: s, e = int(f[1]), int(f[2])
            except ValueError: continue
            if (e - s) < MIN_TE_LEN or (e - s) > MAX_TE_LEN: continue
            feats[sf].append((f[0], s, e))
    return feats


def collect_features_from_gff(gff_p):
    """Returns list of (chrom, start, end, strand_is_plus) for gene features."""
    feats = []
    if not gff_p or not os.path.exists(gff_p): return feats
    gene_feat = "gene"
    with open(gff_p) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.split("\t")
            if len(f) < 3: continue
            if f[2] == "gene": break
            if f[2] in ("exon","CDS"): gene_feat = f[2]; break
    with open(gff_p) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.split("\t")
            if len(f) < 7 or f[2] != gene_feat: continue
            try: s, e = int(f[3]) - 1, int(f[4])
            except ValueError: continue
            if (e - s) < MIN_GENE_LEN or (e - s) > MAX_GENE_LEN: continue
            strand_plus = (f[6].strip() != "-")  # + or . treated as forward
            feats.append((f[0], s, e, strand_plus))
    return feats


def compute_profiles_for_species(entry, classes_te=None, do_pcg=True):
    """
    Returns dict class → (mean_profile, se_profile, n) arrays of length N_TOTAL_BINS.
    """
    if classes_te is None:
        classes_te = ["LTR","LINE","DNA","Helitron","Unknown"]

    sp  = entry["sp"]
    wig = entry["wig"]
    if not os.path.exists(wig):
        print(f"  [SKIP] WIG not found: {wig}"); return {}

    sz   = chrom_sizes_fai(sp)
    arrs = load_nuclpred_wig(wig, sz)
    if sz is None:
        sz = {c: len(a) for c, a in arrs.items()}

    # collect annotations
    bed_p  = te_bed_path(sp)
    gff_p  = gff_path(sp)
    te_feats  = collect_features_from_bed(bed_p, classes_te)
    pcg_feats = collect_features_from_gff(gff_p) if do_pcg else []

    result = {}
    all_classes = [(cls, te_feats.get(cls,[])) for cls in classes_te]
    if do_pcg:
        all_classes.append(("PCG", pcg_feats))

    for cls, feats in all_classes:
        if not feats: continue
        is_pcg = (cls == "PCG")
        rng = np.random.default_rng(42)
        if len(feats) > MAX_FEATURES:
            idx = rng.choice(len(feats), MAX_FEATURES, replace=False)
            feats = [feats[i] for i in idx]

        profiles = []
        for feat in feats:
            chrom = feat[0]; s = feat[1]; e = feat[2]
            occ = arrs.get(chrom)
            if occ is None: continue

            if is_pcg:
                # TSS-centered: feat has (chrom, start, end, strand_plus)
                strand_plus = feat[3] if len(feat) > 3 else True
                tss = s if strand_plus else e
                p = extract_tss_profile(occ, tss, strand_plus)
            else:
                p = extract_profile(occ, s, e, sz.get(chrom, len(occ)))

            if p is not None: profiles.append(p)

        if len(profiles) < 20:
            print(f"    {cls}: only {len(profiles)} profiles — skip")
            continue

        mat = np.vstack(profiles)
        mean_p = mat.mean(axis=0)
        se_p   = mat.std(axis=0) / np.sqrt(mat.shape[0])
        result[cls] = (mean_p, se_p, mat.shape[0])
        print(f"    {cls:12s}: {mat.shape[0]:5d} features")

    return result


# ================================================================ figure ===

def make_figure(profiles_per_species):
    """
    profiles_per_species: list of (entry_dict, result_dict)
    result_dict: class → (mean, se, n)
    For TE classes: profile length = N_TOTAL_TE (200)
    For PCG: profile length = N_TOTAL_TSS (300, TSS-centered)
    """
    all_classes = []
    for _, res in profiles_per_species:
        for cls in CLASSES_SHOW:
            if cls in res and cls not in all_classes:
                all_classes.append(cls)

    n_sp  = len(profiles_per_species)
    n_cls = len(all_classes)

    # TE x-axis (body-scaled)
    x_te    = np.arange(N_TOTAL_TE)
    body_start = N_FLANK_BINS; body_end = N_FLANK_BINS + N_BODY_BINS
    xticks_te  = [0, N_FLANK_BINS//2, body_start,
                  body_start + N_BODY_BINS//4,
                  body_start + N_BODY_BINS//2,
                  body_start + 3*N_BODY_BINS//4,
                  body_end, body_end + N_FLANK_BINS//2, N_TOTAL_TE - 1]
    xlabs_te   = [f"-{FLANK_BP}bp", f"-{FLANK_BP//2}bp", "Start",
                  "25%", "50%", "75%", "End",
                  f"+{FLANK_BP//2}bp", f"+{FLANK_BP}bp"]

    # PCG x-axis (TSS-centered)
    x_tss  = np.arange(N_TOTAL_TSS)
    tss_bin = int(N_TOTAL_TSS * TSS_UP_BP / (TSS_UP_BP + TSS_DN_BP))
    xticks_tss = [0, tss_bin//2, tss_bin,
                  tss_bin + (N_TOTAL_TSS - tss_bin)//4,
                  tss_bin + (N_TOTAL_TSS - tss_bin)//2,
                  N_TOTAL_TSS - 1]
    xlabs_tss  = [f"-{TSS_UP_BP}bp", f"-{TSS_UP_BP//2}bp", "TSS",
                  f"+{TSS_DN_BP//4}bp", f"+{TSS_DN_BP//2}bp", f"+{TSS_DN_BP}bp"]

    fig_h = 2.8 * n_sp + 1.0
    fig, axes = plt.subplots(n_sp, n_cls,
                             figsize=(3.2 * n_cls + 0.8, fig_h),
                             sharey="row")
    fig.patch.set_facecolor("#F9FAFB")

    if n_sp == 1:   axes = axes[np.newaxis, :]
    if n_cls == 1:  axes = axes[:, np.newaxis]

    for row, (entry, res) in enumerate(profiles_per_species):
        lab   = entry["label"]
        assay = entry["assay"]
        for col, cls in enumerate(all_classes):
            ax = axes[row, col]
            ax.set_facecolor("#F9FAFB")
            is_pcg = (cls == "PCG")
            x      = x_tss if is_pcg else x_te
            xticks = xticks_tss if is_pcg else xticks_te
            xlabs  = xlabs_tss  if is_pcg else xlabs_te

            ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.5)

            if is_pcg:
                # mark TSS
                ax.axvline(tss_bin, color="#555", lw=0.9, ls=":", alpha=0.8)
            else:
                ax.axvspan(body_start, body_end, alpha=0.06,
                           color=SF_COLORS.get(cls,"#888888"), zorder=0)
                ax.axvline(body_start, color="#555", lw=0.7, ls=":", alpha=0.7)
                ax.axvline(body_end,   color="#555", lw=0.7, ls=":", alpha=0.7)

            if cls in res:
                mean_p, se_p, n = res[cls]
                color = SF_COLORS.get(cls,"#888888")
                ax.plot(x, mean_p, color=color, lw=1.6, zorder=3)
                ax.fill_between(x, mean_p - 1.96*se_p, mean_p + 1.96*se_p,
                                color=color, alpha=0.22, zorder=2)
                ax.text(0.97, 0.97, f"n={n:,}", transform=ax.transAxes,
                        ha="right", va="top", fontsize=7.5, color="#333333")
                # robust y-limits: mean ± 3 sd, clipped to avoid outlier dominance
                spread = np.nanstd(mean_p)
                mid    = np.nanmedian(mean_p)
                ax.set_ylim(max(0, mid - 3.5*spread), mid + 4.5*spread)
            else:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=8, color="#aaa")

            ax.set_xticks(xticks); ax.set_xticklabels(xlabs, fontsize=6.5,
                                                       rotation=35, ha="right")
            ax.grid(axis="y", alpha=0.2, lw=0.5)
            ax.tick_params(axis="y", labelsize=7)

            if row == 0:
                title_str = cls if not is_pcg else "PCG\n(TSS-centered)"
                ax.set_title(title_str, fontsize=11, fontweight="bold",
                             color=SF_COLORS.get(cls,"#333333"), pad=5)
            if col == 0:
                ax.set_ylabel(f"{lab}\n({assay})\nnorm. signal", fontsize=8)

    fig.suptitle("Average nucleosome-occupancy profiles\n"
                 "TEs: body-scaled ±2kb flanks | PCG: TSS-centered (strand-aware), mean ± 95% CI",
                 fontsize=11, fontweight="bold", y=1.01)
    fig.tight_layout(rect=[0,0,1,1])

    for ext in ("pdf","png"):
        kw = {"bbox_inches":"tight"} if ext=="pdf" else {"dpi":200,"bbox_inches":"tight"}
        fig.savefig(os.path.join(OUT, f"figH1_nucleosome_profiles.{ext}"), **kw)
    plt.close(fig)
    print("  → figH1_nucleosome_profiles saved")


# ================================================================ TSV out ===

def save_tsv(profiles_per_species):
    rows = []
    for entry, res in profiles_per_species:
        for cls, (mean_p, se_p, n) in res.items():
            for i in range(len(mean_p)):
                rows.append({
                    "species": entry["label"], "assay": entry["assay"],
                    "class": cls, "bin": i,
                    "mean": float(mean_p[i]), "se": float(se_p[i]), "n": n,
                    "profile_type": "TSS-centered" if cls == "PCG" else "body-scaled",
                })
    pd.DataFrame(rows).to_csv(
        os.path.join(OUT, "figH1_profiles_summary.tsv"), sep="\t", index=False)
    print("  TSV saved")


# ================================================================== main ===

def main():
    print("\n" + "="*60)
    print("Nucleosome profiles over TE classes and PCG")
    print("="*60)

    profiles_per_species = []
    for entry in IN_VIVO_PANEL:
        print(f"\n[{entry['label']}]")
        res = compute_profiles_for_species(
            entry,
            classes_te=["LTR","LINE","DNA","Helitron","Unknown"],
            do_pcg=True)
        profiles_per_species.append((entry, res))

    make_figure(profiles_per_species)
    save_tsv(profiles_per_species)
    print(f"\n[DONE] outputs in {OUT}")


if __name__ == "__main__":
    main()
