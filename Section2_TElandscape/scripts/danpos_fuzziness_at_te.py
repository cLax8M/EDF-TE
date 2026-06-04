#!/usr/bin/env python3
"""
danpos_fuzziness_at_te.py
=========================
Compare DANPOS per-nucleosome metrics across 3 groups:
  TE    — dyad (smt_pos) falls inside an EarlGrey TE interval
  Gene  — dyad falls inside a protein-coding gene (GFF feature 'gene' or 'exon')
  Inter — neither TE nor gene (intergenic)

Metrics (not affected by multi-mapping coverage bias):
  fuzziness_score  position SD across the read population (lower = better positioned)
  smt_value        smoothed peak height (= sharpness/occupancy of the dyad)
  nrl_bp           distance to the NEXT nucleosome on the same chrom (bp)

Species:
  AS_Ncrass2   MNase in vivo
  AS_Bcine1    MAINE in vivo
  AS_Usnflo1   NucSeq in vitro
  BA_Mereb1    NucSeq in vitro
  BA_Mrafri1   NucSeq in vitro
  BA_Rhoto     NucSeq in vitro
  BL_Catan2    NucSeq in vitro

(MU_Rhimi59 only has a legacy .wig — no positions.xls; excluded here)
(BA_Mereb1_atac is ATAC, no DANPOS — excluded here)

Outputs in
  /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/danpos_fuzziness/
    per_species/<LABEL>_nucleosomes.tsv.gz
    per_species/<LABEL>_summary.tsv
    per_species/<LABEL>_fuzziness.{pdf,png}
    danpos_fuzziness_summary.tsv
    danpos_fuzziness_grid.{pdf,png}
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

MV   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/nucleosomes/nuc-seq"
GFF  = "/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff"

PANEL = [
    dict(label="AS_Ncrass2", sp="AS_Ncrass2", assay="MNase", in_vivo=True,
         xls=f"{MV}/AS_Ncrass2/AS_Ncrass2_danpos/pooled/AS_Ncrass2_dedup.smooth.positions.xls",
         gff=f"{GFF}/AS_Ncrass2.gff"),
    dict(label="AS_Bcine1",  sp="AS_Bcine1",  assay="MAINE", in_vivo=True,
         xls=f"{MV}/AS_Bcine1/AS_Bcine1_danpos/pooled/AS_Bcine1_dedup.smooth.positions.xls",
         gff=f"{GFF}/AS_Bcine1.gff3"),
    dict(label="AS_Usnflo1", sp="AS_Usnflo1", assay="NucSeq", in_vivo=False,
         xls=f"{MV}/AS_Usnflo1/AS_Usnflo1_danpos/pooled/AS_Usnflo1_dedup.smooth.positions.xls",
         gff=f"{GFF}/AS_Usnflo1.gff"),
    dict(label="BA_Mereb1",  sp="BA_Mereb1",  assay="NucSeq", in_vivo=False,
         xls=f"{MV}/BA_Mereb1/BA_Mereb1_danpos/pooled/BA_Mereb1_dedup.smooth.positions.xls",
         gff=f"{GFF}/BA_Mereb1.gff"),
    dict(label="BA_Mrafri1", sp="BA_Mrafri1", assay="NucSeq", in_vivo=False,
         xls=f"{MV}/BA_Mrafri1/BA_Mrafri1_danpos/pooled/BA_Mrafri1_dedup.smooth.positions.xls",
         gff=f"{GFF}/BA_Mrafri1.gff"),
    dict(label="BA_Rhoto",   sp="BA_Rhoto",   assay="NucSeq", in_vivo=False,
         xls=f"{MV}/BA_Rhoto/BA_Rhoto_danpos/pooled/BA_Rhoto_dedup.smooth.positions.xls",
         gff=f"{GFF}/BA_Rhoto.gff"),
    dict(label="BL_Catan2",  sp="BL_Catan2",  assay="NucSeq", in_vivo=False,
         xls=f"{MV}/BL_Catan2/BL_Catan2_danpos/pooled/BL_Catan2_dedup.smooth.positions.xls",
         gff=f"{GFF}/BL_Catan2.gff"),
]

TE_DATA = "/data/carlos/EDF-TE/TE_analysis/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"
TE_MV   = "/legserv/Temp/Carlos/EDF-TE/results/MvFromData/annotated/{sp}_EarlGrey/{sp}_summaryFiles/{sp}.filteredRepeats.bed"

OUT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/danpos_fuzziness"
os.makedirs(os.path.join(OUT, "per_species"), exist_ok=True)


def te_bed(sp):
    a = TE_DATA.format(sp=sp); b = TE_MV.format(sp=sp)
    return a if os.path.exists(a) else b


def load_te_intervals(bed_path):
    """Return dict chrom -> sorted numpy array of (start, end) as Nx2 array."""
    d = {}
    if not os.path.exists(bed_path):
        return d
    with open(bed_path) as fh:
        for ln in fh:
            f = ln.split("\t")
            if len(f) < 3: continue
            c = f[0]
            try: s, e = int(f[1]), int(f[2])
            except ValueError: continue
            d.setdefault(c, []).append((s, e))
    return {c: np.array(sorted(v), dtype=np.int32) for c, v in d.items()}


def load_gene_intervals(gff_path):
    """Load gene body intervals from GFF/GFF3. Uses 'gene' features if present,
    otherwise falls back to 'exon'. Returns dict chrom -> Nx2 array (start,end)."""
    d = {}
    want = set()
    # first pass: detect which feature type to use
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.split("\t")
            if len(f) < 9: continue
            want.add(f[2])
    use = "gene" if "gene" in want else "exon"
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.split("\t")
            if len(f) < 5 or f[2] != use: continue
            c = f[0]
            try: s, e = int(f[3])-1, int(f[4])   # GFF is 1-based inclusive
            except ValueError: continue
            d.setdefault(c, []).append((s, e))
    print(f"  gene intervals ({use}): {sum(len(v) for v in d.values()):,}")
    return {c: np.array(sorted(v), dtype=np.int32) for c, v in d.items()}


def ivs_overlap_vec(df, ivs_dict):
    """Return bool array: True where dyad (smt_pos) falls inside any interval."""
    result = np.zeros(len(df), dtype=bool)
    for c, grp in df.groupby("chr"):
        ivs = ivs_dict.get(c)
        if ivs is None or len(ivs) == 0:
            continue
        pos = grp["smt_pos"].values
        idx = np.searchsorted(ivs[:, 0], pos, side="right") - 1
        in_t = np.where(idx >= 0, pos < ivs[idx.clip(min=0), 1], False)
        result[grp.index.values] = in_t
    return result


def compute_nrl(df):
    """Add nrl_bp = distance to next nucleosome on same chrom (sorted by start)."""
    df = df.sort_values(["chr", "smt_pos"]).copy()
    df["nrl_bp"] = np.nan
    for c, grp in df.groupby("chr"):
        idx = grp.index
        pos = grp["smt_pos"].values
        if len(pos) < 2: continue
        dists = np.diff(pos).astype(float)
        # only keep reasonable NRL (100-400 bp) to exclude gap/assembly breaks
        dists[(dists < 100) | (dists > 400)] = np.nan
        df.loc[idx[:-1], "nrl_bp"] = dists
    return df


def cliffs_d(x, y, nmax=5000):
    x = np.asarray(x); y = np.asarray(y)
    if len(x) == 0 or len(y) == 0: return np.nan
    rng = np.random.RandomState(42)
    if len(x) > nmax: x = rng.choice(x, nmax, replace=False)
    if len(y) > nmax: y = rng.choice(y, nmax, replace=False)
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return float(gt - lt) / (len(x) * len(y))


def run_species(entry):
    label = entry["label"]; sp = entry["sp"]
    print(f"[{label}] loading positions...", flush=True)
    df = pd.read_csv(entry["xls"], sep="\t")
    print(f"  {len(df):,} nucleosomes", flush=True)

    te_ivs   = load_te_intervals(te_bed(sp))
    gene_ivs = load_gene_intervals(entry["gff"])

    df["in_te"]   = ivs_overlap_vec(df, te_ivs)
    df["in_gene"] = ivs_overlap_vec(df, gene_ivs)
    # Priority: TE > Gene > Intergenic (TE wins if overlapping both)
    df["group"] = "Intergenic"
    df.loc[df["in_gene"] & ~df["in_te"], "group"] = "Gene"
    df.loc[df["in_te"],                  "group"] = "TE"

    df = compute_nrl(df)

    # save
    df.to_csv(os.path.join(OUT, "per_species", f"{label}_nucleosomes.tsv.gz"),
              sep="\t", index=False, compression="gzip")

    counts = df["group"].value_counts()
    print(f"  {counts.to_dict()}", flush=True)

    rows = {"species": label, "assay": entry["assay"],
            "in_vivo": entry["in_vivo"],
            "n_nuc_total": len(df),
            "n_te": int((df["group"]=="TE").sum()),
            "n_gene": int((df["group"]=="Gene").sum()),
            "n_inter": int((df["group"]=="Intergenic").sum())}

    # Stats: TE vs Gene (the main comparison)
    for col in ["fuzziness_score", "smt_value", "nrl_bp"]:
        for pair, k in [("TE","Gene"), ("TE","Intergenic"), ("Gene","Intergenic")]:
            x = df.loc[df["group"]==pair, col].dropna().values
            y = df.loc[df["group"]==k,    col].dropna().values
            key = f"{col}__{pair}_vs_{k}"
            if len(x) >= 3 and len(y) >= 3:
                try: _, p = mannwhitneyu(x, y, alternative="two-sided")
                except Exception: p = np.nan
                rows[f"{key}_p"] = float(p)
                rows[f"{key}_d"] = cliffs_d(x, y)
                rows[f"{key}_delta"] = float(np.mean(x) - np.mean(y))
            else:
                rows[f"{key}_p"] = np.nan; rows[f"{key}_d"] = np.nan; rows[f"{key}_delta"] = np.nan
        # group means
        for g in ["TE","Gene","Intergenic"]:
            v = df.loc[df["group"]==g, col].dropna().values
            rows[f"{col}_mean_{g}"] = float(np.mean(v)) if len(v) else np.nan

    # print compact summary
    for col in ["fuzziness_score","smt_value","nrl_bp"]:
        print(f"  {col}: TE={rows[f'{col}_mean_TE']:.2f}  Gene={rows[f'{col}_mean_Gene']:.2f}  "
              f"Inter={rows[f'{col}_mean_Intergenic']:.2f}  "
              f"d(TE-Gene)={rows[f'{col}__TE_vs_Gene_d']:.3f}", flush=True)
    return df, rows


def plot_species(label, df, assay, in_vivo, out_pdf):
    cols = [
        ("fuzziness_score", "Fuzziness (lower = well-positioned)"),
        ("smt_value",       "Sharpness / smt_value (log)"),
        ("nrl_bp",          "NRL — dist. to next nucleosome (bp)"),
    ]
    ORDER  = ["Gene", "Intergenic", "TE"]
    COLORS = {"Gene": "#2CA02C", "Intergenic": "#7F7F7F", "TE": "#D62728"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, (col, ttl) in zip(axes, cols):
        data = [df.loc[df["group"]==g, col].dropna().values for g in ORDER]
        bp = ax.boxplot(data, positions=np.arange(len(ORDER)), widths=0.6,
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=1.5))
        for patch, g in zip(bp["boxes"], ORDER):
            patch.set_facecolor(COLORS[g]); patch.set_edgecolor("black")
        ax.set_xticks(np.arange(len(ORDER)))
        ax.set_xticklabels([f"{g}\n(n={(df['group']==g).sum():,})" for g in ORDER], fontsize=9)
        ax.set_title(ttl, fontsize=10)
        if col == "smt_value": ax.set_yscale("log")
        if col == "nrl_bp":    ax.set_ylim(100, 400)
    tag = f"{'in vivo' if in_vivo else 'in vitro'} {assay}"
    fig.suptitle(f"{label} ({tag}) — per-nucleosome DANPOS metrics: Gene vs TE vs Intergenic", fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.94])
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf",".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def grid_plot(labels, dfs, entries, fname):
    cols = [
        ("fuzziness_score", "Fuzziness score"),
        ("smt_value",       "Sharpness (log)"),
        ("nrl_bp",          "NRL (bp)"),
    ]
    ORDER  = ["Gene", "Intergenic", "TE"]
    COLORS = {"Gene": "#2CA02C", "Intergenic": "#7F7F7F", "TE": "#D62728"}
    fig, axes = plt.subplots(len(cols), len(labels),
                             figsize=(2.6*len(labels), 3.2*len(cols)), squeeze=False)
    for j, lab in enumerate(labels):
        df = dfs[lab]; meta = entries[lab]
        for i, (col, ttl) in enumerate(cols):
            ax = axes[i][j]
            data = [df.loc[df["group"]==g, col].dropna().values for g in ORDER]
            bp = ax.boxplot(data, positions=np.arange(3), widths=0.6,
                            showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.2))
            for patch, g in zip(bp["boxes"], ORDER):
                patch.set_facecolor(COLORS[g]); patch.set_edgecolor("black")
            ax.set_xticks(np.arange(3)); ax.set_xticklabels(ORDER, fontsize=7, rotation=20, ha="right")
            if col == "smt_value": ax.set_yscale("log")
            if col == "nrl_bp":    ax.set_ylim(100, 400)
            if i == 0:
                tag = f"{'in vivo' if meta['in_vivo'] else 'in vitro'} {meta['assay']}"
                ax.set_title(f"{lab}\n{tag}", fontsize=8.5, fontweight="bold")
            if j == 0: ax.set_ylabel(ttl, fontsize=9)
    fig.suptitle("Per-nucleosome DANPOS metrics: Gene vs TE vs Intergenic", fontsize=12)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(os.path.join(OUT, f"{fname}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{fname}.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    all_dfs = {}; all_rows = []; all_entries = {}
    for entry in PANEL:
        if not os.path.exists(entry["xls"]):
            print(f"[SKIP] {entry['label']}: {entry['xls']}"); continue
        df, rows = run_species(entry)
        plot_species(entry["label"], df, entry["assay"], entry["in_vivo"],
                     os.path.join(OUT, "per_species", f"{entry['label']}_fuzziness.pdf"))
        all_dfs[entry["label"]] = df
        all_rows.append(rows)
        all_entries[entry["label"]] = entry

    summ = pd.DataFrame(all_rows)
    cols_first = ["species","assay","in_vivo","n_nuc_total","n_te","n_gene","n_inter"]
    rest = [c for c in summ.columns if c not in cols_first]
    summ = summ[cols_first + rest].sort_values(["in_vivo","species"], ascending=[False,True])
    summ.to_csv(os.path.join(OUT, "danpos_fuzziness_summary.tsv"), sep="\t", index=False)

    print("\n=== TE vs Gene (Cliff's d) ===")
    show_cols = ["species","assay","in_vivo"]
    for col in ["fuzziness_score","smt_value","nrl_bp"]:
        show_cols += [f"{col}__TE_vs_Gene_d", f"{col}__TE_vs_Gene_delta", f"{col}__TE_vs_Gene_p"]
    print(summ[show_cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    labs = list(summ["species"])
    grid_plot(labs, all_dfs, all_entries, "danpos_fuzziness_grid")

    inv = [l for l in labs if all_entries[l]["in_vivo"]]
    if inv:
        grid_plot(inv, {k: all_dfs[k] for k in inv},
                  {k: all_entries[k] for k in inv}, "invivo_only_fuzziness_grid")
    print(f"\n[DONE] {OUT}")


if __name__ == "__main__":
    main()
