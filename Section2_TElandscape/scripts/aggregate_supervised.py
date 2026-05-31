#!/usr/bin/env python3
"""Aggregate already-computed per-species supervised outputs.
Reads cross_species/*_composition.tsv and *_TE_enrichment.tsv;
writes composition_supervised.tsv/.pdf and Het_TE_log2OR_heatmap.pdf.
"""
import os, glob
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from functools import reduce

OUT_ROOT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin"
CS_DIR   = os.path.join(OUT_ROOT, "cross_species")

STATE_ORDER  = ["No_signal","Open","Intermediate","Heterochromatin"]
STATE_COLORS = {"No_signal":"#CCCCCC","Open":"#6BAED6",
                "Intermediate":"#FDAE6B","Heterochromatin":"#E6550D"}

# composition
rows = []
for f in sorted(glob.glob(os.path.join(CS_DIR, "*_composition.tsv"))):
    sp = os.path.basename(f).replace("_composition.tsv","")
    s = pd.read_csv(f, sep="\t", index_col=0).iloc[:,0]
    rows.append([sp] + [float(s.get(st,0.0)) for st in STATE_ORDER])
comp = pd.DataFrame(rows, columns=["species"]+STATE_ORDER)
comp["group"] = comp["species"].str.split("_").str[0]
comp.sort_values(["group","species"], inplace=True)
comp.to_csv(os.path.join(OUT_ROOT,"composition_supervised.tsv"), sep="\t", index=False)
print(f"composition: {len(comp)} species; groups: {comp['group'].value_counts().to_dict()}")

# stacked plot
fig, ax = plt.subplots(figsize=(max(8, 0.06*len(comp)), 4))
bottom = np.zeros(len(comp))
for st in STATE_ORDER:
    ax.bar(comp["species"], comp[st], bottom=bottom,
           color=STATE_COLORS[st], width=0.95, label=st)
    bottom += comp[st].values
# group separators
prev = None
for i,g in enumerate(comp["group"].values):
    if g != prev:
        ax.axvline(i-0.5, color="k", lw=0.4, alpha=0.4)
        ax.text(i, 1.02, g, fontsize=6, rotation=90, va="bottom")
    prev = g
ax.set_xticks([]); ax.set_ylabel("Genome fraction"); ax.set_ylim(0,1)
ax.set_title(f"Supervised chromatin composition (n={len(comp)} species, classifier trained on Rhimi59 H3K9me3)")
ax.legend(loc="upper right", fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT_ROOT,"composition_supervised.pdf"))
plt.close(fig)

# Het TE log2OR matrix
hets = []
for f in sorted(glob.glob(os.path.join(CS_DIR, "*_TE_enrichment.tsv"))):
    sp = os.path.basename(f).replace("_TE_enrichment.tsv","")
    if os.path.getsize(f) == 0:
        continue
    try:
        e = pd.read_csv(f, sep="\t")
    except Exception:
        continue
    if "state" not in e.columns or len(e) == 0:
        continue
    e = e[e["state"]=="Heterochromatin"][["superfam","log2OR","obs_bp"]]
    # drop superfams with 0 obs in this species (they'd be -inf-ish noise)
    e = e[e["obs_bp"]>=50]
    e = e.rename(columns={"log2OR": sp})
    hets.append(e[["superfam", sp]])
if hets:
    het_mat = reduce(lambda a,b: a.merge(b, on="superfam", how="outer"), hets)
    het_mat = het_mat.set_index("superfam")
    print(f"raw matrix: {het_mat.shape}")
    # focus on superfams seen in many species
    good = het_mat.notna().sum(axis=1) >= max(10, len(het_mat.columns)//10)
    het_mat = het_mat.loc[good]
    print(f"filtered (n_species>=20% with signal): {het_mat.shape}")
    het_mat["__mean"] = het_mat.mean(axis=1, skipna=True)
    het_mat = het_mat.sort_values("__mean", ascending=False).drop(columns="__mean")
    het_mat = het_mat.head(30)
    het_mat.to_csv(os.path.join(OUT_ROOT,"Het_TE_log2OR_matrix.tsv"), sep="\t")

    cols_order = [c for c in comp["species"].tolist() if c in het_mat.columns]
    het_plot = het_mat[cols_order]
    fig, ax = plt.subplots(figsize=(max(10, 0.06*het_plot.shape[1]),
                                    0.35*het_plot.shape[0]+1.2))
    m = het_plot.values
    vmax = np.nanpercentile(np.abs(m), 95) or 1.0
    im = ax.imshow(m, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   interpolation="nearest")
    ax.set_yticks(range(het_plot.shape[0])); ax.set_yticklabels(het_plot.index, fontsize=8)
    # group ticks
    groups = comp.set_index("species").loc[het_plot.columns,"group"].values
    prev=None
    for i,g in enumerate(groups):
        if g != prev:
            ax.axvline(i-0.5, color="k", lw=0.3, alpha=0.5)
            ax.text(i, -0.5, g, fontsize=6, rotation=90, va="bottom")
        prev=g
    ax.set_xticks([])
    ax.set_xlabel(f"Species (n={het_plot.shape[1]}, ordered by phylogenetic group)")
    ax.set_title("TE superfamily log2OR in Heterochromatin (supervised; trained on Rhimi59)")
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("log2 odds ratio")
    fig.tight_layout(); fig.savefig(os.path.join(OUT_ROOT,"Het_TE_log2OR_heatmap.pdf"))
    plt.close(fig)

# Per-group composition summary
gcomp = comp.groupby("group")[STATE_ORDER].mean().reset_index()
gcomp.to_csv(os.path.join(OUT_ROOT,"composition_per_group.tsv"), sep="\t", index=False)
fig, ax = plt.subplots(figsize=(max(6, 0.5*len(gcomp)), 4))
bottom = np.zeros(len(gcomp))
for st in STATE_ORDER:
    ax.bar(gcomp["group"], gcomp[st], bottom=bottom,
           color=STATE_COLORS[st], width=0.85, label=st)
    bottom += gcomp[st].values
ax.set_ylabel("Mean genome fraction"); ax.set_ylim(0,1)
ax.set_title("Chromatin composition per phylogenetic group (mean across species)")
ax.legend(loc="upper right", fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT_ROOT,"composition_per_group.pdf"))
plt.close(fig)

print("DONE")
