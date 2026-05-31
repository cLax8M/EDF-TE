#!/usr/bin/env python3
"""
two_state_family_enrichment.py
==============================
Reuses per-species supervised outputs (p_het column) to:
  (1) Collapse 4 states -> 2 states: Open vs Compacted.
      Compacted := p_het >= P_HET_THR AND mean_occ > 1e-3
      Open      := mean_occ > 1e-3 AND p_het < P_HET_THR
      (No_signal windows with mean_occ ~ 0 are EXCLUDED from analysis.)
  (2) Compute TE log2 odds-ratio in Compacted at the FAMILY level
      (Name= field from EarlGrey GFF) AND superfamily (col 4 of BED).
  (3) Emit per-species family tables, group-level rankings, and figures.

Outputs (in OUT_ROOT):
  composition_2state.tsv / .pdf
  composition_2state_per_group.tsv / .pdf
  cross_species_2state/<SP>_family_enrichment.tsv
  cross_species_2state/<SP>_superfam_enrichment.tsv
  Compacted_superfam_log2OR_heatmap.pdf
  Compacted_top_families_per_group.tsv
  Compacted_top_families_per_group.pdf
"""
import os, sys, glob, gzip, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
import json

# -------------------------- paths --------------------------
ROOT      = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin"
CS_DIR    = os.path.join(ROOT, "cross_species")
TE_DIR    = "/data/carlos/EDF-TE/TE_analysis/annotated"
OUT_ROOT  = ROOT
OUT2_DIR  = os.path.join(OUT_ROOT, "cross_species_2state")
os.makedirs(OUT2_DIR, exist_ok=True)

WIN = 500
P_HET_THR = json.load(open(os.path.join(ROOT, "learned_thresholds.json")))["p_het"]
print(f"[INFO] using p_het threshold = {P_HET_THR:.4f}")

STATE_ORDER = ["Open", "Compacted"]
STATE_COLOR = {"Open": "#6BAED6", "Compacted": "#E6550D"}

# -------------------------- helpers --------------------------
NAME_RE = re.compile(r"Name=([^;]+)")


def parse_gff_families(gff_path):
    """Return DataFrame: chrom, start, end, superfam, family."""
    rows = []
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            superfam = f[2]
            attr = f[8]
            m = NAME_RE.search(attr)
            fam = m.group(1) if m else "NA"
            try:
                s = int(f[3]) - 1
                e = int(f[4])
            except ValueError:
                continue
            rows.append((f[0], s, e, superfam, fam))
    if not rows:
        return pd.DataFrame(columns=["chrom","start","end","superfam","family"])
    return pd.DataFrame(rows, columns=["chrom","start","end","superfam","family"])


def assign_2state(df, p_thr=P_HET_THR):
    """Add column 'state2' ∈ {Open, Compacted, NoSig}."""
    state = np.full(len(df), "NoSig", dtype=object)
    info  = df["mean_occ"].values > 1e-3
    state[info & (df["p_het"].values >= p_thr)] = "Compacted"
    state[info & (df["p_het"].values <  p_thr)] = "Open"
    df = df.copy()
    df["state2"] = state
    return df


def state_of_window(df_states):
    """Return per-(chrom,window_idx) state lookup."""
    out = {}
    for c, sub in df_states.groupby("chrom", sort=False):
        out[c] = sub["state2"].values
    return out


def intersect_bp(s, e, win=WIN, state_arr=None):
    """For TE [s,e), return dict state -> bp covered."""
    if state_arr is None:
        return {}
    L = len(state_arr) * win
    if s >= L:
        return {}
    e = min(e, L)
    i0 = s // win
    i1 = (e - 1) // win
    out = {}
    for i in range(i0, i1 + 1):
        ws = i * win
        we = ws + win
        ov = min(e, we) - max(s, ws)
        if ov <= 0: continue
        st = state_arr[i]
        out[st] = out.get(st, 0) + ov
    return out


# -------------------------- per-species worker --------------------------
def process_species(args):
    sp, ts_path, gff_path = args
    try:
        df = pd.read_csv(ts_path, sep="\t")
        if "p_het" not in df.columns:
            return sp, None, None, None
        df = assign_2state(df)

        # composition (Open + Compacted only; NoSig excluded)
        info = df[df["state2"] != "NoSig"]
        tot = len(info)
        if tot == 0:
            return sp, None, None, None
        comp = (info["state2"].value_counts(normalize=True)
                .reindex(STATE_ORDER, fill_value=0.0))

        # state per window per chrom
        sw = state_of_window(df)

        # parse TEs
        if not os.path.exists(gff_path):
            return sp, comp, None, None
        te = parse_gff_families(gff_path)
        if te.empty:
            return sp, comp, None, None

        # accumulate bp per (family, superfam, state)
        bp = {}            # (family, superfam, state) -> bp
        n_copies = {}      # (family, superfam) -> int
        for r in te.itertuples(index=False):
            n_copies[(r.family, r.superfam)] = n_copies.get((r.family, r.superfam), 0) + 1
            arr = sw.get(r.chrom)
            if arr is None: continue
            ov = intersect_bp(r.start, r.end, WIN, arr)
            for st, v in ov.items():
                if st == "NoSig": continue
                k = (r.family, r.superfam, st)
                bp[k] = bp.get(k, 0) + v

        # global state bp (Open + Compacted only)
        bp_state = {st: 0 for st in STATE_ORDER}
        for (_, _, st), v in bp.items():
            if st in bp_state:
                bp_state[st] += v
        # we need state genome bp, not TE-only bp; use windows count
        win_state_bp = {st: int((info["state2"]==st).sum()) * WIN for st in STATE_ORDER}
        tot_bp = sum(win_state_bp.values())
        if tot_bp == 0:
            return sp, comp, None, None

        # family enrichment in Compacted
        fam_rows = []
        fams = sorted({(f, s) for (f, s, _) in bp.keys()})
        for fam, sf in fams:
            fam_total = sum(bp.get((fam, sf, st), 0) for st in STATE_ORDER)
            if fam_total < 200:   # skip tiny families
                continue
            obs = bp.get((fam, sf, "Compacted"), 0)
            # expected if uniformly distributed across Open+Compacted bp
            exp = fam_total * (win_state_bp["Compacted"] / tot_bp)
            log2OR = np.log2((obs + 1) / (exp + 1))
            frac_c = obs / fam_total
            fam_rows.append(dict(
                family=fam, superfam=sf,
                n_copies=n_copies.get((fam, sf), 0),
                bp_total=fam_total,
                bp_compacted=obs,
                frac_compacted=frac_c,
                log2OR_compacted=log2OR,
            ))
        fam_df = pd.DataFrame(fam_rows).sort_values("log2OR_compacted", ascending=False)

        # superfam enrichment (aggregate within superfam)
        sf_bp = {}
        for (fam, sf, st), v in bp.items():
            sf_bp.setdefault(sf, {st: 0 for st in STATE_ORDER})
            if st in STATE_ORDER:
                sf_bp[sf][st] = sf_bp[sf].get(st, 0) + v
        sf_rows = []
        for sf, dd in sf_bp.items():
            tot_sf = sum(dd.values())
            if tot_sf < 500: continue
            obs = dd.get("Compacted", 0)
            exp = tot_sf * (win_state_bp["Compacted"] / tot_bp)
            log2OR = np.log2((obs + 1) / (exp + 1))
            sf_rows.append(dict(
                superfam=sf, bp_total=tot_sf,
                bp_compacted=obs,
                frac_compacted=obs/tot_sf,
                log2OR_compacted=log2OR,
            ))
        sf_df = pd.DataFrame(sf_rows).sort_values("log2OR_compacted", ascending=False)

        return sp, comp, fam_df, sf_df
    except Exception as ex:
        print(f"  [WARN] {sp}: {ex}", file=sys.stderr)
        return sp, None, None, None


# -------------------------- main --------------------------
def main():
    ts_files = sorted(glob.glob(os.path.join(CS_DIR, "*_supervised.tsv.gz")))
    jobs = []
    for f in ts_files:
        sp = os.path.basename(f).replace("_supervised.tsv.gz", "")
        gff = os.path.join(TE_DIR, f"{sp}_EarlGrey", f"{sp}_summaryFiles",
                           f"{sp}.filteredRepeats.gff")
        jobs.append((sp, f, gff))
    print(f"[RUN] {len(jobs)} species")

    comps = {}
    sf_all = {}
    n_done = 0
    with ProcessPoolExecutor(max_workers=16) as ex:
        futs = [ex.submit(process_species, j) for j in jobs]
        for fut in as_completed(futs):
            sp, comp, fam_df, sf_df = fut.result()
            n_done += 1
            if comp is not None:
                comps[sp] = comp
            if fam_df is not None and len(fam_df):
                fam_df.to_csv(os.path.join(OUT2_DIR, f"{sp}_family_enrichment.tsv"),
                              sep="\t", index=False)
            if sf_df is not None and len(sf_df):
                sf_df.to_csv(os.path.join(OUT2_DIR, f"{sp}_superfam_enrichment.tsv"),
                             sep="\t", index=False)
                for r in sf_df.itertuples(index=False):
                    sf_all.setdefault(r.superfam, {})[sp] = r.log2OR_compacted
            if n_done % 25 == 0:
                print(f"  [{n_done}/{len(jobs)}]")

    # ---------- composition ----------
    comp_df = pd.DataFrame(comps).T.fillna(0.0)
    comp_df.index.name = "species"
    comp_df["group"] = [s.split("_", 1)[0] for s in comp_df.index]
    comp_df = comp_df.sort_values(["group", "Compacted"], ascending=[True, False])
    comp_df.to_csv(os.path.join(OUT_ROOT, "composition_2state.tsv"), sep="\t")
    print(f"[COMP] {len(comp_df)} species")

    # stacked plot
    fig, ax = plt.subplots(figsize=(max(8, len(comp_df) * 0.04), 4.0))
    bottom = np.zeros(len(comp_df))
    for st in STATE_ORDER:
        ax.bar(np.arange(len(comp_df)), comp_df[st].values,
               bottom=bottom, color=STATE_COLOR[st], width=1.0,
               edgecolor="none", label=st)
        bottom += comp_df[st].values
    # group separators
    g = comp_df["group"].values
    last = g[0]; xs = []
    for i, gi in enumerate(g):
        if gi != last:
            ax.axvline(i - 0.5, color="black", lw=0.5)
            xs.append((i, gi)); last = gi
    ax.set_xlim(-0.5, len(comp_df) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel("Fraction (Open + Compacted = 1)")
    ax.set_title(f"Two-state chromatin composition (216 fungi) — p_het ≥ {P_HET_THR:.3f}")
    ax.legend(loc="upper right", framealpha=0.9, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "composition_2state.pdf"))
    plt.close(fig)

    # per-group composition
    grp = comp_df.groupby("group")[STATE_ORDER].mean().sort_values("Compacted", ascending=False)
    grp.to_csv(os.path.join(OUT_ROOT, "composition_2state_per_group.tsv"), sep="\t")
    fig, ax = plt.subplots(figsize=(7, 0.35 * len(grp) + 1))
    bottom = np.zeros(len(grp))
    yy = np.arange(len(grp))
    for st in STATE_ORDER:
        ax.barh(yy, grp[st].values, left=bottom, color=STATE_COLOR[st], label=st)
        bottom += grp[st].values
    ax.set_yticks(yy); ax.set_yticklabels(grp.index)
    ax.set_xlim(0, 1); ax.invert_yaxis()
    ax.set_xlabel("Mean fraction (informative bp)")
    ax.set_title("Mean composition per taxonomic group")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "composition_2state_per_group.pdf"))
    plt.close(fig)

    # ---------- superfam heatmap ----------
    sf_mat = pd.DataFrame(sf_all).T   # rows=superfam, cols=species
    sf_mat = sf_mat.reindex(columns=comp_df.index)
    # filter superfams with signal in >=20% species
    keep = sf_mat.notna().sum(axis=1) >= max(20, sf_mat.shape[1] // 5)
    sf_mat = sf_mat.loc[keep]
    sf_mat["__mean"] = sf_mat.mean(axis=1, skipna=True)
    sf_mat = sf_mat.sort_values("__mean", ascending=False).drop(columns="__mean")
    sf_mat.to_csv(os.path.join(OUT_ROOT, "Compacted_superfam_log2OR_matrix.tsv"), sep="\t")

    fig, ax = plt.subplots(figsize=(max(10, sf_mat.shape[1] * 0.05),
                                    max(4, sf_mat.shape[0] * 0.22)))
    vmax = float(np.nanpercentile(np.abs(sf_mat.values), 95)) or 1.0
    im = ax.imshow(sf_mat.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.set_yticks(np.arange(len(sf_mat)))
    ax.set_yticklabels(sf_mat.index, fontsize=7)
    ax.set_xticks([])
    # group separators
    grp_arr = np.array([s.split("_", 1)[0] for s in sf_mat.columns])
    last = grp_arr[0]
    xticks_pos = []; xticks_lab = []; start = 0
    for i, g in enumerate(np.append(grp_arr, ["__END__"])):
        if g != last:
            ax.axvline(i - 0.5, color="black", lw=0.4)
            xticks_pos.append((start + i - 1) / 2); xticks_lab.append(last)
            last = g; start = i
    ax.set_xticks(xticks_pos); ax.set_xticklabels(xticks_lab, fontsize=8, rotation=0)
    ax.set_xlabel("Species (grouped by phylum-level prefix)")
    ax.set_title(f"log2 OR of TE superfamily in Compacted state (p_het ≥ {P_HET_THR:.3f})")
    plt.colorbar(im, ax=ax, shrink=0.6, label="log2 OR")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_ROOT, "Compacted_superfam_log2OR_heatmap.pdf"))
    plt.close(fig)
    print(f"[SF]   matrix {sf_mat.shape}")

    # ---------- top families per group ----------
    fam_files = sorted(glob.glob(os.path.join(OUT2_DIR, "*_family_enrichment.tsv")))
    rows = []
    for f in fam_files:
        sp = os.path.basename(f).replace("_family_enrichment.tsv", "")
        grp = sp.split("_", 1)[0]
        try:
            d = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if d.empty: continue
        d = d[(d["bp_total"] >= 1000) & (d["log2OR_compacted"] > 0.5)]
        d["species"] = sp; d["group"] = grp
        rows.append(d)
    if not rows:
        print("[FAM] no enriched families"); return
    fam_all = pd.concat(rows, ignore_index=True)
    fam_all.to_csv(os.path.join(OUT_ROOT, "Compacted_families_all.tsv"), sep="\t", index=False)

    # top10 most-enriched-on-average per group (by superfam)
    summary = (fam_all.groupby(["group", "superfam"])
               .agg(n_fam=("family", "nunique"),
                    n_sp =("species", "nunique"),
                    mean_log2OR=("log2OR_compacted", "mean"),
                    median_frac_compacted=("frac_compacted", "median"),
                    total_bp=("bp_total", "sum"))
               .reset_index()
               .sort_values(["group", "mean_log2OR"], ascending=[True, False]))
    summary.to_csv(os.path.join(OUT_ROOT, "Compacted_top_families_per_group.tsv"),
                   sep="\t", index=False)

    # plot: per group, top 8 superfams by mean log2OR (only superfams with >=3 fams)
    groups = sorted(summary["group"].unique())
    fig, axes = plt.subplots(nrows=(len(groups) + 3) // 4, ncols=4,
                             figsize=(16, 2.5 * ((len(groups) + 3) // 4)),
                             squeeze=False)
    for ax, g in zip(axes.ravel(), groups):
        d = summary[(summary["group"] == g) & (summary["n_fam"] >= 2)].head(8)
        if d.empty:
            ax.set_axis_off(); ax.set_title(g, fontsize=9); continue
        ax.barh(np.arange(len(d)), d["mean_log2OR"].values,
                color="#E6550D", edgecolor="black", linewidth=0.3)
        ax.set_yticks(np.arange(len(d)))
        ax.set_yticklabels([f"{s}  (n={n})" for s, n in
                            zip(d["superfam"], d["n_fam"])], fontsize=7)
        ax.invert_yaxis()
        ax.set_title(f"{g}  ({d['n_sp'].iloc[0]} sp)", fontsize=9)
        ax.set_xlabel("mean log2 OR (Compacted)", fontsize=7)
        ax.tick_params(axis="x", labelsize=7)
    for ax in axes.ravel()[len(groups):]:
        ax.set_axis_off()
    fig.suptitle("TE superfamilies most enriched in Compacted state (mean over enriched families per group)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(OUT_ROOT, "Compacted_top_families_per_group.pdf"))
    plt.close(fig)
    print("[DONE]")


if __name__ == "__main__":
    main()
