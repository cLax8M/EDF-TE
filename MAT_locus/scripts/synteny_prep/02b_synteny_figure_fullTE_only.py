#!/usr/bin/env python3
"""
Synteny figure: sex locus structure + TEs across representative Mucorales
Groups: Rhizopus | Mucor | Phycomyces | Other Mucorales
"""

import pandas as pd
import numpy  as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot     as plt
import matplotlib.patches    as mpatches
import matplotlib.gridspec   as gridspec
from matplotlib.patches import FancyArrow
from pathlib import Path

WORKDIR  = Path("/data/carlos/EDF-TE/TE_analysis/sex_locus")
SYNDIR   = WORKDIR / "synteny_figure"

# ─── Window definitions (scaffold, win_start_1based, win_end_1based) ─────────
WINDOWS = {
    "Blatri":     ("scaffold_61",                        1,     159599),
    "Rhisto1":    ("scaffold_35",                   151440,     251440),
    "Rhimi59":    ("scaffold_1",                   310000,     420000),
    "Rhior3":     ("supercontig_3.1",             3260000,    3360000),
    "Mucci3":     ("scaffold_1",                   2163094,    2263094),
    "mu402":      ("scaffold_1",                   2175952,    2275952),
    "Mucpir1":    ("scaffold_29",                        1,     368276),
    "Mucmuc1":    ("scaffold_63",                        1,     213167),
    "Mycind1":    ("scaffold_10",                   581915,     681915),
    "PhyblU21":   ("scaffold_17",                   282700,     382700),
    "Phybl4":     ("scaffold_10",                   676371,     776371),
    "PnitS607":   ("scaffold_54",                        1,      92868),
    "PnitSJTU":   ("JAHRBC010000002.1",            6796920,    6896920),
    "Hesve2":     ("scaffold_28",                        1,     341884),
    "Thaele1":    ("scaffold_15",                   460694,     560694),
    "Ellano2":    ("scaffold_13",                   512635,     612635),
    "Chocucu1":   ("scaffold_21",                   200045,     300045),
}

# What portion of each scaffold we display (1-based, inclusive)
VIEW = {
    "Blatri":    (40000,   90000),   # 50 kb window centred on genes
    "Rhisto1":   (151440,  251440),
    "Rhimi59":   (345000,  400000),
    "Rhior3":    (3280000, 3340000),
    "Mucci3":    (2163094, 2263094),
    "mu402":     (2175952, 2275952),
    "Mucpir1":   (155000,  210000),
    "Mucmuc1":   (88000,   140000),
    "Mycind1":   (581915,  681915),
    "PhyblU21":  (282700,  382700),
    "Phybl4":    (676371,  776371),
    "PnitS607":  (1,        92868),
    "PnitSJTU":  (6796920, 6896920),
    "Hesve2":    (88000,   160000),
    "Thaele1":   (460694,  560694),
    "Ellano2":   (512635,  612635),
    "Chocucu1":  (200045,  300045),
}

GROUPS = [
    ("Rhizopus",        ["Blatri",  "Rhisto1",  "Rhimi59",  "Rhior3"]),
    ("Mucor",           ["Mucci3",  "mu402",    "Mucpir1", "Mucmuc1", "Mycind1"]),
    ("Phycomyces",      ["PhyblU21","Phybl4",   "PnitS607","PnitSJTU"]),
    ("Other Mucorales", ["Hesve2",  "Thaele1",  "Ellano2", "Chocucu1"]),
]

TRACK_ORDER = [sp for (_, members) in GROUPS for sp in members]

SPECIES_NAMES = {
    "Blatri":    "Backusella lamprospora",
    "Rhisto1":   "Rhizopus stolonifer",
    "Rhimi59":   "Rh. microsporus",
    "Rhior3":    "Rh. oryzae",
    "Mucci3":    "Mucor circinelloides",
    "mu402":     "Mu. circinelloides",
    "Mucpir1":   "Mu. piriformis",
    "Mucmuc1":   "Mu. mucedo",
    "Mycind1":   "Mycotypha indica",
    "PhyblU21":  "Phycomyces blakesleeanus",
    "Phybl4":    "Ph. blakesleeanus",
    "PnitS607":  "Ph. nitens",
    "PnitSJTU":  "Ph. nitens",
    "Hesve2":    "Hesseltinella vesiculosa",
    "Thaele1":   "Thamnidium elegans",
    "Ellano2":   "Ellisomyces anomalus",
    "Chocucu1":  "Choanephora cucurbitarum",
}

MT_LABELS = {
    "Blatri":   "sexM?", "Rhisto1": "sexM?", "Rhimi59":  "sexP",  "Rhior3":   "sexP",
    "Mucci3":   "sexP",  "mu402":   "sexP",  "Mucpir1":  "sexP",
    "Mucmuc1":  "sexP",  "Mycind1": "sexP",
    "PhyblU21": "sexP",  "Phybl4":  "sexM",
    "PnitS607": "sexP",  "PnitSJTU":"sexM",
    "Hesve2":   "sexM",  "Thaele1": "sexP",  "Ellano2": "sexP", "Chocucu1": "sexP",
}

GENES = {
    "Blatri":   [(47923,  49257,"-","#27AE60","tptA"),(49523, 50149,"+","#C0392B","sexM"),(50344,51654,"+","#2980B9","rnhA"),(57551,59080,"-","#F1C40F","sagA")],
    "Rhisto1":  [(195713,197264,"+","#F1C40F","sagA"),(197696,199942,"-","#2980B9","rnhA"),(201440,202009,"-","#C0392B","sexM"),(204333,204629,"+","#27AE60","arbA?")],
    "Rhimi59":  [(356782,359109,"+","#784212","glrA"),(359294,360854,"+","#F1C40F","sagA"),(365806,367182,"-","#2980B9","rnhA"),(367788,368324,"-","#E8671A","sexP"),(369876,372068,"+","#95A5A6","unk")],
    "Rhior3":   [(3301855,3302544,"+","#27AE60","tptA"),(3309558,3310079,"+","#E8671A","sexP"),(3313526,3314155,"+","#2980B9","rnhA"),(3315110,3316651,"-","#F1C40F","sagA")],
    "Mucci3":   [(2209870,2211150,"+","#27AE60","tptA"),(2212720,2213468,"+","#E8671A","sexP"),(2214300,2215700,"+","#2980B9","rnhA"),(2222333,2223840,"-","#F1C40F","sagA")],
    "mu402":    [(2227187,2230054,"+","#27AE60","tptA"),(2225606,2226298,"+","#E8671A","sexP"),(2223386,2224864,"+","#2980B9","rnhA"),(2235221,2236728,"-","#F1C40F","sagA")],
    "Mucpir1":  [(170087,171592,"+","#F1C40F","sagA"),(175919,177838,"-","#2980B9","rnhA"),(178751,179722,"-","#E8671A","sexP")],
    "Mucmuc1":  [(97178,98663,"+","#F1C40F","sagA"),(101596,102783,"-","#2980B9","rnhA"),(105841,107225,"-","#E8671A","sexP")],
    "Mycind1":  [(625203,626743,"+","#F1C40F","sagA"),(627534,630347,"-","#2980B9","rnhA"),(631915,632844,"-","#E8671A","sexP")],
    "PhyblU21": [(323739,324380,"+","#27AE60","tptA"),(332315,333085,"+","#E8671A","sexP"),(334587,336122,"+","#2980B9","rnhA")],
    "Phybl4":   [(721220,722100,"+","#27AE60","tptA"),(726130,726613,"+","#C0392B","sexM"),(728917,730200,"+","#2980B9","rnhA")],
    "PnitS607": [(29836,31386,"-","#27AE60","tptA"),(42507,43229,"-","#E8671A","sexP"),(46476,49172,"+","#2980B9","rnhA")],
    "PnitSJTU": [(6833122,6834672,"-","#2980B9","rnhA"),(6846720,6847121,"-","#C0392B","sexM"),(6848516,6852002,"+","#27AE60","tptA")],
    "Hesve2":   [(95768,97024,"-","#F1C40F","sagA"),(112706,114064,"-","#2980B9","rnhA"),(117425,118105,"-","#C0392B","sexM"),(119276,119965,"+","#27AE60","tptA")],
    "Thaele1":  [(484077,485522,"+","#F1C40F","sagA"),(509173,510327,"-","#27AE60","tptA"),(510694,511740,"+","#E8671A","sexP"),(513474,515654,"+","#2980B9","rnhA")],
    "Ellano2":  [(552119,553619,"+","#F1C40F","sagA"),(563849,565324,"+","#27AE60","tptA"),(562635,563555,"-","#E8671A","sexP"),(560901,562304,"-","#2980B9","rnhA")],
    "Chocucu1": [(248504,249178,"-","#27AE60","tptA"),(250045,250971,"+","#E8671A","sexP"),(251143,252453,"+","#2980B9","rnhA"),(258177,259694,"-","#F1C40F","sagA")],
}

IDIO = {
    "Blatri":    (49257,  50344),
    "Rhisto1":   (199942, 204333),
    "Rhimi59":   (367182, 368847),
    "Rhior3":    (3302544, 3313526),
    "Mucci3":    (2211150,2214300),
    "mu402":     (2224864,2227187),
    "Mucpir1":   (175919, 179722),
    "Mucmuc1":   (101596, 107225),
    "Mycind1":   (627534, 632844),
    "PhyblU21":  (324380, 334587),
    "Phybl4":    (722100, 728917),
    "PnitS607":  (31386,  46476),
    "PnitSJTU":  (6834672,6848516),
    "Hesve2":    (112706, 119965),
    "Thaele1":   (509173, 515654),
    "Ellano2":   (560901, 565324),
    "Chocucu1":  (248504, 252453),
}

SYNTENY_PAIRS = [
    ("Blatri",   "Rhisto1",   "Blatri_vs_Rhisto1.paf",    "blastn"),
    ("Rhisto1",  "Rhimi59",   "Rhisto1_vs_Rhimi59.paf",   "blastn"),
    ("Rhimi59",  "Rhior3",    "Rhimi59_vs_Rhior3.paf",    "blastn"),
    ("Mucci3",   "mu402",     "Mucci3_vs_mu402.paf",       "minimap2"),
    ("mu402",    "Mucpir1",   "mu402_vs_Mucpir1.paf",      "blastn"),
    ("Mucpir1",  "Mucmuc1",   "Mucpir1_vs_Mucmuc1.paf",   "blastn"),
    ("Mucmuc1",  "Mycind1",   "Mucmuc1_vs_Mycind1.paf",   "blastn"),
    ("PhyblU21", "Phybl4",    "Phybl4_vs_PhyblU21.paf",   "minimap2"),
    ("PnitS607", "PnitSJTU",  "PnitS607_vs_PnitSJTU.paf", "minimap2"),
    # Hesve2 vs Thaele1: only repetitive/microsatellite hits, no real synteny detected
    ("Thaele1",  "Ellano2",   "Thaele1_vs_Ellano2.paf",   "blastn"),
    ("Ellano2",  "Chocucu1",  "Ellano2_vs_Chocucu1.paf",  "blastn"),
]

TE_COLORS = {
    "LTR":     "#8E44AD",
    "LINE":    "#1A5276",
    "DNA":     "#B7950B",
    "RC":      "#0E6655",
    "PLE":     "#117A65",
    "SINE":    "#CB4335",
    "Unknown": "#808B96",
    "Simple_repeat":  "#D5D8DC",
    "Low_complexity": "#D5D8DC",
}

def te_main_class(name):
    if "#" in name:
        cls = name.split("#")[-1].strip().split("/")[0]
        return cls
    return "Unknown"

def te_color(cls):
    return TE_COLORS.get(cls, "#BDC3C7")

def load_TE(sp, min_bp=50):
    scaf, win_s, win_e = WINDOWS[sp]
    is_windowed = (win_s > 1)
    tsv = SYNDIR / "rmask" / f"{sp}_TE.tsv"
    if not tsv.exists():
        return pd.DataFrame(columns=["start","end","cls","strand","full"])
    cols = ["sseqid","sstart","send","pident","length","qseqid","qlen","evalue","bitscore"]
    df = pd.read_csv(tsv, sep="\t", header=None, names=cols)
    if len(df) == 0:
        return pd.DataFrame(columns=["start","end","cls","strand","full"])
    df = df[df["length"] >= min_bp].copy()
    df["strand"] = df.apply(lambda r: "+" if r["sstart"] <= r["send"] else "-", axis=1)
    df["s0"] = df[["sstart","send"]].min(axis=1)
    df["s1"] = df[["sstart","send"]].max(axis=1)
    if is_windowed:
        df["start"] = df["s0"] + win_s - 1
        df["end"]   = df["s1"] + win_s - 1
    else:
        df["start"] = df["s0"]
        df["end"]   = df["s1"]
    df["cls"]  = df["qseqid"].apply(te_main_class)
    df["full"] = (df["length"] / df["qlen"] >= 0.65) & (df["pident"] >= 75)
    # FULL-TE-ONLY version: keep only full-length elements
    df = df[df["full"]].copy()
    return df[["start","end","cls","strand","full"]].copy()

def load_synteny(q_sp, t_sp, paf_file, fmt):
    if paf_file is None:
        return []
    fpath = SYNDIR / paf_file
    if not fpath.exists():
        return []
    q_scaf, q_ws, q_we = WINDOWS[q_sp]
    t_scaf, t_ws, t_we = WINDOWS[t_sp]
    q_windowed = (q_ws > 1)
    t_windowed = (t_ws > 1)
    blocks = []
    if fmt == "minimap2":
        cols = ["qname","qlen","qstart","qend","strand","tname","tlen","tstart","tend","nmatch","alen","mapq"]
        df = pd.read_csv(fpath, sep="\t", header=None, usecols=range(12), names=cols)
        for _, r in df.iterrows():
            if r["alen"] < 500:
                continue
            # PAF 0-based half-open -> 1-based scaffold coord
            qs = int(r["qstart"]) + 1 + (q_ws - 1)
            qe = int(r["qend"])       + (q_ws - 1)
            ts = int(r["tstart"]) + 1 + (t_ws - 1)
            te = int(r["tend"])       + (t_ws - 1)
            blocks.append(dict(q_start=qs, q_end=qe, t_start=ts, t_end=te,
                               strand=r["strand"], alen=int(r["alen"])))
    elif fmt == "blastn":
        cols = ["qid","sid","pident","length","qstart","qend","sstart","send","evalue","bitscore"]
        df = pd.read_csv(fpath, sep="\t", header=None, names=cols)
        for _, r in df.iterrows():
            if int(r["length"]) < 200:
                continue
            qrs = int(r["qstart"]); qre = int(r["qend"])
            srs = int(r["sstart"]); sre = int(r["send"])
            if q_windowed:
                qs = qrs + q_ws - 1; qe = qre + q_ws - 1
            else:
                qs = qrs; qe = qre
            if t_windowed:
                ts = min(srs,sre) + t_ws - 1; te = max(srs,sre) + t_ws - 1
            else:
                ts = min(srs,sre); te = max(srs,sre)
            strand = "-" if srs > sre else "+"
            blocks.append(dict(q_start=qs, q_end=qe, t_start=ts, t_end=te,
                               strand=strand, alen=int(r["length"])))
    return blocks

TRACK_H = 0.065

def draw_track(ax, sp):
    scaf, win_s, win_e = WINDOWS[sp]
    vs, ve = VIEW[sp]
    vlen = ve - vs

    def x(pos):
        return (pos - vs) / vlen

    yc = 0.5
    ax.axhline(yc, color="#444444", lw=0.9, zorder=1)

    # idiomorphic highlight
    idio = IDIO.get(sp)
    if idio:
        ix0 = max(0.0, x(idio[0])); ix1 = min(1.0, x(idio[1]))
        if ix0 < ix1:
            ax.add_patch(mpatches.Rectangle(
                (ix0, yc-0.28), ix1-ix0, 0.56,
                color="#FDEBD0", alpha=0.50, zorder=0, linewidth=0))
            for xi in (ix0, ix1):
                ax.plot([xi,xi],[yc-0.30,yc+0.30], color="#E59866", lw=0.8, ls="--", zorder=2)

    # TE blocks / arrows
    te_df = load_TE(sp)
    for _, tr in te_df.iterrows():
        tx0 = x(tr["start"]); tx1 = x(tr["end"])
        if tx1 <= 0 or tx0 >= 1:
            continue
        tx0 = max(0.0, tx0); tx1 = min(1.0, tx1)
        tw = tx1 - tx0
        if tw < 2e-4:
            continue
        tcol = te_color(tr["cls"])
        if tr["full"]:
            # Full-length TE -> arrow like genes but smaller
            ast = tx0 if tr["strand"] == "+" else tx1
            dx  = tw  if tr["strand"] == "+" else -tw
            hl  = min(abs(dx) * 0.30, 0.009)
            ax.add_patch(FancyArrow(
                ast, yc, dx, 0,
                width=TRACK_H*1.55, head_width=TRACK_H*2.6, head_length=hl,
                length_includes_head=True, facecolor=tcol, edgecolor="white",
                zorder=4, linewidth=0.3, overhang=0.05))
        else:
            ax.add_patch(mpatches.Rectangle(
                (tx0, yc-TRACK_H*0.75), tw, TRACK_H*1.5,
                color=tcol, alpha=0.75, zorder=3, linewidth=0))

    # gene arrows
    for gs, ge, gstrand, gcolor, glabel in GENES.get(sp, []):
        gx0 = x(gs); gx1 = x(ge)
        if gx1 <= 0 or gx0 >= 1:
            continue
        gx0c = max(0.0, gx0); gx1c = min(1.0, gx1)
        gw = gx1c - gx0c
        if gw < 1e-5:
            continue
        ast = gx0c if gstrand == "+" else gx1c
        dx  = gw   if gstrand == "+" else -gw
        hl  = min(abs(dx)*0.35, 0.012)
        ax.add_patch(FancyArrow(
            ast, yc, dx, 0,
            width=TRACK_H*2.0, head_width=TRACK_H*3.0, head_length=hl,
            length_includes_head=True, color=gcolor, zorder=6, linewidth=0, overhang=0.1))
        ax.text((gx0c+gx1c)/2, yc+TRACK_H*3.9, glabel,
                ha="center", va="bottom", fontsize=5.8, color=gcolor, zorder=7,
                clip_on=True, fontstyle="italic")

    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")

    mt = MT_LABELS.get(sp,"?")
    mt_col = "#C0392B" if "sexM" in mt else ("#E8671A" if "sexP" in mt else "#777777")
    sp_name = SPECIES_NAMES.get(sp, sp)
    # Use non-italic + italic rendering; avoid LaTeX space issues
    parts = sp_name.split(" ", 1)
    genus   = parts[0]
    epithet = parts[1] if len(parts) > 1 else ""
    label_str = f"$\\it{{{genus}}}$" + (f" $\\it{{{epithet}}}$" if epithet else "")
    ax.text(-0.01, yc+0.12, label_str,
            ha="right", va="center", fontsize=6.2, transform=ax.transAxes, color="#222222")
    ax.text(-0.01, yc-0.22,
            f"{sp}  {mt}",
            ha="right", va="top", fontsize=5.2, transform=ax.transAxes, color=mt_col)


def draw_synteny_band(fig, ax_q, ax_t, q_sp, t_sp, paf_file, fmt, min_alen=200):
    blocks = load_synteny(q_sp, t_sp, paf_file, fmt)
    if not blocks:
        rev = load_synteny(t_sp, q_sp, paf_file, fmt)
        blocks = [dict(q_start=b["t_start"], q_end=b["t_end"],
                       t_start=b["q_start"], t_end=b["q_end"],
                       strand=b["strand"], alen=b["alen"]) for b in rev]
    if not blocks:
        return

    # Build a dedicated connector axes in the gap between the two tracks
    pos_q = ax_q.get_position()   # upper track
    pos_t = ax_t.get_position()   # lower track
    gap_y0 = pos_t.y1             # top edge of lower axes
    gap_y1 = pos_q.y0             # bottom edge of upper axes
    if gap_y1 <= gap_y0 + 1e-4:
        return
    conn_ax = fig.add_axes([pos_q.x0, gap_y0,
                             pos_q.x1 - pos_q.x0, gap_y1 - gap_y0])
    conn_ax.set_xlim(0, 1); conn_ax.set_ylim(0, 1)
    conn_ax.axis("off"); conn_ax.patch.set_visible(False)

    q_vs, q_ve = VIEW[q_sp]; q_vlen = q_ve - q_vs
    t_vs, t_ve = VIEW[t_sp]; t_vlen = t_ve - t_vs

    for b in blocks:
        if b["alen"] < min_alen:
            continue
        qx0 = (b["q_start"] - q_vs) / q_vlen
        qx1 = (b["q_end"]   - q_vs) / q_vlen
        tx0 = (b["t_start"] - t_vs) / t_vlen
        tx1 = (b["t_end"]   - t_vs) / t_vlen
        qx0, qx1 = sorted([max(0.0, min(1.0, qx0)), max(0.0, min(1.0, qx1))])
        tx0, tx1 = sorted([max(0.0, min(1.0, tx0)), max(0.0, min(1.0, tx1))])
        if qx1 - qx0 < 1e-4 or tx1 - tx0 < 1e-4:
            continue
        col_f = "#5DADE2" if b["strand"] == "+" else "#EC7063"
        col_e = "#2471A3" if b["strand"] == "+" else "#C0392B"
        poly = mpatches.Polygon(
            [(qx0, 1.0), (qx1, 1.0), (tx1, 0.0), (tx0, 0.0)],
            closed=True, facecolor=col_f, edgecolor=col_e,
            alpha=0.60, linewidth=0.6, zorder=5)
        conn_ax.add_patch(poly)


# ─── Build figure ─────────────────────────────────────────────────────────────
heights = []
for gi, (gname, members) in enumerate(GROUPS):
    for _ in members:
        heights.append(3)
    if gi < len(GROUPS)-1:
        heights.append(1.4)

fig_h = 0.75 * len(heights) + 1.2
fig = plt.figure(figsize=(18, fig_h))
fig.patch.set_facecolor("white")

gs_main = gridspec.GridSpec(len(heights), 1, figure=fig,
                             height_ratios=heights, hspace=0.55,
                             left=0.24, right=0.90, top=0.96, bottom=0.06)

axes = {}
row = 0
for gi, (gname, members) in enumerate(GROUPS):
    for sp in members:
        axes[sp] = fig.add_subplot(gs_main[row])
        row += 1
    if gi < len(GROUPS)-1:
        row += 1

for sp in TRACK_ORDER:
    draw_track(axes[sp], sp)

for q_sp, t_sp, paf_file, fmt in SYNTENY_PAIRS:
    draw_synteny_band(fig, axes[q_sp], axes[t_sp], q_sp, t_sp, paf_file, fmt)

# Scale bar on last track
last_sp = TRACK_ORDER[-1]
ax_last = axes[last_sp]
vs, ve = VIEW[last_sp]
bar_kb  = 10000
bar_frac = bar_kb / (ve - vs)
bx1 = 0.92; bx0 = bx1 - bar_frac
ax_last.annotate("", xy=(bx1,0.05), xytext=(bx0,0.05),
                 xycoords="axes fraction",
                 arrowprops=dict(arrowstyle="-", color="black", lw=1.5))
ax_last.text((bx0+bx1)/2, -0.15, "10 kb",
             ha="center", va="top", fontsize=7, transform=ax_last.transAxes)

# Group labels right side
for gname, members in GROUPS:
    ax0 = axes[members[0]]; ax1 = axes[members[-1]]
    p0 = ax0.get_position(); p1 = ax1.get_position()
    xb = p0.x1 + 0.012
    yt = p0.y1; yb = p1.y0; ym = (yt+yb)/2
    ax0.annotate("", xy=(xb,yb), xytext=(xb,yt),
                 xycoords="figure fraction", textcoords="figure fraction",
                 arrowprops=dict(arrowstyle="-", color="#444444", lw=1.2))
    for yy in (yt,yb):
        ax0.annotate("", xy=(xb+0.004,yy), xytext=(xb,yy),
                     xycoords="figure fraction", textcoords="figure fraction",
                     arrowprops=dict(arrowstyle="-", color="#444444", lw=1.0))
    fig.text(xb+0.013, ym, gname,
             ha="left", va="center", fontsize=8.5, fontweight="bold", color="#333333")

# Legend
legend_elems = [
    mpatches.Patch(color="#E8671A",label="sexP gene"),
    mpatches.Patch(color="#C0392B",label="sexM gene"),
    mpatches.Patch(color="#27AE60",label="tptA"),
    mpatches.Patch(color="#2980B9",label="rnhA"),
    mpatches.Patch(color="#F1C40F",label="sagA"),
    mpatches.Patch(color="#784212",label="glrA"),
    mpatches.Patch(color="#95A5A6",label="unk (Rhimi59)"),
    mpatches.Patch(color="#FDEBD0",alpha=0.9,label="Idiomorphic region"),
    mpatches.FancyArrow(0, 0, 0.1, 0, width=0.3, head_width=0.5, head_length=0.02,
                        color="#8E44AD", label="Full TE (LTR)"),
    mpatches.Patch(color="#8E44AD",alpha=0.75,label="LTR remnant"),
    mpatches.FancyArrow(0, 0, 0.1, 0, width=0.3, head_width=0.5, head_length=0.02,
                        color="#1A5276", label="Full TE (LINE)"),
    mpatches.Patch(color="#1A5276",alpha=0.75,label="LINE remnant"),
    mpatches.Patch(color="#B7950B",alpha=0.75,label="DNA TE"),
    mpatches.Patch(color="#0E6655",alpha=0.75,label="RC/Helitron"),
    mpatches.Patch(color="#808B96",alpha=0.75,label="Unknown TE"),
    mpatches.Patch(color="#5DADE2",alpha=0.75,label="Synteny (same strand)"),
    mpatches.Patch(color="#EC7063",alpha=0.75,label="Synteny (inverted)"),
]
fig.legend(handles=legend_elems, loc="lower center", ncol=6,
           fontsize=6.5, frameon=True, framealpha=0.95, edgecolor="lightgray",
           bbox_to_anchor=(0.57, 0.0), bbox_transform=fig.transFigure)

fig.suptitle("Sex determination locus synteny across Mucorales",
             fontsize=11, fontweight="bold", y=0.99)

out_pdf = WORKDIR / "sex_locus_synteny_figure_fullTE_only.pdf"
out_png = WORKDIR / "sex_locus_synteny_figure_fullTE_only.png"
fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
plt.close()
