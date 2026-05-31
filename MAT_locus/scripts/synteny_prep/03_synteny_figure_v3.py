#!/usr/bin/env python3
"""
Sex locus synteny figure v3 — 60 kb windows centred on sexP/sexM
Groups: Mucor | Phycomyces | Rhizopus
All annotated genes shown as grey arrows; known sex-locus genes coloured.
Synteny from pairwise BLASTN of extracted 60-kb windows.
"""

import os, re
import pandas as pd
import numpy  as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot   as plt
import matplotlib.patches  as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrow
from pathlib import Path
from collections import defaultdict

WORKDIR = Path("/data/carlos/EDF-TE/TE_analysis/sex_locus")
SYNDIR  = WORKDIR / "synteny_figure"           # PAF + rmask
WINDIR  = SYNDIR / "windows60"                  # extracted 60kb fastas + new BLASTN PAFs
BGDIR   = SYNDIR / "bg_genes"                   # TBLASTN-derived background gene TSVs
GFF_DIR = Path("/legserv/Temp/Carlos/EDF-TE/data/genes_annot_gff")

# ── 60-kb VIEW windows (scaffold_1-based start, end) ─────────────────────────
VIEW = {
    # Rhizopus group
    "Rhisto1":   ("scaffold_35",      171724, 231724),
    "Rhimi59":   ("scaffold_1",       338056, 398056),
    "ATCC52813": ("scaffold_2",      3186275,3246275),
    "ATCC52814": ("scaffold_1",      3983605,4043605),
    "Rhior3":    ("supercontig_3.1", 3279818,3339818),
    # Mucor group
    "Mucci3":    ("scaffold_1",      2183094,2243094),
    "Muccirc1":  ("KE_123975",         43879, 103879),
    "PS15m":     ("JBOZQX010000005.1",1669161,1729161),
    "PS15p":     ("JBNRQL010000006.1",1709974,1769974),
    "Mucmuc1":   ("scaffold_63",       76533, 136533),
    "Mucpir1":   ("scaffold_29",      149236, 209236),
    # Phycomyces group
    "PhyblU21":  ("scaffold_17",      302700, 362700),
    "Phybl4":    ("scaffold_10",      696371, 756371),
    "PnitS607":  ("scaffold_54",       12868,  72868),
    "PnitSJTU":  ("JAHRBC010000002.1",6816920,6876920),
    # Fenlin pair (CB54 reference vs JGI Fenlin1, 2026-04-28)
    "Fenlin":    ("scaffold_5",      2377539,2437539),
    "Fenlin1":   ("scaffold_45",      286592, 346592),
}

# For TE lookup we keep the original wider WINDOWS for species that have RM data
WINDOWS = {
    "Blatri":    ("scaffold_61",          1,    159599),
    "Rhisto1":   ("scaffold_35",     151440,    251440),
    "Rhimi59":   ("scaffold_1",      310000,    420000),
    "Rhior3":    ("supercontig_3.1",3260000,   3360000),
    "Mucci3":    ("scaffold_1",     2163094,   2263094),
    "Mucmuc1":   ("scaffold_63",          1,    213167),
    "Mucpir1":   ("scaffold_29",          1,    368276),
    "PhyblU21":  ("scaffold_17",     282700,    382700),
    "Phybl4":    ("scaffold_10",     676371,    776371),
    "PnitS607":  ("scaffold_54",          1,     92868),
    "PnitSJTU":  ("JAHRBC010000002.1",6796920,  6896920),
    # new species
    "Muccirc1":  ("KE_123975", 1, 158905),   # full scaffold; TE TSV has absolute coords
    "ATCC52813": VIEW["ATCC52813"],
    "ATCC52814": VIEW["ATCC52814"],
    "PS15m":     VIEW["PS15m"],
    "PS15p":     VIEW["PS15p"],
    # Fenlin pair (TE annotations parsed from Earl Grey filteredRepeats.gff)
    "Fenlin":    VIEW["Fenlin"],
    "Fenlin1":   VIEW["Fenlin1"],
}

GROUPS = [
    ("Rhizopus",   ["Rhisto1","Rhimi59","ATCC52813","ATCC52814","Rhior3"]),
    ("Mucor",      ["Mucci3","Muccirc1","PS15m","PS15p","Mucmuc1","Mucpir1"]),
    ("Phycomyces", ["PhyblU21","Phybl4","PnitS607","PnitSJTU"]),
    ("Fenlin",     ["Fenlin","Fenlin1"]),
]

TRACK_ORDER = [sp for (_, members) in GROUPS for sp in members]

SPECIES_NAMES = {
    "Rhisto1":   "Rhizopus stolonifer",
    "Rhimi59":   "Rh. microsporus",
    "ATCC52813": "Rh. microsporus ATCC52813",
    "ATCC52814": "Rh. microsporus ATCC52814",
    "Rhior3":    "Rh. oryzae",
    "Mucci3":    "Mucor circinelloides CBS",
    "Muccirc1":  "Mu. circinelloides fcirc",
    "PS15m":     "Mu. circinelloides PS15(−)",
    "PS15p":     "Mu. circinelloides PS15(+)",
    "Mucmuc1":   "Mu. mucedo",
    "Mucpir1":   "Mu. piriformis",
    "PhyblU21":  "Phycomyces blakesleeanus",
    "Phybl4":    "Ph. blakesleeanus",
    "PnitS607":  "Ph. nitens S607",
    "PnitSJTU":  "Ph. nitens SJTU",
    "Fenlin":    "Fennellomyces linderi CB54",
    "Fenlin1":   "F. linderi (Fenlin1)",
}

MT_LABELS = {
    "Rhisto1":  "sexM",  "Rhimi59":   "sexP",
    "ATCC52813":"sexM",  "ATCC52814": "sexP",   "Rhior3":   "sexP",
    "Mucci3":   "sexP",  "Muccirc1":  "sexM",
    "PS15m":    "sexM(−)","PS15p":   "sexP(+)", "Mucmuc1":  "sexP",  "Mucpir1":  "sexP",
    "PhyblU21": "sexP",  "Phybl4":   "sexM",    "PnitS607": "sexP",  "PnitSJTU": "sexM",
    "Fenlin":   "sexM",  "Fenlin1":  "sexM",
}

# ── Known sex-locus genes (scaffold coords, 1-based) ─────────────────────────
GENES = {  # (start, end, strand, colour, label)
    "Rhisto1":   [(195713,197264,"+","#F1C40F","sagA"),(197696,199942,"-","#2980B9","rnhA"),
                  (201440,202009,"-","#C0392B","sexM"),(204333,204629,"+","#27AE60","arbA")],
    "Rhimi59":   [(356782,359109,"+","#784212","glrA"),(359294,360854,"+","#F1C40F","sagA"),
                  (365806,367182,"-","#2980B9","rnhA"),(367788,368324,"-","#E8671A","sexP"),
                  (369876,372068,"+","#95A5A6","unk")],
    "ATCC52813": [(3215677,3219700,"-","#2980B9","rnhA")],
    "ATCC52814": [(4009028,4009564,"+","#E8671A","sexP"),(4010180,4014203,"+","#2980B9","rnhA")],
    "Rhior3":    [(3301855,3302544,"+","#27AE60","tptA"),(3309558,3310079,"+","#E8671A","sexP"),
                  (3313526,3314155,"+","#2980B9","rnhA"),(3315110,3316651,"-","#F1C40F","sagA")],
    "Mucci3":    [(2209870,2211150,"+","#27AE60","tptA"),(2212720,2213468,"+","#E8671A","sexP"),
                  (2214300,2215700,"+","#2980B9","rnhA"),(2222333,2223840,"-","#F1C40F","sagA")],
    "Muccirc1":  [(67540,71621,"-","#2980B9","rnhA"),(71964,72932,"-","#C0392B","sexM"),
                  (73156,74602,"+","#27AE60","tptA")],
    "PS15m":     [(1694983,1695705,"-","#2980B9","rnhA"),(1697182,1697886,"+","#C0392B","sexM"),
                  (1699585,1700766,"+","#27AE60","tptA")],
    "PS15p":     [(1737919,1738641,"-","#2980B9","rnhA"),(1739628,1740320,"+","#E8671A","sexP"),
                  (1742380,1743561,"+","#27AE60","tptA")],
    "Mucmuc1":   [(97178,98663,"+","#F1C40F","sagA"),(101596,102783,"-","#2980B9","rnhA"),
                  (105841,107225,"-","#E8671A","sexP")],
    "Mucpir1":   [(170087,171592,"+","#F1C40F","sagA"),(175919,177838,"-","#2980B9","rnhA"),
                  (178751,179722,"-","#E8671A","sexP")],
    "PhyblU21":  [(322900,328244,"-","#2980B9","rnhA"),(332486,333426,"+","#E8671A","sexP"),
                  (334189,336122,"+","#27AE60","tptA")],
    "Phybl4":    [(721220,723470,"-","#27AE60","tptA"),(726130,726613,"+","#C0392B","sexM"),
                  (726831,732097,"+","#2980B9","rnhA")],
    "PnitS607":  [(29836,31386,"-","#27AE60","tptA"),(42507,43229,"-","#E8671A","sexP"),
                  (46476,49172,"+","#2980B9","rnhA")],
    "PnitSJTU":  [(6833122,6834672,"+","#27AE60","tptA"),(6846720,6847121,"-","#C0392B","sexM"),
                  (6848516,6852002,"+","#2980B9","rnhA")],
    # Fenlin pair: sexM + helicase from blastp_proteomes/locus_orthologs_v3.tsv (2026-04-28)
    "Fenlin":    [(2407254,2407823,"+","#C0392B","sexM"),
                  (2413688,2418037,"+","#8E44AD","helicase")],
    "Fenlin1":   [(310225, 314859,"-","#8E44AD","helicase"),
                  (320293, 320898,"-","#C0392B","sexM")],
}

# ── Idiomorphic region (shaded band) ─────────────────────────────────────────
IDIO = {
    "Rhisto1":   (199942, 204333),
    "Rhimi59":   (367182, 369876),
    "ATCC52814": (4009028,4010180),
    "Rhior3":    (3302544,3313526),
    "Mucci3":    (2211150,2214300),
    "Muccirc1":  (71621,  73156),
    "PS15m":     (1695705,1699585),
    "PS15p":     (1738641,1742380),
    "Mucmuc1":   (101596, 107225),
    "Mucpir1":   (175919, 179722),
    "PhyblU21":  (328244, 334189),
    "Phybl4":    (723470, 726831),
    "PnitS607":  (31386,  46476),
    "PnitSJTU":  (6834672,6848516),
    # Fenlin pair (idiomorphic_start/end from sex_locus_master_v3.tsv)
    "Fenlin":    (2402254, 2412823),
    "Fenlin1":   (315293,  325898),
}

# ── Synteny pairs (query → target, filename, fmt) ────────────────────────────
# All new files in WINDIR; old files (fmt="old") use the original PAF format.
# New files are in WINDIR and use window-local coordinates (1-60001).
SYNTENY_PAIRS = [
    # Rhizopus
    ("Rhisto1",   "Rhimi59",   "Rhisto1_vs_Rhimi59_60kb.paf",   "win60"),
    ("Rhimi59",   "ATCC52813", "Rhimi59_vs_ATCC52813.paf",      "win60"),
    ("ATCC52813", "ATCC52814", "ATCC52813_vs_ATCC52814.paf",    "win60"),
    ("ATCC52814", "Rhior3",    "ATCC52814_vs_Rhior3.paf",       "win60"),
    # Mucor
    ("Mucci3",    "Muccirc1",  "Mucci3_vs_Muccirc1.paf",        "win60"),
    ("Muccirc1",  "PS15m",     "Muccirc1_vs_PS15m.paf",         "win60"),
    ("PS15m",     "PS15p",     "PS15m_vs_PS15p.paf",            "win60"),
    ("PS15p",     "Mucmuc1",   "PS15p_vs_Mucmuc1.paf",          "win60"),
    ("Mucmuc1",   "Mucpir1",   "Mucmuc1_vs_Mucpir1.paf",        "win60"),
    # Phycomyces
    ("PhyblU21",  "Phybl4",    "PhyblU21_vs_Phybl4.paf",        "win60"),
    ("Phybl4",    "PnitS607",  "Phybl4_vs_PnitS607.paf",        "win60"),
    ("PnitS607",  "PnitSJTU",  "PnitS607_vs_PnitSJTU_60kb.paf", "win60"),
    # Fenlin pair (BLASTN of 60kb windows; PAF in SYNDIR root, 2026-04-28)
    ("Fenlin",    "Fenlin1",   "Fenlin_vs_Fenlin1.paf",         "win60"),
]

TE_COLORS = None  # populated below from cladogram palette


def _load_palette():
    import sys as _sys
    _sys.path.insert(0, "/data/carlos/EDF-TE/EDF-TE/Section2_TElandscape/scripts")
    from clade_palette import TE_CLASS_COLORS as P
    return {
        "LTR":            P["LTR"],
        "LINE":           P["LINE"],
        "DNA":            P["DNA"],
        "SINE":           P["SINE"],
        "RC":             P["Rolling Circle"],
        "PLE":            P["Penelope"],
        "Penelope":       P["Penelope"],
        "Unknown":        P["Unclassified"],
        "Unclassified":   P["Unclassified"],
        "Simple_repeat":  P["Unclassified"],
        "Low_complexity": P["Unclassified"],
    }

TE_COLORS = _load_palette()

def te_main_class(name):
    if "#" in name:
        return name.split("#")[-1].strip().split("/")[0]
    return "Unknown"

def te_color(cls):
    return TE_COLORS.get(cls, "#bdbdbd")

def load_TE(sp, min_bp=50):
    if sp not in WINDOWS:
        return pd.DataFrame(columns=["start","end","cls","strand","full"])
    # Special path for species annotated with Earl Grey (Fenlin pair, 2026-04-28):
    # parse the filteredRepeats.gff directly instead of the BLAST-vs-RM rmask TSV.
    EG_GFF = {
        "Fenlin":  Path("/data/carlos/EDF-TE/TE_analysis/annotated/MU_Fenlin_EarlGrey/MU_Fenlin_summaryFiles/MU_Fenlin.filteredRepeats.gff"),
        "Fenlin1": Path("/data/carlos/EDF-TE/TE_analysis/annotated/MU_Fenlin1_EarlGrey/MU_Fenlin1_summaryFiles/MU_Fenlin1.filteredRepeats.gff"),
    }
    if sp in EG_GFF:
        gff = EG_GFF[sp]
        if not gff.exists():
            return pd.DataFrame(columns=["start","end","cls","strand","full"])
        vscaf, vs, ve = VIEW[sp]
        rows = []
        with open(gff) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[0] != vscaf:
                    continue
                s, e = int(f[3]), int(f[4])
                if e < vs or s > ve:
                    continue
                if (e - s + 1) < min_bp:
                    continue
                cls    = f[2].split("/")[0].strip() or "Unknown"
                strand = f[6] if f[6] in ("+", "-") else "+"
                # "full" marker: shortte=F flag from Earl Grey AND >= 1 kb hit
                attrs = f[8]
                shortte_f = ("shortte=F" in attrs)
                full = bool(shortte_f and (e - s + 1) >= 1000)
                rows.append((s, e, cls, strand, full))
        return pd.DataFrame(rows, columns=["start","end","cls","strand","full"])
    scaf, win_s, win_e = WINDOWS[sp]
    vscaf, vs, ve = VIEW[sp]
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
    return df[["start","end","cls","strand","full"]].copy()

def load_synteny_win60(q_sp, t_sp, paf_file, min_aln=200):
    """Load BLASTN results from extracted 60kb windows.
    Hit coords are 1-based positions within extracted window.
    Convert to scaffold abs coords: pos_abs = pos_win + view_start - 1.
    """
    # Try WINDIR first, then SYNDIR
    for parent in (WINDIR, SYNDIR):
        fpath = parent / paf_file
        if fpath.exists():
            break
    else:
        return []

    q_vs = VIEW[q_sp][1]; q_ve = VIEW[q_sp][2]; q_vlen = q_ve - q_vs
    t_vs = VIEW[t_sp][1]; t_ve = VIEW[t_sp][2]; t_vlen = t_ve - t_vs

    cols = ["qid","sid","pident","length","qstart","qend","sstart","send","evalue","bitscore"]
    try:
        df = pd.read_csv(fpath, sep="\t", header=None, names=cols)
    except Exception:
        return []
    if len(df) == 0:
        return []

    blocks = []
    for _, r in df.iterrows():
        ln = int(r["length"])
        if ln < min_aln:
            continue
        qs = int(r["qstart"]) + q_vs - 1
        qe = int(r["qend"])   + q_vs - 1
        ss = int(r["sstart"]); se = int(r["send"])
        strand = "-" if ss > se else "+"
        ts = min(ss,se) + t_vs - 1
        te = max(ss,se) + t_vs - 1
        blocks.append(dict(q_start=qs, q_end=qe, t_start=ts, t_end=te,
                           strand=strand, alen=ln))
    return blocks

# ── GFF parsing for background-gene arrows ───────────────────────────────────
GFF_MAP = {
    "Rhisto1":   (GFF_DIR/"MU_Rhisto1.gff",  "scaffold_35",     "jgi"),
    "Rhimi59":   (GFF_DIR/"MU_Rhimi59.gff",  "scaffold_1",      "jgi"),
    "Rhior3":    (GFF_DIR/"MU_Rhior3.gff",   "NW_027137827.1",  "gff3"),
    "Mucci3":    (GFF_DIR/"MU_Mucci3.gff",   "scaffold_1",      "jgi"),
    "Muccirc1":  (GFF_DIR/"MU_Muccirc1.gff", "KE_123975",       "jgi"),
    "Mucmuc1":   (GFF_DIR/"MU_Mucmuc1.gff",  "scaffold_63",     "jgi"),
    "Mucpir1":   (GFF_DIR/"MU_Mucpir1.gff",  "scaffold_29",     "jgi"),
    "PhyblU21":  (Path("/legserv/Temp/Carlos/EDF-TE/data/Extra/PhyblU21_2_GeneCatalog_20191204.gff3"),
                  "scaffold_17", "gff3"),
    "Phybl4":    (GFF_DIR/"MU_Phybl4.gff3",  "scaffold_10",     "gff3"),
    "PnitS607":  (GFF_DIR/"MU_PnitS607.gff", "scaffold_54",     "jgi"),
    # Fenlin pair (2026-04-28)
    "Fenlin":    (GFF_DIR/"MU_Fenlin.gff",  "scaffold_5",      "jgi"),
    "Fenlin1":   (GFF_DIR/"MU_Fenlin1.gff3","scaffold_45",     "gff3"),
    # ATCC52813/52814, PS15m/p, PnitSJTU have no GFF available
}

def _parse_jgi_gff(filepath, scaffold, vs, ve):
    genes = {}
    with open(filepath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f = line.rstrip("\n").split("\t")
            if len(f)<9 or f[0]!=scaffold or f[2]!="CDS": continue
            s,e = int(f[3]),int(f[4])
            if e<vs or s>ve: continue
            strand=f[6]; pid=None
            for tok in f[8].split(";"):
                tok=tok.strip()
                if tok.startswith("proteinId "):
                    pid=tok.split()[1].strip(); break
            if pid is None: continue
            if pid not in genes: genes[pid]=[s,e,strand]
            else: genes[pid][0]=min(genes[pid][0],s); genes[pid][1]=max(genes[pid][1],e)
    return sorted((gs,ge,gst) for gs,ge,gst in genes.values())

def _parse_gff3(filepath, scaffold, vs, ve):
    result=[]
    with open(filepath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<9 or f[0]!=scaffold or f[2]!="gene": continue
            s,e=int(f[3]),int(f[4])
            if e<vs or s>ve: continue
            result.append((s,e,f[6]))
    return sorted(result)

def _load_bg_genes(sp):
    if sp not in GFF_MAP: return []
    filepath, scaffold, fmt = GFF_MAP[sp]
    vs, ve = VIEW[sp][1], VIEW[sp][2]
    try:
        if fmt=="jgi": return _parse_jgi_gff(filepath, scaffold, vs, ve)
        else:          return _parse_gff3(filepath, scaffold, vs, ve)
    except Exception:
        return []

# Species without GFF — use pre-computed TBLASTN clusters (absolute scaffold coords)
_TBLASTN_SPECIES = {"ATCC52813","ATCC52814","PS15m","PS15p","PnitSJTU"}

def _load_bg_tblastn(sp):
    tsv = BGDIR / f"{sp}_bg_genes.tsv"
    if not tsv.exists():
        return []
    vs, ve = VIEW[sp][1], VIEW[sp][2]
    result = []
    with open(tsv) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            s, e, strand = int(parts[0]), int(parts[1]), parts[2]
            if e < vs or s > ve:
                continue
            result.append((s, e, strand))
    return result

ALL_GENES = {}
for sp in TRACK_ORDER:
    if sp in _TBLASTN_SPECIES:
        ALL_GENES[sp] = _load_bg_tblastn(sp)
    else:
        ALL_GENES[sp] = _load_bg_genes(sp)

# ── Drawing ───────────────────────────────────────────────────────────────────
TRACK_H = 0.065

def draw_track(ax, sp):
    _, vs, ve = VIEW[sp]
    vlen = ve - vs

    def x(pos):
        return (pos - vs) / vlen

    yc = 0.5
    ax.axhline(yc, color="#444444", lw=0.9, zorder=1)

    # Idiomorphic region highlight
    idio = IDIO.get(sp)
    if idio:
        ix0 = max(0.0, x(idio[0])); ix1 = min(1.0, x(idio[1]))
        if ix0 < ix1:
            ax.add_patch(mpatches.Rectangle(
                (ix0, yc-0.28), ix1-ix0, 0.56,
                color="#FDEBD0", alpha=0.50, zorder=0, linewidth=0))
            for xi in (ix0, ix1):
                ax.plot([xi,xi],[yc-0.30,yc+0.30], color="#E59866", lw=0.8, ls="--", zorder=2)

    # TE blocks — drawn on a SEPARATE sub-track BELOW the backbone so they
    # never visually overlap the gene arrows (user request 2026-04-22).
    # Anything painted at the backbone (yc) used to collide with gene arrows
    # whenever a RepeatMasker hit was inside / adjacent to a CDS.
    te_yc = yc - TRACK_H * 2.6   # centre of TE sub-track
    # Light divider line so the TE row reads as a track rather than noise
    ax.plot([0, 1], [te_yc, te_yc], color="#bbbbbb", lw=0.4, zorder=1)
    te_df = load_TE(sp)
    for _, tr in te_df.iterrows():
        tx0=x(tr["start"]); tx1=x(tr["end"])
        if tx1<=0 or tx0>=1: continue
        tx0=max(0.0,tx0); tx1=min(1.0,tx1)
        tw=tx1-tx0
        if tw<2e-4: continue
        tcol=te_color(tr["cls"])
        if tr["full"]:
            ast=tx0 if tr["strand"]=="+" else tx1
            dx=tw if tr["strand"]=="+" else -tw
            hl=min(abs(dx)*0.30,0.009)
            ax.add_patch(FancyArrow(ast,te_yc,dx,0,
                width=TRACK_H*1.55,head_width=TRACK_H*2.6,head_length=hl,
                length_includes_head=True,facecolor=tcol,edgecolor="white",
                zorder=4,linewidth=0.3,overhang=0.05))
        else:
            ax.add_patch(mpatches.Rectangle(
                (tx0,te_yc-TRACK_H*0.75),tw,TRACK_H*1.5,
                color=tcol,alpha=0.75,zorder=3,linewidth=0))

    # Background grey gene arrows (all annotated genes from GFF)
    for gs,ge,gstrand in ALL_GENES.get(sp,[]):
        gx0=x(gs); gx1=x(ge)
        if gx1<=0 or gx0>=1: continue
        gx0c=max(0.0,gx0); gx1c=min(1.0,gx1)
        gw=gx1c-gx0c
        if gw<1e-5: continue
        ast=gx0c if gstrand=="+" else gx1c
        dx=gw if gstrand=="+" else -gw
        hl=min(abs(dx)*0.30,0.009)
        ax.add_patch(FancyArrow(ast,yc,dx,0,
            width=TRACK_H*2.0,head_width=TRACK_H*3.0,head_length=hl,
            length_includes_head=True,facecolor="#707B7C",edgecolor="#5D6D7E",
            zorder=5,linewidth=0.2,overhang=0.05))

    # Known sex-locus gene arrows (coloured, labelled). Larger arrowhead +
    # explicit direction glyph + a thick BLACK 5'-end tail bar so the
    # gene direction reads off the position of the black cap, not just the
    # arrowhead — strand inversions become unambiguous (e.g. tptA in
    # PnitS607 vs PnitSJTU). Also stash gene-position metadata on the axis
    # so the group figure can draw cross-track INVERSION connectors.
    ax._sex_gene_positions = []
    for gs,ge,gstrand,gcolor,glabel in GENES.get(sp,[]):
        gx0=x(gs); gx1=x(ge)
        if gx1<=0 or gx0>=1: continue
        gx0c=max(0.0,gx0); gx1c=min(1.0,gx1)
        gw=gx1c-gx0c
        if gw<1e-5: continue
        ast=gx0c if gstrand=="+" else gx1c
        dx=gw if gstrand=="+" else -gw
        hl=min(abs(dx)*0.55, 0.022)            # bigger arrowhead
        ax.add_patch(FancyArrow(ast,yc,dx,0,
            width=TRACK_H*2.0,head_width=TRACK_H*4.4,head_length=hl,
            length_includes_head=True,color=gcolor,zorder=6,linewidth=0,overhang=0.15))
        # 5'-end tail cap (black). + strand → cap on LEFT; − strand → cap on RIGHT.
        cap_w = max(gw * 0.05, 0.0015)
        cap_x0 = gx0c if gstrand == "+" else (gx1c - cap_w)
        ax.add_patch(mpatches.Rectangle(
            (cap_x0, yc - TRACK_H*1.05), cap_w, TRACK_H*2.10,
            facecolor="black", edgecolor="none", zorder=7))
        glyph = "→" if gstrand == "+" else "←"
        ax.text((gx0c+gx1c)/2, yc+TRACK_H*3.9, f"{glabel} {glyph}",
                ha="center",va="bottom",fontsize=6.2,color=gcolor,
                zorder=7,clip_on=True,fontstyle="italic",fontweight="bold")
        ax._sex_gene_positions.append(
            (glabel, (gx0c+gx1c)/2, gstrand, gcolor)
        )

    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")

    mt = MT_LABELS.get(sp,"?")
    mt_col = "#C0392B" if "sexM" in mt else ("#E8671A" if "sexP" in mt else "#777777")
    sp_name = SPECIES_NAMES.get(sp, sp)
    parts = sp_name.split(" ",1)
    genus   = parts[0]
    epithet = parts[1] if len(parts)>1 else ""
    # Shorten for display
    label_str = f"$\\it{{{genus}}}$" + (f" $\\it{{{epithet}}}$" if epithet else "")
    ax.text(-0.01, yc+0.12, label_str,
            ha="right",va="center",fontsize=5.8,transform=ax.transAxes,color="#222222")
    ax.text(-0.01, yc-0.22, f"{sp}  {mt}",
            ha="right",va="top",fontsize=4.8,transform=ax.transAxes,color=mt_col)

def draw_synteny_band(fig, ax_q, ax_t, q_sp, t_sp, paf_file, fmt, min_alen=200):
    blocks = load_synteny_win60(q_sp, t_sp, paf_file, min_alen)
    if not blocks:
        return
    pos_q = ax_q.get_position()
    pos_t = ax_t.get_position()
    gap_y0 = pos_t.y1; gap_y1 = pos_q.y0
    if gap_y1 <= gap_y0 + 1e-4: return
    conn_ax = fig.add_axes([pos_q.x0, gap_y0,
                             pos_q.x1-pos_q.x0, gap_y1-gap_y0])
    conn_ax.set_xlim(0,1); conn_ax.set_ylim(0,1)
    conn_ax.axis("off"); conn_ax.patch.set_visible(False)

    q_vs,q_ve = VIEW[q_sp][1],VIEW[q_sp][2]; q_vlen=q_ve-q_vs
    t_vs,t_ve = VIEW[t_sp][1],VIEW[t_sp][2]; t_vlen=t_ve-t_vs

    for b in blocks:
        if b["alen"] < min_alen: continue
        qx0=(b["q_start"]-q_vs)/q_vlen; qx1=(b["q_end"]-q_vs)/q_vlen
        tx0=(b["t_start"]-t_vs)/t_vlen; tx1=(b["t_end"]-t_vs)/t_vlen
        qx0,qx1=sorted([max(0.0,min(1.0,qx0)),max(0.0,min(1.0,qx1))])
        tx0,tx1=sorted([max(0.0,min(1.0,tx0)),max(0.0,min(1.0,tx1))])
        if qx1-qx0<1e-4 or tx1-tx0<1e-4: continue
        # Uniform palette (2026-04-28, v4): all synteny blocks share the
        # original soft yellow #D5CBA5 regardless of strand.  Inversions are
        # still drawn as a CROSSED trapezoid so the orientation reads from
        # the geometry, not from a colour change.
        col_f = "#D5CBA5"; col_e = "#B8A870"
        if b["strand"] == "+":
            poly = [(qx0,1.0),(qx1,1.0),(tx1,0.0),(tx0,0.0)]
        else:
            poly = [(qx0,1.0),(qx1,1.0),(tx0,0.0),(tx1,0.0)]
        conn_ax.add_patch(mpatches.Polygon(
            poly, closed=True, facecolor=col_f, edgecolor=col_e,
            alpha=0.55, linewidth=0.5, zorder=5))

# ── Build one figure per group ───────────────────────────────────────────────

def draw_group_figure(group_name, members):
    sp_pairs = [(q,t,pf,fm) for (q,t,pf,fm) in SYNTENY_PAIRS if q in members and t in members]
    n = len(members)
    # heights: n track rows + (n-1) synteny gap rows
    heights = []
    for i in range(n):
        heights.append(3)
        if i < n - 1:
            heights.append(1.2)

    fig_h = 0.72 * len(heights) + 1.4
    fig = plt.figure(figsize=(20, fig_h))
    fig.patch.set_facecolor("white")

    gs_main = gridspec.GridSpec(len(heights), 1, figure=fig,
                                height_ratios=heights, hspace=0.55,
                                left=0.26, right=0.91, top=0.97, bottom=0.07)
    axes = {}
    row = 0
    for i, sp in enumerate(members):
        axes[sp] = fig.add_subplot(gs_main[row])
        row += 1
        if i < n - 1:
            row += 1   # skip synteny gap row (filled by draw_synteny_band)

    for sp in members:
        draw_track(axes[sp], sp)

    for q_sp, t_sp, paf_file, fmt in sp_pairs:
        draw_synteny_band(fig, axes[q_sp], axes[t_sp], q_sp, t_sp, paf_file, fmt)

    # ── Cross-track INVERSION connectors ────────────────────────────────────
    # For each coloured sex-locus gene that appears with BOTH + and − strand
    # somewhere in this group, draw a thick figure-level connector between
    # the (+) occurrence and the (−) occurrence with a red "INVERSION" label.
    # Makes the strand flip impossible to miss (e.g. tptA in P. nitens).
    _gene_tbl = defaultdict(list)   # label -> [(sp, x_axes, strand, color)]
    for sp in members:
        for (lbl, xc, strand, gcol) in getattr(axes[sp], "_sex_gene_positions", []):
            _gene_tbl[lbl].append((sp, xc, strand, gcol))
    # Restrict cross-track INVERSION callouts to genes whose strand flip is the
    # biologically interesting feature (tptA): rnhA is colinear with multiple
    # orientations across the panel and would otherwise add a parallel
    # connector that overlaps the tptA one.
    _INVERSION_GENES = {"tptA"}
    for lbl, occs in _gene_tbl.items():
        if lbl not in _INVERSION_GENES:
            continue
        plus  = [o for o in occs if o[2] == "+"]
        minus = [o for o in occs if o[2] == "-"]
        if not plus or not minus:
            continue
        # Pair each (+) with each (−) but emit at most ONE connector per
        # adjacent species pair (avoid overlapping duplicate labels).
        _drawn_pairs = set()
        for (sp_p, xp, _, cp) in plus:
            for (sp_m, xm, _, cm) in minus:
                idx_p = members.index(sp_p); idx_m = members.index(sp_m)
                if abs(idx_p - idx_m) != 1:
                    continue
                pair_key = tuple(sorted((sp_p, sp_m)))
                if pair_key in _drawn_pairs:
                    continue
                _drawn_pairs.add(pair_key)
                ax_p = axes[sp_p]; ax_m = axes[sp_m]
                p_fig = fig.transFigure.inverted().transform(
                    ax_p.transAxes.transform((xp, 0.5)))
                m_fig = fig.transFigure.inverted().transform(
                    ax_m.transAxes.transform((xm, 0.5)))
                # Thick, prominent INVERSION connector (2026-04-28 v3).
                # The BLASTN synteny blocks frequently miss strand-flipped
                # regions because of sequence divergence (e.g. tptA between
                # P. nitens strains), so the gene-position connector is the
                # primary visual cue and must read at-a-glance.
                ln1 = plt.Line2D([p_fig[0], m_fig[0]], [p_fig[1], m_fig[1]],
                                 transform=fig.transFigure,
                                 color="#A04000", lw=2.4, alpha=0.85,
                                 linestyle=(0, (5, 2)), zorder=20,
                                 solid_capstyle="round")
                fig.add_artist(ln1)
                # Add a small filled triangle at each endpoint pointing in the
                # direction of the gene's actual strand on that track.
                for (x_pt, y_pt, strand_pt) in (
                        (p_fig[0], p_fig[1], "+"),
                        (m_fig[0], m_fig[1], "-")):
                    dx = 0.010 if strand_pt == "+" else -0.010
                    tri = mpatches.Polygon(
                        [(x_pt, y_pt - 0.006),
                         (x_pt, y_pt + 0.006),
                         (x_pt + dx, y_pt)],
                        closed=True, transform=fig.transFigure,
                        facecolor="#A04000", edgecolor="#A04000",
                        alpha=0.95, zorder=21)
                    fig.add_artist(tri)
                # 2026-04-28: removed the "INVERSION (tptA)" floating text
                # label per user request.  The dashed red connector and the
                # inverted-block colour gradient already convey the flip.

    # Scale bar on last track
    last_sp = members[-1]
    ax_last = axes[last_sp]
    bar_kb = 10000
    _, vs, ve = VIEW[last_sp]
    bar_frac = bar_kb / (ve - vs)
    bx1 = 0.92; bx0 = bx1 - bar_frac
    ax_last.annotate("", xy=(bx1, 0.05), xytext=(bx0, 0.05),
                     xycoords="axes fraction",
                     arrowprops=dict(arrowstyle="-", color="black", lw=1.5))
    ax_last.text((bx0+bx1)/2, -0.15, "10 kb",
                 ha="center", va="top", fontsize=7, transform=ax_last.transAxes)

    # Group bracket label (right side)
    ax0 = axes[members[0]]; ax1 = axes[members[-1]]
    p0 = ax0.get_position(); p1 = ax1.get_position()
    xb = p0.x1 + 0.012; yt = p0.y1; yb = p1.y0; ym = (yt + yb) / 2
    ax0.annotate("", xy=(xb, yb), xytext=(xb, yt),
                 xycoords="figure fraction", textcoords="figure fraction",
                 arrowprops=dict(arrowstyle="-", color="#444444", lw=1.2))
    for yy in (yt, yb):
        ax0.annotate("", xy=(xb+0.004, yy), xytext=(xb, yy),
                     xycoords="figure fraction", textcoords="figure fraction",
                     arrowprops=dict(arrowstyle="-", color="#444444", lw=1.0))
    fig.text(xb+0.013, ym, group_name,
             ha="left", va="center", fontsize=9, fontweight="bold", color="#333333")

    # Legend
    legend_elems = [
        mpatches.Patch(color="#707B7C",         label="Other annotated genes"),
        mpatches.Patch(color="#E8671A",         label="sexP"),
        mpatches.Patch(color="#C0392B",         label="sexM"),
        mpatches.Patch(color="#27AE60",         label="tptA"),
        mpatches.Patch(color="#2980B9",         label="rnhA"),
        mpatches.Patch(color="#F1C40F",         label="sagA"),
        mpatches.Patch(color="#784212",         label="glrA"),
        mpatches.Patch(color="#8E44AD",         label="helicase"),
        mpatches.Patch(color="#95A5A6",         label="unk (Rhimi59)"),
        mpatches.Patch(color="#FDEBD0", alpha=0.9, label="Idiomorphic region"),
        mpatches.FancyArrow(0, 0, 0.1, 0, width=0.3, head_width=0.5, head_length=0.02,
                            color=TE_COLORS["LTR"], label="Full LTR TE"),
        mpatches.Patch(color=TE_COLORS["LTR"], alpha=0.75, label="LTR remnant"),
        mpatches.FancyArrow(0, 0, 0.1, 0, width=0.3, head_width=0.5, head_length=0.02,
                            color=TE_COLORS["LINE"], label="Full LINE TE"),
        mpatches.Patch(color=TE_COLORS["LINE"], alpha=0.75, label="LINE remnant"),
        mpatches.Patch(color=TE_COLORS["DNA"], alpha=0.75, label="DNA TE"),
        mpatches.Patch(color=TE_COLORS["SINE"], alpha=0.75, label="SINE"),
        mpatches.Patch(color=TE_COLORS["RC"], alpha=0.75, label="RC/Helitron"),
        mpatches.Patch(color=TE_COLORS["PLE"], alpha=0.75, label="Penelope"),
        mpatches.Patch(color=TE_COLORS["Unknown"], alpha=0.75, label="Unclassified TE"),
        mpatches.Patch(color="#D5CBA5", alpha=0.75, label="Synteny block"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", ncol=6,
               fontsize=6.2, frameon=True, framealpha=0.95, edgecolor="lightgray",
               bbox_to_anchor=(0.58, 0.0), bbox_transform=fig.transFigure)

    fig.suptitle(f"Sex determination locus — {group_name} — 60 kb windows",
                 fontsize=11, fontweight="bold", y=0.99)

    out_pdf = WORKDIR / f"sex_locus_synteny_v3_{group_name}.pdf"
    out_png = WORKDIR / f"sex_locus_synteny_v3_{group_name}.png"
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")

    # Also publish to the visible MAT_locus/figures tree (= /Volumes/Temp/...)
    _FIG = Path("/legserv/Temp/Carlos/EDF-TE/results/MAT_locus/figures")
    _FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG / f"sex_locus_synteny_v3_{group_name}.pdf",
                dpi=200, bbox_inches="tight")
    fig.savefig(_FIG / f"sex_locus_synteny_v3_{group_name}.png",
                dpi=150, bbox_inches="tight")
    print(f"Saved: {_FIG}/sex_locus_synteny_v3_{group_name}.{{pdf,png}}")
    plt.close()


for group_name, members in GROUPS:
    draw_group_figure(group_name, members)

