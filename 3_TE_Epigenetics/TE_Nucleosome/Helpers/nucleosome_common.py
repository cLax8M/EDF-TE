#!/usr/bin/env python3
# provides shared helper functions for nucleosome analysis scripts

# handles genomic data parsing, nucleosome metrics, and TE classification

# includes statistical analysis and plot styling utilities

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# TE classification
# ---------------------------------------------------------------------------

TE_CLASS_COLORS = {
    "LTR":            "#e41a1c",
    "LINE":           "#4daf4a",
    "DNA":            "#377eb8",
    "SINE":           "#984ea3",
    "RC/Helitron":    "#ff7f00",
    "Penelope":       "#a65628",
    "Satellite":      "#f781bf",
    "rRNA":           "#a65628",
    "Unknown":        "#999999",
}

TE_CLASSES = ["LTR", "LINE", "DNA", "SINE", "RC/Helitron", "Penelope", "Unknown"]


def classify_te(col4: str) -> str:
    """Map an EarlGrey name/annotation (col4) to a TE class.

    Handles both 'LTR/Gypsy' and 'TE#LTR/Gypsy' style annotations.
    """
    raw = str(col4).strip()
    if "#" in raw:
        raw = raw.split("#", 1)[1]
    c = raw.split("/")[0].rstrip("?").strip().upper()
    if c.startswith("LTR"):
        return "LTR"
    if c.startswith("LINE"):
        return "LINE"
    if c.startswith("DNA"):
        return "DNA"
    if c.startswith("SINE"):
        return "SINE"
    if c.startswith("RC") or c.startswith("HELITRON"):
        return "RC/Helitron"
    if c.startswith("PLE") or c.startswith("PENELOPE"):
        return "Penelope"
    if c.startswith("SAT"):
        return "Satellite"
    if c.startswith("RRNA"):
        return "rRNA"
    return "Unknown"


def load_te_bed(bed_path: str | Path) -> pd.DataFrame:
    """Load an EarlGrey filteredRepeats.bed into a DataFrame with TE class."""
    path = Path(bed_path)
    df = pd.read_csv(path, sep="\t", header=None,
                     usecols=[0, 1, 2, 3],
                     names=["chrom", "start", "end", "name"])
    df["te_class"] = df["name"].apply(classify_te)
    df["length"] = (df["end"] - df["start"]).astype(int)
    return df


def te_intervals_by_chrom(df: pd.DataFrame) -> dict[str, np.ndarray]:
    """Return {chrom: Nx2 int32 array of [start, end)} for fast overlap."""
    out = {}
    for chrom, grp in df.groupby("chrom"):
        out[chrom] = grp[["start", "end"]].to_numpy(np.int32)
    return out


# ---------------------------------------------------------------------------
# WIG / bedGraph loading
# ---------------------------------------------------------------------------

def load_wig(path: str | Path, chrom_sizes: dict[str, int] | None = None):
    """Parse fixedStep/variableStep WIG (optionally gzipped).

    Returns a dict {chrom: np.float32 array}.  If ``chrom_sizes`` is given the
    arrays are pre-allocated to the chromosome lengths; otherwise they are
    grown to the largest position seen.  Chromosome names that carry a
    'SPECIES.' prefix are stripped to match the FAI / BED naming.
    """
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    arrays: dict[str, np.ndarray] = {}
    cur = None
    pos = 0
    step = 1
    span = 1
    mode = None

    def _get(chrom):
        if chrom not in arrays:
            if chrom_sizes is not None and chrom in chrom_sizes:
                arrays[chrom] = np.zeros(chrom_sizes[chrom], dtype=np.float32)
            else:
                arrays[chrom] = np.zeros(1 << 26, dtype=np.float32)  # 64 Mb grow
        return arrays[chrom]

    def _norm(chrom):
        if chrom_sizes is None or chrom in chrom_sizes:
            return chrom
        if "." in chrom:
            cand = chrom.split(".", 1)[1]
            if cand in chrom_sizes:
                return cand
        return chrom

    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith(("track", "browser", "#")):
                continue
            if line.startswith("fixedStep"):
                kv = dict(tok.split("=", 1) for tok in line.split()[1:] if "=" in tok)
                chrom = _norm(kv.get("chrom", ""))
                cur = chrom if (chrom_sizes is None or chrom in chrom_sizes) else None
                pos = int(kv.get("start", 1)) - 1
                step = int(kv.get("step", 1))
                span = int(kv.get("span", 1))
                mode = "fixed"
                continue
            if line.startswith("variableStep"):
                kv = dict(tok.split("=", 1) for tok in line.split()[1:] if "=" in tok)
                chrom = _norm(kv.get("chrom", ""))
                cur = chrom if (chrom_sizes is None or chrom in chrom_sizes) else None
                span = int(kv.get("span", 1))
                mode = "var"
                continue
            if cur is None:
                continue
            arr = _get(cur)
            if mode == "fixed":
                try:
                    v = float(line)
                except ValueError:
                    continue
                if 0 <= pos < len(arr):
                    arr[pos:min(pos + span, len(arr))] = v
                pos += step
            elif mode == "var":
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    p = int(parts[0]) - 1
                    v = float(parts[1])
                except ValueError:
                    continue
                if 0 <= p < len(arr):
                    arr[p:min(p + span, len(arr))] = v
    return arrays


def load_chrom_sizes(fai_path: str | Path) -> dict[str, int]:
    df = pd.read_csv(fai_path, sep="\t", header=None, usecols=[0, 1],
                     names=["chrom", "length"])
    return dict(zip(df["chrom"], df["length"].astype(int)))


# ---------------------------------------------------------------------------
# Window metrics
# ---------------------------------------------------------------------------

def autocorr_peak(sig, lag_min=140, lag_max=240):
    """Autocorrelation peak in the given lag range.

    Returns (ac_lag, ac_val): the lag (bp) of the strongest autocorrelation
    peak in [lag_min, lag_max] and its normalised value.  ac_lag is the
    nucleosome repeat length (NRL) proxy; ac_val is the phasing strength.
    """
    sig = sig - sig.mean()
    if sig.std() < 1e-6 or len(sig) < lag_max + 10:
        return np.nan, 0.0
    n = len(sig)
    f = np.fft.rfft(sig, n=2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac = ac / ac[0]
    seg = ac[lag_min:lag_max + 1]
    if len(seg) == 0:
        return np.nan, 0.0
    i = int(np.argmax(seg))
    return lag_min + i, float(seg[i])


def window_metrics(occ: np.ndarray, win: int):
    """Tiled-window summary of a per-base occupancy array.

    Returns a dict of numpy arrays (one entry per window):
        mean_occ, cv, ac_lag, ac_val, valid
    """
    n = len(occ) // win
    means = np.full(n, np.nan, dtype=np.float32)
    cvs = np.full(n, np.nan, dtype=np.float32)
    aclags = np.full(n, np.nan, dtype=np.float32)
    acvals = np.zeros(n, dtype=np.float32)
    valid = np.zeros(n, dtype=bool)
    for i in range(n):
        s = i * win
        seg = occ[s:s + win]
        m = seg.mean()
        if not np.isfinite(m) or seg.size < win:
            continue
        means[i] = m
        sd = seg.std()
        cvs[i] = sd / (abs(m) + 1e-9)
        lag, val = autocorr_peak(seg)
        aclags[i] = lag
        acvals[i] = val
        valid[i] = True
    return dict(mean_occ=means, cv=cvs, ac_lag=aclags, ac_val=acvals, valid=valid)


# ---------------------------------------------------------------------------
# DANPOS positions.xls
# ---------------------------------------------------------------------------

DANPOS_COLS = ["chr", "start", "end", "smt_pos", "smt_value", "fuzziness_score"]


def load_danpos_positions(path: str | Path) -> pd.DataFrame:
    """Load a DANPOS ``positions.xls`` file (whitespace-separated).

    Keeps the essential per-nucleosome columns:
        chr, start, end, smt_pos (dyad), smt_value (peak height),
        fuzziness_score (positional SD, bp).
    """
    df = pd.read_csv(path, sep=r"\s+", comment="#")
    keep = [c for c in DANPOS_COLS if c in df.columns]
    return df[keep].copy()


def add_nrl(df: pd.DataFrame, lo=100, hi=400) -> pd.DataFrame:
    """Distance to the next dyad on the same chromosome (NRL proxy).

    Only distances in [lo, hi] bp are kept (gap/assembly-break filter).
    """
    df = df.sort_values(["chr", "smt_pos"]).copy()
    df["nrl_bp"] = np.nan
    for chrom, grp in df.groupby("chr"):
        pos = grp["smt_pos"].to_numpy()
        if len(pos) < 2:
            continue
        dists = np.diff(pos).astype(float)
        dists[(dists < lo) | (dists > hi)] = np.nan
        df.loc[grp.index[:-1], "nrl_bp"] = dists
    return df


def classify_nucleosome_context(chrom, pos, te_dict, gene_dict=None):
    """Assign a dyad to TE / gene / intergenic.

    TE takes precedence over gene.  Returns one of "TE", "Gene", "Inter".
    """
    def _inside(ivs):
        if ivs is None or len(ivs) == 0:
            return False
        idx = np.searchsorted(ivs[:, 0], pos, side="right") - 1
        return idx >= 0 and pos < ivs[idx, 1]

    if _inside(te_dict.get(chrom)):
        return "TE"
    if gene_dict is not None and _inside(gene_dict.get(chrom)):
        return "Gene"
    return "Inter"


def load_gene_intervals(gff_path: str | Path) -> dict[str, np.ndarray]:
    """Load gene-body intervals from a GFF/GFF3.

    Uses 'gene' features if present, otherwise falls back to 'exon'.
    Returns {chrom: Nx2 int32 array} with 0-based half-open coords.
    """
    want = set()
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) >= 3:
                want.add(f[2])
    use = "gene" if "gene" in want else "exon"

    d = {}
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 5 or f[2] != use:
                continue
            try:
                s, e = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            d.setdefault(f[0], []).append((s, e))
    return {c: np.array(sorted(v), dtype=np.int32) for c, v in d.items()}


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR-adjusted p-values (5 % default threshold)."""
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = ranked
    return np.clip(out, 0, 1)


def cliffs_delta(x, y, nmax=5000, seed=42):
    """Cliff's delta (rank-biserial effect size) between two groups.

    Returns a value in [-1, 1]; |d| > 0.147, 0.33, 0.474 correspond to small,
    medium and large effects (Vargha & Delaney).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size == 0 or y.size == 0:
        return np.nan
    rng = np.random.RandomState(seed)
    if x.size > nmax:
        x = rng.choice(x, nmax, replace=False)
    if y.size > nmax:
        y = rng.choice(y, nmax, replace=False)
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return float(gt - lt) / (x.size * y.size)


def fmt_p(p: float) -> str:
    return "nan" if not np.isfinite(p) else (f"{p:.2g}" if p < 1e-3 else f"{p:.3f}")


# ---------------------------------------------------------------------------
# Plotting style
# ---------------------------------------------------------------------------

def style_ax(ax, **kw):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=kw.pop("labelsize", 9))
    return ax
