#!/usr/bin/env python3
"""
08_plot_circos.py — Matplotlib-based circos-style synteny plot for species pairs.

Reads BED, genome, and collinear files from the synteny pipeline and generates
a circular chord diagram showing syntenic blocks between two species.

Usage:
  python 08_plot_circos.py --sp_a SP_A --sp_b SP_B --clade CLADE \
      --bed_dir BED_DIR --collinear_dir COLLINEAR_DIR --out_dir OUT_DIR \
      [--max_scaffolds 25] [--min_links 3]

Or with --batch to generate multiple pairs at once:
  python 08_plot_circos.py --batch --bed_dir BED_DIR --collinear_dir COLLINEAR_DIR \
      --out_dir OUT_DIR --ratio_file RATIO_FILE

Author: EDF-TE project
"""

import argparse, sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MPath
import matplotlib.patches as mpatches_patch


# ── Colors ────────────────────────────────────────────────────────────────────
SP_A_COLORS = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
    "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62",
    "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494",
    "#b3b3b3", "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
    "#66a61e", "#e6ab02", "#a6761d", "#666666", "#8dd3c7",
]
SP_B_COLOR = "#c7c7c7"  # neutral grey for species B scaffolds

# ── TE annotation ─────────────────────────────────────────────────────────────
TE_ANNOT_DIR = Path("/data/carlos/EDF-TE/TE_analysis/annotated")
_TE_SKIP = frozenset([
    "Low_complexity", "Simple_repeat", "rRNA", "Satellite",
    "Telomeric_repeat", "scRNA", "snRNA", "srpRNA", "tRNA",
])

# Class I (RNA-based transposons): LTR, LINE, SINE, Retroposons
_CLASS1_PREFIXES = ("LTR", "LINE", "SINE", "Retroposon")
# Class II (DNA-based transposons): DNA TEs, Helitrons, RC elements
_CLASS2_PREFIXES = ("DNA", "RC")

# Colors for Class I / II split track
_COLOR_CLASS1 = "#f0c010"   # amber/yellow  – Class I (LTR, LINE, SINE)
_COLOR_CLASS2 = "#8040b0"   # purple        – Class II (DNA TEs, Helitron)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sp_a", help="Species A prefix")
    p.add_argument("--sp_b", help="Species B prefix")
    p.add_argument("--clade", help="Clade code")
    p.add_argument("--bed_dir", required=True)
    p.add_argument("--collinear_dir", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--max_scaffolds", type=int, default=25)
    p.add_argument("--min_links", type=int, default=3,
                   help="Min collinear pairs on a scaffold to include")
    p.add_argument("--top_by_size", type=int, default=0,
                   help="If >0, pick top N scaffolds per species by SIZE "
                        "instead of by link count (Merenyi-style clean view)")
    p.add_argument("--simple_style", action="store_true",
                   help="Single colour per species + minimalist title 'X/Y (Z%)'")
    p.add_argument("--batch", action="store_true",
                   help="Run multiple pairs from --pairs list")
    p.add_argument("--pairs", nargs="*",
                   help="CLADE:SpA:SpB triplets for batch mode")
    p.add_argument("--class-split", action="store_true", dest="class_split",
                   help="TE track: Class I (yellow) + Class II (purple) stacked bars")
    return p.parse_args()


def load_genome(bed_dir, sp):
    """Load genome sizes."""
    gf = Path(bed_dir) / f"{sp}.genome"
    sizes = {}
    with open(gf) as fh:
        for line in fh:
            c = line.strip().split("\t")
            if len(c) >= 2:
                sizes[c[0]] = int(c[1])
    return sizes


def load_bed(bed_dir, sp):
    """Load BED positions."""
    bf = Path(bed_dir) / f"{sp}.bed"
    pos = {}
    with open(bf) as fh:
        for line in fh:
            c = line.strip().split("\t")
            if len(c) >= 4:
                pos[c[3]] = (c[0], int(c[1]), int(c[2]))
    return pos


def load_collinear(collinear_dir, clade, sp_a, sp_b):
    """Load collinear pairs."""
    d = Path(collinear_dir) / clade
    f1 = d / f"{sp_a}__{sp_b}.collinear.tsv"
    f2 = d / f"{sp_b}__{sp_a}.collinear.tsv"
    target = f1 if f1.exists() else (f2 if f2.exists() else None)
    if not target:
        print(f"ERROR: no collinear file for {sp_a} × {sp_b}", file=sys.stderr)
        return []
    rev = (target == f2)
    pairs = []
    with open(target) as fh:
        next(fh)
        for line in fh:
            c = line.strip().split("\t")
            if len(c) < 6:
                continue
            if rev:
                pairs.append((c[3], c[4], c[5], c[0], c[1], c[2]))
            else:
                pairs.append(tuple(c[:6]))
    return pairs


def load_te_density(sp, scaffolds, genome, n_windows=50, split_classes=False):
    """Compute per-window TE density fraction for selected scaffolds.

    Returns:
      split_classes=False: {scaffold: [(frac_start, frac_end, density), ...]}
      split_classes=True:  {scaffold: [(frac_start, frac_end, dens_I, dens_II), ...]}
    """
    scaf_set = set(scaffolds)
    base = TE_ANNOT_DIR / f"{sp}_EarlGrey"
    bed = base / f"{sp}_mergedRepeats" / "looseMerge" / f"{sp}.filteredRepeats.bed"
    if not bed.exists():
        bed = base / f"{sp}_summaryFiles" / f"{sp}.filteredRepeats.bed"
    if not bed.exists():
        print(f"  [TE] No BED found for {sp}, skipping TE track")
        return {}

    # Store intervals by scaffold; for class split, store (start, end, class_id)
    # class_id: 1 = Class I, 2 = Class II, 0 = unknown/other
    te_ivs: dict = defaultdict(list)
    with open(bed) as fh:
        for line in fh:
            c = line.strip().split("\t")
            if len(c) < 4:
                continue
            scaf, start, end, cls = c[0], int(c[1]), int(c[2]), c[3]
            if scaf not in scaf_set or cls in _TE_SKIP:
                continue
            if split_classes:
                top = cls.split("/")[0]
                if top in _CLASS1_PREFIXES:
                    cid = 1
                elif top in _CLASS2_PREFIXES:
                    cid = 2
                else:
                    cid = 0   # Unknown class — skip for class split view
                if cid > 0:
                    te_ivs[scaf].append((start, end, cid))
            else:
                te_ivs[scaf].append((start, end))

    density = {}
    for scaf in scaffolds:
        sz = genome.get(scaf, 0)
        if sz == 0:
            continue
        wins = []
        for i in range(n_windows):
            ws = int(i * sz / n_windows)
            we = int((i + 1) * sz / n_windows)
            wlen = max(we - ws, 1)
            if split_classes:
                ivs1 = sorted((s, e) for s, e, cid in te_ivs.get(scaf, []) if cid == 1)
                ivs2 = sorted((s, e) for s, e, cid in te_ivs.get(scaf, []) if cid == 2)
                def _cov(ivs):
                    return sum(
                        min(te, we) - max(ts, ws)
                        for ts, te in ivs
                        if ts < we and te > ws
                    )
                dens_I  = min(_cov(ivs1) / wlen, 1.0)
                dens_II = min(_cov(ivs2) / wlen, 1.0 - dens_I)
                wins.append((ws / sz, we / sz, dens_I, dens_II))
            else:
                ivs = sorted(te_ivs.get(scaf, []))
                covered = sum(
                    min(te, we) - max(ts, ws)
                    for ts, te in ivs
                    if ts < we and te > ws
                )
                wins.append((ws / sz, we / sz, min(covered / wlen, 1.0)))
        density[scaf] = wins
    return density


def select_scaffolds(pairs, pos_a, pos_b, max_scaf, min_links,
                     genome_a=None, genome_b=None, top_by_size=0):
    """Select scaffolds either by collinear pair count (default) or by
    scaffold SIZE (top_by_size>0, Merenyi-style clean view)."""
    cnt_a = defaultdict(int)
    cnt_b = defaultdict(int)

    resolved = []
    for scaf_a, g_a1, g_a2, scaf_b, g_b1, g_b2 in pairs:
        pa1 = pos_a.get(g_a1)
        pb1 = pos_b.get(g_b1)
        if pa1 and pb1:
            resolved.append((scaf_a, pa1, scaf_b, pb1))
            cnt_a[scaf_a] += 1
            cnt_b[scaf_b] += 1

    if top_by_size > 0 and genome_a and genome_b:
        # Pick top N scaffolds by genome size, sorted descending
        top_a = [s for s, _ in sorted(genome_a.items(), key=lambda x: -x[1])[:top_by_size]]
        top_b = [s for s, _ in sorted(genome_b.items(), key=lambda x: -x[1])[:top_by_size]]
    else:
        top_a = [s for s, n in sorted(cnt_a.items(), key=lambda x: -x[1])[:max_scaf]
                 if n >= min_links]
        top_b = [s for s, n in sorted(cnt_b.items(), key=lambda x: -x[1])[:max_scaf]
                 if n >= min_links]

    set_a = set(top_a)
    set_b = set(top_b)
    filt = [(sa, pa, sb, pb) for sa, pa, sb, pb in resolved
            if sa in set_a and sb in set_b]

    return filt, top_a, top_b


def make_circos_plot(sp_a, sp_b, pairs, scafs_a, scafs_b,
                     genome_a, genome_b, out_path, synteny_score=None,
                     patristic_dist=None, te_pct_a=None, te_pct_b=None,
                     ri_a=None, ri_b=None,
                     simple_style=False, n_total_pairs=None,
                     n_orthologs=None,
                     te_density_a=None, te_density_b=None):
    """Generate circos-style matplotlib figure."""

    gap_frac = 0.02      # gap between scaffolds (fraction of circle)
    sp_gap_frac = 0.06   # larger gap between species

    # Total genome for selected scaffolds
    tot_a = sum(genome_a.get(s, 0) for s in scafs_a)
    tot_b = sum(genome_b.get(s, 0) for s in scafs_b)
    total = tot_a + tot_b
    if total == 0:
        print(f"  WARNING: no scaffold data for {sp_a} × {sp_b}, skipping")
        return

    n_gaps = len(scafs_a) + len(scafs_b)
    total_gap = gap_frac * n_gaps + sp_gap_frac * 2
    avail = 2 * np.pi * (1 - total_gap / (2 * np.pi + total_gap))
    # Simpler: allocate 2*pi minus gaps
    avail = 2 * np.pi - (n_gaps * gap_frac * 2 * np.pi / (2 * np.pi)) - 2 * sp_gap_frac * 2 * np.pi / (2 * np.pi)
    # Let's just do it properly
    n_gap_small = len(scafs_a) - 1 + len(scafs_b) - 1
    n_gap_large = 2  # between species
    total_gap_angle = n_gap_small * gap_frac + n_gap_large * sp_gap_frac
    avail = 2 * np.pi - total_gap_angle

    # Equal perimeters: each species occupies half the available arc,
    # regardless of genome size.  Within each half the scaffolds are
    # scaled proportionally to their individual lengths.
    avail_a = avail / 2
    avail_b = avail / 2

    # Build angular layout: species A on top (right half), species B on bottom (left half)
    arc_map = {}  # scaffold -> (start_angle, end_angle)
    color_map = {}

    angle = sp_gap_frac / 2  # start with half large gap

    # Species A scaffolds
    sp_a_solid = "#e6772a"  # orange
    sp_b_solid = "#2d8b54"  # green
    for i, s in enumerate(scafs_a):
        sz = genome_a.get(s, 0)
        span = (sz / tot_a) * avail_a if tot_a > 0 else 0
        arc_map[("A", s)] = (angle, angle + span)
        if simple_style:
            color_map[("A", s)] = sp_a_solid
        else:
            color_map[("A", s)] = SP_A_COLORS[i % len(SP_A_COLORS)]
        angle += span + gap_frac

    angle += sp_gap_frac - gap_frac  # large gap between species

    # Species B scaffolds
    for i, s in enumerate(scafs_b):
        sz = genome_b.get(s, 0)
        span = (sz / tot_b) * avail_b if tot_b > 0 else 0
        arc_map[("B", s)] = (angle, angle + span)
        color_map[("B", s)] = sp_b_solid if simple_style else SP_B_COLOR
        angle += span + gap_frac

    # ── Figure ────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})
    ax.set_ylim(0, 1.35)
    ax.set_axis_off()

    ring_inner = 0.95
    ring_outer = 1.05

    # Draw arcs
    for key, (a0, a1) in arc_map.items():
        theta = np.linspace(a0, a1, max(2, int((a1 - a0) * 200)))
        ax.fill_between(theta, ring_inner, ring_outer,
                        color=color_map[key], alpha=0.85, edgecolor="white",
                        linewidth=0.3)

    # Draw scaffold labels (only for larger scaffolds)
    for key, (a0, a1) in arc_map.items():
        span = a1 - a0
        if span > 0.05:  # only label scaffolds with enough room
            mid = (a0 + a1) / 2
            sp_tag, scaf = key
            label = scaf if len(scaf) < 20 else scaf[:15] + "..."
            rot = np.degrees(mid)
            ha = "left" if 0 <= mid <= np.pi else "right"
            if np.pi/2 < mid < 3*np.pi/2:
                rot += 180
            ax.text(mid, ring_outer + 0.03, label, fontsize=4,
                    ha="center", va="center", rotation=rot - 90,
                    rotation_mode="anchor")

    # Draw links as Bezier curves through center
    for scaf_a, pa, scaf_b, pb in pairs:
        key_a = ("A", scaf_a)
        key_b = ("B", scaf_b)
        if key_a not in arc_map or key_b not in arc_map:
            continue

        a0_a, a1_a = arc_map[key_a]
        a0_b, a1_b = arc_map[key_b]
        sz_a = genome_a.get(scaf_a, 1)
        sz_b = genome_b.get(scaf_b, 1)

        # Position of gene in angular coordinates
        _, start_a, end_a = pa
        _, start_b, end_b = pb
        frac_a = (start_a + end_a) / 2 / sz_a
        frac_b = (start_b + end_b) / 2 / sz_b

        theta_a = a0_a + frac_a * (a1_a - a0_a)
        theta_b = a0_b + frac_b * (a1_b - a0_b)

        # Bezier through center (r=0)
        # Convert to cartesian, draw quadratic Bezier
        r = ring_inner - 0.02
        x0 = r * np.cos(theta_a)
        y0 = r * np.sin(theta_a)
        x1 = r * np.cos(theta_b)
        y1 = r * np.sin(theta_b)

        # Control point at center (slightly offset for aesthetics)
        cx = 0.0
        cy = 0.0

        # Draw as simple arc using polar coordinates
        # Use quadratic Bezier in cartesian overlay
        link_color = color_map[key_a]
        ax.plot([theta_a, theta_b], [ring_inner - 0.02, ring_inner - 0.02],
                color=link_color, alpha=0.08, linewidth=0.3,
                solid_capstyle="round")

    # Actually, polar plot lines are straight in polar coords which looks bad.
    # Let's switch to a non-polar approach with manual circle drawing.
    plt.close(fig)

    # ── Better approach: cartesian figure with manual arcs and Bezier links ──
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.set_aspect("equal")
    ax.axis("off")

    R_inner = 0.90
    R_outer = 1.0
    R_TE_BASE = R_outer + 0.02
    R_TE_MAX  = R_outer + 0.16

    # Draw scaffold arcs
    for key, (a0, a1) in arc_map.items():
        n_pts = max(50, int((a1 - a0) * 300))
        theta = np.linspace(a0, a1, n_pts)
        # Filled wedge
        x_outer = R_outer * np.cos(theta)
        y_outer = R_outer * np.sin(theta)
        x_inner = R_inner * np.cos(theta[::-1])
        y_inner = R_inner * np.sin(theta[::-1])

        xs = np.concatenate([x_outer, x_inner])
        ys = np.concatenate([y_outer, y_inner])
        ax.fill(xs, ys, color=color_map[key], alpha=0.85, edgecolor="white",
                linewidth=0.5)

    # Draw TE density track (outer ring histogram)
    if te_density_a or te_density_b:
        # Detect whether we have class-split (4-tuple) or merged (3-tuple) data
        _first_dens = (te_density_a or te_density_b)
        _sample = next(iter(_first_dens.values()), [])
        _is_split = len(_sample) > 0 and len(_sample[0]) == 4

        for key, (a0, a1) in arc_map.items():
            sp_tag, scaf = key
            dens_map = te_density_a if sp_tag == "A" else te_density_b
            if not dens_map or scaf not in dens_map:
                continue
            col = color_map[key]
            for win in dens_map[scaf]:
                if _is_split:
                    frac_s, frac_e, dens_I, dens_II = win
                    if dens_I <= 0 and dens_II <= 0:
                        continue
                    t0 = a0 + frac_s * (a1 - a0)
                    t1 = a0 + frac_e * (a1 - a0)
                    n_pts = max(3, int((t1 - t0) * 300))
                    theta = np.linspace(t0, t1, n_pts)
                    r_base = R_TE_BASE
                    # Draw Class I (yellow) from base to base+dens_I
                    if dens_I > 0:
                        r_top_I = r_base + dens_I * (R_TE_MAX - R_TE_BASE)
                        xs = np.concatenate([r_top_I * np.cos(theta),
                                             r_base * np.cos(theta[::-1])])
                        ys = np.concatenate([r_top_I * np.sin(theta),
                                             r_base * np.sin(theta[::-1])])
                        ax.fill(xs, ys, color=_COLOR_CLASS1, alpha=0.7, linewidth=0)
                    # Draw Class II (purple) stacked on top of Class I
                    if dens_II > 0:
                        r_bot_II = r_base + dens_I * (R_TE_MAX - R_TE_BASE)
                        r_top_II = r_base + (dens_I + dens_II) * (R_TE_MAX - R_TE_BASE)
                        xs = np.concatenate([r_top_II * np.cos(theta),
                                             r_bot_II * np.cos(theta[::-1])])
                        ys = np.concatenate([r_top_II * np.sin(theta),
                                             r_bot_II * np.sin(theta[::-1])])
                        ax.fill(xs, ys, color=_COLOR_CLASS2, alpha=0.7, linewidth=0)
                else:
                    frac_s, frac_e, dens = win
                    if dens <= 0:
                        continue
                    t0 = a0 + frac_s * (a1 - a0)
                    t1 = a0 + frac_e * (a1 - a0)
                    r_top = R_TE_BASE + dens * (R_TE_MAX - R_TE_BASE)
                    n_pts = max(3, int((t1 - t0) * 300))
                    theta = np.linspace(t0, t1, n_pts)
                    xs = np.concatenate([
                        r_top * np.cos(theta),
                        R_TE_BASE * np.cos(theta[::-1])
                    ])
                    ys = np.concatenate([
                        r_top * np.sin(theta),
                        R_TE_BASE * np.sin(theta[::-1])
                    ])
                    ax.fill(xs, ys, color=col, alpha=0.55, linewidth=0)

    # Draw scaffold labels
    for key, (a0, a1) in arc_map.items():
        span = a1 - a0
        if span > 0.08:
            mid = (a0 + a1) / 2
            sp_tag, scaf = key
            label = scaf if len(scaf) < 18 else scaf[:14] + "…"
            lx = (R_TE_MAX + 0.05) * np.cos(mid)
            ly = (R_TE_MAX + 0.05) * np.sin(mid)
            rot = np.degrees(mid) - 90
            if -180 < rot < -90 or 90 < rot < 180:
                rot += 180
            ax.text(lx, ly, label, fontsize=4, ha="center", va="center",
                    rotation=rot, rotation_mode="anchor")

    # Draw Bezier links
    link_r = R_inner - 0.01

    # Pre-count ortholog multiplicity per scaffold pair (used by simple_style
    # to emphasize true syntenic blocks vs scattered orphans). Pairs that link
    # the same scaffold pair many times are real synteny; lone pairs are noise.
    from collections import Counter
    pair_counts = Counter((scaf_a, scaf_b) for scaf_a, _, scaf_b, _ in pairs)
    max_pc = max(pair_counts.values()) if pair_counts else 1
    # For each scaffold, total links and the count of its top partner — used
    # to decide whether a scaffold-pair is the *dominant* (=conserved synteny)
    # partner or just one of many promiscuous matches (=fragmented synteny).
    a_scaf_total = Counter()
    b_scaf_total = Counter()
    for (sa, sb), c in pair_counts.items():
        a_scaf_total[sa] += c
        b_scaf_total[sb] += c

    for scaf_a, pa, scaf_b, pb in pairs:
        key_a = ("A", scaf_a)
        key_b = ("B", scaf_b)
        if key_a not in arc_map or key_b not in arc_map:
            continue

        a0_a, a1_a = arc_map[key_a]
        a0_b, a1_b = arc_map[key_b]
        sz_a = genome_a.get(scaf_a, 1)
        sz_b = genome_b.get(scaf_b, 1)

        _, start_a, end_a = pa
        _, start_b, end_b = pb
        frac_a = np.clip((start_a + end_a) / 2 / sz_a, 0, 1)
        frac_b = np.clip((start_b + end_b) / 2 / sz_b, 0, 1)

        theta_a = a0_a + frac_a * (a1_a - a0_a)
        theta_b = a0_b + frac_b * (a1_b - a0_b)

        x0 = link_r * np.cos(theta_a)
        y0 = link_r * np.sin(theta_a)
        x1 = link_r * np.cos(theta_b)
        y1 = link_r * np.sin(theta_b)

        # Quadratic Bezier through origin
        verts = [(x0, y0), (0, 0), (x1, y1)]
        codes = [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3]
        path = MPath(verts, codes)
        link_col = "#6e6e6e" if simple_style else color_map[key_a]
        link_alpha = 0.35 if simple_style else 0.12
        link_lw = 0.55 if simple_style else 0.3
        if simple_style:
            # Uniform blue palette for all pairs regardless of synteny score
            _hi_col  = "#1f3a8a"   # deep navy   (dominant synteny block)
            _mid_col = "#7bb3e0"   # light blue  (moderate synteny block)
            _lo_col  = "#cfcfcf"   # pale grey   (noise/rearranged)

            # Scale per-link alpha inversely with pair count so that dense
            # (high-synteny) pairs produce a visually denser ribbon while
            # sparse (low-synteny) pairs keep each individual link visible.
            n_links = len(pairs)
            # dens_scale: 1.0 at ≤100 links, 0.25 at ≥3000 links (log-scale)
            dens_scale = max(0.25, 1.0 - 0.3 * np.log10(max(n_links, 1) / 100.0))

            # Visual emphasis = how *dominant* this scaffold-pair is for both
            # scaffolds. Conserved synteny → one scaffold's genes mostly go to
            # one partner. Fragmented/eroded synteny → many small partners.
            pc = pair_counts.get((scaf_a, scaf_b), 1)
            fa = pc / max(a_scaf_total.get(scaf_a, 1), 1)   # share for sp A scaffold
            fb = pc / max(b_scaf_total.get(scaf_b, 1), 1)   # share for sp B scaffold
            dominance = min(fa, fb)                          # 0–1, conservative
            if dominance >= 0.40 and pc >= 5:
                # strong syntenic block
                link_col = _hi_col
                link_alpha = min(0.85, 0.65 * dens_scale * 1.5)
                link_lw = 0.55 + 0.6 * dominance
            elif dominance >= 0.15 and pc >= 2:
                link_col = _mid_col
                link_alpha = 0.30 * dens_scale
                link_lw = 0.35
            else:
                # rearranged / orphan / TE-mediated noise
                link_col = _lo_col
                link_alpha = 0.12 * dens_scale
                link_lw = 0.18
        patch = mpatches_patch.PathPatch(path, facecolor="none",
                                         edgecolor=link_col,
                                         alpha=link_alpha, linewidth=link_lw)
        ax.add_patch(patch)

    # ── Legend and title ──────────────────────────────────────────────────
    sp_a_short = sp_a.split("_", 1)[1] if "_" in sp_a else sp_a
    sp_b_short = sp_b.split("_", 1)[1] if "_" in sp_b else sp_b

    # Species labels
    # Find midpoint of each species' arc span
    if not simple_style:
        r_sp_label = R_TE_MAX + 0.18
        if scafs_a:
            mid_a = (arc_map[("A", scafs_a[0])][0] + arc_map[("A", scafs_a[-1])][1]) / 2
            lx = r_sp_label * np.cos(mid_a)
            ly = r_sp_label * np.sin(mid_a)
            ax.text(lx, ly, sp_a_short, fontsize=10, ha="center", va="center",
                    fontweight="bold", style="italic",
                    color=SP_A_COLORS[0])
        if scafs_b:
            mid_b = (arc_map[("B", scafs_b[0])][0] + arc_map[("B", scafs_b[-1])][1]) / 2
            lx = r_sp_label * np.cos(mid_b)
            ly = r_sp_label * np.sin(mid_b)
            ax.text(lx, ly, sp_b_short, fontsize=10, ha="center", va="center",
                    fontweight="bold", style="italic", color="#555")

    title = f"Synteny: {sp_a} × {sp_b}\n({len(pairs)} collinear pairs"
    if synteny_score is not None:
        title += f", r = {synteny_score:.3f}"
    if patristic_dist is not None:
        title += f", d = {patristic_dist:.3f}"
    title += ")"
    extra = []
    if te_pct_a is not None and te_pct_b is not None:
        extra.append(f"TE: {te_pct_a:.1f}% / {te_pct_b:.1f}%")
    if ri_a is not None and ri_b is not None:
        extra.append(f"RI: {ri_a:.2f} / {ri_b:.2f}")
    if extra:
        title += "\n" + "   |   ".join(extra)

    if simple_style:
        # Merenyi-style: "Sp1-Sp2" / "N/Total (X.X%)" + r and d
        nb = n_total_pairs if n_total_pairs is not None else len(pairs)
        denom = n_orthologs if n_orthologs is not None else nb
        pct = (synteny_score * 100) if synteny_score is not None else (100.0 * nb / denom if denom else 0.0)
        sub_bits = [f"{nb}/{denom} ({pct:.1f}%)"]
        if synteny_score is not None:
            sub_bits.append(f"r = {synteny_score:.3f}")
        if patristic_dist is not None:
            sub_bits.append(f"d = {patristic_dist:.3f}")
        title = (f"$\\it{{{sp_a_short}}}$ \u2013 $\\it{{{sp_b_short}}}$\n"
                 + "   |   ".join(sub_bits))
        if te_pct_a is not None or te_pct_b is not None:
            te_parts = []
            if te_pct_a is not None:
                te_parts.append(f"{sp_a_short}: {te_pct_a:.1f}%")
            if te_pct_b is not None:
                te_parts.append(f"{sp_b_short}: {te_pct_b:.1f}%")
            title += "\nTE: " + "   |   ".join(te_parts)
        ax.set_title(title, fontsize=11, pad=12, fontweight="bold")
    else:
        ax.set_title(title, fontsize=12, pad=10)

    # Save
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=200, bbox_inches="tight")
    fig.savefig(str(out_path).replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.batch and args.pairs:
        jobs = []
        for triplet in args.pairs:
            parts = triplet.split(":")
            if len(parts) != 3:
                print(f"ERROR: bad triplet '{triplet}', expected CLADE:SpA:SpB",
                      file=sys.stderr)
                continue
            jobs.append(tuple(parts))
    elif args.sp_a and args.sp_b and args.clade:
        jobs = [(args.clade, args.sp_a, args.sp_b)]
    else:
        print("ERROR: provide --sp_a/--sp_b/--clade or --batch --pairs", file=sys.stderr)
        sys.exit(1)

    # Load TE content per species
    te_data_file = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                        "figure_species_data.tsv")
    te_pct = {}
    if te_data_file.exists():
        with open(te_data_file) as fh:
            next(fh)
            for line in fh:
                c = line.strip().split("\t")
                if len(c) >= 2:
                    try:
                        te_pct[c[0]] = float(c[1])
                    except ValueError:
                        pass

    # Load recency index per species
    ki_file = Path("/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/"
                   "C01_kimura_species.tsv")
    ri_vals = {}
    if ki_file.exists():
        with open(ki_file) as fh:
            for line in fh:
                c = line.strip().split("\t")
                if len(c) >= 7 and c[0] != "recency_index":
                    try:
                        ri_vals[c[6]] = float(c[0])
                    except ValueError:
                        pass

    # Load patristic distances AND corrected r from corrected file
    corrected_file = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny_guide_corrected/"
                          "corrected_synteny_ratios.tsv")
    distances = {}
    scores = {}   # use r_corrected_norm (col 12, 0-based) as the displayed r
    if corrected_file.exists():
        with open(corrected_file) as fh:
            next(fh)
            for line in fh:
                c = line.strip().split("\t")
                if len(c) >= 13:
                    key = tuple(sorted([c[1], c[2]]))
                    try:
                        distances[key] = float(c[9])   # d_ij
                        scores[key]    = float(c[12])  # r_corrected_norm
                    except ValueError:
                        pass

    for clade, sp_a, sp_b in jobs:
        print(f"\n[Circos] {sp_a} × {sp_b} (clade {clade})")

        genome_a = load_genome(args.bed_dir, sp_a)
        genome_b = load_genome(args.bed_dir, sp_b)
        pos_a = load_bed(args.bed_dir, sp_a)
        pos_b = load_bed(args.bed_dir, sp_b)

        raw_pairs = load_collinear(args.collinear_dir, clade, sp_a, sp_b)
        if not raw_pairs:
            print(f"  No collinear pairs found, skipping")
            continue

        filt, scafs_a, scafs_b = select_scaffolds(
            raw_pairs, pos_a, pos_b, args.max_scaffolds, args.min_links,
            genome_a=genome_a, genome_b=genome_b,
            top_by_size=args.top_by_size)

        if not filt:
            print(f"  No links after filtering, skipping")
            continue

        print(f"  {len(raw_pairs)} raw pairs → {len(filt)} filtered")
        print(f"  {len(scafs_a)} scaffolds A, {len(scafs_b)} scaffolds B")

        score_key = tuple(sorted([sp_a, sp_b]))
        score = scores.get(score_key)
        dist = distances.get(score_key)

        # Load TE density for outer track
        print(f"  Loading TE density...")
        te_dens_a = load_te_density(sp_a, scafs_a, genome_a,
                                    split_classes=args.class_split)
        te_dens_b = load_te_density(sp_b, scafs_b, genome_b,
                                    split_classes=args.class_split)

        out_path = out_dir / f"circos_{sp_a}__{sp_b}.pdf"
        make_circos_plot(sp_a, sp_b, filt, scafs_a, scafs_b,
                        genome_a, genome_b, out_path, synteny_score=score,
                        patristic_dist=dist,
                        te_pct_a=te_pct.get(sp_a),
                        te_pct_b=te_pct.get(sp_b),
                        ri_a=ri_vals.get(sp_a),
                        ri_b=ri_vals.get(sp_b),
                        simple_style=args.simple_style,
                        n_total_pairs=len(raw_pairs),
                        te_density_a=te_dens_a,
                        te_density_b=te_dens_b)

    print("\n[DONE]")


if __name__ == "__main__":
    main()
