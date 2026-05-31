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
                     n_orthologs=None):
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

    # Build angular layout: species A on top (right half), species B on bottom (left half)
    arc_map = {}  # scaffold -> (start_angle, end_angle)
    color_map = {}

    angle = sp_gap_frac / 2  # start with half large gap

    # Species A scaffolds
    sp_a_solid = "#e6772a"  # orange
    sp_b_solid = "#2d8b54"  # green
    for i, s in enumerate(scafs_a):
        sz = genome_a.get(s, 0)
        span = (sz / total) * avail
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
        span = (sz / total) * avail
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
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")

    R_inner = 0.90
    R_outer = 1.0

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

    # Draw scaffold labels
    for key, (a0, a1) in arc_map.items():
        span = a1 - a0
        if span > 0.08:
            mid = (a0 + a1) / 2
            sp_tag, scaf = key
            label = scaf if len(scaf) < 18 else scaf[:14] + "…"
            lx = (R_outer + 0.04) * np.cos(mid)
            ly = (R_outer + 0.04) * np.sin(mid)
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
            # Visual emphasis = how *dominant* this scaffold-pair is for both
            # scaffolds. Conserved synteny → one scaffold's genes mostly go to
            # one partner. Fragmented/eroded synteny → many small partners.
            pc = pair_counts.get((scaf_a, scaf_b), 1)
            fa = pc / max(a_scaf_total.get(scaf_a, 1), 1)   # share for sp A scaffold
            fb = pc / max(b_scaf_total.get(scaf_b, 1), 1)   # share for sp B scaffold
            dominance = min(fa, fb)                          # 0–1, conservative
            if dominance >= 0.40 and pc >= 5:
                # strong syntenic block — deep blue
                link_col = "#1f3a8a"
                link_alpha = 0.65
                link_lw = 0.55 + 0.6 * dominance
            elif dominance >= 0.15 and pc >= 2:
                link_col = "#7bb3e0"
                link_alpha = 0.30
                link_lw = 0.35
            else:
                # rearranged / orphan / TE-mediated noise
                link_col = "#cfcfcf"
                link_alpha = 0.12
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
        if scafs_a:
            mid_a = (arc_map[("A", scafs_a[0])][0] + arc_map[("A", scafs_a[-1])][1]) / 2
            lx = 1.2 * np.cos(mid_a)
            ly = 1.2 * np.sin(mid_a)
            ax.text(lx, ly, sp_a_short, fontsize=10, ha="center", va="center",
                    fontweight="bold", style="italic",
                    color=SP_A_COLORS[0])
        if scafs_b:
            mid_b = (arc_map[("B", scafs_b[0])][0] + arc_map[("B", scafs_b[-1])][1]) / 2
            lx = 1.2 * np.cos(mid_b)
            ly = 1.2 * np.sin(mid_b)
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
        ax.set_title(title, fontsize=13, pad=12, fontweight="bold")
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

    # Load synteny ratios and patristic distances for labels
    ratio_file = Path(args.bed_dir).parent / "ratios" / "all_clades_synteny_ratios.tsv"
    scores = {}
    if ratio_file.exists():
        with open(ratio_file) as fh:
            next(fh)
            for line in fh:
                c = line.strip().split("\t")
                if len(c) >= 9:
                    key = tuple(sorted([c[1], c[2]]))
                    scores[key] = float(c[8])

    # Load patristic distances from corrected file
    corrected_file = Path("/legserv/Temp/Carlos/EDF-TE/results/Synteny_guide_corrected/"
                          "corrected_synteny_ratios.tsv")
    distances = {}
    if corrected_file.exists():
        with open(corrected_file) as fh:
            next(fh)
            for line in fh:
                c = line.strip().split("\t")
                if len(c) >= 10:
                    key = tuple(sorted([c[1], c[2]]))
                    try:
                        distances[key] = float(c[9])  # d_ij column
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

        out_path = out_dir / f"circos_{sp_a}__{sp_b}.pdf"
        make_circos_plot(sp_a, sp_b, filt, scafs_a, scafs_b,
                        genome_a, genome_b, out_path, synteny_score=score,
                        patristic_dist=dist,
                        te_pct_a=te_pct.get(sp_a),
                        te_pct_b=te_pct.get(sp_b),
                        ri_a=ri_vals.get(sp_a),
                        ri_b=ri_vals.get(sp_b),
                        simple_style=args.simple_style,
                        n_total_pairs=len(raw_pairs))

    print("\n[DONE]")


if __name__ == "__main__":
    main()
