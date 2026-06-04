#!/usr/bin/env python3
"""
chromatin_schematic.py
======================
Conceptual diagram explaining NRL, CV (rugosity) and Phasing
for a typical PCG region vs a typical TE region.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

OUT = "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/supervised_chromatin/invivo_landscape"

rng = np.random.default_rng(42)

# ── nucleosome positions ──────────────────────────────────────────────────
def make_positions(n, mean_gap, cv, start=0):
    """Generate n nucleosome centre positions with given mean gap and CV."""
    sd = cv * mean_gap
    gaps = rng.normal(mean_gap, sd, size=n).clip(min=50)
    return np.cumsum(np.concatenate([[start], gaps]))

NUC_W = 147   # nucleosome footprint bp
NUC_H = 0.55  # oval height (arbitrary units)

def draw_chromatin_row(ax, positions, color, label, y=0, linker_col="#cccccc"):
    """Draw nucleosomes as ovals on a DNA line."""
    x_min = positions[0] - NUC_W / 2 - 10
    x_max = positions[-1] + NUC_W / 2 + 10
    ax.plot([x_min, x_max], [y, y], color=linker_col, lw=2.5, zorder=1)
    for xc in positions:
        ellipse = mpatches.Ellipse(
            (xc, y), width=NUC_W, height=NUC_H,
            facecolor=color, edgecolor="white", linewidth=1.5, zorder=2)
        ax.add_patch(ellipse)
    ax.text(x_min - 30, y, label, ha="right", va="center",
            fontsize=10, fontweight="bold", color=color)

def draw_spacing_arrows(ax, positions, y_base, y_arr, color):
    """Draw double-headed arrows for each inter-nucleosome gap."""
    for i in range(len(positions) - 1):
        x0 = positions[i] + NUC_W / 2
        x1 = positions[i + 1] - NUC_W / 2
        if x1 <= x0:
            continue
        ax.annotate("", xy=(x1, y_arr), xytext=(x0, y_arr),
                    arrowprops=dict(arrowstyle="<->", color=color, lw=1.4))
        mid = (x0 + x1) / 2
        gap = positions[i + 1] - positions[i]
        ax.text(mid, y_arr + 0.06, f"{gap:.0f} bp",
                ha="center", va="bottom", fontsize=7, color=color)

def autocorrelation_profile(positions, total_len=2000, lag_max=500):
    """Compute nucleosome occupancy autocorrelation."""
    occ = np.zeros(total_len)
    for xc in positions:
        lo = max(0, int(xc - NUC_W // 2))
        hi = min(total_len, int(xc + NUC_W // 2))
        occ[lo:hi] = 1.0
    lags = np.arange(0, lag_max)
    ac = np.array([np.corrcoef(occ[:-l], occ[l:])[0, 1] if l > 0 else 1.0
                   for l in lags])
    return lags, ac

# ── figure layout ─────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10))
fig.patch.set_facecolor("white")

# GridSpec: 3 rows (chromatin | spacing zoom | autocorrelation)
# × 2 cols (PCG | TE)
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
gs_outer = GridSpec(3, 2, figure=fig,
                    height_ratios=[1.6, 1.0, 1.4],
                    hspace=0.55, wspace=0.35,
                    left=0.12, right=0.97, top=0.93, bottom=0.06)

COL_PCG = "#2a7fb5"
COL_TE  = "#c0392b"

N_NUC   = 9
MEAN_NRL_PCG = 175   # bp between centres
MEAN_NRL_TE  = 210
CV_PCG  = 0.05       # very regular
CV_TE   = 0.20       # irregular

pos_pcg = make_positions(N_NUC, MEAN_NRL_PCG, CV_PCG, start=100)
pos_te  = make_positions(N_NUC, MEAN_NRL_TE,  CV_TE,  start=100)

for col_i, (pos, color, title_lbl, cv_val, nrl_val) in enumerate([
    (pos_pcg, COL_PCG, "Protein-coding gene (PCG)", CV_PCG, MEAN_NRL_PCG),
    (pos_te,  COL_TE,  "Transposable element (TE)",  CV_TE,  MEAN_NRL_TE),
]):
    # ── Row 0: nucleosome diagram ────────────────────────────────────
    ax0 = fig.add_subplot(gs_outer[0, col_i])
    ax0.set_facecolor("white")
    draw_chromatin_row(ax0, pos, color, "", y=0.3)
    draw_chromatin_row(ax0, pos, color, "", y=0.3)   # just once

    # NRL annotation: bracket between first two centres
    nrl_obs = pos[1] - pos[0]
    y_brk = 0.68
    ax0.annotate("", xy=(pos[1], y_brk), xytext=(pos[0], y_brk),
                 arrowprops=dict(arrowstyle="<->", color="#333333", lw=1.5))
    ax0.text((pos[0] + pos[1]) / 2, y_brk + 0.06,
             f"NRL ≈ {nrl_val} bp", ha="center", va="bottom",
             fontsize=9, fontweight="bold", color="#333333")

    # Spacing heterogeneity markers
    if col_i == 1:   # TE: highlight variable gaps with colored brackets
        for i in range(len(pos) - 1):
            gap = pos[i + 1] - pos[i]
            shade = min(1.0, abs(gap - MEAN_NRL_TE) / 80)
            y_b = -0.12 - (i % 2) * 0.08
            ax0.annotate("", xy=(pos[i+1] - NUC_W/2, y_b),
                         xytext=(pos[i] + NUC_W/2, y_b),
                         arrowprops=dict(arrowstyle="<->",
                                         color=plt.cm.Reds(0.4 + shade * 0.5),
                                         lw=1.2))

    ax0.set_xlim(pos[0] - NUC_W, pos[-1] + NUC_W)
    ax0.set_ylim(-0.35, 0.92)
    ax0.axis("off")
    ax0.set_title(title_lbl, fontsize=12, fontweight="bold", color=color, pad=8)

    # ── Row 1: metric summary box ────────────────────────────────────
    ax1 = fig.add_subplot(gs_outer[1, col_i])
    ax1.set_facecolor("#f8f8f8")
    ax1.set_xlim(0, 1); ax1.set_ylim(0, 1)
    ax1.axis("off")

    gaps_val  = np.diff(pos)
    cv_real   = gaps_val.std() / gaps_val.mean()
    nrl_real  = gaps_val.mean() + NUC_W   # centre-to-centre
    # fake phasing for display
    phase_str = "HIGH  →  regular arrays" if col_i == 0 \
                else "LOW  →  disordered arrays"
    cv_str    = f"LOW ({cv_real:.2f})  →  uniform spacing" if col_i == 0 \
                else f"HIGH ({cv_real:.2f})  →  irregular spacing"

    txt = (f"NRL  =  {nrl_real:.0f} bp\n\n"
           f"CV  =  {cv_str}\n\n"
           f"Phasing  =  {phase_str}")
    ax1.text(0.5, 0.52, txt, ha="center", va="center",
             fontsize=9.5, color="#333333",
             bbox=dict(boxstyle="round,pad=0.5",
                       fc="white", ec=color, lw=1.8))

    # ── Row 2: autocorrelation profile ──────────────────────────────
    ax2 = fig.add_subplot(gs_outer[2, col_i])
    ax2.set_facecolor("white")

    total = int(pos[-1] + NUC_W + 200)
    lags, ac = autocorrelation_profile(pos, total_len=total, lag_max=450)
    ax2.plot(lags, ac, color=color, lw=2.0)
    ax2.axhline(0, color="#888888", lw=0.8, ls="--")
    ax2.axvline(0, color="#888888", lw=0.8)

    # annotate first positive peak (NRL)
    peak_idx = np.argmax(ac[50:]) + 50
    ax2.annotate(f"1st peak\n≈ NRL = {lags[peak_idx]} bp",
                 xy=(lags[peak_idx], ac[peak_idx]),
                 xytext=(lags[peak_idx] + 30, ac[peak_idx] + 0.08),
                 arrowprops=dict(arrowstyle="->", color=color, lw=1.2),
                 fontsize=8, color=color)
    # annotate peak height = phasing
    ax2.annotate(f"Peak height\n= Phasing",
                 xy=(lags[peak_idx], ac[peak_idx]),
                 xytext=(lags[peak_idx] + 30, ac[peak_idx] - 0.15),
                 arrowprops=dict(arrowstyle="->", color="#555555", lw=1.0),
                 fontsize=8, color="#555555")

    ax2.set_xlabel("Lag (bp)", fontsize=9)
    ax2.set_ylabel("Autocorrelation", fontsize=9)
    ax2.set_title("Nucleosome occupancy autocorrelation",
                  fontsize=9, color="#555555")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.tick_params(labelsize=8)

# ── global title ─────────────────────────────────────────────────────────
fig.suptitle(
    "Chromatin metrics explained: NRL, CV (rugosity) and Phasing",
    fontsize=13, fontweight="bold", y=0.975)

for ext in ("pdf", "png"):
    kw = {"bbox_inches": "tight"} if ext == "pdf" \
         else {"dpi": 200, "bbox_inches": "tight"}
    fig.savefig(f"{OUT}/chromatin_schematic.{ext}", **kw)
print("Saved chromatin_schematic.pdf / .png")
