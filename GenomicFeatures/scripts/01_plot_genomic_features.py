#!/usr/bin/env python3
"""
Plot genomic features across 238 species grouped by clade.

Usage:
  python3 01_plot_genomic_features.py <ALL_SPECIES_ALL_SUMMARY.txt> <out_dir>

Outputs (PDF + PNG):
  - boxplot per metric (Total_Genes, Avg_Exons_Per_Gene, Avg_Introns_Multiexon,
    Avg_Isoforms_Per_Gene, Pct_Exonic, Pct_Intergenic)
  - stacked bar: mean Pct_Exonic / Pct_Intronic / Pct_Intergenic per clade
  - heatmap: 238 species x 6 metrics (z-score, ordered by clade)
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Phylogenetic order (canonical, includes NE)
CLADE_ORDER = ['HO', 'NU', 'AP', 'CR', 'MI', 'BL', 'OL', 'CH', 'NE',
               'BS', 'EN', 'ZO', 'KX', 'GL', 'MO', 'MU', 'BA', 'AS']

# Canonical clade palette (single source of truth)
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                   '..', '..', 'Section2_TElandscape', 'scripts'))
from clade_palette import CLADE_COLORS, CLADE_NAMES, get_clade  # noqa: E402

# (column, y-axis label, log_scale)
METRICS_BOX = [
    ('Total_Genes',            'Total genes',                   True),
    ('Avg_Exons_Per_Gene',     'Avg exons per gene',            False),
    ('Avg_Introns_Multiexon',  'Avg introns / gene (multiexon)', False),
    ('Avg_Isoforms_Per_Gene',  'Avg isoforms per gene',         False),
    ('Pct_Exonic',             'Exonic (%)',                    False),
    ('Pct_Intergenic',         'Intergenic (%)',                False),
]

HEATMAP_COLS = ['Total_Genes', 'Avg_Exons_Per_Gene', 'Avg_Introns_Multiexon',
                'Avg_Isoforms_Per_Gene', 'Pct_Exonic', 'Pct_Intergenic']
HEATMAP_LABELS = ['Total genes', 'Exons/gene', 'Introns/gene\n(multiexon)',
                  'Isoforms/gene', 'Exonic %', 'Intergenic %']


def load_data(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    df['Species'] = df['Species'].str.replace(r'\.aa$', '', regex=True)
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df['Clade'] = df['Species'].apply(get_clade)
    df = df[df['Clade'].isin(CLADE_ORDER)].copy()
    df['Clade'] = pd.Categorical(df['Clade'], categories=CLADE_ORDER, ordered=True)
    return df.sort_values('Clade').reset_index(drop=True)


def boxplot_metric(df, col, label, logscale, out_dir):
    present = [c for c in CLADE_ORDER if c in df['Clade'].values]
    data    = [df[df['Clade'] == c][col].dropna().values for c in present]
    colors  = [CLADE_COLORS.get(c, '#ccc') for c in present]

    fig, ax = plt.subplots(figsize=(14, 5))
    bp = ax.boxplot(data, patch_artist=True,
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(linewidth=0.8),
                    flierprops=dict(marker='o', markersize=3, alpha=0.5))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)

    ax.set_xticks(range(1, len(present) + 1))
    ax.set_xticklabels(present, fontsize=9)
    ax.set_xlabel('Clade', fontsize=11)
    ax.set_ylabel(label, fontsize=11)
    ax.set_title(f'{label}  by clade  (n={len(df)} species)', fontsize=12)
    if logscale:
        ax.set_yscale('log')
    ax.grid(axis='y', alpha=0.3)
    sns.despine(ax=ax)

    ymin = ax.get_ylim()[0]
    for i, (c, v) in enumerate(zip(present, data)):
        ax.annotate(f'n={len(v)}', xy=(i + 1, ymin),
                    ha='center', va='bottom', fontsize=7, color='#666')

    plt.tight_layout()
    out = os.path.join(out_dir, f'boxplot_{col}.pdf')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    fig.savefig(out.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  OK  boxplot_{col}")


def plot_stacked_bar(df, out_dir):
    cols    = ['Pct_Exonic', 'Pct_Intronic', 'Pct_Intergenic']
    labels  = ['Exonic', 'Intronic', 'Intergenic']
    bcolors = ['#2166ac', '#92c5de', '#d1e5f0']
    present = [c for c in CLADE_ORDER if c in df['Clade'].values]
    means   = df.groupby('Clade', observed=True)[cols].mean().reindex(present)

    fig, ax = plt.subplots(figsize=(14, 5))
    bottom  = np.zeros(len(present))
    for col, label, color in zip(cols, labels, bcolors):
        vals = means[col].values
        ax.bar(present, vals, bottom=bottom, label=label,
               color=color, edgecolor='white', linewidth=0.5)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 4:
                ax.text(i, b + v / 2, f'{v:.1f}%', ha='center', va='center',
                        fontsize=7, color='white' if v > 8 else 'black')
        bottom += vals

    ax.set_xlabel('Clade', fontsize=11)
    ax.set_ylabel('Mean genomic composition (%)', fontsize=11)
    ax.set_title('Genomic composition by clade  (mean per species)', fontsize=12)
    ax.set_ylim(0, 105)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    sns.despine(ax=ax)
    plt.tight_layout()

    out = os.path.join(out_dir, 'stacked_bar_genomic_composition.pdf')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    fig.savefig(out.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  OK  stacked_bar_genomic_composition")


def plot_heatmap(df, out_dir):
    mat = df.set_index('Species')[HEATMAP_COLS].dropna()
    z   = (mat - mat.mean()) / mat.std()

    clade_map  = df.set_index('Species')['Clade'].reindex(z.index)
    row_colors = clade_map.map(CLADE_COLORS).fillna('#ccc')

    g = sns.clustermap(
        z,
        row_cluster=False, col_cluster=False,
        row_colors=row_colors,
        cmap='RdBu_r', center=0, vmin=-3, vmax=3,
        yticklabels=True, xticklabels=HEATMAP_LABELS,
        linewidths=0,
        figsize=(10, max(12, len(z) * 0.065)),
        cbar_kws={'label': 'Z-score', 'shrink': 0.3}
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=4)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                  fontsize=9, rotation=30, ha='right')
    g.fig.suptitle('Genomic features — 238 Holomycota species', y=1.001, fontsize=12)

    handles = [mpatches.Patch(color=CLADE_COLORS.get(c, '#ccc'), label=c)
               for c in CLADE_ORDER if c in clade_map.values]
    g.fig.legend(handles=handles, title='Clade', loc='lower left',
                 bbox_to_anchor=(0.02, 0.01), fontsize=7, title_fontsize=8,
                 ncol=3, framealpha=0.8)

    out = os.path.join(out_dir, 'heatmap_genomic_features.pdf')
    g.fig.savefig(out, dpi=150, bbox_inches='tight')
    g.fig.savefig(out.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
    plt.close(g.fig)
    print(f"  OK  heatmap_genomic_features")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    tsv_path, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    print(f"\nLoading {tsv_path}")
    df = load_data(tsv_path)
    print(f"  {len(df)} species | {df['Clade'].nunique()} clades")

    print("\nBoxplots...")
    for col, label, logscale in METRICS_BOX:
        if col in df.columns:
            boxplot_metric(df, col, label, logscale, out_dir)

    print("\nStacked bar...")
    plot_stacked_bar(df, out_dir)

    print("\nHeatmap...")
    plot_heatmap(df, out_dir)

    print(f"\nDone — plots in {out_dir}")


if __name__ == '__main__':
    main()
