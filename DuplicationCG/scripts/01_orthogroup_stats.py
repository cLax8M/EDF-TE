#!/usr/bin/env python3
"""
Per-species OrthoFinder duplication / paralog statistics.

Computes for each species:
  N_genes_in_OGs           total genes assigned to orthogroups
  N_paralog_genes          genes in OGs where the species has > 1 copy
  Pct_paralog_genes        N_paralog_genes / N_genes_in_OGs * 100
  N_OGs_present            OGs where species has >= 1 gene
  N_OGs_multicopy          OGs where species has > 1 gene (multi-copy families)
  Pct_OGs_multicopy        N_OGs_multicopy / N_OGs_present * 100
  N_OGs_above_clade_avg    OGs where copy# > mean of own clade
  N_OGs_comparable_clade   OGs used for clade comparison (species present + >=1 clade-mate)
  Pct_OGs_above_clade_avg  N_OGs_above_clade_avg / N_OGs_comparable_clade * 100
  N_OGs_above_fungal_avg   OGs where copy# > mean across all Fungi (all clades except HO)
  N_OGs_comparable_fungal  OGs used for kingdom comparison
  Pct_OGs_above_fungal_avg N_OGs_above_fungal_avg / N_OGs_comparable_fungal * 100
                           [NaN for HO species — they are the non-fungal outgroup]

Usage:
  python3 01_orthogroup_stats.py <Orthogroups.GeneCount.tsv> <out_tsv>
"""

import sys
import os
import pandas as pd
import numpy as np

CLADE_ORDER = ['CR', 'OL', 'ZO', 'GL', 'BL', 'EN', 'KX', 'NU',
               'BS', 'BA', 'AS', 'CH', 'MI', 'MO', 'MU', 'HO', 'AP']

NON_FUNGAL = {'HO'}   # clades excluded from the fungal kingdom reference


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    in_tsv, out_tsv = sys.argv[1], sys.argv[2]
    print(f"\nLoading {in_tsv} ...")
    df = pd.read_csv(in_tsv, sep='\t', index_col=0)

    # Drop summary row / Total column if present
    if 'Total' in df.columns:
        df = df.drop(columns=['Total'])
    df.columns = [c.replace('.aa', '') for c in df.columns]
    print(f"  {len(df):,} orthogroups x {len(df.columns)} species")

    species_list  = list(df.columns)
    clade_of      = {sp: sp[:2] for sp in species_list}

    # Group species by clade
    clade_species = {}
    for sp, cl in clade_of.items():
        clade_species.setdefault(cl, []).append(sp)

    # Fungal reference: all species NOT in NON_FUNGAL clades
    fungal_sp  = [sp for sp in species_list if clade_of[sp] not in NON_FUNGAL]
    fungal_df  = df[fungal_sp]
    # Pre-compute row-wise sum and presence count across all fungi
    fungal_sum = fungal_df.sum(axis=1)
    fungal_n   = (fungal_df > 0).sum(axis=1)
    print(f"  Fungal kingdom reference: {len(fungal_sp)} species "
          f"(excluded: {NON_FUNGAL})")

    records = []
    for i, sp in enumerate(species_list):
        clade        = clade_of[sp]
        col          = df[sp]
        present_mask = col > 0
        n_present    = int(present_mask.sum())
        n_genes      = int(col[present_mask].sum())

        # --- Multi-copy / paralog genes ---
        mc_mask         = col > 1
        n_mc            = int(mc_mask.sum())
        n_paralog_genes = int(col[mc_mask].sum()) if n_mc > 0 else 0
        pct_mc          = round(n_mc / n_present * 100, 3)           if n_present > 0 else 0
        pct_para        = round(n_paralog_genes / n_genes * 100, 3)  if n_genes   > 0 else 0

        # --- Clade-level above-average ---
        clade_others = [s for s in clade_species.get(clade, []) if s != sp]
        if clade_others:
            cdf           = df[clade_others]
            mask_c        = (cdf > 0).any(axis=1) & present_mask
            clade_mean    = cdf[mask_c].mean(axis=1)
            n_above_clade = int((col[mask_c] > clade_mean).sum())
            n_comp_clade  = int(mask_c.sum())
            pct_above_cl  = round(n_above_clade / n_comp_clade * 100, 3) if n_comp_clade else 0
        else:
            n_above_clade = n_comp_clade = 0
            pct_above_cl  = float('nan')

        # --- Fungal kingdom above-average (only for fungal species) ---
        if clade not in NON_FUNGAL:
            # Subtract self from precomputed fungal totals to avoid self-comparison
            if sp in fungal_sp:
                adj_sum = fungal_sum - col
                adj_n   = fungal_n   - present_mask.astype(int)
            else:
                adj_sum = fungal_sum
                adj_n   = fungal_n
            mask_f      = (adj_n > 0) & present_mask
            fungal_mean = adj_sum[mask_f] / adj_n[mask_f]
            n_above_f   = int((col[mask_f] > fungal_mean).sum())
            n_comp_f    = int(mask_f.sum())
            pct_above_f = round(n_above_f / n_comp_f * 100, 3) if n_comp_f else 0
        else:
            n_above_f = n_comp_f = 0
            pct_above_f = float('nan')

        records.append({
            'Species':                  sp,
            'Clade':                    clade,
            'N_genes_in_OGs':           n_genes,
            'N_paralog_genes':          n_paralog_genes,
            'Pct_paralog_genes':        pct_para,
            'N_OGs_present':            n_present,
            'N_OGs_multicopy':          n_mc,
            'Pct_OGs_multicopy':        pct_mc,
            'N_OGs_above_clade_avg':    n_above_clade,
            'N_OGs_comparable_clade':   n_comp_clade,
            'Pct_OGs_above_clade_avg':  pct_above_cl,
            'N_OGs_above_fungal_avg':   n_above_f,
            'N_OGs_comparable_fungal':  n_comp_f,
            'Pct_OGs_above_fungal_avg': pct_above_f,
        })

        if (i + 1) % 25 == 0 or (i + 1) == len(species_list):
            print(f"  [{i+1}/{len(species_list)}] {sp}: "
                  f"para={pct_para:.1f}%  mc={pct_mc:.1f}%  "
                  f">clade={pct_above_cl:.1f}%  "
                  f">fungi={'N/A' if np.isnan(pct_above_f) else f'{pct_above_f:.1f}%'}")

    out_df = pd.DataFrame(records)
    out_df['Clade'] = pd.Categorical(out_df['Clade'],
                                     categories=CLADE_ORDER, ordered=True)
    out_df = out_df.sort_values(['Clade', 'Species']).reset_index(drop=True)

    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    out_df.to_csv(out_tsv, sep='\t', index=False)
    print(f"\nDone — {len(out_df)} species written to {out_tsv}")

    # Quick summary
    print("\nPreview:")
    print(out_df[['Species', 'Pct_paralog_genes', 'Pct_OGs_multicopy',
                  'Pct_OGs_above_clade_avg', 'Pct_OGs_above_fungal_avg']
                ].head(10).to_string(index=False))


if __name__ == '__main__':
    main()
