#!/usr/bin/env python3
"""
04_build_pa_matrix.py — Build cross-species TE family presence/absence matrix
                        from MMseqs2/cd-hit cluster output + per-species FL TSVs.

Inputs:
  --cluster_tsv    MMseqs2 *_cluster.tsv (rep_id<TAB>member_id) where IDs are
                   formatted "{SPECIES}|{family_name}" (from step 02).
  --fl_dir         per-species *_fl_families.tsv (from step 01)
  --whitelist      species whitelist (one species_prefix per line; pairs/species
                   not in it are dropped)  [optional]

Output:
  {out_dir}/te_family_pa_matrix.tsv     rows=cluster_id, cols=species, vals=0/1
  {out_dir}/te_family_pa_matrix.count   Count v9.1 input format (taxa table)
  {out_dir}/cluster_members.tsv         cluster_id<TAB>member1,member2,...
  {out_dir}/summary.txt                 stats: n_clusters, n_species, fill rate

The Count format is:
  Family  sp1  sp2  sp3 ...
  fam1    0    1    1
  ...

Usage:
  python 04_build_pa_matrix.py \
      --cluster_tsv /path/to/cluster.tsv \
      --fl_dir /path/to/fl_families \
      --whitelist /path/to/species_whitelist.txt \
      --out_dir /path/to/out

Author: EDF-TE project
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict

# Species names in cluster/FL files that differ from tree tip names
_NAME_MAP = {
    "AS_AspfumUD": "AS_AspfumUD1",
    "NU_Parvatl":  "NU_Parvatl_c",
    "HO_Magisepte": "HO_Magisepte_c",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cluster_tsv", required=True)
    p.add_argument("--fl_dir", required=True)
    p.add_argument("--whitelist", default=None,
                   help="(deprecated; ignored if --tree given)")
    p.add_argument("--tree", default=None,
                   help="Newick tree file; all tip names become columns")
    p.add_argument("--out_dir", required=True)
    return p.parse_args()


def load_whitelist(path):
    if not path:
        return None
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def load_tree_tips(tree_path):
    """Extract tip names from a Newick tree file."""
    with open(tree_path) as fh:
        nwk = fh.read()
    return set(re.findall(r'[\w.]+(?=:)', nwk))


def load_fl_set(fl_dir, species):
    tsv = Path(fl_dir) / f"{species}_fl_families.tsv"
    if not tsv.exists():
        return set()
    out = set()
    with open(tsv) as fh:
        next(fh)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 4:
                try:
                    if int(cols[3]) >= 1:
                        out.add(cols[0])
                except ValueError:
                    pass
    return out


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Determine the full species set from tree (preferred) or whitelist
    tree_tips = load_tree_tips(args.tree) if args.tree else None
    whitelist = load_whitelist(args.whitelist) if not tree_tips else None
    # Build reverse name map for cluster species → tree tip name
    rev_map = {v: v for v in (tree_tips or set())}
    for old, new in _NAME_MAP.items():
        if tree_tips and new in tree_tips:
            rev_map[old] = new

    # cluster_id -> list of (species, family)
    cl_members = defaultdict(list)
    members_seen = set()
    species_present = set()

    with open(args.cluster_tsv) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            rep, mem = cols[0], cols[1]
            if mem in members_seen:
                continue
            members_seen.add(mem)
            try:
                sp_raw, fam = mem.split("|", 1)
            except ValueError:
                continue
            sp = rev_map.get(sp_raw, sp_raw)
            if whitelist is not None and sp not in whitelist:
                continue
            if tree_tips is not None and sp not in tree_tips:
                continue
            cl_members[rep].append((sp, fam))
            species_present.add(sp)

    # The final species list = all tree tips (including those with no clusters)
    if tree_tips:
        all_species = sorted(tree_tips)
    else:
        all_species = sorted(species_present)
    print(f"Loaded {len(cl_members)} clusters across {len(species_present)} species "
          f"(matrix will have {len(all_species)} columns)")

    # For each species, FL-present family set
    fl_sets = {}
    for sp in all_species:
        # Try direct name first, then reverse-mapped name
        fl_sets[sp] = load_fl_set(args.fl_dir, sp)
        if not fl_sets[sp]:
            for old, new in _NAME_MAP.items():
                if new == sp:
                    fl_sets[sp] = load_fl_set(args.fl_dir, old)
                    break

    # Build matrix
    sp_idx = {sp: i for i, sp in enumerate(all_species)}
    matrix = {}
    for cl, mems in cl_members.items():
        row = [0] * len(all_species)
        for sp, fam in mems:
            if fam in fl_sets.get(sp, set()):
                row[sp_idx[sp]] = 1
        matrix[cl] = row

    # Filter clusters present in at least 1 species (should be all by construction)
    informative_cls = [cl for cl, row in matrix.items() if sum(row) >= 1]
    print(f"Informative clusters (≥1 species): {len(informative_cls)}")

    # Write standard TSV matrix
    pa_file = out_dir / "te_family_pa_matrix.tsv"
    with open(pa_file, "w") as fh:
        fh.write("cluster_id\t" + "\t".join(all_species) + "\n")
        for cl in informative_cls:
            fh.write(cl + "\t" + "\t".join(str(v) for v in matrix[cl]) + "\n")

    # Write Count v9.1 input format (same as TSV with "Family" header)
    count_file = out_dir / "te_family_pa_matrix.count"
    with open(count_file, "w") as fh:
        fh.write("Family\t" + "\t".join(all_species) + "\n")
        for cl in informative_cls:
            fh.write(cl + "\t" + "\t".join(str(v) for v in matrix[cl]) + "\n")

    # Write cluster members
    members_file = out_dir / "cluster_members.tsv"
    with open(members_file, "w") as fh:
        fh.write("cluster_id\tn_members\tmembers\n")
        for cl in informative_cls:
            mems = cl_members[cl]
            fh.write(f"{cl}\t{len(mems)}\t" +
                     ",".join(f"{s}|{f}" for s, f in mems) + "\n")

    # Summary
    n_cells = len(informative_cls) * len(all_species)
    n_present = sum(sum(matrix[cl]) for cl in informative_cls)
    fill = n_present / n_cells if n_cells else 0
    summary_file = out_dir / "summary.txt"
    with open(summary_file, "w") as fh:
        fh.write(f"clusters\t{len(informative_cls)}\n")
        fh.write(f"species\t{len(all_species)}\n")
        fh.write(f"presence_cells\t{n_present}\n")
        fh.write(f"fill_rate\t{fill:.4f}\n")
    print(f"Matrix:  {pa_file}")
    print(f"Count:   {count_file}")
    print(f"Members: {members_file}")
    print(f"Summary: {summary_file} (fill rate {fill:.2%})")


if __name__ == "__main__":
    main()
