#!/usr/bin/env python3
"""
fig_H_gypsy_chromo_phylotree.py
Phylogenetic tree of Gypsy/ERV chromovirus co-opted genes.

Step 1 (this script): Find genes with Chromo + RT domains across species,
  extract protein sequences, write multi-species FASTA.

Step 2 (cluster): mafft + iqtree (via SLURM).

Step 3 (visualize): ggtree R script (fig_H_plot_phylo.R).

Outputs:
  figH_gypsy_chromo_candidates.tsv   — per-gene table
  figH_gypsy_chromo_seqs.faa         — multi-FASTA for phylogeny
  figH_run_phylo.sh                  — SLURM script for alignment + tree
"""
import os, sys, glob, re, collections
import pandas as pd

# ─── Paths ──────────────────────────────────────────────────────────────────
COOPTION_DIR  = "/legserv/Temp/Carlos/EDF-TE/results/Co-option"
TE_MASTER     = f"{COOPTION_DIR}/config/TE_domains_master.tsv"
FUSIONS_GLOB  = f"{COOPTION_DIR}/per_species/fusions/*.TE_host_fusions.tsv"
CLASS_GLOB    = f"{COOPTION_DIR}/per_species/classification/*.true_TE_coopted_genes.tsv"
PROTEIN_DIR   = "/legserv/Temp/Carlos/EDF-TE/data/proteins"
OUT_DIR       = f"{COOPTION_DIR}/figures"
CLADE_FILE    = "/legserv/Temp/Carlos/EDF-TE/results/HTT/clade_assignment.tsv"

CLADE_LABELS = {
    "AP":"Aphelida","AS":"Ascomycota","BA":"Basidiomycota",
    "BL":"Blastocladiomycota","BS":"Basidiobolomycota","CH":"Chytridiomycota",
    "CR":"Cryptomycota","EN":"Entomophthoromycota","GL":"Glomeromycota",
    "HO":"Holozoa","KX":"Kickxellomycotina","MI":"Microsporidia",
    "MO":"Mortierellomycota","MU":"Mucoromycotina","NU":"Nuclearida",
    "OL":"Olpidiomycota","ZO":"Zoopagomycota",
}

# ─── TE master: which domains define Gypsy/ERV and Chromo ────────────────────
te_master = pd.read_csv(TE_MASTER, sep="\t")

# Chromo domain IDs
CHROMO_PATTERN = re.compile(r"[Cc]hromo|CHD|chromovir", re.IGNORECASE)
chromo_ids = set(te_master[
    te_master["domain_name"].astype(str).str.contains(CHROMO_PATTERN, regex=True) |
    te_master["domain_id"].astype(str).str.contains(CHROMO_PATTERN, regex=True)
]["domain_id"].astype(str))
# Also add known Chromo Pfam IDs
chromo_ids.update({"PF00385", "PS50013", "PS00537"})  # Chromo, chromo shadow

# RT / Gypsy domain IDs
RT_PATTERN = re.compile(r"[Rr]everse_?[Tt]rans|[Gg]ypsy|RVT|RNH|RdRp", re.IGNORECASE)
rt_ids = set(te_master[
    te_master["domain_name"].astype(str).str.contains(RT_PATTERN, regex=True) |
    te_master["domain_id"].astype(str).str.contains(RT_PATTERN, regex=True) |
    te_master["TE_superfamily"].astype(str).str.contains(r"Gypsy|Reverse_trans|LTR/Gypsy", regex=True)
]["domain_id"].astype(str))
# Add known RT Pfam IDs
rt_ids.update({"PF00078","PF07727","PF13456","PF17917","PF00075"})  # RVT_1/2, RNH

print(f"[INFO] Chromo domain IDs: {len(chromo_ids)} — {sorted(chromo_ids)[:8]}")
print(f"[INFO] RT/Gypsy domain IDs: {len(rt_ids)} — {sorted(rt_ids)[:8]}")

# ─── Clade info ──────────────────────────────────────────────────────────────
clades = pd.read_csv(CLADE_FILE, sep="\t", header=None, names=["species","clade"])
sp2clade = dict(zip(clades["species"], clades["clade"]))

# ─── Load true co-opted genes ────────────────────────────────────────────────
true_genes = collections.defaultdict(set)
for f in sorted(glob.glob(CLASS_GLOB)):
    sp = os.path.basename(f).split(".")[0]
    df = pd.read_csv(f, sep="\t", usecols=["gene_id"])
    true_genes[sp].update(df["gene_id"].astype(str))

# ─── Find genes with Chromo + RT in true co-opted set ────────────────────────
candidates = []
print("[INFO] Scanning fusion files for Chromo+RT co-opted genes...")

for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    true_sp = true_genes.get(sp, set())
    if not true_sp:
        continue
    clade = sp2clade.get(sp, sp.split("_")[0])

    df = pd.read_csv(fpath, sep="\t", dtype=str).fillna("")
    df_true = df[df["gene_id"].isin(true_sp)]

    for gid, gdf in df_true.groupby("gene_id"):
        domain_ids = set(gdf["domain_id"].astype(str))
        has_chromo = bool(domain_ids & chromo_ids)
        has_rt     = bool(domain_ids & rt_ids)
        # Also match by domain name
        domain_names = " ".join(gdf["domain_name"].astype(str))
        if not has_chromo:
            has_chromo = bool(CHROMO_PATTERN.search(domain_names))
        if not has_rt:
            has_rt = bool(RT_PATTERN.search(domain_names))

        if has_chromo and has_rt:
            candidates.append({
                "species": sp, "clade": clade,
                "clade_label": CLADE_LABELS.get(clade, clade),
                "gene_id": gid,
                "domains": ";".join(sorted(domain_ids)),
                "has_chromo": has_chromo,
                "has_rt": has_rt,
            })

df_cand = pd.DataFrame(candidates)
if len(df_cand) == 0:
    print("[WARN] No genes with Chromo+RT found!")
else:
    print(f"[INFO] Found {len(df_cand)} genes with Chromo+RT in {df_cand['species'].nunique()} species")
    print(df_cand.groupby("clade_label")["gene_id"].count().sort_values(ascending=False).to_string())

    out_tsv = f"{OUT_DIR}/figH_gypsy_chromo_candidates.tsv"
    df_cand.to_csv(out_tsv, sep="\t", index=False)
    print(f"[OK] {out_tsv}")

# ─── Try to extract protein sequences ────────────────────────────────────────
# Look for protein FASTA files per species
print(f"\n[INFO] Looking for protein FASTAs in {PROTEIN_DIR} ...")

def find_protein_fasta(sp, search_dirs):
    """Find protein FASTA for species sp."""
    for d in search_dirs:
        for ext in (".faa", ".fasta", ".fa", ".pep", ".aa.fasta"):
            for prefix in (sp, sp.split("_")[1] if "_" in sp else sp):
                p = os.path.join(d, f"{prefix}{ext}")
                if os.path.exists(p):
                    return p
                p2 = os.path.join(d, sp, f"{prefix}{ext}")
                if os.path.exists(p2):
                    return p2
    return None

search_dirs = [
    PROTEIN_DIR,
    "/legserv/Temp/Carlos/EDF-TE/data/proteins",
    "/data/carlos/EDF-TE/TE_analysis/annotated",
    "/legserv/Temp/Carlos/EDF-TE/data",
]

# Sample a few species to check
if len(df_cand) > 0:
    test_sps = df_cand["species"].unique()[:3]
    for tsp in test_sps:
        fp = find_protein_fasta(tsp, search_dirs)
        print(f"  {tsp}: {fp or 'NOT FOUND'}")

# ─── Build the FASTA (if protein files exist) ─────────────────────────────────
import gzip as gz

def parse_fasta(fasta_path):
    """Return dict gene_id → sequence."""
    seqs = {}
    cur_id = None
    cur_seq = []
    opener = gz.open if str(fasta_path).endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)
    return seqs


out_faa = f"{OUT_DIR}/figH_gypsy_chromo_seqs.faa"
n_written = 0

if len(df_cand) > 0:
    with open(out_faa, "w") as fout:
        for sp, sp_grp in df_cand.groupby("species"):
            clade = sp_grp["clade"].iloc[0]
            fasta_path = find_protein_fasta(sp, search_dirs)
            if not fasta_path:
                continue
            try:
                seqs = parse_fasta(fasta_path)
                for _, row in sp_grp.iterrows():
                    gid = str(row["gene_id"])
                    seq = seqs.get(gid)
                    if seq:
                        # Header: species|clade|gene_id for tree labeling
                        fout.write(f">{sp}|{clade}|{gid}\n{seq}\n")
                        n_written += 1
            except Exception as e:
                print(f"  [WARN] {sp}: {e}")

    print(f"[INFO] Sequences written to {out_faa}: {n_written}")
else:
    n_written = 0

# ─── Write SLURM script for phylogenetic analysis ────────────────────────────
slurm_script = f"""\
#!/bin/bash
#SBATCH --job-name=figH_phylo
#SBATCH --partition=normal.168h
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/data/carlos/EDF-TE/logs/figH_phylo_%j.out
#SBATCH --error=/data/carlos/EDF-TE/logs/figH_phylo_%j.err

set -e
source ~/.bashrc
micromamba activate esptreeEDFTE

OUTDIR="{OUT_DIR}"
FAA="{out_faa}"
PREFIX="${{OUTDIR}}/figH_gypsy_chromo"

echo "[1] MAFFT alignment..."
mafft --auto --thread 16 "${{FAA}}" > "${{PREFIX}}_aligned.faa"

echo "[2] Trimming with trimAl..."
if command -v trimal &>/dev/null; then
    trimal -in "${{PREFIX}}_aligned.faa" -out "${{PREFIX}}_trimmed.faa" -automated1
    ALN_FILE="${{PREFIX}}_trimmed.faa"
else
    ALN_FILE="${{PREFIX}}_aligned.faa"
fi

echo "[3] IQ-TREE..."
iqtree2 -s "${{ALN_FILE}}" \\
  -m TEST -B 1000 -T 16 \\
  --prefix "${{PREFIX}}_iqtree" \\
  --seed 42

echo "[DONE] Tree: ${{PREFIX}}_iqtree.treefile"
"""

slurm_path = f"{OUT_DIR}/figH_run_phylo.sh"
with open(slurm_path, "w") as fh:
    fh.write(slurm_script)
os.chmod(slurm_path, 0o755)
print(f"[OK] SLURM script: {slurm_path}")
print("\n[DONE] Fig H step 1 complete.")
print(f"  If n_written={n_written}: sbatch {slurm_path}")
print("  Then run fig_H_plot_phylo.R after tree is built.")
