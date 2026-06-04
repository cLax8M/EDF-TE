#!/usr/bin/env python3
"""
fig_H_widespread_cooption.py
Find co-option architectures that are highly frequent across fungi/all species,
then for each selected architecture extract protein sequences for phylogenetic analysis.

Logic:
  1. From figC_domain_tepair_summary.tsv (host-domain × TE-superfam pairs with
     n_genes and n_species), select architectures present in ≥MIN_SPECIES species
     AND spanning ≥MIN_CLADES different clades.
  2. For selected architectures, scan per_species classification to get gene IDs.
  3. Find protein FASTA files and extract sequences.
  4. Write one FASTA per selected architecture for phylogenetic analysis.
  5. Also write a SETMAR check (SET domain + Transposase/MULE).

Selected architectures (based on data):
  A) Chromo + Helitron  : 151 species, present in ALL clades (most widespread)
  B) Chromo + RT        : 15 species, across 5 clades (Gypsy-ERV-like, example shown)
  C) RVP (retroviral protease) alone : check n clades
  D) DUF4217 + Helitron : 156 species, highly abundant
  Plus: user-defined via --arch argument

Outputs:
  figH_arch_summary.tsv                — candidate architectures with n_sp, n_clades
  figH_{ARCH_LABEL}_candidates.tsv     — gene table per architecture
  figH_{ARCH_LABEL}_seqs.faa           — sequences for phylogeny
"""
import os, sys, glob, re, collections
import pandas as pd

# ─── Config ──────────────────────────────────────────────────────────────────
COOPTION_DIR  = "/legserv/Temp/Carlos/EDF-TE/results/Co-option"
TE_MASTER     = f"{COOPTION_DIR}/config/TE_domains_master.tsv"
FUSIONS_GLOB  = f"{COOPTION_DIR}/per_species/fusions/*.TE_host_fusions.tsv"
CLASS_GLOB    = f"{COOPTION_DIR}/per_species/classification/*.true_TE_coopted_genes.tsv"
OUT_DIR       = f"{COOPTION_DIR}/figures"
CLADE_FILE    = "/legserv/Temp/Carlos/EDF-TE/results/HTT/clade_assignment.tsv"
PAIR_SUMMARY  = f"{OUT_DIR}/figC_domain_tepair_summary.tsv"

# Protein FASTA search paths (in order of priority)
PROTEIN_DIRS = [
    "/legserv/Temp/Carlos/EDF-TE/data/proteins_annot_faa",
    "/legserv/Temp/Carlos/EDF-TE/data/proteins",
    "/data/carlos/EDF-TE/TE_analysis/annotated",
    "/legserv/Temp/Carlos/EDF-TE/data",
]

MIN_SPECIES = 20      # minimum species for "widespread"
MIN_CLADES  = 3       # minimum clades for "widespread"
MAX_PER_CLADE = 65    # max representative sequences per clade in output FASTA

# Clades ordered by phylogenetic interest (Dikarya first, then early-diverging)
PRIORITY_CLADES = [
    "AS", "BA",             # Dikarya
    "HO", "NU",             # non-fungi outgroups
    "MI", "ZO", "KX",      # early-diverging fungi
    "MO", "MU", "GL",      # Mucoromycota
    "EN", "CH", "BL",      # other early
    "AP", "BS", "CR", "OL", "NE",
]

CLADE_LABELS = {
    "AP":"Aphelida","AS":"Ascomycota","BA":"Basidiomycota",
    "BL":"Blastocladiomycota","BS":"Basidiobolomycota","CH":"Chytridiomycota",
    "CR":"Cryptomycota","EN":"Entomophthoromycota","GL":"Glomeromycota",
    "HO":"Holozoa","KX":"Kickxellomycotina","MI":"Microsporidia",
    "MO":"Mortierellomycota","MU":"Mucoromycotina","NU":"Nuclearida",
    "OL":"Olpidiomycota","ZO":"Zoopagomycota",
}

# ─── TE master ────────────────────────────────────────────────────────────────
te_master = pd.read_csv(TE_MASTER, sep="\t")
te_core_ids  = set(te_master.loc[te_master["TE_status"]=="TE_core","domain_id"].astype(str))
te_assoc_ids = set(te_master.loc[te_master["TE_status"]=="TE_associated","domain_id"].astype(str))
te_all_ids   = te_core_ids | te_assoc_ids
id2name = dict(zip(te_master["domain_id"].astype(str), te_master["domain_name"].astype(str)))
id2sf   = dict(zip(te_master["domain_id"].astype(str), te_master["TE_superfamily"].astype(str)))

# Domain name → Pfam IDs
# Key TE domain patterns:
HELITRON_NAMES = re.compile(r"[Hh]elitron|SNF2|DEAD.helic|Helicase", re.I)
RT_NAMES       = re.compile(r"[Rr]everse.?[Tt]rans|RVT|[Gg]ypsy", re.I)
CHROMO_NAMES   = re.compile(r"[Cc]hromo|CHD|CHROMO", re.I)
RVP_NAMES      = re.compile(r"[Rr]etroviral.*[Pp]rote|[Pp]epsin.*RVP|^RVP$", re.I)
SET_NAMES      = re.compile(r"^SET domain|Su\(var\)|[Ll]ysine.*methyltrans.*SET|SET and RING", re.I)
TRANSPOSASE_NAMES = re.compile(r"[Tt]ransposase|MULE|MuDR|DDE.*transpos", re.I)
# Top-3 empirical architectures (by n_species in EDF dataset):
HELICASE_NAMES = re.compile(r"[Hh]elicas|DEAD|DEAH|SNF2|DExH|HA2|OB.fold|Oligonucl", re.I)
DUF_NAMES      = re.compile(r"DUF4217|[Ii]nositol.poly|[Tt]ype III restrict", re.I)
APE_NAMES      = re.compile(r"AP.endonuclease|AP_endonuclease|Endo.IV|XthA|exodeoxyribonuclease III|ExoIII|xth", re.I)
EXODNASE_NAMES = re.compile(r"exodeoxyribonuclease III|ExoDNase|Xth|xthA|exonuclease III|endonuclease", re.I)
LRR_NAMES      = re.compile(r"[Ll]eucine.rich.repeat|LRR.domain|LRR[_ ]|Leu-rich", re.I)
CHROMO_RT_NAMES= re.compile(r"[Cc]hromodomain|chromo.domain", re.I)
RHOGAP_NAMES   = re.compile(r"RhoGAP|Rho.GTPase.activat|GTPase.activat.*protein|RhoGAP.domain", re.I)
DUF4208_NAMES  = re.compile(r"DUF4208", re.I)

# ─── Clade info ──────────────────────────────────────────────────────────────
clades_df = pd.read_csv(CLADE_FILE, sep="\t", header=None, names=["species","clade"])
sp2clade = dict(zip(clades_df["species"], clades_df["clade"]))

# ─── Load pair summary (from figC) to pick architectures ─────────────────────
pairs = pd.read_csv(PAIR_SUMMARY, sep="\t")
pairs["clade_list"] = None  # will fill from scan below

# Empirical architectures ranked by n_species across EDF dataset:
#   1. Helicase+Helitron  — 199 sp (pure Helitron helicase domestication)
#   2. Chromo+Helitron    — 194 sp (chromodomain + Helitron helicase)
#   3. DUF4217+LINE       — 154 sp (LINE AP-endo/DUF domain)
#   4. APE+LINE           — 263 sp (AP endonuclease from LINE EN domain)
#   5. ExoDNase+LINE      — 240 sp (exodeoxyribonuclease III from LINE)
#   6. LRR+LINE           — 204 sp (leucine-rich repeat + LINE)
#   7. Chromodomain+LINE  — 102 sp (chromodomain + LINE RT, distinct from #2)
#   8. RhoGAP+LINE        —  96 sp (RhoGAP domain + LINE)
ARCH_DEFINITIONS = [
    # (host_re, te_sf_str_match, label, min_sp_override)
    (HELICASE_NAMES,  "Helitron",    "Helicase_Helitron",  30),  # #1: 199 sp
    (CHROMO_NAMES,    "Helitron",    "Chromo_Helitron",    30),  # #2: 194 sp
    (DUF_NAMES,       "LINE",        "DUF4217_LINE",        20),  # #3: 154 sp
    (APE_NAMES,       "LINE",        "APE_LINE",            20),  # #4: 263 sp
    (EXODNASE_NAMES,  "LINE",        "ExoDNase_LINE",       20),  # #5: 240 sp
    (LRR_NAMES,       "LINE",        "LRR_LINE",            20),  # #6: 204 sp
    (CHROMO_RT_NAMES, "LINE",        "Chromo_LINE",         20),  # #7: 102 sp
    (RHOGAP_NAMES,    "LINE",        "RhoGAP_LINE",         15),  # #8:  96 sp
    (DUF4208_NAMES,   "Helitron",    "DUF4208_Helitron",     3),  # fungal-specific: 118 sp, 11 clades
]

# From pairs, compute per-architecture n_species and infer n_clades (approximate)
arch_summary = []
for host_re, te_sf, label, min_sp in ARCH_DEFINITIONS:
    if host_re is not None and te_sf is not None:
        mask = (
            pairs["host_domain"].str.contains(host_re, regex=True, na=False) &
            pairs["te_superfam"].str.contains(te_sf, na=False)
        )
    elif host_re is not None:
        mask = pairs["host_domain"].str.contains(host_re, regex=True, na=False)
    elif te_sf is not None:
        mask = pairs["te_superfam"].str.contains(te_sf, na=False)
    else:
        mask = pd.Series([False]*len(pairs), index=pairs.index)

    sub = pairs[mask]
    total_sp = sub["n_species"].sum()
    total_genes = sub["n_genes"].sum()
    arch_summary.append({
        "label": label, "n_species_est": total_sp, "n_genes_est": total_genes,
        "min_sp": min_sp, "host_re": str(host_re.pattern), "te_sf": te_sf or "any"
    })

df_arch = pd.DataFrame(arch_summary)
print("[INFO] Architecture candidates:")
print(df_arch[["label","n_species_est","n_genes_est"]].to_string(index=False))
df_arch.to_csv(f"{OUT_DIR}/figH_arch_summary.tsv", sep="\t", index=False)

# Filter to architectures meeting threshold
selected = [
    r for _, r in df_arch.iterrows()
    if r["n_species_est"] >= r["min_sp"]
]
print(f"\n[INFO] Selected architectures: {[r['label'] for r in selected]}")

# ─── Load true co-opted genes ─────────────────────────────────────────────────
true_genes = collections.defaultdict(set)
for f in sorted(glob.glob(CLASS_GLOB)):
    sp = os.path.basename(f).split(".")[0]
    df = pd.read_csv(f, sep="\t", usecols=["gene_id"])
    true_genes[sp].update(df["gene_id"].astype(str))

# ─── Scan fusions: gene-level domain matches per architecture ─────────────────
arch_gene_dicts = {r["label"]: [] for r in selected}  # label → list of row dicts

host_re_map = {r["label"]: re.compile(r["host_re"], re.I) for r in selected}
te_sf_map   = {r["label"]: r["te_sf"] for r in selected}

print("[INFO] Scanning per-species fusions...")
for fpath in sorted(glob.glob(FUSIONS_GLOB)):
    sp = os.path.basename(fpath).split(".")[0]
    true_sp = true_genes.get(sp, set())
    if not true_sp:
        continue
    clade = sp2clade.get(sp, sp.split("_")[0])

    df = pd.read_csv(fpath, sep="\t", dtype=str).fillna("")
    df_true = df[df["gene_id"].isin(true_sp)]

    for gid, gdf in df_true.groupby("gene_id"):
        domain_ids   = set(gdf["domain_id"].astype(str))
        domain_names = " ".join(gdf["domain_name"].astype(str).str.strip())
        te_sfs_gene  = set()
        host_names_gene = set()

        for _, row in gdf.iterrows():
            did   = str(row["domain_id"])
            dname = str(row["domain_name"]).strip()
            dclass = str(row.get("domain_class",""))
            sf = id2sf.get(did, "")
            if did in te_core_ids or dclass == "TE_core":
                te_sfs_gene.add(sf if sf else "Unknown")
            elif did not in te_assoc_ids:
                if dname not in ("nan","-",""):
                    host_names_gene.add(dname)

        for label, h_re in host_re_map.items():
            te_sf_required = te_sf_map[label]
            # Check host domain match
            has_host = any(h_re.search(n) for n in host_names_gene) if host_names_gene else \
                       h_re.search(domain_names)
            # Check TE superfam match
            has_te = True
            if te_sf_required != "any" and te_sf_required is not None:
                has_te = any(te_sf_required.lower() in s.lower() for s in te_sfs_gene) or \
                         bool(re.search(te_sf_required, domain_names, re.I))

            if has_host and has_te:
                arch_gene_dicts[label].append({
                    "species": sp, "clade": clade,
                    "clade_label": CLADE_LABELS.get(clade, clade),
                    "gene_id": gid,
                    "host_domains": ";".join(sorted(host_names_gene))[:100],
                    "te_superfamilies": ";".join(sorted(te_sfs_gene)),
                })

# ─── Save per-architecture gene tables ───────────────────────────────────────
for label, rows in arch_gene_dicts.items():
    if not rows:
        print(f"[SKIP] {label}: no genes found")
        continue
    df_out = pd.DataFrame(rows)
    n_sp  = df_out["species"].nunique()
    n_cl  = df_out["clade"].nunique()
    print(f"[INFO] {label}: {len(df_out)} genes, {n_sp} species, {n_cl} clades")
    print("  Clades:", sorted(df_out["clade"].unique()))
    tsv_out = f"{OUT_DIR}/figH_{label}_candidates.tsv"
    df_out.to_csv(tsv_out, sep="\t", index=False)
    print(f"  → {tsv_out}")

# ─── Extract protein sequences ────────────────────────────────────────────────
import gzip as gz

def find_protein_fasta(sp):
    sp2 = sp.split("_")[1] if "_" in sp else sp
    for d in PROTEIN_DIRS:
        for sub in ("", sp, sp2):
            base = os.path.join(d, sub)
            for fname in (f"{sp}.faa", f"{sp}.fasta", f"{sp}.fa", f"{sp}.pep",
                          f"{sp}.aa.fasta", f"{sp}.aa.fa",
                          f"{sp2}.faa", f"{sp2}.fasta"):
                p = os.path.join(base, fname) if sub else os.path.join(d, fname)
                if os.path.exists(p):
                    return p
            # glob
            for pat in (f"{sp}*.faa", f"{sp}*.fasta", f"{sp2}*.faa"):
                hits = glob.glob(os.path.join(base if sub else d, pat))
                if hits:
                    return hits[0]
    return None


def read_fasta(path):
    """Parse FASTA and index sequences by multiple possible IDs.

    Handles JGI-format headers (>jgi|SPECIES|PROTEINID|...) by also indexing
    the protein ID (3rd pipe field) so integer gene_ids from the candidates
    TSV can be found even when the full header is longer.
    Also indexes by 'proteinId=N' if present in the header line.
    """
    seqs = {}
    cur_ids, cur_seq = [], []
    opener = gz.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_ids:
                    seq = "".join(cur_seq)
                    for k in cur_ids:
                        seqs[k] = seq
                header = line[1:]
                primary = header.split()[0]
                ids = [primary]
                # JGI: >jgi|PORTAL|PROTEINID|ANNOTATION
                parts = primary.split("|")
                if len(parts) >= 3:
                    ids.append(parts[2])   # numeric protein ID
                # proteinId=N anywhere in the header
                import re as _re
                m = _re.search(r"proteinId=(\d+)", header)
                if m:
                    ids.append(m.group(1))
                cur_ids = list(dict.fromkeys(ids))  # deduplicate, preserve order
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_ids:
        seq = "".join(cur_seq)
        for k in cur_ids:
            seqs[k] = seq
    return seqs


# Sample test for protein FASTA availability
test_sp = list(true_genes.keys())[:5]
found_count = 0
for tsp in test_sp:
    fp = find_protein_fasta(tsp)
    status = "OK" if fp else "NOT FOUND"
    print(f"  [FASTA check] {tsp}: {status}{' -> '+fp if fp else ''}")
    if fp:
        found_count += 1

if found_count == 0:
    print("\n[WARN] No protein FASTAs found in default directories.")
    print("  Sequences will not be extracted. Check PROTEIN_DIRS paths.")
    print("  Candidate gene tables have been saved for manual sequence extraction.")
else:
    # Extract sequences per architecture — stratified sampling by clade
    for label, rows in arch_gene_dicts.items():
        if not rows:
            continue
        df_cand = pd.DataFrame(rows)
        faa_out = f"{OUT_DIR}/figH_{label}_seqs.faa"
        n_written = 0

        # Stratified sampling: up to MAX_PER_CLADE sequences per clade,
        # one per species within each clade, following PRIORITY_CLADES order.
        # Clades not in PRIORITY_CLADES are appended at the end.
        seen_clades = set(df_cand["clade"].unique())
        clade_order = [c for c in PRIORITY_CLADES if c in seen_clades] + \
                      [c for c in seen_clades if c not in PRIORITY_CLADES]

        sampled = []  # list of (sp, gene_id, clade) tuples
        for clade in clade_order:
            clade_df = df_cand[df_cand["clade"] == clade]
            n_clade  = 0
            for sp in clade_df["species"].unique():
                if n_clade >= MAX_PER_CLADE:
                    break
                sp_row = clade_df[clade_df["species"] == sp].iloc[0]
                sampled.append((sp, str(sp_row["gene_id"]), clade))
                n_clade += 1

        print(f"  [{label}] {len(sampled)} seqs selected across "
              f"{len(set(c for _,_,c in sampled))} clades")

        with open(faa_out, "w") as fout:
            for sp, gid, clade in sampled:
                fasta_path = find_protein_fasta(sp)
                if not fasta_path:
                    continue
                try:
                    seqs = read_fasta(fasta_path)
                    seq  = seqs.get(gid)
                    if seq and len(seq) > 50:
                        safe_gid = re.sub(r"[^A-Za-z0-9._-]", "_", gid)
                        fout.write(f">{sp}|{clade}|{safe_gid}\n{seq}\n")
                        n_written += 1
                except Exception:
                    pass

        if n_written > 0:
            print(f"[OK] {label}: {n_written} sequences → {faa_out}")
        else:
            print(f"[WARN] {label}: 0 sequences extracted (FASTAs may not match gene IDs)")
            os.remove(faa_out)

print("\n[DONE] Fig H prep complete.")
