#!/usr/bin/env python3

#Summarize OrthoFinder gene duplication events per node of the species tree.


import argparse
from pathlib import Path

import pandas as pd
from ete3 import Tree  # requires ete3 installed in the environment


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize OrthoFinder duplications per node")
    parser.add_argument("--orthofinder-dir", required=True, help="Root directory of OrthoFinder run")
    parser.add_argument("--out-dir", required=True, help="Output directory for summary tables")
    parser.add_argument(
        "--gene-trees-dir",
        required=False,
        help=(
            "Optional directory with gene trees. If not provided, the script will try "
            "<orthofinder-dir>/Gene_Trees and <orthofinder-dir>/Resolved_Gene_Trees in this order."
        ),
    )
    return parser.parse_args()


def load_species_tree(orthofinder_dir: Path) -> Tree:
    tree_path = orthofinder_dir / "Species_Tree" / "SpeciesTree_rooted.txt"
    if not tree_path.is_file():
        raise FileNotFoundError(f"Species tree not found at {tree_path}")
    with tree_path.open() as fh:
        newick = fh.read().strip()
    return Tree(newick, format=1)


def load_duplications(orthofinder_dir: Path) -> pd.DataFrame:
    dup_path = orthofinder_dir / "Gene_Duplication_Events" / "Duplications.tsv"
    if not dup_path.is_file():
        raise FileNotFoundError(f"Duplications.tsv not found at {dup_path}")

    df = pd.read_csv(dup_path, sep="\t")

    # Try to standardize key columns if present; we avoid depending on
    # any previous custom parser and work directly with OrthoFinder output.
    expected_cols = {"Node", "Species Tree Node", "Gene Tree", "Orthogroup"}
    missing = expected_cols.intersection(df.columns.symmetric_difference(expected_cols))
    # We do not fail hard on column-name variations, but we record them.

    return df


def parse_species_from_gene_token(token: str) -> str:
    """Extract the species identifier from an OrthoFinder gene token.

    OrthoFinder duplication tables usually store genes as strings like
    'CH_Clapol1_aa_jgi|Clapol1_1|655623|...'. We normalize the leading
    species prefix so the output uses the same identifiers as the rest
    of the pipeline.
    """

    token = str(token).strip()
    if not token:
        return ""

    if "_aa_" in token:
        return token.split("_aa_", 1)[0] + ".aa"
    if ".aa_" in token:
        return token.split(".aa_", 1)[0] + ".aa"

    prefix = token.split("_jgi|", 1)[0]
    if prefix.endswith("_aa"):
        prefix = prefix[:-3] + ".aa"
    elif "|" in prefix:
        prefix = prefix.split("|", 1)[0]
    return prefix


def parse_species_set(cell: str) -> set[str]:
    cell = str(cell).strip()
    if not cell or cell.lower() == "nan":
        return set()
    species: set[str] = set()
    for token in cell.split(","):
        token = token.strip()
        if not token:
            continue
        species_id = parse_species_from_gene_token(token)
        if species_id:
            species.add(species_id)
    return species


def count_gene_tokens(cell: str) -> int:
    cell = str(cell).strip()
    if not cell or cell.lower() == "nan":
        return 0
    return sum(1 for token in cell.split(",") if token.strip())


def classify_node(node_id: str) -> str:
    node_id = str(node_id)
    return "internal" if node_id.startswith("N") else "terminal"


def summarize_per_node(df_dup: pd.DataFrame) -> pd.DataFrame:
    """Summarize number of duplications per species-tree node.

    We rely on the column that encodes the species-tree node; in most
    recent OrthoFinder versions this is 'Species Tree Node' or 'Node'.
    """

    node_col = None
    for cand in ("Species Tree Node", "Node", "Species_Tree_Node"):
        if cand in df_dup.columns:
            node_col = cand
            break
    if node_col is None:
        raise ValueError(
            "Could not find a species-tree node column in Duplications.tsv; "
            "looked for 'Species Tree Node', 'Node', 'Species_Tree_Node'."
        )

    group = df_dup.groupby(node_col)
    records = []
    for node_id, sub in group:
        node_type = classify_node(node_id)

        species_set_1: set[str] = set()
        species_set_2: set[str] = set()
        n_genes_1 = 0
        n_genes_2 = 0
        for cell in sub.get("Genes 1", pd.Series(dtype=str)).dropna().astype(str):
            n_genes_1 += count_gene_tokens(cell)
            species_set_1.update(parse_species_set(cell))
        for cell in sub.get("Genes 2", pd.Series(dtype=str)).dropna().astype(str):
            n_genes_2 += count_gene_tokens(cell)
            species_set_2.update(parse_species_set(cell))

        rec = {
            "species_tree_node_id": node_id,
            "node_type": node_type,
            "n_duplications": len(sub),
            "n_genes_side1": n_genes_1,
            "n_genes_side2": n_genes_2,
            "n_species_side1": len(species_set_1),
            "n_species_side2": len(species_set_2),
            "n_species_union": len(species_set_1 | species_set_2),
        }
        # crude normalization: duplications per distinct gene tree / orthogroup if available
        if "Gene Tree" in sub.columns:
            rec["n_gene_trees"] = sub["Gene Tree"].nunique()
            rec["dup_per_gene_tree"] = rec["n_duplications"] / rec["n_gene_trees"]
        if "Orthogroup" in sub.columns:
            rec["n_orthogroups"] = sub["Orthogroup"].nunique()
            rec["dup_per_orthogroup"] = rec["n_duplications"] / rec["n_orthogroups"]

        records.append(rec)

    return pd.DataFrame.from_records(records).sort_values("n_duplications", ascending=False)


def summarize_per_species_branch(df_dup: pd.DataFrame, tree: Tree) -> pd.DataFrame:
    """Summarize duplication events per species found in the duplicated genes.

    OrthoFinder's duplication table contains the reconciliation node plus two
    gene lists. The raw node ids are not always species-tree leaf names, so we
    derive a species-level summary directly from Genes 1/2.
    """

    node_col = None
    for cand in ("Species Tree Node", "Node", "Species_Tree_Node"):
        if cand in df_dup.columns:
            node_col = cand
            break
    if node_col is None:
        raise ValueError(
            "Could not find a species-tree node column in Duplications.tsv; "
            "looked for 'Species Tree Node', 'Node', 'Species_Tree_Node'."
        )

    records = []
    for node_id, sub in df_dup.groupby(node_col):
        species_set: set[str] = set()
        for gene_col in ("Genes 1", "Genes 2"):
            if gene_col in sub.columns:
                for cell in sub[gene_col].dropna().astype(str):
                    species_set.update(parse_species_set(cell))

        for sp in sorted(species_set):
            records.append({"species": sp, "node_id": node_id, "n_duplications": len(sub)})

    if not records:
        return pd.DataFrame(columns=["species", "node_id", "n_duplications"])

    df = pd.DataFrame.from_records(records)
    df = (
        df.groupby(["species", "node_id"], as_index=False)["n_duplications"]
        .sum()
        .sort_values(["species", "n_duplications"], ascending=[True, False])
    )
    return df


def detect_wgd(
    orthofinder_dir: Path,
    tree: Tree,
    df_dup: pd.DataFrame,
    out_dir: Path,
    gene_trees_dir: Path | None = None,
) -> pd.DataFrame:
    """Detect WGD signals using 2:1 triplets and OrthoFinder duplications.

    Esta versión sigue el enfoque de tripletes 2:1 pero usa la reconciliación
    de OrthoFinder (Duplications.tsv + árbol de especies) para validar que la
    duplicación ocurrió en la rama de la especie con 2 copias después del
    split con la especie con 1 copia. Los árboles de genes se usan de forma
    indirecta a través de Duplications.tsv, que ya incorpora esa información.

    Salidas:
      - wgd_triplets_per_species_pairs.tsv: conteo por par de especies (A,B).
      - wgd_triplets_per_species.tsv: resumen por especie A (con 2 copias).
    """

    og_path = orthofinder_dir / "Orthogroups" / "Orthogroups.tsv"
    if not og_path.is_file():
        raise FileNotFoundError(f"Orthogroups.tsv not found at {og_path}")

    # Tabla de orthogrupos estándar de OrthoFinder.
    df_og = pd.read_csv(og_path, sep="\t").fillna("")
    if "Orthogroup" not in df_og.columns:
        raise ValueError("Expected an 'Orthogroup' column in Orthogroups.tsv")

    species_cols = [c for c in df_og.columns if c != "Orthogroup"]
    if not species_cols:
        raise ValueError("No species columns found in Orthogroups.tsv")

    def parse_gene_list(cell: str) -> list[str]:
        cell = str(cell).strip()
        if not cell or cell.lower() == "nan":
            return []
        tokens = cell.replace(";", ",").replace(" ", ",").split(",")
        return [g for g in (t.strip() for t in tokens) if g]

    # Conteo de genes por orthogrupo y especie.
    counts: dict[str, pd.Series] = {}
    for sp_col in species_cols:
        counts[sp_col] = df_og[sp_col].apply(parse_gene_list).apply(len)

    df_counts = pd.DataFrame(counts)
    total_genes = df_counts.sum(axis=1)

    # Orthogrupos con exactamente 3 genes en total.
    mask_total3 = total_genes == 3
    df_og3 = df_og[mask_total3].reset_index(drop=True)
    df_counts3 = df_counts[mask_total3].reset_index(drop=True)

    # Número total de orthogrupos donde cada especie tiene >=1 gen.
    n_og_per_species: dict[str, int] = {}
    for sp_col in species_cols:
        n_og_per_species[sp_col] = int((df_counts[sp_col] > 0).sum())

    # Información de duplicaciones por orthogrupo y nodo del árbol de especies.
    node_col = None
    for cand in ("Species Tree Node", "Node", "Species_Tree_Node"):
        if cand in df_dup.columns:
            node_col = cand
            break

    has_og_col = "Orthogroup" in df_dup.columns
    dup_by_og = df_dup.groupby("Orthogroup") if has_og_col else None

    # Índice de nodos del árbol de especies para poder comprobar qué especies
    # están bajo cada nodo.
    name_to_node = {n.name: n for n in tree.traverse() if n.name}

    triplet_records: list[dict[str, str]] = []

    for idx, row in df_og3.iterrows():
        og_id = str(row["Orthogroup"])
        counts_row = df_counts3.loc[idx]

        species_with2 = [sp for sp in species_cols if counts_row[sp] == 2]
        species_with1 = [sp for sp in species_cols if counts_row[sp] == 1]
        if not species_with2 or not species_with1:
            continue

        for sp2 in species_with2:
            for sp1 in species_with1:
                is_valid = False

                if has_og_col and node_col is not None and dup_by_og is not None:
                    try:
                        sub_dup = dup_by_og.get_group(og_id)
                    except KeyError:
                        sub_dup = None

                    if sub_dup is not None and not sub_dup.empty:
                        try:
                            mrca = tree.get_common_ancestor(sp2, sp1)
                        except Exception:
                            mrca = None

                        for _, dup_row in sub_dup.iterrows():
                            node_id = str(dup_row[node_col])
                            tnode = name_to_node.get(node_id)
                            if tnode is None:
                                continue

                            leaves = {leaf.name for leaf in tnode.get_leaves()}
                            if sp2 in leaves and sp1 not in leaves:
                                if mrca is not None:
                                    mrca_leaves = {leaf.name for leaf in mrca.get_leaves()}
                                    if not leaves.issubset(mrca_leaves):
                                        continue
                                is_valid = True
                                break
                else:
                    # Sin información de duplicaciones por nodo/orthogrupo,
                    # aceptamos el triplete 2:1 como candidato.
                    is_valid = True

                if not is_valid:
                    continue

                triplet_records.append(
                    {
                        "orthogroup": og_id,
                        "species_two": sp2,
                        "species_one": sp1,
                    }
                )

    # Si no hay ningún triplete válido, escribimos ficheros vacíos.
    if not triplet_records:
        df_pairs_empty = pd.DataFrame(
            columns=[
                "species_two",
                "species_one",
                "n_triplets_2_1",
                "n_orthogroups",
                "ratio_triplets_per_orthogroup",
            ]
        )
        df_species_empty = pd.DataFrame(
            columns=["species", "n_triplets_2_1", "n_orthogroups", "ratio_triplets_per_orthogroup"]
        )
        df_pairs_empty.to_csv(out_dir / "wgd_triplets_per_species_pairs.tsv", sep="\t", index=False)
        df_species_empty.to_csv(out_dir / "wgd_triplets_per_species.tsv", sep="\t", index=False)
        return df_species_empty

    df_trip = pd.DataFrame.from_records(triplet_records)

    # Conteo por par de especies (A con 2 copias, B con 1 copia).
    df_pairs = (
        df_trip.groupby(["species_two", "species_one"])["orthogroup"]
        .nunique()
        .reset_index(name="n_triplets_2_1")
    )

    df_pairs["n_orthogroups"] = df_pairs["species_two"].map(n_og_per_species).fillna(0).astype(int)
    df_pairs["ratio_triplets_per_orthogroup"] = df_pairs.apply(
        lambda r: r["n_triplets_2_1"] / r["n_orthogroups"] if r["n_orthogroups"] > 0 else 0.0,
        axis=1,
    )
    df_pairs = df_pairs.sort_values("ratio_triplets_per_orthogroup", ascending=False)

    df_pairs.to_csv(out_dir / "wgd_triplets_per_species_pairs.tsv", sep="\t", index=False)

    # Resumen por especie (sumando sobre todos los B posibles).
    df_species = (
        df_pairs.groupby("species_two")["n_triplets_2_1"]
        .sum()
        .reset_index(name="n_triplets_2_1")
        .rename(columns={"species_two": "species"})
    )
    df_species["n_orthogroups"] = df_species["species"].map(n_og_per_species).fillna(0).astype(int)
    df_species["ratio_triplets_per_orthogroup"] = df_species.apply(
        lambda r: r["n_triplets_2_1"] / r["n_orthogroups"] if r["n_orthogroups"] > 0 else 0.0,
        axis=1,
    )
    df_species = df_species.sort_values("ratio_triplets_per_orthogroup", ascending=False)

    df_species.to_csv(out_dir / "wgd_triplets_per_species.tsv", sep="\t", index=False)

    return df_species


def main() -> None:
    args = parse_args()

    orthofinder_dir = Path(args.orthofinder_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    tree = load_species_tree(orthofinder_dir)
    df_dup = load_duplications(orthofinder_dir)

    df_node = summarize_per_node(df_dup)
    df_node.to_csv(out_dir / "duplications_per_node.tsv", sep="\t", index=False)

    df_branch = summarize_per_species_branch(df_dup, tree)
    df_branch.to_csv(out_dir / "duplications_per_species_branch.tsv", sep="\t", index=False)

    # Nuevo análisis: detección de señales de WGD mediante tripletes 2:1.
    gene_trees_dir = Path(args.gene_trees_dir).resolve() if getattr(args, "gene_trees_dir", None) else None
    df_wgd = detect_wgd(orthofinder_dir, tree, df_dup, out_dir, gene_trees_dir)
    print(f"[WGD] Wrote 2:1 triplet summary to {out_dir / 'wgd_triplets_per_species.tsv'}")
    if not df_wgd.empty:
        print("[WGD] Top species by 2:1 triplet ratio:")
        # Mostramos un pequeño resumen en stdout para inspección rápida.
        print(df_wgd.head().to_string(index=False))

    print(f"[WGD] Wrote per-node summary to {out_dir / 'duplications_per_node.tsv'}")
    print(f"[WGD] Wrote per-species-branch summary to {out_dir / 'duplications_per_species_branch.tsv'}")


if __name__ == "__main__":
    main()
