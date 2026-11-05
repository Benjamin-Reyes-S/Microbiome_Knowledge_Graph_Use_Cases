"""Run the queries of kg_queries.py and write one table per entity type: KG_human_{reactions,genes,metabolites}.csv.
Reactions and metabolites: one row per KG node. Genes: one row per KG Gene node for the model's gene ids only (the expanded KG has
~8.4M Gene nodes), so this script needs Human_GEM_genes.csv from 02 (organizer.sh runs 02 first).
One boolean column per route plus its n_support count.
New KG identifiers (BiGG, MetaNetX, ChEBI): every node row keeps source_db and the arrays bigg_ids / metanetx_ids (';'-joined) and, for
compounds, kegg_id / chebi_id.  Because the arrays can hold several ids and an id can sit on several nodes, the exploded one-row-per-
(node, id) tables KG_xref_reactions.csv and KG_xref_metabolites.csv are written as well (03 should match model ids against these)."""
import pandas as pd
import common as c
import kg_queries as q

H = c.HUMAN


def run(section, name, **p):
    rows = c.cypher(section[name], **p)
    print(f"{name:28s} {len(rows):>8} rows", flush=True)
    return rows


def add_routes(df, section, routes, **extra):
    """routes: {column_prefix: query_name}; adds <prefix> (bool) and <prefix>_n (support count)."""
    for col, name in routes.items():
        hit = pd.DataFrame(run(section, name, h=H, **extra))
        if hit.empty:
            hit = pd.DataFrame({"id": [], "n_support": []})
        n = dict(zip(hit.id, hit.n_support)) if "n_support" in hit else {i: 1 for i in hit.id}
        df[col] = df.id.isin(n)
        df[col + "_n"] = df.id.map(n).fillna(0).astype(int)
    return df


def join(x):
    return ";".join(x) if isinstance(x, list) else ""


def write_xref(section, names, fname):
    """Exploded cross-reference table (id, namespace, xref) + n_nodes_with_xref = how many KG nodes carry the same xref (1:n mapping)."""
    x = pd.concat([pd.DataFrame(run(section, n), columns=["id", "namespace", "xref"]) for n in names], ignore_index=True)
    x["n_nodes_with_xref"] = x.groupby(["namespace", "xref"]).id.transform("nunique")
    x.to_csv(c.RESULTS / fname, index=False)
    print(f"{fname}: {len(x)} rows;", x.groupby("namespace").agg(nodes=("id", "nunique"), xrefs=("xref", "nunique"),
          xrefs_on_several_nodes=("n_nodes_with_xref", lambda v: int((v > 1).sum()))).to_dict("index"), flush=True)


# ---- reactions -----------------------------------------------------------------------------------------------------
r = pd.DataFrame(run(q.REACTIONS, "R_NODES"))
r["n_ec"], r["n_ko"] = r.ecs.map(len), r.kos.map(len)
r["ecs"], r["kos"] = r.ecs.map(join), r.kos.map(join)
r["n_bigg"], r["n_metanetx"] = r.bigg_ids.map(lambda x: len(x) if isinstance(x, list) else 0), r.metanetx_ids.map(lambda x: len(x) if isinstance(x, list) else 0)
r["bigg_ids"], r["metanetx_ids"] = r.bigg_ids.map(join), r.metanetx_ids.map(join)
r = add_routes(r, q.REACTIONS, {"R1_enzyme_protein": "R1_enzyme_protein", "R2_ko_protein": "R2_ko_protein",
                                "R3_module_pathway": "R3_module_pathway", "R3b_pathway_direct": "R3b_pathway_direct",
                                "R4_ko_hsa_gene": "R4_ko_hsa_gene"})
r["strong_human"] = r.R1_enzyme_protein | r.R2_ko_protein
r.to_csv(c.RESULTS / "KG_human_reactions.csv", index=False)
write_xref(q.REACTIONS, ["R_XREF_BIGG", "R_XREF_METANETX"], "KG_xref_reactions.csv")

# ---- genes ---------------------------------------------------------------------------------------------------------
mg = c.RESULTS / "Human_GEM_genes.csv"
if not mg.exists():
    raise SystemExit(f"{mg} missing: run 02_extract_human_gem.py first (organizer.sh does)")
gene_ids = sorted({x for s in pd.read_csv(mg).kg_gene.dropna() for x in s.split(";") if x})
print(f"model gene ids to look up in the KG: {len(gene_ids)}")
g = pd.DataFrame(run(q.GENES, "G_NODES", ids=gene_ids))
g["symbols"], g["xrefs"], g["kos"] = g.symbols.map(join), g.xrefs.map(join), g.kos.map(join)
g = add_routes(g, q.GENES, {"G1_uniprot_protein": "G1_uniprot_protein", "G2_ko_protein": "G2_ko_protein",
                            "G3_hsa_prefix": "G3_hsa_prefix"}, ids=gene_ids)
g["edge_type"] = "none (no taxon edge on Gene)"
g.loc[g.G1_uniprot_protein, "edge_type"] = "Protein-BELONGS_TO (via UniProt xref)"
g["strong_human"] = g.G1_uniprot_protein
g.to_csv(c.RESULTS / "KG_human_genes.csv", index=False)
print("direct Gene-Taxon edges (expect none):", run(q.GENES, "G_direct_taxon_edge_check", h=H, ids=gene_ids))

# ---- metabolites ---------------------------------------------------------------------------------------------------
m = pd.DataFrame(run(q.METABOLITES, "M_NODES"))
m["synonyms"] = m.synonyms.map(join)
m["n_bigg"], m["n_metanetx"] = m.bigg_ids.map(lambda x: len(x) if isinstance(x, list) else 0), m.metanetx_ids.map(lambda x: len(x) if isinstance(x, list) else 0)
m["bigg_ids"], m["metanetx_ids"] = m.bigg_ids.map(join), m.metanetx_ids.map(join)
m[["formula", "kegg_id", "chebi_id"]] = m[["formula", "kegg_id", "chebi_id"]].fillna("")
m["in_any_reaction"] = ~m.id.isin({x["id"] for x in run(q.METABOLITES, "M_no_reaction")})
m = add_routes(m, q.METABOLITES, {"M1a_via_R1": "M1a_via_R1", "M1b_via_R2": "M1b_via_R2", "M2_pathway": "M2_pathway",
                                  "M3_found_in": "M3_found_in"})
m["strong_human"] = m.M1a_via_R1 | m.M1b_via_R2 | m.M3_found_in
m.to_csv(c.RESULTS / "KG_human_metabolites.csv", index=False)
write_xref(q.METABOLITES, ["M_XREF_CHEBI", "M_XREF_BIGG", "M_XREF_METANETX"], "KG_xref_metabolites.csv")

for name, d in (("reactions", r), ("genes", g), ("metabolites", m)):
    cols = [x for x in d.columns if d[x].dtype == bool]
    print(f"\nKG {name}: {len(d)} nodes\n", d[cols].sum().to_string())
