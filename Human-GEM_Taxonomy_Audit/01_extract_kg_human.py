"""Run the queries of kg_queries.py and write one table per entity type: KG_human_{reactions,genes,metabolites}.csv.
Reactions and metabolites: one row per KG node. Genes: one row per KG Gene node for the model's gene ids only (the expanded KG has
~8.4M Gene nodes), so this script needs Human_GEM_genes.csv from 02 (organizer.sh runs 02 first).
One boolean column per route plus its n_support count."""
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


# ---- reactions -----------------------------------------------------------------------------------------------------
r = pd.DataFrame(run(q.REACTIONS, "R_NODES"))
r["n_ec"], r["n_ko"] = r.ecs.map(len), r.kos.map(len)
r["ecs"], r["kos"] = r.ecs.map(join), r.kos.map(join)
r = add_routes(r, q.REACTIONS, {"R1_enzyme_protein": "R1_enzyme_protein", "R2_ko_protein": "R2_ko_protein",
                                "R3_module_pathway": "R3_module_pathway", "R4_ko_hsa_gene": "R4_ko_hsa_gene"})
r["strong_human"] = r.R1_enzyme_protein | r.R2_ko_protein
r.to_csv(c.RESULTS / "KG_human_reactions.csv", index=False)

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
m["in_any_reaction"] = ~m.id.isin({x["id"] for x in run(q.METABOLITES, "M_no_reaction")})
m = add_routes(m, q.METABOLITES, {"M1a_via_R1": "M1a_via_R1", "M1b_via_R2": "M1b_via_R2", "M2_pathway": "M2_pathway"})
m["strong_human"] = m.M1a_via_R1 | m.M1b_via_R2
m.to_csv(c.RESULTS / "KG_human_metabolites.csv", index=False)

for name, d in (("reactions", r), ("genes", g), ("metabolites", m)):
    cols = [x for x in d.columns if d[x].dtype == bool]
    print(f"\nKG {name}: {len(d)} nodes\n", d[cols].sum().to_string())
