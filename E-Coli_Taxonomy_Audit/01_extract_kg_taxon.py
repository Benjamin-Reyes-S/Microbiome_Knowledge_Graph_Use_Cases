"""Run the queries of kg_queries.py and write the KG tables for the model's target taxon (E. coli K-12 MG1655, scope in common.py):
KG_taxon_{reactions,genes,proteins,metabolites}.csv.
Reactions and metabolites: one row per KG node. Genes: KG Gene nodes for the model's Entrez ids only, plus KG Protein nodes for the
model's UniProt accessions (the expanded KG has ~8.4M Gene nodes), so this script needs Model_genes.csv from 02 (organizer.sh runs 02 first).
One boolean column per route (strict scope: strain + K-12) plus <route>_n (support count) and <route>_species (wider scope: + species)."""
import pandas as pd
import common as c
import kg_queries as q

P = dict(taxa=c.TAXA, taxa_sp=c.TAXA_SPECIES)


def run(section, name, **p):
    rows = c.cypher(section[name], **P, **p)
    print(f"{name:28s} {len(rows):>8} rows", flush=True)
    return rows


def add_routes(df, section, routes, key="id", **extra):
    """routes: {column: query_name}; adds <col> (bool, strict scope), <col>_n (strict support count), <col>_species (bool, wider scope)."""
    for col, name in routes.items():
        hit = pd.DataFrame(run(section, name, **extra), columns=["id", "n_support", "n_support_species"])
        n, ns = dict(zip(hit.id, hit.n_support)), dict(zip(hit.id, hit.n_support_species))
        df[col] = df[key].map(n).fillna(0).astype(int) > 0
        df[col + "_n"] = df[key].map(n).fillna(0).astype(int)
        df[col + "_species"] = df[key].map(ns).fillna(0).astype(int) > 0
    return df


def join(x):
    return ";".join(x) if isinstance(x, list) else ""


# ---- reactions -----------------------------------------------------------------------------------------------------
r = pd.DataFrame(run(q.REACTIONS, "R_NODES"))
r["n_ec"], r["n_ko"] = r.ecs.map(len), r.kos.map(len)
r["ecs"], r["kos"] = r.ecs.map(join), r.kos.map(join)
r = add_routes(r, q.REACTIONS, {"R1_enzyme_protein": "R1_enzyme_protein", "R2_ko_protein": "R2_ko_protein",
                                "R3_module_pathway": "R3_module_pathway"})
r["strong_taxon"] = r.R1_enzyme_protein | r.R2_ko_protein
r["strong_taxon_species"] = r.R1_enzyme_protein_species | r.R2_ko_protein_species
r.to_csv(c.RESULTS / "KG_taxon_reactions.csv", index=False)

# ---- genes ---------------------------------------------------------------------------------------------------------
mg = c.RESULTS / "Model_genes.csv"
if not mg.exists():
    raise SystemExit(f"{mg} missing: run 02_extract_model.py first (organizer.sh does)")
mg = pd.read_csv(mg).fillna("")
split = lambda s: [x for x in s.split(";") if x]
gene_ids = sorted({x for s in mg.kg_gene for x in split(s)})
uniprots = sorted({x for s in mg.uniprot for x in split(s)})
locus = sorted({x for s in mg.locus for x in split(s)})
print(f"model ids to look up: {len(gene_ids)} NCBIGene, {len(uniprots)} UniProt, {len(locus)} locus tags")

g = pd.DataFrame(run(q.GENES, "G_NODES", gene_ids=gene_ids), columns=["id", "symbols", "xrefs"])
g["symbols"], g["xrefs"] = g.symbols.map(join), g.xrefs.map(join)
g = add_routes(g, q.GENES, {"G2_encodes_protein": "G2_encodes_protein"}, gene_ids=gene_ids)
g.to_csv(c.RESULTS / "KG_taxon_genes.csv", index=False)

# proteins (UniProt accessions of the model): G1 is the main gene route, G3 (transcript by locus tag) is supporting
p = pd.DataFrame({"id": uniprots})
p = add_routes(p, q.GENES, {"G1_uniprot_protein": "G1_uniprot_protein"}, uniprots=uniprots)
p.to_csv(c.RESULTS / "KG_taxon_proteins.csv", index=False)
lt = pd.DataFrame({"id": locus})
lt = add_routes(lt, q.GENES, {"G3_transcript_locus": "G3_transcript_locus"}, locus=locus)
lt.to_csv(c.RESULTS / "KG_taxon_locus_transcripts.csv", index=False)

# ---- metabolites ---------------------------------------------------------------------------------------------------
m = pd.DataFrame(run(q.METABOLITES, "M_NODES"))
m["synonyms"] = m.synonyms.map(join)
m["in_any_reaction"] = ~m.id.isin({x["id"] for x in run(q.METABOLITES, "M_no_reaction")})
m = add_routes(m, q.METABOLITES, {"M1a_via_R1": "M1a_via_R1", "M1b_via_R2": "M1b_via_R2", "M2_pathway": "M2_pathway"})
m["strong_taxon"] = m.M1a_via_R1 | m.M1b_via_R2
m["strong_taxon_species"] = m.M1a_via_R1_species | m.M1b_via_R2_species
m.to_csv(c.RESULTS / "KG_taxon_metabolites.csv", index=False)

for name, d in (("reactions", r), ("genes", g), ("proteins", p), ("locus transcripts", lt), ("metabolites", m)):
    cols = [x for x in d.columns if d[x].dtype == bool]
    print(f"\nKG {name}: {len(d)} rows\n", d[cols].sum().to_string())
