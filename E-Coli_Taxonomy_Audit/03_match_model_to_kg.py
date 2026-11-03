"""Match iML1515 genes, reactions and metabolites to the components of the KG that belong to the target taxon (E. coli K-12 MG1655:
strain + K-12 node, see common.py), then show what does NOT match and which reaction types those are.

Match levels (same wording for all three entity types):
  matched_taxon      the entity maps to a KG node that has a STRONG route to the taxon (strict scope: strain + K-12)
  mapped_not_taxon   the entity maps to a KG node, but the KG has no strong taxon route for it
  unmappable         no identifier of the entity exists in the KG
Strong routes: reactions R1/R2 (KEGG id) or EC of a taxon protein (exact id), genes G1 (UniProt Protein in the taxon) or G2 (Gene
ENCODES such a Protein), metabolites M1a/M1b. Column `species_only` marks entities not matched in the strict scope but matched when
the species node (E. coli, 562) is added (wider scope); they are NOT counted as matched.
Reaction match here is strict (exact KEGG R id or exact EC); rescue attempts are in 04.
"""
import pandas as pd
import common as c
import kg_queries as q
import plotstyle as ps

R = c.RESULTS
mr, mg, mm = (pd.read_csv(R / f"Model_{n}.csv").fillna("") for n in ("reactions", "genes", "metabolites"))
kr, kg, km = (pd.read_csv(R / f"KG_taxon_{n}.csv").fillna("") for n in ("reactions", "genes", "metabolites"))
kp = pd.read_csv(R / "KG_taxon_proteins.csv")
split = lambda s: [x for x in s.split(";") if x]
STATUS = ["matched_taxon", "mapped_not_taxon", "unmappable"]

# ---------------- genes ----------------
g2 = dict(zip(kg.id, kg.G2_encodes_protein.astype(bool)))                 # KG Gene node (NCBIGene:<entrez>) -> ENCODES a taxon protein
g2s = dict(zip(kg.id, kg.G2_encodes_protein_species.astype(bool)))
g1 = dict(zip(kp.id, kp.G1_uniprot_protein.astype(bool)))                 # model UniProt accession -> Protein in the taxon
g1s = dict(zip(kp.id, kp.G1_uniprot_protein_species.astype(bool)))
accs = sorted({u.split(":")[1] for s in mg.uniprot for u in split(s)})
prot_in_kg = set()
for i in range(0, len(accs), 500):
    prot_in_kg |= {"UniProtKB:" + row["uniprot"] for row in c.cypher(q.GENES["P_by_uniprot"], ids=accs[i:i + 500])}
anyin = lambda s, d: any(d.get(x, False) for x in split(s))
mg["in_kg_gene"] = mg.kg_gene.map(lambda s: any(x in g2 for x in split(s)))
mg["uniprot_in_kg"] = mg.uniprot.map(lambda s: any(x in prot_in_kg for x in split(s)))
mg["G2_taxon_gene"] = mg.kg_gene.map(lambda s: anyin(s, g2))
mg["G1_taxon_protein"] = mg.uniprot.map(lambda s: anyin(s, g1))
mg["species_hit"] = mg.kg_gene.map(lambda s: anyin(s, g2s)) | mg.uniprot.map(lambda s: anyin(s, g1s))
mg["status"] = "unmappable"
mg.loc[mg.in_kg_gene | mg.uniprot_in_kg, "status"] = "mapped_not_taxon"
mg.loc[mg.G2_taxon_gene | mg.G1_taxon_protein, "status"] = "matched_taxon"
mg["species_only"] = (mg.status != "matched_taxon") & mg.species_hit
mg.to_csv(R / "match_genes.csv", index=False)
gene_ok = set(mg.id[mg.status == "matched_taxon"])

# ---------------- metabolites (one row per unique compound; compartment copies collapsed) ----------------
kc = km.set_index("id")
u = mm.drop_duplicates("unique_id").copy()
u["in_kg"] = u.kegg_compound.map(lambda s: any(x in kc.index for x in split(s)))
u["M1_taxon"] = u.kegg_compound.map(lambda s: any(bool(kc.strong_taxon.get(x, False)) for x in split(s)))
u["species_hit"] = u.kegg_compound.map(lambda s: any(bool(kc.strong_taxon_species.get(x, False)) for x in split(s)))
u["M2_pathway"] = u.kegg_compound.map(lambda s: any(bool(kc.M2_pathway.get(x, False)) for x in split(s)))
u["has_kegg_id"] = u.kegg_compound != ""
u["status"] = "unmappable"
u.loc[u.in_kg, "status"] = "mapped_not_taxon"
u.loc[u.M1_taxon, "status"] = "matched_taxon"
u["species_only"] = (u.status != "matched_taxon") & u.species_hit
u.to_csv(R / "match_metabolites.csv", index=False)

# ---------------- reactions ----------------
kr_idx = kr.set_index("id")
enz = pd.DataFrame(c.cypher(q.REACTIONS["E_NODES"]))
ec_in_kg = set(enz.id)
ec_rows = c.cypher(q.REACTIONS["E_taxon"], taxa=c.TAXA, taxa_sp=c.TAXA_SPECIES)
ec_taxon = {x["id"] for x in ec_rows if x["n_support"] > 0}
ec_species = {x["id"] for x in ec_rows}
mr["kegg_in_kg"] = mr.kegg_reaction.map(lambda s: any(x in kr_idx.index for x in split(s)))
mr["kegg_taxon"] = mr.kegg_reaction.map(lambda s: any(bool(kr_idx.strong_taxon.get(x, False)) for x in split(s)))
mr["kegg_species"] = mr.kegg_reaction.map(lambda s: any(bool(kr_idx.strong_taxon_species.get(x, False)) for x in split(s)))
mr["ec_in_kg"] = mr.ec.map(lambda s: any(x in ec_in_kg for x in split(s)))
mr["ec_taxon"] = mr.ec.map(lambda s: any(x in ec_taxon for x in split(s)))
mr["ec_species"] = mr.ec.map(lambda s: any(x in ec_species for x in split(s)))
mr["gene_taxon"] = mr.genes.map(lambda s: any(x in gene_ok for x in split(s)))      # supporting only, not an identity match
mr["status"] = "unmappable"
mr.loc[mr.kegg_in_kg | mr.ec_in_kg, "status"] = "mapped_not_taxon"
mr.loc[mr.kegg_taxon | mr.ec_taxon, "status"] = "matched_taxon"
mr["species_only"] = (mr.status != "matched_taxon") & (mr.kegg_species | mr.ec_species)
mr["route"] = ""
mr.loc[mr.kegg_taxon & ~mr.ec_taxon, "route"] = "KEGG"
mr.loc[~mr.kegg_taxon & mr.ec_taxon, "route"] = "EC"
mr.loc[mr.kegg_taxon & mr.ec_taxon, "route"] = "KEGG+EC"
mr.to_csv(R / "match_reactions.csv", index=False)

# ---------------- summary ----------------
summ = []
for ent, d in (("genes", mg), ("reactions", mr), ("metabolites (unique)", u)):
    vc = d.status.value_counts()
    summ.append(dict(entity=ent, total=len(d), **{s: int(vc.get(s, 0)) for s in STATUS},
                     species_only=int(d.species_only.sum()),
                     pct_matched_taxon=round(100 * vc.get("matched_taxon", 0) / len(d), 1)))
summ = pd.DataFrame(summ)
summ.to_csv(R / "match_summary.csv", index=False)
print(summ.to_string(index=False))
print("\ngenes: Gene node in KG", int(mg.in_kg_gene.sum()), "| taxon via Gene-ENCODES", int(mg.G2_taxon_gene.sum()),
      "| UniProt Protein in taxon", int(mg.G1_taxon_protein.sum()))
print("metabolites: with KEGG id", int(u.has_kegg_id.sum()), "of", len(u), "| KEGG id in KG", int(u.in_kg.sum()))
print("reactions by route:", mr[mr.status == "matched_taxon"].route.value_counts().to_dict(), "| gene-supported among unmatched:",
      int(mr[mr.status != "matched_taxon"].gene_taxon.sum()))

# ---------------- what does not match: reaction type ----------------
# (iML1515 has no subsystem labels, so there is no per-subsystem table / figure here)
cls = pd.crosstab(mr.reaction_class, mr.status).reindex(columns=STATUS, fill_value=0)
cls["total"] = cls.sum(axis=1)
cls.to_csv(R / "match_by_reaction_class.csv")
print("\n", cls.to_string())

# ---------------- plots ----------------
plt = ps.plt
LABEL = c.TAXON_LABEL
fig, ax = plt.subplots(figsize=(8, 3.2))
ps.stacked_h(ax, summ.entity.tolist(), summ[STATUS].values, STATUS, totals=summ.total.tolist())
fig.suptitle(f"iML1515 entities matched to {LABEL} components of the KG ({c.VERSION_KG})", x=0.01, ha="left", fontsize=12)
ps.finish(fig, "fig1_match_overview.png")

cl = cls.reindex(["metabolic", "transport", "exchange_demand", "artificial"]).dropna()
fig, ax = plt.subplots(figsize=(8, 3.4))
ps.stacked_h(ax, cl.index.tolist(), cl[STATUS].values, STATUS, totals=cl.total.tolist())
fig.suptitle(f"Reaction type vs. match to {LABEL} KG reactions / enzymes ({c.VERSION_KG})", x=0.01, ha="left", fontsize=12)
ps.finish(fig, "fig2_reaction_class_match.png")
