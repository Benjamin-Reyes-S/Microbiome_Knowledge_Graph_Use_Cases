"""Match MitoCore genes, reactions and metabolites to the human components of the KG (the six CSVs of 01/02), then show what
does NOT match and which reaction types / subsystems those are.

Match levels (same wording for all three entity types):
  matched_human      the entity maps to a KG node that has a STRONG human route
  mapped_not_human   the entity maps to a KG node, but the KG has no strong human route for it
  unmappable         no identifier of the entity exists in the KG
Strong routes: reactions R1/R2 (KEGG id) or human EC (exact id), genes G1 (UniProt Protein is human), metabolites M1a/M1b/M3.
Identifiers used for the match (all exact; an id counts when ANY KG node carrying it has a strong human route):
  reactions    KEGG R id, EC, BiGG reaction id, MetaNetX reaction id   (BiGG / MetaNetX are ARRAYS in the KG: one model id can hit several
                nodes, KEGG:R.. and/or BIGG.R:..; see KG_xref_reactions.csv; BIGG.R nodes have no human route of their own)
  metabolites  KEGG C id, ChEBI id, BiGG metabolite id, MetaNetX chemical id   (union over the compartment copies of a unique compound)
  `status_kegg_ec` keeps the old KEGG+EC-only status of the reactions, to show what the new identifiers add.
Rescue attempts (weaker evidence) are in 04.
"""
import pandas as pd
import common as c
import kg_queries as q
import plotstyle as ps

R = c.RESULTS
mr, mg, mm = (pd.read_csv(R / f"Model_{n}.csv").fillna("") for n in ("reactions", "genes", "metabolites"))
kr, kg, km = (pd.read_csv(R / f"KG_human_{n}.csv", low_memory=False).fillna("") for n in ("reactions", "genes", "metabolites"))
xr, xm = (pd.read_csv(R / f"KG_xref_{n}.csv").fillna("") for n in ("reactions", "metabolites"))
xmap = lambda x, ns: x[x.namespace == ns].groupby("xref").id.apply(list).to_dict()      # xref -> all KG nodes carrying it (1:n)
hit = lambda ids, d: sorted({n for i in ids for n in d.get(i, [])})
split = lambda s: [x for x in s.split(";") if x]
STATUS = ["matched_human", "mapped_not_human", "unmappable"]

# ---------------- genes ----------------
kg_gene_strong = dict(zip(kg.id, kg.strong_human.astype(bool)))
accs = sorted({u.split(":")[1] for s in mg.uniprot for u in split(s)})
prot = {}
for i in range(0, len(accs), 500):
    for row in c.cypher(q.GENES["P_by_uniprot"], ids=accs[i:i + 500]):
        prot[row["uniprot"]] = c.HUMAN in row["taxa"]
mg["in_kg_gene"] = mg.kg_gene.map(lambda s: any(x in kg_gene_strong for x in split(s)))
mg["G1_human_gene"] = mg.kg_gene.map(lambda s: any(kg_gene_strong.get(x, False) for x in split(s)))
mg["uniprot_in_kg"] = mg.uniprot.map(lambda s: any(u.split(":")[1] in prot for u in split(s)))
mg["uniprot_human"] = mg.uniprot.map(lambda s: any(prot.get(u.split(":")[1], False) for u in split(s)))
mg["status"] = "unmappable"
mg.loc[mg.in_kg_gene | mg.uniprot_in_kg, "status"] = "mapped_not_human"
mg.loc[mg.G1_human_gene | mg.uniprot_human, "status"] = "matched_human"
mg.to_csv(R / "match_genes.csv", index=False)
gene_ok = set(mg.id[mg.status == "matched_human"])

# ---------------- metabolites (one row per unique compound; compartment copies collapsed, ids = union over the copies) ----------------
kc = km.set_index("id")
u = mm.drop_duplicates("unique_id").copy()
NS = {"kegg": ("kegg_compound", None), "chebi": ("chebi", "chebi"), "bigg": ("bigg_metabolite", "bigg.metabolite"),
      "mnx": ("metanetx_chemical", "metanetx.chemical")}
for tag, (col, ns) in NS.items():
    u[col] = u.unique_id.map(mm.groupby("unique_id")[col].agg(lambda s: ";".join(sorted({x for v in s for x in split(v)}))))
    d = None if ns is None else xmap(xm, ns)
    nodes = u[col].map(lambda s: [x for x in split(s) if x in kc.index] if d is None else hit(split(s), d))
    u[f"{tag}_nodes"] = nodes.map(";".join)
    u[f"{tag}_in_kg"] = nodes.map(bool)
    u[f"{tag}_human"] = nodes.map(lambda l: any(bool(kc.strong_human.get(x, False)) for x in l))
u["in_kg"] = u[[f"{t}_in_kg" for t in NS]].any(axis=1)
u["M1_human"] = u[[f"{t}_human" for t in NS]].any(axis=1)
u["M2_pathway"] = u[[f"{t}_nodes" for t in NS]].apply(lambda r: any(bool(kc.M2_pathway.get(x, False)) for v in r for x in split(v)), axis=1)
u["has_kegg_id"] = u.kegg_compound != ""
u["has_any_kg_id"] = u[[v[0] for v in NS.values()]].ne("").any(axis=1)
u["status_kegg"] = "unmappable"
u.loc[u.kegg_in_kg, "status_kegg"] = "mapped_not_human"
u.loc[u.kegg_human, "status_kegg"] = "matched_human"
u["status"] = "unmappable"
u.loc[u.in_kg, "status"] = "mapped_not_human"
u.loc[u.M1_human, "status"] = "matched_human"
u.to_csv(R / "match_metabolites.csv", index=False)

# ---------------- reactions ----------------
kr_idx = kr.set_index("id")
enz = pd.DataFrame(c.cypher(q.REACTIONS["E_NODES"]))
ec_in_kg = set(enz.id)
ec_human = {x["id"] for x in c.cypher(q.REACTIONS["E_human"], h=c.HUMAN)}
mr["kegg_in_kg"] = mr.kegg_reaction.map(lambda s: any(x in kr_idx.index for x in split(s)))
mr["kegg_human"] = mr.kegg_reaction.map(lambda s: any(bool(kr_idx.strong_human.get(x, False)) for x in split(s)))
mr["ec_in_kg"] = mr.ec.map(lambda s: any(x in ec_in_kg for x in split(s)))
mr["ec_human"] = mr.ec.map(lambda s: any(x in ec_human for x in split(s)))
mr["gene_human"] = mr.genes.map(lambda s: any(x in gene_ok for x in split(s)))      # supporting only, not an identity match
strong_r = dict(zip(kr.id, kr.strong_human.astype(bool)))
for tag, col, ns in (("bigg", "bigg", "bigg.reaction"), ("mnx", "metanetx", "metanetx.reaction")):
    d = xmap(xr, ns)
    nodes = mr[col].map(lambda s: hit(split(s), d))
    mr[f"{tag}_nodes"] = nodes.map(";".join)                                   # every KG node the model ids point to (arrays -> 1:n)
    mr[f"{tag}_in_kg"] = nodes.map(bool)
    mr[f"{tag}_human"] = nodes.map(lambda l: any(strong_r.get(n, False) for n in l))
mr["status_kegg_ec"] = "unmappable"
mr.loc[mr.kegg_in_kg | mr.ec_in_kg, "status_kegg_ec"] = "mapped_not_human"
mr.loc[mr.kegg_human | mr.ec_human, "status_kegg_ec"] = "matched_human"
mr["status"] = "unmappable"
mr.loc[mr.kegg_in_kg | mr.ec_in_kg | mr.bigg_in_kg | mr.mnx_in_kg, "status"] = "mapped_not_human"
mr.loc[mr.kegg_human | mr.ec_human | mr.bigg_human | mr.mnx_human, "status"] = "matched_human"
ROUTES = (("KEGG", "kegg_human"), ("EC", "ec_human"), ("BiGG", "bigg_human"), ("MetaNetX", "mnx_human"))
mr["route"] = mr.apply(lambda r: "+".join(n for n, f in ROUTES if r[f]), axis=1)
mr.to_csv(R / "match_reactions.csv", index=False)

# ---------------- summary ----------------
summ = []
for ent, d in (("genes", mg), ("reactions", mr), ("metabolites (unique)", u)):
    vc = d.status.value_counts()
    summ.append(dict(entity=ent, total=len(d), **{s: int(vc.get(s, 0)) for s in STATUS},
                     pct_matched_human=round(100 * vc.get("matched_human", 0) / len(d), 1)))
summ = pd.DataFrame(summ)
summ.to_csv(R / "match_summary.csv", index=False)
print(summ.to_string(index=False))
gain = pd.DataFrame([
    dict(entity="reactions", ids="KEGG+EC only", matched_human=int((mr.status_kegg_ec == "matched_human").sum()), mapped=int((mr.status_kegg_ec != "unmappable").sum())),
    dict(entity="reactions", ids="+ BiGG + MetaNetX", matched_human=int((mr.status == "matched_human").sum()), mapped=int((mr.status != "unmappable").sum())),
    dict(entity="metabolites", ids="KEGG only", matched_human=int((u.status_kegg == "matched_human").sum()), mapped=int((u.status_kegg != "unmappable").sum())),
    dict(entity="metabolites", ids="+ ChEBI + BiGG + MetaNetX", matched_human=int((u.status == "matched_human").sum()), mapped=int((u.status != "unmappable").sum()))])
gain.to_csv(R / "match_id_gain.csv", index=False)
print("\nwhat the new identifiers add:\n", gain.to_string(index=False))
print("reactions: model ids hitting >1 KG node  BiGG", int(mr.bigg_nodes.str.contains(";").sum()), "| MetaNetX", int(mr.mnx_nodes.str.contains(";").sum()),
      "| mapped only to a BIGG.R node (no human route possible):", int((mr.bigg_in_kg & ~mr.kegg_in_kg & ~mr.ec_in_kg & ~mr.bigg_human & ~mr.mnx_human).sum()))
print("metabolites human per identifier:", {t: int(u[f"{t}_human"].sum()) for t in NS}, "| in KG:", {t: int(u[f"{t}_in_kg"].sum()) for t in NS})
print("\ngenes: in KG as gene node", int(mg.in_kg_gene.sum()), "| human via KG gene", int(mg.G1_human_gene.sum()),
      "| UniProt is human Protein", int(mg.uniprot_human.sum()))
print("metabolites: with KEGG id", int(u.has_kegg_id.sum()), "of", len(u), "| with any KG-usable id", int(u.has_any_kg_id.sum()), "| in KG (any id)", int(u.in_kg.sum()))
print("reactions by route:", mr[mr.status == "matched_human"].route.value_counts().to_dict(), "| gene-supported among unmatched:",
      int(mr[mr.status != "matched_human"].gene_human.sum()))

# ---------------- what does not match: reaction type, subsystem ----------------
cls = pd.crosstab(mr.reaction_class, mr.status).reindex(columns=STATUS, fill_value=0)
cls["total"] = cls.sum(axis=1)
cls.to_csv(R / "match_by_reaction_class.csv")
print("\n", cls.to_string())

sub = pd.crosstab(mr.subsystem, mr.status).reindex(columns=STATUS, fill_value=0)
sub["total"] = sub.sum(axis=1)
sub["unmatched"] = sub.mapped_not_human + sub.unmappable
sub["pct_unmatched"] = (100 * sub.unmatched / sub.total).round(1)
sub = sub.sort_values("unmatched", ascending=False)
sub.to_csv(R / "unmatched_by_subsystem.csv")


print("\nsubsystems with most unmatched reactions:\n", sub.head(12).to_string())

# ---------------- plots ----------------
plt = ps.plt
# fig1: overall match per entity type
fig, ax = plt.subplots(figsize=(8, 3.2))
ps.stacked_h(ax, summ.entity.tolist(), summ[STATUS].values, STATUS, totals=summ.total.tolist())
fig.suptitle("MitoCore entities matched to human components of the KG", x=0.01, ha="left", fontsize=12)
ps.finish(fig, "fig1_match_overview.png")

# fig2: reaction class x match status
cl = cls.reindex(["metabolic", "transport", "exchange_demand", "isolated", "artificial"]).dropna()
fig, ax = plt.subplots(figsize=(8, 3.6))
ps.stacked_h(ax, cl.index.tolist(), cl[STATUS].values, STATUS, totals=cl.total.tolist())
fig.suptitle("Reaction type vs. match to human KG reactions / enzymes", x=0.01, ha="left", fontsize=12)
ps.finish(fig, "fig2_reaction_class_match.png")

# fig3: subsystems with most unmatched reactions
top = sub.head(20)[::-1]
fig, ax = plt.subplots(figsize=(8.5, 7))
ps.stacked_h(ax, top.index.tolist(), top[STATUS].values, STATUS, totals=top.total.tolist(), labels=False)
fig.suptitle("Subsystems with the most reactions NOT matched to human KG evidence (top 20)", x=0.01, ha="left", fontsize=12)
ps.finish(fig, "fig3_unmatched_subsystems.png")
