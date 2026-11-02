"""Taxonomy of the 4,102 unique model metabolites. A compound has no taxon of its own in any of these databases, so every route goes
through enzymes/reactions/models that use it:
  MX1 KEGG C -> ENZYME list (KEGG) -> UniProt          MX2 ChEBI -> Rhea reactions (+EC) -> UniProt
  MX3 BiGG metabolite -> organisms of BiGG models      MX4 HMDB id -> implied human (human metabolome database, no taxon field, no query)
No taxonomy route in scope for PubChem, InChI, LIPID MAPS, MetaNetX or SEED ids.
Hub compounds (ATP, water...) list thousands of ECs: ECs are checked in order and the search stops at the first human hit,
at most 60 ECs per compound (so a hub with 60 non-human ECs is called 'other' - noted as a caveat).
Output: results/api_metabolites.csv
"""
import pandas as pd
import common as c
import api_queries as q

m = pd.read_csv(c.RESULTS / "model_xrefs_metabolites.csv").fillna("")
chunks = lambda xs, n: [tuple(xs[i:i + n]) for i in range(0, len(xs), n)]
CAP = 60

kc = sorted({x for s in m["kegg.compound"] for x in c.split(s)})
kegg = {}
for res in c.run_parallel(q.kegg_get, chunks(kc, 10), workers=3, label="kegg C").values():
    kegg.update(res)
chebi = sorted({x for s in m.chebi for x in c.split(s)})
rh = c.run_parallel(q.rhea_ecs_of_chebi, chebi, label="rhea chebi")
bg = c.run_parallel(lambda b: q.bigg("metabolites", b), sorted({x for s in m["bigg.metabolite"] for x in c.split(s)}), label="bigg")

# every EC that could be needed (first CAP per compound), checked in parallel once; the per-compound verdict then reads the cache
need = set()
for _, x in m.iterrows():
    for key, src in (("kegg.compound", lambda i: kegg.get(i, {}).get("ecs", [])), ("chebi", lambda i: rh.get(i) or [])):
        need.update(sorted({e for i in c.split(x[key]) for e in src(i)})[:CAP])
print("distinct ECs to check:", len(need))
ec_cache = c.run_parallel(q.uniprot_ec, sorted(need), label="uniprot EC")


def ec_verdict(ecs):
    """first human EC wins; otherwise 'other' if any EC has any organism; 'none' if no EC at all / unknown to UniProt"""
    best = "none"
    for e in sorted(set(ecs))[:CAP]:
        if e not in ec_cache:
            ec_cache[e] = q.uniprot_ec(e)
        if ec_cache[e] == "human":
            return "human"
        if ec_cache[e] == "other":
            best = "other"
    return best


todo = [(i, x) for i, x in m.iterrows()]


def one(item):
    _, x = item
    kg, ch, bi, hm = (c.split(x[k]) for k in ("kegg.compound", "chebi", "bigg.metabolite", "hmdb"))
    k_ecs = [e for i in kg for e in kegg.get(i, {}).get("ecs", [])]
    c_ecs = [e for i in ch for e in (rh.get(i) or [])]
    return {
        "MX1_kegg_enzyme_uniprot": "noid" if not kg else ec_verdict(k_ecs),
        "MX2_chebi_rhea_uniprot": "noid" if not ch else ec_verdict(c_ecs),
        "MX3_bigg_models": "noid" if not bi else ("none" if not any(bg.get(i) for i in bi) else
                           ("human" if any(bg.get(i) and "Homo sapiens" in bg[i]["orgs"] for i in bi) else "other")),
        "MX4_hmdb_implied_human": "noid" if not hm else "human",
    }


rows = [one(t) for t in todo]                      # EC checks are cached; 
out = pd.concat([m[["unique_id", "name"]],
                 pd.DataFrame({"has_kegg": m["kegg.compound"] != "", "has_chebi": m.chebi != "", "has_bigg": m["bigg.metabolite"] != "",
                               "has_hmdb": m.hmdb != "", "has_pubchem": m["pubchem.compound"] != "", "has_inchi": m.inchi != "",
                               "has_lipidmaps": m.lipidmaps != "", "has_metanetx": m["metanetx.chemical"] != "",
                               "has_seed": m["seed.compound"] != ""}), pd.DataFrame(rows)], axis=1)
out.to_csv(c.RESULTS / "api_metabolites.csv", index=False)
c.log_calls("04_parse_metabolites.py", "KEGG, Rhea, UniProt, BiGG")
print(out.filter(like="MX").apply(pd.Series.value_counts).fillna(0).astype(int).to_string())
