"""Taxonomy of the 12,748 model reactions through seven parsing routes (identifier -> API chain -> organisms):
  RX1 KEGG R -> EC (KEGG) -> UniProt                 RX2 KEGG R -> KO (KEGG) -> organisms listed under the KO (HSA = human)
  RX3 model EC -> UniProt                            RX4 Rhea id (+ its 3 sibling ids) -> UniProt
  RX5 BiGG reaction -> organisms of BiGG models      RX6 BiGG reaction -> EC (database_links) -> UniProt
  RX7 Reactome id -> species
No API route exists in scope for MetaNetX, SEED or TCDB ids (no taxonomy service / flat-file download only).
Output: results/api_reactions.csv  (id, reaction_class, has-id flags, one verdict column per route: human / other / none / noid)
"""
import pandas as pd
import common as c
import api_queries as q

r = pd.read_csv(c.RESULTS / "model_xrefs_reactions.csv").fillna("")
chunks = lambda xs, n: [tuple(xs[i:i + n]) for i in range(0, len(xs), n)]

# --- KEGG: reaction entries -> ECs and KOs; KO entries -> organism codes
kr = sorted({x for s in r["kegg.reaction"] for x in c.split(s)})
kegg = {}
for res in c.run_parallel(q.kegg_get, chunks(kr, 10), workers=3, label="kegg R").values():
    kegg.update(res)
kos = sorted({k for v in kegg.values() for k in v["kos"]})
ko_org = {}
for res in c.run_parallel(q.kegg_get, chunks(kos, 10), workers=3, label="kegg K").values():
    ko_org.update({k: v["orgs"] for k, v in res.items()})

# --- BiGG reactions (one call each) and Reactome events
bg = c.run_parallel(lambda b: q.bigg("reactions", b), sorted({x for s in r["bigg.reaction"] for x in c.split(s)}), label="bigg")
rc = c.run_parallel(q.reactome_species, sorted({x for s in r.reactome for x in c.split(s)}), label="reactome")

# --- UniProt: every EC seen anywhere (model, KEGG, BiGG) and every Rhea id
ecs = sorted({e for s in r["ec-code"] for e in c.split(s)} | {e for v in kegg.values() for e in v["ecs"]}
             | {e for v in bg.values() if v for e in v["ecs"]})
ec_v = c.run_parallel(q.uniprot_ec, ecs, label="uniprot EC")
rh_v = c.run_parallel(q.uniprot_rhea, sorted({x for s in r.rhea for x in c.split(s)}), label="uniprot Rhea")
print("distinct ECs queried:", len(ecs))


def combine(vs):
    """vs: list of per-id verdicts (None = id not answered) -> route verdict"""
    vs = [v for v in vs if v]
    return "none" if not vs else ("human" if "human" in vs else ("other" if "other" in vs else "none"))


def ec_route(eclist):
    return combine([ec_v.get(e) for e in eclist])


out = r[["id", "reaction_class"]].copy()
for col, key in (("has_kegg", "kegg.reaction"), ("has_ec", "ec-code"), ("has_rhea", "rhea"), ("has_bigg", "bigg.reaction"),
                 ("has_reactome", "reactome"), ("has_metanetx", "metanetx.reaction"), ("has_seed", "seed.reaction"), ("has_tcdb", "tcdb")):
    out[col] = r[key] != ""
rows = []
for _, x in r.iterrows():
    kg_ids, ec_ids, rh_ids, bg_ids, rc_ids = (c.split(x[k]) for k in ("kegg.reaction", "ec-code", "rhea", "bigg.reaction", "reactome"))
    k_ec = sorted({e for i in kg_ids for e in kegg.get(i, {}).get("ecs", [])})
    k_ko = sorted({k for i in kg_ids for k in kegg.get(i, {}).get("kos", [])})
    b_ec = sorted({e for i in bg_ids if bg.get(i) for e in bg[i]["ecs"]})
    rows.append({
        "RX1_kegg_ec_uniprot": "noid" if not kg_ids else ec_route(k_ec),
        "RX2_kegg_ko_genes": "noid" if not kg_ids else combine([("human" if "HSA" in ko_org[k] else "other") if ko_org.get(k) else None for k in k_ko]),
        "RX3_ec_uniprot": "noid" if not ec_ids else ec_route(ec_ids),
        "RX4_rhea_uniprot": "noid" if not rh_ids else combine([rh_v.get(i) for i in rh_ids]),
        "RX5_bigg_models": "noid" if not bg_ids else combine([(("human" if "Homo sapiens" in bg[i]["orgs"] else "other") if bg.get(i) and bg[i]["orgs"] else None) for i in bg_ids]),
        "RX6_bigg_ec_uniprot": "noid" if not bg_ids else ec_route(b_ec),
        "RX7_reactome_species": "noid" if not rc_ids else combine([(("human" if "Homo sapiens" in rc[i] else "other") if rc.get(i) else None) for i in rc_ids]),
    })
out = pd.concat([out, pd.DataFrame(rows)], axis=1)
out.to_csv(c.RESULTS / "api_reactions.csv", index=False)
c.log_calls("03_parse_reactions.py", "KEGG, UniProt, Rhea, BiGG, Reactome")
print(out.filter(like="RX").apply(pd.Series.value_counts).fillna(0).astype(int).to_string())
