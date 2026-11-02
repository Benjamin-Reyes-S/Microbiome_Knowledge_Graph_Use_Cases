"""Read the three KEGG modules from the knowledge graph (read-only Cypher) and diff them with the paper's pathway files."""
import pandas as pd
import common as c

Q = """
MATCH (m:Module {id: $module})-[:HAS_REACTION]->(r:Reaction)
OPTIONAL MATCH (r)-[:CATALYZED_BY]->(k:KeggOrthology)
WITH m, r, collect(DISTINCT k.id) AS kos
OPTIONAL MATCH (r)-[:CATALYZED_BY]->(e:Enzyme)
RETURN m.name AS module_name, r.id AS rid, r.name AS reaction_name, r.equation AS equation, kos, collect(DISTINCT e.id) AS ecs
"""
Q_MOD_KO = "MATCH (m:Module {id:$module})-[:HAS_ORTHOLOG]->(k:KeggOrthology) RETURN collect(k.id) AS kos"

rows, diff = [], []
for pw, d in c.PATHWAYS.items():
    kg = c.cypher(Q, module=d["module"])
    mod_kos = set(c.cypher(Q_MOD_KO, module=d["module"])[0]["kos"])
    paper = {s["rid"]: s for s in c.parse_pathway_csv(c.DATA / "pathways" / d["kegg_csv"])}
    kg_ids = {r["rid"] for r in kg}
    for r in kg:
        rows.append(dict(pathway=pw, module=d["module"], rid=r["rid"], name=r["reaction_name"], equation=r["equation"],
                         kg_ko=",".join(sorted(r["kos"])), kg_ec=",".join(sorted(r["ecs"])), in_paper=r["rid"] in paper))
        p = paper.get(r["rid"])
        pk = {"KEGG:" + k for k in p["ko"]} if p else set()
        pe = {"EC:" + e for e in p["ec"]} if p else set()
        diff.append(dict(pathway=pw, rid=r["rid"], status="both" if p else "KG only",
                         ko_only_kg=",".join(sorted(set(r["kos"]) - pk)), ko_only_paper=",".join(sorted(pk - set(r["kos"]))),
                         ec_only_kg=",".join(sorted(set(r["ecs"]) - pe)), ec_only_paper=",".join(sorted(pe - set(r["ecs"])))))
    for rid, p in paper.items():
        if rid not in kg_ids:
            diff.append(dict(pathway=pw, rid=rid, status="paper only", ko_only_paper=",".join(p["ko"]), ec_only_paper=",".join(p["ec"])))
    print(f"{pw}: KG {len(kg_ids)} reactions, paper {len(paper)}, shared {len(kg_ids & set(paper))}; module-level KOs in KG: {len(mod_kos)}")
pd.DataFrame(rows).to_csv(c.RESULTS / "kg_module_definitions.csv", index=False)
dd = pd.DataFrame(diff).fillna("")
dd.to_csv(c.RESULTS / "definition_diff_kg_vs_paper.csv", index=False)
print(dd.to_string(max_colwidth=60))
