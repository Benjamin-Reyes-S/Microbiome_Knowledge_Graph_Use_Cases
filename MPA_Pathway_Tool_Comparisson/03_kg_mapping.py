"""Map the Supp. 9 metaproteomics data onto the KEGG modules *stored in the knowledge graph*.
Matching happens inside Neo4j (read-only Cypher): module -> reactions -> CATALYZED_BY KO / EC, intersected with each protein's KO/EC
and the paper's taxonomic requirement. Spectral counts are then summed per sample."""
import json
import numpy as np, pandas as pd
import common as c

df, samples = c.load_proteins()
prots = c.protein_params(df)
rows, mod_rows, matches_out = [], [], {}
for pw, d in c.PATHWAYS.items():
    res = c.cypher(c.MATCH_KG_Q, module=d["module"], prots=prots, taxa=d["taxa"])
    matches = {r["rid"]: r["idx"] for r in res}
    kg_rids = [r["rid"] for r in c.cypher("MATCH (:Module {id:$m})-[:HAS_REACTION]->(r:Reaction) RETURN r.id AS rid", m=d["module"])]
    rx, multi = c.counts_from_matches(matches, df, samples, kg_rids)
    paper_ids = set(c.load_supp_reactions(d["supp_rx"]))
    shared = {k: v for k, v in rx.items() if k in paper_ids}
    for rid, v in rx.items():
        rows.append(dict(pathway=pw, rid=rid, n_proteins=len(matches.get(rid, [])), **dict(zip(samples, v))))
        matches_out[rid + "|" + pw] = [df.iloc[i, 0] for i in matches.get(rid, [])]
    for variant, vec in (("KG all reactions (summed over reactions)", multi),
                         ("KG reactions shared with paper (summed over reactions)", sum(shared.values()))):
        mod_rows.append(dict(pathway=pw, variant=variant, steps_found=sum(1 for v in matches.values() if v), steps_total=len(kg_rids),
                             **dict(zip(samples, vec))))
    print(pw, d["module"], "reactions", len(kg_rids), "found", sum(1 for v in matches.values() if v))
pd.DataFrame(rows).to_csv(c.RESULTS / "kg_reaction_counts.csv", index=False)
pd.DataFrame(mod_rows).to_csv(c.RESULTS / "kg_module_counts.csv", index=False)
json.dump(matches_out, open(c.RESULTS / "kg_matched_proteins.json", "w"))
