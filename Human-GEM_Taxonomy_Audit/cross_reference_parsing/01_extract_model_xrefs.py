"""Read Human-GEM.xml (read-only, in the parent folder) and write every cross-reference of reactions, genes and metabolites.
One row per entity, one column per identifier namespace (several ids joined by ';'). Metabolites are collapsed over compartments
(unique_id, union of the ids of the copies); reactions carry the reaction_class computed by ../02_extract_human_gem.py.
Outputs: results/model_xrefs_{reactions,genes,metabolites}.csv, results/model_xref_namespaces.csv (how many entities carry each id type).
`sbo` (ontology term) and `vmh*` (model's own ids) are not external identifiers and are dropped.
"""
import pandas as pd
import common as c

DROP = {"sbo", "vmhreaction", "vmhmetabolite"}


def ann(o):
    return {k: ([v] if isinstance(v, str) else list(v)) for k, v in o.annotation.items() if k not in DROP}


import cobra
model = cobra.io.read_sbml_model(str(c.MODEL_XML))

rows = []
for r in model.reactions:
    rows.append({"id": r.id, "name": r.name, "subsystem": r.subsystem, **{k: ";".join(v) for k, v in ann(r).items()}})
rx = pd.DataFrame(rows).fillna("")
cls = pd.read_csv(c.PARENT_RESULTS / "Human_GEM_reactions.csv", usecols=["id", "reaction_class"])
rx = rx.merge(cls, on="id", how="left")

gn = pd.DataFrame([{"id": g.id, "symbol": g.name, **{k: ";".join(v) for k, v in ann(g).items()}} for g in model.genes]).fillna("")
gn["ensembl"] = gn.ensembl.map(lambda s: ";".join(x for x in s.split(";") if x.startswith("ENSG")))      # gene-level ids only (not ENST/ENSP)

uni = {}
for m in model.metabolites:
    uid = m.id[:-1] if m.compartment and m.id.endswith(m.compartment) else m.id
    d = uni.setdefault(uid, {"unique_id": uid, "name": m.name})
    for k, v in ann(m).items():
        d.setdefault(k, set()).update(v)
mt = pd.DataFrame([{k: (";".join(sorted(v)) if isinstance(v, set) else v) for k, v in d.items()} for d in uni.values()]).fillna("")

for name, d in (("reactions", rx), ("genes", gn), ("metabolites", mt)):
    d.to_csv(c.RESULTS / f"model_xrefs_{name}.csv", index=False)
    print(name, len(d))
ns = []
for name, d, skip in (("reactions", rx, {"id", "name", "subsystem", "reaction_class"}), ("genes", gn, {"id", "symbol"}),
                      ("metabolites", mt, {"unique_id", "name"})):
    for col in d.columns:
        if col not in skip:
            ns.append(dict(component=name, namespace=col, n=int((d[col] != "").sum()), total=len(d), pct=round(100 * (d[col] != "").mean(), 1)))
pd.DataFrame(ns).to_csv(c.RESULTS / "model_xref_namespaces.csv", index=False)
print(pd.DataFrame(ns).to_string(index=False))
