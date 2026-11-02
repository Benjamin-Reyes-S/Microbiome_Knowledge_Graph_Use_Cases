"""Which proteins explain the only count difference (reaction R04347, methyl-CoM:CoM methyltransferase)?"""
import json
import common as c
df, samples = c.load_proteins()
kg = set(json.load(open(c.RESULTS / "kg_matched_proteins.json"))["KEGG:R04347|Hydrogenotrophic methanogenesis"])
paper_ids = {"EC:2.1.1.86"} | {"KEGG:" + k for k in "K00577 K00578 K00579 K00580 K00581 K00582 K00583 K00584".split()}
rows = []
for i, (pid, ko, ec, sk) in enumerate(zip(df.iloc[:, 0], df.iloc[:, 1], df.iloc[:, 2], df["superkingdom"])):
    ids = {"KEGG:" + k for k in c._split(ko)} | {"EC:" + e for e in c._split(ec)}
    if sk == "Archaea" and ids & paper_ids and pid not in kg:
        rows.append(dict(protein=pid[:70], ko=ko, ec=ec, total_counts=df.iloc[i][samples].sum()))
import pandas as pd
out = pd.DataFrame(rows); out.to_csv(c.RESULTS / "R04347_proteins_lost_in_kg.csv", index=False)
print(out.to_string(), "\nsum of counts:", out.total_counts.sum() if len(out) else 0)
