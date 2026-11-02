"""Taxonomy of the 2,870 model genes through four parsing routes (each needs its own identifier AND its own API):
  GX1 UniProt accession -> UniProt organism            GX2 NCBI Gene id -> NCBI esummary taxid
  GX3 Ensembl ENSG id   -> Ensembl lookup species      GX4 HGNC symbol  -> HGNC (human-only database: existence = human)
Output: results/api_genes.csv (id, one verdict column per route: human / other / none / noid).
"""
import pandas as pd
import common as c
import api_queries as q

g = pd.read_csv(c.RESULTS / "model_xrefs_genes.csv").fillna("")
chunks = lambda xs, n: [xs[i:i + n] for i in range(0, len(xs), n)]

uni_ids = sorted({a for s in g.uniprot for a in c.split(s)})
org = {}
for res in c.run_parallel(q.uniprot_organisms, [tuple(x) for x in chunks(uni_ids, 300)], label="uniprot").values():
    org.update(res)
ent_ids = sorted({a for s in g.ncbigene for a in c.split(s)})
tax = {}
for res in c.run_parallel(q.ncbi_gene_taxids, [tuple(x) for x in chunks(ent_ids, 200)], workers=2, label="ncbi").values():
    tax.update(res)
ens_ids = sorted({a for s in g.ensembl for a in c.split(s)})
sp = {}
for res in c.run_parallel(q.ensembl_species, [tuple(x) for x in chunks(ens_ids, 500)], workers=2, label="ensembl").values():
    sp.update(res)
sym = sorted({a for s in g["hgnc.symbol"] for a in c.split(s)})
hg = c.run_parallel(q.hgnc_symbol_exists, sym, workers=4, label="hgnc")


def route(ids, lookup):
    if not ids:
        return "noid"
    vals = [lookup(i) for i in ids]
    vals = [v for v in vals if v is not None]
    return "none" if not vals else ("human" if any(vals) else "other")


out = pd.DataFrame({"id": g.id, "symbol": g.symbol})
out["GX1_uniprot"] = [route(c.split(s), lambda a: org[a] == c.HUMAN if a in org else None) for s in g.uniprot]
out["GX2_ncbigene"] = [route(c.split(s), lambda a: tax[a] == c.HUMAN if a in tax else None) for s in g.ncbigene]
out["GX3_ensembl"] = [route(c.split(s), lambda a: sp[a] == "homo_sapiens" if a in sp else None) for s in g.ensembl]
out["GX4_hgnc"] = [route(c.split(s), lambda a: True if hg.get(a) else None) for s in g["hgnc.symbol"]]
out.to_csv(c.RESULTS / "api_genes.csv", index=False)
c.log_calls("02_parse_genes.py", "UniProt, NCBI Gene, Ensembl, HGNC")
print(out.iloc[:, 2:].apply(pd.Series.value_counts).fillna(0).astype(int).to_string())
