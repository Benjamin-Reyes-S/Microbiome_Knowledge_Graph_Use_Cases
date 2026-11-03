"""Extract MitoCore into three tables: Model_{reactions,genes,metabolites}.csv (one row per entity, ids and xrefs only).
Identifiers are written in KG format (KEGG:R#####, EC:x.x.x.x, KEGG:C#####, UniProtKB:<acc>, KEGG:hsa:<entrez>) so 03 can join directly."""
import re
import pandas as pd
import common as c

m = c.load_model()


def ann(o, key):
    return [str(x) for x in c.xref_ns(o).get(key, [])]


def tail(x):  # identifiers.org URI or prefixed id -> bare id
    return re.split(r"[/:]", x.strip())[-1] if x.startswith("http") else x.strip()


j = lambda xs: ";".join(sorted(set(xs)))

def classify(r):
    """Reaction type from id, structure and MitoCore's subsystem label (it has no 'Isolated' class).
    artificial: biomass / objective-function reactions and the 'diffusion / artificial' transporters;
    exchange_demand: boundary reaction (EX_/DM_/sink, subsystem 'Boundary ...', or a single metabolite species);
    transport: subsystem mentions transport, or same compounds on both sides in different compartments; else metabolic."""
    left = {c.uid(x.id) for x, co in r.metabolites.items() if co < 0}
    right = {c.uid(x.id) for x, co in r.metabolites.items() if co > 0}
    ncomp = len({x.compartment for x in r.metabolites})
    sub = (r.subsystem or "").lower()
    if r.id.lower().startswith(("biomass", "of_")) or sub.startswith("objective function") or "artificial" in sub: return "artificial"
    if len(r.metabolites) == 1 or r.id.startswith(("EX_", "DM_", "sink_")) or sub.startswith("boundary"): return "exchange_demand"
    if "transport" in sub or (ncomp > 1 and left == right): return "transport"
    return "metabolic"


def signature(r):
    """Compartment-free stoichiometry, direction-independent: identical for the same reaction copied into several compartments."""
    s = sorted(f"{c.uid(x.id)}:{co:g}" for x, co in r.metabolites.items())
    t = sorted(f"{c.uid(x.id)}:{-co:g}" for x, co in r.metabolites.items())
    return "|".join(min(s, t))


rows = []
for r in m.reactions:
    rows.append(dict(
        id=r.id, name=r.name, subsystem=r.subsystem or "(none)", reaction_class=classify(r),
        compartments=j(x.compartment for x in r.metabolites), signature=signature(r),
        n_genes=len(r.genes), gpr=r.gene_reaction_rule,
        genes=j(g.id for g in r.genes), lower=r.lower_bound, upper=r.upper_bound,
        n_metabolites=len({c.uid(x.id) for x in r.metabolites}), uids=j(c.uid(x.id) for x in r.metabolites),
        compounds=j("KEGG:" + k for x in r.metabolites for k in ann(x, "kegg.compound")),
        kegg_reaction=j("KEGG:" + tail(x) for x in ann(r, "kegg.reaction")),
        ec=j("EC:" + tail(x) for x in ann(r, "ec-code")),
        rhea=j(ann(r, "rhea")), bigg=j(ann(r, "bigg.reaction")), metanetx=j(ann(r, "metanetx.reaction"))))
pd.DataFrame(rows).to_csv(c.RESULTS / "Model_reactions.csv", index=False)

rows = []
for g in m.genes:
    ent = [tail(x) for x in ann(g, "ncbigene")]
    up = [tail(x) for x in ann(g, "uniprot")]
    rows.append(dict(id=g.id, symbol=j(ann(g, "hgnc.symbol")), kg_gene=j("KEGG:hsa:" + e for e in ent),
                     entrez=j(ent), uniprot=j("UniProtKB:" + u for u in up),
                     n_reactions=len(g.reactions)))
pd.DataFrame(rows).to_csv(c.RESULTS / "Model_genes.csv", index=False)

rows = []
for x in m.metabolites:
    rows.append(dict(id=x.id, unique_id=c.uid(x.id), name=x.name, compartment=x.compartment, formula=x.formula,
                     kegg_compound=j("KEGG:" + tail(k) for k in ann(x, "kegg.compound")),
                     chebi=j(ann(x, "chebi")), hmdb=j(ann(x, "hmdb")), n_reactions=len(x.reactions)))
pd.DataFrame(rows).to_csv(c.RESULTS / "Model_metabolites.csv", index=False)

for n in ("reactions", "genes", "metabolites"):
    d = pd.read_csv(c.RESULTS / f"Model_{n}.csv")
    print(f"{n}: {len(d)} rows; non-empty per column:", {k: int(d[k].notna().sum()) for k in d.columns if k not in ("id", "name")})

rc = pd.read_csv(c.RESULTS / "Model_reactions.csv")
print(rc.reaction_class.value_counts().to_string())
print("unique compartment-free reaction signatures:", rc.signature.nunique(), "of", len(rc),
      "| metabolic only:", rc[rc.reaction_class == "metabolic"].signature.nunique(), "of", int((rc.reaction_class == "metabolic").sum()))
