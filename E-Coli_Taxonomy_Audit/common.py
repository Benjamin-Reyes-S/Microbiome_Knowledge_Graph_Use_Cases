"""Shared helpers for the E. coli iML1515 taxonomy audit against the knowledge graph.

All knowledge-graph access goes through `cypher()`, which only allows read-only
queries against the Neo4j HTTP transaction endpoint (user neo4j, empty password).
Taxon scope (TAXA / TAXA_SPECIES) is defined below; queries receive it as $taxa / $taxa_sp.
"""
import base64, json, re, urllib.request
from pathlib import Path
import os, sys

import pandas as pd

BASE = Path(__file__).resolve().parent

# The KG endpoint is mandatory and comes from the environment (set by organizer.sh from its port argument).
NEO4J_URL = os.getenv("NEO4J_URL")
_m = re.match(r"^https?://[^:/]+:(\d+)/", NEO4J_URL or "")
if not _m:
    sys.exit("NEO4J_URL is not set or has no port. Run the pipeline via: ./organizer.sh <port>")
VERSION_KG = f"KG_port{_m.group(1)}"   # one output folder per KG version

OUT = BASE / VERSION_KG
DATA, RESULTS, FIGURES = OUT / "data", OUT / "results", OUT / "figures"
for p in (DATA, RESULTS, FIGURES):
    p.mkdir(parents=True, exist_ok=True)

MODEL_XML = BASE / "iML1515.xml"
# Target taxon of the model: E. coli str. K-12 substr. MG1655 (ncbi_taxid 511145). In the KG almost all of its proteins are
# attached to the parent K-12 node (83333, the UniProt reference-proteome taxon) or to the species (562), not to the strain node:
# Protein-BELONGS_TO counts 511145: 50 | 83333: 4,579 | 562: 18,735. Matching on the strain node alone would miss nearly everything.
TAXON = "NCBITaxon:511145"                              # the strain
TAXA = [TAXON, "NCBITaxon:83333"]                       # STRICT scope: strain + K-12 (the strain's own proteome)
TAXA_SPECIES = TAXA + ["NCBITaxon:562"]                 # WIDER scope: + species E. coli (reported separately as *_species)
TAXON_LABEL = "E. coli K-12 MG1655"
_WRITE = re.compile(r"\b(CREATE|MERGE|SET|DELETE|REMOVE|DROP|LOAD\s+CSV|FOREACH|CALL\s*\{|apoc\.(create|merge|periodic))\b", re.I)


def cypher(query, **params):
    """Run a READ-ONLY Cypher query, return list of row dicts."""
    if _WRITE.search(query):
        raise ValueError("write keyword in query - this analysis is read-only")
    body = json.dumps({"statements": [{"statement": query, "parameters": params}]}).encode()
    req = urllib.request.Request(NEO4J_URL, body, {
        "Content-Type": "application/json",
        "Authorization": "Basic " + base64.b64encode(b"neo4j:").decode()})
    with urllib.request.urlopen(req, timeout=600) as r:
        out = json.load(r)
    if out["errors"]:
        raise RuntimeError(out["errors"])
    res = out["results"][0]
    return [dict(zip(res["columns"], d["row"])) for d in res["data"]]


def load_model():
    import cobra
    return cobra.io.read_sbml_model(str(MODEL_XML))


def uid(met_id):
    """iML1515 metabolite id 'atp_c' -> compartment-free 'atp' (BiGG: compartment after the last underscore)."""
    return met_id.rsplit("_", 1)[0]


def xref_ns(obj):
    """COBRApy annotation dict -> {namespace: [ids]} (always lists)."""
    return {k: ([v] if isinstance(v, str) else list(v)) for k, v in obj.annotation.items()}
