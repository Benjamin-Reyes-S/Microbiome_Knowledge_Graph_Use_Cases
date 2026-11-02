"""Shared helpers for the Human-GEM taxonomy audit against the knowledge graph.

All knowledge-graph access goes through `cypher()`, which only allows read-only
queries against the Neo4j HTTP transaction endpoint (user neo4j, empty password).
"""
import base64, json, re, urllib.request
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent
DATA = BASE / "data"
RESULTS = BASE / "results"
FIGURES = BASE / "figures"
for p in (DATA, RESULTS, FIGURES):
    p.mkdir(exist_ok=True)

MODEL_XML = BASE / "Human-GEM.xml"
HUMAN = "NCBITaxon:9606"
NEO4J_URL = "http://localhost:8004/db/neo4j/tx/commit"
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


def xref_ns(obj):
    """COBRApy annotation dict -> {namespace: [ids]} (always lists)."""
    return {k: ([v] if isinstance(v, str) else list(v)) for k, v in obj.annotation.items()}
