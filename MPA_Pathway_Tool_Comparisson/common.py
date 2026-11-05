"""Shared helpers for the MPA_Pathway_Tool vs. knowledge-graph comparison.

All knowledge-graph access goes through `cypher()`, which only allows read-only
queries against the Neo4j HTTP transaction endpoint (user neo4j, empty password).
"""
import base64, csv, json, os, re, sys, urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
DATA = BASE / "data"                      # inputs (paper supplements, pathway files): shared by all KG versions

# The KG endpoint is mandatory and comes from the environment (set by organizer.sh from its port argument).
NEO4J_URL = os.getenv("NEO4J_URL")
_m = re.match(r"^https?://[^:/]+:(\d+)/", NEO4J_URL or "")
if not _m:
    sys.exit("NEO4J_URL is not set or has no port. Run the pipeline via: ./organizer.sh <port>")
VERSION_KG = f"KG_port{_m.group(1)}"      # one output folder per KG version
OUT = BASE / VERSION_KG
RESULTS, FIGURES = OUT / "results", OUT / "figures"
for p in (RESULTS, FIGURES):
    p.mkdir(parents=True, exist_ok=True)

_WRITE = re.compile(r"\b(CREATE|MERGE|SET|DELETE|REMOVE|DROP|LOAD\s+CSV|FOREACH|CALL\s*\{|apoc\.(create|merge|periodic))\b", re.I)

# pathway -> (KEGG module in KG, taxonomic requirement used in the paper, supp. file prefix of paper per-reaction table)
PATHWAYS = {
    "Wood-Ljungdahl":              dict(module="KEGG:M00377", taxa=["Viruses", "Unknown", "Bacteria"],
                                        kegg_csv="KEGG-WoodLjungdahl pathway.csv", supp1="WoodLjungdahlKEGGModule.csv",
                                        supp_rx="Supp. file 13_"),
    "Hydrogenotrophic methanogenesis": dict(module="KEGG:M00567", taxa=["Archaea"],
                                        kegg_csv="KEGG-hydrogenotrophic methanogenesis.csv", supp1="hydrogenotrophicMethanogenesisKEGGModule.csv",
                                        supp_rx="Supp. file 14_"),
    "Acetoclastic methanogenesis": dict(module="KEGG:M00357", taxa=["Archaea"],
                                        kegg_csv="KEGG-acetoclastic methanogenesis.csv", supp1="acetoclasticMethanogenesisKEGGModule.csv",
                                        supp_rx="Supp. file 15_"),
}


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


def load_proteins():
    """Supp. 9 metaproteomics table -> DataFrame (+ list of the 44 sample names)."""
    path = next(DATA.glob("Supp. file 9_*"))
    df = pd.read_csv(path, sep="\t", dtype=str, quoting=csv.QUOTE_NONE, encoding="utf-8", encoding_errors="replace")
    samples = list(df.columns[12:])
    df[samples] = df[samples].astype(float)
    return df, samples


def _split(s):
    return [] if s in (None, "", "-") or pd.isna(s) else [x.strip() for x in str(s).split("|") if x.strip() and x.strip() != "-"]


def protein_params(df):
    """Protein records for Cypher: i = row index, ids = ['KEGG:Kxxxxx', 'EC:x.x.x.x'], sk = superkingdom."""
    recs = []
    for i, (ko, ec, sk) in enumerate(zip(df.iloc[:, 1], df.iloc[:, 2], df["superkingdom"])):
        ids = ["KEGG:" + k for k in _split(ko)] + ["EC:" + e for e in _split(ec)]
        if ids:
            recs.append(dict(i=i, ids=ids, sk=sk))
    return recs


MATCH_PARAM_Q = """
UNWIND $reactions AS r
UNWIND $prots AS p
WITH r, p WHERE p.sk IN $taxa AND any(x IN p.ids WHERE x IN r.ids)
RETURN r.rid AS rid, collect(p.i) AS idx
"""

MATCH_KG_Q = """
MATCH (m:Module {id: $module})-[:HAS_REACTION]->(r:Reaction)
OPTIONAL MATCH (r)-[:CATALYZED_BY]->(x) WHERE x:KeggOrthology OR x:Enzyme
WITH r, collect(DISTINCT x.id) AS rids
UNWIND $prots AS p
WITH r, rids, p WHERE p.sk IN $taxa AND any(x IN p.ids WHERE x IN rids)
RETURN r.id AS rid, collect(p.i) AS idx
"""


def counts_from_matches(matches, df, samples, rids_all=()):
    """matches: {rid: [protein row idx]} -> per-reaction sample sums + module sums.
    Module count = sum over reactions (a protein matching several reactions counts once per reaction, as in the paper)."""
    X = df[samples].to_numpy()
    rx = {rid: X[idx].sum(axis=0) if idx else np.zeros(len(samples)) for rid, idx in matches.items()}
    for rid in rids_all:
        rx.setdefault(rid, np.zeros(len(samples)))
    mod_multi = sum(rx.values()) if rx else np.zeros(len(samples))
    return rx, mod_multi


def parse_pathway_csv(path):
    """Pathway-Creator CSV -> list of dict(rid, ko[], ec[], taxa[])  (one per step)."""
    rows = list(csv.reader(open(path, encoding="utf-8", errors="replace"), delimiter=";"))[1:]
    steps = {}
    for r in rows:
        if r[0] in steps:
            continue
        m = re.search(r"\b(R\d{5})\s*$", r[1])
        steps[r[0]] = dict(rid="KEGG:" + m.group(1) if m else r[1],
                           ko=[k for k in r[2].split(",") if k], ec=[e for e in r[3].split(",") if e],
                           taxa=[t.split(":")[1] for t in r[8].split("&&") if ":" in t])
    return list(steps.values())


def load_supp1():
    """Supp. 1 (Calculator result table): rows = pathways; 4 info columns + 44 sample sums. Header in file is broken."""
    path = next(DATA.glob("Supp. file 1_*"))
    rows = list(csv.reader(open(path)))[1:]
    return {r[0]: dict(found=float(r[2]), total=float(r[3]), counts=np.array(r[4:], dtype=float)) for r in rows}


def load_supp_reactions(prefix):
    """Supp. 10-15 per-reaction mapping table -> {KEGG reaction id: counts[44]} + KO/EC lists as given by the paper."""
    path = next(DATA.glob(prefix + "*"))
    lines = path.read_text(encoding="utf-8", errors="replace").split("\n")[1:]
    out = {}
    for ln in lines:
        f = ln.split(";")
        if len(f) < 48:
            continue
        m = re.search(r"\b(R\d{5})\s*$", f[0].strip()) or re.search(r"\b(R\d{5})\b", f[0])
        rid = "KEGG:" + m.group(1) if m else f[0].strip()
        out[rid] = dict(ko=[k for k in f[1].split(",") if k], ec=[e for e in f[2].split(",") if e],
                        counts=np.array(f[-44:], dtype=float), name=f[0].strip())
    return out
