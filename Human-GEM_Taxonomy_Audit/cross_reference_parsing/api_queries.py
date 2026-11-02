"""Every web-API call used to read the taxonomy of a Human-GEM entity, in one place (counterpart of ../kg_queries.py).

One function per API endpoint, all going through common.fetch (throttled, retried, cached). Sections:
  KEGG     rest.kegg.jp          get R/C/K entries (batch of 10)  -> EC numbers, KO ids, organisms listed under a KO
  UNIPROT  rest.uniprot.org      organism of accessions; does an EC / Rhea id occur in organism 9606 (or in any organism)?
  RHEA     www.rhea-db.org       reactions (+EC) that use a ChEBI compound
  BIGG     bigg.ucsd.edu         organisms of the BiGG models that contain a reaction / metabolite, EC links of a reaction
  NCBI     eutils.ncbi.nlm.nih.gov  taxid of NCBI Gene ids      ENSEMBL  rest.ensembl.org  species of ENSG ids
  HGNC     rest.genenames.org    existence of a gene symbol (HGNC = human genes only)      REACTOME  reactome.org  species of an event
Outcome vocabulary of every route: 'human' (9606 among the taxa) | 'other' (taxa found, 9606 not among them) | 'none' (id known to
nobody / no taxon information) | 'noid' (the entity has no id of that type).
"""
import json, re, urllib.parse
import common as c

EC_RE = re.compile(r"\b\d+\.\d+\.\d+\.(?:\d+|n\d+|-)\b")
UNI = "https://rest.uniprot.org/uniprotkb"


def verdict(taxa_human, taxa_any):
    return "human" if taxa_human else ("other" if taxa_any else "none")


# ------------------------------------------------------------------ KEGG
def kegg_get(ids):
    """GET /get/<id>+<id>... (max 10). -> {id: {'ecs': [...], 'kos': [...], 'orgs': [KEGG organism codes under GENES]}}"""
    s, _, body = c.fetch("https://rest.kegg.jp/get/" + "+".join(ids))
    out = {}
    if s != 200:
        return out
    for entry in body.split("\n///"):
        if not entry.strip():
            continue
        sec, cur, eid = {}, None, None
        for line in entry.split("\n"):
            if line[:1].strip():                                   # new section header
                cur = line.split()[0]
                sec[cur] = [line[12:]]
                if cur == "ENTRY":
                    eid = line.split()[1]
            elif cur:
                sec[cur].append(line[12:])
        if eid:
            text = lambda k: "\n".join(sec.get(k, []))
            out[eid] = {"ecs": sorted(set(EC_RE.findall(text("ENZYME")))),
                        "kos": re.findall(r"\bK\d{5}\b", text("ORTHOLOGY")),
                        "orgs": sorted({l.split(":")[0].strip() for l in sec.get("GENES", []) if re.match(r"\s*[A-Z0-9]{2,6}:", l)})}
    return out


# ------------------------------------------------------------------ UNIPROT
def uniprot_count(query):
    """Number of UniProtKB entries matching a query (X-Total-Results header, size=1); None if the API kept failing."""
    s, h, _ = c.fetch(f"{UNI}/search?query={urllib.parse.quote(query)}&size=1&fields=accession")
    return int(h["x-total-results"]) if s == 200 and "x-total-results" in h else None


def uniprot_route(term):
    """term e.g. 'ec:2.7.1.40' or '(rhea:1 OR rhea:2)': human if an entry from organism 9606 exists, other if any organism has one."""
    n_h = uniprot_count(f"{term} AND organism_id:{c.HUMAN}")
    if n_h:
        return "human"
    if n_h is None:
        return "none"
    n_all = uniprot_count(term)
    return "other" if n_all else "none"


def uniprot_ec(ec):
    return uniprot_route(f"ec:{ec}")


def rhea_block(rid):
    """Rhea ids come in groups of four (master, left->right, right->left, bidirectional); UniProt cites one of the directional ones,
    models/KEGG cite whichever they like, so the whole block is tried."""
    base = int(rid) - int(rid) % 4
    return sorted({int(rid)} | set(range(base, base + 4)))


def uniprot_rhea(rid):
    return uniprot_route("(" + " OR ".join(f"rhea:{x}" for x in rhea_block(rid)) + ")")


def uniprot_organisms(accessions):
    """Batch (<=300) /accessions -> {accession: organism taxid}. Accessions UniProt does not know are absent."""
    s, _, body = c.fetch(f"{UNI}/accessions?accessions={','.join(accessions)}&fields=accession,organism_id&format=tsv")
    if s != 200:
        return {}
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in body.strip().split("\n")[1:] if "\t" in l and l.split("\t")[1].isdigit()}


# ------------------------------------------------------------------ RHEA
def rhea_ecs_of_chebi(chebi):
    """All Rhea reactions that have this ChEBI compound as participant -> sorted EC numbers of those reactions (None on failure)."""
    num = chebi.split(":")[-1]
    s, _, body = c.fetch(f"https://www.rhea-db.org/rhea?query=chebi:{num}&columns=rhea-id,ec&format=tsv&limit=500")
    if s != 200:
        return None
    ecs = set()
    for l in body.strip().split("\n")[1:]:
        ecs.update(EC_RE.findall(l))
    return sorted(ecs)


# ------------------------------------------------------------------ BIGG
def bigg(kind, bid):
    """kind 'reactions' | 'metabolites' -> {'orgs': [organism names of models containing it], 'ecs': [EC from database_links]} or None if unknown."""
    s, _, body = c.fetch(f"https://bigg.ucsd.edu/api/v2/universal/{kind}/{urllib.parse.quote(bid)}")
    if s != 200:
        return None
    j = json.loads(body)
    models = j.get("models_containing_reaction") or j.get("compartments_in_models") or []
    links = (j.get("database_links") or {}).get("EC Number", [])
    return {"orgs": sorted({m["organism"] for m in models}), "ecs": sorted({l["id"] for l in links})}


# ------------------------------------------------------------------ NCBI / ENSEMBL / HGNC / REACTOME
def ncbi_gene_taxids(entrez_ids):
    """Batch (<=200) esummary -> {entrez: taxid}."""
    s, _, body = c.fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&retmode=json&id=" + ",".join(entrez_ids))
    if s != 200:
        return {}
    res = json.loads(body).get("result", {})
    return {u: int(res[u]["organism"]["taxid"]) for u in res.get("uids", []) if res.get(u, {}).get("organism")}


def ensembl_species(ensg_ids):
    """Batch POST /lookup/id (<=1000) -> {ENSG: species name}."""
    s, _, body = c.fetch("https://rest.ensembl.org/lookup/id", data=json.dumps({"ids": ensg_ids}),
                         headers={"Content-Type": "application/json", "Accept": "application/json"})
    return {k: v["species"] for k, v in json.loads(body).items() if v} if s == 200 else {}


def hgnc_symbol_exists(symbol):
    s, _, body = c.fetch(f"https://rest.genenames.org/fetch/symbol/{urllib.parse.quote(symbol)}", headers={"Accept": "application/json"})
    return s == 200 and json.loads(body)["response"]["numFound"] > 0


def reactome_species(stid):
    s, _, body = c.fetch(f"https://reactome.org/ContentService/data/query/{stid}/speciesName")
    return body.strip() if s == 200 else None
