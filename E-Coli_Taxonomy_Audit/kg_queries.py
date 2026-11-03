"""Every Cypher query used against the knowledge graph, in one place, in three sections:
REACTIONS, GENES, METABOLITES.  Scripts import these strings and run them unchanged (read-only, see common.cypher).

Taxon scope (defined in common.py): $taxa = strain + K-12 (STRICT: 511145, 83333), $taxa_sp = + species (562, WIDER).
Every route returns n_support (strict scope) and n_support_species (wider scope); the strict count drives the match.
Target model: E. coli K-12 MG1655 (iML1515).
Pattern of each section:  *_NODES  -> every KG node of the type (+ matching ids/xrefs), one row per node
                          route R*/G*/M*  -> nodes connected to the target taxon by that route (id [, n_support])
Routes are never merged; strong = R1, R2, G1, G2, M1a, M1b.  Protein is the only node type with a direct taxon edge.
"""

# =====================================================================================================================
# SECTION 1: REACTIONS   (Reaction nodes have no taxon edge; taxon evidence comes via Enzyme / KO / Module-Pathway)
# =====================================================================================================================
REACTIONS = {}

# every Reaction node with its EC and KO annotations (collected lists)
REACTIONS["R_NODES"] = """
MATCH (r:Reaction)
OPTIONAL MATCH (r)-[:CATALYZED_BY]->(e:Enzyme)
WITH r, collect(DISTINCT e.id) AS ecs
OPTIONAL MATCH (r)-[:CATALYZED_BY]->(k:KeggOrthology)
RETURN r.id AS id, r.name AS name, r.equation AS equation, ecs, collect(DISTINCT k.id) AS kos
"""

# R1: reaction -> enzyme (EC) <- protein of the taxon          [strong]  n = distinct proteins
REACTIONS["R1_enzyme_protein"] = """
MATCH (r:Reaction)-[:CATALYZED_BY]->(:Enzyme)<-[:HAS_EC]-(p:Protein)-[:BELONGS_TO]->(t:Taxon) WHERE t.id IN $taxa_sp
RETURN r.id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN p END) AS n_support, count(DISTINCT p) AS n_support_species
"""

# R2: reaction -> ortholog group (KO) <- protein of the taxon    [strong]
REACTIONS["R2_ko_protein"] = """
MATCH (r:Reaction)-[:CATALYZED_BY]->(:KeggOrthology)<-[:MEMBER_OF]-(p:Protein)-[:BELONGS_TO]->(t:Taxon) WHERE t.id IN $taxa_sp
RETURN r.id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN p END) AS n_support, count(DISTINCT p) AS n_support_species
"""

# R3: reaction <- module -> pathway <- taxon         [weak, supporting only]
REACTIONS["R3_module_pathway"] = """
MATCH (r:Reaction)<-[:HAS_REACTION]-(:Module)-[:PART_OF]->(pw:Pathway)<-[:HAS_PATHWAY]-(t:Taxon) WHERE t.id IN $taxa_sp
RETURN r.id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN pw END) AS n_support, count(DISTINCT pw) AS n_support_species
"""
# (the human route R4 'KO <- KEGG hsa gene' has no E. coli counterpart: no Gene node carries a KEGG:eco: id in this KG)

# --- matching of model reactions (script 03) -------------------------------------------------------------------------
# every Enzyme node with its parent in the EC hierarchy (EC:1.1.1.1 -> EC:1.1.1.-)
REACTIONS["E_NODES"] = """
MATCH (e:Enzyme)
OPTIONAL MATCH (e)-[:CHILD_OF]->(pa:Enzyme)
RETURN e.id AS id, pa.id AS parent
"""

# enzymes (EC) with at least one protein of the taxon          n = distinct proteins
REACTIONS["E_taxon"] = """
MATCH (e:Enzyme)<-[:HAS_EC]-(p:Protein)-[:BELONGS_TO]->(t:Taxon) WHERE t.id IN $taxa_sp
RETURN e.id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN p END) AS n_support, count(DISTINCT p) AS n_support_species
"""

# --- evidence for rescue / other taxa (script 04); batch with $ids ----------------------------------------------------
# compounds of every Reaction node (to match a model reaction by its participants)
REACTIONS["R_compounds"] = """
MATCH (r:Reaction)
OPTIONAL MATCH (r)-[:HAS_REACTANT|HAS_PRODUCT]->(x:Compound)
RETURN r.id AS id, collect(DISTINCT x.id) AS compounds
"""

# EC ids -> taxon groups of the proteins annotated with them, EXCLUDING the target taxa (other-taxa evidence)
REACTIONS["E_other_taxa"] = """
UNWIND $ids AS i
MATCH (e:Enzyme {id:i})<-[:HAS_EC]-(p:Protein)-[:BELONGS_TO]->(t:Taxon)
WHERE NOT t.id IN $taxa_sp
RETURN i AS id, p.taxon_group AS taxon_group, count(DISTINCT p) AS n_proteins, count(DISTINCT t) AS n_taxa
"""

# KEGG reaction ids -> taxon groups of proteins via KO, EXCLUDING the target taxa
REACTIONS["R_ko_other_taxa"] = """
UNWIND $ids AS i
MATCH (r:Reaction {id:i})-[:CATALYZED_BY]->(:KeggOrthology)<-[:MEMBER_OF]-(p:Protein)-[:BELONGS_TO]->(t:Taxon)
WHERE NOT t.id IN $taxa_sp
RETURN i AS id, p.taxon_group AS taxon_group, count(DISTINCT p) AS n_proteins, count(DISTINCT t) AS n_taxa
"""

# =====================================================================================================================
# SECTION 2: GENES   (iML1515 genes: locus tag bNNNN, Entrez id, UniProt accession)
# KG Gene nodes are NCBIGene:<entrez> (with KEGG:eco:/EnsemblGenome: locus tags in xrefs) and have NO taxon edge;
# Gene -ENCODES-> Protein -BELONGS_TO-> Taxon, and Transcript -BELONGS_TO-> Taxon (Transcript.gene_id = locus tag).
# The expanded KG has ~8.4M Gene nodes, so all gene queries are restricted to the model's ids.
# =====================================================================================================================
GENES = {}

# the Gene nodes for the model's Entrez ids ($gene_ids = ['NCBIGene:945001', ...]) with symbols, xrefs, locus tags
GENES["G_NODES"] = """
MATCH (g:Gene) WHERE g.id IN $gene_ids
RETURN g.id AS id, g.symbols AS symbols, g.xrefs AS xrefs
"""

# G1: model UniProt accession -> Protein node with taxon edge ($uniprots = ['UniProtKB:P0A825', ...])        [strong]
GENES["G1_uniprot_protein"] = """
UNWIND $uniprots AS u
MATCH (p:Protein {id:u})-[:BELONGS_TO]->(t:Taxon) WHERE t.id IN $taxa_sp
RETURN u AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN t END) AS n_support, count(DISTINCT t) AS n_support_species
"""

# G2: Gene node -ENCODES-> Protein with taxon edge                                                           [strong]
GENES["G2_encodes_protein"] = """
MATCH (g:Gene) WHERE g.id IN $gene_ids
MATCH (g)-[:ENCODES]->(p:Protein)-[:BELONGS_TO]->(t:Taxon) WHERE t.id IN $taxa_sp
RETURN g.id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN p END) AS n_support, count(DISTINCT p) AS n_support_species
"""

# G3: Transcript of the taxon whose gene_id is the model's locus tag ($locus = ['b0001', ...])               [supporting]
GENES["G3_transcript_locus"] = """
MATCH (tr:Transcript)-[:BELONGS_TO]->(t:Taxon) WHERE t.id IN $taxa_sp AND tr.gene_id IN $locus
RETURN tr.gene_id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN tr END) AS n_support, count(DISTINCT tr) AS n_support_species
"""

# Proteins by UniProt accession, in the taxon or not (used to map model uniprot ids); batch with $ids
GENES["P_by_uniprot"] = """
UNWIND $ids AS u
MATCH (p:Protein {id:'UniProtKB:'+u})
OPTIONAL MATCH (p)-[:BELONGS_TO]->(t:Taxon)
RETURN u AS uniprot, p.id AS id, collect(DISTINCT t.id) AS taxa
"""

# =====================================================================================================================
# SECTION 3: METABOLITES   (Compound has only id, name, synonyms, sources: NO ChEBI, NO taxon edge -> evidence is derived)
# =====================================================================================================================
METABOLITES = {}

# every Compound node
METABOLITES["M_NODES"] = """
MATCH (x:Compound)
RETURN x.id AS id, x.name AS name, x.synonyms AS synonyms
"""

# M1a: compound in a reaction that is in the taxon via R1 (enzyme route)   n = distinct supporting reactions   [strong, derived]
METABOLITES["M1a_via_R1"] = """
MATCH (x:Compound)<-[:HAS_REACTANT|HAS_PRODUCT]-(r:Reaction)-[:CATALYZED_BY]->(:Enzyme)<-[:HAS_EC]-(p:Protein)-[:BELONGS_TO]->(t:Taxon) WHERE t.id IN $taxa_sp
RETURN x.id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN r END) AS n_support, count(DISTINCT r) AS n_support_species
"""

# M1b: compound in a reaction that is in the taxon via R2 (KO route)       [strong, derived]
METABOLITES["M1b_via_R2"] = """
MATCH (x:Compound)<-[:HAS_REACTANT|HAS_PRODUCT]-(r:Reaction)-[:CATALYZED_BY]->(:KeggOrthology)<-[:MEMBER_OF]-(p:Protein)-[:BELONGS_TO]->(t:Taxon) WHERE t.id IN $taxa_sp
RETURN x.id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN r END) AS n_support, count(DISTINCT r) AS n_support_species
"""

# M2: compound participates in a pathway linked to the taxon            [weak, supporting only]
METABOLITES["M2_pathway"] = """
MATCH (x:Compound)-[:PARTICIPATES_IN]->(pw:Pathway)<-[:HAS_PATHWAY]-(t:Taxon) WHERE t.id IN $taxa_sp
RETURN x.id AS id, count(DISTINCT CASE WHEN t.id IN $taxa THEN pw END) AS n_support, count(DISTINCT pw) AS n_support_species
"""

# compounds in no reaction at all (no evidence possible through M1)
METABOLITES["M_no_reaction"] = """
MATCH (x:Compound) WHERE NOT (x)<-[:HAS_REACTANT|HAS_PRODUCT]-(:Reaction)
RETURN x.id AS id
"""
