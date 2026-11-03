"""Every Cypher query used against the knowledge graph, in one place, in three sections:
REACTIONS, GENES, METABOLITES.  Scripts import these strings and run them unchanged (read-only, see common.cypher).

Parameter  $h = 'NCBITaxon:9606' (Homo sapiens).
Pattern of each section:  *_NODES  -> every KG node of the type (+ matching ids/xrefs), one row per node
                          route R*/G*/M*  -> nodes connected to human by that route (id [, n_support])
Routes are never merged; strong = R1, R2, G1, M1a, M1b.  Protein is the only node type with a direct taxon edge.
"""

# =====================================================================================================================
# SECTION 1: REACTIONS   (Reaction nodes have no taxon edge; human evidence comes via Enzyme / KO / Module-Pathway)
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

# R1: reaction -> enzyme (EC) <- human protein            [strong]  n = distinct human proteins
REACTIONS["R1_enzyme_protein"] = """
MATCH (r:Reaction)-[:CATALYZED_BY]->(:Enzyme)<-[:HAS_EC]-(p:Protein)-[:BELONGS_TO]->(:Taxon {id:$h})
RETURN r.id AS id, count(DISTINCT p) AS n_support
"""

# R2: reaction -> ortholog group (KO) <- human protein    [strong]
REACTIONS["R2_ko_protein"] = """
MATCH (r:Reaction)-[:CATALYZED_BY]->(:KeggOrthology)<-[:MEMBER_OF]-(p:Protein)-[:BELONGS_TO]->(:Taxon {id:$h})
RETURN r.id AS id, count(DISTINCT p) AS n_support
"""

# R3: reaction <- module -> pathway <- human taxon         [weak, supporting only]
REACTIONS["R3_module_pathway"] = """
MATCH (r:Reaction)<-[:HAS_REACTION]-(:Module)-[:PART_OF]->(pw:Pathway)<-[:HAS_PATHWAY]-(:Taxon {id:$h})
RETURN r.id AS id, count(DISTINCT pw) AS n_support
"""

# R4: reaction -> KO <- KEGG hsa gene (human by id prefix, no taxon edge)   [supporting]
REACTIONS["R4_ko_hsa_gene"] = """
MATCH (r:Reaction)-[:CATALYZED_BY]->(:KeggOrthology)<-[:MEMBER_OF]-(g:Gene)
WHERE g.id STARTS WITH 'KEGG:hsa:'
RETURN r.id AS id, count(DISTINCT g) AS n_support
"""

# --- matching of model reactions (script 03) -------------------------------------------------------------------------
# every Enzyme node with its parent in the EC hierarchy (EC:1.1.1.1 -> EC:1.1.1.-)
REACTIONS["E_NODES"] = """
MATCH (e:Enzyme)
OPTIONAL MATCH (e)-[:CHILD_OF]->(pa:Enzyme)
RETURN e.id AS id, pa.id AS parent
"""

# enzymes (EC) with at least one human protein          n = distinct human proteins
REACTIONS["E_human"] = """
MATCH (e:Enzyme)<-[:HAS_EC]-(p:Protein)-[:BELONGS_TO]->(:Taxon {id:$h})
RETURN e.id AS id, count(DISTINCT p) AS n_support
"""

# --- evidence for rescue / other taxa (script 04); batch with $ids ----------------------------------------------------
# compounds of every Reaction node (to match a model reaction by its participants)
REACTIONS["R_compounds"] = """
MATCH (r:Reaction)
OPTIONAL MATCH (r)-[:HAS_REACTANT|HAS_PRODUCT]->(x:Compound)
RETURN r.id AS id, collect(DISTINCT x.id) AS compounds
"""

# EC ids -> taxon groups of the proteins annotated with them, EXCLUDING human (other-taxa evidence)
REACTIONS["E_other_taxa"] = """
UNWIND $ids AS i
MATCH (e:Enzyme {id:i})<-[:HAS_EC]-(p:Protein)-[:BELONGS_TO]->(t:Taxon)
WHERE t.id <> $h
RETURN i AS id, p.taxon_group AS taxon_group, count(DISTINCT p) AS n_proteins, count(DISTINCT t) AS n_taxa
"""

# KEGG reaction ids -> taxon groups of proteins via KO, EXCLUDING human
REACTIONS["R_ko_other_taxa"] = """
UNWIND $ids AS i
MATCH (r:Reaction {id:i})-[:CATALYZED_BY]->(:KeggOrthology)<-[:MEMBER_OF]-(p:Protein)-[:BELONGS_TO]->(t:Taxon)
WHERE t.id <> $h
RETURN i AS id, p.taxon_group AS taxon_group, count(DISTINCT p) AS n_proteins, count(DISTINCT t) AS n_taxa
"""

# =====================================================================================================================
# SECTION 2: GENES   (Gene = KEGG:hsa:<entrez>; no taxon edge and no ENCODES edge; link to Protein goes by UniProt xref)
# =====================================================================================================================
GENES = {}

# NOTE: the expanded KG has ~8.4M Gene nodes (LocusTag, NCBIGene, EMBLCDS, KEGG), so every gene query is restricted
# to the model's gene ids ($ids = the KEGG:hsa:<entrez> ids in Human_GEM_genes.csv, written by 02).

# the Gene nodes for the model's ids, with symbols, UniProt/RefSeq xrefs and KO memberships
GENES["G_NODES"] = """
MATCH (g:Gene) WHERE g.id IN $ids
OPTIONAL MATCH (g)-[:MEMBER_OF]->(k:KeggOrthology)
RETURN g.id AS id, g.symbols AS symbols, g.xrefs AS xrefs, collect(DISTINCT k.id) AS kos
"""

# G1: gene xref UniProtKB:<acc> -> Protein with taxon edge to human                [strong]
GENES["G1_uniprot_protein"] = """
MATCH (g:Gene) WHERE g.id IN $ids
UNWIND [x IN g.xrefs WHERE x STARTS WITH 'UniProtKB:'] AS u
MATCH (p:Protein {id:u})-[:BELONGS_TO]->(:Taxon {id:$h})
RETURN g.id AS id, count(DISTINCT p) AS n_support
"""

# G2: gene -> KO <- human protein                                                    [supporting; KO membership is broad]
GENES["G2_ko_protein"] = """
MATCH (g:Gene) WHERE g.id IN $ids
MATCH (g)-[:MEMBER_OF]->(:KeggOrthology)<-[:MEMBER_OF]-(p:Protein)-[:BELONGS_TO]->(:Taxon {id:$h})
RETURN g.id AS id, count(DISTINCT p) AS n_support
"""

# G3: Gene id has the KEGG human prefix hsa (organism code, not a taxon edge)        [supporting, edge_type = id-prefix]
GENES["G3_hsa_prefix"] = """
MATCH (g:Gene) WHERE g.id IN $ids AND g.id STARTS WITH 'KEGG:hsa:'
RETURN g.id AS id
"""

# check, expected 0 in this KG: direct taxon edge on Gene
GENES["G_direct_taxon_edge_check"] = """
MATCH (g:Gene) WHERE g.id IN $ids
MATCH (g)-[r]-(:Taxon {id:$h})
RETURN type(r) AS edge_type, count(DISTINCT g) AS n
"""

# Proteins by UniProt accession, human or not (used to map model uniprot ids); batch with $ids
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

# M1a: compound in a reaction that is human via R1 (enzyme route)   n = distinct supporting reactions   [strong, derived]
METABOLITES["M1a_via_R1"] = """
MATCH (x:Compound)<-[:HAS_REACTANT|HAS_PRODUCT]-(r:Reaction)-[:CATALYZED_BY]->(:Enzyme)<-[:HAS_EC]-(p:Protein)-[:BELONGS_TO]->(:Taxon {id:$h})
RETURN x.id AS id, count(DISTINCT r) AS n_support
"""

# M1b: compound in a reaction that is human via R2 (KO route)       [strong, derived]
METABOLITES["M1b_via_R2"] = """
MATCH (x:Compound)<-[:HAS_REACTANT|HAS_PRODUCT]-(r:Reaction)-[:CATALYZED_BY]->(:KeggOrthology)<-[:MEMBER_OF]-(p:Protein)-[:BELONGS_TO]->(:Taxon {id:$h})
RETURN x.id AS id, count(DISTINCT r) AS n_support
"""

# M2: compound participates in a pathway linked to human            [weak, supporting only]
METABOLITES["M2_pathway"] = """
MATCH (x:Compound)-[:PARTICIPATES_IN]->(pw:Pathway)<-[:HAS_PATHWAY]-(:Taxon {id:$h})
RETURN x.id AS id, count(DISTINCT pw) AS n_support
"""

# compounds in no reaction at all (no evidence possible through M1)
METABOLITES["M_no_reaction"] = """
MATCH (x:Compound) WHERE NOT (x)<-[:HAS_REACTANT|HAS_PRODUCT]-(:Reaction)
RETURN x.id AS id
"""
