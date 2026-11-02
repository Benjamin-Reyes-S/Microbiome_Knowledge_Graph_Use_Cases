# MicrobiomeKG use cases

Three read-only analyses that test the MicrobiomeKG (Neo4j) against external references, or against what you get without it.

| Folder | Question |
|---|---|
| `MPA_Pathway_Tool_Comparisson/` | Do the MPA_Pathway_Tool results (Walke et al. 2021) agree with the KG? |
| `Human-GEM_Taxonomy_Audit/` | Which Human-GEM reactions, genes and metabolites does the KG connect to *Homo sapiens* (`NCBITaxon:9606`)? |
| `Human-GEM_Taxonomy_Audit/cross_reference_parsing/` | How far do chained public web APIs (KEGG, UniProt, Rhea, BiGG, ...) get on the same question, without the KG? |

## Requirements
- **KG** (pipelines 1 and 2): Neo4j HTTP endpoint `http://localhost:8004/db/neo4j/tx/commit`, user `neo4j`, empty password. Queries are read-only
  (`common.cypher` rejects write keywords). Check it is up before running: `curl -s localhost:8004`.
- **Python:** the MPA pipeline needs pandas, numpy, matplotlib, scipy. The Human-GEM pipelines need COBRApy as well and use the virtual
  environment `Human-GEM_Taxonomy_Audit/.env/` (`Human-GEM_Taxonomy_Audit/.env/bin/python`).
- **Internet** (pipeline 3 only): access to rest.kegg.jp, rest.uniprot.org, www.rhea-db.org, bigg.ucsd.edu, eutils.ncbi.nlm.nih.gov,
  rest.ensembl.org, rest.genenames.org and reactome.org.
- Original inputs are never modified. Each pipeline writes only to its own `results/` and `figures/`.

## 1. MPA pathway tool vs. KG (`MPA_Pathway_Tool_Comparisson/`)
Inputs are copies in `data/`. Run from inside the folder, in order:
1. `01_kg_module_definitions.py` - KEGG modules M00377 / M00567 / M00357 from the KG, diffed against the paper's pathway files
2. `02_calibrate_replication.py` - reproduces the paper's Calculator numbers (Supp. 1, 10-15) from the paper's own pathway files
3. `03_kg_mapping.py` - maps Supp. 9 proteins onto the KG modules (matching done in Cypher)
4. `04_compare_and_plot.py` - module-level, then reaction-level comparison and figures
5. `05_explain_R04347_difference.py` - traces the only count difference to a protein

```bash
cd MPA_Pathway_Tool_Comparisson
for s in 0[1-5]_*.py; do python "$s"; done
```
Outputs: `results/*.csv`, `figures/*.png`. The method write-up is `WALKTHROUGH.md`.

## 2. Human-GEM taxonomy audit against the KG (`Human-GEM_Taxonomy_Audit/`)
Reads `Human-GEM.xml` (Human-GEM 2.1.1). All Cypher is in `kg_queries.py` (sections REACTIONS / GENES / METABOLITES); the scripts only run it.
Run from inside the folder, in order:
1. `01_extract_kg_human.py` - route queries -> `results/KG_human_{reactions,genes,metabolites}.csv` (one row per KG node, one column per route)
2. `02_extract_human_gem.py` - the model -> `results/Human_GEM_{reactions,genes,metabolites}.csv` (reaction class, compartment-free signature, ids)
3. `03_match_model_to_kg.py` - strict match (`matched_human` / `mapped_not_human` / `unmappable`), unmatched by class and subsystem, fig1-3
4. `04_rescue_unmatched_reactions.py` - rescue by compartment copy, EC hierarchy, participants, GPR gene; other-taxa evidence; fig5-6
5. `05_xref_profile_unmatched.py` - which identifier types the unmatched entities carry (reactions fig7, metabolites fig9); needs only the model and the 03 tables, not the KG

```bash
cd Human-GEM_Taxonomy_Audit
for s in 0[1-5]_*.py; do .env/bin/python -W ignore -u "$s"; done
```
Strong human routes are R1 (EC), R2 (KO), G1 (UniProt), M1 (via human reaction); weak or supporting routes are always reported separately.
Context and the verified KG schema are in `SKILL.md`, the work plan in `PLAN.md`.

## 3. Cross-reference parsing without the KG (`Human-GEM_Taxonomy_Audit/cross_reference_parsing/`)
Classifies the same entities to human by chaining public APIs, one parser per identifier type. Every API call is in `api_queries.py`; all requests
go through the throttled, retrying `common.fetch()` and are cached in `data/api_cache.sqlite`, so re-runs are free and a failed run can be restarted.
Run from inside the folder, in order (the full run sends about 15,000 requests and took roughly two hours, mostly UniProt):
1. `01_extract_model_xrefs.py` - every identifier of the model -> `results/model_xrefs_*.csv`
2. `02_parse_genes.py` - routes GX1-GX4 (UniProt, NCBI Gene, Ensembl, HGNC) -> `results/api_genes.csv`
3. `03_parse_reactions.py` - routes RX1-RX7 (KEGG, EC, Rhea, BiGG, Reactome) -> `results/api_reactions.csv`
4. `04_parse_metabolites.py` - routes MX1-MX4 (KEGG, ChEBI via Rhea, BiGG, HMDB) -> `results/api_metabolites.csv`
5. `05_patterns_and_figures.py` - classification, mapping patterns, comparison with the KG tables of pipeline 2, fig1-5

```bash
cd Human-GEM_Taxonomy_Audit/cross_reference_parsing
for s in 0[1-5]_*.py; do ../.env/bin/python -W ignore -u "$s"; done
```
Pipeline 2 must have been run first: step 05 reads its `match_*.csv`, and step 01 reads its `Human_GEM_reactions.csv`.
An entity is human if any route reaches 9606. The BiGG-models and HMDB routes only state membership in a human model/database, so results are
reported both with and without them (`results/taxonomy_summary.csv`). Give-ups after retries are counted in `results/api_call_log.csv`.

## Reading the results
- Results are CSVs in each `results/`, figures in each `figures/`.
- "Not mapped" and "no KG evidence" are findings, not errors: transport, exchange and pseudo-reactions of a metabolic model are expected to lack KEGG ids.
- Do not mix numbers across pipelines unless the folder says they are comparable (`Human-GEM_Taxonomy_Audit/legacy/` was produced against a wrong KG).
