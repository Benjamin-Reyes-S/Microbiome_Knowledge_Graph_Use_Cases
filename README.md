# MicrobiomeKG use cases

Read-only analyses that test the MicrobiomeKG (Neo4j) against external references, or against what you get without it. The same audit
is run for three metabolic models and can be repeated against **different KG versions** (one Neo4j instance per version, selected by port),
to see what changes when the KG is expanded.

| Folder | Question |
|---|---|
| `MPA_Pathway_Tool_Comparisson/` | Do the MPA_Pathway_Tool results (Walke et al. 2021) agree with the KG? |
| `Human-GEM_Taxonomy_Audit/` | Which Human-GEM reactions, genes and metabolites does the KG connect to *Homo sapiens* (`NCBITaxon:9606`)? |
| `E-Coli_Taxonomy_Audit/` | Same audit for the gold-standard *E. coli* model iML1515 against *E. coli* K-12 MG1655 (`NCBITaxon:511145`) |
| `MitoCore-Human-GEM_Audit/` | Same audit for the small curated human mitochondrial model MitoCore against *Homo sapiens* |
| `Human-GEM_Taxonomy_Audit/cross_reference_parsing/` | How far do chained public web APIs (KEGG, UniProt, Rhea, BiGG, ...) get on the same question, without the KG? |

## Requirements
- **KG** (pipelines 1, 2, 4, 5): Neo4j HTTP endpoint `http://localhost:<port>/db/neo4j/tx/commit`, user `neo4j`, empty password. The port is
  chosen at run time (see "Running against a KG version"). Queries are read-only (`common.cypher` rejects write keywords).
  Check the instance answers before running: `curl -s -m 10 -u neo4j: -H 'Content-Type: application/json' -d '{"statements":[{"statement":"RETURN 1"}]}' localhost:<port>/db/neo4j/tx/commit`
  (a listening port can still hang while Neo4j is starting or the machine is under load).
- **Python:** the MPA pipeline needs pandas, numpy, matplotlib, scipy. The Human-GEM pipelines need COBRApy as well and use the virtual
  environment `Human-GEM_Taxonomy_Audit/.env/` (`Human-GEM_Taxonomy_Audit/.env/bin/python`).
- **Internet** (pipeline 3 only): access to rest.kegg.jp, rest.uniprot.org, www.rhea-db.org, bigg.ucsd.edu, eutils.ncbi.nlm.nih.gov,
  rest.ensembl.org, rest.genenames.org and reactome.org.
- Original inputs are never modified. The MPA pipeline writes to its own `results/` and `figures/`; the model audits write to one folder per KG version (below).

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

## Running against a KG version (pipelines 2, 4, 5)
Each model folder has `organizer.sh`, which takes the Neo4j HTTP port as its only argument, activates the shared virtual environment, exports
`NEO4J_URL` and runs the scripts in the right order:
```bash
cd Human-GEM_Taxonomy_Audit      # or E-Coli_Taxonomy_Audit / MitoCore-Human-GEM_Audit
./organizer.sh 8034              # e.g. the expanded KG; run again with another port for another KG version
```
- `common.py` requires `NEO4J_URL` (it exits with a message if it is missing) and derives `VERSION_KG = KG_port<port>`.
- All outputs go to `<model folder>/KG_port<port>/{data,results,figures}`, so runs against different KG versions never overwrite each other.
- `Human-GEM_Taxonomy_Audit/legacy/` holds the first run (before versioning, original port 8004) and is kept as is.
- Order is `02 -> 01 -> 03 -> 04 -> 05`: the expanded KG has about 8.4 million `Gene` nodes (LocusTag, NCBIGene, EMBLCDS, KEGG), so `01` looks up
  only the model's own gene ids, which `02` writes. A full-table gene query exceeded Neo4j's 8 GiB transaction limit.

## 2. Human-GEM taxonomy audit against the KG (`Human-GEM_Taxonomy_Audit/`)
Reads `Human-GEM.xml` (Human-GEM 2.1.1). All Cypher is in `kg_queries.py` (sections REACTIONS / GENES / METABOLITES); the scripts only run it.
`./organizer.sh <port>` runs, in order:
1. `02_extract_human_gem.py` - the model -> `results/Human_GEM_{reactions,genes,metabolites}.csv` (reaction class, compartment-free signature, ids)
2. `01_extract_kg_human.py` - route queries -> `results/KG_human_{reactions,genes,metabolites}.csv` (reactions and metabolites: one row per KG node; genes: the model's gene ids only)
3. `03_match_model_to_kg.py` - strict match (`matched_human` / `mapped_not_human` / `unmappable`), unmatched by class and subsystem, fig1-3
4. `04_rescue_unmatched_reactions.py` - rescue by compartment copy, EC hierarchy, participants, GPR gene; other-taxa evidence; fig5-6
5. `05_xref_profile_unmatched.py` - which identifier types the unmatched entities carry (reactions fig7, metabolites fig9); needs only the model and the 03 tables, not the KG

Strong human routes are R1 (EC), R2 (KO), G1 (UniProt), M1 (via human reaction); weak or supporting routes are always reported separately.
Context and the verified KG schema are in `SKILL.md`, the work plan in `PLAN.md`.

## 2b. E. coli iML1515 audit (`E-Coli_Taxonomy_Audit/`)
A copy of pipeline 2 adapted to the E. coli K-12 MG1655 model (`iML1515.xml`, 2,712 reactions, 1,516 genes, 1,877 metabolites). Run `./organizer.sh <port>`.
- **Taxon scope** (`common.py`): in the KG most MG1655 proteins sit on the parent K-12 node (83333) or the species (562), not on the strain node (511145;
  only 50 proteins). Strict scope = strain + K-12 (`TAXA`); species hits are reported separately (`*_species` columns, `species_only` flag) and not counted.
- **Genes:** G1 = model UniProt accession -> Protein in the taxon (strong); G2 = `NCBIGene` node -ENCODES-> such a Protein (strong); G3 = Transcript by locus tag (supporting).
  The human route R4 (KO <- `KEGG:hsa:` gene) has no counterpart and is dropped.
- **Model parsing:** BiGG ids with `_c/_p/_e` suffixes (`common.uid()`), no subsystem labels, so reaction classes come from ids and structure; no subsystem table or fig3.
- Status names: `matched_taxon` / `mapped_not_taxon` / `unmappable`. Scripts: `01_extract_kg_taxon.py`, `02_extract_model.py`, `03`-`05` as in pipeline 2.
- Known gap in the KG tested (port 8034): the KO route R2 returns nothing (the strain's proteins have no KO links), so all reaction matches come from EC.

## 2c. MitoCore audit (`MitoCore-Human-GEM_Audit/`)
Pipeline 2 for `Mitocore_Human_GEM.xml` (560 reactions, 593 genes, 445 metabolites; compartments cytosol / mitochondrion / external). Target taxon is human (9606).
Differences from Human-GEM: BiGG-style ids (`common.uid()`), reaction classes from the model's own subsystems (biomass / `OF_*` / "artificial" -> artificial;
`EX_` / `DM_` / "Boundary ..." -> exchange_demand; "transport" subsystems -> transport), no `isolated` class. Run `./organizer.sh <port>`.

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
**Not updated for versioned outputs:** this pipeline still reads `../results/`, which was moved to `../legacy/results/`. It is finished and left untouched; point
`PARENT_RESULTS` in its `common.py` at `legacy/results` (or a `KG_port*/results` folder) before re-running it.
An entity is human if any route reaches 9606. The BiGG-models and HMDB routes only state membership in a human model/database, so results are
reported both with and without them (`results/taxonomy_summary.csv`). Give-ups after retries are counted in `results/api_call_log.csv`.

## Reading the results
- Results are CSVs in each `results/`, figures in each `figures/`.
- "Not mapped" and "no KG evidence" are findings, not errors: transport, exchange and pseudo-reactions of a metabolic model are expected to lack KEGG ids.
- Do not mix numbers across pipelines or KG versions unless the folder says they are comparable (`Human-GEM_Taxonomy_Audit/legacy/` was produced against a wrong KG).
- Results of the model audits are per KG version: compare `KG_port<a>/results/match_summary.csv` with `KG_port<b>/results/match_summary.csv`.
