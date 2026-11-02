"""Classify every reaction, gene and metabolite by taxonomy (human / other taxa only / not mappable) from the API routes of 02-04,
define the MAPPING PATTERN of each entity (= the set of parsing routes that reach human) and compare with the knowledge graph.
Rule: an entity counts as human if ANY route reaches NCBITaxon 9606, even when the same route also returns other taxa.
Not mappable = no route returned any taxon; split into 'no id of a queryable type' vs 'ids present but every API came back empty'.
Outputs (results/): patterns_{reactions,genes,metabolites}.csv, route_coverage.csv, taxonomy_summary.csv, kg_vs_api.csv, api_by_reaction_class.csv
Figures: fig1-3 patterns per component (A: % per mapping pattern, B: what each single route achieves), fig4 API parsing vs KG,
         fig5 how many independent routes support each human call.
"""
import pandas as pd
import numpy as np
import common as c
import plotstyle as ps

plt = ps.plt
R, P = c.RESULTS, c.PARENT_RESULTS
LABEL = {"GX1_uniprot": "UniProt", "GX2_ncbigene": "NCBI Gene", "GX3_ensembl": "Ensembl", "GX4_hgnc": "HGNC",
         "RX1_kegg_ec_uniprot": "KEGG>EC>UniProt", "RX2_kegg_ko_genes": "KEGG>KO>genes", "RX3_ec_uniprot": "EC>UniProt",
         "RX4_rhea_uniprot": "Rhea>UniProt", "RX5_bigg_models": "BiGG>models", "RX6_bigg_ec_uniprot": "BiGG>EC>UniProt",
         "RX7_reactome_species": "Reactome",
         "MX1_kegg_enzyme_uniprot": "KEGG>ENZYME>UniProt", "MX2_chebi_rhea_uniprot": "ChEBI>Rhea>EC>UniProt",
         "MX3_bigg_models": "BiGG>models", "MX4_hmdb_implied_human": "HMDB (implied)"}
# routes that only say "this id belongs to a human model / human database" (membership, not enzyme or protein evidence):
# Human-GEM's own HMR ids are registered in BiGG's Recon3D, HMDB is the human metabolome. Reported as a separate, stricter tier.
WEAK = {"RX5_bigg_models", "MX3_bigg_models", "MX4_hmdb_implied_human"}
# an entity with none of these ids cannot even be queried
QUERYABLE = {"reactions": ["has_kegg", "has_ec", "has_rhea", "has_bigg", "has_reactome"],
             "genes": [], "metabolites": ["has_kegg", "has_chebi", "has_bigg", "has_hmdb"]}
COMP = {"reactions": ("api_reactions.csv", "id", "match_reactions.csv", "id"),
        "genes": ("api_genes.csv", "id", "match_genes.csv", "id"),
        "metabolites": ("api_metabolites.csv", "unique_id", "match_metabolites.csv", "unique_id")}
TITLE = {"reactions": "reactions", "genes": "genes", "metabolites": "metabolites (unique compounds)"}

summary, cover, kgcmp, data = [], [], [], {}
for comp, (f, key, kgf, kgkey) in COMP.items():
    d = pd.read_csv(R / f).fillna("")
    routes = [x for x in d.columns if x[:2] in ("GX", "RX", "MX")]
    H = d[routes].eq("human")
    O = d[routes].eq("other")
    strong = [r for r in routes if r not in WEAK]
    d["n_human_routes"] = H.sum(axis=1)
    d["tax_strict"] = np.where(H[strong].any(axis=1), "human", np.where(O[strong].any(axis=1), "other", "unmappable"))
    d["tax"] = np.where(H.any(axis=1), "human", np.where(O.any(axis=1), "other", "unmappable"))
    q = QUERYABLE[comp]
    can = d[q].any(axis=1) if q else pd.Series(True, index=d.index)
    d["pattern"] = [" + ".join(LABEL[r] for r in routes if H.loc[i, r]) if t == "human" else
                    ("other taxa only" if t == "other" else
                     ("no id of a queryable type" if not can[i] else "ids present, no taxon returned"))
                    for i, t in zip(d.index, d.tax)]
    data[comp] = d
    pat = d.pattern.value_counts().rename_axis("pattern").reset_index(name="n")
    pat["pct"] = (100 * pat.n / len(d)).round(2)
    pat.to_csv(R / f"patterns_{comp}.csv", index=False)
    for r in routes:
        vc = d[r].value_counts()
        cover.append(dict(component=comp, route=r, label=LABEL[r], total=len(d), **{k: int(vc.get(k, 0)) for k in ("human", "other", "none", "noid")}))
    vc = d.tax.value_counts()
    summary.append(dict(component=comp, total=len(d), human=int(vc.get("human", 0)), other_only=int(vc.get("other", 0)),
                        unmappable=int(vc.get("unmappable", 0)), pct_human=round(100 * vc.get("human", 0) / len(d), 1),
                        pct_other_only=round(100 * vc.get("other", 0) / len(d), 1), pct_unmappable=round(100 * vc.get("unmappable", 0) / len(d), 1),
                        unmappable_no_queryable_id=int((d.pattern == "no id of a queryable type").sum()),
                        unmappable_ids_but_empty=int((d.pattern == "ids present, no taxon returned").sum()),
                        pct_human_strict=round(100 * (d.tax_strict == "human").mean(), 1),
                        pct_other_only_strict=round(100 * (d.tax_strict == "other").mean(), 1),
                        pct_unmappable_strict=round(100 * (d.tax_strict == "unmappable").mean(), 1),
                        n_patterns=int((d.tax == "human").any() and d[d.tax == "human"].pattern.nunique())))
    kg = pd.read_csv(P / kgf)[[kgkey, "status"]].rename(columns={kgkey: key, "status": "kg_status"})
    j = d.merge(kg, on=key, how="left")
    a, k = j.tax == "human", j.kg_status == "matched_human"
    kgcmp.append(dict(component=comp, total=len(j), kg_human=int(k.sum()), api_human=int(a.sum()), both=int((a & k).sum()),
                      only_api=int((a & ~k).sum()), only_kg=int((~a & k).sum()), neither=int((~a & ~k).sum()),
                      kg_mapped_not_human=int((j.kg_status == "mapped_not_human").sum()), kg_unmappable=int((j.kg_status == "unmappable").sum())))
    if comp == "reactions":
        cl = pd.crosstab(d.reaction_class, d.tax).reindex(columns=["human", "other", "unmappable"], fill_value=0)
        cl["total"] = cl.sum(axis=1)
        cl["pct_human"] = (100 * cl.human / cl.total).round(1)
        cl.to_csv(R / "api_by_reaction_class.csv")
        print(cl.to_string())

summary, cover, kgcmp = pd.DataFrame(summary), pd.DataFrame(cover), pd.DataFrame(kgcmp)
summary.to_csv(R / "taxonomy_summary.csv", index=False); cover.to_csv(R / "route_coverage.csv", index=False); kgcmp.to_csv(R / "kg_vs_api.csv", index=False)
print(summary.to_string(index=False)); print(kgcmp.to_string(index=False))

# ---------------- fig1-3: mapping patterns per component ----------------
TOPN = 10
OUT_COLOR = {"human": ps.BLUE, "other taxa only": ps.ORANGE, "no id of a queryable type": ps.GRAY, "ids present, no taxon returned": "#c4c3be"}
for n, comp in enumerate(COMP, 1):
    d, tot = data[comp], len(data[comp])
    pat = pd.read_csv(R / f"patterns_{comp}.csv")
    hum = pat[~pat.pattern.isin(OUT_COLOR)]
    top = hum.head(TOPN)
    rows = [(p, k, "human") for p, k in zip(top.pattern, top.n)]
    if len(hum) > TOPN:
        rows.append((f"{len(hum) - TOPN} other combinations", int(hum.n.iloc[TOPN:].sum()), "human"))
    for p in ("other taxa only", "ids present, no taxon returned", "no id of a queryable type"):
        k = int(pat.n[pat.pattern == p].sum())
        if k:
            rows.append((p, k, p))
    cv = cover[cover.component == comp]
    fig, (a, b) = plt.subplots(1, 2, figsize=(15, max(4.2, 0.45 * max(len(rows), len(cv)) + 1.8)), gridspec_kw={"width_ratios": [1.25, 1]})
    for i, (p, k, kind) in enumerate(rows[::-1]):
        a.barh(i, 100 * k / tot, color=OUT_COLOR.get(kind, ps.BLUE), height=0.65)
        a.text(100 * k / tot, i, f"  {100 * k / tot:.1f}%  ({k:,})", va="center", fontsize=8.5, color=ps.INK2)
    a.set_yticks(range(len(rows))); a.set_yticklabels([r[0] for r in rows[::-1]], fontsize=8.5)
    a.set_xlim(0, max(100 * r[1] / tot for r in rows) * 1.3); a.set_xlabel(f"% of {TITLE[comp]} (n={tot:,})")
    a.set_title("A  mapping pattern = routes that reach human\n   (blue = human, orange = other taxa only, gray = not mappable)", loc="left", fontsize=10, color=ps.INK)
    # panel B: what each single route achieves
    cols = [("human", ps.BLUE, "human"), ("other", ps.ORANGE, "other taxa only"), ("none", "#c4c3be", "id known to nobody / empty"),
            ("noid", ps.RULE, "entity has no such id")]
    left = np.zeros(len(cv))
    for col, colr, nm in cols:
        v = 100 * cv[col].values / tot
        b.barh(range(len(cv)), v, left=left, color=colr, height=0.62, edgecolor=ps.SURFACE, linewidth=1.5, label=nm)
        for i, x in enumerate(v):
            if x >= 6:
                b.text(left[i] + x / 2, i, f"{x:.0f}%", ha="center", va="center", fontsize=8.5, color="white" if col in ("human", "other") else ps.INK2)
        left += v
    b.set_yticks(range(len(cv))); b.set_yticklabels([l + ("  [membership]" if r in WEAK else "") for l, r in zip(cv.label, cv.route)], fontsize=9); b.invert_yaxis(); b.set_xlim(0, 100)
    b.set_xlabel("% of entities, per single parsing route")
    b.set_title("B  outcome of each route on its own", loc="left", fontsize=10, color=ps.INK)
    b.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0, -0.13), ncol=2)
    for ax in (a, b):
        ps.style(ax)
    fig.suptitle(f"Taxonomy of Human-GEM {TITLE[comp]} by parsing public APIs: {summary.pct_human[n - 1]:.1f}% human, "
                 f"{summary.pct_unmappable[n - 1]:.1f}% not mappable", x=0.01, ha="left", fontsize=12.5)
    ps.finish(fig, f"fig{n}_patterns_{comp}.png")

# ---------------- fig4: API parsing vs knowledge graph ----------------
fig, ax = plt.subplots(figsize=(10.5, 5.6))
names, vals = [], []
for comp, s, k in zip(COMP, summary.itertuples(), kgcmp.itertuples()):
    nm, d_ = TITLE[comp].split(" (")[0], data[comp]
    ns = len(d_.filter(like=("GX" if comp == "genes" else "RX" if comp == "reactions" else "MX")).columns)
    names += [f"{nm}: KG", f"{nm}: API, enzyme/protein evidence only", f"{nm}: API, all {ns} routes (incl. model membership)"]
    st = d_.tax_strict.value_counts()
    vals += [[k.kg_human, k.kg_mapped_not_human, k.kg_unmappable], [st.get("human", 0), st.get("other", 0), st.get("unmappable", 0)],
             [s.human, s.other_only, s.unmappable]]
v = np.array(vals, float); pct = 100 * v / v.sum(axis=1, keepdims=True); left = np.zeros(len(v))
y = np.array([0, 1, 2, 3.4, 4.4, 5.4, 6.8, 7.8, 8.8])[::-1]
for j, (col, nm) in enumerate(((ps.BLUE, "human"), (ps.ORANGE, "other taxa only / in KG without human evidence"), (ps.GRAY, "not mappable / no id in KG"))):
    ax.barh(y, pct[:, j], left=left, color=col, height=0.8, edgecolor=ps.SURFACE, linewidth=2, label=nm)
    for yy, x, l, n_ in zip(y, pct[:, j], left, v[:, j]):
        if x >= 6:
            ax.text(l + x / 2, yy, f"{x:.0f}%", ha="center", va="center", color="white", fontsize=9)
    left += pct[:, j]
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=9); ax.set_xlim(0, 100); ax.set_xlabel("% of entities")
ax.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0, -0.14), ncol=3); ps.style(ax)
ax.set_title("Taxonomy assignment: integrated KG vs. chains of per-database API parsers", loc="left", fontsize=12)
ps.finish(fig, "fig4_api_vs_kg.png")

# ---------------- fig5: number of independent routes behind each human call ----------------
fig, ax = plt.subplots(figsize=(8.5, 3.6))
cat = ["not mappable / other only", "1 route", "2 routes", "3+ routes"]
colors = [ps.GRAY, "#9cc3ef", "#5a9be0", ps.BLUE]
rows = []
for comp in COMP:
    d = data[comp]; nr = d.n_human_routes
    rows.append([(nr == 0).sum(), (nr == 1).sum(), (nr == 2).sum(), (nr >= 3).sum()])
rows = np.array(rows, float); pc = 100 * rows / rows.sum(axis=1, keepdims=True); left = np.zeros(3)
for j in range(4):
    ax.barh(range(3), pc[:, j], left=left, color=colors[j], height=0.62, edgecolor=ps.SURFACE, linewidth=2, label=cat[j])
    for i, x in enumerate(pc[:, j]):
        if x >= 5:
            ax.text(left[i] + x / 2, i, f"{x:.0f}%", ha="center", va="center", color="white", fontsize=9)
    left += pc[:, j]
ax.set_yticks(range(3)); ax.set_yticklabels([f"{TITLE[k].split(' (')[0]} (n={len(data[k]):,})" for k in COMP]); ax.invert_yaxis(); ax.set_xlim(0, 100)
ax.set_xlabel("% of entities, by number of parsing routes that reach human"); ps.style(ax)
ax.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0, -0.2), ncol=4)
ax.set_title("A human call often rests on a single database chain", loc="left", fontsize=12)
ps.finish(fig, "fig5_route_support.png")
