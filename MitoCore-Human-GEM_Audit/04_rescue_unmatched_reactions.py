"""Take the MitoCore reactions that 03 could not match to human KG evidence, try alternative matches, and check whether the KG
links them to OTHER taxa.

Rescue attempts (each is its own column, none is merged into the strict match of 03):
  copy          same compartment-free stoichiometry as a reaction that 03 matched (a compartment copy)      [proxy for human]
  ec_hier       partial EC (e.g. EC:2.4.1.-) with a human enzyme below it in the EC hierarchy (Enzyme-CHILD_OF->Enzyme)
  participants  all non-ubiquitous compounds of the reaction (KEGG ids) equal the compound set of one KG Reaction node,
                which has a strong human route (R1/R2)                                                        [proxy for human]
  gene          at least one GPR gene is human in the KG (G1): supports the reaction but does not identify it
Other taxa: proteins of non-human taxa annotated with the reaction's EC (Enzyme) or KEGG reaction (via KO), grouped by Protein.taxon_group.
"""
import pandas as pd
import common as c
import kg_queries as q
import plotstyle as ps

R = c.RESULTS
mr = pd.read_csv(R / "match_reactions.csv").fillna("")
mm = pd.read_csv(R / "match_metabolites.csv").fillna("")
kr = pd.read_csv(R / "KG_human_reactions.csv").fillna("")
split = lambda s: [x for x in s.split(";") if x] if isinstance(s, str) else []
UBIQ = {"KEGG:C00001", "KEGG:C00080"}                       # water, proton (not decisive for reaction identity)

un = mr[mr.status != "matched_human"].copy()
print(f"unmatched reactions: {len(un)} of {len(mr)}")

# ---- copy: signature shared with a matched reaction ----
matched_sigs = set(mr.signature[mr.status == "matched_human"])
un["copy"] = un.signature.isin(matched_sigs)

# ---- ec_hier: partial EC with a human descendant ----
enz = pd.DataFrame(c.cypher(q.REACTIONS["E_NODES"]))
parent = dict(zip(enz.id, enz.parent))
human_ec = {x["id"] for x in c.cypher(q.REACTIONS["E_human"], h=c.HUMAN)}
anc = set()
for e in human_ec:
    while isinstance(parent.get(e), str) and parent[e] not in anc:
        e = parent[e]
        anc.add(e)
un["ec_hier"] = un.ec.map(lambda s: any(e.endswith("-") and e in anc for e in split(s)))

# ---- participants: compound-set equality with a KG reaction ----
kgc = {x["id"]: set(x["compounds"]) - UBIQ for x in c.cypher(q.REACTIONS["R_compounds"])}
by_set = {}
for rid, s in kgc.items():
    if len(s) >= 2:
        by_set.setdefault(frozenset(s), []).append(rid)
uid2k = {r.unique_id: set(split(r.kegg_compound)) for r in mm.itertuples()}
strong = set(kr.id[kr.strong_human.astype(bool)])


def participants(r):
    if r.reaction_class not in ("metabolic", "isolated"):
        return []
    ks = set()
    for uid in split(r.uids):
        k = uid2k.get(uid, set())
        if k and k <= UBIQ:
            continue
        if not k:
            return []                                        # a non-ubiquitous metabolite has no KEGG id -> set incomplete
        ks |= k
    return by_set.get(frozenset(ks - UBIQ), []) if len(ks - UBIQ) >= 2 else []


un["part_cands"] = un.apply(participants, axis=1)
un["participants"] = un.part_cands.map(lambda l: any(x in strong for x in l))
un["part_cands"] = un.part_cands.map(";".join)
un["gene"] = un.gene_human.astype(bool)

# ---- other taxa ----
ecs = sorted({e for s in un.ec for e in split(s)} & set(enz.id))
kegg_ids = sorted({k for s in un.kegg_reaction for k in split(s)} | {k for s in un.part_cands for k in split(s)})


def batched(name, ids, n=100):
    out = []
    for i in range(0, len(ids), n):
        out += c.cypher(q.REACTIONS[name], ids=ids[i:i + n], h=c.HUMAN)
        print(f"  {name}: {min(i + n, len(ids))}/{len(ids)}", flush=True)
    return pd.DataFrame(out, columns=["id", "taxon_group", "n_proteins", "n_taxa"])


e_tax = batched("E_other_taxa", ecs)
k_tax = batched("R_ko_other_taxa", kegg_ids)
e_tax.to_csv(R / "other_taxa_by_ec.csv", index=False)
k_tax.to_csv(R / "other_taxa_by_kegg_reaction.csv", index=False)
groups_ec = e_tax.groupby("id").taxon_group.agg(lambda x: set(x.dropna()))
groups_k = k_tax.groupby("id").taxon_group.agg(lambda x: set(x.dropna()))
prot_ec = e_tax.groupby("id").n_proteins.sum()
prot_k = k_tax.groupby("id").n_proteins.sum()


def other(r):
    g, n = set(), 0
    for e in split(r.ec):
        g |= groups_ec.get(e, set()); n += int(prot_ec.get(e, 0))
    for k in split(r.kegg_reaction) + split(r.part_cands):
        g |= groups_k.get(k, set()); n += int(prot_k.get(k, 0))
    return pd.Series(dict(other_taxa_groups=";".join(sorted(g)), other_taxa_proteins=n))


un[["other_taxa_groups", "other_taxa_proteins"]] = un.apply(other, axis=1)
un["other_taxa"] = un.other_taxa_groups != ""

# ---- final category (first applicable) ----
CAT = ["rescued_human_proxy", "other_taxa_only", "gene_supported_only", "no_evidence"]
NICE = {"rescued_human_proxy": "human by proxy (copy / EC hierarchy / participants)", "other_taxa_only": "KG knows it in other taxa only",
        "gene_supported_only": "only a human GPR gene supports it", "no_evidence": "no KG evidence"}
COL = {"rescued_human_proxy": "#2a78d6", "other_taxa_only": "#eb6834", "gene_supported_only": "#1baf7a", "no_evidence": "#9a9a95"}


def final(r):
    if r["copy"] or r.ec_hier or r.participants: return "rescued_human_proxy"
    if r.other_taxa: return "other_taxa_only"
    if r.gene: return "gene_supported_only"
    return "no_evidence"


un["final_category"] = un.apply(final, axis=1)
cols = ["id", "name", "subsystem", "reaction_class", "compartments", "status", "kegg_reaction", "ec", "genes", "copy", "ec_hier",
        "participants", "part_cands", "gene", "other_taxa", "other_taxa_groups", "other_taxa_proteins", "final_category"]
un[cols].to_csv(R / "unmatched_reactions.csv", index=False)

# ---- summaries ----
ct = pd.crosstab(un.reaction_class, un.final_category).reindex(columns=CAT, fill_value=0)
ct["total"] = ct.sum(axis=1)
ct.to_csv(R / "rescue_summary.csv")
print("\nrescue flags among unmatched:", {k: int(un[k].sum()) for k in ("copy", "ec_hier", "participants", "gene", "other_taxa")})
print("\n", ct.to_string())
tg = un[un.other_taxa].other_taxa_groups.map(split).explode().value_counts()
tg.rename_axis("taxon_group").reset_index(name="n_reactions").to_csv(R / "other_taxa_groups_summary.csv", index=False)
print("\nreactions with other-taxa evidence by taxon group:\n", tg.to_string())
ms = un[un.reaction_class.isin(["metabolic", "isolated"]) & (un.final_category == "no_evidence")]
print("\nmetabolic/isolated reactions left with no evidence:", len(ms),
      "| of which no KEGG id and no EC:", int(((ms.kegg_reaction == "") & (ms.ec == "")).sum()))
print("top subsystems of the residual no-evidence reactions:\n", un[un.final_category == "no_evidence"].subsystem.value_counts().head(8).to_string())

# ---- plots ----
plt = ps.plt
order = [x for x in ["metabolic", "isolated", "transport", "exchange_demand", "artificial"] if x in ct.index]
c2 = ct.loc[order]
fig, ax = plt.subplots(figsize=(8.5, 3.8))
ps.stacked_h(ax, c2.index.tolist(), c2[CAT].values, CAT, totals=c2.total.tolist(), colors=COL, nice=NICE)
ax.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=2, bbox_to_anchor=(0, -0.2))
fig.suptitle("What the KG can say about the MitoCore reactions that did not match human evidence", x=0.01, ha="left", fontsize=12)
ps.finish(fig, "fig5_unmatched_rescue.png")

fig, ax = plt.subplots(figsize=(7.5, 3.8))
t = tg.sort_values()
ax.barh(range(len(t)), t.values, color="#2a78d6", height=0.62)
for i, v in enumerate(t.values):
    ax.text(v, i, f" {v:,}", va="center", fontsize=9, color=ps.INK)
ax.set_yticks(range(len(t)))
ax.set_yticklabels(t.index)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.set_xlabel("unmatched reactions with KG evidence in this taxon group (a reaction can appear in several)")
fig.suptitle("Unmatched MitoCore reactions: where else does the KG place them?", x=0.01, ha="left", fontsize=12)
ps.finish(fig, "fig6_unmatched_other_taxa.png")
