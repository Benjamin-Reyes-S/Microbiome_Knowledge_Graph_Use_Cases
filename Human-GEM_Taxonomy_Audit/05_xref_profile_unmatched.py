"""Cross-reference profile of the Human-GEM entities that did NOT match human KG evidence (output of 03).

Question: are the unmatched reactions / genes / metabolites unmatched only because they lack a KEGG (or EC) identifier,
and which OTHER identifier types do they carry that the KG could be linked to (improve the KG) or used to match them?

Needs no KG access: it reads the model annotations (identifiers.org namespaces) and the match tables written by 03.
Outputs (results/):  xref_coverage_{reactions,metabolites}.csv   one row per namespace x group (% carrying it)
                     xref_patterns_{reactions,metabolites}.csv   one row per combination of namespaces (unmatched only)
                     xref_unmatched_{reactions,metabolites}.csv  one row per unmatched entity with its namespaces + rescue_potential
                     xref_summary.csv                                  headline: unmatched without KG-usable id but with another id
Figures: fig7_xref_reactions.png, fig9_xref_metabolites.png  (A: coverage per namespace, B: top id combinations)
KG-usable namespaces = ids that the strict match in 03 can use today (reactions: kegg.reaction, ec-code; metabolites: kegg.compound;
genes are not profiled: 2,869/2,870 match). Everything else is a candidate for a new KG cross-reference / mapping table.
"""
import numpy as np
import pandas as pd
import common as c
import plotstyle as ps

R = c.RESULTS
plt = ps.plt
KG_USABLE = {"reactions": {"kegg.reaction", "ec-code"}, "metabolites": {"kegg.compound"}}

# ---------------- namespaces carried by each model entity ----------------
model = c.load_model()
# not external identifiers: `sbo` is an ontology term present on every entity, `vmh*` is the model's own id namespace (always present)
NOT_XREF = {"sbo", "vmhreaction", "vmhmetabolite"}
ns = lambda o: set(c.xref_ns(o)) - NOT_XREF
ns_r = {r.id: ns(r) for r in model.reactions}
ns_m = {}                                           # unique metabolite = compartment-collapsed id (union over the copies)
for m in model.metabolites:
    ns_m.setdefault(m.id[:-1] if m.compartment and m.id.endswith(m.compartment) else m.id, set()).update(ns(m))

mr = pd.read_csv(R / "match_reactions.csv").fillna("")
mm = pd.read_csv(R / "match_metabolites.csv").fillna("")            # one row per unique metabolite

# non-model-annotation ids that the model carries as plain fields are already inside the annotation; name/formula are not ids
mr["ns"] = mr.id.map(ns_r); mm["ns"] = mm.unique_id.map(ns_m)

# groups compared in panel A: matched (reference) vs. unmatched, reactions split by type
mr["group"] = np.where(mr.status == "matched_human", "matched to human",
              np.where(mr.reaction_class.isin(["metabolic", "isolated"]), "unmatched: metabolic / isolated",
              np.where(mr.reaction_class == "transport", "unmatched: transport", "unmatched: exchange / artificial")))
mm["group"] = np.where(mm.status == "matched_human", "matched to human", "unmatched")
GROUP_COLOR = {"matched to human": "#2a78d6", "unmatched": "#eb6834", "unmatched: metabolic / isolated": "#eb6834",
               "unmatched: transport": "#7a5cc4", "unmatched: exchange / artificial": "#9a9a95"}

summary = []


def profile(name, d, ns_col="ns"):
    """Write the three CSVs and return (coverage table, pattern table) for one component."""
    usable = KG_USABLE[name]
    d = d.copy()
    allns = sorted(set().union(*d[ns_col]))
    unm = d[d.status != "matched_human"].copy()

    # coverage: % of the group's entities that carry each namespace
    cov = []
    for g, sub in d.groupby("group"):
        for n in allns:
            k = sum(n in s for s in sub[ns_col])
            cov.append(dict(component=name, group=g, namespace=n, n=k, group_size=len(sub), pct=round(100 * k / len(sub), 1),
                            kg_usable=n in usable))
    cov = pd.DataFrame(cov)
    cov.to_csv(R / f"xref_coverage_{name}.csv", index=False)

    # per unmatched entity: which ids it has, and what that means for fixing it
    def potential(s):
        if s & usable:
            return "has KG-usable id (mapped, but no human route in KG)"
        if s:
            return "no KG-usable id, but another id type (new KG xref / mapping would rescue)"
        return "no identifier at all (only name / stoichiometry / structure)"
    unm["namespaces"] = unm[ns_col].map(lambda s: ";".join(sorted(s)))
    unm["rescue_potential"] = unm[ns_col].map(potential)
    keep = [x for x in ("id", "unique_id", "name", "symbol", "subsystem", "reaction_class", "status") if x in unm.columns]
    unm[keep + ["namespaces", "rescue_potential"]].to_csv(R / f"xref_unmatched_{name}.csv", index=False)

    # combinations of namespaces among the unmatched
    pat = unm.namespaces.replace("", "(no identifier)").value_counts().rename_axis("pattern").reset_index(name="n")
    pat["pct"] = (100 * pat.n / max(len(unm), 1)).round(1)
    pat.to_csv(R / f"xref_patterns_{name}.csv", index=False)

    vc = unm.rescue_potential.value_counts()
    for k, v in vc.items():
        summary.append(dict(component=name, unmatched=len(unm), rescue_potential=k, n=int(v), pct=round(100 * v / len(unm), 1)))
    if name == "reactions":                       # the same split for the reactions that are NOT exchange/transport/artificial
        core = unm[unm.reaction_class.isin(["metabolic", "isolated"])]
        for k, v in core.rescue_potential.value_counts().items():
            summary.append(dict(component="reactions (metabolic/isolated only)", unmatched=len(core), rescue_potential=k,
                                n=int(v), pct=round(100 * v / len(core), 1)))
    return cov, pat, len(unm)


def plot(name, cov, pat, n_unm, title, fname, topn=12):
    order = (cov[cov.group != "matched to human"].groupby("namespace").n.sum().sort_values().index.tolist())   # most common unmatched last
    groups = [g for g in GROUP_COLOR if g in set(cov.group)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, max(4, 0.42 * len(order) + 1.6)), gridspec_kw={"width_ratios": [1.05, 1]})
    h = 0.8 / len(groups)
    for j, g in enumerate(groups):
        v = cov[cov.group == g].set_index("namespace").pct.reindex(order).fillna(0)
        size = int(cov[cov.group == g].group_size.iloc[0])
        a.barh(np.arange(len(order)) + (j - (len(groups) - 1) / 2) * h, v.values, height=h * 0.9, color=GROUP_COLOR[g],
               label=f"{g} (n={size:,})")
    a.set_yticks(range(len(order)))
    a.set_yticklabels([o + ("  *" if o in KG_USABLE[name] else "") for o in order])
    a.set_xlim(0, 100); a.set_xlabel("% of entities carrying the identifier type")
    a.set_title("A  identifier types carried (* = usable by the KG match today)", loc="left", fontsize=10.5, color=ps.INK)
    a.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0, -0.16 if len(order) > 6 else -0.2))
    p = pat.head(topn)[::-1]
    short = [x if len(x) < 58 else x[:55] + "..." for x in p.pattern]
    b.barh(range(len(p)), p.n, color="#eb6834", height=0.65)
    for i, (n, pc) in enumerate(zip(p.n, p.pct)):
        b.text(n, i, f"  {n:,} ({pc:g}%)", va="center", fontsize=8.5, color=ps.INK2)
    b.set_yticks(range(len(p))); b.set_yticklabels(short, fontsize=8.5)
    b.set_xlim(0, p.n.max() * 1.25)
    b.set_xlabel("unmatched entities")
    b.set_title(f"B  identifier combinations among the {n_unm:,} unmatched (top {len(p)})", loc="left", fontsize=10.5, color=ps.INK)
    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)
        ax.xaxis.grid(True, color=ps.RULE, linewidth=0.5); ax.set_axisbelow(True)
    fig.suptitle(title, x=0.01, ha="left", fontsize=12.5)
    ps.finish(fig, fname)


cov, pat, n = profile("reactions", mr)
plot("reactions", cov, pat, n, "Human-GEM reactions: identifier types of those NOT matched to human KG evidence", "fig7_xref_reactions.png")
cov, pat, n = profile("metabolites", mm)
plot("metabolites", cov, pat, n, "Human-GEM metabolites (unique compounds): identifier types of those NOT matched", "fig9_xref_metabolites.png")

summary = pd.DataFrame(summary)
summary.to_csv(R / "xref_summary.csv", index=False)
print(summary.to_string(index=False))
