"""Compare module-level and reaction-level spectral counts: paper's Pathway-Calculator (Supp. 1, 13-15) vs. knowledge graph."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import common as c

# dataviz palette slots 1-3 (light surface), text/ink tokens
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "text.color": INK,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": GRID,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left"})

_, samples = c.load_proteins()
supp1 = c.load_supp1()
kgm = pd.read_csv(c.RESULTS / "kg_module_counts.csv")
kgr = pd.read_csv(c.RESULTS / "kg_reaction_counts.csv")
V_ALL, V_SHARED, V_UNIQ = (kgm.variant.unique()[i] for i in (0, 1, 2))


def stats(a, b):
    return dict(pearson_r=pearsonr(a, b)[0] if a.std() and b.std() else np.nan, spearman=spearmanr(a, b)[0] if a.std() and b.std() else np.nan,
                mean_rel_dev=np.mean(np.abs(b - a) / np.maximum(a, 1)), total_paper=a.sum(), total_kg=b.sum(), ratio_kg_to_paper=b.sum() / a.sum() if a.sum() else np.nan)

# ---------------- module level ----------------
mod_rows, M = [], {}
for pw, d in c.PATHWAYS.items():
    paper = supp1[d["supp1"]]
    M[pw] = dict(paper=paper["counts"])
    for tag, v in (("kg", V_ALL), ("kg_shared", V_SHARED), ("kg_unique", V_UNIQ)):
        r = kgm[(kgm.pathway == pw) & (kgm.variant == v)].iloc[0]
        M[pw][tag] = r[samples].to_numpy(float)
        mod_rows.append(dict(pathway=pw, module=d["module"], comparison=f"paper vs {v}", steps_found_paper=paper["found"], steps_total_paper=paper["total"],
                             steps_found_kg=r.steps_found, steps_total_kg=r.steps_total, **stats(paper["counts"], M[pw][tag])))
mod = pd.DataFrame(mod_rows)
mod.to_csv(c.RESULTS / "module_comparison.csv", index=False)
pd.DataFrame([dict(pathway=pw, sample=s, paper_calculator=M[pw]["paper"][i], kg=M[pw]["kg"][i], kg_shared_reactions=M[pw]["kg_shared"][i], kg_unique_proteins=M[pw]["kg_unique"][i])
              for pw in M for i, s in enumerate(samples)]).to_csv(c.RESULTS / "module_comparison_per_sample.csv", index=False)

# ---------------- reaction level ----------------
rx_rows, R = [], {}
defs = pd.read_csv(c.RESULTS / "kg_module_definitions.csv").set_index(["pathway", "rid"])
diff = pd.read_csv(c.RESULTS / "definition_diff_kg_vs_paper.csv").fillna("").set_index(["pathway", "rid"])
for pw, d in c.PATHWAYS.items():
    paper = c.load_supp_reactions(d["supp_rx"])
    for rid in sorted(set(paper) | set(kgr[kgr.pathway == pw].rid)):
        a = paper[rid]["counts"] if rid in paper else np.zeros(44)
        b = kgr[(kgr.pathway == pw) & (kgr.rid == rid)][samples].to_numpy(float)
        b = b[0] if len(b) else np.zeros(44)
        R[(pw, rid)] = (a, b)
        dd = diff.loc[(pw, rid)] if (pw, rid) in diff.index else None
        rx_rows.append(dict(pathway=pw, rid=rid, status="both" if rid in paper and (pw, rid) in defs.index else ("paper only" if rid in paper else "KG only"),
                            name=(defs.loc[(pw, rid)]["name"] if (pw, rid) in defs.index else paper[rid]["name"].split("\t")[0])[:70],
                            **stats(a, b), identical=bool(np.allclose(a, b)),
                            kg_added_ko=dd["ko_only_kg"] if dd is not None else "", kg_added_ec=dd["ec_only_kg"] if dd is not None else "",
                            paper_only_ec=dd["ec_only_paper"] if dd is not None else ""))
rx = pd.DataFrame(rx_rows)
rx.to_csv(c.RESULTS / "reaction_comparison.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 40)
print("MODULE LEVEL\n", mod.drop(columns=["module"]).round(3).to_string())
print("\nREACTION LEVEL (non-identical or KG/paper-only)\n", rx[~rx.identical | (rx.status != "both")].drop(columns=["name"]).round(3).to_string())
print("\nidentical reactions:", rx.identical.sum(), "of", len(rx))

# ---------------- figures ----------------
PW = list(c.PATHWAYS)
x = np.arange(44)
short = [s.replace("BGP_", "") for s in samples]

fig, axs = plt.subplots(3, 1, figsize=(12, 9.5), sharex=True)
for ax, pw in zip(axs, PW):
    ax.plot(x, M[pw]["paper"], color=BLUE, lw=2, marker="o", ms=4, label="Paper: Pathway-Calculator (KEGG module, 2020 KEGG)")
    ax.plot(x, M[pw]["kg"], color=ORANGE, lw=2, marker="s", ms=4, label="Knowledge graph (current KEGG module)")
    ax.plot(x, M[pw]["kg_unique"], color=AQUA, lw=1.6, ls="--", label="Knowledge graph, each protein counted once")
    ax.set_title(f"{pw}  ({c.PATHWAYS[pw]['module'].split(':')[1]})"); ax.set_ylabel("summed spectral counts"); ax.set_ylim(bottom=0)
    ax.grid(axis="x", visible=False)
axs[0].legend(frameon=False, loc="upper right", ncol=1, fontsize=8)
axs[-1].set_xticks(x); axs[-1].set_xticklabels(short, rotation=90, fontsize=7); axs[-1].set_xlabel("biogas plant sample")
fig.suptitle("Module level: summed spectral counts per sample, paper vs. knowledge graph", x=0.01, ha="left", fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(c.FIGURES / "fig1_module_per_sample.png", dpi=160); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(13, 4.4))
for ax, pw in zip(axs, PW):
    a, b = M[pw]["paper"], M[pw]["kg"]
    m = max(a.max(), b.max()) * 1.05
    ax.plot([0, m], [0, m], color=INK2, lw=1, ls=":")
    ax.scatter(a, b, s=28, color=ORANGE, edgecolor=SURF, linewidth=0.8, zorder=3)
    r = mod[(mod.pathway == pw) & (mod.comparison.str.contains("all reactions \\(summed"))].iloc[0]
    ax.set_title(pw, fontsize=10); ax.set_xlabel("paper Calculator (spectral counts)"); ax.set_ylabel("knowledge graph (spectral counts)")
    ax.text(0.04, 0.95, f"r = {r.pearson_r:.3f}\nKG/paper total = {r.ratio_kg_to_paper:.2f}", transform=ax.transAxes, va="top", color=INK, fontsize=9)
    ax.set_xlim(0, m); ax.set_ylim(0, m); ax.set_aspect("equal")
fig.suptitle("Module level: one point per sample (dotted line = identity)", x=0.01, ha="left", fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(c.FIGURES / "fig2_module_scatter.png", dpi=160); plt.close(fig)

fig, axs = plt.subplots(3, 1, figsize=(11, 12), gridspec_kw=dict(height_ratios=[10, 13, 11]))
for ax, pw in zip(axs, PW):
    sub = rx[rx.pathway == pw].sort_values("total_paper")
    y = np.arange(len(sub)); h = 0.36
    ax.barh(y + h / 2, sub.total_paper, h - 0.04, color=BLUE, label="Paper: Pathway-Calculator")
    ax.barh(y - h / 2, sub.total_kg, h - 0.04, color=ORANGE, label="Knowledge graph")
    ax.set_yticks(y); ax.set_yticklabels([f"{r.rid.split(':')[1]}  {r['name'][:38]}" + ("  [KG only]" if r.status == "KG only" else "") for _, r in sub.iterrows()], fontsize=8)
    ax.set_title(pw); ax.set_xlabel("spectral counts summed over 44 samples"); ax.grid(axis="y", visible=False)
axs[0].legend(frameon=False, loc="lower right", fontsize=8)
fig.suptitle("Reaction level: total spectral counts per KEGG reaction", x=0.01, ha="left", fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(c.FIGURES / "fig3_reaction_totals.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.2, 6))
cols = dict(zip(PW, (BLUE, ORANGE, AQUA)))
for pw in PW:
    pts = np.array([(R[k][0][i], R[k][1][i]) for k in R if k[0] == pw for i in range(44)])
    ax.scatter(pts[:, 0], pts[:, 1], s=14, color=cols[pw], edgecolor=SURF, linewidth=0.4, alpha=0.85, label=pw)
mx = max(max(v[0].max(), v[1].max()) for v in R.values()) * 1.1
ax.plot([0, mx], [0, mx], color=INK2, lw=1, ls=":"); ax.set_xscale("symlog", linthresh=10); ax.set_yscale("symlog", linthresh=10)
ax.set_xlim(0, mx); ax.set_ylim(0, mx)
allp = np.array([(R[k][0][i], R[k][1][i]) for k in R for i in range(44)])
ax.set_xlabel("paper Calculator (spectral counts, per reaction and sample)"); ax.set_ylabel("knowledge graph (spectral counts)")
ax.set_title(f"Reaction level: r = {pearsonr(allp[:,0], allp[:,1])[0]:.3f}, {(allp[:,0]==allp[:,1]).mean()*100:.0f}% of points identical", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig(c.FIGURES / "fig4_reaction_scatter.png", dpi=160); plt.close(fig)
