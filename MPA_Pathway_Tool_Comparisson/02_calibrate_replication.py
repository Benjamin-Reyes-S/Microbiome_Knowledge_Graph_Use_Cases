"""Calibration: reproduce the paper's Calculator results (Supp. 1, 13-15) from the paper's own pathway files,
using the same Cypher matching query (parameter-supplied reactions, no graph content), to confirm the matching rule."""
import numpy as np, pandas as pd
import common as c

df, samples = c.load_proteins()
prots = c.protein_params(df)
supp1 = c.load_supp1()
out = []
for pw, d in c.PATHWAYS.items():
    steps = c.parse_pathway_csv(c.DATA / "pathways" / d["kegg_csv"])
    reactions = [dict(rid=s["rid"], ids=["KEGG:" + k for k in s["ko"]] + ["EC:" + e for e in s["ec"]]) for s in steps]
    taxa = steps[0]["taxa"]
    res = c.cypher(c.MATCH_PARAM_Q, reactions=reactions, prots=prots, taxa=taxa)
    matches = {r["rid"]: r["idx"] for r in res}
    rx, multi = c.counts_from_matches(matches, df, samples, [r["rid"] for r in reactions])
    ref = supp1[d["supp1"]]
    out.append(dict(pathway=pw, kind="KEGG-module", variant="summed over reactions",
                    rel_abs_dev_vs_supp1=np.abs(multi - ref["counts"]).sum() / ref["counts"].sum(),
                    max_abs_diff=np.abs(multi - ref["counts"]).max(), steps_found=sum(1 for v in matches.values() if v),
                    steps_found_supp1=ref["found"], total_steps=len(steps), total_supp1=ref["total"]))
    for rid, p in c.load_supp_reactions(d["supp_rx"]).items():
        mine = rx.get(rid, np.zeros(44))
        out.append(dict(pathway=pw, kind="KEGG-reaction", variant=rid,
                        rel_abs_dev_vs_supp1=np.abs(mine - p["counts"]).sum() / max(p["counts"].sum(), 1), max_abs_diff=np.abs(mine - p["counts"]).max()))
r = pd.DataFrame(out)
r.to_csv(c.RESULTS / "calibration.csv", index=False)
pd.set_option("display.width", 250)
print(r.round(4).to_string())
