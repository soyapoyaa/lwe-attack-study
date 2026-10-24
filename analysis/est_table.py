"""Table 8 / Section 5.7: compare primal_usvp's returned point with the best beta>=40 point and the global minimum
of cost_gsa, using the exhaustive scans in results/scan_n*.csv (made by src/extend_scan_beta.sage with the real estimator).
RET = beta returned by primal_usvp, as recorded in the study protocol; replay_local_minimum.py reproduces these values."""
import os, json, pandas as pd
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
RET = {32: 65, 36: 73, 40: 81, 44: 54, 48: 63}
out = {}
for n in (32, 36, 40, 44, 48):
    s = pd.read_csv(os.path.join(RES, f"scan_n{n}.csv")).dropna(subset=["log2rop"]).set_index("beta")
    g = s.log2rop.idxmin(); b40 = s.loc[40:].log2rop.idxmin(); r = RET[n]
    out[str(n)] = dict(g=[int(g), int(s.loc[g].d), float(s.loc[g].log2rop)], b40=[int(b40), int(s.loc[b40].d), float(s.loc[b40].log2rop)],
        ret=[r, int(s.loc[r].d), float(s.loc[r].log2rop)], floor_bits=float(s.loc[b40].log2rop - s.loc[g].log2rop),
        search_bits=float(s.loc[r].log2rop - s.loc[b40].log2rop), clamp_start=int(s[s.d == s.index + 1].index.min()))
    o = out[str(n)]; print(n, "global", o["g"], "| best>=40", o["b40"], "| returned", o["ret"], "| floor %.2f | search %.2f | clamp starts at beta=%d" % (o["floor_bits"], o["search_bits"], o["clamp_start"]))
json.dump(out, open("est_table.json", "w"), indent=1)
