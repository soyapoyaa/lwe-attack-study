"""Analyse the long-tour control run (results/raw/long.csv, from `run_sweep.py long`).

Questions it answers, in order:
  1. Does the run reproduce the original sweep? For every (n, beta, seed) in
     both files, "solved within 8 tours" in long.csv must equal "solved" in
     prob.csv, and the tour of success must match. If this fails, the derived
     budgets below are not exact and must not be reported as such.
  2. Where does the success curve settle as the budget grows? beta_50, beta_90
     (logistic, bootstrap over instances) and success at the prediction for
     budgets 4, 8, 12, 16, 24, 32, and for a call that stops the way a default
     fpylll BKZ call would (first auto-abort or clean tour).
  3. How many tours does a default fpylll call actually run here?
  4. Does a larger budget remove the non-monotone dips (Section 5.4)?

Run from analysis/:  python3 tour_convergence.py  -> tour_convergence.json
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd

from sweep_stats import ests, wilson

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "results", "raw")
GSA = {32: 20, 36: 34, 40: 44, 44: 54}
BUDGETS = (4, 8, 12, 16, 24, 32)
B = 500
LN9 = math.log(9)


def load():
    L = pd.read_csv(os.path.join(RAW, "long.csv"))
    P = pd.read_csv(os.path.join(RAW, "prob.csv"))
    for c in ("solved_tour", "abort_tour", "clean_tour"):
        L[c] = pd.to_numeric(L[c], errors="coerce")
    return L, P


def outcome(L, policy):
    """0/1 outcome of each row under a stopping policy; NaN if unknown (timed out)."""
    st = L.solved_tour
    if policy == "fpylll_stop":
        stop = L[["abort_tour", "clean_tour"]].min(axis=1)       # NaN = never stopped
        y = (st.notna() & (stop.isna() | (st <= stop))).astype(float)
        # unsolved and timed out before any stop: unknown
        unknown = (L.timed_out == 1) & st.isna() & stop.isna()
    else:
        T = int(policy)
        y = (st.notna() & (st <= T)).astype(float)
        unknown = (L.timed_out == 1) & st.isna() & (L.tours_run < T)
    y[unknown] = np.nan
    return y


def fit_boot(g, rng):
    g = g.dropna(subset=["solved"])
    e = ests(g)
    e["b90"] = e["logit"] + e["k"] * LN9
    seeds = g.seed.unique()
    by = {s: g[g.seed == s] for s in seeds}
    dr = {k: [] for k in ("logit", "b90")}
    for _ in range(B):
        gg = pd.concat([by[s] for s in rng.choice(seeds, len(seeds))])
        eb = ests(gg)
        dr["logit"].append(eb["logit"])
        dr["b90"].append(eb["logit"] + eb["k"] * LN9)
    ci = {k: [float(x) for x in np.nanpercentile(v, [2.5, 97.5])] for k, v in dr.items()}
    return dict(b50=float(e["logit"]), b50_ci=ci["logit"], b90=float(e["b90"]), b90_ci=ci["b90"],
                first=float(e["first"]), last=float(e["last"]), iso=float(e["iso"]))


def main():
    L, P = load()
    ns = sorted(n for n in L.n.unique() if n in GSA)
    res = {"n": [int(n) for n in ns], "replication": {}, "budgets": {}, "stop_tours": {},
           "down_steps": {}, "cumulative_at_pred": {}}

    # 1. replication of the original 8-tour sweep
    m = L.merge(P, on=["n", "beta", "seed", "m"], suffixes=("_long", "_prob"))
    m = m[(m.timed_out_prob == 0) & ~((m.timed_out_long == 1) & m.solved_tour.isna() & (m.tours_run < 8))]
    s8 = m.solved_tour.notna() & (m.solved_tour <= 8)
    agree = int((s8.astype(int) == m.solved_prob).sum())
    same_tour = int(((m.solved_prob == 1) & (m.solved_tour == m.tours_used)).sum())
    res["replication"] = dict(pairs=int(len(m)), agree_8tour=agree,
                              solved_in_prob=int(m.solved_prob.sum()), same_tour=same_tour)
    print(f"1. Replication of prob.csv: {agree}/{len(m)} outcomes agree within 8 tours; "
          f"{same_tour}/{int(m.solved_prob.sum())} successes on the same tour")
    if agree != len(m):
        print("   WARNING: outcomes are not reproducible. Report budgets as measured on this run,"
              " not as exact extensions of prob.csv.", file=sys.stderr)

    # 2. budgets
    rng = np.random.default_rng(7)
    print("\n2. Success curve by budget (logistic, bootstrap over instances)")
    for pol in [str(T) for T in BUDGETS] + ["fpylll_stop"]:
        d = L.copy()
        d["solved"] = outcome(L, pol)
        res["budgets"][pol] = {}
        kk = NN = 0
        for n in ns:
            g = d[d.n == n]
            r = fit_boot(g, rng)
            at = g[(g.beta == GSA[n])].solved.dropna()
            r["at_pred"] = [int(at.sum()), int(len(at))]
            kk += int(at.sum()); NN += int(len(at))
            res["budgets"][pol][int(n)] = r
        lo, hi = wilson(kk, NN) if NN else (float("nan"), float("nan"))
        res["budgets"][pol]["pooled_at_pred"] = [kk, NN, float(lo), float(hi)]
        res["down_steps"][pol] = int(sum(
            np.sum((s[:-1] == 1) & (s[1:] == 0))
            for _, gg in d.dropna(subset=["solved"]).groupby(["n", "seed"])
            for s in [gg.sort_values("beta").solved.values]))
        label = f"{pol} tours" if pol != "fpylll_stop" else "default stop"
        print(f"  {label:13s} pooled at pred {kk}/{NN} | downward steps {res['down_steps'][pol]} | "
              + " | ".join(f"n={n}: b50 {res['budgets'][pol][int(n)]['b50']:.1f} "
                           f"(x{res['budgets'][pol][int(n)]['b50']/GSA[n]:.2f})" for n in ns))

    # 3. where a default fpylll call would stop
    print("\n3. Tour at which a default fpylll call would stop (unsolved rows)")
    for n in ns:
        g = L[(L.n == n) & L.solved_tour.isna()]
        stop = g[["abort_tour", "clean_tour"]].min(axis=1)
        q = stop.describe(percentiles=[.1, .5, .9])
        res["stop_tours"][int(n)] = dict(count=int(stop.notna().sum()), never=int(stop.isna().sum()),
                                         p10=float(q.get("10%", np.nan)), median=float(q.get("50%", np.nan)),
                                         p90=float(q.get("90%", np.nan)))
        print(f"  n={n}: stopped {stop.notna().sum()} / never within budget {stop.isna().sum()} | "
              f"median stop tour {res['stop_tours'][int(n)]['median']}")

    # cumulative solved at the prediction, by tour
    print("\n4. Solved at the predicted beta, by tour (cumulative)")
    for n in ns:
        h = L[(L.n == n) & (L.beta == GSA[n])]
        cum = [int((h.solved_tour <= t).sum()) for t in range(1, int(L.tours_max.max()) + 1)]
        res["cumulative_at_pred"][int(n)] = cum
        print(f"  n={n} (of {len(h)}): " + " ".join(str(c) for c in cum))

    json.dump(res, open(os.path.join(HERE, "tour_convergence.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
