"""Tour budget: how the answer depends on how many BKZ tours a call may run.

Replaces tour4_stats.py and extends it to every budget from 1 to 8 tours.

Why the derived curves are exact, not an approximation
------------------------------------------------------
src/attack.py::solve_single runs `bkz.tour(par)` once per loop iteration and
checks the basis for the secret after every tour. Nothing in a tour depends on
the total budget, and fpylll's AUTO_ABORT flag is only read by
BKZReduction.__call__, which this loop never calls (confirmed in the data: all
552 uncensored failures ran all 8 tours). A call with budget T therefore does
exactly what the first T tours of an 8-tour call do. Combined with the fact that
re-running a call reproduces its outcome (replay_ladder.py, check 1), "solved
within T tours" in the 8-tour data *is* the outcome of a T-tour call.

Outputs (tour_budget.json):
  - solved by tour t at the predicted beta, per n (cumulative)
  - pooled success at the prediction for each budget, with Wilson interval
  - beta_50 and beta_90 (logistic) with bootstrap intervals over instances,
    and the three other beta_50 estimators, for budgets 2, 4, 6, 8
  - downward steps (solved at one beta, unsolved at the next) per budget

Run from analysis/:  python3 tour_budget.py   (about 5-10 minutes)
"""
import json
import math
import os

import numpy as np
import pandas as pd

from sweep_stats import ests, p as P, wilson

HERE = os.path.dirname(os.path.abspath(__file__))
GSA = {32: 20, 36: 34, 40: 44, 44: 54}
LN9 = math.log(9)
B = 1000
PC = P[P.timed_out == 0].copy()


def with_budget(df, T):
    d = df.copy()
    d.loc[d.tours_used > T, "solved"] = 0
    return d


def bootstrap(g, rng):
    seeds = g.seed.unique()
    by_seed = {s: g[g.seed == s] for s in seeds}
    e = ests(g)
    e["b90"] = e["logit"] + e["k"] * LN9
    draws = {k: [] for k in e}
    for _ in range(B):
        gg = pd.concat([by_seed[s] for s in rng.choice(seeds, len(seeds))])
        eb = ests(gg)
        eb["b90"] = eb["logit"] + eb["k"] * LN9
        for k in draws:
            draws[k].append(eb[k])
    return {k: (float(e[k]), [float(x) for x in np.nanpercentile(draws[k], [2.5, 97.5])])
            for k in e}


def down_steps(df):
    out = 0
    for _, g in df.groupby(["n", "seed"]):
        s = g.sort_values("beta").solved.values
        out += int(np.sum((s[:-1] == 1) & (s[1:] == 0)))
    return out


def main():
    res = {"cumulative_at_pred": {}, "pooled_at_pred": {}, "estimates": {}, "down_steps": {}}

    print("Solved at the predicted beta, by tour (cumulative, of 20):")
    for n, b in GSA.items():
        h = PC[(PC.n == n) & (PC.beta == b)]
        cum = [int(((h.solved == 1) & (h.tours_used <= t)).sum()) for t in range(1, 9)]
        res["cumulative_at_pred"][n] = cum
        print(f"  n={n}: {cum}")

    print("\nPooled at the prediction (80 instances):")
    for T in range(1, 9):
        d = with_budget(PC, T)
        k = int(sum(d[(d.n == n) & (d.beta == b)].solved.sum() for n, b in GSA.items()))
        lo, hi = wilson(k, 80)
        res["pooled_at_pred"][T] = [k, 80, float(lo), float(hi)]
        res["down_steps"][T] = down_steps(d)
        print(f"  {T} tours: {k}/80 = {k/80:.3f}  Wilson [{lo:.3f}, {hi:.3f}]   "
              f"downward steps {res['down_steps'][T]}")

    rng = np.random.default_rng(2026)
    for T in (2, 4, 6, 8):
        d = with_budget(PC, T)
        res["estimates"][T] = {}
        print(f"\n{T} tours:")
        for n, b in GSA.items():
            r = bootstrap(d[d.n == n], rng)
            res["estimates"][T][n] = r
            print(f"  n={n} pred {b}: b50 logit {r['logit'][0]:.1f} [{r['logit'][1][0]:.1f}, "
                  f"{r['logit'][1][1]:.1f}] (x{r['logit'][0]/b:.2f}) | b90 {r['b90'][0]:.1f} "
                  f"[{r['b90'][1][0]:.1f}, {r['b90'][1][1]:.1f}] | first/last/iso "
                  f"{r['first'][0]:.1f}/{r['last'][0]:.1f}/{r['iso'][0]:.1f}")
    json.dump(res, open(os.path.join(HERE, "tour_budget.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
