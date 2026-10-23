"""Does instance difficulty track the length of the planted vector t?

Regenerates every instance from its seed (src/lwe_instance.gen_lwe, exactly as
the runs did), computes ||t|| = sqrt(||e||^2 + ||s||^2 + tau^2), and:

  1. checks the regeneration against ladder.csv's recorded target_norm;
  2. per n, Spearman correlation between ||t|| and the fraction of block sizes
     at which the instance was solved (8-tour sweep and 32-tour run);
  3. compares the mean ||t|| of the 20 instances with the expected length the
     estimator uses, sqrt((m + n) sigma^2 + tau^2), and with the slightly
     larger expectation for a rounded Gaussian (variance sigma^2 + 1/12).

Run from analysis/:  python3 target_norm.py   -> target_norm.json
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from lwe_instance import gen_lwe  # noqa: E402

RAW = os.path.join(HERE, "..", "results", "raw")
Q, SIGMA, TAU = 257, 4.0, 4
M = {32: 71, 36: 84, 40: 88, 44: 91}


def tnorm(n, m, seed):
    i = gen_lwe(n, Q, m, seed, "dg", SIGMA, "dg", SIGMA)
    return float(math.sqrt(float((i.e ** 2).sum() + (i.s ** 2).sum()) + TAU ** 2))


def main():
    res = {"regeneration": {}, "spearman": {}, "mean_vs_expected": {}}
    L = pd.read_csv(os.path.join(RAW, "ladder.csv"))
    ok = [abs(tnorm(int(r.n), int(r.m), int(r.seed)) - r.target_norm) < 1e-2 for _, r in L.iterrows()]
    res["regeneration"] = [int(sum(ok)), len(ok)]
    print(f"Regenerated ||t|| matches ladder.csv for {sum(ok)}/{len(ok)} instances")

    runs = {"8-tour sweep": pd.read_csv(os.path.join(RAW, "prob.csv"))}
    if os.path.exists(os.path.join(RAW, "long.csv")):
        runs["32-tour run"] = pd.read_csv(os.path.join(RAW, "long.csv"))
    print("\nSpearman correlation between ||t|| and the fraction of block sizes solved (negative = shorter t, easier)")
    for name, D in runs.items():
        D = D[D.timed_out == 0]
        res["spearman"][name] = {}
        for n, m in M.items():
            rate = D[D.n == n].groupby("seed").solved.mean()
            t = pd.Series({s: tnorm(n, m, int(s)) for s in rate.index})
            rho, p = spearmanr(t, rate)
            res["spearman"][name][n] = dict(rho=float(rho), p=float(p), t_min=t.min(), t_max=t.max())
            print(f"  {name:12s} n={n}: rho = {rho:+.2f} (p = {p:.3f}), ||t|| from {t.min():.1f} to {t.max():.1f}")

    print("\nMean ||t|| of the 20 instances against the expected length")
    for n, m in M.items():
        t = np.array([tnorm(n, m, s) for s in range(1, 21)])
        exp = math.sqrt((m + n) * SIGMA ** 2 + TAU ** 2)
        exp_r = math.sqrt((m + n) * (SIGMA ** 2 + 1 / 12) + TAU ** 2)
        res["mean_vs_expected"][n] = dict(mean=float(t.mean()), expected=exp, expected_rounded=exp_r,
                                          ratio=float(t.mean() / exp), below_expected=int((t < exp).sum()))
        print(f"  n={n}: mean {t.mean():.2f}, estimator's expectation {exp:.2f} (ratio {t.mean() / exp:.3f}), "
              f"rounded-Gaussian expectation {exp_r:.2f}, {int((t < exp).sum())}/20 below")
    json.dump(res, open(os.path.join(HERE, "target_norm.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
