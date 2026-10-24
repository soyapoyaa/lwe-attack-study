"""Measured Gram-Schmidt profiles against the GSA and CN11 shape models.

Reads results/raw/profiles.csv (from `run_sweep.py profile`). For each n and
beta (prediction - 8 and prediction) and each recorded tour t, it averages
log2 ||b_i*|| over the calls still unsolved after t tours, and compares, at the
index the success condition uses (i = d - beta, 0-based):

  measured      mean log2 ||b*_{d-beta}|| of the reduced bases
  GSA           the estimator's GSA formula (estimator/simulator.py::GSA)
  CN11 (q-ary)  fpylll's CN11 simulator after t tours, started from the q-ary
                basis, as the estimator's CN11 model does
  CN11 (LLL)    the same simulator started from the measured post-LLL profile
  threshold     log2 sqrt(sigma^2 (beta - 1) + tau^2), the expected length of
                the target's projection that the estimator compares against

A positive (value - threshold) means the condition predicts success. If the
measured margin is larger than the GSA margin, the real basis is "better than
GSA" where it matters, which would explain success below the predicted beta.

Also checks that the run reproduced prob.csv (solved within 8 tours).

Run from analysis/ in an environment with fpylll:
    python3 gs_profiles.py
Writes gs_profiles.json and figs/fig7_gs_profile.{png,pdf}.
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from gsa_estimate import delta_f  # noqa: E402  (same delta(beta) table/formula as the estimator)

RAW = os.path.join(HERE, "..", "results", "raw")
Q, SIGMA, TAU = 257, 4.0, 4.0
GSA_PRED = {32: 20, 36: 34, 40: 44, 44: 54, 48: 63}
SNAPS = (1, 2, 4, 8, 16)


def gsa_log2(d, n, beta):
    log_vol = math.log2(Q) * (d - n - 1) + math.log2(TAU)          # xi = 1
    ld = math.log2(delta_f(beta))
    return [(d - 1 - 2 * i) * ld + log_vol / d for i in range(d)]


def cn11_log2(start_log2, beta, tours):
    from fpylll import BKZ
    from fpylll.tools.bkz_simulator import simulate
    r = [2.0 ** (2 * x) for x in start_log2]
    out, _ = simulate(r, BKZ.Param(block_size=beta, max_loops=tours))
    return [0.5 * math.log2(x) for x in out]


def main():
    df = pd.read_csv(os.path.join(RAW, "profiles.csv"))
    df["P"] = df.profiles.apply(json.loads)
    df["solved_tour"] = pd.to_numeric(df.solved_tour, errors="coerce")
    res = {"replication": {}, "rows": [], "notes": []}

    # replication against the original 8-tour sweep
    P = pd.read_csv(os.path.join(RAW, "prob.csv"))
    m = df.merge(P, on=["n", "beta", "seed", "m"], suffixes=("_p", "_s"))
    m = m[m.timed_out_s == 0]
    agree = int(((m.solved_tour <= 8).astype(int) == m.solved_s).sum())
    res["replication"] = dict(pairs=int(len(m)), agree_within_8_tours=agree)
    print(f"Replication of prob.csv (solved within 8 tours): {agree}/{len(m)}")

    print("\nlog2 ||b*_{d-beta}||, minus the threshold log2 sqrt(sigma^2(beta-1)+tau^2); positive = condition met")
    print(" n   beta  tour  #unsolved  measured   GSA   CN11(q-ary)  CN11(LLL)")
    for (n, beta), g in df.groupby(["n", "beta"]):
        d = int(g.d.iloc[0]); k = d - beta
        thr = 0.5 * math.log2(SIGMA ** 2 * (beta - 1) + TAU ** 2)
        gsa = gsa_log2(d, n, beta)
        lll_mean = np.mean(np.array([p["0"] for p in g.P]), axis=0)
        qary = [math.log2(Q)] * (d - n - 1) + [0.0] * n + [math.log2(TAU)]
        for t in SNAPS:
            profs = [p[str(t)] for p in g.P if str(t) in p]
            if not profs:
                continue
            A = np.array(profs)
            meas = float(A[:, k].mean())
            cq = cn11_log2(qary, beta, t)[k]
            cl = cn11_log2(list(lll_mean), beta, t)[k]
            row = dict(n=int(n), beta=int(beta), tour=t, unsolved=len(profs), d=d, index=k,
                       threshold=thr, measured=meas, measured_sd=float(A[:, k].std()),
                       gsa=gsa[k], cn11_qary=cq, cn11_lll=cl,
                       margin_measured=meas - thr, margin_gsa=gsa[k] - thr,
                       margin_cn11_qary=cq - thr, margin_cn11_lll=cl - thr,
                       mean_profile=[round(float(x), 4) for x in A.mean(axis=0)])
            res["rows"].append(row)
            print(f"{n:3d} {beta:5d} {t:5d} {len(profs):9d}  {meas - thr:+8.3f} {gsa[k] - thr:+7.3f} "
                  f"{cq - thr:+10.3f} {cl - thr:+10.3f}")

    json.dump(res, open(os.path.join(HERE, "gs_profiles.json"), "w"), indent=1)

    # figure: n = 40 and 44 at the predicted beta, after 8 tours
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLU, VER, GRN, INK, GRY = "#0072B2", "#D55E00", "#009E73", "#222222", "#7f7f7f"
    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.alpha": .25,
                         "figure.dpi": 200, "savefig.dpi": 200, "legend.frameon": False})
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 3.4), sharey=False)
    for ax, n in zip(axs, (40, 44)):
        beta = GSA_PRED[n]
        rows = [r for r in res["rows"] if r["n"] == n and r["beta"] == beta and r["tour"] == 8]
        if not rows:
            ax.set_visible(False); continue
        r = rows[0]; d, k = r["d"], r["index"]
        g = df[(df.n == n) & (df.beta == beta)]
        A = np.array([p["8"] for p in g.P if "8" in p])
        x = np.arange(d)
        ax.fill_between(x, np.percentile(A, 10, axis=0), np.percentile(A, 90, axis=0), color=BLU, alpha=.18, lw=0)
        ax.plot(x, A.mean(axis=0), color=BLU, lw=1.4, label="measured (mean, 10–90%)")
        ax.plot(x, gsa_log2(d, n, beta), color=INK, ls="--", lw=1.1, label="GSA")
        lll_mean = np.mean(np.array([p["0"] for p in g.P]), axis=0)
        ax.plot(x, cn11_log2(list(lll_mean), beta, 8), color=VER, lw=1.1, label="CN11 from measured LLL")
        ax.axvline(k, color=GRY, ls=":", lw=1)
        ax.plot([k], [r["threshold"]], marker="_", ms=14, mew=2, color=GRN, label="threshold at d − β")
        ax.set_title(f"n = {n}, β = {beta}, after 8 tours ({len(A)} unsolved)", fontsize=9)
        ax.set_xlabel("index i"); ax.set_ylabel("log₂ ‖b*ᵢ‖")
    axs[0].legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, "figs", f"fig7_gs_profile.{ext}"))
    print("\nwrote gs_profiles.json and figs/fig7_gs_profile.png/.pdf")


if __name__ == "__main__":
    main()
