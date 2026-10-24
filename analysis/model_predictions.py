"""Model predictions for our exact instances, to compare with the measurements.

Two models, both taking the number of BKZ tours into account:

  (A) Postlethwaite-Virdia [PKC 2021] BKZ success-probability simulator.
      Their code (github.com/fvirdia/usvp-simulation) has no licence file, so it
      is NOT copied here: pass the path of a local clone (commit 333416d) with --pv21. We import
      their CN11 BKZ simulator, LLL profile simulator and success accumulator
      unchanged, and replace only their Sage chi-squared CDF with scipy's
      (the same function). Embedding: their "baigal" with nu = 1, which is our
      lattice (q*I_m, I_n, tau) with d = m + n + 1.
      Output: P(success within tau tours) for every beta on our grids, for
      tau in {4, 8, 12, 16, 24, 32}.

  (B) The estimator's own success condition with the CN11 shape model, at the
      lattice dimension we attack: smallest beta with
          r[d - beta] > sigma^2 (beta - 1) + tau^2,
      r simulated by fpylll's CN11 simulator from the q-ary basis (exactly as
      estimator/simulator.py does), after T tours, T in {4, 8, 16, 32} and
      "converged" (the estimator's default: run until no progress).
      The GSA version of the same condition gives our beta_GSA (20, 34, 44, 54, 63).

Then compares (A) with the measured success curves (prob.csv for 4 and 8 tours,
long.csv for 8 to 32 tours): 50% points and the curves themselves.

Run from analysis/ in an environment with fpylll and scipy:
    python3 model_predictions.py --pv21 ../../usvp-simulation
Writes model_predictions.json and figs/fig6_model_vs_measured.{png,pdf}.
"""
import argparse
import contextlib
import io
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "results", "raw")
Q, SIGMA, TAU = 257, 4.0, 4.0
NS = {32: (71, 20), 36: (84, 34), 40: (88, 44), 44: (91, 54), 48: (96, 63)}   # n: (m used, beta_GSA)
TOURS = (4, 8, 12, 16, 24, 32)


# ---------------------------------------------------------------- (A) PV21
def pv21_curve(pv21, n, m, betas, tours):
    """P(success within `tours` tours of BKZ-beta), PV21's simulate_bkz with CN11, embedding baigal, nu=1."""
    CN11, lll_simulator, SuccessAccumulator, BKZ = pv21
    d = n + m + 1
    with contextlib.redirect_stdout(io.StringIO()):   # their LLL simulator prints a diagnostic in small dimensions
        lll = lll_simulator(d, m, Q, scale=1)       # = their LLLProfile(n, q, m, nu=1, embedding="baigal")
    out = {}
    for beta in betas:
        sa = SuccessAccumulator()
        for t in range(1, tours + 1):
            prof, _ = CN11.simulate(lll[:], BKZ.Param(block_size=beta, max_loops=t))
            p = float(chi2.cdf(prof[d - beta] / SIGMA ** 2, beta))   # their cmf_chi_sq(profile[d-beta]/sd^2, beta)
            sa.accumulate(min(max(p, 0.0), 1.0))
        out[int(beta)] = sa.success_probability()
    return out


def load_pv21(path):
    sys.path.insert(0, os.path.abspath(path))
    import bkz_simulators.CN11 as CN11                 # their copy of the CN11 simulator
    from lll_sim import lll_simulator                  # LLL "Z-shape" profile
    from success_accumulator import SuccessAccumulator
    from fpylll import BKZ
    return CN11, lll_simulator, SuccessAccumulator, BKZ


# ---------------------------------------------------------------- (B) estimator condition with CN11
def cn11_threshold(n, m, tours):
    """Smallest beta passing the estimator's uSVP condition at d = m+n+1 under CN11 after `tours` (None = converge)."""
    from fpylll import BKZ
    from fpylll.tools.bkz_simulator import simulate
    d = m + n + 1
    r0 = [Q ** 2] * m + [1.0] * n + [TAU ** 2]          # estimator.simulator.qary_simulator, xi = 1
    for beta in range(3, d):
        par = BKZ.Param(block_size=beta, max_loops=0 if tours is None else tours)
        r, _ = simulate(r0[:], par)
        if r[d - beta] > SIGMA ** 2 * (beta - 1) + TAU ** 2:
            return beta
    return None


# ---------------------------------------------------------------- measured
def measured_curves():
    P = pd.read_csv(os.path.join(RAW, "prob.csv")); P = P[P.timed_out == 0]
    out = {}
    if os.path.exists(os.path.join(RAW, "long.csv")):
        L = pd.read_csv(os.path.join(RAW, "long.csv"))
        L["solved_tour"] = pd.to_numeric(L.solved_tour, errors="coerce")
        for T in TOURS:
            s = L.assign(ok=(L.solved_tour <= T).astype(float))
            s.loc[(L.timed_out == 1) & L.solved_tour.isna() & (L.tours_run < T), "ok"] = np.nan
            for n, g in s.dropna(subset=["ok"]).groupby("n"):
                out[(int(n), T)] = g.groupby("beta").ok.agg(["sum", "count"])
    for T in (4, 8):   # from the original sweep (its grid), if long.csv lacks them
        s = P.assign(ok=((P.solved == 1) & (P.tours_used <= T)).astype(float))
        for n, g in s.groupby("n"):
            out.setdefault((int(n), T), g.groupby("beta").ok.agg(["sum", "count"]))
    return out


def cross50(curve):
    """Linear interpolation of the first 0.5 crossing of a (beta -> p) dict."""
    b = sorted(curve)
    for b0, b1 in zip(b, b[1:]):
        if curve[b0] < 0.5 <= curve[b1]:
            return b0 + (0.5 - curve[b0]) * (b1 - b0) / (curve[b1] - curve[b0])
    return float("nan")


def logistic50(cnt):
    from scipy.optimize import minimize
    x = np.repeat(cnt.index.values.astype(float), cnt["count"].astype(int))
    y = np.concatenate([[1] * int(k) + [0] * int(c - k) for k, c in zip(cnt["sum"], cnt["count"])])
    def nll(t):
        p = np.clip(1 / (1 + np.exp(-(x - t[0]) / np.exp(t[1]))), 1e-12, 1 - 1e-12)
        return -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
    return float(minimize(nll, [x.mean(), 1.0], method="Nelder-Mead").x[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pv21", required=True, help="path to a clone of github.com/fvirdia/usvp-simulation")
    ap.add_argument("--no-fig", action="store_true")
    args = ap.parse_args()
    pv21 = load_pv21(args.pv21)
    meas = measured_curves()
    res = {"pv21": {}, "cn11_threshold": {}, "comparison": {}}

    print("(B) estimator success condition at the attacked d (smallest beta that passes)")
    for n, (m, bg) in NS.items():
        row = {"gsa": bg}
        for T in (4, 8, 16, 32, None):
            row["converged" if T is None else str(T)] = cn11_threshold(n, m, T)
        res["cn11_threshold"][n] = row
        print(f"  n={n}: GSA {bg} | CN11 after 4/8/16/32 tours / converged: "
              + " / ".join(str(row[k]) for k in ("4", "8", "16", "32", "converged")))

    print("\n(A) PV21 BKZ simulator: predicted 50% point vs measured (logistic), by tour budget")
    for n, (m, bg) in NS.items():
        betas = list(range(max(4, bg - 24), bg + 17, 1))
        res["pv21"][n] = {}
        for T in TOURS:
            curve = pv21_curve(pv21, n, m, betas, T)
            res["pv21"][n][T] = curve
            p50 = cross50(curve)
            mc = meas.get((n, T))
            m50 = logistic50(mc) if mc is not None else float("nan")
            # mean absolute difference between predicted and measured success fraction on the measured grid
            mad = (float(np.mean([abs(curve.get(int(b), 1.0) - k / c) for b, (k, c) in mc.iterrows()]))
                   if mc is not None else float("nan"))
            # The simulator adds up per-tour success chances as if tours were independent; at very small
            # beta and many tours this saturates towards 1 even where no instance is ever solved. Flag it.
            low = min(curve)
            reliable = curve[low] < 0.10
            res["comparison"][f"{n}_{T}"] = dict(n=n, tours=T, pv21_b50=p50 if reliable else None, measured_b50=m50,
                                                 mean_abs_diff=mad, pv21_p_at_pred=curve.get(bg),
                                                 pv21_p_at_lowest_beta=[low, curve[low]], reliable=reliable)
            flag = "" if reliable else f"  <- UNRELIABLE: simulator gives P={curve[low]:.2f} already at beta={low}"
            print(f"  n={n} tours={T:2d}: PV21 50% at {p50 if reliable else float('nan'):5.1f} | measured {m50:5.1f} | "
                  f"PV21 P(at prediction {bg}) = {curve.get(bg, float('nan')):.2f} | mean |pred - meas| = {mad:.2f}{flag}")

    json.dump(res, open(os.path.join(HERE, "model_predictions.json"), "w"), indent=1,
              default=lambda o: None if (isinstance(o, float) and math.isnan(o)) else o)

    if not args.no_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        BLU, VER, INK = "#0072B2", "#D55E00", "#222222"
        plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9, "axes.spines.top": False,
                             "axes.spines.right": False, "axes.grid": True, "grid.alpha": .25,
                             "figure.dpi": 200, "savefig.dpi": 200, "legend.frameon": False})
        fig, axs = plt.subplots(2, 2, figsize=(6.6, 5.2), sharey=True)
        for ax, n in zip(axs.flat, (32, 36, 40, 44)):
            m, bg = NS[n]
            for T, col, mk in ((8, BLU, "o"), (32, VER, "s")):
                c = res["pv21"][n][T]; b = sorted(c)
                ax.plot(b, [c[x] for x in b], color=col, lw=1.4, label=f"PV21 simulator, {T} tours")
                mc = meas.get((n, T))
                if mc is not None:
                    ax.plot(mc.index, mc["sum"] / mc["count"], mk, color=col, ms=4, mec="white", mew=.5,
                            label=f"measured, {T} tours")
            ax.axvline(bg, color=INK, ls="--", lw=1)
            ax.set_title(f"n = {n}  (prediction β = {bg})", fontsize=9)
            ax.set_xlabel("block size β")
        axs[0, 0].set_ylabel("fraction solved"); axs[1, 0].set_ylabel("fraction solved")
        h, l = axs[0, 0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=2, fontsize=7.4, bbox_to_anchor=(0.5, -0.005))
        fig.tight_layout(rect=[0, 0.07, 1, 1])
        os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(HERE, "figs", f"fig6_model_vs_measured.{ext}"))
        print("\nwrote model_predictions.json and figs/fig6_model_vs_measured.png/.pdf")


if __name__ == "__main__":
    main()
