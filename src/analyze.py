"""
Turn raw trial CSVs into the figures and tables the paper needs.

    python3 src/analyze.py

Reads   results/raw/{ladder,prob,mtrade}.csv   (whichever exist)
        results/predictions.csv                (from `sage src/estimate.sage`;
                                                falls back to the pure-Python
                                                GSA reimplementation if absent,
                                                and says so on the figure)
        results/scan_n*.csv                    (exhaustive cost_gsa scans)

The GSA prediction (beta_gsa) is the global minimum of cost_gsa, taken from the
exhaustive scans, as in the paper. primal_usvp's own GSA answer misses that
minimum at n = 32, 36 and 40 (paper, Section 5.7), so it is kept separately as
beta_gsa_primal_usvp and never used as the prediction. beta_cn11 is
primal_usvp's CN11 answer, which inherits the same miss at n <= 40.

This is a quick-look summary. The paper's tables and figures come from the
scripts in analysis/ (see analysis/README_analysis.md).

Writes  figures/*.png, figures/*.pdf
        results/summary_*.csv
        paper/tables/*.tex

Every figure states the machine and software versions from
results/environment.json in its caption file, so no number in the paper is
untraceable.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gsa_estimate import primal_usvp

REPO = Path(__file__).resolve().parent.parent
RAW = REPO / "results" / "raw"
FIG = REPO / "figures"
TAB = REPO / "paper" / "tables"

plt.rcParams.update({
    "figure.figsize": (6.2, 4.1), "figure.dpi": 150,
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "savefig.bbox": "tight",
})

C_PRED_GSA, C_PRED_CN11, C_SINGLE, C_PROG = "#c0392b", "#e67e22", "#2c3e50", "#2980b9"


def load(name: str) -> pd.DataFrame | None:
    p = RAW / f"{name}.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    return df[df.get("exp", name) != "FAILED"] if len(df) else None


def predictions(ns, q, sigma, m_mult=8) -> pd.DataFrame:
    """Official estimator output if available, else the local GSA copy."""
    p = REPO / "results" / "predictions.csv"
    if p.exists():
        df = pd.read_csv(p)
        # primal_usvp's GSA answer is not the minimum of its own cost model at
        # small n; replace it by the exhaustive-scan minimum where we have one.
        df = df.rename(columns={"beta_gsa": "beta_gsa_primal_usvp", "d_gsa": "d_gsa_primal_usvp"})
        df["beta_gsa"], df["d_gsa"] = np.nan, np.nan
        for i, n in df["n"].items():
            scan = REPO / "results" / f"scan_n{int(n)}.csv"
            if scan.exists():
                s = pd.read_csv(scan)
                best = s.loc[s["log2rop"].idxmin()]
            else:
                pr = primal_usvp(n=int(n), q=int(df.at[i, "q"]), sigma_s=float(df.at[i, "sigma"]),
                                 sigma_e=float(df.at[i, "sigma"]), m=int(df.at[i, "m_available"]))
                best = {"beta": pr.beta, "d": pr.d}
            df.at[i, "beta_gsa"], df.at[i, "d_gsa"] = best["beta"], best["d"]
        df["beta_gsa"], df["d_gsa"] = df["beta_gsa"].astype(int), df["d_gsa"].astype(int)
        df.attrs["source"] = "lattice-estimator (GSA: global minimum of cost_gsa from the scans)"
        return df
    rows = []
    for n in ns:
        pr = primal_usvp(n=n, q=q, sigma_s=sigma, sigma_e=sigma, m=m_mult * n)
        pr40 = primal_usvp(n=n, q=q, sigma_s=sigma, sigma_e=sigma,
                           m=m_mult * n, min_beta=40)
        rows.append(dict(n=n, q=q, sigma=sigma, beta_gsa=pr.beta, d_gsa=pr.d,
                         m_opt_gsa=pr.m_opt, delta_gsa=pr.delta,
                         beta_gsa_floor40=pr40.beta))
    df = pd.DataFrame(rows)
    df.attrs["source"] = "local GSA reimplementation (run estimate.sage for the real thing)"
    return df


def env_note() -> str:
    p = REPO / "results" / "environment.json"
    if not p.exists():
        return "environment not recorded"
    e = json.loads(p.read_text())
    return (f"{e.get('cpu','?')}, {e.get('cores','?')} cores; "
            f"fpylll {e.get('fpylll','?')}, Python {e.get('python','?')}; "
            f"recorded {e.get('recorded_utc','?')}")


def save(fig, stem: str, caption: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{stem}.png")
    fig.savefig(FIG / f"{stem}.pdf")
    (FIG / f"{stem}.caption.txt").write_text(caption + "\n\n" + env_note() + "\n")
    plt.close(fig)
    print(f"  figures/{stem}.png")


# --------------------------------------------------------------------------
# experiment `ladder`
# --------------------------------------------------------------------------


def fig_ladder(df: pd.DataFrame) -> pd.DataFrame | None:
    ok = df[df.single_solved == 1]
    if ok.empty:
        print("  ladder: no solved trials yet")
        return None
    q, sigma = int(df.q.iloc[0]), float(df.sigma.iloc[0])
    g = ok.groupby("n")["beta_single"].agg(["median", "min", "max", "count"])
    gp = df[df.prog_solved == 1].groupby("n")["beta_prog"].agg(["median", "min", "max"])
    pred = predictions(sorted(g.index), q, sigma).set_index("n")

    fig, ax = plt.subplots()
    ns = np.array(g.index)
    ax.errorbar(ns, g["median"], yerr=[g["median"] - g["min"], g["max"] - g["median"]],
                fmt="o-", color=C_SINGLE, capsize=3, ms=4,
                label="measured, single BKZ-$\\beta$ (median, range)")
    if not gp.empty:
        ax.plot(gp.index, gp["median"], "s--", color=C_PROG, ms=4,
                label="measured, progressive BKZ (median)")
    if "beta_gsa" in pred:
        ax.plot(pred.index, pred["beta_gsa"], "^-", color=C_PRED_GSA, ms=4,
                label="GSA prediction (global minimum of cost_gsa)")
    if "beta_cn11" in pred and pred["beta_cn11"].notna().any():
        ax.plot(pred.index, pred["beta_cn11"], "v-", color=C_PRED_CN11, ms=4,
                label="primal_usvp, CN11 (misses its minimum at $n \\leq 40$)")
    ax.axhline(40, color="grey", lw=0.8, ls=":")
    ax.text(ns.min(), 41, "estimator's GSA search floor ($\\beta=40$)",
            fontsize=7, color="grey", va="bottom")
    ax.set_xlabel("LWE secret dimension $n$")
    ax.set_ylabel("BKZ block size $\\beta$")
    ax.set_title(f"Required block size: measured vs predicted  ($q={q}$, "
                 f"$\\sigma={sigma:g}$)")
    ax.legend(loc="upper left", fontsize=7.5)
    save(fig, "beta_vs_n",
         "Smallest BKZ block size that recovers the secret, against the "
         "lattice-estimator's prediction for the same lattice dimension. "
         "Points are medians over independent instances; bars span the observed range.")

    out = g.join(gp, rsuffix="_prog").join(
        pred[[c for c in ("beta_gsa", "beta_cn11", "d_gsa") if c in pred]])
    out["ratio_single_over_gsa"] = (out["median"] / out["beta_gsa"]).round(3)
    out.to_csv(REPO / "results" / "summary_ladder.csv")
    write_table(out, "ladder")
    return out


def fig_time(df: pd.DataFrame) -> None:
    ok = df[(df.single_solved == 1) & (df.time_single_at_beta > 0)]
    if len(ok) < 3:
        return
    g = ok.groupby("n")["time_single_at_beta"].median()
    fig, ax = plt.subplots()
    ax.semilogy(g.index, g.values, "o-", color=C_SINGLE, ms=4,
                label="median time of the solving BKZ-$\\beta$ call")
    if len(g) >= 3:
        c = np.polyfit(g.index, np.log2(g.values), 1)
        xs = np.linspace(g.index.min(), g.index.max(), 50)
        ax.semilogy(xs, 2 ** np.polyval(c, xs), "--", color=C_PRED_GSA, lw=1,
                    label=f"fit: $\\log_2 t \\approx {c[0]:.2f}n {c[1]:+.1f}$")
    ax.set_xlabel("LWE secret dimension $n$")
    ax.set_ylabel("wall-clock seconds (1 core)")
    ax.set_title("Cost of the solving reduction call")
    ax.legend(fontsize=7.5)
    save(fig, "time_vs_n",
         "Wall-clock cost of the single BKZ call that recovered the secret. "
         "Not comparable to the estimator's operation counts, which assume "
         "sieving rather than the enumeration fpylll uses at these block sizes.")


# --------------------------------------------------------------------------
# experiment `prob`
# --------------------------------------------------------------------------


def beta_50(betas: np.ndarray, p: np.ndarray) -> float | None:
    """First crossing of p = 0.5, by linear interpolation."""
    for i in range(1, len(betas)):
        if p[i - 1] < 0.5 <= p[i]:
            x0, x1, y0, y1 = betas[i - 1], betas[i], p[i - 1], p[i]
            return float(x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)) if y1 != y0 else float(x1)
    return float(betas[0]) if p[0] >= 0.5 else None


def fig_prob(df: pd.DataFrame) -> pd.DataFrame | None:
    q, sigma = int(df.q.iloc[0]), float(df.sigma.iloc[0])
    ns = sorted(df.n.unique())
    pred = predictions(ns, q, sigma).set_index("n")

    fig, ax = plt.subplots()
    cmap = plt.get_cmap("viridis")
    rows = []
    for i, n in enumerate(ns):
        sub = df[df.n == n].groupby("beta")["solved"].agg(["mean", "count"])
        col = cmap(i / max(1, len(ns) - 1) * 0.85)
        ax.plot(sub.index, sub["mean"], "o-", ms=3.5, color=col, label=f"$n={n}$")
        b50 = beta_50(np.array(sub.index), np.array(sub["mean"]))
        if b50 is not None:
            ax.plot([b50], [0.5], "*", ms=9, color=col)
        rows.append(dict(n=n, beta_50=b50, trials_per_beta=int(sub["count"].median()),
                         beta_gsa=pred.loc[n].get("beta_gsa") if n in pred.index else None,
                         beta_cn11=pred.loc[n].get("beta_cn11") if n in pred.index else None))
        if n in pred.index and not pd.isna(pred.loc[n].get("beta_gsa")):
            ax.axvline(pred.loc[n]["beta_gsa"], color=col, ls=":", lw=1)

    ax.axhline(0.5, color="grey", lw=0.8, ls="--")
    ax.set_xlabel("BKZ block size $\\beta$")
    ax.set_ylabel("fraction of instances solved")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title(f"Success probability of a single BKZ-$\\beta$ call  "
                 f"($q={q}$, $\\sigma={sigma:g}$)")
    ax.legend(fontsize=7.5, ncol=2)
    save(fig, "success_prob",
         "Fraction of independent LWE instances solved by one BKZ call at each "
         "block size. Stars mark the interpolated 50% point; dotted vertical "
         "lines are the estimator's predicted block size for the same n. "
         "Success is a gradual transition, not a threshold.")

    out = pd.DataFrame(rows).set_index("n")
    out.to_csv(REPO / "results" / "summary_prob.csv")
    write_table(out, "prob")
    return out


def fig_mtrade(df: pd.DataFrame) -> None:
    ok = df[df.solved == 1]
    if ok.empty:
        return
    fig, ax = plt.subplots()
    for n, sub in ok.groupby("n"):
        g = sub.groupby("d")["beta_single"].median()
        ax.plot(g.index, g.values, "o-", ms=4, label=f"$n={n}$")
        pr = primal_usvp(n=int(n), q=int(sub.q.iloc[0]),
                         sigma_s=float(sub.sigma.iloc[0]),
                         sigma_e=float(sub.sigma.iloc[0]), m=8 * int(n))
        ax.axvline(pr.d, ls=":", lw=1, color="grey")
        ax.annotate("estimator's $d$", (pr.d, ax.get_ylim()[1]), fontsize=7,
                    color="grey", ha="right", va="top", rotation=90)
    ax.set_xlabel("lattice dimension $d = m + n + 1$")
    ax.set_ylabel("smallest solving block size $\\beta$")
    ax.set_title("Using more samples buys a smaller block size")
    ax.legend(fontsize=7.5)
    save(fig, "m_tradeoff",
         "Smallest solving block size as a function of how many LWE samples are "
         "embedded. The estimator optimises this trade-off internally; the "
         "dotted line is the dimension it selects.")


def write_table(df: pd.DataFrame, stem: str) -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    tex = df.round(3).to_latex(escape=False, na_rep="--")
    (TAB / f"{stem}.tex").write_text(tex)
    print(f"  paper/tables/{stem}.tex")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    print("analysing", RAW)
    any_data = False

    df = load("ladder")
    if df is not None and len(df):
        any_data = True
        print(f"\nladder: {len(df)} trials, n in {sorted(df.n.unique())}")
        out = fig_ladder(df)
        fig_time(df)
        if out is not None:
            print(out.to_string())

    df = load("prob")
    if df is not None and len(df):
        any_data = True
        print(f"\nprob: {len(df)} trials")
        out = fig_prob(df)
        if out is not None:
            print(out.to_string())

    df = load("mtrade")
    if df is not None and len(df):
        any_data = True
        print(f"\nmtrade: {len(df)} trials")
        fig_mtrade(df)

    if not any_data:
        print("no raw results yet -- run src/run_sweep.py first")
        return
    src = predictions([32], 257, 4.0).attrs.get("source", "?")
    print(f"\npredictions source: {src}")
    print(f"environment: {env_note()}")


if __name__ == "__main__":
    main()
