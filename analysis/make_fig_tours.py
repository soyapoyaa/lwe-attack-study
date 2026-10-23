"""Figure: 50% point relative to the prediction, against the tour budget (long-tour run).

Reads tour_convergence.json (from tour_convergence.py, i.e. results/raw/long.csv).
Writes figs/fig5_tour_budget.png and .pdf. Same style and colours as make_figs.py.
Run from analysis/:  python3 make_fig_tours.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BLU, ORG, GRN, VER, GRY, INK = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#7f7f7f", "#222222"
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": .25, "grid.linewidth": .5,
                     "figure.dpi": 200, "savefig.dpi": 200, "legend.frameon": False})

J = json.load(open(os.path.join(HERE, "tour_convergence.json")))["budgets"]
GSA = {"32": 20, "36": 34, "40": 44, "44": 54}
BUD = [4, 8, 12, 16, 24, 32]
STYLE = {"32": (BLU, "o"), "36": (ORG, "s"), "40": (GRN, "D"), "44": (VER, "^")}
OFF = {"32": -0.45, "36": -0.15, "40": 0.15, "44": 0.45}   # small x-offsets so intervals don't overlap

fig, ax = plt.subplots(figsize=(6.4, 3.9))
ax.axhline(1.0, color=INK, lw=1.2, ls="--", zorder=1)
ax.text(32.6, 1.012, "prediction", ha="right", va="bottom", fontsize=7.5, color=INK)
ax.axvline(8, color=GRY, lw=.8, ls=":", zorder=0)
ax.text(8.3, 1.235, "8 tours: work charged by several\nestimator cost models (8·d SVP calls)",
        fontsize=7, color=GRY, va="top")
for n, (col, mk) in STYLE.items():
    g = GSA[n]
    y = [J[str(t)][n]["b50"] / g for t in BUD]
    lo = [y[i] - J[str(t)][n]["b50_ci"][0] / g for i, t in enumerate(BUD)]
    hi = [J[str(t)][n]["b50_ci"][1] / g - y[i] for i, t in enumerate(BUD)]
    x = [t + OFF[n] for t in BUD]
    ax.errorbar(x, y, yerr=[lo, hi], fmt=mk + "-", color=col, lw=1.4, ms=5.5, capsize=2.5,
                mec="white", mew=.6, label=f"n = {n}  (prediction β = {g})", zorder=3)
    ax.text(33.0, {"44": .815, "32": .775, "40": .735, "36": y[-1]}[n], f"n = {n}", fontsize=7.5, color=INK, va="center", ha="left")
ax.set_xticks(BUD)
ax.set_xlim(2.5, 35.5)
ax.set_xlabel("tours per BKZ call (budget)")
ax.set_ylabel("50% point ÷ predicted block size")
ax.set_title("The 50% point falls as the tour budget grows (q = 257, σ = 4)", fontsize=10)
ax.legend(loc="upper right", fontsize=7.4)
fig.tight_layout()
os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, "figs", f"fig5_tour_budget.{ext}"))
print("wrote figs/fig5_tour_budget.png/.pdf")
