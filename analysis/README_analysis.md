# Analysis scripts

These scripts compute every table, figure and number in the paper from the
committed data (`results/raw/*.csv` and `results/scan_n*.csv`). None of them
runs BKZ. Run them from this directory (`analysis/`); they expect `results/`
and `src/` one level up.

## Requirements

All scripts need numpy, pandas, scipy and matplotlib. `model_predictions.py`
and `gs_profiles.py` also need `fpylll`, and `model_predictions.py` needs a
local clone of the simulator of Postlethwaite and Virdia
(<https://github.com/fvirdia/usvp-simulation>, commit `333416d`), whose path is
passed with `--pv21`.

## Scripts, in order

### Main sweep and ladder (8 tours)

| Script | Output | Used in |
|---|---|---|
| `ladder_stats.py` | Ladder statistics; paired comparison of progressive and single-β BKZ | Sections 5.1 and 5.5, Tables 3 and 7 |
| `sweep_stats.py` | Success fraction p(β), four β₅₀ estimators, bootstrap intervals (slow) | Section 5.2, Table 4 |
| `timing_stats.py` | Cost scaling | Section 5.6, Figure 3 |
| `replay_ladder.py` | Exact replay of the ladder search on the sweep outcomes; reproducibility check | Section 5.3, Table 5 |
| `est_table.py` | `primal_usvp` output, best point with β ≥ 40 and global minimum, from the scans | Section 5.7, Table 8 |
| `make_numbers.py` | `final_numbers.json`, the single source for the tables | All tables |
| `extra_stats.py` | Stratified pooled interval, estimator ranges, reversals, tour counts, hand check of the success condition | Sections 2.4 and 5.2–5.4 |
| `replay_local_minimum.py` | Replay of the estimator's `local_minimum` search over the scans | Section 5.7 |

### Tour budget

| Script | Output | Used in |
|---|---|---|
| `tour4_stats.py` | 4-tour success curves (input to `make_figs.py`) | Figure 2 |
| `tour_budget.py` | β₅₀, β₉₀ and successes at the prediction for every budget from 1 to 8 tours (5–10 minutes) | Section 5.4 |
| `tour_convergence.py` | The same for budgets up to 32 tours, from `results/raw/long.csv` (see `../RUN_LONG.md`) | Section 5.4, Table 6 |
| `make_fig_tours.py` | Figure 5, from `tour_convergence.json` | Figure 5 |

### Explanations of success below the prediction

| Script | Output | Used in |
|---|---|---|
| `target_norm.py` | Length of the planted vector against success; regenerates instances from their seeds | Section 5.8 |
| `model_predictions.py --pv21 PATH` | Postlethwaite–Virdia simulator and the CN11 condition against the measurements | Section 5.8, Figure 6 |
| `gs_profiles.py` | Measured Gram–Schmidt profiles against the GSA and CN11 | Section 5.8, Figure 7 |

### Figures

| Script | Output |
|---|---|
| `make_figs.py` | Figures 1–4 (requires `tour4_stats.json`). File names differ from the paper's numbering: `fig1_beta_vs_n.png` is Figure 2, `fig2_success_prob.png` is Figure 1, `fig3_time_vs_n.png` is Figure 3 and `fig4_cost_curve_n40.png` is Figure 4. |

`anchor_sim.py` is superseded by `replay_ladder.py`. It assumed independent
probes, whereas the data show that outcomes are fixed for each instance and
block size. It is kept for reference only.

## Provenance

These scripts were written with AI assistance (Claude) and run by the author,
who is responsible for the code and its results. The statistical steps
(logistic fit, isotonic regression, bootstrap over instances, search replay)
have not been independently reviewed.
