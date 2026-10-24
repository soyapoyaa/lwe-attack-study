# Analysis scripts

Every number and figure in the paper is computed by these scripts from `results/raw/*.csv` and `results/scan_n*.csv`.
Run from this directory (`analysis/`), with `results/` and `src/` one level up.

    python3 ladder_stats.py          # Table 2 inputs, paired progressive-vs-single comparison
    python3 sweep_stats.py           # p(beta), four beta_50 estimators, bootstrap intervals (slow)
    python3 timing_stats.py          # cost scaling (Section 5.5)
    python3 replay_ladder.py         # exact replay of the ladder search on sweep outcomes (Section 5.3, Table 5); reproducibility check
    python3 anchor_sim.py            # SUPERSEDED by replay_ladder.py (assumed independent probes; the data show outcomes are fixed per instance)
    python3 est_table.py             # Table 6 from the scans
    python3 make_numbers.py          # -> final_numbers.json (single source of truth for tables)
    python3 extra_stats.py           # stratified pooled interval, estimator ranges, reversals, tour counts, hand check of the success condition
    python3 tour4_stats.py           # 4-tour success curves (kept for make_figs.py; tour_budget.py gives the same numbers)
    python3 tour_budget.py           # every budget 1-8 tours: beta_50, beta_90, solved-by-tour at the prediction (Section 5.4); about 5-10 minutes
    python3 tour_convergence.py      # analyses results/raw/long.csv from `run_sweep.py long` (see ../RUN_LONG.md)
    python3 make_fig_tours.py        # Figure 5 (50% point against the tour budget), from tour_convergence.json
    python3 target_norm.py           # length of the planted vector against success (Section 5.8); regenerates instances from their seeds
    python3 model_predictions.py --pv21 PATH   # Postlethwaite-Virdia simulator and CN11 condition vs measurements, Figure 6; needs fpylll and a clone of github.com/fvirdia/usvp-simulation (see ../RUN_V5.md)
    python3 gs_profiles.py           # measured Gram-Schmidt profiles vs GSA and CN11, Figure 7; needs fpylll
    python3 make_figs.py             # -> figs/*.png (needs tour4_stats.json)
    python3 replay_local_minimum.py  # replay of the estimator's local_minimum over the scans

Requires numpy, pandas, scipy, matplotlib; the last two scripts also need fpylll.

## Provenance and responsibility
These scripts were written with AI assistance (Claude, in the conversation that produced the draft). They have been run, but
not independently reviewed. Before relying on them, read them, and test the statistical steps (logistic fit, isotonic
regression, bootstrap over instances, search simulation) yourself. The author is responsible for the code and its results.
`replay_local_minimum.py` comes from earlier work in the same project; check its provenance before describing it in the paper.
