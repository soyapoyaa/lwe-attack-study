# Analysis scripts

Every number and figure in the paper is computed by these scripts from `results/raw/*.csv` and `results/scan_n*.csv`.
Run from this directory (`analysis/`), with `results/` and `src/` one level up.

    python3 ladder_stats.py          # Table 2 inputs, paired progressive-vs-single comparison
    python3 sweep_stats.py           # p(beta), four beta_50 estimators, bootstrap intervals (slow)
    python3 timing_stats.py          # cost scaling (Section 5.5)
    python3 anchor_sim.py            # simulation of the ladder search (Section 5.3, Table 4)
    python3 est_table.py             # Table 6 from the scans
    python3 make_numbers.py          # -> final_numbers.json (single source of truth for tables)
    python3 extra_stats.py           # stratified pooled interval, estimator ranges, reversals, tour counts, hand check of the success condition
    python3 tour4_stats.py           # derived 4-tour success curves (assumes a 4-tour call repeats the first 4 tours of an 8-tour call)
    python3 make_figs.py             # -> figs/*.png (needs tour4_stats.json)
    python3 replay_local_minimum.py  # replay of the estimator's local_minimum over the scans

Requires numpy, pandas, scipy, matplotlib.

## Provenance and responsibility
These scripts were written with AI assistance (Claude, in the conversation that produced the draft). They have been run, but
not independently reviewed. Before relying on them, read them, and test the statistical steps (logistic fit, isotonic
regression, bootstrap over instances, search simulation) yourself. The author is responsible for the code and its results.
`replay_local_minimum.py` comes from earlier work in the same project; check its provenance before describing it in the paper.
