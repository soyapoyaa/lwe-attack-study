# Long-tour control run

**What it tests.** The paper's headline ("the 50% point is below the prediction") holds with 8 BKZ tours but reverses with 4, and success at the predicted β was still rising at the eighth tour. This run attacks the same 20 instances per n, on the same lattices and β grid, with up to **32 tours** per call. It also records when fpylll's own stopping rules (`AUTO_ABORT`, clean tour) would have ended each call. The estimator's δ₀ table for β ≤ 40 was measured with that stopping rule. One row per trial gives the outcome for every budget from 1 to 32 tours.

It writes `results/raw/long.csv` and `results/environment_long.json`. It never touches `prob.csv`, `ladder.csv` or `environment.json`.

## Before you start

1. Use the same Python environment as the original sweep (Python 3.14, fpylll 0.6.4) and the same strategies file (`data/bkz_strategies.json`). The analysis checks that the first 8 tours reproduce `prob.csv`. A different fpylll or strategies file could break that check.
2. Commit the current code so the run records a real commit. (The recorded run used `ab5a2b9` with uncommitted changes; see the README, "Provenance of the runs".)
3. Plug the laptop in and stop it from sleeping during the run.

## Run

From the repository root:

```bash
# n = 32, 36, 40 on the full prob grid (prediction ± 16): 1,020 calls, about 1 hour
python3 src/run_sweep.py long --ns 36 32 40 --trials 20 --jobs 7

# n = 44 up to prediction + 6 (its high-β calls are very slow): 240 calls, about 2 hours
python3 src/run_sweep.py long --ns 44 --trials 20 --jobs 7 --beta-above 6
```

After the first analysis, the 32-tour 50% point at n = 36 was at the bottom edge of the grid, so the n = 36 grid was extended down to β = 14 and 16 (40 more calls; the 8-tour value in Figure 5 differs from Table 4 for this reason):

```bash
python3 src/run_sweep.py long --ns 36 --trials 20 --jobs 7 --beta-window 20
```

Add `--dry-run` to any command to see the job count without running anything. If a run is interrupted, run the same command again. It skips finished trials.

Time estimates are upper bounds from `prob.csv`, assuming every 8-tour failure runs all 32 tours. Each call's limit is 7,200 s, checked between tours.

## Analyse

```bash
cd analysis
python3 tour_convergence.py      # about 10-15 minutes (bootstrap)
```

**Check line 1 first.** It must say that all outcomes agree within 8 tours. If any disagree, outcomes are not reproducible on your setup, and the derived budgets are not exact. Don't report them as exact. Stop and look at which rows differ.

The rest reports:
- β₅₀ and β₉₀ (with intervals) and success at the prediction, for 4, 8, 12, 16, 24 and 32 tours and for the default fpylll stopping rule
- the tour at which fpylll would have stopped
- the cumulative number solved at the prediction, by tour
- the number of downward steps (the non-monotone dips) at each budget

## Then

The results are in Section 5.4 and Figure 5 of the paper (`analysis/make_fig_tours.py`).
