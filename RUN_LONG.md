# Long-tour control run

## Purpose

With 8 BKZ tours per call, the block size at which half of the instances are
solved lies below the prediction; with 4 tours it lies above the prediction at
three of the four sizes, and success at the predicted β was still rising at the
eighth tour. This run attacks the same 20
instances per n, on the same lattices and β grid, with up to **32 tours** per
call. It also records when fpylll's own stopping rules (`AUTO_ABORT` and a
clean tour) would have ended each call; the estimator's δ₀ table for β ≤ 40 was
measured with that rule. Each row therefore gives the outcome for every budget
from 1 to 32 tours.

The run writes `results/raw/long.csv` and `results/environment_long.json`. It
does not modify `prob.csv`, `ladder.csv` or `environment.json`.

## Requirements

- The same environment as the main sweep (Python 3.14, fpylll 0.6.4) and the
  same strategies file (`data/bkz_strategies.json`). The analysis checks that
  the first 8 tours reproduce `prob.csv`; a different fpylll version or
  strategies file may break that check.
- A clean working tree, so that the run records a real commit. The recorded run
  used `ab5a2b9` with uncommitted changes (README, Section 6).
- A machine that will not sleep during the run.

## Commands

From the repository root:

```bash
# n = 32, 36, 40 on the full grid (prediction ± 16): 1,020 calls, about 1 hour
python3 src/run_sweep.py long --ns 36 32 40 --trials 20 --jobs 7

# n = 44 up to prediction + 6 (its high-β calls are very slow): 240 calls, about 2 hours
python3 src/run_sweep.py long --ns 44 --trials 20 --jobs 7 --beta-above 6
```

After the first analysis, the 32-tour 50% point at n = 36 lay at the bottom
edge of the grid, so the n = 36 grid was extended down to β = 14 and 16
(40 further calls). For this reason the 8-tour value at n = 36 in Figure 5
differs slightly from Table 4.

```bash
python3 src/run_sweep.py long --ns 36 --trials 20 --jobs 7 --beta-window 20
```

`--dry-run` prints the number of jobs without running them. An interrupted run
resumes when the same command is repeated. The time estimates above are upper
bounds from `prob.csv`, assuming that every call failing at 8 tours runs all
32. Each call is limited to 7,200 s, checked between tours.

## Analysis

```bash
cd analysis
python3 tour_convergence.py      # about 10–15 minutes (bootstrap)
python3 make_fig_tours.py        # Figure 5
```

The first line of the output must state that all outcomes agree with
`prob.csv` within 8 tours. If any disagree, outcomes are not reproducible in
that environment and the derived budgets are not exact.

The script reports, for budgets of 4, 8, 12, 16, 24 and 32 tours and for
fpylll's default stopping rule: β₅₀ and β₉₀ with intervals; success at the
prediction; the tour at which fpylll would have stopped; the cumulative number
of instances solved at the prediction, by tour; and the number of non-monotone
steps in each success curve. The results are reported in Section 5.4 and
Figure 5 of the paper.
