# Empirical hardness of toy LWE: BKZ attacks vs. the LWE Estimator

Small, deliberately weak LWE instances are generated, embedded into lattices,
attacked with BKZ via `fpylll`, and the block size that actually breaks them is
compared against the block size the
[lattice-estimator](https://github.com/malb/lattice-estimator) predicts.

Everything here is a controlled experiment on instances we generate ourselves.
No deployed scheme or third-party system is involved.

---

## 1. Install

### Option A — conda/mamba (recommended: one environment, everything works)

SageMath ships `fpylll` **and** is what the estimator needs, so this avoids the
pip problems described in §2.

```bash
# install micromamba first if you don't have conda: https://mamba.readthedocs.io
mamba create -n lwe -c conda-forge sage python=3.11 numpy pandas matplotlib scipy
mamba activate lwe

git clone https://github.com/malb/lattice-estimator.git   # into the repo root
git -C lattice-estimator checkout 53da598   # version used in the paper
python3 src/strategies.py        # should print a real path, not "none"
```

### Option B — pip only (no Sage; predictions come from the local copy)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./scripts/fetch_strategies.sh    # required, see §2
```

Option B cannot run `src/estimate.sage`. Use it to build and debug the
experiment, but the paper's headline comparison must use the real estimator —
either switch to Option A, or run the estimator in
[CoCalc](https://cocalc.com) or the `sagemath/sagemath` Docker image and copy
`results/predictions.csv` back.

### Verify before running anything

```bash
python3 src/lwe_instance.py    # instance generation + self-check
python3 src/embedding.py       # must print "target in lattice: True"
python3 src/gsa_estimate.py    # must print "MATCH" for Kyber512
python3 src/strategies.py      # must print a path, not "none"
```

`gsa_estimate.py` printing `MATCH` means our reimplementation of the
estimator's GSA predicate reproduces its published Kyber-512 output
(β=406, d=998, δ=1.003941). If that breaks, stop and fix it before trusting
any prediction in this repository.

---

## 2. Two installation traps

**`pip install fpylll` gives you a BKZ that cannot run.** The wheels ship with
`BKZ.DEFAULT_STRATEGY` pointing at the *build machine's* filesystem
(`/project/local/share/fplll/strategies/default.json`), which does not exist on
yours. Any `BKZ.Param(strategies=BKZ.DEFAULT_STRATEGY)` then raises
`RuntimeError: Cannot open strategies file`. `scripts/fetch_strategies.sh`
downloads the file from the fplll repository and `src/strategies.py` resolves
it. Without a strategies file BKZ still runs correctly but with no
preprocessing or pruning, which is several times slower — so **never mix timing
measurements taken with and without it.** Every run records which was used in
`results/environment.json`.

**`pip install fpylll` alone fails at import.** It needs `cysignals`, which is
not pulled in automatically. It is in `requirements.txt`.

---

## 3. Run the experiments

```bash
# predictions for the grid (needs Sage)
sage src/estimate.sage                       # -> results/predictions.csv

# headline: success probability vs block size
python3 src/run_sweep.py prob   --ns 32 36 40 44 --trials 20 --jobs 8

# block-size ladder + wall-clock scaling
python3 src/run_sweep.py ladder --ns 32 36 40 44 48 --trials 10 --jobs 8

# sample/dimension trade-off
python3 src/run_sweep.py mtrade --ns 40 --trials 5 --jobs 8

# tour-budget control: same instances, one BKZ-beta call with up to 32 tours,
# recording the tour of success (see RUN_LONG.md)
python3 src/run_sweep.py long --ns 36 32 40 --trials 20 --jobs 7
python3 src/run_sweep.py long --ns 44 --trials 20 --jobs 7 --beta-above 6
python3 src/run_sweep.py long --ns 36 --trials 20 --jobs 7 --beta-window 20   # adds beta = 14, 16

# basis profiles: log2 ||b*_i|| after LLL and after 1, 2, 4, 8, 16 tours,
# at beta = prediction and prediction - 8 (see RUN_V5.md)
python3 src/run_sweep.py profile --ns 32 36 40 44 --trials 20 --jobs 7
```

### Analysis

Every table, figure and number in the paper is produced by the scripts in
`analysis/` from the committed CSVs, without running BKZ. The order and what
each script produces are in [`analysis/README_analysis.md`](analysis/README_analysis.md).
Figures are written to `analysis/figs/`.

`python3 src/analyze.py` is an older quick-look summary (`figures/`, ignored by
git, and `results/summary_*.csv`). Its GSA column is the global minimum of
`cost_gsa` from the scans, not `primal_usvp`'s answer (see §6).

Sweeps are **resumable**: re-running the same command runs only the trials
missing from the CSV, so interrupt them freely. Use `--dry-run` to see the job
count before committing compute. `--jobs` should be your physical core count
minus one.

---

## 4. Parameters, and why these

Default family: `q = 257`, secret and error both discrete Gaussian with
`σ = 4`, `n` swept from 32 upward, `m` set to the number of samples the
estimator's own optimum uses.

- **Secret and error share a distribution** so the estimator's rescaling factor
  `ξ = σ_e/σ_s` is exactly 1. The embedding basis then stays integral and is
  *the same lattice* the estimator models. Changing one distribution and not
  the other breaks that correspondence — read `src/embedding.py` first.
- **`σ` is large relative to `q`** compared to a deployed scheme (our noise rate
  `α = σ√(2π)/q ≈ 0.039` versus Kyber-512's `≈ 0.00092`). Toy dimensions are
  only hard if the noise rate is raised. State this as a limitation: the regime
  being probed is not the deployed one.
- **`n ≥ 32`** because below that the estimator's optimal `m` yields instances
  that are not well-posed (see §5), and the run aborts those points with a
  printed reason.

`src/gsa_estimate.py` run as a script prints a design table so you can choose a
different `(q, σ)` family: aim for predicted `β` in 40–80, lattice dimension
`d ≤ ~180`, and `gap ≤ 0.5`.

---

## 5. Two ways this experiment silently produces nonsense

Both are guarded in code. Do not remove the guards.

**Instances with too few samples are not uniquely decodable.** With `m ≈ n`,
an instance has astronomically many secrets whose error is small, LLL alone
finds one, and every dimension looks "broken" instantly. The guard is
`usvp_gap = ‖t‖ / λ₁(Gaussian heuristic)`, which must be well below 1;
`run_sweep.py` refuses any parameter point above `GAP_LIMIT = 0.6` and says
which and why.

**A loose success test accepts wrong answers.** `embedding.py` counts a trial as
solved only if the recovered secret is *exactly* the planted one, or if its
error is no larger in ℓ₂ than the planted error. In every recorded success the
recovered secret was the planted one (column `exact`).

---

## 6. Three things the estimator does not claim

Carry these into the write-up; each is checkable in the estimator's source.

1. **It cannot report `β < 40`.** Under the GSA shape model, `primal_usvp`
   searches block sizes from 40 upward (`local_minimum(40, max_beta, ...)` in
   `PrimalUSVP.__call__`). A measured `β` below 40 at such a grid point is
   outside the tool's intended range, not evidence against it. Separately,
   at n = 32, 36 and 40 its search also misses the cheapest point *inside*
   that range (paper, Section 5.7; lattice-estimator issue
   [#225](https://github.com/malb/lattice-estimator/issues/225)). The
   predictions used throughout are therefore the global minima of `cost_gsa`
   from exhaustive scans (`results/scan_n*.csv`), not `primal_usvp`'s output.
2. **Its operation counts assume sieving** (MATZOV, BDGL16, ADPS16 …), while
   `fpylll` uses enumeration with pruning at these block sizes. Wall-clock
   times are therefore *not* comparable to `rop`. The comparable quantity is
   the **block size**.
3. **It predicts a single BKZ-β call**, not a progressive ladder. Progressive
   BKZ is strictly stronger, so comparing progressive results against the
   estimator's β understates the estimator. `run_sweep.py` measures both and
   labels them separately.

---

## 7. Layout

```
src/lwe_instance.py        seeded instance generation, fingerprinting
src/embedding.py           Kannan embedding, lattice-membership test, success test
src/attack.py              progressive BKZ, single-β mode, β search, long and profile modes
src/gsa_estimate.py        Sage-free reimplementation of the GSA predicate
src/run_sweep.py           experiment driver: resumable, parallel, records provenance
src/estimate.sage          official estimator over the grid -> results/predictions.csv
src/confirm_primal_usvp.sage, src/check_kyber_usvp.sage   checks for Section 5.7
src/analyze.py             quick-look summary (not the paper's figures)
analysis/                  scripts that produce the paper's tables, figures and numbers
analysis/figs/             the paper's figures (PNG, PDF)
results/raw/*.csv          one row per trial, append-only (ladder, prob, long, profiles)
results/scan_n*.csv        exhaustive cost_gsa scans
results/environment*.json  CPU, core count, library versions, git commit, per run
MANIFEST.sha256            checksums of code and data: sha256sum -c MANIFEST.sha256;
                           rebuild with scripts/make_manifest.sh after any change
PROTOCOL.md                the protocol, and the deviations from it
```

Raw CSVs, scans and the paper's figures are committed deliberately: a reviewer
can regenerate every figure and number without running a single BKZ call.

### Provenance of the runs

Each `results/environment*.json` records the commit the run used. The main
sweep and ladder (`environment.json`, 2026-09-23) were run from commit
`0f3cce2`, with uncommitted changes. That commit is not in this repository's
history; the first commit here is `99117b0` (2026-10-13). The long-tour run
used `ab5a2b9`, also with uncommitted changes, and reproduces the 8-tour
outcome of all 1,260 instance/block-size pairs the two runs share (the 10
sweep calls that timed out are left out). The profile run used
`581bfca`, clean.

---

## 8. Reproducing a single number

Each result row carries `seed` and a `fingerprint` (a hash of `A` and `b`), so
any trial can be replayed exactly:

```python
from lwe_instance import gen_lwe
from attack import solve_single
inst = gen_lwe(n=40, q=257, m=88, seed=2, dist_s="dg", param_s=4,
               dist_e="dg", param_e=4)
print(inst.fingerprint())          # must match the CSV row
print(solve_single(inst, beta=46, tours=8))
```

The recorded outcome of this call is solved, at tour 5. In our runs the outcome
of a call is reproducible for a given instance and block size: the long-tour
run repeats the 8-tour calls of the sweep and agrees on all 1,260 shared pairs
(timeouts left out),
and on the tour of success in 694 of 708 solved cases. Success probability is
therefore a fraction of instances, not a chance per call; the paper reports it
that way and resamples instances for its intervals.
