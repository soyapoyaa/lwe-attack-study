# Protocol

Pre-registered description of the experiment, written before interpreting
results. It describes what the code in `src/` actually does; if you change the
code, change this file in the same commit, and note the change in the paper.

Frozen on: [DATE]  ·  Commit: [HASH]
[AUTHOR: fill in when and in which commit this file was frozen. The first commit
of this repository that contains it is 99117b0 (2026-10-13); the main sweep ran
on 2026-09-23 from commit 0f3cce2 of an earlier local repository. If the protocol
was not frozen before the sweep, say so here rather than giving a date.]

---

## 1. Instance family

| quantity | value |
|---|---|
| modulus `q` | 257 (prime) |
| secret distribution `Xs` | rounded Gaussian, σ = 4, coefficient-wise i.i.d. |
| error distribution `Xe` | rounded Gaussian, σ = 4, coefficient-wise i.i.d. |
| secret dimension `n` | 32, 36, 40, 44, 48 (52 if compute allows) |
| samples available | 8n |
| samples used `m` | the estimator's own optimum, `d_gsa − n − 1` |
| `A` | uniform in `Z_q^{m×n}` |

`Xs = Xe` is deliberate: it makes the estimator's rescaling factor
`ξ = σ_e/σ_s` exactly 1, so the embedding basis is integral and is the same
lattice the estimator models. A rounded Gaussian is used rather than an exact
discrete Gaussian; this is reported as such.

Each instance is determined entirely by `(n, q, m, seed)` and is identified in
the results by a SHA-256 prefix of `(A, b)`.

## 2. Reduction to uSVP

Kannan embedding, dimension `d = m + n + 1`, row basis

```
[  q·I_m      0        0  ]
[  −Aᵀ      ξ·I_n      0  ]
[   bᵀ        0        τ  ]
```

with `ξ = max(1, σ_e/σ_s) = 1` and `τ = round(σ_e) = 4`, matching
`PrimalUSVP._xi_factor` and `cost_gsa` in the estimator. The target is
`t = (e, ξ·s, τ)`; `src/embedding.py` run as a script verifies by explicit
integer combination that `t` lies in the lattice.

## 3. Success criterion

A trial is solved if and only if some row of the reduced basis has last
coordinate `±τ` and yields a candidate secret that is **either** exactly the
planted secret **or** has error no larger in ℓ₂ norm than the planted error.

Nothing weaker is accepted. A bound of the form "error is small" admits wrong
answers on instances with many valid secrets.

## 4. Well-posedness guard

Define `gap = ‖t‖ / λ₁(Gaussian heuristic)` where
`λ₁ ≈ √(d/2πe) · vol(L)^{1/d}` and `vol(L) = q^m · ξ^n · τ`.

Parameter points with `gap > 0.6` are **not run**, and the rejection is printed
with its reason. Above that threshold the planted vector is not the unique
shortest vector, the instance has many small-error secrets, and LLL alone
"solves" it — a meaningless measurement.

`n = 28` is rejected under this rule at the estimator's own optimal `m`.

## 5. Reduction algorithm

LLL (fplll defaults), then BKZ 2.0 (`fpylll.algorithms.bkz2.BKZReduction`) with
fplll's default pruning and preprocessing strategies, flags
`GH_BND | AUTO_ABORT`, at most 8 tours per block size, checking for a solution
after every tour. Floating-point type by dimension: `d` (d ≤ 90), `ld`
(≤ 180), `dd` (above); MPFR precision 160.

Two modes:

- **single-β** - fresh LLL-reduced basis, one block size only. This is what the
  estimator models and is the mode used for every comparison against it.
- **progressive** - one basis, ascending ladder of block sizes. This is what a
  real attacker runs. Reported separately and never compared against the
  estimator's β as though equivalent.

## 6. Locating the smallest solving β

β grid: even values, 4 to 80. Exponential search outward from the estimator's
predicted β to bracket the transition, then binary search within the bracket.

Success is expected to be monotone in β but is not guaranteed, pruning is
randomised so the reported per-instance value is "the smallest β located by
this search", and the headline observable is a success *probability*, not a
per-instance threshold.

## 7. Observables

| symbol | definition |
|---|---|
| `p(β)` | fraction of independent instances solved by one BKZ-β call |
| `β₅₀` | interpolated β where `p = 0.5`; **the quantity compared to the estimator** |
| `β_single` | smallest solving β for one instance, single-β mode |
| `β_prog` | smallest solving β for one instance, progressive mode |
| `t_β` | wall-clock seconds of the solving call; reported for scaling only |

Wall-clock time is never compared against the estimator's operation counts:
those assume sieving, fpylll enumerates.

## 8. Trial counts

- `prob`: 20 independent instances per `(n, β)`, β within ±16 of the prediction.
- `ladder`: 10 independent instances per `n`.
- `mtrade`: 5 independent instances per `(n, m)`.

Minimum acceptable for `prob` is 10; below that the spread is not measurable
and the headline observable is not meaningful.

## 9. Predictions

From `lattice-estimator` via `sage src/estimate.sage`, recording
`primal_usvp` under the GSA, CN11 and ZGSA basis-shape models, plus
`primal_bdd` and `dual` for context, and `rop` under ADPS16, BDGL16, MATZOV and
CheNgu12 cost models.

Known constraint: under GSA, `primal_usvp` searches block sizes from 40 upward,
so it cannot report β < 40. Grid points whose true β falls below 40 are outside
the tool's range and a measured β below 40 there is **not** evidence against it.

## 10. Provenance

Every sweep writes `results/environment.json` with CPU model, core count,
Python/numpy/fpylll versions, the BKZ strategies file in use, and the git
commit plus dirty flag. One CSV row per trial, append-only; sweeps are
resumable and never overwrite. All figures regenerate from committed CSVs with
no lattice reduction.

## 11. Deviations from this protocol

Record each one here with a date and a reason. An empty section is the goal;
a populated one is honest.

- 2026-09-25: No deviations to the (n, m) grid, success criteria, or trial
  counts above. One methodological finding surfaced during verification of
  the estimator predictions in §9, recorded here for traceability.

  `LWE.primal_usvp()`'s bounded search over β does not always find the global
  minimum of its own `cost_gsa` cost function. Verified by exhaustive scan
  (β = 4 .. m+n, calling `cost_gsa` directly) for all five grid points:

      n     primal_usvp() output      true global minimum (exhaustive scan)
      32    β=65 (d=66)               β=20 (d=104)   -- below the β≥40 floor (§9, expected)
      36    β=73 (d=74)               β=34 (d=121)   -- below the β≥40 floor (§9, expected)
      40    β=81 (d=82)               β=44 (d=129)   -- NOT below the floor: search bug
      44    β=54 (d=136)              β=54 (d=136)   -- matches
      48    β=63 (d=145)              β=63 (d=145)   -- matches

  n=40's true minimum is above the documented floor, and the bounded search
  still misses it. This is a `primal_usvp()` search limitation, not a floor
  case. All GSA predictions used elsewhere in this paper are the
  exhaustively-verified minima above, not raw `primal_usvp()` output.
  Reported upstream: https://github.com/malb/lattice-estimator/issues/225

- 2026-10-16: Added, after the 8-tour results had been seen: a long-tour
  control run (`run_sweep.py long`, `results/raw/long.csv`, RUN_LONG.md). Same
  instances, lattices and β grid as `prob`, one BKZ-β call with up to 32 tours,
  recording the tour of success and when fpylll's own stopping rules would have
  stopped. Reason: success at the predicted β was still rising at the eighth
  tour, and the 8-tour budget is not fixed by the estimator. Deviations from
  the `prob` grid: n = 44 only up to the prediction + 6 (its high-β calls are
  very slow); n = 36 extended down to β = 14 because the 32-tour 50% point was
  at the edge of the grid. The run reproduces the 8-tour outcome of all 1,260
  shared (instance, β) pairs (the 10 `prob` calls that timed out excluded).

- 2026-10-23: Recorded now, decided while analysing the 8-tour sweep: three
  changes to how §6 and §7 are applied.
  1. β₅₀. §7 defines it by linear interpolation at the first crossing of
     p = 0.5. The measured curves are not monotone, so the first crossing is
     unstable. The paper reports four estimators (first crossing, last
     crossing, isotonic regression, logistic maximum likelihood) and uses the
     logistic fit as the main one, with bootstrap intervals over instances.
  2. Timeouts. The 10 `prob` calls that hit the time limit are treated as
     unknown and left out, not counted as failures; counting them as failures
     creates a false dip in p(β).
  3. 4-tour curve. Not a separate run: a trial counts as solved with 4 tours
     only if it was solved within its first 4 tours, read off the recorded
     tour of success. This is valid because the code checks for the secret
     after every tour and nothing in a tour depends on the total budget.
  Also, §6 says pruning is randomised and outcomes may vary per call. In these
  runs the outcome reproduces for a given instance and β (2026-10-16 entry), so
  p(β) is a fraction of the 20 instances.

- 2026-10-23: Correction to the 2026-09-25 entry. At n = 32 and 36 the miss is not
  only the β ≥ 40 floor. Exhaustive scans show `primal_usvp()` also misses the
  cheapest point inside its own range, by 2.0 bits at n = 32 (best β ≥ 40:
  β = 40, d = 78) and 3.6 bits at n = 36 (β = 40, d = 99); the floor accounts
  for a further 11.6 and 1.4 bits (paper, Table 8). That entry's table is left as
  recorded. The CN11 mode inherits the miss at n = 32, 36 and 40, because it
  searches only a window around the GSA answer (`results/predictions.csv`:
  β = 59, 66, 72). Every prediction in the paper is the global minimum of
  `cost_gsa` from the scans, as stated in that entry.

- 2026-10-23: Clarification of §8. The 20 instances per n are the same at every
  β (seeds 1–20), so the success curve is a fraction of one set of 20
  instances, and outcomes are reproducible per (instance, β) (see the 2026-10-16
  entry). The paper's intervals resample instances, not calls.

- 2026-10-23: Added, after the results above had been seen: a basis-profile
  run (`run_sweep.py profile`, `results/raw/profiles.csv`, 160 calls at β = the
  prediction and the prediction − 8, up to 16 tours) and comparisons with the
  success-probability simulator of Postlethwaite and Virdia and with the CN11
  shape model (RUN_V5.md). Reason: to test explanations for success below the
  prediction. These are exploratory analyses.
