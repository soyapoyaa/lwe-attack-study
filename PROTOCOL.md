# Protocol

Pre-registered description of the experiment, written before interpreting
results. It describes what the code in `src/` actually does; if you change the
code, change this file in the same commit, and note the change in the paper.

Frozen on: [DATE]  ·  Commit: [HASH]

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

- **single-β** — fresh LLL-reduced basis, one block size only. This is what the
  estimator models and is the mode used for every comparison against it.
- **progressive** — one basis, ascending ladder of block sizes. This is what a
  real attacker runs. Reported separately and never compared against the
  estimator's β as though equivalent.

## 6. Locating the smallest solving β

β grid: even values, 4 to 80. Exponential search outward from the estimator's
predicted β to bracket the transition, then binary search within the bracket.

Success is expected to be monotone in β but is not guaranteed — pruning is
randomised — so the reported per-instance value is "the smallest β located by
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
  Reported upstream: [GitHub issue URL once filed].
