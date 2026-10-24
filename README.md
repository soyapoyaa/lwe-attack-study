# Empirical hardness of toy LWE: BKZ attacks vs. the LWE Estimator

This repository contains the code, raw data and analysis scripts for the paper
*Measuring the Block Size That Breaks Toy LWE: An Empirical Check of the LWE
Estimator's GSA Prediction*.

The study generates small, deliberately weak Learning With Errors (LWE)
instances, embeds each into a lattice, and attacks it with BKZ through
`fpylll`. The smallest block size that recovers the secret is then compared
with the block size predicted by the
[lattice-estimator](https://github.com/malb/lattice-estimator) for the same
lattice. All instances are generated within the study; no deployed scheme or
third-party system is attacked.

## Contents

1. [Repository layout](#1-repository-layout)
2. [Installation](#2-installation)
3. [Reproducing the results](#3-reproducing-the-results)
4. [Experimental design](#4-experimental-design)
5. [Relation to the estimator](#5-relation-to-the-estimator)
6. [Provenance and reproducibility](#6-provenance-and-reproducibility)
7. [Licence](#7-licence)

---

## 1. Repository layout

| Path | Contents |
|---|---|
| `src/lwe_instance.py` | Seeded instance generation and fingerprinting |
| `src/embedding.py` | Kannan embedding, lattice-membership test and success test |
| `src/attack.py` | Single-β BKZ, progressive BKZ, β search, and the long-tour and profile modes |
| `src/gsa_estimate.py` | Sage-free reimplementation of the estimator's GSA success condition |
| `src/run_sweep.py` | Experiment driver (resumable, parallel, records provenance) |
| `src/estimate.sage` | Estimator predictions for the parameter grid |
| `src/confirm_primal_usvp.sage`, `src/check_kyber_usvp.sage` | Checks reported in Section 5.7 of the paper |
| `src/analyze.py` | Quick-look summary; not used for the paper's figures |
| `analysis/` | Scripts that produce every table, figure and number in the paper ([`analysis/README_analysis.md`](analysis/README_analysis.md)) |
| `analysis/figs/` | The paper's figures (PNG and PDF) |
| `results/raw/` | One row per trial: `ladder.csv`, `prob.csv`, `long.csv`, `profiles.csv` |
| `results/scan_n*.csv` | Exhaustive scans of the estimator's `cost_gsa` over β |
| `results/predictions.csv` | Estimator output for the grid |
| `results/environment*.json` | Hardware, library versions and git commit recorded by each run |
| `PROTOCOL.md` | The experimental protocol and a dated record of deviations from it |
| `RUN_LONG.md`, `RUN_V5.md` | Instructions for the long-tour run and the basis-profile run |
| `PLAN.md` | Original execution plan and timing estimates |
| `MANIFEST.sha256` | SHA-256 checksums of the code, data and analysis outputs |
| `scripts/` | `fetch_strategies.sh` (BKZ strategies file) and `make_manifest.sh` |

The raw data, scans and figures are committed so that every result in the
paper can be regenerated without running BKZ.

---

## 2. Installation

Two environments are supported. The recorded runs used Option B for BKZ
(Python 3.14.4, fpylll 0.6.4) and SageMath 10.7 for the estimator.

### Option A: conda or mamba (single environment)

SageMath includes `fpylll` and is required by the estimator, so one
environment covers every step.

```bash
mamba create -n lwe -c conda-forge sage python=3.11 numpy pandas matplotlib scipy
mamba activate lwe

git clone https://github.com/malb/lattice-estimator.git    # into the repository root
git -C lattice-estimator checkout 53da598                   # version used in the paper
python3 src/strategies.py                                   # must print a path, not "none"
```

### Option B: pip only (no SageMath)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./scripts/fetch_strategies.sh                               # required; see below
```

Option B cannot run `src/estimate.sage`. Estimator predictions must then be
produced separately, with Option A, CoCalc or the `sagemath/sagemath` Docker
image, and copied to `results/predictions.csv`.

### Known installation issues

- **Missing strategies file.** The `fpylll` wheels on PyPI set
  `BKZ.DEFAULT_STRATEGY` to a path on the build machine
  (`/project/local/share/fplll/strategies/default.json`). BKZ with that
  strategy then fails with `RuntimeError: Cannot open strategies file`.
  `scripts/fetch_strategies.sh` downloads the file from the fplll repository
  into `data/`, and `src/strategies.py` locates it. Without the file, BKZ runs
  without preprocessing or pruning and is several times slower. Timings taken
  with and without it must not be compared; every run records which file it
  used in its `environment*.json`.
- **Missing `cysignals`.** `fpylll` fails at import without `cysignals`, which
  pip does not install automatically. It is listed in `requirements.txt`.

### Verification

Run these checks before any experiment:

```bash
python3 src/lwe_instance.py    # instance generation and self-check
python3 src/embedding.py       # must print "target in lattice: True"
python3 src/gsa_estimate.py    # must print "MATCH" for Kyber-512
python3 src/strategies.py      # must print a path, not "none"
```

`MATCH` confirms that the reimplementation of the GSA condition reproduces the
estimator's published Kyber-512 output (β = 406, d = 998, δ = 1.003941). If it
does not, no prediction in this repository should be trusted.

---

## 3. Reproducing the results

### 3.1 From the committed data (no BKZ)

Every table, figure and number in the paper is computed by the scripts in
`analysis/` from the committed CSV files. The order of the scripts and the
output of each are listed in
[`analysis/README_analysis.md`](analysis/README_analysis.md). Figures are
written to `analysis/figs/`.

To verify that the data and code match the published versions:

```bash
sha256sum -c MANIFEST.sha256
```

After any change to code or data, rebuild the manifest with
`scripts/make_manifest.sh` and commit it with the change.

### 3.2 Re-running the experiments

```bash
# Estimator predictions for the grid (requires SageMath)
sage src/estimate.sage                                       # -> results/predictions.csv

# Main sweep: success fraction against block size, 8 tours per call
python3 src/run_sweep.py prob   --ns 32 36 40 44 --trials 20 --jobs 8

# Block-size ladder and wall-clock scaling
python3 src/run_sweep.py ladder --ns 32 36 40 44 48 --trials 10 --jobs 8

# Sample-count trade-off
python3 src/run_sweep.py mtrade --ns 40 --trials 5 --jobs 8

# Tour-budget control: same instances, up to 32 tours per call (RUN_LONG.md)
python3 src/run_sweep.py long --ns 36 32 40 --trials 20 --jobs 7
python3 src/run_sweep.py long --ns 44 --trials 20 --jobs 7 --beta-above 6
python3 src/run_sweep.py long --ns 36 --trials 20 --jobs 7 --beta-window 20   # adds β = 14, 16

# Basis profiles after LLL and after 1, 2, 4, 8 and 16 tours (RUN_V5.md)
python3 src/run_sweep.py profile --ns 32 36 40 44 --trials 20 --jobs 7
```

All sweeps are resumable: re-running a command executes only the trials
missing from its CSV file. `--dry-run` prints the number of jobs without
running them. `--jobs` should be set to the number of physical cores minus one.

`src/analyze.py` produces an earlier quick-look summary (`figures/`, not
tracked, and `results/summary_*.csv`). Its GSA column is the global minimum of
`cost_gsa` from the scans, not the output of `primal_usvp` (see Section 5).

### 3.3 Replaying a single trial

Each row of the raw data records the instance seed and a fingerprint (a hash
of `A` and `b`), so any trial can be replayed:

```python
from lwe_instance import gen_lwe
from attack import solve_single

inst = gen_lwe(n=40, q=257, m=88, seed=2, dist_s="dg", param_s=4,
               dist_e="dg", param_e=4)
print(inst.fingerprint())                    # must match the CSV row
print(solve_single(inst, beta=46, tours=8))  # recorded outcome: solved at tour 5
```

---

## 4. Experimental design

### 4.1 Parameters

The instance family uses `q = 257`, with secret and error both drawn from a
discrete Gaussian distribution with `σ = 4`. The secret dimension `n` runs from
32 upward, and the number of samples `m` is the number the estimator's own
optimum uses.

- **Equal secret and error distributions.** The estimator's rescaling factor
  `ξ = σ_e/σ_s` is then exactly 1, so the embedding basis is integral and is
  the lattice the estimator models. Changing one distribution without the
  other breaks this correspondence (see `src/embedding.py`).
- **High noise rate.** The noise rate `α = σ√(2π)/q ≈ 0.039` is about 42 times
  that of Kyber-512 (`≈ 0.00092`). Instances of this size are only non-trivial
  at such a rate, so the results describe a regime different from deployed
  parameters.
- **Lower limit `n ≥ 32`.** Below this, the estimator's optimal `m` yields
  instances that are not well posed (Section 4.2), and the driver refuses them.

Run as a script, `src/gsa_estimate.py` prints a design table for choosing a
different `(q, σ)` family. Suitable families have a predicted β between 40 and
80, lattice dimension `d` up to about 180, and a gap of at most 0.5.

### 4.2 Safeguards

Two failure modes would invalidate the results without any visible error.
Both are guarded in the code.

- **Instances that are not uniquely decodable.** With `m ≈ n`, many secrets
  have a small error, LLL alone finds one of them, and every dimension appears
  broken. The guard is the ratio `‖t‖ / λ₁`, where `λ₁` is the Gaussian
  heuristic. `run_sweep.py` refuses any parameter point with a ratio above
  `GAP_LIMIT = 0.6` and reports the reason.
- **A permissive success test.** `embedding.py` counts a trial as solved only
  if the recovered secret equals the planted one, or if its error is no larger
  in ℓ₂ norm than the planted error. In every recorded success the recovered
  secret was the planted one (column `exact`).

---

## 5. Relation to the estimator

Three properties of the estimator govern how the comparison is made. Each can
be checked in its source code.

1. **Range of block sizes.** Under the GSA shape model, `primal_usvp` searches
   block sizes from 40 upward (`local_minimum(40, max_beta, ...)` in
   `PrimalUSVP.__call__`), so it cannot report β < 40. A measured β below 40
   where the prediction is below 40 is therefore outside the tool's intended
   range. In addition, at n = 32, 36 and 40 the search misses the cheapest
   point within its own range (paper, Section 5.7; lattice-estimator issue
   [#225](https://github.com/malb/lattice-estimator/issues/225)). All
   predictions used in this study are therefore the global minima of
   `cost_gsa` from exhaustive scans (`results/scan_n*.csv`), not the output of
   `primal_usvp`.
2. **Cost model.** The estimator's operation counts assume sieving (MATZOV,
   BDGL16, ADPS16 and others), whereas `fpylll` uses pruned enumeration at
   these block sizes. Wall-clock times are therefore not comparable with the
   estimator's `rop`; the quantity compared is the block size.
3. **Single BKZ call.** The estimator models one BKZ-β call, not a progressive
   ladder of block sizes. Progressive BKZ is stronger, so comparing it with the
   estimator's β would understate the estimate. Both modes are measured and
   reported separately.

---

## 6. Provenance and reproducibility

Each `results/environment*.json` file records the commit used by the run.

| Run | Date | Commit |
|---|---|---|
| Main sweep and ladder (`environment.json`) | 2026-09-23 | `0f3cce2`, with uncommitted changes |
| Long-tour run (`environment_long.json`) | 2026-10-16 | `ab5a2b9`, with uncommitted changes |
| Basis-profile run (`environment_profile.json`) | 2026-10-23 | `581bfca`, clean |

Commit `0f3cce2` is not part of this repository's history; the first commit
here is `99117b0` (2026-10-13).

The outcome of a BKZ call is reproducible for a given instance and block size.
The long-tour run repeats the 8-tour calls of the main sweep and agrees with it
on all 1,260 shared instance–block-size pairs (the 10 sweep calls that timed
out are excluded), and on the tour of success in 694 of the 708 solved cases.
Success probability is therefore reported as a fraction of instances, and
confidence intervals resample instances rather than calls.

The experiment follows `PROTOCOL.md`. Changes made after data collection began
are recorded there with dates and reasons.

---

## 7. Licence

The code is released under the MIT Licence; see [`LICENSE`](LICENSE). The
`lattice-estimator`, the fplll strategies file and the simulator of
Postlethwaite and Virdia are not included and remain under their authors'
terms.
