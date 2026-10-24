# Basis-profile run and model comparison

## Purpose

These runs and analyses test explanations for success below the predicted block
size (Section 5.8 of the paper). Only the basis-profile run executes BKZ, and it
takes about 15–20 minutes.

| Step | Question | Environment | Time |
|---|---|---|---|
| 2 | What does the reduced basis look like after 1–16 tours? | fpylll | 15–20 min |
| 3a | Does the difficulty of an instance track the length of its planted vector? | fpylll | < 1 min |
| 3b | How do the Postlethwaite–Virdia simulator and the CN11 condition compare with the measurements? | fpylll | about 1 min |
| 3c | How does the measured basis compare with the GSA and CN11 at the index the success condition uses? | fpylll | about 1 min |
| 4 | Estimator predictions, including CN11 (`results/predictions.csv`) | SageMath | a few minutes |
| 5 | Does the search miss of issue #225 occur for Kyber? | SageMath | a few minutes |

The code for these runs is commit `581bfca`, which the profile run recorded with
a clean working tree (`results/environment_profile.json`). The CN11 setting in
`analysis/model_predictions.py` was corrected afterwards in `830091e`;
`model_predictions.json` was produced by the corrected script.

## 1. Simulator code

The simulator of Postlethwaite and Virdia has no licence file, so it is not
included in this repository. The analysis imports it from a separate clone:

```bash
git clone https://github.com/fvirdia/usvp-simulation ../usvp-simulation   # commit 333416d
```

## 2. Basis-profile run

From the repository root, in the environment with fpylll:

```bash
python3 src/run_sweep.py profile --ns 32 36 40 44 --trials 20 --jobs 7 --dry-run   # 160 jobs
python3 src/run_sweep.py profile --ns 32 36 40 44 --trials 20 --jobs 7
```

## 3. Analyses

```bash
cd analysis
python3 target_norm.py                                   | tee target_norm.log
python3 model_predictions.py --pv21 ../../usvp-simulation | tee model_predictions.log
python3 gs_profiles.py                                   | tee gs_profiles.log
cd ..
```

The first line of `gs_profiles.log` must state that every outcome agrees with
`prob.csv` within 8 tours.

## 4. Estimator predictions

In the SageMath environment:

```bash
sage src/estimate.sage | tee results/estimate.log         # writes results/predictions.csv
```

## 5. Kyber check for issue #225

```bash
sage src/check_kyber_usvp.sage | tee results/check_kyber_usvp.txt
```

The recorded output is `same point: True` for Kyber-512, Kyber-768 and
Kyber-1024: the search miss reported in issue #225 does not occur there.
