# Basis-profile run and model comparison

The runs and analyses behind Section 5.8 of the paper. Only step 3 needs real BKZ runs, and it takes about 15–20 minutes. Everything else takes a few minutes.

| Step | What it answers | Environment | Time |
|---|---|---|---|
| 3 | Basis-shape run: what the reduced basis actually looks like after 1–16 tours | lwe-venv | ~15–20 min |
| 4a | Does instance difficulty track the planted vector's length? | lwe-venv | < 1 min |
| 4b | Postlethwaite–Virdia simulator and CN11 condition against the measurements | lwe-venv | ~1 min |
| 4c | Measured basis shape against the GSA and CN11 at the index the success condition uses | lwe-venv | ~1 min |
| 5 | Estimator predictions incl. CN11 (`results/predictions.csv`) | sage | a few min |
| 6 | Optional: does issue #225 affect Kyber? | sage | a few min |

## 1. Code

The code for these runs is in commit `581bfca`; the profile run recorded that commit, clean (`results/environment_profile.json`). The CN11 setting in `analysis/model_predictions.py` was corrected afterwards in `830091e`, and `model_predictions.json` comes from the corrected script.

## 2. Get the Postlethwaite–Virdia simulator code

Their repository has no licence file, so it is not copied into this one. The analysis imports it from your clone.

```bash
git clone https://github.com/fvirdia/usvp-simulation ~/usvp-simulation
```

## 3. Basis-shape run (lwe-venv)

```bash
conda deactivate 2>/dev/null; source ~/lwe-venv/bin/activate
python3 -c "import fpylll, scipy, pandas, matplotlib; print('ok')" || pip install matplotlib
python3 src/run_sweep.py profile --ns 32 36 40 44 --trials 20 --jobs 7 --dry-run   # expect 160 jobs
python3 src/run_sweep.py profile --ns 32 36 40 44 --trials 20 --jobs 7
```

## 4. Analyses (lwe-venv, from analysis/)

```bash
cd analysis
python3 target_norm.py                             | tee target_norm.log
python3 model_predictions.py --pv21 ~/usvp-simulation | tee model_predictions.log
python3 gs_profiles.py                             | tee gs_profiles.log
cd ..
```

Check the first line of `gs_profiles.log`: it should say every outcome agrees with `prob.csv` within 8 tours.

## 5. Estimator predictions (Sage)

```bash
deactivate; conda activate sage
sage src/estimate.sage | tee results/estimate.log            # writes results/predictions.csv
```

## 6. Optional: Kyber check for issue #225 (Sage)

```bash
sage src/check_kyber_usvp.sage | tee results/check_kyber_usvp.txt
```

It printed `same point: True` for Kyber-512, -768 and -1024: the search miss of issue #225 does not occur there.
