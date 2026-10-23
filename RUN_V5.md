# Runs and analyses for v5

Four analyses for the v5 draft. Only step 3 needs real BKZ runs, and it takes about 15–20 minutes. Everything else takes a few minutes.

| Step | What it answers | Environment | Time |
|---|---|---|---|
| 3 | Basis-shape run: what the reduced basis actually looks like after 1–16 tours | lwe-venv | ~15–20 min |
| 4a | Does instance difficulty track the planted vector's length? | lwe-venv | < 1 min |
| 4b | Postlethwaite–Virdia simulator and CN11 condition against your measurements | lwe-venv | ~1 min |
| 4c | Measured basis shape against the GSA and CN11 at the index the success condition uses | lwe-venv | ~1 min |
| 5 | Estimator predictions incl. CN11 (`results/predictions.csv`, missing so far) | sage | a few min |
| 6 | Optional: does issue #225 affect Kyber? | sage | a few min |

## 1. Apply the changes and commit

```bash
SRC=/mnt/c/Users/maina/OneDrive/Desktop/lwe-v5-scripts      # wherever you unzipped the zip
REPO=/mnt/c/Users/maina/OneDrive/Desktop/quant/lwe-attack-study
cd "$REPO"
git status                                   # should list nothing modified
git apply --check "$SRC/v5-scripts.patch" && git apply "$SRC/v5-scripts.patch"
git add -A && git commit -m "Add profile run, model comparison and target-length analyses"
```

## 2. Get the Postlethwaite–Virdia simulator code

Their repository has no licence file, so it isn't copied into yours. The analysis imports it from your clone.

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

If it prints `same point: True` for all three, you can write in the issue that the second-basin miss does not affect Kyber.

## 7. Commit and push

```bash
conda deactivate
git add -A results analysis
git commit -m "v5 runs: basis profiles, model comparison, target lengths, estimator predictions"
git push
```

Your repository is public, so once it's pushed I can read the results directly. You don't need to attach anything.
