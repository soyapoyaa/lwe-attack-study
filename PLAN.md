# Seven-day execution plan

Written against measured timings from this repository on **one** core. Times
below are single-core; divide by your core count minus one.

| parameter point | one full β-search | 20×9 single-BKZ calls (`prob`) |
|---|---|---|
| n=32 (d=104) | ~3 s | ~3 min |
| n=36 (d=121) | ~7 s | ~6 min |
| n=40 (d=129) | ~30 s | ~30 min |
| n=44 (d=136) | ~3 min | ~2 h |
| n=48 (d=145) | ~10 min (est.) | ~7 h |
| n=52 (d=164) | ~40 min (est.) | ~1.5 days |

So: **n ≤ 48 is a comfortable week. n=52 is an overnight stretch goal. n≥56 is
not a one-week target — do not plan on it.**

The single most important scheduling fact: compute runs unattended, writing
does not. Launch every sweep before you start writing, never after.

---

## Day 1 — Environment and protocol freeze

**Goal: every self-test passes, and the protocol is written down and not
changed again.**

1. Install (README §1). Prefer conda/Sage — it gives `fpylll` and the estimator
   in one environment and avoids the pip strategies bug.
2. Run all four verification commands in README §1. Do not proceed until
   `embedding.py` prints `target in lattice: True` and `gsa_estimate.py` prints
   `MATCH`. The `MATCH` line is your evidence that you understand what the
   estimator claims; without it every later comparison is guesswork.
3. `git init`, commit. From here on, every sweep records the commit hash in
   `results/environment.json`. **Never run a sweep with uncommitted changes to
   `src/` —** `git_dirty: true` in that file makes the run unciteable.
4. Smoke test, ~1 minute:
   ```bash
   python3 src/run_sweep.py ladder --ns 32 --trials 2 --jobs 2
   python3 src/analyze.py
   ```
   Open `figures/beta_vs_n.png`. If a figure appears, the whole pipeline works.
5. **Write `PROTOCOL.md` now, before seeing results.** One page:
   - parameter family (q, σ, distributions, n range) and why;
   - `m` = the estimator's own optimal number of samples;
   - success = exact recovery of the planted secret;
   - well-posedness guard gap ≤ 0.6;
   - single-β mode is the estimator comparison; progressive is reported
     separately;
   - trial count per point;
   - the β grid and the number of BKZ tours.

   Committing this *before* looking at the data is what separates a measurement
   from a story fitted to noise. If you later change any of it, say so in the
   paper and say why.

**Deliverable:** passing self-tests, `PROTOCOL.md`, first commit.

---

## Day 2 — Predictions, then launch the long sweep

1. Clone `lattice-estimator` into the repo root, then:
   ```bash
   sage src/estimate.sage
   ```
   Check the last lines reproduce the estimator's own Kyber-512 docstring
   output. Commit `results/predictions.csv`.
2. Read the printed table. For each `n`, note `beta_gsa`, `beta_cn11`, `d_gsa`.
   If `beta_gsa` is 40 for several small `n`, that is the search floor, not a
   prediction — those points cannot be used for the headline claim.
3. **Launch the headline sweep and leave it running:**
   ```bash
   python3 src/run_sweep.py prob --ns 32 36 40 44 --trials 20 --jobs 8 \
       > logs/prob.log 2>&1 &
   ```
   Use `--dry-run` first to see the job count. On 8 cores this is roughly half
   a day for n ≤ 44.
4. While it runs, on the *same machine but a different core budget* — or just
   wait — do not start a second sweep. Contending sweeps corrupt your timing
   measurements, which are a reported quantity.

**Deliverable:** `results/predictions.csv` committed; `prob` sweep running.

---

## Day 3 — Ladder, trade-off, first look at the data

1. When `prob` finishes, launch:
   ```bash
   python3 src/run_sweep.py ladder --ns 32 36 40 44 48 --trials 10 --jobs 8
   python3 src/run_sweep.py mtrade --ns 40 --trials 5 --jobs 8
   ```
2. `python3 src/analyze.py`, then look hard at `figures/success_prob.png`.
   Three questions:
   - Does each curve rise from 0 to 1 across the β window? If a curve is flat
     at 0 or 1, the window is wrong — widen `--beta-window` and resume.
   - Is β₅₀ defined (a 0.5 crossing exists) for every n?
   - How far is β₅₀ from `beta_gsa` and from `beta_cn11`?
3. Decide the paper's claim from the data, in one sentence. Likely one of:
   - the estimator is conservative here, by a factor you can quote;
   - the estimator is accurate within the seed spread;
   - CN11 tracks the measurements and GSA does not.

   All three are perfectly publishable. Write the sentence down.
4. If n=48 is cheaper than feared, launch n=52 overnight as a stretch point.

**Deliverable:** all three experiments done for n ≤ 44, figures, and one
sentence stating the finding.

---

## Day 4 — Methods and Results (the data is in hand; write it)

Write these two sections in full. They are the ones that require the data and
nothing else, so they are the ones to write while the stretch runs finish.

**Methods** — in this order: the LWE instance; the Kannan embedding with the
explicit basis and why `(e, ξs, τ)` lies in it; ξ and τ chosen to match the
estimator; what success means; the well-posedness guard with the gap
definition; the β search; single vs progressive; trial counts; hardware and
versions from `results/environment.json`.

**Results** — one subsection per figure. For each: what was measured, what the
figure shows, and the number. State the seed spread every single time you quote
a β. A reader who sees "β=34 at n=40" and later learns the three seeds gave
34, 46, 46 will stop trusting the paper.

Also write the limitations paragraph today, while you are still close to the
mechanics:
- toy dimensions and an elevated noise rate `α ≈ 0.039` versus Kyber's
  `≈ 0.00092` — a different regime from deployed parameters;
- the estimator's GSA search floor at β=40;
- enumeration (fpylll) versus the sieving assumed by the cost models, so no
  wall-clock/`rop` comparison is made;
- pruning is randomised, so a single trial's β is not a threshold;
- a rounded Gaussian is used, not an exact discrete Gaussian;
- single-machine timings, one core per trial.

**Deliverable:** Methods and Results drafted; limitations paragraph written.

---

## Day 5 — Introduction, related work, discussion

1. **Related work, and be precise here.** Your experiment is close to
   Postlethwaite and Virdia, *On the Success Probability of Solving Unique SVP
   via BKZ* (PKC 2021), which studied exactly the tightness of the 2016
   estimate, and to Albrecht, Göpfert, Virdia, Wunderer and the
   Dachman-Soled–Ducas–Gong–Rossi line of work on refined predictions. Cite
   them, state what they did, and state what yours adds — the very small-`n`
   regime, the quantified single-versus-progressive gap, and the GSA-versus-CN11
   comparison on the same instances.

   Position the work as an independent replication plus a small extension.
   Reviewers in this area know these papers. A replication study is a
   legitimate and welcome SIURO contribution; an overclaimed novelty is the one
   thing that will sink the submission.

2. **Introduction:** LWE underpins ML-KEM/Kyber; concrete parameter choices rest
   on the estimator; the estimator rests on heuristics (GSA, the 2016 success
   condition); heuristics deserve measurement; measuring at deployed sizes is
   infeasible, so measure at toy sizes and report honestly what that does and
   does not tell you.

3. **Discussion:** why the direction of the discrepancy you found makes sense.
   If measured β is below predicted, candidate explanations to discuss: the
   2016 condition is sufficient rather than necessary; GSA misdescribes basis
   shape in small dimension (test this — compare against your CN11 column);
   pruning plus multiple tours effectively exceeds one clean BKZ-β call.

**Deliverable:** complete draft.

---

## Day 6 — Reproducibility pass

The highest-value day, and the one most often skipped.

1. **Fresh clone, clean environment, rerun a subset.**
   ```bash
   git clone <your repo> /tmp/check && cd /tmp/check
   # install per README §1
   python3 src/run_sweep.py ladder --ns 32 36 --trials 3 --jobs 4
   python3 src/analyze.py
   ```
   Do the n=32 and n=36 numbers land inside the range in your paper? If not,
   find out why before submitting, not after.
2. Check every number in the text against `results/summary_*.csv`. Every quoted
   figure should be traceable to a CSV column.
3. Confirm `results/environment.json` shows `git_dirty: false`.
4. Regenerate all figures from committed CSVs alone, with no BKZ calls. That is
   the property a reviewer will actually exercise.
5. Proofread the abstract last. It should contain the finding and its size, not
   a promise of one.

**Deliverable:** verified-reproducible repository; final PDF.

---

## Day 7 — Submit

1. **IACR ePrint.** Submit as a PDF at <https://eprint.iacr.org/submit>.
   Category: public-key cryptography / cryptanalysis. Include the repository
   URL in the abstract and in an artifact-availability paragraph. ePrint is not
   peer-reviewed and posting does not conflict with a later SIURO submission.
2. **SIURO.** Read the author guidelines at
   <https://www.siam.org/publications/siuro> before formatting — SIURO requires
   a faculty advisor's endorsement letter confirming the work is
   undergraduate-authored, so **ask your advisor for that on Day 1, not Day 7.**
   Use SIAM's LaTeX macros. SIURO expects the exposition to be readable by an
   undergraduate audience, so keep the lattice background self-contained.
3. Tag the repository (`git tag v1.0-preprint`) and archive it to Zenodo for a
   citable DOI.

---

## Cut scope in this order if you fall behind

1. Drop `mtrade` (the sample/dimension trade-off figure). Nice, not load-bearing.
2. Drop the largest `n`. A clean four-point ladder beats a ragged six-point one.
3. Drop the progressive-BKZ comparison and report single-β only. Halves the
   compute.
4. Reduce `prob` trials from 20 to 10 — but **never** below 10. The whole point
   is quantifying spread, and spread from 5 samples is not a measurement.

Do **not** cut, under any circumstances: the well-posedness guard, exact-secret
success checking, the seed spread in reported numbers, or the limitations
section.

---

## If the results look boring

"The estimator predicted β and we measured β" is a *result*, and a useful one —
it is the first independent check of that heuristic in this regime. Do not
manufacture a discrepancy by loosening the success test or changing parameters
until something disagrees. If you do change parameters mid-study, report every
family you ran, including the ones that agreed.
