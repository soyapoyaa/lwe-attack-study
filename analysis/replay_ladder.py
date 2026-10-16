"""Exact replay of the ladder search on the sweep's own outcomes (Section 5.3).

Replaces anchor_sim.py. The simulation drew every probe independently from one
fitted curve. The data show that is the wrong model: a single BKZ-beta call on
a given instance gives the same outcome every time it is run (check 1 below), so
the only randomness is *between* instances. We can therefore run the ladder's
search procedure directly on each instance's recorded sweep outcomes.

  1. Determinism check. For every ladder instance (seeds 1-10, n = 32-44) the
     ladder's search must have seen a success at beta_single and, unless
     beta_single is the bottom of the grid, a failure at beta_single - 2. Both
     points were also attacked in the sweep. Count agreements.
  2. Replay. Run the ladder's search (same code path as attack.find_beta_single)
     on the sweep outcomes of seeds 1-10 and compare with the observed ladder,
     instance by instance.
  3. Counterfactual starts. Replay on all 20 sweep instances, starting the
     search at the prediction, 8 below it, and 8 above it.

Probes outside the sweep's grid (beta_GSA +- 16) are reported, never guessed
silently. Timed-out sweep trials (n = 44, beta 64-70) are treated as unknown;
the script reports whether any replay touched one.

Run from analysis/:  python3 replay_ladder.py   -> replay_ladder.json
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "results", "raw")
P = pd.read_csv(os.path.join(RAW, "prob.csv"))
L = pd.read_csv(os.path.join(RAW, "ladder.csv"))
GSA = {32: 20, 36: 34, 40: 44, 44: 54}
GRID = list(range(4, 82, 2))          # the ladder's grid (run_sweep defaults)


class OffGrid(Exception):
    pass


def search(ok, hint):
    """Same control flow as attack.find_beta_single, on a predicate ok(index)."""
    start = min(range(len(GRID)), key=lambda i: abs(GRID[i] - hint))
    if ok(start):
        lo_fail, hi_ok, step = -1, start, 1
        while hi_ok - step > lo_fail:
            idx = hi_ok - step
            if ok(idx):
                hi_ok, step = idx, 1
            else:
                lo_fail = idx
                break
    else:
        lo_fail, hi_ok, step, idx = start, None, 1, start
        while hi_ok is None:
            idx = idx + step
            if idx > len(GRID) - 1:
                if lo_fail >= len(GRID) - 1:
                    return None
                idx = len(GRID) - 1
            if ok(idx):
                hi_ok = idx
            elif idx == len(GRID) - 1:
                return None
            else:
                lo_fail = idx
                step *= 2
    while hi_ok - lo_fail > 1:
        mid = (lo_fail + hi_ok) // 2
        if ok(mid):
            hi_ok = mid
        else:
            lo_fail = mid
    return GRID[hi_ok]


def outcomes(n, seed):
    g = P[(P.n == n) & (P.seed == seed)].set_index("beta")
    return {int(b): (None if r.timed_out else bool(r.solved)) for b, r in g.iterrows()}


def replay(n, seed, hint):
    out = outcomes(n, seed)
    probed = []

    def ok(i):
        b = GRID[i]
        probed.append(b)
        if b not in out:
            raise OffGrid(b)
        if out[b] is None:
            raise OffGrid(f"censored@{b}")
        return out[b]

    try:
        return search(ok, hint), probed, None
    except OffGrid as exc:
        return None, probed, str(exc)


def main():
    res = {"determinism": {}, "replay_seeds_1_10": {}, "counterfactual": {}}

    # 1. determinism
    tot_s = ok_s = tot_f = ok_f = 0
    for _, r in L[L.n.isin(GSA)].iterrows():
        out = outcomes(int(r.n), int(r.seed))
        b = int(r.beta_single)
        if b in out and out[b] is not None:
            tot_s += 1
            ok_s += out[b] is True
        if b > GRID[0] and b - 2 in out and out[b - 2] is not None:
            tot_f += 1
            ok_f += out[b - 2] is False
    res["determinism"] = dict(success_points=[ok_s, tot_s], failure_points=[ok_f, tot_f])
    print(f"1. Same (instance, beta) in ladder and sweep: successes agree {ok_s}/{tot_s}, "
          f"failures agree {ok_f}/{tot_f}")

    # 2. replay on the ladder's own instances
    print("\n2. Replay of the ladder search on sweep outcomes, seeds 1-10")
    for n in GSA:
        obs = L[L.n == n].set_index("seed").beta_single.astype(int)
        rep, notes = {}, []
        for sd in range(1, 11):
            v, _, err = replay(n, sd, GSA[n])
            rep[sd] = v
            if err:
                notes.append((sd, err))
        match = sum(rep[s] == obs[s] for s in obs.index)
        res["replay_seeds_1_10"][n] = dict(
            observed=[int(obs[s]) for s in range(1, 11)],
            replayed=[rep[s] for s in range(1, 11)],
            per_instance_match=int(match), observed_median=float(obs.median()),
            replayed_median=float(np.median([v for v in rep.values() if v is not None])),
            off_grid=notes)
        print(f"   n={n}: per-instance match {match}/10 | median observed "
              f"{obs.median():.0f}, replayed {res['replay_seeds_1_10'][n]['replayed_median']:.0f}"
              + (f" | off-grid: {notes}" if notes else ""))

    # 3. counterfactual starting points on all 20 instances
    print("\n3. Same search on all 20 sweep instances, from different starting points")
    for n in GSA:
        row = {}
        for label, shift in (("pred-8", -8), ("pred", 0), ("pred+8", 8)):
            vals, notes, lower = [], [], []
            for sd in range(1, 21):
                v, probed, err = replay(n, sd, GSA[n] + shift)
                if err:
                    # The search left the sweep grid while stepping UP after failures,
                    # so the answer is above every failed probe: a lower bound only.
                    out = outcomes(n, sd)
                    fails = [b for b in probed if out.get(b) is False]
                    notes.append((sd, err))
                    lower.append(max(fails) + 2)
                else:
                    vals.append(v)
            # median of all 20, with each off-grid instance placed at +infinity
            # (valid because its true value exceeds its lower bound, which is
            # above the median in every case; checked below)
            full = sorted(vals + [np.inf] * len(lower))
            med_all = float(np.median(full[:20])) if len(full) == 20 else float("nan")
            assert all(lb > med_all for lb in lower), "lower bound not above median"
            row[label] = dict(median=med_all, n_ok=len(vals), off_grid=notes,
                              off_grid_lower_bounds=lower, values=vals)
        res["counterfactual"][n] = row
        print(f"   n={n} (pred {GSA[n]}): " + " | ".join(
            f"start {k}: median {v['median']:.0f} ({v['n_ok']}/20)"
            + (f" off-grid {v['off_grid']}" if v["off_grid"] else "")
            for k, v in row.items()))

    json.dump(res, open(os.path.join(HERE, "replay_ladder.json"), "w"), indent=1, default=int)


if __name__ == "__main__":
    main()
