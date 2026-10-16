"""
Experiment driver. Resumable, parallel, one CSV row per trial.

Experiments
-----------
prob    Success probability vs block size. For each (n, beta), attack N fresh
        instances with a single BKZ-beta call. Gives p(beta) and hence
        beta_50 -- the observable that is actually comparable to the
        estimator's predicted beta. This is the headline experiment.

ladder  Smallest solving beta per instance, by bisection, for each n. Both
        single-beta and progressive mode. Cheaper than `prob`, and gives the
        beta-vs-n ladder plus wall-clock scaling.

long    Convergence control for `prob`. Same instances (seeds), same lattice,
        one BKZ-beta call per (n, beta, seed), but up to --tours tours
        (default 32) and a record of when fpylll's AutoAbort / clean-tour
        tests would have stopped the call. Writes results/raw/long.csv and
        results/environment_long.json; never touches prob.csv. Every fixed
        budget up to --tours is recovered from one row (solved_tour <= T).

mtrade  Sample/dimension trade-off: fix n, vary the number of samples used
        (hence the lattice dimension d), find the solving beta for each.
        Shows the optimisation over d that the estimator performs internally.

Usage
-----
    python3 src/run_sweep.py prob   --ns 32 36 40 --trials 20 --jobs 4
    python3 src/run_sweep.py ladder --ns 28 32 36 40 44 --trials 5 --jobs 4
    python3 src/run_sweep.py mtrade --ns 40 --trials 3 --jobs 4
    python3 src/run_sweep.py long   --ns 36 32 40 44 --trials 20 --jobs 7

Re-running the same command appends only the missing trials, so the sweep can
be interrupted and resumed. Nothing is ever overwritten.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from attack import (beta_schedule, find_beta_single, run_attack, solve_long,
                    solve_single)
from embedding import xi_tau
from gsa_estimate import primal_usvp, usvp_gap
from lwe_instance import gen_lwe
from strategies import describe as strategies_describe

REPO = Path(__file__).resolve().parent.parent
RAW = REPO / "results" / "raw"

# Default parameter family. Chosen so that the estimator's predicted block size
# spans roughly 20-80 while the lattice dimension stays under ~170, i.e. the
# experiment is feasible on a laptop and the uSVP instance stays well-posed
# (usvp_gap well below 1). See PLAN.md for how these were selected.
DEFAULTS = dict(q=257, dist="dg", sigma=4.0)

GAP_LIMIT = 0.6  # refuse to run instances that are not honest uSVP instances


# --------------------------------------------------------------------------
# bookkeeping
# --------------------------------------------------------------------------


def environment() -> dict:
    def _cpu():
        try:
            for line in open("/proc/cpuinfo"):
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
        except OSError:
            pass
        return platform.processor() or "unknown"

    import fpylll

    return {
        "recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "cpu": _cpu(),
        "cores": os.cpu_count(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "fpylll": fpylll.__version__,
        "bkz_strategies": strategies_describe(),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
    }


def _git(*args) -> str:
    try:
        return subprocess.run(["git", *args], cwd=REPO, capture_output=True,
                              text=True, timeout=5).stdout.strip()
    except Exception:
        return ""


def append_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists()
    keys = list(rows[0].keys())
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        if new:
            w.writeheader()
        for r in rows:
            w.writerow(r)


def done_keys(path: Path, key_fields: list[str]) -> set[tuple]:
    if not path.exists():
        return set()
    with open(path) as fh:
        return {tuple(row[k] for k in key_fields) for row in csv.DictReader(fh)}


def _wellposed(n: int, m: int, q: int, sigma: float) -> bool:
    """Print and skip parameter points that are not honest uSVP instances."""
    try:
        guard(n, m, q, sigma)
        return True
    except ValueError as exc:
        print(f"SKIP n={n}: {exc}")
        return False


def make_instance(n: int, m: int, seed: int, q: int, dist: str, sigma: float):
    return gen_lwe(n=n, q=q, m=m, seed=seed,
                   dist_s=dist, param_s=sigma, dist_e=dist, param_e=sigma)


def guard(n: int, m: int, q: int, sigma: float) -> float:
    """Well-posedness check. Raises if the planted vector is not unique-shortest."""
    gap = usvp_gap(n, q, m, sigma, sigma)
    if gap > GAP_LIMIT:
        raise ValueError(
            f"n={n}, m={m}: usvp gap ||t||/GH = {gap:.2f} > {GAP_LIMIT}. The "
            "planted vector is not the unique shortest vector, so the instance "
            "has many valid secrets and 'success' would be meaningless. "
            "Increase m, increase q, or decrease sigma."
        )
    return gap


# --------------------------------------------------------------------------
# one trial per experiment
# --------------------------------------------------------------------------


def trial_prob(n, beta, seed, q, dist, sigma, m, tours, timeout) -> dict:
    inst = make_instance(n, m, seed, q, dist, sigma)
    xi, tau = xi_tau(inst)
    out = solve_single(inst, beta, tours=tours, timeout=timeout)
    return {
        "exp": "prob", "n": n, "q": q, "m": m, "d": out["d"], "sigma": sigma,
        "dist": dist, "xi": xi, "tau": tau, "seed": seed, "beta": beta,
        "solved": int(out["solved"]), "exact": out["exact"],
        "tours_used": out["tours_used"], "time": round(out["time"], 4),
        "timed_out": int(out["timed_out"]),
        "gap": round(usvp_gap(n, q, m, sigma, sigma), 4),
        "fingerprint": inst.fingerprint(), "error": out["error"] or "",
    }


def trial_long(n, beta, seed, q, dist, sigma, m, tours_max, timeout) -> dict:
    inst = make_instance(n, m, seed, q, dist, sigma)
    xi, tau = xi_tau(inst)
    out = solve_long(inst, beta, tours=tours_max, timeout=timeout)
    blank = lambda v: "" if v is None else v
    return {
        "exp": "long", "n": n, "q": q, "m": m, "d": out["d"], "sigma": sigma,
        "dist": dist, "xi": xi, "tau": tau, "seed": seed, "beta": beta,
        "tours_max": tours_max, "solved": int(out["solved"]), "exact": out["exact"],
        "solved_tour": blank(out["solved_tour"]), "tours_run": out["tours_run"],
        "abort_tour": blank(out["abort_tour"]), "clean_tour": blank(out["clean_tour"]),
        "time": round(out["time"], 4), "timed_out": int(out["timed_out"]),
        "gap": round(usvp_gap(n, q, m, sigma, sigma), 4),
        "fingerprint": inst.fingerprint(), "error": out["error"] or "",
    }


def trial_ladder(n, seed, q, dist, sigma, m, betas, tours, timeout,
                 beta_hint=None) -> dict:
    inst = make_instance(n, m, seed, q, dist, sigma)
    xi, tau = xi_tau(inst)

    single = find_beta_single(inst, betas, tours=tours, timeout=timeout,
                              beta_hint=beta_hint)
    prog, meta = run_attack(inst, betas=betas, tours=4, timeout=timeout)

    return {
        "exp": "ladder", "n": n, "q": q, "m": m, "d": meta["d"], "sigma": sigma,
        "dist": dist, "xi": xi, "tau": tau, "seed": seed,
        "beta_single": single.get("beta_single"),
        "beta_prog": prog.beta_success,
        "single_solved": int(single["solved"]),
        "single_probes": single.get("probes"),
        "prog_solved": int(prog.solved),
        "time_single_at_beta": round(single.get("time_at_beta") or 0.0, 4),
        "time_single_total": round(single.get("time_total") or 0.0, 4),
        "time_prog_to_success": round(prog.time_to_success, 4)
        if prog.solved else "",
        "time_prog_total": round(prog.time_total, 4),
        "time_lll": round(prog.time_lll, 4),
        "delta_0_final": round(prog.delta_0, 6),
        "exact": prog.exact_secret,
        "timed_out": int(bool(single.get("timed_out")) or prog.timed_out),
        "gap": round(usvp_gap(n, q, m, sigma, sigma), 4),
        "target_norm": round(meta["target_norm"], 3),
        "gh_norm": round(meta["gh_norm"], 3),
        "fingerprint": inst.fingerprint(),
        "error": (single.get("error") or "") or (prog.error or ""),
    }


def trial_mtrade(n, m, seed, q, dist, sigma, betas, tours, timeout,
                 beta_hint=None) -> dict:
    inst = make_instance(n, m, seed, q, dist, sigma)
    single = find_beta_single(inst, betas, tours=tours, timeout=timeout,
                              beta_hint=beta_hint)
    return {
        "exp": "mtrade", "n": n, "q": q, "m": m, "d": m + n + 1, "sigma": sigma,
        "dist": dist, "seed": seed,
        "beta_single": single.get("beta_single"),
        "solved": int(single["solved"]),
        "time_total": round(single.get("time_total") or 0.0, 4),
        "gap": round(usvp_gap(n, q, m, sigma, sigma), 4),
        "timed_out": int(bool(single.get("timed_out"))),
        "fingerprint": inst.fingerprint(), "error": single.get("error") or "",
    }


# --------------------------------------------------------------------------
# job construction
# --------------------------------------------------------------------------


def build_jobs(args) -> tuple[Path, list, list[str], callable]:
    q, dist, sigma = args.q, args.dist, args.sigma
    betas = beta_schedule(args.beta_min, args.beta_max, args.beta_step)
    jobs = []

    if args.exp == "prob":
        path = RAW / "prob.csv"
        keys = ["n", "beta", "seed", "m"]
        for n in args.ns:
            pred = primal_usvp(n=n, q=q, sigma_s=sigma, sigma_e=sigma, m=8 * n)
            m = args.m if args.m else pred.m_opt
            if not _wellposed(n, m, q, sigma):
                continue
            lo = max(args.beta_min, pred.beta - args.beta_window)
            hi = min(args.beta_max, pred.beta + args.beta_window)
            for beta in range(lo, hi + 1, args.beta_step):
                for seed in range(1, args.trials + 1):
                    jobs.append((trial_prob,
                                 dict(n=n, beta=beta, seed=seed, q=q, dist=dist,
                                      sigma=sigma, m=m, tours=args.tours,
                                      timeout=args.timeout)))
        return path, jobs, keys, None

    if args.exp == "long":
        path = RAW / "long.csv"
        keys = ["n", "beta", "seed", "m", "tours_max"]
        for n in args.ns:
            pred = primal_usvp(n=n, q=q, sigma_s=sigma, sigma_e=sigma, m=8 * n)
            m = args.m if args.m else pred.m_opt
            if not _wellposed(n, m, q, sigma):
                continue
            lo = max(args.beta_min, pred.beta - args.beta_window)
            hi = min(args.beta_max, pred.beta + args.beta_above)
            # same even grid as `prob`, so every row has a partner in prob.csv
            for beta in range(lo, hi + 1, args.beta_step):
                for seed in range(1, args.trials + 1):
                    jobs.append((trial_long,
                                 dict(n=n, beta=beta, seed=seed, q=q, dist=dist,
                                      sigma=sigma, m=m, tours_max=args.tours,
                                      timeout=args.timeout)))
        return path, jobs, keys, None

    if args.exp == "ladder":
        path = RAW / "ladder.csv"
        keys = ["n", "seed", "m"]
        for n in args.ns:
            pred = primal_usvp(n=n, q=q, sigma_s=sigma, sigma_e=sigma, m=8 * n)
            m = args.m if args.m else pred.m_opt
            if not _wellposed(n, m, q, sigma):
                continue
            for seed in range(1, args.trials + 1):
                jobs.append((trial_ladder,
                             dict(n=n, seed=seed, q=q, dist=dist, sigma=sigma,
                                  m=m, betas=betas, tours=args.tours,
                                  timeout=args.timeout, beta_hint=pred.beta)))
        return path, jobs, keys, None

    if args.exp == "mtrade":
        path = RAW / "mtrade.csv"
        keys = ["n", "m", "seed"]
        for n in args.ns:
            pred = primal_usvp(n=n, q=q, sigma_s=sigma, sigma_e=sigma, m=8 * n)
            centre = pred.m_opt
            ms = sorted({max(n + 5, int(centre * f)) for f in
                         (0.6, 0.7, 0.8, 0.9, 1.0, 1.15, 1.3, 1.5, 1.8)})
            for m in ms:
                if usvp_gap(n, q, m, sigma, sigma) > GAP_LIMIT:
                    continue
                for seed in range(1, args.trials + 1):
                    jobs.append((trial_mtrade,
                                 dict(n=n, m=m, seed=seed, q=q, dist=dist,
                                      sigma=sigma, betas=betas,
                                      tours=args.tours, timeout=args.timeout,
                                      beta_hint=pred.beta)))
        return path, jobs, keys, None

    raise ValueError(args.exp)


def _call(fn_kw):
    fn, kw = fn_kw
    try:
        return fn(**kw)
    except Exception as exc:
        return {"exp": "FAILED", "error": f"{type(exc).__name__}: {exc}", **kw}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("exp", choices=["prob", "ladder", "mtrade", "long"])
    ap.add_argument("--ns", type=int, nargs="+", required=True)
    ap.add_argument("--trials", type=int, default=10)
    ap.add_argument("--q", type=int, default=DEFAULTS["q"])
    ap.add_argument("--sigma", type=float, default=DEFAULTS["sigma"])
    ap.add_argument("--dist", default=DEFAULTS["dist"], choices=["dg", "cb"])
    ap.add_argument("--m", type=int, default=None,
                    help="samples to use; default = estimator-optimal m")
    ap.add_argument("--beta-min", type=int, default=4)
    ap.add_argument("--beta-max", type=int, default=80)
    ap.add_argument("--beta-step", type=int, default=2)
    ap.add_argument("--beta-window", type=int, default=16,
                    help="prob mode: scan predicted beta +- this; "
                         "long mode: this far below the prediction")
    ap.add_argument("--beta-above", type=int, default=16,
                    help="long mode: scan up to this far above the prediction "
                         "(default 16 = the prob grid; use 6 for n = 44, whose "
                         "high-beta trials are very slow)")
    ap.add_argument("--tours", type=int, default=None,
                    help="tours per BKZ call (default 8; 32 for `long`)")
    ap.add_argument("--timeout", type=float, default=None,
                    help="seconds per call, checked at tour boundaries "
                         "(default 900; 7200 for `long`)")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    if args.tours is None:
        args.tours = 32 if args.exp == "long" else 8
    if args.timeout is None:
        args.timeout = 7200.0 if args.exp == "long" else 900.0

    path, jobs, keys, _ = build_jobs(args)

    already = done_keys(path, keys)
    todo = [j for j in jobs
            if tuple(str(j[1][k]) for k in keys) not in already]

    print(f"experiment : {args.exp}")
    print(f"output     : {path}")
    print(f"parameters : q={args.q} dist={args.dist} sigma={args.sigma} "
          f"ns={args.ns} trials={args.trials} tours={args.tours} "
          f"timeout={args.timeout:.0f}s")
    print(f"jobs       : {len(jobs)} total, {len(already)} already done, "
          f"{len(todo)} to run, {args.jobs} workers")
    if args.dry_run:
        return

    # `long` keeps its own record so the original sweep's environment.json survives
    env_name = "environment_long.json" if args.exp == "long" else "environment.json"
    env_path = REPO / "results" / env_name
    env_path.parent.mkdir(parents=True, exist_ok=True)
    env_path.write_text(json.dumps(environment(), indent=2))
    print(f"environment: {env_path}\n")

    t0 = time.time()
    buf: list[dict] = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(_call, j): j for j in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            row = fut.result()
            buf.append(row)
            if len(buf) >= 1:   # flush every row: an interrupted sweep loses nothing
                append_rows(path, buf)
                buf = []
            el = time.time() - t0
            eta = el / i * (len(todo) - i)
            tag = " ".join(f"{k}={row.get(k)}" for k in
                           ("n", "beta", "m", "seed") if k in row)
            got = row.get("solved", row.get("beta_single", "?"))
            print(f"[{i:4d}/{len(todo)}] {tag:38s} -> {got}   "
                  f"elapsed {el/60:5.1f}m  eta {eta/60:5.1f}m", flush=True)
    append_rows(path, buf)
    print(f"\ndone in {(time.time()-t0)/60:.1f} min -> {path}")


if __name__ == "__main__":
    main()
