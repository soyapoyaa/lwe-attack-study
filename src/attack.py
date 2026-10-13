"""
The attack: progressive BKZ on the primal embedding lattice, instrumented.

Protocol (fixed once, for every trial):
  1. LLL-reduce the embedding basis.  Check for a solution.
  2. For beta in the schedule: run up to `tours` BKZ-2.0 tours at block size
     beta, checking for a solution after every tour.
  3. Stop at the first solution, at the beta budget, or at the wall-clock
     timeout.

Recorded per trial: the smallest beta that worked (`beta_success`), the time
spent up to and including that beta, the achieved root-Hermite factor, and the
log Gram-Schmidt profile. The block size -- not the wall-clock time -- is the
quantity comparable to the estimator, because fpylll uses enumeration with
pruning at these sizes whereas the estimator's cost models assume sieving.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from fpylll import BKZ, GSO, LLL, FPLLL, IntegerMatrix
from fpylll.algorithms.bkz2 import BKZReduction

from embedding import primal_basis, scan_basis
from strategies import bkz_param
from lwe_instance import LWEInstance


def float_type_for(d: int) -> str:
    """Enough precision without paying for it. Bump one level if fpylll warns."""
    if d <= 90:
        return "d"
    if d <= 180:
        return "ld"
    return "dd"


def beta_schedule(start: int = 2, stop: int = 60, step: int = 2) -> list[int]:
    """Block sizes to try, in increasing order. beta=2 is plain LLL."""
    return [b for b in range(start, stop + 1, step)]


def log_gs_profile(M: GSO.Mat) -> list[float]:
    """log2 ||b_i*|| for the current basis."""
    M.update_gso()
    return [0.5 * np.log2(M.get_r(i, i)) for i in range(M.d)]


def root_hermite(M: GSO.Mat, log_vol: float) -> float:
    """delta_0 = (||b_0|| / vol^(1/d))^(1/d), computed in logs."""
    M.update_gso()
    d = M.d
    log2_b0 = 0.5 * np.log2(M.get_r(0, 0))
    return float(2.0 ** ((log2_b0 - log_vol / d) / d))


@dataclass
class TrialResult:
    solved: bool = False
    beta_success: int | None = None
    tours_at_success: int | None = None
    time_to_success: float = float("nan")
    time_total: float = 0.0
    time_lll: float = 0.0
    exact_secret: bool | None = None
    err_inf: int | None = None
    timed_out: bool = False
    delta_0: float = float("nan")
    b0_norm: float = float("nan")
    per_beta: list[dict] = field(default_factory=list)
    gs_profile: list[float] = field(default_factory=list)
    error: str | None = None


def run_attack(
    inst: LWEInstance,
    betas: list[int] | None = None,
    tours: int = 8,
    timeout: float = 600.0,
    verbose: bool = False,
) -> tuple[TrialResult, dict]:
    """Attack one instance. Returns (result, embedding metadata)."""
    B, target, meta = primal_basis(inst)
    d = meta["d"]
    betas = betas or beta_schedule(stop=min(60, d - 1))
    betas = [b for b in betas if b <= d]
    log_vol = (
        (d - inst.n - 1) * np.log2(inst.q)
        + inst.n * np.log2(meta["xi"])
        + np.log2(meta["tau"])
    )

    res = TrialResult()
    t_start = time.perf_counter()

    try:
        ft = float_type_for(d)
        FPLLL.set_precision(160)

        t0 = time.perf_counter()
        LLL.reduction(B)
        res.time_lll = time.perf_counter() - t0

        M = GSO.Mat(B, float_type=ft, update=True)
        bkz = BKZReduction(M)

        hit = scan_basis(inst, B, meta)
        if hit is not None:
            res.solved, res.beta_success, res.tours_at_success = True, 2, 0
            res.time_to_success = res.time_lll
            res.exact_secret = hit["exact"]
            res.err_inf = hit["err_inf"]
        res.per_beta.append(
            {"beta": 2, "tours": 0, "time": res.time_lll,
             "delta_0": root_hermite(M, log_vol), "solved": res.solved}
        )
        if verbose:
            print(f"  LLL      t={res.time_lll:7.2f}s  "
                  f"delta_0={res.per_beta[-1]['delta_0']:.5f}  solved={res.solved}")

        for beta in betas:
            if res.solved:
                break
            if beta <= 2:
                continue
            par = bkz_param(beta)
            t_beta0 = time.perf_counter()
            for tour in range(1, tours + 1):
                if time.perf_counter() - t_start > timeout:
                    res.timed_out = True
                    break
                bkz.tour(par)
                hit = scan_basis(inst, B, meta)
                if hit is not None:
                    res.solved = True
                    res.beta_success = beta
                    res.tours_at_success = tour
                    res.time_to_success = time.perf_counter() - t_start
                    res.exact_secret = hit["exact"]
                    res.err_inf = hit["err_inf"]
                    break
            dt = time.perf_counter() - t_beta0
            rec = {"beta": beta, "tours": tour, "time": dt,
                   "delta_0": root_hermite(M, log_vol), "solved": res.solved}
            res.per_beta.append(rec)
            if verbose:
                print(f"  beta={beta:3d}  t={dt:7.2f}s  "
                      f"delta_0={rec['delta_0']:.5f}  solved={res.solved}")
            if res.timed_out:
                break

        M.update_gso()
        res.delta_0 = root_hermite(M, log_vol)
        res.b0_norm = float(np.sqrt(M.get_r(0, 0)))
        res.gs_profile = log_gs_profile(M)

    except Exception as exc:  # keep a sweep alive when one trial explodes
        res.error = f"{type(exc).__name__}: {exc}"

    res.time_total = time.perf_counter() - t_start
    return res, meta


# --------------------------------------------------------------------------
# single-block-size mode: the estimator's actual model
# --------------------------------------------------------------------------


def solve_single(inst: LWEInstance, beta: int, tours: int = 8,
                 timeout: float = 600.0) -> dict:
    """LLL, then BKZ at ONE block size, from a fresh basis.

    This is what the estimator models: a single BKZ-beta call on a freshly
    reduced basis. `run_attack` in progressive mode instead climbs a ladder of
    block sizes and is therefore strictly stronger -- comparing progressive
    results to the estimator's beta understates the estimator.
    """
    B, _, meta = primal_basis(inst)
    d = meta["d"]
    out = {"beta": beta, "d": d, "solved": False, "tours_used": 0,
           "time": 0.0, "timed_out": False, "exact": None, "error": None}
    t0 = time.perf_counter()
    try:
        FPLLL.set_precision(160)
        LLL.reduction(B)
        M = GSO.Mat(B, float_type=float_type_for(d), update=True)
        bkz = BKZReduction(M)
        hit = scan_basis(inst, B, meta)
        if hit is not None:
            out.update(solved=True, exact=hit["exact"], tours_used=0)
        else:
            par = bkz_param(min(beta, d))
            for tour in range(1, tours + 1):
                if time.perf_counter() - t0 > timeout:
                    out["timed_out"] = True
                    break
                bkz.tour(par)
                out["tours_used"] = tour
                hit = scan_basis(inst, B, meta)
                if hit is not None:
                    out.update(solved=True, exact=hit["exact"])
                    break
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    out["time"] = time.perf_counter() - t0
    return out


def find_beta_single(inst: LWEInstance, betas: list[int], tours: int = 8,
                     timeout: float = 600.0, beta_hint: int | None = None,
                     verbose: bool = False) -> dict:
    """Smallest beta in `betas` for which a single BKZ-beta call succeeds.

    Exponential search outward from `beta_hint` (pass the estimator's predicted
    beta) to find a bracket, then binary search inside it. Starting from the
    hint rather than from max(betas) matters a lot: one wasted probe at a large
    block size can cost more than the entire rest of the search.

    Success is *usually* monotone in beta but not guaranteed -- pruning is
    probabilistic and BKZ is randomised -- so the return value is "the smallest
    beta located by this search". That is why the headline result uses success
    *probability* over many instances (experiment `prob`) rather than one
    threshold per instance.
    """
    betas = sorted(betas)
    trace: list[dict] = []

    def attempt(idx: int) -> dict:
        r = solve_single(inst, betas[idx], tours=tours, timeout=timeout)
        trace.append(r)
        if verbose:
            print(f"    beta={betas[idx]:3d} solved={r['solved']} t={r['time']:6.1f}s")
        return r

    def fail(reason=None):
        return {"beta_single": None, "solved": False, "trace": trace,
                "time_total": sum(x["time"] for x in trace),
                "probes": len(trace),
                "timed_out": any(x["timed_out"] for x in trace),
                "error": reason}

    if beta_hint is None:
        start = len(betas) // 2
    else:
        start = min(range(len(betas)), key=lambda i: abs(betas[i] - beta_hint))

    r = attempt(start)

    if r["solved"]:
        lo_fail, hi_ok, step = -1, start, 1
        while hi_ok - step > lo_fail:
            idx = hi_ok - step
            rr = attempt(idx)
            if rr["solved"]:
                hi_ok, step = idx, 1
            else:
                lo_fail = idx
                break
        else:
            pass
    else:
        lo_fail, hi_ok, step, idx = start, None, 1, start
        while hi_ok is None:
            idx = idx + step
            if idx > len(betas) - 1:
                if lo_fail >= len(betas) - 1:
                    return fail()
                idx = len(betas) - 1
            rr = attempt(idx)
            if rr["solved"]:
                hi_ok = idx
            elif idx == len(betas) - 1:
                return fail()
            else:
                lo_fail = idx
                step *= 2

    while hi_ok - lo_fail > 1:
        mid = (lo_fail + hi_ok) // 2
        rr = attempt(mid)
        if rr["solved"]:
            hi_ok = mid
        else:
            lo_fail = mid

    best = next(x for x in trace if x["beta"] == betas[hi_ok] and x["solved"])
    return {"beta_single": betas[hi_ok], "solved": True, "trace": trace,
            "exact": best["exact"], "time_at_beta": best["time"],
            "time_total": sum(x["time"] for x in trace), "probes": len(trace),
            "timed_out": any(x["timed_out"] for x in trace), "error": None}


if __name__ == "__main__":
    from lwe_instance import gen_lwe

    from gsa_estimate import primal_usvp

    for n in (32, 40):
        p = primal_usvp(n=n, q=257, sigma_s=4, sigma_e=4, m=8 * n)
        inst = gen_lwe(n=n, q=257, m=p.m_opt, seed=1,
                       dist_s="dg", param_s=4, dist_e="dg", param_e=4)
        print(f"n={n} d={p.d} predicted beta={p.beta}")
        res, meta = run_attack(inst, betas=beta_schedule(2, 60, 2),
                               tours=4, timeout=600)
        print(f"  progressive: beta*={res.beta_success} t={res.time_to_success:.1f}s")
        out = find_beta_single(inst, beta_schedule(4, 60, 2), tours=8,
                               timeout=600, verbose=True)
        print(f"  single-beta: beta*={out['beta_single']}\n")
