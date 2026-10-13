"""
A pure-Python copy of the lattice-estimator's GSA-based primal-uSVP estimate.

Why this file exists
--------------------
1. It runs without SageMath, so parameter selection, well-posedness checks and
   the experiment driver need no heavy dependency.
2. It is an *independent reimplementation* of the same formulas, so agreement
   with `LWE.primal_usvp(..., red_shape_model="gsa")` is a cheap correctness
   check on our understanding of what the estimator claims. The paper's
   headline comparison must still use the real estimator (`estimate.sage`) --
   this is scaffolding, not a substitute.

Formulas follow estimator/lwe_primal.py (PrimalUSVP.cost_gsa, _solve_for_d)
and estimator/reduction.py (ReductionCost._delta) as of 2024-2025.

Validation target, straight from the estimator's own docstring:
    LWE.primal_usvp(schemes.Kyber512)
    -> rop: 2^143.8, delta: 1.003941, beta: 406, d: 998
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# --------------------------------------------------------------------------
# root-Hermite factor delta(beta)
# --------------------------------------------------------------------------

_SMALL_DELTA = (
    (2, 1.02190), (5, 1.01862), (10, 1.01616), (15, 1.01485),
    (20, 1.01420), (25, 1.01342), (28, 1.01331), (40, 1.01295),
)


def delta_f(beta: float) -> float:
    """Root-Hermite factor achieved by BKZ-beta (estimator convention)."""
    beta = round(beta)
    if beta <= 2:
        return 1.0219
    if beta < 40:
        for i in range(1, len(_SMALL_DELTA)):
            if _SMALL_DELTA[i][0] > beta:
                return _SMALL_DELTA[i - 1][1]
    if beta == 40:
        return _SMALL_DELTA[-1][1]
    return (beta / (2 * math.pi * math.e) * (math.pi * beta) ** (1 / beta)) ** (
        1 / (2 * (beta - 1))
    )


# --------------------------------------------------------------------------
# core-SVP style cost models (log2 of gate/operation count)
# --------------------------------------------------------------------------


def log2_rop_adps16(beta: int, d: int) -> float:
    """Core-SVP sieving, classical: 2^(0.292 beta)."""
    return 0.2920 * beta


def log2_rop_bdgl16(beta: int, d: int) -> float:
    """BDGL16 sieving with the estimator's polynomial factors."""
    return 0.2920 * beta + math.log2(8 * d) + 16.4


# --------------------------------------------------------------------------
# the uSVP predicate
# --------------------------------------------------------------------------


def xi_factor(sigma_s: float, sigma_e: float) -> float:
    return sigma_e / sigma_s if sigma_s < sigma_e else 1.0


def usvp_predicate(beta: int, d: int, n: int, q: int,
                   sigma_s: float, sigma_e: float, tau: float | None = None) -> bool:
    """True if BKZ-beta in dimension d is predicted to solve the uSVP instance.

    lhs = log ||pi_{d-beta}(t)||   (expected projected target length)
    rhs = log ||b*_{d-beta}||      (GSA-predicted basis profile)
    """
    delta = delta_f(beta)
    xi = xi_factor(sigma_s, sigma_e)
    tau = sigma_e if tau is None else tau
    lhs = math.log(math.sqrt(sigma_e**2 * (beta - 1) + tau**2))
    rhs = (
        math.log(delta) * (2 * beta - d - 1)
        + (math.log(tau) + math.log(xi) * n + math.log(q) * (d - n - 1)) / d
    )
    return lhs <= rhs


def solve_for_d(beta: int, n: int, q: int, sigma_s: float, sigma_e: float,
                m: float = math.inf, tau: float | None = None) -> float:
    """Smallest lattice dimension d satisfying the predicate, or inf."""
    delta = delta_f(beta)
    xi = xi_factor(sigma_s, sigma_e)
    tau = sigma_e if tau is None else tau

    a = -math.log(delta)
    C = math.log(sigma_e**2 * (beta - 1) + tau**2) / 2.0
    c = math.log(tau) + n * math.log(xi) - (n + 1) * math.log(q)
    b = math.log(delta) * (2 * beta - 1) + math.log(q) - C

    if a * n * n + b * n + c >= 0:
        return n
    disc = b * b - 4 * a * c
    if disc < 0:
        return math.inf
    d1 = (-b + math.sqrt(disc)) / (2 * a)
    d2 = (-b - math.sqrt(disc)) / (2 * a)
    if a > 0:
        return min(m, math.ceil(d2))
    if n <= d1:
        return min(m, math.ceil(d1))
    return math.inf


@dataclass
class Prediction:
    beta: int
    d: int
    delta: float
    m_opt: int
    log2_rop_adps16: float
    log2_rop_bdgl16: float
    feasible: bool = True


def primal_usvp(n: int, q: int, sigma_s: float, sigma_e: float,
                m: float = math.inf, tau: float | None = None,
                max_beta: int = 1754, min_beta: int = 2) -> Prediction:
    """Smallest feasible beta and the matching minimal d (GSA shape model).

    `m` is the number of LWE *samples* available. Following Bai-Galbraith (and
    PrimalUSVP.__call__), when sigma_s <= sigma_e the secret block contributes
    n extra usable dimensions, so the lattice dimension is capped at
    d <= m + n + 1; otherwise at d <= m + 1.

    `min_beta`: set to 40 to mirror the real estimator, whose GSA search starts
    at beta = 40 and therefore never reports a smaller block size. Left at 2
    here so the function is usable for designing toy parameters.
    """
    if m == math.inf:
        m_cap = math.inf
    else:
        m_lat = m + n if sigma_s <= sigma_e else m
        m_cap = m_lat + 1
    for beta in range(max(2, min_beta), max_beta + 1):
        d = solve_for_d(beta, n, q, sigma_s, sigma_e, m=m_cap, tau=tau)
        if d == math.inf or d > m_cap:
            continue
        d = int(max(d, beta))
        if d == beta and d < m_cap:
            d += 1
        if not usvp_predicate(beta, d, n, q, sigma_s, sigma_e, tau=tau):
            continue
        return Prediction(
            beta=beta, d=d, delta=delta_f(beta), m_opt=d - n - 1,
            log2_rop_adps16=log2_rop_adps16(beta, d),
            log2_rop_bdgl16=log2_rop_bdgl16(beta, d),
        )
    return Prediction(beta=-1, d=-1, delta=float("nan"), m_opt=-1,
                      log2_rop_adps16=float("inf"),
                      log2_rop_bdgl16=float("inf"), feasible=False)


# --------------------------------------------------------------------------
# well-posedness of the experiment
# --------------------------------------------------------------------------


def gaussian_heuristic_norm(d: int, n: int, q: int, xi: float, tau: float) -> float:
    log_vol = (d - n - 1) * math.log(q) + n * math.log(xi) + math.log(tau)
    return math.sqrt(d / (2 * math.pi * math.e)) * math.exp(log_vol / d)


def usvp_gap(n: int, q: int, m: int, sigma_s: float, sigma_e: float,
             tau: float | None = None) -> float:
    """||t|| / lambda_1(GH).

    Must be well below 1, otherwise the planted vector is not the unique
    shortest vector, the instance has many valid secrets, and LLL alone will
    "solve" it -- a meaningless measurement. Aim for <= 0.5.
    """
    tau = sigma_e if tau is None else tau
    xi = xi_factor(sigma_s, sigma_e)
    d = m + n + 1
    t = math.sqrt(m * sigma_e**2 + n * (xi * sigma_s) ** 2 + tau**2)
    return t / gaussian_heuristic_norm(d, n, q, xi, tau)


if __name__ == "__main__":
    print("== validation against the estimator's published numbers ==")
    cases = [
        ("Kyber512", dict(n=512, q=3329, sigma_s=1.2247, sigma_e=1.2247, m=512),
         dict(beta=406, d=998, delta=1.003941)),
        ("Kyber768", dict(n=768, q=3329, sigma_s=1.0, sigma_e=1.0, m=768), None),
        ("Kyber1024", dict(n=1024, q=3329, sigma_s=1.0, sigma_e=1.0, m=1024), None),
    ]
    for name, kw, want in cases:
        p = primal_usvp(**kw)
        line = (f"{name:10s} beta={p.beta:5d}  d={p.d:5d}  delta={p.delta:.6f}  "
                f"log2 rop (BDGL16)={p.log2_rop_bdgl16:6.1f}")
        if want:
            ok = (p.beta == want["beta"] and p.d == want["d"]
                  and abs(p.delta - want["delta"]) < 1e-6)
            line += f"   expected beta={want['beta']} d={want['d']} -> {'MATCH' if ok else 'MISMATCH'}"
        print(line)

    print("\n== toy sweep: q=257, Xs=Xe=DG(sigma), m = 8n samples ==")
    hdr = f"{'n':>4} {'sig':>4} {'beta':>5} {'b>=40':>6} {'d':>5} {'m_opt':>6} {'gap':>7}"
    print(hdr)
    for sigma in (2.0, 4.0, 8.0):
        for n in range(20, 73, 8):
            p = primal_usvp(n=n, q=257, sigma_s=sigma, sigma_e=sigma, m=8 * n)
            p40 = primal_usvp(n=n, q=257, sigma_s=sigma, sigma_e=sigma, m=8 * n,
                              min_beta=40)
            g = usvp_gap(n, 257, p.m_opt, sigma, sigma) if p.feasible else float("nan")
            print(f"{n:4d} {sigma:4.0f} {p.beta:5d} {p40.beta:6d} {p.d:5d} "
                  f"{p.m_opt:6d} {g:7.3f}")
