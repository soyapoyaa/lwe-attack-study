"""
Primal (Kannan) embedding of an LWE instance into a uSVP lattice instance.

With A in Z^{m x n}, b = A s + e mod q, the row lattice of

        [  q*I_m     0        0  ]   m rows
    B = [  -A^T    xi*I_n     0  ]   n rows
        [   b^T      0       tau ]   1 row

contains the target vector

    t = ( e | xi*s | tau ),      dim d = m + n + 1

because the integer combination with coefficients (k, s, 1) gives
q*k - A s + b = q*k + e = e for the right k.

This is exactly the lattice the lattice-estimator's `primal_usvp` models, with
    xi  = sigma_e / sigma_s        (rescales the secret block)
    tau = sigma_e                  (embedding coordinate)
so keeping Xs = Xe gives xi = 1 and an integral basis. See
`estimator/lwe_primal.py::PrimalUSVP._xi_factor`.

The attack succeeds when lattice reduction finds +-t. Since ||t|| is much
shorter than the Gaussian-heuristic length of the lattice, t is (with
overwhelming probability) the unique shortest vector -- hence "unique SVP".
"""

from __future__ import annotations

import numpy as np
from fpylll import IntegerMatrix

from lwe_instance import LWEInstance, centered


def xi_tau(inst: LWEInstance) -> tuple[int, int]:
    """Integral (xi, tau) following the estimator's choice.

    We require xi to be a positive integer; that holds whenever
    sigma_e / sigma_s is an integer, in particular for Xs = Xe (xi = 1).
    """
    ratio = inst.sigma_e / inst.sigma_s
    xi = ratio if ratio > 1 else 1.0
    if abs(xi - round(xi)) > 1e-9:
        raise ValueError(
            f"xi = sigma_e/sigma_s = {xi:.4f} is not an integer; the integral "
            "embedding would be an approximation. Use Xs = Xe (xi = 1) or "
            "choose stddevs with an integer ratio."
        )
    tau = max(1, int(round(inst.sigma_e)))
    return int(round(xi)), tau


def primal_basis(inst: LWEInstance) -> tuple[IntegerMatrix, np.ndarray, dict]:
    """Return (basis as fpylll IntegerMatrix, target vector, metadata)."""
    n, m, q = inst.n, inst.m, inst.q
    xi, tau = xi_tau(inst)
    d = m + n + 1

    M = np.zeros((d, d), dtype=object)

    # block 1: q * I_m
    for i in range(m):
        M[i, i] = int(q)

    # block 2: [ -A^T | xi*I_n | 0 ]   (row j of this block is -A[:, j])
    Ac = centered(inst.A, q)                      # smaller entries -> better numerics
    for j in range(n):
        for i in range(m):
            M[m + j, i] = int(-Ac[i, j])
        M[m + j, m + j] = int(xi)

    # block 3: [ b | 0 | tau ]
    bc = centered(inst.b, q)
    for i in range(m):
        M[d - 1, i] = int(bc[i])
    M[d - 1, d - 1] = int(tau)

    B = IntegerMatrix.from_matrix([[int(v) for v in row] for row in M])

    target = np.concatenate([inst.e, xi * inst.s, [tau]]).astype(np.int64)
    meta = {
        "d": d, "xi": xi, "tau": tau,
        "target_norm": float(np.linalg.norm(target)),
        # sqrt(m*sigma_e^2 + n*(xi*sigma_s)^2 + tau^2): expected target length
        "target_norm_expected": float(
            np.sqrt(m * inst.sigma_e**2 + n * (xi * inst.sigma_s) ** 2 + tau**2)
        ),
        # log2 of the Gaussian-heuristic length lambda_1 ~ sqrt(d/(2 pi e)) * det^(1/d)
        "gh_norm": gaussian_heuristic_norm(d, n, q, xi, tau),
    }
    return B, target, meta


def gaussian_heuristic_norm(d: int, n: int, q: int, xi: int, tau: int) -> float:
    """sqrt(d/(2 pi e)) * vol^(1/d) with vol = q^(d-n-1) * xi^n * tau."""
    log_vol = (d - n - 1) * np.log(q) + n * np.log(xi) + np.log(tau)
    return float(np.sqrt(d / (2 * np.pi * np.e)) * np.exp(log_vol / d))


def verify_candidate(inst: LWEInstance, row, meta: dict):
    """Is this reduced-basis row a solution?

    Accepts +-t. Recovers a candidate secret from the middle block and requires
    either exact recovery of the planted secret, or a candidate whose error is
    no larger (in l2) than the planted error -- i.e. an equally good or better
    solution to the LWE instance.

    Do NOT loosen this to "error is smallish". For m close to n the instance has
    astronomically many secrets with small error, so a loose bound reports
    success on instances that are not uniquely decodable, and LLL alone will
    "break" them. The well-posedness guard in run_sweep.py (usvp_gap) is the
    other half of that defence.
    """
    n, m, q = inst.n, inst.m, inst.q
    xi, tau = meta["xi"], meta["tau"]
    v = np.array([int(x) for x in row], dtype=np.int64)

    if abs(v[-1]) != tau:
        return None
    if v[-1] < 0:
        v = -v

    s_cand, r = np.divmod(v[m : m + n], xi)
    if np.any(r != 0):
        return None

    exact = bool(np.array_equal(s_cand, inst.s))
    e_cand = centered(inst.b - inst.A @ s_cand, q)
    if not exact:
        if np.linalg.norm(e_cand) > np.linalg.norm(inst.e) * (1 + 1e-9):
            return None

    return {
        "secret": s_cand,
        "exact": exact,
        "err_inf": int(np.abs(e_cand).max()),
        "err_l2": float(np.linalg.norm(e_cand)),
    }


def scan_basis(inst: LWEInstance, B: IntegerMatrix, meta: dict):
    """First solving row of B, or None."""
    for i in range(B.nrows):
        hit = verify_candidate(inst, B[i], meta)
        if hit is not None:
            hit["row"] = i
            return hit
    return None


if __name__ == "__main__":
    from lwe_instance import gen_lwe

    inst = gen_lwe(n=10, q=257, m=20, seed=7)
    B, t, meta = primal_basis(inst)
    print("dim d =", meta["d"], "| xi =", meta["xi"], "| tau =", meta["tau"])
    print("||t|| =", round(meta["target_norm"], 2),
          "| E||t|| =", round(meta["target_norm_expected"], 2),
          "| GH =", round(meta["gh_norm"], 2))
    print("gap ||t||/GH =", round(meta["target_norm"] / meta["gh_norm"], 4))

    # the target must actually lie in the lattice: check via the known coeffs
    Bm = np.array([[int(B[i, j]) for j in range(B.ncols)] for i in range(B.nrows)],
                  dtype=object)
    Ac, bc = centered(inst.A, inst.q), centered(inst.b, inst.q)
    k = (Ac @ inst.s - bc + inst.e) // inst.q
    coeffs = np.concatenate([k, inst.s, [1]]).astype(object)
    got = np.array([int(x) for x in coeffs @ Bm])
    print("target in lattice:", np.array_equal(got, t))
