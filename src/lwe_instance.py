"""
Generation of toy LWE instances.

Convention (matches the lattice-estimator):
    A in Z_q^{m x n},  s <- Xs^n,  e <- Xe^m,  b = A s + e  mod q.

Secret and error are drawn from the *same* centered-binomial distribution by
default, which makes the estimator's rescaling factor xi = sigma_e/sigma_s
exactly 1 and therefore keeps our embedding lattice integral and identical to
the one the estimator models. Do not change this without reading
`embedding.py` first.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field, asdict

import numpy as np

# --------------------------------------------------------------------------
# noise distributions
# --------------------------------------------------------------------------


def centered_binomial(eta: int, size, rng: np.random.Generator) -> np.ndarray:
    """Centered binomial with parameter eta: variance eta/2, stddev sqrt(eta/2).

    eta = 2  ->  sigma = 1.0        (used by Kyber for the error term)
    eta = 3  ->  sigma = 1.2247...  (used by Kyber-512 for the secret)
    """
    a = rng.integers(0, 2, size=(eta,) + (size if isinstance(size, tuple) else (size,)))
    b = rng.integers(0, 2, size=(eta,) + (size if isinstance(size, tuple) else (size,)))
    return (a.sum(axis=0) - b.sum(axis=0)).astype(np.int64)


def discrete_gaussian(sigma: float, size, rng: np.random.Generator) -> np.ndarray:
    """Rounded Gaussian. Close enough to a discrete Gaussian for sigma >= 1.

    The paper describes it as a rounded Gaussian, since it is not an exact
    discrete Gaussian.
    """
    n = size if isinstance(size, tuple) else (size,)
    return np.rint(rng.normal(0.0, sigma, size=n)).astype(np.int64)


def stddev_of(kind: str, param: float) -> float:
    if kind == "cb":
        return float(np.sqrt(param / 2.0))
    if kind == "dg":
        return float(param)
    raise ValueError(f"unknown distribution {kind!r}")


def sample(kind: str, param: float, size, rng) -> np.ndarray:
    if kind == "cb":
        return centered_binomial(int(param), size, rng)
    if kind == "dg":
        return discrete_gaussian(float(param), size, rng)
    raise ValueError(f"unknown distribution {kind!r}")


# --------------------------------------------------------------------------
# instances
# --------------------------------------------------------------------------


def centered(x: np.ndarray, q: int) -> np.ndarray:
    """Representatives in (-q/2, q/2]."""
    y = x % q
    y[y > q // 2] -= q
    return y


@dataclass
class LWEInstance:
    n: int
    q: int
    m: int
    dist_s: str
    param_s: float
    dist_e: str
    param_e: float
    seed: int
    A: np.ndarray = field(repr=False)
    b: np.ndarray = field(repr=False)
    s: np.ndarray = field(repr=False)
    e: np.ndarray = field(repr=False)

    @property
    def sigma_s(self) -> float:
        return stddev_of(self.dist_s, self.param_s)

    @property
    def sigma_e(self) -> float:
        return stddev_of(self.dist_e, self.param_e)

    def truncate(self, m_use: int) -> "LWEInstance":
        """Use only the first m_use samples (the estimator optimises over this)."""
        if m_use > self.m:
            raise ValueError(f"asked for {m_use} samples, instance only has {self.m}")
        d = asdict(self)
        d.update(
            m=m_use,
            A=self.A[:m_use].copy(),
            b=self.b[:m_use].copy(),
            e=self.e[:m_use].copy(),
            s=self.s.copy(),
        )
        return LWEInstance(**d)

    def check(self) -> bool:
        """Sanity: b == A s + e mod q."""
        lhs = self.b % self.q
        rhs = (self.A @ self.s + self.e) % self.q
        return bool(np.array_equal(lhs, rhs))

    def fingerprint(self) -> str:
        """Stable hash of the instance, so a result row can be traced back."""
        h = hashlib.sha256()
        for arr in (self.A, self.b):
            h.update(np.ascontiguousarray(arr % self.q, dtype=np.int64).tobytes())
        return h.hexdigest()[:16]


def gen_lwe(
    n: int,
    q: int,
    m: int,
    seed: int,
    dist_s: str = "cb",
    param_s: float = 2,
    dist_e: str = "cb",
    param_e: float = 2,
) -> LWEInstance:
    """Fresh LWE instance. `seed` fully determines it."""
    rng = np.random.default_rng([seed, n, q, m])
    A = rng.integers(0, q, size=(m, n), dtype=np.int64)
    s = sample(dist_s, param_s, n, rng)
    e = sample(dist_e, param_e, m, rng)
    b = (A @ s + e) % q
    inst = LWEInstance(
        n=n, q=q, m=m,
        dist_s=dist_s, param_s=param_s, dist_e=dist_e, param_e=param_e,
        seed=seed, A=A, b=b, s=s, e=e,
    )
    assert inst.check(), "instance generation is broken"
    return inst


if __name__ == "__main__":
    inst = gen_lwe(n=20, q=257, m=40, seed=1)
    print(inst)
    print("valid:", inst.check(), "| fingerprint:", inst.fingerprint())
    print("sigma_s:", inst.sigma_s, "sigma_e:", inst.sigma_e)
    print("s:", inst.s)
    print("|e|_inf:", int(np.abs(inst.e).max()))
    t = inst.truncate(30)
    print("truncated:", t.m, t.check())
