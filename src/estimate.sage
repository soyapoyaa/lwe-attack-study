## -*- coding: utf-8 -*-
"""
Official lattice-estimator predictions for the experiment grid.

Run under SageMath, from the repository root, with lattice-estimator cloned
into ./lattice-estimator (see README):

    sage src/estimate.sage

Writes results/predictions.csv with, for every n in the grid:

  beta_gsa,  d_gsa    primal uSVP, GSA basis-shape model  (the "2016 estimate")
  beta_cn11, d_cn11   primal uSVP, CN11 BKZ simulator     (more accurate, esp.
                                                           in small dimension)
  beta_zgsa, d_zgsa   primal uSVP, ZGSA (GSA + q-vectors)
  beta_bdd,  d_bdd    primal BDD (uSVP + a final CVP call)
  beta_dual, d_dual   dual attack
  rop_*               log2 operation counts under several cost models

Three properties of the estimator, all verifiable in its source:

  1. Under the GSA shape model, `primal_usvp` searches block sizes starting at
     beta = 40 (`local_minimum(40, max_beta, ...)` in PrimalUSVP.__call__), so
     it can never report beta < 40. Any grid point whose true block size is
     below 40 is therefore outside the tool's intended range, and a measured
     beta < 40 there is not evidence against the estimator. Separately, at
     n = 32, 36 and 40 the search misses the cheapest point within its own
     range (lattice-estimator issue #225); the study therefore uses the global
     minima of cost_gsa from exhaustive scans (results/scan_n*.csv).
  2. Cost models (MATZOV, BDGL16, ADPS16, ...) assume *sieving*. fpylll uses
     enumeration with pruning at these block sizes. So rop is not comparable to
     our wall-clock times; the comparable quantity is beta.
  3. The estimator reports the block size at which a *single* BKZ-beta call is
     predicted to succeed. Compare it against our single-beta measurements, not
     against progressive-BKZ measurements.
"""

import csv
import os
import sys

sys.path.insert(0, os.path.join(os.getcwd(), "lattice-estimator"))

from estimator import *            # noqa: F401,F403
from estimator import LWE, ND, RC, Simulator

# ---------------------------------------------------------------- grid
Q = 257
SIGMA = 4.0
NS = [32, 36, 40, 44, 48, 52, 56]
SAMPLE_MULT = 8                    # m = SAMPLE_MULT * n samples available

OUT = os.path.join("results", "predictions.csv")


def params_for(n):
    X = ND.DiscreteGaussian(SIGMA)
    return LWE.Parameters(n=n, q=Q, Xs=X, Xe=X, m=SAMPLE_MULT * n,
                          tag="toy-n%d" % n)


def grab(cost, keys=("beta", "d", "rop", "delta")):
    out = {}
    for k in keys:
        try:
            v = cost[k]
        except (KeyError, TypeError):
            v = None
        out[k] = v
    return out


def log2_or_blank(v):
    from sage.all import log, RR, oo
    if v is None or v == oo:
        return ""
    try:
        return float(RR(log(v, 2)))
    except Exception:
        return ""


rows = []
for n in NS:
    p = params_for(n)
    print("=" * 66)
    print(p)

    rec = {"n": n, "q": Q, "sigma": SIGMA, "m_available": SAMPLE_MULT * n}

    variants = [
        ("gsa",  lambda p: LWE.primal_usvp(p, red_shape_model="gsa")),
        ("cn11", lambda p: LWE.primal_usvp(p, red_shape_model=Simulator.CN11)),
        ("zgsa", lambda p: LWE.primal_usvp(p, red_shape_model=Simulator.ZGSA)),
        ("bdd",  lambda p: LWE.primal_bdd(p)),
        ("dual", lambda p: LWE.dual(p)),
    ]
    for tag, fn in variants:
        try:
            c = fn(p)
            print("  %-5s %s" % (tag, repr(c)))
            g = grab(c)
            rec["beta_%s" % tag] = g["beta"]
            rec["d_%s" % tag] = g["d"]
            rec["delta_%s" % tag] = float(g["delta"]) if g["delta"] else ""
            rec["log2rop_%s" % tag] = log2_or_blank(g["rop"])
        except Exception as exc:
            print("  %-5s FAILED: %s: %s" % (tag, type(exc).__name__, exc))
            for suf in ("beta", "d", "delta", "log2rop"):
                rec["%s_%s" % (suf, tag)] = ""

    # same attack, alternative reduction cost models: changes rop, not beta
    for name, model in (("adps16", RC.ADPS16), ("bdgl16", RC.BDGL16),
                        ("matzov", RC.MATZOV), ("chengu12", RC.CheNgu12)):
        try:
            c = LWE.primal_usvp(p, red_shape_model="gsa", red_cost_model=model)
            rec["log2rop_gsa_%s" % name] = log2_or_blank(c["rop"])
        except Exception:
            rec["log2rop_gsa_%s" % name] = ""

    # m the estimator would actually use, for the GSA variant
    if rec.get("d_gsa"):
        rec["m_opt_gsa"] = int(rec["d_gsa"]) - n - 1
    rows.append(rec)

keys = sorted({k for r in rows for k in r})
keys = ["n", "q", "sigma", "m_available", "m_opt_gsa"] + \
       [k for k in keys if k not in ("n", "q", "sigma", "m_available", "m_opt_gsa")]

os.makedirs("results", exist_ok=True)
with open(OUT, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in keys})

print("=" * 66)
print("wrote %s (%d rows)" % (OUT, len(rows)))
print("\nSanity check -- this must reproduce the estimator's own docstring:")
print("  LWE.primal_usvp(schemes.Kyber512) ->", repr(LWE.primal_usvp(schemes.Kyber512)))
print("  expected: rop: 2^143.8, delta: 1.003941, beta: 406, d: 998")
