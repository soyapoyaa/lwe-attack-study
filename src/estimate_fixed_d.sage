import csv, os, sys
sys.path.insert(0, os.path.join(os.getcwd(), "lattice-estimator"))
from estimator import LWE, ND
from estimator.lwe_primal import PrimalUSVP
from sage.all import oo

Q = 257
SIGMA = 4.0

used = {}
with open("results/raw/ladder.csv") as fh:
    for row in csv.DictReader(fh):
        n = int(row["n"]); m = int(row["m"])
        used.setdefault(n, m)  # first row per n -- m is constant per n here

rows = []
for n, m in sorted(used.items()):
    d = m + n + 1
    X = ND.DiscreteGaussian(SIGMA)
    p = LWE.Parameters(n=n, q=Q, Xs=X, Xe=X, m=m, tag="fixed-d-n%d" % n)

    best = None
    for beta in range(4, 82, 2):
        c = PrimalUSVP.cost_gsa(beta=beta, params=p, d=d)
        if c["rop"] != oo:
            best = (beta, c)
            break

    if best is None:
        print("n=%d m=%d d=%d: no beta in [4,80] solves it under GSA" % (n, m, d))
        rows.append({"n": n, "m_used": m, "d_used": d, "beta_gsa_fixed_d": ""})
        continue

    beta, c = best
    print("n=%d m=%d d=%d -> beta_gsa (fixed d) = %d" % (n, m, d, beta))
    rows.append({"n": n, "m_used": m, "d_used": d, "beta_gsa_fixed_d": beta})

out = os.path.join("results", "predictions_fixed_d.csv")
with open(out, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=["n", "m_used", "d_used", "beta_gsa_fixed_d"])
    w.writeheader()
    for r in rows:
        w.writerow(r)
print("wrote", out)
