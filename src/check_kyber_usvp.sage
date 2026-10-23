# Does primal_usvp's GSA search find the cheapest beta for real parameter sets?
# (Supports lattice-estimator issue #225: is the second-basin miss limited to small n?)
#
# Run from the repository root, under SageMath, with ./lattice-estimator present:
#     sage src/check_kyber_usvp.sage | tee results/check_kyber_usvp.txt
import os, sys
sys.path.insert(0, os.path.join(os.getcwd(), "lattice-estimator"))
from estimator import LWE, schemes
from estimator.lwe_primal import PrimalUSVP

for name in ("Kyber512", "Kyber768", "Kyber1024"):
    p = getattr(schemes, name).normalize()
    m = p.m + p.n if p.Xs <= p.Xe else p.m
    ret = LWE.primal_usvp(p, red_shape_model="gsa")
    top = min(int(m), 1200)
    best = min((PrimalUSVP.cost_gsa(beta=b, params=p, m=m) for b in range(40, top)),
               key=lambda c: c["rop"])
    print(f"{name}: primal_usvp beta={ret['beta']} d={ret['d']} log2rop={float(log(ret['rop'], 2)):.2f} | "
          f"cheapest scanned beta={best['beta']} d={best['d']} log2rop={float(log(best['rop'], 2)):.2f} | "
          f"same point: {ret['beta'] == best['beta']}")
