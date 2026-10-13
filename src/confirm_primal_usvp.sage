import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), "lattice-estimator"))
from estimator import LWE, ND
from estimator.lwe_primal import PrimalUSVP
from sage.all import log, RR

Q, SIGMA = 257, 4.0
GLOBAL = {32: 20, 36: 34, 40: 44, 44: 54, 48: 63}   # from results/scan_n*.csv

def l2(c):
    return float(RR(log(c["rop"], 2)))

for n in (32, 36, 40, 44, 48):
    m = 8 * n
    X = ND.DiscreteGaussian(SIGMA)
    p = LWE.Parameters(n=n, q=Q, Xs=X, Xe=X, m=m)
    r = LWE.primal_usvp(p, red_shape_model="gsa")
    print("n=%d  primal_usvp ->" % n, r)
    best40 = 40 if n in (32, 36) else GLOBAL[n]
    for label, beta in (("best beta>=40", best40), ("global minimum", GLOBAL[n])):
        c = PrimalUSVP.cost_gsa(beta=beta, params=p, m=m + n)
        print("    %-15s beta=%d d=%s log2rop=%.2f" % (label, beta, c["d"], l2(c)))
