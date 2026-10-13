import sys, os, csv

sys.path.insert(0, os.path.join(os.getcwd(), "lattice-estimator"))

from estimator import LWE, ND
from estimator.lwe_primal import PrimalUSVP
from sage.all import oo, log, RR

Q, SIGMA = 257, 4.0

# n=40 was already exhaustively scanned.
for n in (32, 36, 44, 48):

    m_avail = 8 * n
    m_ext = m_avail + n

    X = ND.DiscreteGaussian(SIGMA)

    p = LWE.Parameters(
        n=n,
        q=Q,
        Xs=X,
        Xe=X,
        m=m_avail,
        tag="scan-n%d" % n
    )

    path = "results/scan_n%d.csv" % n

    # Read any existing results, if present.
    existing = []
    last_beta = 3

    if os.path.exists(path):
        with open(path, "r") as fh:
            r = csv.DictReader(fh)
            for row in r:
                existing.append(
                    (
                        int(row["beta"]),
                        int(row["d"]),
                        float(row["log2rop"])
                    )
                )

        if existing:
            last_beta = max(x[0] for x in existing)

    # Create CSV if needed.
    if not os.path.exists(path):
        with open(path, "w") as fh:
            w = csv.writer(fh)
            w.writerow(["beta", "d", "log2rop"])

    # Current best result.
    best = None

    for row in existing:
        if best is None or row[2] < best[2]:
            best = row

    start_beta = max(4, last_beta + 1)

    print("=== n=%d ===" % n)
    print("  existing data through beta=%d" % last_beta)
    print("  extending through beta=%d" % m_ext)

    if start_beta > m_ext:
        print("  already complete")
        print(
            "  global minimum: beta=%d d=%s log2rop=%.6f"
            % best
        )
        continue

    with open(path, "a") as fh:
        w = csv.writer(fh)

        for beta in range(start_beta, m_ext + 1):

            c = PrimalUSVP.cost_gsa(
                beta=beta,
                params=p,
                m=m_ext
            )

            rop = c["rop"]

            if rop == oo:
                continue

            log2rop = float(RR(log(rop, 2)))
            d = c["d"]

            w.writerow([beta, d, log2rop])
            fh.flush()

            row = (beta, d, log2rop)

            if best is None or log2rop < best[2]:
                best = row
                print(
                    "  NEW BEST: beta=%d d=%s log2rop=%.6f"
                    % best
                )

            if beta % 10 == 0:
                print("  reached beta=%d" % beta)

    print(
        "  global minimum: beta=%d d=%s log2rop=%.6f"
        % best
    )
