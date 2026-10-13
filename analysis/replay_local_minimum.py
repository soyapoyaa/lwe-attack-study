"""Replay of lattice-estimator's local_minimum search (estimator/util.py, commit 53da598) over the exhaustive cost_gsa scans in
results/scan_n*.csv. Reproduces the beta values primal_usvp returned (65, 73, 81, 54, 63 for n = 32..48) and the first probe
(beta = 355 at n = 40 on a grid of step 5). Run from this directory: python3 replay_local_minimum.py"""
import os
"""Replay estimator.util.local_minimum (faithful port) over real cost_gsa scan data."""
import math, pandas as pd
INF = float("inf")
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results") + "/"

def ceil(x): return math.ceil(x)
def floor(x): return math.floor(x)

class LocalMin:
    def __init__(s, start, stop, precision=1):
        s.prec = precision; s.orig = (start, stop)
        start = ceil(start / precision); stop = floor(stop / precision)
        s._start = start; s._stop = stop - 1
        s.ib = (start, stop - 1)
        s.dir = -1; s.last = None; s.next = s._stop
        s.best_low = None; s.best_high = None; s.all = set()
    def __iter__(s): return s
    def __next__(s):
        if s.next is not None and s.next not in s.all and s.ib[0] <= s.next <= s.ib[1]:
            s.last = s.next; s.next = None
            return s.last * s.prec
        raise StopIteration
    @property
    def x(s): return s.best_low * s.prec
    @property
    def neighborhood(s):
        a, b = s.orig
        return range(max(a, s.x - s.prec), min(b, s.x + s.prec))
    def update(s, res):
        s.all.add(s.last)
        if s.best_low is None:
            s.best_low, s.best_high = s.last, res
        if res <= s.best_high:            # smallerf: x <= best
            s.best_low, s.best_high = s.last, res
            if abs(s.dir) != 1:
                s.dir = -1; s.next = s.last - 1
            elif s.dir == -1:
                s.dir = -2; s._stop = s.last; s.next = ceil((s._start + s._stop) / 2)
            elif s.dir == 1:
                s.dir = 2; s._start = s.last; s.next = floor((s._start + s._stop) / 2)
        else:
            if s.dir == -1: s.dir = 1; s.next = s.last + 2
            elif s.dir == 1: s.next = None
            elif s.dir == -2: s._start = s.last; s.next = ceil((s._start + s._stop) / 2)
            elif s.dir == 2: s._stop = s.last; s.next = floor((s._start + s._stop) / 2)
        if s.next == s.last: s.next = None

def replay(n, lo=40, precision=5):
    s = pd.read_csv(R + f"scan_n{n}.csv")
    tab = {int(r.beta): (int(r.d), float(r.log2rop)) for r in s.itertuples()}
    f = lambda b: tab[b][1] if b in tab else INF
    m = 8 * n + n                                   # Xs<=Xe -> m = params.m + n
    max_beta = max(min(1754, m), lo + 5)
    it = LocalMin(lo, max_beta, precision)
    probes = []
    for b in it:
        probes.append(b); it.update(f(b))
    best_after_main = (it.x, f(it.x))
    best_y = it.best_high
    for b in it.neighborhood:
        probes.append(b)
        v = f(b)
        if v <= it.best_high: it.best_high = v; ybeta = b
    # y is the Cost object whose beta field is the beta that produced it
    # find beta whose cost equals best_high
    chosen = [b for b in tab if abs(tab[b][1] - it.best_high) < 1e-12]
    gb = min(tab, key=lambda b: tab[b][1])
    return dict(n=n, m=m, main_min=best_after_main, chosen_beta=chosen, chosen_d=[tab[b][0] for b in chosen],
                chosen_log2rop=round(it.best_high, 2), global_beta=gb, global_d=tab[gb][0],
                global_log2rop=round(tab[gb][1], 2), probes_main=len(probes))

if __name__ == "__main__":
    for n in (32, 36, 40, 44, 48):
        print(replay(n))
