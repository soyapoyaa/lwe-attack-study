import os
import pandas as pd, numpy as np, json
from sweep_stats import p, analyse, gsa
pc=p[p.timed_out==0].copy()
# Assumption: a call with a 4-tour budget behaves like the first 4 tours of an 8-tour call, so a trial solved after more than 4 tours counts as unsolved.
pc.loc[pc.tours_used>4,"solved"]=0
A4=analyse(pc,"4tours",B=500,seed=3); out={}
for n in (32,36,40,44):
    g=pc[(pc.n==n)&(pc.beta==gsa[n])]; k=int(g.solved.sum())
    out[str(n)]=dict(logit=A4[n]["logit"][0],logit_ci=A4[n]["logit"][1],k_at=k,ratio=A4[n]["logit"][0]/gsa[n])
    print(n,"4-tour beta50 (logistic) = %.1f [%.1f, %.1f]  ratio to GSA %.2f | solved at beta_GSA with 4 tours: %d/20"%(A4[n]["logit"][0],*A4[n]["logit"][1],A4[n]["logit"][0]/gsa[n],k))
tot=sum(out[n]["k_at"] for n in out); print("pooled at beta_GSA, 4 tours: %d/80"%tot); out["pooled_k"]=tot
json.dump(out,open("tour4_stats.json","w"),indent=1)
