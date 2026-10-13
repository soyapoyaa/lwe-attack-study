import os
import numpy as np, pandas as pd, json
from sweep_stats import p as P, logit, gsa
rng=np.random.default_rng(2026)
GRID=list(range(4,82,2))
def find(prob, hint):
    start=min(range(len(GRID)),key=lambda i:abs(GRID[i]-hint)); ok=lambda i: prob(GRID[i])
    if ok(start):
        lo_fail,hi_ok,step=-1,start,1
        while hi_ok-step>lo_fail:
            idx=hi_ok-step
            if ok(idx): hi_ok,step=idx,1
            else: lo_fail=idx; break
    else:
        lo_fail,hi_ok,step,idx=start,None,1,start
        while hi_ok is None:
            idx=idx+step
            if idx>len(GRID)-1: idx=len(GRID)-1
            if ok(idx): hi_ok=idx
            elif idx==len(GRID)-1: return None
            else: lo_fail=idx; step*=2
    while hi_ok-lo_fail>1:
        mid=(lo_fail+hi_ok)//2
        if ok(mid): hi_ok=mid
        else: lo_fail=mid
    return GRID[hi_ok]
l=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results") + "/raw/ladder.csv"); Pc=P[P.timed_out==0]
out={}
print("n | observed median (range) | sim single 5-95 | sim MEDIAN-OF-10: p5,p50,p95 | obs percentile in sim median-of-10 | obs range width vs sim 10-sample range p95")
for n in (32,36,40,44):
    t=logit(Pc[Pc.n==n]); b0,k=t[0],np.exp(t[1]); pf=lambda b: 1/(1+np.exp(-(b-b0)/k))
    B=10000; singles=[]; meds=[]; rngs=[]
    for _ in range(B):
        v=[]
        while len(v)<10:
            r=find(lambda b: rng.random()<pf(b), gsa[n])
            if r is not None: v.append(r)
        singles+=v; meds.append(np.median(v)); rngs.append(max(v)-min(v))
    singles=np.array(singles); meds=np.array(meds); rngs=np.array(rngs)
    o=l[l.n==n].beta_single; om=float(o.median()); orng=int(o.max()-o.min())
    pr=float((meds<om).mean()+0.5*(meds==om).mean())*100
    out[n]=dict(obs_med=om,obs_range=[int(o.min()),int(o.max())],obs_width=orng,
        single_p5_95=[float(np.percentile(singles,5)),float(np.percentile(singles,95))],
        med10_p5_50_95=[float(np.percentile(meds,q)) for q in (5,50,95)],obs_pctile_in_med10=pr,
        width10_p95=float(np.percentile(rngs,95)),width10_med=float(np.median(rngs)),
        frac_single_le_obs_min=float((singles<=o.min()).mean()))
    print(n,"|",om,(int(o.min()),int(o.max())),"|",out[n]['single_p5_95'],"|",out[n]['med10_p5_50_95'],"|",round(pr,1),"|",orng,"vs",out[n]['width10_p95'],"| P(single<=obs min)=%.4f"%out[n]['frac_single_le_obs_min'])
json.dump(out,open("anchor_sim.json","w"),indent=1)
