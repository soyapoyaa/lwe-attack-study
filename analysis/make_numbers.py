import os
import json, numpy as np, pandas as pd
from sweep_stats import p, logit, ests, analyse, gsa, wilson, cross
R=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results") + "/raw/"
l=pd.read_csv(R+"ladder.csv")
N={}
pc=p[p.timed_out==0]
# --- sweep, primary = censored excluded
A=analyse(pc,"primary",B=1000,seed=11)
Aall=analyse(p,"asrecorded",B=400,seed=11)
sw={}
for n,g in pc.groupby('n'):
    n=int(n); e=A[n]; gg=p[p.n==n]
    k_at=int(gg[gg.beta==gsa[n]].solved.sum()); lo,hi=wilson(k_at,20)
    s=g.groupby('beta').solved.mean()
    b90=int(s[s>=0.9].index.min()); 
    counts={int(b):[int(v.solved.sum()),int(len(v))] for b,v in g.groupby('beta')}
    rec={int(b):[int(v.solved.sum()),int(len(v))] for b,v in gg.groupby('beta')}
    ncens=int((gg.timed_out==1).sum())
    sw[n]=dict(gsa=gsa[n],k_at=k_at,p_at=k_at/20,p_lo=lo,p_hi=hi,
        logit=e['logit'][0],logit_ci=e['logit'][1],k=e['k'][0],w1090=2*e['k'][0]*np.log(9),
        first=e['first'][0],first_ci=e['first'][1],last=e['last'][0],last_ci=e['last'][1],iso=e['iso'][0],iso_ci=e['iso'][1],
        ratio_logit=e['logit'][0]/gsa[n],ratio_first=e['first'][0]/gsa[n],ratio_last=e['last'][0]/gsa[n],ratio_iso=e['iso'][0]/gsa[n],
        ratio_logit_ci=[e['logit'][1][0]/gsa[n],e['logit'][1][1]/gsa[n]],
        b90=b90,counts=counts,recorded=rec,ncens=ncens,
        logit_asrecorded=Aall[n]['logit'][0], k_asrecorded=Aall[n]['k'][0], w1090_asrecorded=2*Aall[n]['k'][0]*np.log(9))
N['sweep']=sw
tot=sum(sw[n]['k_at'] for n in sw); lo,hi=wilson(tot,80); N['pooled']=dict(k=tot,n=80,p=tot/80,lo=lo,hi=hi)
N['n_trials']=int(len(p)); N['n_solved']=int(p.solved.sum()); N['n_cens']=int(p.timed_out.sum())
# --- ladder
lad={}
for n,g in l.groupby('n'):
    n=int(n); s=g.beta_single; q=g.beta_prog.dropna(); h=g.dropna(subset=['beta_prog'])
    lad[n]=dict(m=int(g.m.iloc[0]),d=int(g.d.iloc[0]),gsa=gsa[n],gap=float(g.gap.iloc[0]),
       s_med=float(s.median()),s_min=int(s.min()),s_max=int(s.max()),s_all=sorted(map(int,s)),ratio=float(s.median()/gsa[n]),
       p_n=int(len(q)),p_med=float(q.median()),p_min=int(q.min()),p_max=int(q.max()),p_all=sorted(map(int,q)),p_ratio=float(q.median()/gsa[n]),
       med_gap=float(1-q.median()/s.median()),
       pair_lt=int((h.beta_prog<h.beta_single).sum()),pair_eq=int((h.beta_prog==h.beta_single).sum()),pair_gt=int((h.beta_prog>h.beta_single).sum()),
       t_med=float(g.time_single_at_beta.median()),t_min=float(g.time_single_at_beta.min()),t_max=float(g.time_single_at_beta.max()),
       n_timed_out_rows=int(g.timed_out.sum()))
N['ladder']=lad
N['pairs_total']=dict(n=sum(v['p_n'] for v in lad.values()),lt=sum(v['pair_lt'] for v in lad.values()),eq=sum(v['pair_eq'] for v in lad.values()),gt=sum(v['pair_gt'] for v in lad.values()))
# --- timing fit
x=l.n.values.astype(float); y=np.log2(l.time_single_at_beta.values)
A_=np.vstack([x,np.ones_like(x)]).T; b,a=np.linalg.lstsq(A_,y,rcond=None)[0]
rng=np.random.default_rng(3); bs=[np.linalg.lstsq(A_[i],y[i],rcond=None)[0][0] for i in (rng.integers(0,len(x),len(x)) for _ in range(4000))]
yh=a+b*x; N['timefit']=dict(slope=float(b),intercept=float(a),ci=[float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))],r2=float(1-((y-yh)**2).sum()/((y-y.mean())**2).sum()))
# --- anchor sim
N['anchor2']=json.load(open("anchor_sim.json")); N['est_table']=json.load(open("est_table.json"))
N['floor_best']={k:v['b40'][:2] for k,v in N['est_table'].items()}
# --- misc totals
N['cpu_hours']=dict(prob=float(p.time.sum()/3600),ladder=float((l.time_single_total.sum()+l.time_prog_total.sum())/3600))
N['tours']=p.tours_used.value_counts().sort_index().to_dict()
json.dump(N,open("final_numbers.json","w"),indent=1,default=lambda o: o.item() if hasattr(o,'item') else str(o))
for n,v in sw.items(): print(n,{k:(round(x,3) if isinstance(x,float) else x) for k,x in v.items() if k in('logit','logit_ci','w1090','first','last','iso','ratio_logit','b90','k_at','p_lo','p_hi','ncens','logit_asrecorded','w1090_asrecorded')})
print(N['pooled'],N['pairs_total'],N['timefit'],N['cpu_hours'],N['tours'])
