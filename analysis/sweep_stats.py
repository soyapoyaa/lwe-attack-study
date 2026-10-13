import os
import pandas as pd, numpy as np, json
from scipy.optimize import minimize
R=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results") + "/raw/"
p=pd.read_csv(R+"prob.csv")
gsa={32:20,36:34,40:44,44:54,48:63}

def wilson(k,n,z=1.96):
    ph=k/n; d=1+z*z/n; c=(ph+z*z/(2*n))/d; h=z*np.sqrt(ph*(1-ph)/n+z*z/(4*n*n))/d; return c-h,c+h
def cross(x,y,first=True):
    idx=[i for i in range(1,len(x)) if y[i-1]<0.5<=y[i]]
    if not idx: return np.nan
    i=idx[0] if first else idx[-1]
    return x[i-1]+(0.5-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1])
def pava(y):
    out=[]
    for v in map(float,y):
        out.append([v,1,1])
        while len(out)>1 and out[-2][0]>out[-1][0]:
            b2=out.pop(); b1=out.pop(); w=b1[1]+b2[1]; out.append([(b1[0]*b1[1]+b2[0]*b2[1])/w,w,b1[2]+b2[2]])
    r=[]
    for v,w,c in out: r+=[v]*c
    return np.array(r)
def nll(t,X,Y):
    pr=1/(1+np.exp(-(X-t[0])/np.exp(t[1]))); pr=np.clip(pr,1e-12,1-1e-12)
    return -np.sum(Y*np.log(pr)+(1-Y)*np.log(1-pr))
def logit(g,init=None):
    X=g.beta.values.astype(float); Y=g.solved.values.astype(float)
    s=g.groupby('beta').solved.mean(); x0=float(s.index[np.argmin(abs(s.values-.5))])
    o=minimize(nll,init or [x0,0.5],args=(X,Y),method='Nelder-Mead',options=dict(xatol=1e-4,fatol=1e-6,maxiter=4000)); return o.x
def ests(g):
    s=g.groupby('beta').solved.mean(); x=np.array(s.index,float); y=s.values
    t=logit(g)
    return dict(first=cross(x,y,True),last=cross(x,y,False),iso=cross(x,pava(y),True),logit=t[0],k=np.exp(t[1]))

def analyse(df,label,B=600,seed=0):
    rng=np.random.default_rng(seed); out={}
    for n,g in df.groupby('n'):
        e=ests(g); seeds=g.seed.unique(); bs={k:[] for k in e}
        for _ in range(B):
            pick=rng.choice(seeds,len(seeds))
            gg=pd.concat([g[g.seed==s_] for s_ in pick])
            eb=ests(gg)
            for k in bs: bs[k].append(eb[k])
        out[int(n)]={k:(float(e[k]),[float(v) for v in np.nanpercentile(bs[k],[2.5,97.5])]) for k in e}
    return out

if __name__=="__main__":
    print("== successes /20 by beta ==")
    for n,g in p.groupby('n'):
        s=g.groupby('beta').solved.sum().astype(int); print(n," ".join(f"{b}:{v}" for b,v in s.items()))
    print("\n== non-monotone drops in empirical p ==")
    for n,g in p.groupby('n'):
        s=g.groupby('beta').solved.sum(); d=s.diff().dropna(); print(n,{int(b):f"{int(s.shift(1)[b])}->{int(s[b])}" for b in d[d<0].index})

    print("\n== p at the GSA-predicted beta (grid has even beta; gsa values are even) ==")
    tot=0
    pa={}
    for n,g in p.groupby('n'):
        k=int(g[g.beta==gsa[n]].solved.sum()); N=int((g.beta==gsa[n]).sum()); tot+=k
        lo,hi=wilson(k,N); pa[int(n)]=(k,N,lo,hi); print(f"n={n} beta={gsa[n]}: {k}/{N} = {k/N:.2f}  Wilson95=[{lo:.2f},{hi:.2f}]")
    lo,hi=wilson(tot,80); print(f"pooled n=32..44: {tot}/80 = {tot/80:.3f}  Wilson95=[{lo:.3f},{hi:.3f}]")

    print("\n== beta_50 (ALL 1360 trials, timeouts counted as recorded=unsolved) ==")
    A=analyse(p,"all")
    for n,d in A.items():
        print(n,"gsa",gsa[n],{k:(round(v[0],1),[round(z,1) for z in v[1]]) for k,v in d.items() if k!='k'},"| logistic scale k=%.2f -> 10-90%% width=%.1f"%(d['k'][0],2*d['k'][0]*np.log(9)))
    print("\n== SENSITIVITY: drop the 10 censored (timed-out) trials ==")
    p2=p[p.timed_out==0]
    B2=analyse(p2,"drop-censored")
    for n,d in B2.items():
        print(n,{k:round(v[0],1) for k,v in d.items() if k!='k'})
    json.dump(dict(all=A,dropcens=B2,at_pred={k:list(v) for k,v in pa.items()},pooled=[tot,80,*wilson(tot,80)]),open("sweep_stats.json","w"),indent=1,default=float)
