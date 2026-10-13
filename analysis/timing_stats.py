import os
import pandas as pd, numpy as np
R=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results") + "/raw/"
p=pd.read_csv(R+"prob.csv"); l=pd.read_csv(R+"ladder.csv")
gsa={32:20,36:34,40:44,44:54,48:63}
rng=np.random.default_rng(1)

print("A) LADDER: wall time of the solving single-BKZ call at the instance's smallest solving beta")
rows=[]
for n,g in l.groupby('n'):
    t=g.time_single_at_beta; rows.append((n,len(t),t.median(),t.min(),t.max()))
    print(f"n={n}: N={len(t)} median={t.median():.1f}s  min={t.min():.1f}  max={t.max():.1f}   log2(median)={np.log2(t.median()):.2f}")
def fit(x,y):
    A=np.vstack([x,np.ones_like(x)]).T; return np.linalg.lstsq(A,y,rcond=None)[0]
# fit on per-instance log2 time vs n
x=l.n.values.astype(float); y=np.log2(l.time_single_at_beta.values)
b,a=fit(x,y)
bs=[]
for _ in range(2000):
    i=rng.integers(0,len(x),len(x)); bs.append(fit(x[i],y[i])[0])
print(f"   per-instance OLS log2(time) = {a:.2f} + {b:.3f} n   bootstrap95 slope=[{np.percentile(bs,2.5):.3f},{np.percentile(bs,97.5):.3f}]  -> time x2 every {1/b:.1f} units of n; x{2**(4*b):.1f} per +4 in n")
# fit on medians
xm=np.array([r[0] for r in rows],float); ym=np.log2([r[2] for r in rows]); bm,am=fit(xm,ym)
print(f"   OLS on 5 medians: slope={bm:.3f} (x2 per {1/bm:.1f} n)")
# R^2
yhat=a+b*x; r2=1-((y-yhat)**2).sum()/((y-y.mean())**2).sum(); print("   R^2 per-instance:",round(r2,3))
# mixed contribution: solving beta also grows with n -> partial fit log2 t = a + b*beta_single
x2=l.beta_single.values.astype(float); b2,a2=fit(x2,y); print(f"   log2(time) vs solving beta (pooled, all n): slope={b2:.3f} per unit beta (x2 per {1/b2:.1f} beta)")
# 2-var fit
A=np.vstack([l.n.values,l.beta_single.values,np.ones(len(l))]).T; co=np.linalg.lstsq(A,y,rcond=None)[0]; print("   2-var OLS log2t = %.2f + %.3f*n + %.3f*beta"%(co[2],co[0],co[1]))

print("\nB) SWEEP: time of SOLVED single-BKZ calls at the GSA-predicted beta")
for n,g in p.groupby('n'):
    s=g[(g.beta==gsa[n])&(g.solved==1)]
    print(f"n={n} beta={gsa[n]}: solved={len(s)} median time={s.time.median():.1f}s (min {s.time.min():.1f}, max {s.time.max():.1f})")
print("\nC) SWEEP: within fixed n, log2(time of solved calls) vs beta")
for n,g in p.groupby('n'):
    s=g[g.solved==1]; s=s[s.time>0]
    b,a=fit(s.beta.values.astype(float),np.log2(s.time.values)); print(f"n={n}: slope={b:.3f} per unit beta (x2 per {1/b:.1f} beta), N={len(s)}")
print("\nD) enumerated total CPU: prob sweep total seconds =",round(p.time.sum()/3600,1),"h ; ladder single_total+prog_total =",round((l.time_single_total.sum()+l.time_prog_total.sum())/3600,1),"h")
print("   log file tail:")
