import os
import pandas as pd, numpy as np
R=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results") + "/raw/"
l=pd.read_csv(R+"ladder.csv")
gsa={32:20,36:34,40:44,44:54,48:63}; dgsa={32:104,36:121,40:129,44:136,48:145}
print("LADDER RESULTS (paper Table 3)")
rows=[]
for n,g in l.groupby("n"):
    s=g.beta_single; q=g.beta_prog.dropna()
    rows.append(dict(n=n,m=int(g.m.iloc[0]),d=int(g.d.iloc[0]),gsa=gsa[n],
        s_med=s.median(),s_min=s.min(),s_max=s.max(),s_mean=round(s.mean(),1),
        ratio=round(s.median()/gsa[n],3),
        p_n=len(q),p_med=q.median(),p_min=q.min(),p_max=q.max(),
        p_ratio=round(q.median()/gsa[n],3),
        gap_vs_single_med=round(1-q.median()/s.median(),3)))
T=pd.DataFrame(rows); print(T.to_string(index=False))
print("\nPER-INSTANCE paired comparison (prog solved only)")
for n,g in l.groupby("n"):
    h=g.dropna(subset=["beta_prog"])
    lt=(h.beta_prog<h.beta_single).sum(); eq=(h.beta_prog==h.beta_single).sum(); gt=(h.beta_prog>h.beta_single).sum()
    dif=(h.beta_single-h.beta_prog)
    print(f"n={n}: paired={len(h)}  prog<single:{lt}  equal:{eq}  prog>single:{gt}   mean(single-prog)={dif.mean():.1f}  median={dif.median():.1f}")
    if gt: print("      prog>single instances:",h[h.beta_prog>h.beta_single][['seed','beta_single','beta_prog']].values.tolist())
allh=l.dropna(subset=["beta_prog"])
print("ALL:",len(allh),"paired; prog<single",(allh.beta_prog<allh.beta_single).sum(),"eq",(allh.beta_prog==allh.beta_single).sum(),"prog>single",(allh.beta_prog>allh.beta_single).sum())
print("\nsingle-BKZ instance beta relative to GSA (per instance)")
for n,g in l.groupby("n"):
    r=g.beta_single/gsa[n]; print(n,"ratio min/med/max",round(r.min(),2),round(r.median(),3),round(r.max(),2),
      "| frac instances with beta_single<=gsa:",round((g.beta_single<=gsa[n]).mean(),2), " <gsa:",round((g.beta_single<gsa[n]).mean(),2))
print("\nsingle_probes (search cost) per n:",l.groupby('n').single_probes.mean().round(1).to_dict())
print("first probe start = grid point nearest gsa; beta parity single:",(l.beta_single%2).value_counts().to_dict())
T.to_csv("table1_fresh.csv",index=False)
