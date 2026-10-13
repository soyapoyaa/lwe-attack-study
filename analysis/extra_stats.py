import os
import pandas as pd, numpy as np, json, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import gsa_estimate as g
p=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results") + "/raw/prob.csv")
N=json.load(open("final_numbers.json")); S=N["sweep"]; gsa={32:20,36:34,40:44,44:54,48:63}
out={}; rng=np.random.default_rng(5)
# (a) pooled success at beta_GSA: stratified bootstrap (resample instances within each n)
rows={n:p[(p.n==n)&(p.beta==gsa[n])].sort_values("seed").solved.values for n in (32,36,40,44)}
est=[sum(rng.choice(v,len(v)).sum() for v in rows.values())/80 for _ in range(20000)]
lo,hi=np.percentile(est,[2.5,97.5]); out["pooled_strat"]=[float(np.mean(est)),float(lo),float(hi)]
print("(a) pooled 54/80: stratified-bootstrap 95%% interval = %.3f-%.3f  (Wilson was 0.566-0.768)"%(lo,hi))
# (b) beta50 relative to prediction under all four estimators
print("(b) ratio beta50/beta_GSA by estimator:")
rr={}
for k in ("first","last","iso","logit"):
    r=[S[str(n)][k]/gsa[n] for n in (32,36,40,44)]; rr[k]=r; print("   %-6s"%k,[round(x,3) for x in r],"-> %.1f%% to %.1f%% below"%(100*(1-max(r)),100*(1-min(r))))
out["ratios"]=rr
# (c) reversals within instances (same instance solved at beta, unsolved at next beta), censored trials dropped
print("(c) instance-level reversals (solved at one beta, unsolved at the next grid value):")
rev={}
for n in (32,36,40,44):
    g_=p[(p.n==n)&(p.timed_out==0)]; piv=g_.pivot(index="seed",columns="beta",values="solved")
    betas=list(piv.columns); down=up=0; inst=0
    for sd,row in piv.iterrows():
        d=0
        for a,b in zip(betas[:-1],betas[1:]):
            x,y=row[a],row[b]
            if pd.notna(x) and pd.notna(y):
                if x==1 and y==0: down+=1; d+=1
                if x==0 and y==1: up+=1
        inst+= d>0
    rev[n]=dict(instances_with_reversal=int(inst),down=int(down),up=int(up)); print("   n=%d: %d of 20 instances have at least one 1->0 step; 1->0 steps=%d, 0->1 steps=%d"%(n,inst,down,up))
out["reversals"]=rev
# n=36: solved at 30, unsolved at 32
g36=p[(p.n==36)&(p.timed_out==0)]; a=g36[g36.beta==30].set_index("seed"); b=g36[g36.beta==32].set_index("seed")
flip=[s for s in a.index if a.loc[s,"solved"]==1 and b.loc[s,"solved"]==0]
print("   n=36 instances solved at beta=30 but not at 32:",flip)
print("   tours used at beta=30 (solved):",[int(a.loc[s,'tours_used']) for s in flip],"| tours used at beta=32 (failed):",[int(b.loc[s,'tours_used']) for s in flip])
out["flip36"]=dict(seeds=[int(s) for s in flip],tours30=[int(a.loc[s,'tours_used']) for s in flip],tours32=[int(b.loc[s,'tours_used']) for s in flip])
# (d) tours: do failures stop early (auto-abort) and when do successes happen?
f=p[(p.solved==0)&(p.timed_out==0)]; s_=p[p.solved==1]
print("(d) failed (uncensored) trials: %d; of these used all 8 tours: %d; fewer than 8: %d"%(len(f),(f.tours_used==8).sum(),(f.tours_used<8).sum()))
print("    solved trials: %d; solved by tour 1: %.1f%%, by tour 2: %.1f%%, by tour 4: %.1f%%, needed more than 4: %d"%(len(s_),100*(s_.tours_used<=1).mean(),100*(s_.tours_used<=2).mean(),100*(s_.tours_used<=4).mean(),(s_.tours_used>4).sum()))
for n in (32,36,40,44):
    t=s_[s_.n==n]; print("      n=%d: solved=%d, needed >4 tours: %d, max tours used: %d"%(n,len(t),(t.tours_used>4).sum(),t.tours_used.max()))
out["tours"]=dict(fail_total=int(len(f)),fail_all8=int((f.tours_used==8).sum()),fail_lt8=int((f.tours_used<8).sum()),solved_total=int(len(s_)),solved_le4=int((s_.tours_used<=4).sum()),solved_gt4=int((s_.tours_used>4).sum()),solved_le4_frac=float((s_.tours_used<=4).mean()))
# (e) hand check of the success condition at the attacked dimension
print("(e) hand check: left side sqrt(sigma^2 (beta-1)+tau^2) vs right side delta^(2beta-d-1) * (tau * xi^n * q^(d-n-1))^(1/d)")
q,sig,tau,xi=257,4.0,4.0,1.0; hand=[]
for n in (32,36,40,44,48):
    d=N["ladder"][str(n)]["d"]
    for beta in (gsa[n]-1,gsa[n]):
        dl=g.delta_f(beta); lhs=math.sqrt(sig**2*(beta-1)+tau**2)
        rhs=math.exp((2*beta-d-1)*math.log(dl)+(math.log(tau)+n*math.log(xi)+(d-n-1)*math.log(q))/d)
        hand.append(dict(n=n,d=d,beta=beta,delta=dl,lhs=lhs,rhs=rhs,holds=bool(lhs<=rhs)))
        print("   n=%d d=%d beta=%d delta0=%.5f left=%.3f right=%.3f holds=%s"%(n,d,beta,dl,lhs,rhs,lhs<=rhs))
out["hand"]=hand
json.dump(out,open("extra_stats.json","w"),indent=1)
