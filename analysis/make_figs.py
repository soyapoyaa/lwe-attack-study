import os
os.makedirs('figs', exist_ok=True)
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
N=json.load(open("final_numbers.json"))
R=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results") + "/"
p=pd.read_csv(R+"raw/prob.csv"); l=pd.read_csv(R+"raw/ladder.csv")
gsa={32:20,36:34,40:44,44:54,48:63}
# Okabe-Ito
BLK="#000000"; BLU="#0072B2"; ORG="#E69F00"; GRN="#009E73"; VER="#D55E00"; GRY="#7f7f7f"
plt.rcParams.update({"font.family":"DejaVu Serif","font.size":9,"axes.spines.top":False,"axes.spines.right":False,
    "axes.grid":True,"grid.alpha":.25,"grid.linewidth":.5,"figure.dpi":200,"savefig.dpi":200,"legend.frameon":False})
ns=[32,36,40,44,48]

# ---------- Fig 1: beta vs n
fig,ax=plt.subplots(figsize=(6.4,4.1))
ax.axhline(40,color=GRY,ls=":",lw=1); ax.text(31.0,40.7,"estimator GSA search floor (β = 40)",ha="left",va="bottom",fontsize=7.5,color=GRY)
ax.plot(ns,[gsa[n] for n in ns],"-s",color=BLK,lw=1.8,ms=6,label="GSA prediction (global minimum of cost_gsa)",zorder=5)
x=np.array(ns)-0.35
med=[N['ladder'][str(n)]['s_med'] for n in ns]; lo=[med[i]-N['ladder'][str(n)]['s_min'] for i,n in enumerate(ns)]; hi=[N['ladder'][str(n)]['s_max']-med[i] for i,n in enumerate(ns)]
ax.errorbar(x,med,yerr=[lo,hi],fmt="o-",color=BLU,lw=1.3,ms=5,capsize=3,label="ladder: single-BKZ smallest solving β (median, range; 10/n)")
xp=np.array(ns)+0.35
pm=[N['ladder'][str(n)]['p_med'] for n in ns]; plo=[pm[i]-N['ladder'][str(n)]['p_min'] for i,n in enumerate(ns)]; phi=[N['ladder'][str(n)]['p_max']-pm[i] for i,n in enumerate(ns)]
ax.errorbar(xp[:4],pm[:4],yerr=[plo[:4],phi[:4]],fmt="^--",color=ORG,lw=1.3,ms=5,capsize=3,label="ladder: progressive BKZ (median, range; 10/n)")
ax.errorbar(xp[4:],pm[4:],yerr=[plo[4:],phi[4:]],fmt="^",mfc="white",color=ORG,ms=6,capsize=3)
ax.plot([xp[3],xp[4]],[pm[3],pm[4]],"--",color=ORG,lw=1.3)
ax.annotate("6 of 10 solved;\n4 hit the 900 s cap",(xp[4],pm[4]),xytext=(44.9,33.0),fontsize=7.5,color=ORG,arrowprops=dict(arrowstyle="-",color=ORG,lw=.8))
n4=[32,36,40,44]
b=[N['sweep'][str(n)]['logit'] for n in n4]; ci=[N['sweep'][str(n)]['logit_ci'] for n in n4]
ax.errorbar(n4,b,yerr=[[b[i]-ci[i][0] for i in range(4)],[ci[i][1]-b[i] for i in range(4)]],fmt="D-",color=GRN,lw=1.6,ms=5.5,capsize=3,
    label="sweep: β₅₀, 8 tours, logistic fit (95% CI)")
T4=json.load(open("tour4_stats.json")); b4=[T4[str(n)]["logit"] for n in n4]; c4=[T4[str(n)]["logit_ci"] for n in n4]
ax.errorbar(np.array(n4)-0.15,b4,yerr=[[b4[i]-c4[i][0] for i in range(4)],[c4[i][1]-b4[i] for i in range(4)]],fmt="s:",color=VER,lw=1.4,ms=5,capsize=3,label="sweep: β₅₀, 4 tours (derived)")
ax.set_xticks(ns); ax.set_xlabel("LWE secret dimension n"); ax.set_ylabel("BKZ block size β"); ax.set_ylim(8,72); ax.set_xlim(30.8,49.3)
ax.set_title("Measured vs. predicted block size (q = 257, σ = 4)",fontsize=10)
ax.legend(loc="upper left",fontsize=7.4)
fig.tight_layout(); fig.savefig("figs/fig1_beta_vs_n.png"); plt.close(fig)

# ---------- Fig 2: success probability
def wil(k,n,z=1.96):
    ph=k/n; d=1+z*z/n; c=(ph+z*z/(2*n))/d; h=z*np.sqrt(ph*(1-ph)/n+z*z/(4*n*n))/d; return c-h,c+h
fig,axs=plt.subplots(2,2,figsize=(6.6,5.2),sharey=True)
for ax,n in zip(axs.flat,[32,36,40,44]):
    S=N['sweep'][str(n)]; g=p[p.n==n]
    bs=sorted(g.beta.unique()); 
    rec=[g[g.beta==b].solved.mean() for b in bs]
    ex=[]; 
    for b in bs:
        h=g[(g.beta==b)&(g.timed_out==0)]; ex.append(h.solved.mean())
    xs=np.linspace(min(bs),max(bs),300); pf=1/(1+np.exp(-(xs-S['logit'])/S['k']))
    ax.plot(xs,pf,color=GRN,lw=1.6,zorder=2)
    ax.axvline(S['gsa'],color=BLK,ls="--",lw=1.2,zorder=1); ax.axvline(S['logit'],color=GRN,ls=":",lw=1.2,zorder=1)
    ax.axhline(.5,color=GRY,lw=.6)
    ax.plot(bs,ex,"o",color=BLU,ms=4.2,zorder=4)
    cens=[(b,r) for b,r,e in zip(bs,rec,ex) if abs(r-e)>1e-9]
    if cens:
        ax.plot([c[0] for c in cens],[c[1] for c in cens],"o",mfc="white",mec=VER,ms=5.2,zorder=5)
    ax.set_title(f"n = {n}  (m = {N['ladder'][str(n)]['m']}, d = {N['ladder'][str(n)]['d']})",fontsize=9)
    ax.text(0.97,0.06,f"solved at β_GSA = {S['gsa']}: {S['k_at']}/20\nfitted β₅₀ = {S['logit']:.1f}",transform=ax.transAxes,fontsize=7.5,va="bottom",ha="right")
    ax.set_ylim(-0.04,1.06); ax.set_xlabel("BKZ block size β")
for ax in axs[:,0]: ax.set_ylabel("fraction solved (of 20)")
h=[Line2D([],[],marker="o",color=BLU,ls="",ms=4.5,label="p̂(β), timed-out trials excluded"),
   Line2D([],[],marker="o",mfc="white",mec=VER,color=VER,ls="",ms=5.5,label="p̂(β) if timeouts counted as failures (n=44 only)"),
   Line2D([],[],color=GRN,lw=1.6,label="logistic fit"),Line2D([],[],color=BLK,ls="--",label="GSA prediction"),Line2D([],[],color=GRN,ls=":",label="fitted β₅₀")]
fig.legend(handles=h,loc="lower center",ncol=3,fontsize=7.3,bbox_to_anchor=(0.5,-0.005))
fig.tight_layout(rect=[0,0.085,1,1]); fig.savefig("figs/fig2_success_prob.png"); plt.close(fig)

# ---------- Fig 3: time vs n
fig,ax=plt.subplots(figsize=(5.6,3.7)); rng=np.random.default_rng(0)
xj=l.n.values+rng.uniform(-.45,.45,len(l)); ax.scatter(xj,l.time_single_at_beta,s=16,color=BLU,alpha=.75,label="one instance (solving call)")
mt=[N['ladder'][str(n)]['t_med'] for n in ns]; ax.plot(ns,mt,"D",color=VER,ms=6,label="median per n",zorder=5)
T=N['timefit']; xs=np.linspace(31,49,50); ax.plot(xs,2**(T['intercept']+T['slope']*xs),color=BLK,lw=1.2,label=f"fit: log₂ t = {T['intercept']:.2f} + {T['slope']:.3f}·n  (R² = {T['r2']:.2f})")
ax.set_yscale("log",base=2); ax.set_xticks(ns); ax.set_xlabel("LWE secret dimension n"); ax.set_ylabel("wall time of solving call (s)")
ax.set_yticks([1,2,4,8,16,32,64,128,256,512,1024]); ax.set_yticklabels(["1","2","4","8","16","32","64","128","256","512","1024"])
ax.set_title("Cost of the solving single-BKZ call (fpylll enumeration)",fontsize=10); ax.legend(fontsize=7.4,loc="upper left")
fig.tight_layout(); fig.savefig("figs/fig3_time_vs_n.png"); plt.close(fig)

# ---------- Fig 4: cost curve n=40
s=pd.read_csv(R+"scan_n40.csv"); s=s[(s.beta>=40)&(s.beta<=130)]
fig,ax=plt.subplots(figsize=(6.0,3.7))
ax.plot(s.beta,s.log2rop,color=BLK,lw=1.4)
gl=s.loc[s.beta==44].iloc[0]; pk=s.loc[s.beta==77].iloc[0]; l81=s.loc[s.beta==81].iloc[0]
ax.plot(44,gl.log2rop,"o",color=GRN,ms=8,zorder=5); ax.annotate("β = 44, d = 129\nglobal minimum: 2^40.7",(44,gl.log2rop),xytext=(50,38.9),fontsize=8,color=GRN,arrowprops=dict(arrowstyle="-",color=GRN,lw=.8))
ax.plot(77,pk.log2rop,"s",color=GRY,ms=6,zorder=5); ax.annotate("cost peaks at β = 77\n(2^46.5)",(77,pk.log2rop),xytext=(58,49.6),fontsize=8,color=GRY,arrowprops=dict(arrowstyle="-",color=GRY,lw=.8))
ax.plot(81,l81.log2rop,"X",color=VER,ms=9,zorder=6); ax.annotate("β = 81, d = 82: value returned by\nprimal_usvp (2^44.8) — a second basin",(81,l81.log2rop),xytext=(84,41.4),fontsize=8,color=VER,arrowprops=dict(arrowstyle="-",color=VER,lw=.8))
ax.set_xlabel("BKZ block size β"); ax.set_ylabel("log₂ rop  (cost_gsa, best d per β)"); ax.set_xlim(40,130); ax.set_ylim(38,55)
ax.set_title("n = 40: cost_gsa is not unimodal in β",fontsize=10)
fig.tight_layout(); fig.savefig("figs/fig4_cost_curve_n40.png"); plt.close(fig)
print("ok")
