"""Is there resolvable structure beyond a two-level telegraph signal?
Checks: (1) within-state peak shape vs Gaussian, (2) within-state autocorrelation,
(3) noise variance vs intensity (shot-noise scaling) across traces/states,
(4) excess of intermediate-intensity frames beyond transition straddling,
(5) dwell-time survival curves.  Output: figures/within_state.png + printed stats.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
dt = 0.1

def load(k):
    x = np.loadtxt(f"data/{k}Thermal.txt"); return x, x/x.mean()

def two_level(I):
    # crude 2-Gaussian EM on 1D
    mu = np.array([np.percentile(I,20), np.percentile(I,95)]); sg = np.array([0.02,0.02]); w = np.array([.7,.3])
    for _ in range(50):
        p = w*np.exp(-0.5*((I[:,None]-mu)/sg)**2)/sg; r = p/p.sum(1,keepdims=True)
        w = r.mean(0); mu = (r*I[:,None]).sum(0)/r.sum(0); sg = np.sqrt((r*(I[:,None]-mu)**2).sum(0)/r.sum(0))
    return mu, sg, w

def acf(res, maxlag):
    res = res - res.mean(); v = (res*res).mean()
    return np.array([ (res[:-l]*res[l:]).mean()/v if l else 1.0 for l in range(maxlag+1)])

fig, axs = plt.subplots(2, 3, figsize=(15, 8.5)); axs = axs.ravel()
noise_pts = []; summary = {}
for k in "1234":
    raw, I = load(k); mu, sg, w = two_level(I); thr = mu.mean()
    s = (I > thr).astype(int)                      # 0 = open (low), 1 = closed (high)
    jumps = np.flatnonzero(np.diff(s)) + 1         # index of first frame of a new dwell
    bounds = np.r_[0, jumps, s.size]
    dwl = np.diff(bounds); dst = s[bounds[:-1]]
    # frames not adjacent to a transition
    adj = np.zeros(s.size, bool); adj[jumps] = True; adj[jumps-1] = True
    stats = {}
    for st, name in [(0,"open"),(1,"closed")]:
        m = (s==st) & ~adj; r = I[m]-I[m].mean()
        sk = (r**3).mean()/r.std()**3; ku = (r**4).mean()/r.std()**4 - 3
        stats[name] = dict(mu=I[m].mean(), sd=I[m].std(), skew=sk, exkurt=ku, n=m.sum())
        noise_pts.append((raw[m].mean(), raw[m].var(), k, name))
    # within-state ACF on long dwells (>=100 frames), per-dwell mean removed, lags<=20
    L=100; maxlag=20; ac = {0:[],1:[]}
    for b0,b1,st in zip(bounds[:-1],bounds[1:],dst):
        if b1-b0 >= L:
            seg = I[b0+1:b1-1]; ac[st].append(acf(seg, maxlag)*len(seg))
    for st,name,c in [(0,"open","C0"),(1,"closed","C3")]:
        if ac[st]:
            a = np.sum(ac[st],0)/sum(len(I[b0+1:b1-1]) for b0,b1,s_ in zip(bounds[:-1],bounds[1:],dst) if b1-b0>=L and s_==st)
            axs[1].plot(np.arange(1,maxlag+1)*dt, a[1:], "o-" if k=="1" else ".-", ms=4, color=c, alpha=0.4+0.15*int(k), label=f"{name}, trace {k}")
            stats[name]["acf1"] = a[1]; stats[name]["acf5"] = a[5]
    # intermediate frames: between mu_open+2sd and mu_closed-2sd
    lo = stats["open"]["mu"]+2*stats["open"]["sd"]; hi = stats["closed"]["mu"]-2*stats["closed"]["sd"]
    mid = (I>lo)&(I<hi); n_mid = mid.sum(); n_tr = jumps.size
    exp_straddle = n_tr*(hi-lo)/(stats["closed"]["mu"]-stats["open"]["mu"])  # 1 uniformly-placed straddle frame per transition
    mid_at_tr = (mid & adj).sum()
    summary[k] = dict(w_closed=w[1], n_tr=n_tr, n_mid=n_mid, mid_at_tr=mid_at_tr, exp_straddle=exp_straddle, stats=stats,
                      mean_dwell_open=dwl[dst==0].mean()*dt, mean_dwell_closed=dwl[dst==1].mean()*dt,
                      frac1frame_open=(dwl[dst==0]==1).mean(), frac1frame_closed=(dwl[dst==1]==1).mean())
    # fine histogram (log) with Gaussian overlays, trace 1 only
    if k=="1":
        ax=axs[0]; h,e=np.histogram(I,bins=300,density=True); c=(e[1:]+e[:-1])/2
        ax.semilogy(c,h,"k-",lw=1,label="data")
        for st,name,col in [(0,"open","C0"),(1,"closed","C3")]:
            g=w[st]*np.exp(-0.5*((c-mu[st])/sg[st])**2)/(sg[st]*np.sqrt(2*np.pi)); ax.semilogy(c,g,"--",color=col,label=f"Gaussian {name}")
        ax.axvspan(lo,hi,color="0.9",label="intermediate zone"); ax.set_ylim(1e-3,60); ax.set_xlabel("I/<I>"); ax.set_ylabel("density"); ax.set_title("trace 1: histogram (log) vs 2-Gaussian"); ax.legend(fontsize=8)
    # survival curves
    for st,name,ax in [(0,"open",axs[3]),(1,"closed",axs[4])]:
        d = np.sort(dwl[dst==st])*dt; S = 1-np.arange(d.size)/d.size
        ax.loglog(d, S, lw=1, label=f"trace {k}"); ax.set_xlabel("dwell time (s)"); ax.set_ylabel("P(t >= tau)"); ax.set_title(f"{name} dwell survival")
    # intensity at frames around transitions: trace 1
    if k=="1":
        ax=axs[5]
        for off,lab in [(-2,"t-2"),(-1,"t-1"),(0,"t (first of new dwell)"),(1,"t+1")]:
            idx = jumps+off; idx=idx[(idx>=0)&(idx<I.size)]
            ax.hist(I[idx],bins=120,range=(0.9,1.18),histtype="step",lw=1.2,label=lab)
        ax.axvspan(lo,hi,color="0.9"); ax.set_xlabel("I/<I>"); ax.set_ylabel("count"); ax.set_title("trace 1: intensity of frames around open<->closed jumps"); ax.legend(fontsize=8)

# noise scaling
ax=axs[2]; m=np.array([p[0] for p in noise_pts]); v=np.array([p[1] for p in noise_pts])
for (mm,vv,k,name) in noise_pts: ax.plot(mm,vv,"o" if name=="open" else "s",color=f"C{int(k)-1}",label=f"trace {k} {name}")
A=np.vstack([m,np.ones_like(m)]).T; a,b=np.linalg.lstsq(A,v,rcond=None)[0]; xx=np.linspace(0,m.max()*1.05,10)
ax.plot(xx,a*xx+b,"k--",lw=1,label=f"fit var = {a:.1f}*I + {b:.0f}")
ax.set_xlabel("mean counts in state"); ax.set_ylabel("variance of counts in state"); ax.set_title("noise variance vs intensity (shot-noise test)"); ax.legend(fontsize=7)
axs[1].axhline(0,color="k",lw=0.5); axs[1].set_xlabel("lag (s)"); axs[1].set_ylabel("ACF of within-state residual"); axs[1].set_title("within-state autocorrelation (dwells >= 10 s)"); axs[1].legend(fontsize=7,ncol=2)
axs[3].legend(fontsize=8); axs[4].legend(fontsize=8)
fig.tight_layout(); fig.savefig("figures/within_state.png",dpi=140); print("wrote figures/within_state.png")
print(f"noise fit: var = {a:.2f}*I + {b:.0f}  (pure Poisson with gain g: var = g*I)")
for k,S in summary.items():
    print(f"\ntrace {k}: closed fraction {S['w_closed']:.3f}, transitions {S['n_tr']}, mean dwell open {S['mean_dwell_open']:.1f}s closed {S['mean_dwell_closed']:.1f}s, "
          f"1-frame dwells: open {S['frac1frame_open']:.2f} closed {S['frac1frame_closed']:.2f}")
    for name,st in S['stats'].items():
        print(f"   {name:6s}: mu={st['mu']:.4f} sd={st['sd']:.4f} skew={st['skew']:+.2f} exkurt={st['exkurt']:+.2f} acf(0.1s)={st.get('acf1',np.nan):+.3f} acf(0.5s)={st.get('acf5',np.nan):+.3f} n={st['n']}")
    print(f"   intermediate frames: {S['n_mid']} total, {S['mid_at_tr']} adjacent to a jump; expected from straddling ~{S['exp_straddle']:.0f}")
