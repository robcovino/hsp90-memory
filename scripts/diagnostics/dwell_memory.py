"""Model-free diagnostics on the dwell sequence: (1) correlations between successive dwells vs lag,
(2) dependency of closed dwell on preceding open dwell (Magleby-style), (3) time-reversal symmetry
of the joint (open, closed) dwell distribution (detailed-balance test), (4) stationarity of kinetics."""
import numpy as np
dt=0.1; rng=np.random.default_rng(0)
def dwells(k):
    raw=np.loadtxt(f"data/{k}Thermal.txt"); I=raw/raw.mean()
    thr=(np.percentile(I,20)+np.percentile(I,95))/2; s=(I>thr).astype(int)
    j=np.flatnonzero(np.diff(s))+1; b=np.r_[0,j,s.size]; d=np.diff(b)*dt; st=s[b[:-1]]
    return d[1:-1], st[1:-1]   # drop censored first/last
def spearman(a,b):
    ra=np.argsort(np.argsort(a)); rb=np.argsort(np.argsort(b)); return np.corrcoef(ra,rb)[0,1]
allO=[];allC=[];allOC=[];allCO=[]
print("=== (1) Spearman correlation between dwells, by lag (in dwells of the same state); shuffle 95% band in brackets")
for k in "1234":
    d,st=dwells(k); O=d[st==0]; C=d[st==1]
    # pair up: sequence alternates; make arrays of (O_i, C_i following) and (C_i, O_{i+1} following)
    if st[0]==1: d,st=d[1:],st[1:]
    n=min((st==0).sum(),(st==1).sum()); O=d[st==0][:n]; C=d[st==1][:n]
    # O_i is followed by C_i; C_i is followed by O_{i+1}
    rOC=spearman(O,C); rCO=spearman(C[:-1],O[1:])
    lags=[1,2,3,5,10,20,50,100]
    rOO=[spearman(O[:-l],O[l:]) for l in lags]; rCC=[spearman(C[:-l],C[l:]) for l in lags]
    band=1.96/np.sqrt(n)
    print(f"trace {k} (n={n} cycles): O->next C {rOC:+.3f}, C->next O {rCO:+.3f}  [±{band:.3f}]")
    print("   O-O by lag:", " ".join(f"{l}:{r:+.3f}" for l,r in zip(lags,rOO)))
    print("   C-C by lag:", " ".join(f"{l}:{r:+.3f}" for l,r in zip(lags,rCC)))
    allO.append(O);allC.append(C)
O=np.concatenate(allO); C=np.concatenate(allC); n=O.size
print(f"\npooled n={n}: O->next C {spearman(O,C):+.3f}, C->next O {spearman(C[:-1],O[1:]):+.3f}  [±{1.96/np.sqrt(n):.3f}]")
lags=[1,2,3,5,7,10,15,20,30,50,70,100,150,200]
print("pooled O-O by lag:", " ".join(f"{l}:{spearman(O[:-l],O[l:]):+.3f}" for l in lags))
print("pooled C-C by lag:", " ".join(f"{l}:{spearman(C[:-l],C[l:]):+.3f}" for l in lags))

print("\n=== (2) Dependency: mean closed dwell (s) conditional on preceding open dwell class; and mean open dwell conditional on preceding closed dwell")
eO=np.array([0,0.3,1,3,10,30,1e4]); eC=np.array([0,0.3,1,3,10,1e4])
for i in range(len(eO)-1):
    m=(O>=eO[i])&(O<eO[i+1]); print(f"  prev open in [{eO[i]:g},{eO[i+1]:g}) s: N={m.sum():5d}  <next closed>={C[m].mean():.2f} s  median={np.median(C[m]):.2f}")
for i in range(len(eC)-1):
    m=(C[:-1]>=eC[i])&(C[:-1]<eC[i+1]); print(f"  prev closed in [{eC[i]:g},{eC[i+1]:g}) s: N={m.sum():5d}  <next open>={O[1:][m].mean():.2f} s  median={np.median(O[1:][m]):.2f}")

print("\n=== (3) Time-reversal test: joint density of (open, FOLLOWING closed) vs (open, PRECEDING closed), log-binned; chi2 vs shuffle null")
lb=np.array([0,0.2,0.5,1,2,5,10,30,1e4])
Hf,_,_=np.histogram2d(O[1:],C[1:],bins=[lb,lb]); Hp,_,_=np.histogram2d(O[1:],C[:-1],bins=[lb,lb])
chi2=np.nansum((Hf-Hp)**2/(Hf+Hp+1e-9)); 
null=[]
for _ in range(200):
    # null: randomly reverse direction of time within each cycle pair assignment
    flip=rng.random(O.size-1)<0.5; Cf=np.where(flip,C[:-1],C[1:]); Cp=np.where(flip,C[1:],C[:-1])
    a,_,_=np.histogram2d(O[1:],Cf,bins=[lb,lb]); b_,_,_=np.histogram2d(O[1:],Cp,bins=[lb,lb]); null.append(np.nansum((a-b_)**2/(a+b_+1e-9)))
print(f"  chi2 observed = {chi2:.1f}; shuffle null mean {np.mean(null):.1f}, 95th pct {np.percentile(null,95):.1f}")

print("\n=== (4) Stationarity of kinetics: per-hour mean dwells and closed fraction")
for k in "1234":
    raw=np.loadtxt(f"data/{k}Thermal.txt"); I=raw/raw.mean(); thr=(np.percentile(I,20)+np.percentile(I,95))/2; s=(I>thr).astype(int)
    out=[]
    for h in range(6):
        ss=s[h*36000:(h+1)*36000]; j=np.flatnonzero(np.diff(ss))+1; b=np.r_[0,j,ss.size]; d=np.diff(b)*dt; st=ss[b[:-1]]
        out.append(f"h{h}: fC={ss.mean():.2f} nTr={j.size} <O>={d[st==0].mean():.1f} <C>={d[st==1].mean():.1f}")
    print(f"trace {k}: "+" | ".join(out))
