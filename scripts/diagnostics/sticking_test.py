"""Does the gating mode correlate with nanoparticle mobility?  Compare within-state widths and levels
between slow and fast gating episodes (2-min classification)."""
import numpy as np
dt=0.1; W=1200
print(f"{'trace':5s} {'mode':4s} {'nframes':>8s} {'open level':>10s} {'closed lvl':>10s} {'sd_rel open':>11s} {'sd_rel clsd':>11s} {'ratio C/O':>9s} {'shot exp':>8s} {'excess open':>11s}")
for k in "1234":
    raw=np.loadtxt(f"data/{k}Thermal.txt"); I=raw/raw.mean(); thr=(np.percentile(I,20)+np.percentile(I,95))/2; s=(I>thr).astype(int)
    tr=np.r_[np.abs(np.diff(s)),0]; r=tr.reshape(-1,W).sum(1)/(W*dt); mode=np.repeat((r>0.3).astype(int),W)
    j=np.flatnonzero(np.diff(s))+1; adj=np.zeros(s.size,bool); adj[j]=True; adj[j-1]=True
    for m,nm in [(0,"slow"),(1,"fast")]:
        mo=(s==0)&~adj&(mode==m); mc=(s==1)&~adj&(mode==m)
        lo,lc=I[mo].mean(),I[mc].mean(); so,sc=I[mo].std()/lo, I[mc].std()/lc
        shot=np.sqrt(lo/lc); exc=np.sqrt(max(so**2-(sc/shot)**2,0))
        print(f"{k:5s} {nm:4s} {mo.sum()+mc.sum():8d} {lo:10.4f} {lc:10.4f} {so:11.4f} {sc:11.4f} {sc/so:9.3f} {shot:8.3f} {exc:11.4f}")
