"""Posterior-predictive checks for the fitted aggregated Markov models.
For each model: simulate 4 x 6 h at 10 Hz from the fitted generator and compare with the data on
  (a) dwell survival curves, (b) Spearman dwell-lag correlations, (c) 2-min local-rate histogram (modal gating),
  (d) Chapman-Kolmogorov divergence profile D_tau(t) (Vollmar et al. 2024) out to 600 s.
Usage: python3 analysis/predictive_checks.py modal1 modal2 modal3 ...   -> figures/predictive_checks.png, results/predictive_checks.json"""
import sys, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from hsp90memory import amm
DT = 0.1; NF = 216000

def dwells(s):
    j = np.flatnonzero(np.diff(s)) + 1; b = np.r_[0, j, s.size]; return np.diff(b) * DT, s[b[:-1]]

def spearman(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b)); return np.corrcoef(ra, rb)[0, 1]

LAGS = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100]
def lagcorr(series_list):
    O = []; C = []
    for s in series_list:
        d, st = dwells(s); d, st = d[1:-1], st[1:-1]
        if st[0] == 1: d, st = d[1:], st[1:]
        n = min((st == 0).sum(), (st == 1).sum()); O.append(d[st == 0][:n]); C.append(d[st == 1][:n])
    O = np.concatenate(O); C = np.concatenate(C)
    return np.array([spearman(O[:-l], O[l:]) for l in LAGS]), np.array([spearman(C[:-l], C[l:]) for l in LAGS]), spearman(O, C)

def survival(series_list, cls):
    d = np.concatenate([dwells(s)[0][dwells(s)[1] == cls] for s in series_list]); d = np.sort(d)
    return d, 1 - np.arange(d.size) / d.size

def rate_hist(series_list, W=1200):
    r = np.concatenate([np.r_[np.abs(np.diff(s)), 0].reshape(-1, W).sum(1) / (W * DT) for s in series_list])
    h, e = np.histogram(r, bins=np.arange(0, 1.1, 0.05)); return h / h.sum(), e

TLAGS = np.unique(np.round(np.logspace(0, np.log10(6000), 40)).astype(int))   # 0.1 s .. 600 s in frames
def propagator(s, L):
    a = s[:-L]; b = s[L:]; G = np.zeros((2, 2))
    for c0 in (0, 1):
        m = a == c0; G[c0, 0] = (b[m] == 0).mean(); G[c0, 1] = (b[m] == 1).mean()
    return G

def ck_divergence(s, tau_frames=10):
    """D_{tau,q0}(t) for both q0; returns (t in s, D[q0=open], D[q0=closed])."""
    Gtau = propagator(s, tau_frames); out = []
    for L in TLAGS:
        if L <= tau_frames: out.append((np.nan, np.nan)); continue
        G = propagator(s, L); Gck = Gtau @ propagator(s, L - tau_frames)
        D = [np.sum(G[q0] * np.log(np.clip(G[q0], 1e-12, None) / np.clip(Gck[q0], 1e-12, None))) for q0 in (0, 1)]
        out.append(D)
    out = np.array(out); return TLAGS * DT, out[:, 0], out[:, 1]

def summaries(series_list):
    oo, cc, oc = lagcorr(series_list)
    sO = survival(series_list, 0); sC = survival(series_list, 1)
    h, e = rate_hist(series_list)
    ck = np.array([np.stack(ck_divergence(s)[1:], 0) for s in series_list])    # (ntraces, 2, nlags)
    return dict(oo=oo, cc=cc, oc=oc, sO=sO, sC=sC, h=h, e=e, ck_mean=np.nanmean(ck, 0), ck_sd=np.nanstd(ck, 0))

keys = sys.argv[1:]
data = [amm.load_dwells(f"data/{k}Thermal.txt")[0] for k in "1234"]
S = {"data": summaries(data)}
models = {"modal1": amm.modal(1), "modal2": amm.modal(2), "modal2_noDB": amm.modal(2, db=False), "modal3": amm.modal(3),
          "modal4": amm.modal(4), "general2o2c": amm.general(2, 2), "general3o2c": amm.general(3, 2), "general3o3c": amm.general(3, 3)}
for key in keys:
    r = json.load(open(f"results/{key}.json")); m = models[key]; th = np.array(r["theta"])
    sims = [amm.simulate(m, th, NF, seed=100 + i) for i in range(4)]
    S[key] = summaries(sims); print("done", key, flush=True)

fig, axs = plt.subplots(2, 3, figsize=(16, 9)); axs = axs.ravel()
cols = {"data": "k"}; palette = ["C3", "C0", "C2", "C1", "C4", "C5", "C6"]
for i, k in enumerate(keys): cols[k] = palette[i % len(palette)]
for k, v in S.items():
    lw = 2.2 if k == "data" else 1.3; ls = "-" if k == "data" else "--"
    axs[0].loglog(*v["sO"], color=cols[k], lw=lw, ls=ls, label=k); axs[1].loglog(*v["sC"], color=cols[k], lw=lw, ls=ls, label=k)
    axs[2].semilogx(LAGS, v["oo"], "o-" if k == "data" else "s--", color=cols[k], lw=lw, ms=4, label=f"{k} O-O")
    axs[2].semilogx(LAGS, v["cc"], "o:" if k == "data" else "s:", color=cols[k], lw=lw, ms=3, alpha=0.6)
    c = 0.5 * (v["e"][1:] + v["e"][:-1]); axs[3].step(c, v["h"], where="mid", color=cols[k], lw=lw, ls=ls, label=k)
    t = TLAGS * DT
    axs[4].semilogx(t, v["ck_mean"][0], color=cols[k], lw=lw, ls=ls, label=k)
    axs[5].semilogx(t, v["ck_mean"][1], color=cols[k], lw=lw, ls=ls, label=k)
    if k == "data":
        for ax, q in ((axs[4], 0), (axs[5], 1)):
            ax.fill_between(t, v["ck_mean"][q] - v["ck_sd"][q], v["ck_mean"][q] + v["ck_sd"][q], color="k", alpha=0.12)
axs[0].set(xlabel="open dwell (s)", ylabel="P(t >= tau)", title="open dwell survival"); axs[0].legend(fontsize=8)
axs[1].set(xlabel="closed dwell (s)", ylabel="P(t >= tau)", title="closed dwell survival")
axs[2].axhline(0, color="k", lw=0.5); axs[2].set(xlabel="lag (dwells)", ylabel="Spearman r", title="dwell-lag correlation (solid O-O, dotted C-C)"); axs[2].legend(fontsize=7)
axs[3].set(xlabel="transitions / s in 2-min window", ylabel="fraction of windows", title="local rate histogram (modal gating)")
axs[4].set(xlabel="t (s)", ylabel="D (nats)", title="Chapman-Kolmogorov divergence, q0 = open, tau = 1 s"); axs[4].axhline(0, color="k", lw=0.5)
axs[5].set(xlabel="t (s)", ylabel="D (nats)", title="Chapman-Kolmogorov divergence, q0 = closed, tau = 1 s"); axs[5].axhline(0, color="k", lw=0.5); axs[5].legend(fontsize=8)
fig.suptitle("Predictive checks: data (black, band = sd over 4 traces) vs simulations from fitted aggregated Markov models (4 x 6 h each)")
fig.tight_layout(); fig.savefig("figures/predictive_checks.png", dpi=140); print("wrote figures/predictive_checks.png")
out = {k: dict(oo=list(map(float, v["oo"])), cc=list(map(float, v["cc"])), oc=float(v["oc"]), rate_hist=list(map(float, v["h"])),
               ck_open=list(map(float, v["ck_mean"][0])), ck_closed=list(map(float, v["ck_mean"][1]))) for k, v in S.items()}
out["lags"] = LAGS; out["t_ck"] = list(map(float, TLAGS * DT))
json.dump(out, open("results/predictive_checks.json", "w"), indent=1)
