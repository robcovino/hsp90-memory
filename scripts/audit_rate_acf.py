"""Audit: a direct look at the slow modulation. Local closing rate k_w = (# O->C transitions) / (open time) in 1-min
windows; Spearman autocorrelation of k_w at lags 1-40 min (pooled over traces), and the variance of log local rate
(windows with >= 3 closings). Data vs 40 simulated 4 x 6 h data sets from ringF1 and from ringF1 + OU eq. (grid 15),
each at its ML parameters (plug-in predictive check).
Output: results/audit/rate_acf.json.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
from pathlib import Path
import numpy as np, numba
from scipy.linalg import expm
from scipy.stats import spearmanr
from hsp90memory import amm, io

W = 600; LAGS = [1, 2, 3, 5, 8, 12, 20, 30, 40]; NREP = 40


@numba.njit(cache=True)
def sim(C, x0, u):
    out = np.empty(u.size, np.int64); x = x0
    for t in range(u.size):
        out[t] = x; x = np.searchsorted(C[x], u[t])
        if x >= C.shape[0]: x = C.shape[0] - 1
    return out


def window_rates(s):
    n = s.size // W; s = s[:n * W].reshape(n, W).astype(int)
    closings = np.sum(np.diff(s, axis=1) == 1, 1) + np.r_[0, (s[1:, 0] - s[:-1, -1]) == 1]
    t_open = np.sum(s == 0, 1) * io.DT
    return closings, t_open


def stats(series):
    acf = {L: [] for L in LAGS}; lv = []
    for s in series:
        c, t = window_rates(s); k = np.where(t > 0, c / np.maximum(t, 1e-9), np.nan)
        for L in LAGS:
            a, b = k[:-L], k[L:]; ok = np.isfinite(a) & np.isfinite(b); acf[L].append((a[ok], b[ok]))
        lv.append(np.log(c[c >= 3] / t[c >= 3]))
    out = {f"acf_{L}min": float(spearmanr(np.concatenate([p[0] for p in acf[L]]), np.concatenate([p[1] for p in acf[L]]))[0])
           for L in LAGS}
    out["var_log_k"] = float(np.var(np.concatenate(lv)))
    return out


def main():
    S = [io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in "1234"]
    dat = stats(S); print("data", {k: round(v, 3) for k, v in dat.items()}, flush=True)
    th1 = np.array(json.load(open("results/exact/fit_ringF1.json"))["theta"])
    th2 = np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"])
    models = {"ringF1": (amm.named_graph("ringF1"), th1),
              "ringF1_oueq": (amm.modal_ou(5, n_grid=15, couple="eq", topology="ringF1"), th2)}
    res = {"settings": dict(window_s=W * io.DT, lags_min=LAGS, n_rep=NREP), "data": dat, "sims": {}}
    rng = np.random.default_rng(23)
    for name, (m, th) in models.items():
        Q = m.Q(th); P = expm(Q * io.DT); C = np.cumsum(P, 1); C[:, -1] = 1.0; pi = amm.stationary(Q); cls = np.asarray(m.classes)
        R = []
        for r in range(NREP):
            ser = [cls[sim(C, rng.choice(pi.size, p=pi), rng.random(216000))] for _ in range(4)]
            R.append(stats(ser))
        res["sims"][name] = R
        print(name)
        for k in dat:
            v = np.array([x[k] for x in R])
            print(f"  {k:10s} data {dat[k]:7.3f}  sims {np.median(v):7.3f} [{np.percentile(v, 2.5):.3f}, {np.percentile(v, 97.5):.3f}]"
                  f"  z {(dat[k] - v.mean()) / v.std():5.2f}", flush=True)
    json.dump(res, open("results/audit/rate_acf.json", "w"), indent=1)


if __name__ == "__main__":
    main()
