"""Posterior-predictive checks of the best model (ringF1 + equilibrium OU on kOC) and of its discrete skeleton (ringF1)
for the extended-abstract figure.

For each model: N draws from the exact-likelihood MCMC posterior (results/exact/thin/bayes_<tag>_thin.npz), one
simulated measurement (4 x 6 h at 10 Hz, exact sampling of P = expm(Q dt), grid-OU on 15 points) per draw. Statistics,
computed identically for data and simulations:
  - survival of interior open and closed dwells on a log time grid;
  - Spearman correlation between open dwells i and i + lag (summaries.LAGS);
  - histogram of the local transition rate in 1-min windows (summaries.local_rate_hist);
  - Chapman-Kolmogorov divergence D_{1 s, q0}(t) for q0 = open, closed (Vollmar et al. 2024; summaries.ck_profile);
  - autocorrelation (Spearman, pooled over traces) of the local closing rate k = (# O->C) / (open time) in 1-min
    windows, lags 1-20 min.
Output: results/audit/ppc_paper.npz (settings embedded). Usage: python scripts/ppc_paper_figure.py [--n 200]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, json, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np, numba
from scipy.linalg import expm
from scipy.stats import spearmanr
from hsp90memory import amm, io, summaries as sm

T_GRID = np.logspace(-1, 3.3, 44)
ACF_LAGS_MIN = np.array([1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20])
W = 600
MODELS = {"ringF1_oueq": "ringF1 + OU (eq.)", "ringF1": "ringF1"}


@numba.njit(cache=True)
def _sim(C, x0, u):
    out = np.empty(u.size, np.int64); x = x0
    for t in range(u.size):
        out[t] = x; x = np.searchsorted(C[x], u[t])
        if x >= C.shape[0]: x = C.shape[0] - 1
    return out


def survival(series, cls):
    d = np.concatenate([(lambda dw, dc: dw[1:-1][dc[1:-1] == cls])(*io.dwells(s)) for s in series]) * io.DT
    d = np.sort(d)
    return 1.0 - np.searchsorted(d, T_GRID - 1e-9, side="left") / d.size


def closing_rate_acf(series):
    ks = []
    for s in series:
        n = s.size // W; x = s[:n * W].reshape(n, W).astype(int)
        c = np.sum(np.diff(x, axis=1) == 1, 1) + np.r_[0, (x[1:, 0] - x[:-1, -1]) == 1]
        t = np.sum(x == 0, 1) * io.DT
        ks.append(np.where(t > 0, c / np.maximum(t, 1e-9), np.nan))
    out = []
    for L in ACF_LAGS_MIN:
        a = np.concatenate([k[:-L] for k in ks]); b = np.concatenate([k[L:] for k in ks]); ok = np.isfinite(a) & np.isfinite(b)
        out.append(spearmanr(a[ok], b[ok])[0])
    return np.array(out)


def stats(series):
    O, C = sm.dwell_pairs(series); oo = sm.lag_correlations(O, C)[0]
    return dict(sO=survival(series, 0), sC=survival(series, 1), oo=oo, rate=sm.local_rate_hist(series, W),
                ck=sm.ck_profile(series), kacf=closing_rate_acf(series))


def model_of(tag):
    if tag.endswith("_oueq"):
        m = amm.modal_ou(5, n_grid=15, couple="eq", topology="ringF1"); nb = m.base.n_par
        return m, (lambda phi: np.r_[phi[:nb], np.log(phi[nb]), phi[nb + 1]])
    m = amm.named_graph("ringF1"); return m, (lambda phi: np.asarray(phi, float))


def one(args):
    tag, phi, seed = args
    m, to_theta = model_of(tag); Q = m.Q(to_theta(phi))
    P = expm(Q * io.DT); C = np.cumsum(P, 1); C[:, -1] = 1.0; pi = amm.stationary(Q); cls = np.asarray(m.classes)
    rng = np.random.default_rng(seed)
    series = [cls[_sim(C, rng.choice(pi.size, p=pi), rng.random(216000))].astype(np.int8) for _ in range(4)]
    return stats(series)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=200); ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args(); t0 = time.time()
    data = stats([io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in "1234"])
    out = {f"data_{k}": v for k, v in data.items()}
    rng = np.random.default_rng(2026)
    for tag in MODELS:
        phi = np.load(f"results/exact/thin/bayes_{tag}_thin.npz")["phi"]
        idx = rng.choice(len(phi), a.n, replace=False)
        jobs = [(tag, phi[i], 1000 + j) for j, i in enumerate(idx)]
        with get_context("fork").Pool(a.workers) as pool:
            R = pool.map(one, jobs, chunksize=4)
        for k in R[0]:
            out[f"{tag}_{k}"] = np.array([r[k] for r in R])
        out[f"{tag}_draws"] = idx
        print(tag, f"{time.time() - t0:.0f} s", flush=True)
    out["t_grid"] = T_GRID; out["oo_lags"] = sm.LAGS; out["rate_edges"] = sm.RATE_EDGES; out["ck_lags_s"] = sm.CK_LAGS * io.DT
    out["kacf_lags_min"] = ACF_LAGS_MIN
    out["settings"] = json.dumps(dict(n=a.n, models=MODELS, grid=15, window_s=W * io.DT, seed=2026,
                                      posterior="results/exact/thin/bayes_<tag>_thin.npz"))
    np.savez_compressed("results/audit/ppc_paper.npz", **out)
    print("done", f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
