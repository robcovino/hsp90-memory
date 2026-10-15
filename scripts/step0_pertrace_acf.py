"""Audit step 0, part 4: does molecule-to-molecule heterogeneity of the modulation explain the 5-10 min misfit of the
local closing rate? Gating rates and tau_w fixed at the final model (ringF1 + equilibrium OU, 15-point grid); sigma_OC
re-fitted per trace (1D maximum likelihood); variants
  shared   : the final model for every trace
  sigma4   : each trace simulated with its own sigma_OC
  off1     : trace 1 without modulation (its discrete ringF1 rates = the final-model rates), traces 2-4 shared sigma
Statistic: as scripts/audit_rate_acf.py (1-min windows, Spearman ACF of the local closing rate), pooled over the four
traces and per trace; 40 simulated measurements per variant.
Output: results/step0/pertrace_acf.json
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, sys
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_rate_acf import sim, stats, window_rates, LAGS                # noqa: E402
from hsp90memory import amm, io                                          # noqa: E402

NREP = 40
TH = np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"])
M = amm.modal_ou(5, n_grid=15, couple="eq", topology="ringF1"); G = amm.named_graph("ringF1")


def Q_sigma(s):
    return G.Q(TH[:24]) if s == 0 else M.Q(np.r_[TH[:24], np.log(s), TH[25]])


def per_trace_acf(s):
    c, t = window_rates(s); k = np.where(t > 0, c / np.maximum(t, 1e-9), np.nan); out = {}
    for L in LAGS:
        a, b = k[:-L], k[L:]; ok = np.isfinite(a) & np.isfinite(b); out[f"acf_{L}min"] = float(spearmanr(a[ok], b[ok])[0])
    return out


def one(args):
    sig, r = args
    rng = np.random.default_rng(5000 + r); ser = []
    for s in sig:
        Q = Q_sigma(s); cls = np.asarray(G.classes if s == 0 else M.classes)
        P = expm(Q * io.DT); C = np.cumsum(P, 1); C[:, -1] = 1.0; pi = amm.stationary(Q)
        ser.append(cls[sim(C, rng.choice(pi.size, p=pi), rng.random(216000))])
    return stats(ser), [per_trace_acf(x) for x in ser]


def main():
    data = [io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in "1234"]
    dl = [io.dwells(x) for x in data]
    sig_k = []
    for k in range(4):
        f = lambda ls: -amm.loglik_Q(Q_sigma(np.exp(ls)), M.classes, [dl[k]])
        r = minimize_scalar(f, bounds=(np.log(0.01), np.log(1.0)), method="bounded", options=dict(xatol=1e-3))
        ll0 = amm.loglik_Q(G.Q(TH[:24]), G.classes, [dl[k]])
        sig_k.append(float(np.exp(r.x))); print(f"trace {k + 1}: sigma {np.exp(r.x):.3f}, gain vs no modulation {-r.fun - ll0:+.1f}", flush=True)
    variants = {"shared": [float(np.exp(TH[24]))] * 4, "sigma4": sig_k, "off1": [0.0] + [float(np.exp(TH[24]))] * 3}
    dat_pool = stats(data); dat_tr = [per_trace_acf(x) for x in data]
    res = {"sigma_per_trace": sig_k, "data_pooled": dat_pool, "data_per_trace": dat_tr, "variants": {}}
    with get_context("fork").Pool(4) as pool:
        for name, sig in variants.items():
            R = pool.map(one, [(sig, r) for r in range(NREP)], chunksize=2)
            pooled = {k: [x[0][k] for x in R] for k in dat_pool}
            zp = {k: float((dat_pool[k] - np.mean(v)) / np.std(v)) for k, v in pooled.items()}
            ztr = [{k: float((dat_tr[i][k] - np.mean([x[1][i][k] for x in R])) / np.std([x[1][i][k] for x in R]))
                    for k in dat_tr[i]} for i in range(4)]
            med = {k: float(np.median(v)) for k, v in pooled.items()}
            res["variants"][name] = {"sigma": sig, "pooled_median": med, "pooled_z": zp, "per_trace_z": ztr}
            print(f"{name:7s} pooled median acf 5/8 min {med['acf_5min']:.3f} {med['acf_8min']:.3f}  z 5/8 min "
                  f"{zp['acf_5min']:+.1f} {zp['acf_8min']:+.1f}  | per-trace z(8 min) {[round(z['acf_8min'], 1) for z in ztr]}", flush=True)
    print("data pooled acf 5/8 min", round(dat_pool["acf_5min"], 3), round(dat_pool["acf_8min"], 3),
          "| per trace 5 min", [round(d["acf_5min"], 2) for d in dat_tr], "8 min", [round(d["acf_8min"], 2) for d in dat_tr])
    json.dump(res, open("results/step0/pertrace_acf.json", "w"), indent=1)


if __name__ == "__main__":
    main()
