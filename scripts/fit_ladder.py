"""Global maximum-likelihood search for discrete modal-K AMMs (K >= 5), where L-BFGS from a single family of starts
can stop in a local optimum (modal6 from modal5 splits: -98087.8; MCMC later found -98060.5).

Candidates for modal-K:
  - existing fits (results/modalK.json, results/exact/fit_modalK.json) and the best MCMC sample
    (results/exact/bayes_modalK_samples.npz);
  - every split of the best (K-1)-mode solution (rates x/÷1.5, copies coupled at 1/300 s);
  - every deletion of one mode from the best (K+1)-mode solution, if available;
  - graded-kOC ladders: modal4 ML plus K-4 modes appended after the fast mode with the fast kCO, kOC graded
    (open times 2-15 s), coupled from the open state.
Each candidate: L-BFGS (exact likelihood amm.loglik_Q), then --hops basin-hopping rounds (log-rate perturbation sd 0.3,
re-optimise, keep if better). Best -> results/exact/fit_modalK_global.json (all candidates listed).

Usage: python scripts/fit_ladder.py 5 6 7 [--hops 3] [--workers 4]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, sys, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from hsp90memory import amm

E = Path("results/exact")


def split_theta(th, K):
    modes = [th[2 * i:2 * i + 2] for i in range(K)]; coup = [th[2 * K + 3 * i:2 * K + 3 * i + 3] for i in range(K - 1)]
    return modes, coup


def join(modes, coup):
    return np.concatenate(list(modes) + list(coup))


def best_theta(K):
    """Best available solution for modal-K: *_global, then fits, then MCMC samples."""
    cands = []
    for p in (E / f"fit_modal{K}_global.json", Path(f"results/modal{K}.json"), E / f"fit_modal{K}.json"):
        if p.exists():
            d = json.load(open(p)); cands.append((d["ll"], np.array(d["theta"]), str(p)))
    p = E / f"bayes_modal{K}_samples.npz"
    if p.exists():
        z = np.load(p); i = np.unravel_index(np.argmax(z["loglik"]), z["loglik"].shape)
        cands.append((float(z["loglik"][i]), np.array(z["phi"][i]), str(p)))
    return max(cands, key=lambda c: c[0]) if cands else None


def candidates(K):
    out = []
    for p in (Path(f"results/modal{K}.json"), E / f"fit_modal{K}.json"):
        if p.exists(): out.append((str(p), np.array(json.load(open(p))["theta"])))
    p = E / f"bayes_modal{K}_samples.npz"
    if p.exists():
        z = np.load(p); L = z["loglik"]
        for c in range(L.shape[0]):
            i = np.argmax(L[c]); out.append((f"mcmc chain {c}", np.array(z["phi"][c, i])))
    b = best_theta(K - 1)
    if b is not None:
        modes, coup = split_theta(b[1], K - 1); f = np.log(1.5)
        for j in range(K - 1):
            out.append((f"split{j} of modal{K - 1}", join(modes[:j] + [modes[j] + [f, -f], modes[j] + [-f, f]] + modes[j + 1:],
                                                     coup[:j] + [np.full(3, np.log(1 / 300.))] + coup[j:])))
    b = best_theta(K + 1)
    if b is not None:
        modes, coup = split_theta(b[1], K + 1)
        for j in range(K + 1):
            m2 = modes[:j] + modes[j + 1:]
            c2 = coup[:j] + coup[j + 1:] if j < K else coup[:j - 1] + coup[j:]
            out.append((f"delete{j} of modal{K + 1}", join(m2, c2[:K - 1])))
    th4 = np.array(json.load(open("results/modal4.json"))["theta"]); modes, coup = split_theta(th4, 4)
    kco_fast = modes[3][1]
    for opens in ([2.0, 6.0, 15.0], [3.0, 10.0, 30.0], [1.8, 4.0, 10.0]):
        extra = [np.array([np.log(1 / o), kco_fast]) for o in opens[:K - 4]]
        c_extra = [np.log([1 / 150., 1 / 50., 1e-4]) for _ in range(K - 4)]
        out.append((f"ladder {opens[:K - 4]}", join(modes + extra, coup + c_extra)))
    return out


def run(args):
    K, name, x0, hops, seed = args
    m = amm.modal(K); dl = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
    f = lambda th: (lambda ll: -ll if np.isfinite(ll) else 1e12)(amm.loglik_Q(m.Q(th), m.classes, dl))
    t0 = time.time(); rng = np.random.default_rng(seed)
    x0 = np.clip(x0, -20, 6)
    r = minimize(f, x0, method="L-BFGS-B", options=dict(maxiter=800)); best = r; hist = [float(-r.fun)]
    for h in range(hops):
        r2 = minimize(f, best.x + rng.normal(0, 0.3, len(x0)), method="L-BFGS-B", options=dict(maxiter=800))
        hist.append(float(-r2.fun))
        if r2.fun < best.fun: best = r2
    print(f"modal{K} {name}: start ll {-f(x0):.2f} -> {hist} ({time.time() - t0:.0f} s)", flush=True)
    return dict(K=K, start=name, ll=float(-best.fun), theta=best.x.tolist(), hops=hist, time=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("K", type=int, nargs="+")
    ap.add_argument("--hops", type=int, default=3); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    n_dw = sum(len(amm.load_dwells(f"data/{k}Thermal.txt")[1]) for k in "1234")
    for K in a.K:
        out = E / f"fit_modal{K}_global.json"
        if out.exists():
            print(f"{out} exists, skipped"); continue
        C = candidates(K); t0 = time.time()
        with get_context("fork").Pool(a.workers) as pool:
            res = pool.map(run, [(K, n, x, a.hops, 100 + i) for i, (n, x) in enumerate(C)], chunksize=1)
        best = max(res, key=lambda r: r["ll"]); m = amm.modal(K); th = np.array(best["theta"]); Q = m.Q(th)
        js = dict(model=f"modal{K}", n_par=m.n_par, theta=th.tolist(), rates=np.exp(th).tolist(), ll=best["ll"],
                  n_dwells=n_dw, bic=float(amm.bic(best["ll"], m.n_par, n_dw)), best_start=best["start"],
                  candidates=[{k: r[k] for k in ("start", "ll", "hops", "time")} for r in sorted(res, key=lambda r: -r["ll"])],
                  classes=list(map(int, m.classes)), Q=Q.tolist(), pi=amm.stationary(Q).tolist(), time=time.time() - t0)
        json.dump(js, open(out, "w"), indent=1)
        print(f"modal{K} global best {best['ll']:.3f} from {best['start']}; BIC {js['bic']:.1f}", flush=True)


if __name__ == "__main__":
    main()
