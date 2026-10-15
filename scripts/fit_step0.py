"""Audit step 0, part 1: non-exponential mode lifetimes. Each variant replaces one gating mode of the final skeleton
(ringF1: 0 hub, 1 slow, 2 open-biased, 3 fast, 4 long frozen X) by k identical sub-states (amm.split_mode):
  eq  : equilibrium trap chain (symmetric hopping; lifetime a heavier-tailed mixture of exponentials)
  drv : driven forward chain (enter first copy, leave from last; Erlang-like, more regular lifetime; no detailed balance)
with or without the equilibrium OU modulation of kOC. ML on a 9-point grid from the final-model optimum and several
hopping-rate starts, best refined on 15 points. Question: does a different lifetime shape fit better, and does it absorb
the OU gain (sigma -> smaller)?

Usage: python scripts/fit_step0.py --variants 1:eq 1:drv 2:eq ... [--k 2] [--no-ou] [--workers 4]
Output: results/step0/fit_<name>.json, results/step0/fit_step0.log
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from hsp90memory import amm

OUT = Path("results/step0")


def make(mode, k, driven, ou, grid):
    sm = amm.split_mode(amm.named_graph("ringF1"), mode, k, driven)
    return amm.modal_ou(0, n_grid=grid, couple="eq", base=sm) if ou else sm


def x0(mode, k, driven, ou, logh):
    if ou:
        th = np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"])
        return np.r_[th[:24], logh, th[24:]]
    th = np.array(json.load(open("results/exact/fit_ringF1.json"))["theta"])
    return np.r_[th, logh]


def run(args):
    mode, k, driven, ou, logh, grid, maxiter = args
    m = make(mode, k, driven, ou, grid)
    dl = [amm.load_dwells(f"data/{i}Thermal.txt")[1:] for i in "1234"]
    f = lambda th: (lambda ll: -ll if np.isfinite(ll) else 1e12)(amm.loglik_Q(m.Q(th), m.classes, dl))
    t0 = time.time(); start = x0(mode, k, driven, ou, logh) if np.isscalar(logh) else np.asarray(logh)
    r = minimize(f, start, method="L-BFGS-B", options=dict(maxiter=maxiter))
    return dict(mode=mode, k=k, driven=driven, ou=ou, grid=grid, start_logh=float(start[24]), ll=float(-r.fun),
                theta=r.x.tolist(), nit=int(r.nit), time=time.time() - t0)


def name(mode, k, driven, ou):
    return f"ringF1{'_oueq' if ou else ''}_split{mode}x{k}{'drv' if driven else 'eq'}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="+", required=True, help="mode:eq|drv")
    ap.add_argument("--k", type=int, default=2); ap.add_argument("--no-ou", action="store_true")
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--maxiter", type=int, default=300)
    a = ap.parse_args(); ou = not a.no_ou
    OUT.mkdir(parents=True, exist_ok=True); logf = open(OUT / "fit_step0.log", "a")
    def log(msg): print(msg, flush=True); logf.write(msg + "\n"); logf.flush()
    specs = [(int(v.split(":")[0]), v.split(":")[1] == "drv") for v in a.variants]
    specs = [s for s in specs if not (OUT / f"fit_{name(s[0], a.k, s[1], ou)}.json").exists()]
    hstarts = {False: [np.log(1 / 300.), np.log(1 / 50.)], True: [np.log(1 / 100.), np.log(1 / 20.)]}
    jobs = [(m, a.k, d, ou, h, 9, a.maxiter) for m, d in specs for h in hstarts[d]]
    t0 = time.time()
    with get_context("fork").Pool(a.workers) as pool:
        res = pool.map(run, jobs, chunksize=1)
        for r in res:
            log(f"{name(r['mode'], r['k'], r['driven'], r['ou'])} g9 start h=1/{np.exp(-r['start_logh']):.0f}s: ll {r['ll']:.3f} nit {r['nit']} ({r['time']:.0f} s)")
        best = {}
        for r in res:
            key = (r["mode"], r["driven"])
            if key not in best or r["ll"] > best[key]["ll"]: best[key] = r
        fine = pool.map(run, [(m, a.k, d, ou, np.array(b["theta"]), 15 if ou else 9, a.maxiter) for (m, d), b in best.items()], chunksize=1)
    for (key, b), r in zip(best.items(), fine):
        mm = make(r["mode"], a.k, r["driven"], ou, 15 if ou else 9); th = np.array(r["theta"]); nb = 25
        js = dict(name=name(r["mode"], a.k, r["driven"], ou), mode=r["mode"], k=a.k, driven=r["driven"], ou=ou,
                  n_par=mm.n_par, ll=r["ll"], grid=r["grid"], ll_g9=b["ll"], theta=th.tolist(), h=float(np.exp(th[24])),
                  sigma_OC=float(np.exp(th[nb])) if ou else 0.0, tau_w=float(np.exp(th[nb + 1])) if ou else None,
                  starts=[{k_: x[k_] for k_ in ("start_logh", "ll", "nit")} for x in res if (x["mode"], x["driven"]) == key])
        json.dump(js, open(OUT / f"fit_{js['name']}.json", "w"), indent=1)
        log(f"done {js['name']}: ll {js['ll']:.3f} (grid {js['grid']}), h 1/{1 / js['h']:.0f} s"
            + (f", sigma {js['sigma_OC']:.3f}, tau {js['tau_w']:.0f} s" if ou else ""))
    log(f"total {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
