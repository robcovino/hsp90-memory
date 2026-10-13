"""Leave-one-trace-out cross-validation of the exact-likelihood models: fit by ML on three traces (started from the
pooled ML parameters), evaluate the log-likelihood of the held-out trace. Out-of-sample predictive comparison,
no penalty for parameter count needed.

Usage: python scripts/exact_cv.py --models modal4 modal5 modal6 modal4_ou modal4_oueq [--grid 9] [--workers 4]
Output: results/exact/cv.json
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, re, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from hsp90memory import amm
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent)); from _fits import best_fit_file  # noqa: E402


def model_and_start(name, grid):
    m = re.match(r"(modal|ring)(\d+)(_noDB)?(_ou(eq|2)?F?)?$", name); K = int(m.group(2))
    topo = "chain" if m.group(1) == "modal" else ("ring_noDB" if m.group(3) else "ring")
    if m.group(4):
        couple = {None: "OC", "eq": "eq", "2": "both"}[m.group(5)]
        return amm.modal_ou(K, n_grid=grid, couple=couple, topology=topo), np.array(json.load(open(best_fit_file(name)))["theta"])
    return (amm.modal(K) if topo == "chain" else amm.ring(K, db=(topo == "ring"))), np.array(json.load(open(best_fit_file(name)))["theta"])


def fold(args):
    name, grid, k = args
    dl = [amm.load_dwells(f"data/{i}Thermal.txt")[1:] for i in "1234"]
    train = [d for i, d in enumerate(dl) if i != k]; test = [dl[k]]
    model, x0 = model_and_start(name, grid); t0 = time.time()
    f = lambda th: (lambda ll: -ll if np.isfinite(ll) else 1e12)(amm.loglik_Q(model.Q(th), model.classes, train))
    r = minimize(f, x0, method="L-BFGS-B", options=dict(maxiter=300))
    ll_test = amm.loglik_Q(model.Q(r.x), model.classes, test)
    return dict(model=name, fold=k, ll_train=float(-r.fun), ll_test=float(ll_test), nit=int(r.nit), time=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=["modal4", "modal5", "modal6", "modal4_ou", "modal4_oueq"])
    ap.add_argument("--grid", type=int, default=9); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default="results/exact/cv.json")
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists")
    jobs = [(m, a.grid, k) for m in a.models for k in range(4)]
    with get_context("fork").Pool(a.workers) as pool:
        res = pool.map(fold, jobs, chunksize=1)
    summ = {}
    for m in a.models:
        r = [x for x in res if x["model"] == m]
        summ[m] = {"ll_test_total": sum(x["ll_test"] for x in r), "ll_test_per_trace": [x["ll_test"] for x in sorted(r, key=lambda x: x["fold"])]}
    ref = "modal5" if "modal5" in summ else a.models[0]
    for m in summ:
        summ[m]["delta_vs_" + ref] = summ[m]["ll_test_total"] - summ[ref]["ll_test_total"]
        summ[m]["delta_per_trace_vs_" + ref] = list(np.array(summ[m]["ll_test_per_trace"]) - np.array(summ[ref]["ll_test_per_trace"]))
    json.dump({"grid": a.grid, "models": a.models, "folds": res, "summary": summ}, open(out, "w"), indent=1)
    for m, s in summ.items():
        print(m, round(s["ll_test_total"], 2), "delta", round(s["delta_vs_" + ref], 2), np.round(s["delta_per_trace_vs_" + ref], 2))


if __name__ == "__main__":
    main()
