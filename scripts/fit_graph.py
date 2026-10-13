"""ML fits of the discrete ring-skeleton competitors (amm.named_graph: ringL1, ringL2, ringS2), i.e. the ring of four
gating modes with a discrete kOC ladder branching off the fast mode. Starts: the ring4 optimum
(results/exact/fit_ring4.json) plus ladder modes (open times 2-15 s, closed time of the fast mode, entered from the open
state), several open-time choices, each followed by --hops basin-hopping rounds.
Output: results/exact/fit_<name>.json. Usage: python scripts/fit_graph.py ringL1 ringL2 ringS2 [--hops 3]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from hsp90memory import amm


def starts(name):
    t4 = np.array(json.load(open("results/exact/fit_ring4.json"))["theta"]); modes = t4[:8]; sw4 = t4[8:]
    n, edges = amm.NAMED_GRAPHS[name]; kco_fast = t4[7]; out = []
    for opens in ([2.0, 8.0], [3.0, 12.0], [1.8, 5.0], [2.5, 20.0]):
        lm = np.concatenate([[np.log(1 / o), kco_fast] for o in opens[:n - 4]])
        extra = np.concatenate([np.log([1 / 200., 1 / 60., 1e-4]) for _ in range(n - 4)])
        out.append((f"ring4 + ladder {opens[:n - 4]}", np.r_[modes, lm, sw4, extra]))
    return out


def run(args):
    name, label, x0, hops, seed = args
    m = amm.named_graph(name); dl = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
    f = lambda th: (lambda ll: -ll if np.isfinite(ll) else 1e12)(amm.loglik_Q(m.Q(th), m.classes, dl))
    rng = np.random.default_rng(seed); t0 = time.time()
    best = minimize(f, x0, method="L-BFGS-B", options=dict(maxiter=800)); hist = [float(-best.fun)]
    for h in range(hops):
        r = minimize(f, best.x + rng.normal(0, 0.3, len(x0)), method="L-BFGS-B", options=dict(maxiter=800))
        hist.append(float(-r.fun)); best = r if r.fun < best.fun else best
    print(f"{name} {label}: {hist} ({time.time() - t0:.0f} s)", flush=True)
    return dict(start=label, ll=float(-best.fun), theta=best.x.tolist(), hops=hist)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("names", nargs="+"); ap.add_argument("--hops", type=int, default=3)
    ap.add_argument("--workers", type=int, default=4); a = ap.parse_args()
    n_dw = sum(len(amm.load_dwells(f"data/{k}Thermal.txt")[1]) for k in "1234")
    for name in a.names:
        out = Path(f"results/exact/fit_{name}.json")
        if out.exists(): print(out, "exists"); continue
        S = starts(name)
        with get_context("fork").Pool(a.workers) as pool:
            res = pool.map(run, [(name, l, x, a.hops, 7 + i) for i, (l, x) in enumerate(S)], chunksize=1)
        b = max(res, key=lambda r: r["ll"]); m = amm.named_graph(name); Q = m.Q(np.array(b["theta"]))
        json.dump(dict(model=name, edges=m.edges, n_par=m.n_par, theta=b["theta"], ll=b["ll"], n_dwells=n_dw,
                       bic=float(amm.bic(b["ll"], m.n_par, n_dw)), candidates=sorted(res, key=lambda r: -r["ll"]),
                       Q=Q.tolist(), pi=amm.stationary(Q).tolist(), classes=list(map(int, m.classes))),
                  open(out, "w"), indent=1)
        print(f"{name}: best {b['ll']:.3f} ({b['start']}), BIC {amm.bic(b['ll'], m.n_par, n_dw):.1f}", flush=True)


if __name__ == "__main__":
    main()
