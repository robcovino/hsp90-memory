"""Maximum-likelihood fit of modal5 (5 gating modes, detailed balance) to the four pooled traces.
Seeds: the default initialisation, and modal4 (results/modal4.json) with each of its modes split in two
(rates x/÷ 1.5 in opposite directions, the two copies coupled at 1/300 s). L-BFGS from every seed in parallel; the best
fit goes to results/modal5.json in the format of scripts/run_fits.py, all seeds to results/modal5_seeds.json.

Usage: python scripts/fit_modal5.py [--workers 4]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, sys, time
from multiprocessing import get_context
import numpy as np
from scipy.optimize import minimize
from hsp90memory import amm

M5 = amm.modal(5)


def seeds_from_modal4():
    th4 = np.array(json.load(open("results/modal4.json"))["theta"])
    modes = [th4[2 * m:2 * m + 2] for m in range(4)]
    coup = [th4[8 + 3 * m:8 + 3 * m + 3] for m in range(3)]
    f = np.log(1.5); out = {}
    for j in range(4):
        new_modes = modes[:j] + [modes[j] + [f, -f], modes[j] + [-f, f]] + modes[j + 1:]
        new_coup = coup[:j] + [np.full(3, np.log(1 / 300.))] + coup[j:]
        out[f"split{j}"] = np.concatenate(new_modes + new_coup)
    return out


def run(args):
    name, x0, dl, maxiter = args
    t0 = time.time()
    res = minimize(lambda th: -amm.loglik(M5, th, dl), x0, method="L-BFGS-B", options=dict(maxiter=maxiter))
    print(f"{name}: ll={-res.fun:.2f} nit={res.nit} {time.time() - t0:.0f} s", flush=True)
    return name, res.x.tolist(), float(-res.fun), int(res.nit), time.time() - t0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--maxiter", type=int, default=400); a = ap.parse_args()
    if os.path.exists("results/modal5.json"):
        sys.exit("results/modal5.json exists")
    dl = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
    seeds = {"default": M5.init, **seeds_from_modal4()}
    t0 = time.time()
    with get_context("fork").Pool(a.workers) as pool:
        res = pool.map(run, [(k, v, dl, a.maxiter) for k, v in seeds.items()], chunksize=1)
    json.dump([dict(seed=r[0], theta=r[1], ll=r[2], nit=r[3], time=r[4]) for r in res],
              open("results/modal5_seeds.json", "w"), indent=1)
    best = max(res, key=lambda r: r[2]); th = np.array(best[1]); n = sum(len(x[0]) for x in dl)
    Q = M5.Q(th)
    out = dict(model="modal5", n_par=M5.n_par, theta=th.tolist(), rates=np.exp(th).tolist(), ll=best[2], n_dwells=n,
               bic=float(amm.bic(best[2], M5.n_par, n)), time=time.time() - t0, classes=list(map(int, M5.classes)),
               seed=best[0], Q=Q.tolist(), pi=amm.stationary(Q).tolist())
    json.dump(out, open("results/modal5.json", "w"), indent=1)
    print("modal5 best", best[0], "ll", best[2], "BIC", out["bic"], f"time {out['time']:.0f} s")
    print("Q (1/s):\n", np.round(Q, 5)); print("pi:", np.round(out["pi"], 4))


if __name__ == "__main__":
    main()
