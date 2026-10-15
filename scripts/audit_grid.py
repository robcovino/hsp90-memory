"""Audit: is the hidden modulator of the closing rate resolved as continuous? Refit ringF1 + OU (eq coupling) with the
z coordinate on n_grid = 2, 3, 5 points (same 26 parameters; n_grid = 2 is a symmetric two-level switch that multiplies
every closing rate by exp(+-sigma) and flips with relaxation time tau_w). Starts: the grid-15 optimum and two
(sigma, tau_w) variants. Compare max log L with grid 15 (-97918.86) and ringF1 (-97996.32).
Output: results/audit/grid_switch.json (settings embedded), log results/audit/grid_switch.log.
Usage: python scripts/audit_grid.py [--grids 2 3 5]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from hsp90memory import amm

OUT = Path("results/audit")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--grids", type=int, nargs="*", default=[2, 3, 5])
    ap.add_argument("--maxiter", type=int, default=300); a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); logf = open(OUT / "grid_switch.log", "a")
    def log(s): print(s, flush=True); logf.write(s + "\n"); logf.flush()
    dl = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
    th0 = np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"]); nb = 24
    bounds = [(-20, 5)] * nb + [(np.log(1e-3), np.log(3.0)), (0.0, np.log(1e5))]
    starts = {"g15 optimum": th0,
              "sigma x1.6": np.r_[th0[:nb], th0[nb] + np.log(1.6), th0[nb + 1]],
              "sigma x0.7, tau x0.4": np.r_[th0[:nb], th0[nb] + np.log(0.7), th0[nb + 1] + np.log(0.4)]}
    res = {"settings": vars(a), "ref": {"ringF1_oueq_g15": -97918.856, "ringF1": -97996.319}, "grids": {}}
    for g in a.grids:
        m = amm.modal_ou(5, n_grid=g, couple="eq", topology="ringF1")
        f = lambda th: (lambda v: -v if np.isfinite(v) else 1e12)(amm.loglik_Q(m.Q(th), m.classes, dl))
        out = []
        for lab, x0 in starts.items():
            t0 = time.time()
            r = minimize(f, x0, method="L-BFGS-B", bounds=bounds, options=dict(maxiter=a.maxiter))
            out.append(dict(start=lab, ll=float(-r.fun), sigma=float(np.exp(r.x[nb])), tau_w=float(np.exp(r.x[nb + 1])),
                            nit=int(r.nit), theta=r.x.tolist(), time_s=time.time() - t0))
            log(f"grid {g} start '{lab}': ll={-r.fun:.3f} sigma={np.exp(r.x[nb]):.3f} tau={np.exp(r.x[nb+1]):.1f} "
                f"nit={r.nit} {time.time()-t0:.0f}s")
        res["grids"][str(g)] = dict(best=max(out, key=lambda d: d["ll"]), starts=out, z=m.z.tolist())
        json.dump(res, open(OUT / "grid_switch.json", "w"), indent=1)


if __name__ == "__main__":
    main()
