"""Audit: is the OU gain carried by the shortest events? Remove every interior 1-frame dwell (the frame takes the class
of its predecessor, merging it into the surrounding dwell), then refit ringF1 (from its optimum) and ringF1 + OU eq.
(grid 9, from the grid-15 optimum) and compare the gain with the unmodified data (grid 9: -97918.32 vs -97996.32).
Both models are equally misspecified for the filtered data (they still predict 1-frame dwells), so only the
difference is of interest. Output: results/audit/short_merge.json, log results/audit/short_merge.log.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from hsp90memory import amm, io

OUT = Path("results/audit")


def merge_singletons(s):
    s = s.copy(); dw, dc = io.dwells(s); end = np.cumsum(dw); start = end - dw
    for i in range(1, len(dw) - 1):
        if dw[i] == 1: s[start[i]] = s[start[i] - 1]
    return s


def main():
    logf = open(OUT / "short_merge.log", "a")
    def log(m): print(m, flush=True); logf.write(m + "\n"); logf.flush()
    S = [merge_singletons(io.idealise(io.load_trace(f"data/{k}Thermal.txt"))) for k in "1234"]
    dl = [io.dwells(s) for s in S]
    log(f"dwells after merging 1-frame events: {sum(len(d[0]) for d in dl)} (min interior dwell "
        f"{min(int(d[0][1:-1].min()) for d in dl)} frames)")
    th1 = np.array(json.load(open("results/exact/fit_ringF1.json"))["theta"])
    th2 = np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"])
    m1 = amm.named_graph("ringF1"); m2 = amm.modal_ou(5, n_grid=9, couple="eq", topology="ringF1")
    res = {}
    for name, m, x0, bounds in (("ringF1", m1, th1, None),
                                ("ringF1_oueq_g9", m2, th2, [(-20, 5)] * 24 + [(np.log(1e-3), np.log(3.0)), (0.0, np.log(1e5))])):
        f = lambda th: (lambda v: -v if np.isfinite(v) else 1e12)(amm.loglik_Q(m.Q(th), m.classes, dl))
        t0 = time.time(); r = minimize(f, x0, method="L-BFGS-B", bounds=bounds, options=dict(maxiter=400))
        res[name] = dict(ll_start=float(-f(x0)), ll=float(-r.fun), nit=int(r.nit), theta=r.x.tolist(), time_s=time.time() - t0)
        extra = f" sigma={np.exp(r.x[24]):.3f} tau={np.exp(r.x[25]):.1f}" if "oueq" in name else ""
        log(f"{name}: ll start {-f(x0):.2f} -> {-r.fun:.3f}{extra} nit={r.nit} {time.time()-t0:.0f}s")
    res["gain"] = res["ringF1_oueq_g9"]["ll"] - res["ringF1"]["ll"]; res["gain_unfiltered_g9"] = -97918.32 + 97996.32
    log(f"gain filtered {res['gain']:.1f} vs unfiltered (grid 9) {res['gain_unfiltered_g9']:.1f}")
    json.dump(res, open(OUT / "short_merge.json", "w"), indent=1)


if __name__ == "__main__":
    main()
