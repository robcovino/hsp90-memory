"""Audit step 0, part 3: the 5-10 min misfit of the local closing rate. Same statistic and simulator as
scripts/audit_rate_acf.py (1-min windows, pooled Spearman ACF at 1-40 min, variance of log local rate; 40 simulated
4 x 6 h data sets at the ML parameters, plug-in), evaluated for the step-0 variants (results/step0/fit_*.json) next
to the final model. Output: results/step0/rate_acf.json.

Usage: python scripts/step0_rate_acf.py [--names ringF1_oueq_split3x2drv ...] [--nrep 40] [--workers 4]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, sys
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.linalg import expm
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_rate_acf import sim, stats                                  # noqa: E402
from hsp90memory import amm, io                                         # noqa: E402


def model_Q(name):
    if name == "ringF1_oueq":
        m = amm.modal_ou(5, n_grid=15, couple="eq", topology="ringF1")
        return m.Q(np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"])), m.classes
    d = json.load(open(f"results/step0/fit_{name}.json"))
    sm = amm.split_mode(amm.named_graph("ringF1"), d["mode"], d["k"], d["driven"])
    m = amm.modal_ou(0, n_grid=15, couple="eq", base=sm) if d["ou"] else sm
    return m.Q(np.array(d["theta"])), m.classes


def one(args):
    name, r = args
    Q, cls = model_Q(name); cls = np.asarray(cls)
    P = expm(Q * io.DT); C = np.cumsum(P, 1); C[:, -1] = 1.0; pi = amm.stationary(Q)
    rng = np.random.default_rng(1000 + r)
    return stats([cls[sim(C, rng.choice(pi.size, p=pi), rng.random(216000))] for _ in range(4)])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--names", nargs="*", default=None)
    ap.add_argument("--nrep", type=int, default=40); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    names = a.names or ["ringF1_oueq"] + sorted(p.stem[4:] for p in Path("results/step0").glob("fit_ringF1*.json"))
    dat = stats([io.idealise(io.load_trace(f"data/{k}Thermal.txt")) for k in "1234"])
    res = {"data": dat, "nrep": a.nrep, "models": {}}
    with get_context("fork").Pool(a.workers) as pool:
        for n in names:
            R = pool.map(one, [(n, r) for r in range(a.nrep)], chunksize=2)
            z = {k: float((dat[k] - np.mean([x[k] for x in R])) / np.std([x[k] for x in R])) for k in dat}
            med = {k: float(np.median([x[k] for x in R])) for k in dat}
            res["models"][n] = {"median": med, "z": z}
            print(n, "| median acf 1/5/8/12 min", [round(med[f"acf_{L}min"], 3) for L in (1, 5, 8, 12)],
                  "| z", [round(z[f"acf_{L}min"], 1) for L in (1, 2, 3, 5, 8, 12)], "| var_log_k z", round(z["var_log_k"], 1), flush=True)
    print("data acf 1/5/8/12 min", [round(dat[f"acf_{L}min"], 3) for L in (1, 5, 8, 12)])
    json.dump(res, open("results/step0/rate_acf.json", "w"), indent=1)


if __name__ == "__main__":
    main()
