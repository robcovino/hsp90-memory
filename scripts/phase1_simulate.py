"""Phase 1: checkpointed simulation driver for NPE over the nested modal-OU model.

Prior (settings embedded in every batch file and in <out>/settings.json):
  - the 12 free log-rates of modal3 (results/modal3.json, theta): independent normal around the ML values with sd ln 2
    (log-normal, factor 2); Q rebuilt with amm.modal(3), so detailed balance holds for every draw;
  - sigma_OC, sigma_CO ~ U[0, 1]; tau_w ~ log-U[1, 300] s; rho ~ U[-1, 1].
Each simulation is one measurement of 4 x 6 h at 10 Hz (simulators.simulate_modal_ou_nb), reduced to
summaries.summary_vector. Batches of --batch-size simulations go to <out>/batch_XXXX.npz (written atomically).
Parameters and simulator seeds of batch b depend only on (--seed, b), so a rerun skips every complete batch file and
produces exactly the missing ones: kill and restart at any time.

Usage: python scripts/phase1_simulate.py --n-batches 10 [--workers 4] [--out results/phase1]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, json, logging, resource, sys, time, warnings
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from hsp90memory import amm, simulators as sim, summaries as sm

NF, NT = 216000, 4
PARAM_NAMES = [f"log_{n}" for n in ("kOC0", "kCO0", "kOC1", "kCO1", "kOC2", "kCO2",
                                    "a0", "b0", "c0", "a1", "b1", "c1")] + ["sigma_OC", "sigma_CO", "tau_w", "rho"]


def settings_from_args(a):
    r = json.load(open(a.amm))
    return {"amm": a.amm, "theta_hat": r["theta"], "theta_sd": float(np.log(2.0)), "model": r["model"],
            "sigma_OC": [0.0, 1.0], "sigma_CO": [0.0, 1.0], "tau_w_loguniform": [1.0, 300.0], "rho": [-1.0, 1.0],
            "batch_size": a.batch_size, "seed": a.seed, "n_frames": NF, "n_traces": NT, "dt": 0.1,
            "simulator": "simulate_modal_ou_nb", "summaries": "summary_vector", "param_names": PARAM_NAMES}


def sample_prior(S, b):
    rng = np.random.default_rng([S["seed"], b])
    n = S["batch_size"]; th = np.asarray(S["theta_hat"])
    theta = th + S["theta_sd"] * rng.standard_normal((n, th.size))
    sOC = rng.uniform(*S["sigma_OC"], n); sCO = rng.uniform(*S["sigma_CO"], n)
    tau = np.exp(rng.uniform(*np.log(S["tau_w_loguniform"]), n)); rho = rng.uniform(*S["rho"], n)
    seeds = rng.integers(0, 2 ** 32, n, dtype=np.uint64)
    return np.column_stack([theta, sOC, sCO, tau, rho]), seeds


def batch_path(out, b):
    return Path(out) / f"batch_{b:04d}.npz"


def is_complete(path, n):
    try:
        with np.load(path) as f:
            return f["x"].shape[0] == n and f["params"].shape[0] == n
    except Exception:
        return False


def run_batch(job):
    b, S, out = job
    log = logging.getLogger("phase1")
    model = amm.modal(3)
    P, seeds = sample_prior(S, b); n = len(P)
    X, T = [], np.empty(n)
    t_batch = time.time()
    for i in range(n):
        t0 = time.perf_counter()
        theta, (sOC, sCO, tau, rho) = P[i, :12], P[i, 12:]
        kOC, kCO, SO, SC = sim.modal_params_from_Q(model.Q(theta))
        Y = sim.simulate_modal_ou_nb(kOC, kCO, SO, SC, (sOC, sCO), tau, rho, NF, n_rep=NT, seed=int(seeds[i]))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            v, names = sm.summary_vector(list(Y))
        X.append(v); T[i] = time.perf_counter() - t0
        if (i + 1) % 200 == 0:
            log.info(f"batch {b:04d}: {i + 1}/{n}, {T[:i + 1].mean():.3f} s/sim")
    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    path = batch_path(out, b); tmp = path.with_suffix(".tmp")
    with open(tmp, "wb") as f:
        np.savez(f, params=P, param_names=np.array(PARAM_NAMES), x=np.array(X), names=np.array(names),
                 seeds=seeds, sim_time=T, batch=b, peak_rss_mb=rss_mb, settings=json.dumps(S))
    os.replace(tmp, path)
    dt = time.time() - t_batch
    log.info(f"batch {b:04d} done: {dt:.0f} s, {T.mean():.3f} s/sim, peak RSS {rss_mb:.0f} MB, "
             f"non-finite summaries {int((~np.isfinite(np.array(X))).sum())}")
    return b, dt, rss_mb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-batches", type=int, required=True, help="total batches wanted (existing ones are kept)")
    ap.add_argument("--batch-size", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--out", default="results/phase1")
    ap.add_argument("--amm", default="results/modal3.json")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(processName)s %(message)s",
                        handlers=[logging.FileHandler(out / "phase1.log"), logging.StreamHandler(sys.stdout)])
    log = logging.getLogger("phase1")
    S = settings_from_args(a); sfile = out / "settings.json"
    if sfile.exists():
        old = json.load(open(sfile))
        if old != S:
            sys.exit(f"{sfile} differs from the requested settings; use another --out (results are never overwritten)")
    else:
        json.dump(S, open(sfile, "w"), indent=1)
    for p in out.glob("batch_*.tmp"):
        p.unlink()                                   # partial writes of an interrupted run
    todo = [b for b in range(a.n_batches) if not is_complete(batch_path(out, b), a.batch_size)]
    log.info(f"{a.n_batches - len(todo)} complete batches found, {len(todo)} to run on {a.workers} workers")
    if not todo:
        return
    sim.simulate_modal_ou_nb(np.ones(2), np.ones(2), np.zeros((2, 2)), np.zeros((2, 2)), 0.5, 10.0, 0.0, 10)  # compile
    t0 = time.time()
    with get_context("fork").Pool(a.workers, maxtasksperchild=1) as pool:
        for b, dt, rss in pool.imap_unordered(run_batch, [(b, S, str(out)) for b in todo]):
            pass
    log.info(f"all {len(todo)} batches done in {time.time() - t0:.0f} s wall")


if __name__ == "__main__":
    main()
