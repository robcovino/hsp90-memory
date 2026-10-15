"""Calibration with Laplace intervals (no per-data-set MCMC). Same synthetic data sets as exact_coverage.py (parameter
draws from the data posterior with --seed, simulation seeds 5000 + i, simulation grid --sim-grid). For each data set:
maximum a posteriori in the logit sampling space (bayes_exact.Target, same priors and likelihood; L-BFGS from the true
value), numerical Hessian there, posterior approximated by N(u_map, (-H)^-1). Reported per parameter: coverage of
central intervals and the uniformity of Phi((u_true - u_map) / sd). Also the Laplace approximation for the real data,
compared with the MCMC posterior (bayes_<tag>.json), so that the test is interpretable for the MCMC results.

Usage: python scripts/exact_coverage_laplace.py ringF1_ou --couple eq --n-sets 40 --fit-grid 7 --sim-grid 15
Output: results/exact/coverage_laplace_<tag>.{json,npz}
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, sys, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm, kstest
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bayes_exact as be                                                     # noqa: E402
from hsp90memory import amm, io, simulators as sim                          # noqa: E402

LEVELS = np.array([0.5, 0.8, 0.9, 0.95])


def laplace(T, u0):
    f = lambda u: -T.logpost(u)[0] if np.isfinite(T.logpost(u)[0]) else 1e12
    r = minimize(f, u0, method="L-BFGS-B", options=dict(maxiter=600))
    H = be.hessian(lambda u: T.logpost(u)[0], r.x)
    w, V = np.linalg.eigh(-H); w = np.clip(w, 1e-6, None)
    cov = (V / w) @ V.T
    return r.x, np.sqrt(np.diag(cov)), float(-r.fun), int(r.nit)


def one(args):
    name, couple, sim_grid, fit_grid, u_true, seed = args
    t0 = time.time()
    Ts = be.Target(name, sim_grid, couple)
    X = sim.simulate_amm(Ts.model.Q(Ts.to_theta(Ts.phi(u_true))), Ts.model.classes, 216000, n_rep=4, seed=seed)
    be.DL = [io.dwells(x) for x in X]
    T = be.Target(name, fit_grid, couple)
    u_map, sd, lp, nit = laplace(T, u_true)
    z = (u_true - u_map) / sd
    return dict(z=z, sd=sd, nit=nit, time=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name"); ap.add_argument("--couple", default="OC"); ap.add_argument("--suffix", default="")
    ap.add_argument("--n-sets", type=int, default=40); ap.add_argument("--fit-grid", type=int, default=7)
    ap.add_argument("--sim-grid", type=int, default=15); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    is_ou = "_ou" in a.name
    tag = a.name + ({"OC": "", "both": "2", "eq": "eq"}[a.couple] if is_ou else "") + a.suffix
    out = Path(f"results/exact/coverage_laplace_{tag}.json")
    if out.exists(): sys.exit(f"{out} exists")
    z = np.load(f"results/exact/bayes_{tag}_samples.npz"); U = z["u"].reshape(-1, z["u"].shape[-1])
    rng = np.random.default_rng(a.seed); picks = U[rng.choice(len(U), a.n_sets, replace=False)]   # as exact_coverage.py
    t0 = time.time()
    # real data: Laplace vs MCMC
    be.DL = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
    T = be.Target(a.name, 9, a.couple)
    u_start = U[np.argmax(z["logpost"].reshape(-1))]
    u_map, sd, lp, nit = laplace(T, u_start)
    mc_mean, mc_sd = U.mean(0), U.std(0)
    data_cmp = {"laplace_mean": u_map.tolist(), "laplace_sd": sd.tolist(), "mcmc_mean": mc_mean.tolist(),
                "mcmc_sd": mc_sd.tolist(), "sd_ratio": (sd / mc_sd).tolist(), "mean_shift_in_mcmc_sd": ((u_map - mc_mean) / mc_sd).tolist()}
    print("real data: Laplace/MCMC sd ratio", np.round(sd / mc_sd, 2), flush=True)
    with get_context("fork").Pool(a.workers) as pool:
        res = pool.map(one, [(a.name, a.couple, a.sim_grid, a.fit_grid, picks[i], 5000 + i) for i in range(a.n_sets)], chunksize=1)
    Z = np.array([r["z"] for r in res]); d = Z.shape[1]
    names = [f"log_rate_{i}" for i in range(d - (2 if is_ou else 0))] + (["sigma_OC", "log_tau_w"] if is_ou else [])
    q = norm.ppf((1 + LEVELS) / 2)
    cov = {n: [float(np.mean(np.abs(Z[:, j]) <= qq)) for qq in q] for j, n in enumerate(names)}
    ks = {n: float(kstest(norm.cdf(Z[:, j]), "uniform").pvalue) for j, n in enumerate(names)}
    rep = {"tag": tag, "args": vars(a), "levels": LEVELS.tolist(), "coverage": cov, "ks_uniform_p": ks,
           "coverage_mean_over_params": np.mean(list(cov.values()), 0).tolist(), "real_data_laplace_vs_mcmc": data_cmp,
           "names": names, "time_s": time.time() - t0, "nit": [r["nit"] for r in res]}
    json.dump(rep, open(out, "w"), indent=1)
    np.savez(f"results/exact/coverage_laplace_{tag}.npz", z=Z, picks=picks, names=np.array(names))
    print(json.dumps({k: rep[k] for k in ("coverage_mean_over_params", "time_s")}))
    for n in names: print(f"  {n:12s} {np.round(cov[n], 2)}  KS p {ks[n]:.3f}")


if __name__ == "__main__":
    main()
