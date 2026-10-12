"""Calibration of the exact-likelihood posteriors (coverage / rank test in the region of the data).

Parameter sets are drawn from the posterior of the real data (results/exact/bayes_<tag>_samples.npz); for each, a
synthetic measurement (4 traces x 216,000 frames) is simulated exactly from the generator (simulators.simulate_amm;
grid-OU models on --sim-grid points) and the posterior is sampled again with the same priors and likelihood
(MCMC, one chain, random-walk proposal with the covariance of the real-data posterior in logit space, started at the
true value plus jitter, first half discarded; grid-OU models on --fit-grid points, which also tests the grid
approximation). Output: rank of the true value among thinned posterior samples, per parameter, and the empirical
coverage of central credible intervals. Under correct calibration ranks are uniform and coverage is nominal.

Usage: python scripts/exact_coverage.py modal5 --n-sets 100 --steps 8000
       python scripts/exact_coverage.py modal4_ou --couple eq --n-sets 60 --steps 8000 --fit-grid 7 --sim-grid 15
Outputs: results/exact/coverage_<tag>.json, .npz; figures/exact_coverage_<tag>.png
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, sys, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bayes_exact as be                                                     # noqa: E402
from hsp90memory import amm, io, simulators as sim                          # noqa: E402

LEVELS = np.array([0.5, 0.8, 0.9, 0.95])


def one(args):
    name, couple, sim_grid, fit_grid, u_true, cov, steps, seed = args
    rng = np.random.default_rng(seed)
    Ts = be.Target(name, sim_grid, couple)
    Q = Ts.model.Q(Ts.to_theta(Ts.phi(u_true)))
    X = sim.simulate_amm(Q, Ts.model.classes, 216000, n_rep=4, seed=seed)
    be.DL = [io.dwells(x) for x in X]
    T = be.Target(name, fit_grid, couple); d = T.d
    L = np.linalg.cholesky(2.38 ** 2 / d * cov + 1e-10 * np.eye(d))
    x = u_true + L @ rng.standard_normal(d) * 0.3; lp = T.logpost(x)[0]
    burn = steps // 2; S = []; acc = 0; t0 = time.time()
    for t in range(steps):
        y = x + L @ rng.standard_normal(d); ly = T.logpost(y)[0]
        if np.log(rng.random()) < ly - lp:
            x, lp = y, ly; acc += t >= burn
        if t >= burn: S.append(x.copy())
    S = np.array(S)
    import emcee
    tau = np.array([emcee.autocorr.integrated_time(S[:, j], quiet=True)[0] for j in range(d)])
    thin = max(1, int(np.ceil(np.median(tau)))); St = S[::thin]
    rank = (St < u_true).mean(0)                                # fraction of posterior samples below the truth
    return dict(rank=rank, n_eff=len(St), acc=acc / (steps - burn), tau=tau, time=time.time() - t0,
                n_dwells=int(sum(len(d_[0]) for d_ in be.DL)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name"); ap.add_argument("--couple", default="OC")
    ap.add_argument("--n-sets", type=int, default=100); ap.add_argument("--steps", type=int, default=8000)
    ap.add_argument("--fit-grid", type=int, default=7); ap.add_argument("--sim-grid", type=int, default=15)
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    is_ou = "_ou" in a.name
    tag = a.name + ({"OC": "", "both": "2", "eq": "eq"}[a.couple] if is_ou else "")
    out = Path(f"results/exact/coverage_{tag}.json")
    if out.exists():
        sys.exit(f"{out} exists")
    z = np.load(f"results/exact/bayes_{tag}_samples.npz"); U = z["u"].reshape(-1, z["u"].shape[-1])
    cov = np.cov(U.T); rng = np.random.default_rng(a.seed)
    picks = U[rng.choice(len(U), a.n_sets, replace=False)]
    t0 = time.time()
    with get_context("fork").Pool(a.workers) as pool:
        res = pool.map(one, [(a.name, a.couple, a.sim_grid, a.fit_grid, picks[i], cov, a.steps, 5000 + i)
                             for i in range(a.n_sets)], chunksize=1)
    R = np.array([r["rank"] for r in res]); d = R.shape[1]
    T = be.Target(a.name, a.fit_grid, a.couple)
    names = [f"log_rate_{i}" for i in range(T.model.base.n_par if is_ou else T.model.n_par)] + (["sigma_OC", "log_tau_w"] if is_ou else [])
    cov_emp = {n: [float(np.mean(np.abs(R[:, j] - 0.5) <= l / 2)) for l in LEVELS] for j, n in enumerate(names)}
    from scipy.stats import kstest
    ks = {n: float(kstest(R[:, j], "uniform").pvalue) for j, n in enumerate(names)}
    rep = {"tag": tag, "args": vars(a), "n_sets": a.n_sets, "levels": LEVELS.tolist(), "coverage": cov_emp,
           "ks_uniform_p": ks, "acceptance": [r["acc"] for r in res], "n_eff": [r["n_eff"] for r in res],
           "n_dwells": [r["n_dwells"] for r in res], "time_s": time.time() - t0,
           "coverage_mean_over_params": np.mean([v for v in cov_emp.values()], 0).tolist()}
    json.dump(rep, open(out, "w"), indent=1)
    np.savez(f"results/exact/coverage_{tag}.npz", ranks=R, picks=picks, names=np.array(names))
    nc = 5; nr = int(np.ceil(d / nc)) + 1
    fig, axs = plt.subplots(nr, nc, figsize=(16, 2.4 * nr))
    for j, ax in enumerate(axs.ravel()):
        if j < d:
            ax.hist(R[:, j], bins=10, range=(0, 1), color="C0"); ax.set_title(f"{names[j]} KS p={ks[names[j]]:.2f}", fontsize=8)
            ax.axhline(len(R) / 10, color="k", lw=0.6); ax.set_yticks([])
        elif j == d:
            for jj, n in enumerate(names):
                ax.plot(LEVELS, cov_emp[n], "-", color="C3" if n in ("sigma_OC", "log_tau_w") else "0.6", lw=1)
            ax.plot([0.45, 1], [0.45, 1], "k:"); ax.set_title("coverage (red: sigma, tau)", fontsize=8)
        else:
            ax.axis("off")
    fig.suptitle(f"{tag}: ranks of true values, {a.n_sets} synthetic measurements from the data posterior", fontsize=10)
    fig.tight_layout(); fig.savefig(f"figures/exact_coverage_{tag}.png", dpi=100)
    print(json.dumps({k: rep[k] for k in ("coverage_mean_over_params", "time_s")}), flush=True)
    print({n: [round(x, 2) for x in v] for n, v in cov_emp.items()})


if __name__ == "__main__":
    main()
