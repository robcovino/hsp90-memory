"""Bayesian inference with the exact likelihood for discrete modal-K AMMs and grid-OU modal-K models:
posterior by adaptive-Metropolis MCMC, marginal likelihood (evidence) by bridge sampling.

Priors (data-independent; box-uniform in the sampling parameters phi):
  every free rate: log-uniform on [1e-7, 1e2] s^-1;   sigma_OC: uniform on [0, 2];   tau_w: log-uniform on [1, 1e4] s.
Sampling space: u = logit((phi - lo) / (hi - lo)), prior density in u = prod s(u)(1 - s(u)).
MCMC: --chains independent chains in parallel processes, each started at the ML point + small jitter, Gaussian random-
walk proposal starting from the numerical Hessian at the mode, adapted during the first half (burn-in, discarded;
Andrieu & Thoms 2008: running covariance with clipped eigenvalues + Robbins-Monro global scale, target acceptance
0.234); the second half uses the frozen proposal. (Runs before 2026-10-13 03:00 UTC used plain Haario adaptation.) Diagnostics: split-R-hat, integrated
autocorrelation time and ESS (emcee.autocorr), acceptance.
Evidence: bridge sampling (Meng & Wong 1996; iterative scheme, Gronau et al. 2017) with a multivariate normal proposal
fitted to the first half of the post-burn-in samples and evaluated against the second half; reported per chain
(spread = error estimate) and pooled. log Z includes the conditioning on the first-frame class (as the likelihood).
Grid-OU models are sampled on a coarse grid (--grid); a posterior-averaged correction to finer grids
(--grid-correct) is applied to log Z: dlogZ = log E_post[exp(ll_fine - ll_coarse)].

Usage: python scripts/bayes_exact.py modal5 --steps 40000
       python scripts/bayes_exact.py modal4_ou --grid 9 --grid-correct 15 21 --steps 40000 [--couple eq]
Outputs: results/exact/bayes_<name>.json, bayes_<name>_samples.npz, bayes_<name>.log
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, re, sys, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.special import expit, logit, logsumexp
from hsp90memory import amm

OUT = Path("results/exact")
RATE_BOX = (np.log(1e-7), np.log(1e2)); SIGMA_BOX = (0.0, 2.0); LOGTAU_BOX = (0.0, np.log(1e4))
DL = None


def build(name, grid, couple):
    m = re.match(r"modal(\d+)(_ou)?$", name); K = int(m.group(1))
    if m.group(2):
        model = amm.modal_ou(K, n_grid=grid, couple=couple); nb = model.base.n_par
        lo = np.r_[np.full(nb, RATE_BOX[0]), SIGMA_BOX[0], LOGTAU_BOX[0]]
        hi = np.r_[np.full(nb, RATE_BOX[1]), SIGMA_BOX[1], LOGTAU_BOX[1]]
        def to_theta(phi):        # phi: log-rates, sigma, log tau -> model theta: log-rates, log sigma, log tau
            return np.r_[phi[:nb], np.log(max(phi[nb], 1e-300)), phi[nb + 1]]
        def to_phi(theta):
            return np.r_[theta[:nb], np.exp(theta[nb]), theta[nb + 1]]
    else:
        model = amm.modal(K); nb = model.n_par
        lo, hi = np.full(nb, RATE_BOX[0]), np.full(nb, RATE_BOX[1])
        to_theta = to_phi = lambda x: np.asarray(x, float)
    return model, lo, hi, to_theta, to_phi


class Target:
    def __init__(self, name, grid, couple):
        self.model, self.lo, self.hi, self.to_theta, self.to_phi = build(name, grid, couple)
        self.d = len(self.lo)

    def phi(self, u):
        return self.lo + (self.hi - self.lo) * expit(u)

    def u(self, phi):
        return logit(np.clip((phi - self.lo) / (self.hi - self.lo), 1e-12, 1 - 1e-12))

    def loglik(self, u):
        ll = amm.loglik_Q(self.model.Q(self.to_theta(self.phi(u))), self.model.classes, DL)
        return ll if np.isfinite(ll) else -np.inf

    def logprior(self, u):
        return float(np.sum(-np.logaddexp(0, -u) - np.logaddexp(0, u)))     # log s(u) + log(1 - s(u))

    def logpost(self, u):
        lp = self.logprior(u)
        return lp + self.loglik(u), lp


def hessian(f, x, h=1e-3):
    d = len(x); H = np.zeros((d, d)); f0 = f(x)
    for i in range(d):
        for j in range(i, d):
            e_i = np.zeros(d); e_i[i] = h; e_j = np.zeros(d); e_j[j] = h
            if i == j:
                H[i, i] = (f(x + e_i) - 2 * f0 + f(x - e_i)) / h ** 2
            else:
                H[i, j] = H[j, i] = (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4 * h * h)
    return H


def _clip_cov(C, lo=1e-6, hi=4.0):
    w, V = np.linalg.eigh((C + C.T) / 2); return (V * np.clip(w, lo, hi)) @ V.T


def chain(args):
    """Adaptive Metropolis (Andrieu & Thoms 2008, alg. 4): during burn-in the proposal covariance is the running
    covariance of the chain (eigenvalues clipped to [1e-6, 4] in logit units, i.e. sd <= 2) and a global log-scale is
    tuned by Robbins-Monro towards acceptance 0.234; both are frozen for the sampling half."""
    name, grid, couple, u0, cov0, steps, seed = args
    T = Target(name, grid, couple); rng = np.random.default_rng(seed); d = T.d
    x = u0 + rng.multivariate_normal(np.zeros(d), 0.1 * cov0); lp, lpr = T.logpost(x)
    sc = 2.38 ** 2 / d; C = _clip_cov(cov0); L = np.linalg.cholesky(sc * C + 1e-10 * np.eye(d)); log_s = 0.0
    burn = steps // 2; X = np.empty((steps - burn, d)); LP = np.empty(steps - burn); LL = np.empty(steps - burn)
    acc = 0; mean = x.copy(); M2 = np.zeros((d, d)); n_ad = 0; t0 = time.time()
    for t in range(steps):
        y = x + np.exp(log_s) * (L @ rng.standard_normal(d))
        ly, lyr = T.logpost(y)
        a_prob = np.exp(min(0.0, ly - lp)) if np.isfinite(ly) else 0.0
        if rng.random() < a_prob:
            x, lp, lpr = y, ly, lyr; acc += t >= burn
        if t < burn:
            log_s += (t + 1) ** -0.6 * (a_prob - 0.234)
            n_ad += 1; dx = x - mean; mean += dx / n_ad; M2 += np.outer(dx, x - mean)
            if t >= 1000 and t % 200 == 0:
                C = _clip_cov(M2 / (n_ad - 1)); L = np.linalg.cholesky(sc * C + 1e-10 * np.eye(d))
        else:
            X[t - burn] = x; LP[t - burn] = lp; LL[t - burn] = lp - lpr
    return dict(X=X, LP=LP, LL=LL, acc=acc / (steps - burn), log_scale=log_s, time=time.time() - t0)


def bridge(lq1, lg1, lq2, lg2, tol=1e-10, maxit=1000):
    """Meng & Wong iterative bridge sampling. lq: log unnormalised posterior, lg: log proposal density;
    1 = posterior samples, 2 = proposal samples. Returns log Z."""
    n1, n2 = len(lq1), len(lq2); s1, s2 = n1 / (n1 + n2), n2 / (n1 + n2)
    l1 = lq1 - lg1; l2 = lq2 - lg2; lstar = np.median(l1)
    r = 0.0
    for _ in range(maxit):
        num = logsumexp(l2 - lstar - np.logaddexp(np.log(s1) + l2 - lstar, np.log(s2) + r)) - np.log(n2)
        den = logsumexp(-np.logaddexp(np.log(s1) + l1 - lstar, np.log(s2) + r)) - np.log(n1)
        r_new = num - den
        if abs(r_new - r) < tol: r = r_new; break
        r = r_new
    return r + lstar


def mvn_logpdf(X, mu, L):
    z = np.linalg.solve(L, (X - mu).T); d = len(mu)
    return -0.5 * np.sum(z * z, 0) - np.sum(np.log(np.diag(L))) - 0.5 * d * np.log(2 * np.pi)


def eval_logpost(args):
    name, grid, couple, U = args
    T = Target(name, grid, couple)
    return np.array([T.logpost(u)[0] for u in U])


def eval_ll_grids(args):
    name, couple, grids, U = args
    out = []
    for g in grids:
        T = Target(name, g, couple); out.append([T.loglik(u) for u in U])
    return np.array(out)


def main():
    global DL
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--couple", default="OC")
    ap.add_argument("--grid", type=int, default=9)
    ap.add_argument("--grid-correct", type=int, nargs="*", default=[15, 21])
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--n-bridge", type=int, default=4000, help="proposal draws per chain for bridge sampling")
    ap.add_argument("--n-grid-corr", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--start", default=None, help="JSON with 'theta' (model parametrisation) for the start point")
    ap.add_argument("--suffix", default="", help="appended to the output tag (e.g. _long for a longer rerun)")
    a = ap.parse_args()
    is_ou = "_ou" in a.name
    base_tag = a.name + ({"OC": "", "both": "2", "eq": "eq"}[a.couple] if is_ou else "")
    tag = base_tag + a.suffix
    OUT.mkdir(parents=True, exist_ok=True); out = OUT / f"bayes_{tag}.json"
    if out.exists():
        sys.exit(f"{out} exists")
    logf = open(OUT / f"bayes_{tag}.log", "a")
    def log(msg):
        print(msg, flush=True); logf.write(msg + "\n"); logf.flush()
    DL = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
    t0 = time.time()
    T = Target(a.name, a.grid, a.couple)
    if a.start:
        start = a.start
    elif is_ou:
        start = f"results/exact/fit_{base_tag}.json"
    else:
        start = str(next(q for q in (Path(f"results/exact/fit_{a.name}_global.json"), Path(f"results/{a.name}.json"),
                                     Path(f"results/exact/fit_{a.name}.json")) if q.exists()))
    th = np.array(json.load(open(start))["theta"]); u_hat = T.u(T.to_phi(th))
    f = lambda u: T.logpost(u)[0]
    log(f"{tag}: d={T.d}, start {start}, logpost at start {f(u_hat):.3f}")
    H = hessian(f, u_hat); w, V = np.linalg.eigh(-H)
    w = np.clip(w, 1.0 / 4.0, None)                    # cap proposal sd at 2 in logit units for flat directions
    cov0 = (V / w) @ V.T
    log(f"  Hessian done ({time.time() - t0:.0f} s); smallest curvature {np.linalg.eigvalsh(-H).min():.3g}")
    with get_context("fork").Pool(a.chains) as pool:
        res = pool.map(chain, [(a.name, a.grid, a.couple, u_hat, cov0, a.steps, a.seed * 100 + c) for c in range(a.chains)])
    log(f"  MCMC done ({time.time() - t0:.0f} s); acceptance {[round(r['acc'], 3) for r in res]}")
    Xs = np.array([r["X"] for r in res]); LPs = np.array([r["LP"] for r in res]); LLs = np.array([r["LL"] for r in res])
    import emcee
    n = Xs.shape[1]; half = n // 2
    split = np.concatenate([Xs[:, :half], Xs[:, half:2 * half]], 0)              # split-R-hat
    W = split.var(1, ddof=1).mean(0); B = half * split.mean(1).var(0, ddof=1)
    rhat = np.sqrt(((half - 1) / half * W + B / half) / W)
    tau = np.array([[emcee.autocorr.integrated_time(Xs[c][:, j], quiet=True)[0] for j in range(T.d)]
                    for c in range(a.chains)])
    ess = (n / tau).sum(0)
    log(f"  max R-hat {rhat.max():.3f}; min ESS {ess.min():.0f}; max tau {tau.max():.0f}")

    # bridge sampling per chain: proposal from the first half (thinned), evaluated with the second half
    thin = max(1, int(np.ceil(np.median(tau))))
    jobs, props = [], []
    rng = np.random.default_rng(a.seed + 7)
    for c in range(a.chains):
        fit = Xs[c, :half]; mu = fit.mean(0); L = np.linalg.cholesky(np.cov(fit.T) + 1e-10 * np.eye(T.d))
        U2 = mu + rng.standard_normal((a.n_bridge, T.d)) @ L.T
        props.append((mu, L, U2)); jobs.append((a.name, a.grid, a.couple, U2))
    with get_context("fork").Pool(a.chains) as pool:
        lq2s = pool.map(eval_logpost, jobs)
    logZ_chain = []
    for c in range(a.chains):
        mu, L, U2 = props[c]
        U1 = Xs[c, half::thin]; lq1 = LPs[c, half::thin]
        logZ_chain.append(bridge(lq1, mvn_logpdf(U1, mu, L), lq2s[c], mvn_logpdf(U2, mu, L)))
    logZ_chain = np.array(logZ_chain)
    logZ = float(logsumexp(logZ_chain) - np.log(a.chains)); logZ_se = float(logZ_chain.std(ddof=1) / np.sqrt(a.chains))
    log(f"  bridge log Z per chain {np.round(logZ_chain, 3)}; pooled {logZ:.3f} +- {logZ_se:.3f} ({time.time() - t0:.0f} s)")

    corr = {}
    if is_ou and a.grid_correct:
        flat = Xs.reshape(-1, T.d); sel = flat[np.random.default_rng(3).choice(len(flat), a.n_grid_corr, replace=False)]
        chunks = np.array_split(sel, a.chains)
        with get_context("fork").Pool(a.chains) as pool:
            parts = pool.map(eval_ll_grids, [(a.name, a.couple, [a.grid] + a.grid_correct, ch) for ch in chunks])
        LLg = np.concatenate(parts, 1)
        for i, g in enumerate(a.grid_correct):
            dl_ = LLg[i + 1] - LLg[0]
            corr[str(g)] = {"dlogZ": float(logsumexp(dl_) - np.log(len(dl_))), "dll_mean": float(dl_.mean()),
                            "dll_sd": float(dl_.std())}
        log(f"  grid correction: {corr}")

    flat = Xs.reshape(-1, T.d); PH = np.array([T.phi(u) for u in flat])
    names = [f"p{i}" for i in range(T.d)]
    q = np.quantile(PH, [0.05, 0.5, 0.95], 0)
    res_js = {"name": tag, "model": T.model.name, "d": T.d, "args": vars(a),
              "priors": {"rates_log_uniform": [1e-7, 1e2], "sigma_uniform": list(SIGMA_BOX), "tau_log_uniform": [1, 1e4]},
              "start": start, "acceptance": [r["acc"] for r in res], "rhat": rhat.tolist(), "ess": ess.tolist(),
              "tau_int": tau.tolist(), "logZ_chain": logZ_chain.tolist(), "logZ": logZ, "logZ_se": logZ_se,
              "grid_correction": corr,
              "logZ_corrected": {g: logZ + v["dlogZ"] for g, v in corr.items()},
              "ll_max_sampled": float(LLs.max()), "ll_post_mean": float(LLs.mean()),
              "phi_q05": q[0].tolist(), "phi_median": q[1].tolist(), "phi_q95": q[2].tolist(),
              "phi_lo": T.lo.tolist(), "phi_hi": T.hi.tolist(), "time_s": time.time() - t0}
    json.dump(res_js, open(out, "w"), indent=1)
    np.savez(OUT / f"bayes_{tag}_samples.npz", u=Xs, logpost=LPs, loglik=LLs, phi=PH.reshape(Xs.shape))
    log(f"done {tag}: log Z = {logZ:.3f} +- {logZ_se:.3f}; time {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
