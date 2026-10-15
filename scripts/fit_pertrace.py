"""Audit step 0, part 2: is the modulation the same in every molecule? Final model (ringF1 + equilibrium OU on kOC)
with the gating rates shared and
  shared  : one sigma_OC, one tau_w (the final model, 26 parameters)
  sigma4  : sigma_OC per trace, tau_w shared (29)
  st4     : sigma_OC and tau_w per trace (32)
  off1    : no modulation in trace 1, one sigma_OC and tau_w for traces 2-4 (26)
ML (L-BFGS from the final optimum) and Laplace evidence in the logit space of bayes_exact.py (same priors: rates
log-uniform on [1e-7, 1e2] s^-1, sigma U[0, 2], tau_w log-U[1, 1e4] s); the shared model's Laplace value is compared
with its bridge-sampling value. Per-trace log-likelihoods are reported.

Usage: python scripts/fit_pertrace.py [--grid 9] [--workers 4]
Output: results/step0/pertrace.json
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import argparse, json, sys, time
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bayes_exact import RATE_BOX, SIGMA_BOX, LOGTAU_BOX, hessian      # noqa: E402
from hsp90memory import amm                                            # noqa: E402

NB = 24
VARIANTS = {"shared": ([0, 0, 0, 0], [0, 0, 0, 0]), "sigma4": ([0, 1, 2, 3], [0, 0, 0, 0]),
            "st4": ([0, 1, 2, 3], [0, 1, 2, 3]), "off1": ([None, 0, 0, 0], [None, 0, 0, 0])}


class PerTrace:
    def __init__(self, variant, grid):
        self.m = amm.modal_ou(0, n_grid=grid, couple="eq", topology="ringF1")
        self.dl = [amm.load_dwells(f"data/{k}Thermal.txt")[1:] for k in "1234"]
        self.si, self.ti = VARIANTS[variant]
        self.ns = len({i for i in self.si if i is not None}); self.nt = len({i for i in self.ti if i is not None})
        self.d = NB + self.ns + self.nt
        self.lo = np.r_[np.full(NB, RATE_BOX[0]), np.full(self.ns, SIGMA_BOX[0]), np.full(self.nt, LOGTAU_BOX[0])]
        self.hi = np.r_[np.full(NB, RATE_BOX[1]), np.full(self.ns, SIGMA_BOX[1]), np.full(self.nt, LOGTAU_BOX[1])]
        self.base = amm.named_graph("ringF1")

    def phi(self, u): return self.lo + (self.hi - self.lo) * expit(u)
    def u(self, phi): return logit(np.clip((phi - self.lo) / (self.hi - self.lo), 1e-12, 1 - 1e-12))

    def per_trace_ll(self, u):
        ph = self.phi(u); rates = ph[:NB]; sig = ph[NB:NB + self.ns]; ltau = ph[NB + self.ns:]
        out = []
        for k in range(4):
            if self.si[k] is None:
                Q = self.base.Q(rates)
            else:
                Q = self.m.Q(np.r_[rates, np.log(max(sig[self.si[k]], 1e-300)), ltau[self.ti[k]]])
            out.append(amm.loglik_Q(Q, self.base.classes if self.si[k] is None else self.m.classes, [self.dl[k]]))
        return np.array(out)

    def logprior(self, u): return float(np.sum(-np.logaddexp(0, -u) - np.logaddexp(0, u)))

    def logpost(self, u):
        ll = self.per_trace_ll(u).sum(); return ll + self.logprior(u) if np.isfinite(ll) else -np.inf


def fit(args):
    variant, grid = args
    T = PerTrace(variant, grid); t0 = time.time()
    th = np.array(json.load(open("results/exact/fit_ringF1_oueq_global.json"))["theta"])
    phi0 = np.r_[th[:NB], np.full(T.ns, np.exp(th[NB])), np.full(T.nt, th[NB + 1])]
    f = lambda u: (lambda v: -v if np.isfinite(v) else 1e12)(T.logpost(u))
    r = minimize(f, T.u(phi0), method="L-BFGS-B", options=dict(maxiter=400))
    u_hat = r.x; lp = -r.fun; ll_tr = T.per_trace_ll(u_hat)
    H = hessian(T.logpost, u_hat); w = np.linalg.eigvalsh(-H)
    w_c = np.clip(w, 1.0 / 4.0, None)                     # flat directions capped (proposal sd 2 in logit units)
    logZ = lp + 0.5 * T.d * np.log(2 * np.pi) - 0.5 * np.sum(np.log(w_c))
    ph = T.phi(u_hat)
    return dict(variant=variant, grid=grid, d=T.d, logpost_max=lp, ll=float(ll_tr.sum()), ll_per_trace=ll_tr.tolist(),
                logZ_laplace=float(logZ), n_capped=int(np.sum(w < 0.25)), nit=int(r.nit),
                sigma=ph[NB:NB + T.ns].tolist(), tau_w=np.exp(ph[NB + T.ns:]).tolist(), time=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--grid", type=int, default=9); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(); out = Path("results/step0/pertrace.json"); out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists(): sys.exit(f"{out} exists")
    with get_context("fork").Pool(a.workers) as pool:
        res = pool.map(fit, [(v, a.grid) for v in VARIANTS], chunksize=1)
    ref = next(r for r in res if r["variant"] == "shared")
    bridge = json.load(open("results/exact/bayes_ringF1_oueq.json"))["logZ"]
    for r in res:
        r["dlogZ_vs_shared"] = r["logZ_laplace"] - ref["logZ_laplace"]; r["dll_vs_shared"] = r["ll"] - ref["ll"]
        print(f"{r['variant']:7s} d={r['d']} ll {r['ll']:.2f} (d {r['dll_vs_shared']:+.2f})  logZ_L {r['logZ_laplace']:.2f} "
              f"(d {r['dlogZ_vs_shared']:+.2f})  sigma {np.round(r['sigma'], 3)}  tau {np.round(r['tau_w'], 0)}  "
              f"per-trace ll {np.round(r['ll_per_trace'], 1)}", flush=True)
    json.dump({"grid": a.grid, "bridge_logZ_shared_g9": bridge, "results": res}, open(out, "w"), indent=1)
    print("shared: Laplace", round(ref["logZ_laplace"], 2), "vs bridge (g9)", round(bridge, 2))


if __name__ == "__main__":
    main()
